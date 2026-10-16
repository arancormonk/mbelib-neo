#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the DVSI scoreboard's measurements and gates on synthetic data (no DVSI material)."""

import importlib.util
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile

SKIP = 77  # CTest SKIP_RETURN_CODE: the pitch measurements need NumPy


def load_scoreboard():
    path = Path(__file__).resolve().parents[1] / "tools" / "quality" / "dvsi_scoreboard.py"
    spec = importlib.util.spec_from_file_location("dvsi_scoreboard", path)
    assert spec is not None and spec.loader is not None, path
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def harmonic_signal(f0, seconds=1.0, harmonics_to=3600.0, seed=1):
    rng = random.Random(seed)
    count = int(harmonics_to / f0)
    phases = [rng.uniform(0, 2 * math.pi) for _ in range(count)]
    return [
        sum(3000.0 / k * math.cos(2 * math.pi * k * f0 * n / 8000.0 + phases[k - 1]) for k in range(1, count + 1))
        + rng.gauss(0, 20.0)
        for n in range(int(8000 * seconds))
    ]


def check_pitch(board):
    for f0 in (83.7, 123.4, 212.9, 331.0):
        signal = harmonic_signal(f0)
        # A guess off by up to 4% (the size of the D-STAR mismatch) still finds the series.
        for guess_ratio in (1.0, 0.96, 1.04):
            measured = board.harmonic_f0(signal, 4000, f0 * guess_ratio)
            assert measured is not None, (f0, guess_ratio)
            assert abs(measured / f0 - 1.0) < 3e-4, (f0, guess_ratio, measured)
    rng = random.Random(7)
    noise = [rng.gauss(0, 3000.0) for _ in range(8000)]
    assert board.harmonic_f0(noise, 4000, 150.0) is None
    quiet = [0.001 * x for x in harmonic_signal(150.0)]
    assert board.harmonic_f0(quiet, 4000, 150.0) is None
    assert board.harmonic_f0(harmonic_signal(150.0), 100, 150.0) is None  # window leaves the signal


def record(w0_hz, L=20, voiced=None, flags=0, status=0):
    w0 = 2 * math.pi * w0_hz / 8000.0
    model = {"status": status, "w0": w0, "L": L, "Vl": voiced or [1] * L, "Ml": [1.0] * L, "log2Ml": [0.0] * L}
    return {"f": 0, "flags": flags, "errors": 0, "decoded": model, "history": model, "synth": model}


def check_frame_selection(board):
    # Status flags (C0/C4 valid, soft input) mark ordinary voice frames.
    records = [record(150.0, flags=0x02 | 0x04) for _ in range(9)]
    assert board.stable_voiced_frames(records) == [2, 3, 4, 5, 6]
    for handling in (0x10, 0x20, 0x40, 0x80, 0x100):  # tone, erasure, repeat, mute, silence
        records[4] = record(150.0, flags=0x02 | handling)  # breaks the +-2 frame window around it
        assert board.stable_voiced_frames(records) == [], hex(handling)
    records = [record(150.0) for _ in range(9)]
    records[6] = record(160.0)  # a 6.7% jump is not stable
    assert board.stable_voiced_frames(records) == [2, 3]
    records = [record(150.0, voiced=[1] * 10 + [0] * 10) for _ in range(9)]
    assert board.stable_voiced_frames(records) == []

    # A decoder that plays 3% sharp shows up as a 1.03 ratio against the coded pitch.
    records = [record(140.0, L=25) for _ in range(30)]
    sharp = harmonic_signal(140.0 * 1.03, seconds=30 * 160 / 8000.0)
    pairs = board.pitch_ratios(records, sharp, 0)
    assert len(pairs) >= 20
    summary = board.summarize_ratios(pairs)
    assert abs(summary["100_150"]["median"] - 1.03) < 1e-3, summary
    assert abs(board.max_ratio_deviation(summary) - 0.03) < 1e-3


def check_voicing(board):
    # 125 Hz: harmonics 1-7 below 1 kHz, 8-15 in 1-2 kHz, 16-23 in 2-3 kHz, 24-31 in 3-4 kHz.
    voiced = [1] * 15 + [0] * 16
    model = record(125.0, L=31, voiced=voiced)["decoded"]
    assert board.band_voicing(model) == [True, True, False, False]
    a = [record(125.0, L=31, voiced=voiced) for _ in range(4)]
    b = [record(125.0, L=31, voiced=[1] * 31) for _ in range(4)]
    assert board.voicing_agreement(a, b, 0) == [1.0, 1.0, 0.0, 0.0]
    assert board.voicing_agreement(a, b, 5) == [None] * 4


def vector(name, partition, lsd, env, rail=0):
    ours = {key: 0.0 for key in sorted(set(("level_offset_db", "lsd_db", "env_corr", "crest_delta_db",
                                            "pcm_rail_samples")))}
    ours.update({"lsd_db": lsd, "env_corr": env, "pcm_rail_samples": rail})
    pitch = {"all": {"n": 0, "median": None}}
    return {"mode": "dstar", "name": name, "partition": partition, "ours": ours, "dvsi": {}, "vs_dvsi": {},
            "pitch": {"ours": pitch, "dvsi": pitch, "input": pitch}}


def check_gates(board):
    names = [f"v{i}" for i in range(6)]
    baseline = {"vectors": [vector(n, "validation", 7.0, 0.88) for n in names]}
    better = {"vectors": [vector(n, "validation", 6.5 + 0.01 * i, 0.89, rail=0) for i, n in enumerate(names)]}
    mixed = {"vectors": [vector(n, "validation", 7.0 + (0.3 if i % 2 else -0.3), 0.88) for i, n in enumerate(names)]}
    gates = {"partition": "validation", "checks": [
        {"mode": "dstar", "metric": "ours.lsd_db", "rule": "decrease"},
        {"mode": "dstar", "metric": "ours.env_corr", "rule": "increase"},
        {"mode": "dstar", "metric": "ours.pcm_rail_samples", "rule": "no_increase"},
    ]}
    passed, lines = board.check_gates(gates, board.aggregate(better["vectors"]), board.compare(better, baseline))
    assert passed, lines
    passed, lines = board.check_gates(gates, board.aggregate(mixed["vectors"]), board.compare(mixed, baseline))
    assert not passed and any(line.startswith("FAIL dstar/validation ours.lsd_db") for line in lines), lines
    absolute = {"partition": "validation", "checks": [
        {"mode": "dstar", "metric": "ours.lsd_db", "rule": "abs_delta_max", "value": 0.0}]}
    assert board.check_gates(absolute, {}, board.compare(baseline, baseline))[0]
    assert not board.check_gates(absolute, {}, board.compare(better, baseline))[0]
    # Per-vector changes that cancel in the mean pass on the mean alone; "worst" catches them.
    cancelling = {"vectors": [vector(n, "validation", 7.0 + (0.5 if i == 0 else -0.5 if i == 1 else 0.0), 0.88)
                              for i, n in enumerate(names)]}
    assert board.check_gates(absolute, {}, board.compare(cancelling, baseline))[0]
    absolute["checks"][0]["worst"] = 0.0
    assert not board.check_gates(absolute, {}, board.compare(cancelling, baseline))[0]
    assert board.check_gates(absolute, {}, board.compare(baseline, baseline))[0]


def check_encoder_metrics(board):
    # Our encoding against DVSI's, both decoded by this library: the level gap is unsigned.
    result = vector("v0", "validation", 7.0, 0.88)
    result["encoder"] = {
        "speech": {"lsd_db": 5.0, "level_offset_db": 0.4},
        "vs_dvsi_bits": {"lsd_db": 4.2, "level_offset_db": -2.5, "env_corr": None},
        "rail_excess_samples": 3, "silence_frame_rate": 0.25, "dvsi_silence_frame_rate": 0.0,
        "pitch_cents_median": 6.0, "dvsi_pitch_cents_median": 9.0, "octave_error_rate": 0.0,
        "dvsi_octave_error_rate": 0.01, "voicing_agreement": [0.9, None, 0.7, 0.8],
    }
    flat = board.flatten(result)
    assert flat["encoder.vs_dvsi_bits.lsd_db"] == 4.2 and "encoder.vs_dvsi_bits.env_corr" not in flat, flat
    assert abs(flat["encoder.level_gap_db"] - 2.5) < 1e-12, flat
    assert flat["encoder.rail_excess_samples"] == 3 and flat["encoder.silence_frame_rate"] == 0.25, flat
    assert abs(flat["encoder.voicing_agreement_mean"] - 0.8) < 1e-12, flat
    del result["encoder"]["vs_dvsi_bits"], result["encoder"]["rail_excess_samples"]
    flat = board.flatten(result)  # scoreboards from before these metrics still flatten
    assert "encoder.level_gap_db" not in flat and "encoder.rail_excess_samples" not in flat, flat
    # A comparison pinned at the evaluator's lag search limit, or with no alignment at all (a silent
    # decode), is rejected, not scored.
    good = {"alignment_at_limit": False, "alignment_corr": 0.97, "lsd_db": 4.0}
    assert board.aligned(good, "x")["lsd_db"] == 4.0
    for bad in ({"alignment_at_limit": True, "alignment_corr": 0.9, "auto_lag_samples": -160},
                {"alignment_at_limit": False, "alignment_corr": None},
                {"alignment_at_limit": False, "alignment_corr": float("nan")}, {}):
        assert not board.alignment_ok(bad), bad
        try:
            board.aligned(bad, "dstar_t03_enc_vs_dvsi")
        except RuntimeError as error:
            assert "dstar_t03_enc_vs_dvsi" in str(error), error
        else:
            raise AssertionError(f"an unusable alignment must fail: {bad}")
    # Our encoding, delayed so that DVSI's later one can stay the reference inside the lag search.
    with tempfile.TemporaryDirectory() as work:
        delayed = Path(work) / "delayed.raw"
        board.write_delayed(delayed, [1, -2, 3], 2)
        assert delayed.read_bytes() == struct.pack("<5h", 0, 0, 1, -2, 3)
    # An encoding that cannot be aligned with the input is excluded from the means and pairs, counted, and
    # fails any gate on its partition.
    failed = vector("v1", "validation", 7.0, 0.88)
    failed["encoder"] = {"alignment_failures": ["v1_enc"]}
    flat = board.flatten(failed)
    assert flat["encoder.alignment_failures"] == 1 and "encoder.lsd_db" not in flat, flat
    clean = vector("v2", "validation", 7.0, 0.88)
    clean["encoder"] = dict(result["encoder"], alignment_failures=[])
    assert board.flatten(clean)["encoder.alignment_failures"] == 0
    gates = {"partition": "validation", "checks": [{"mode": "dstar", "metric": "ours.lsd_db", "rule": "no_increase"}]}
    # Both runs' encodings decoded by one decoder, as compare() requires for encoder deltas.
    candidate, baseline = ({"decoder_identity": "d", "vectors": [failed, clean]},
                           {"decoder_identity": "d", "vectors": [failed, clean]})
    passed, lines = board.check_gates(gates, board.aggregate(candidate["vectors"]), board.compare(candidate, baseline))
    assert not passed and any("alignment" in line for line in lines), lines
    # A metric the candidate has but the baseline lacks (or the reverse) is not a pair to skip: the
    # check fails rather than pass on the remaining vectors.
    encoder_gate = {"partition": "validation", "checks": [
        {"mode": "dstar", "metric": "encoder.lsd_db", "rule": "no_increase"}]}
    worse = dict(clean, encoder=dict(clean["encoder"], speech={"lsd_db": 100.0}))
    candidate, baseline = ({"decoder_identity": "d", "vectors": [worse, result]},
                           {"decoder_identity": "d", "vectors": [dict(failed, name="v2"), result]})
    passed, lines = board.check_gates(encoder_gate, board.aggregate(candidate["vectors"]),
                                      board.compare(candidate, baseline))
    assert not passed and any("unpaired" in line for line in lines), lines
    # So is a baseline vector the candidate run left out.
    candidate, baseline = ({"decoder_identity": "d", "vectors": [result]},
                           {"decoder_identity": "d", "vectors": [clean, result]})
    passed, lines = board.check_gates(encoder_gate, board.aggregate(candidate["vectors"]),
                                      board.compare(candidate, baseline))
    assert not passed and any("unpaired" in line and "v2" in line for line in lines), lines
    tone, voice = {"flags": 0x0010}, {"flags": 0x0002}
    assert board.tone_frame_rate([tone, voice, voice, voice]) == 0.25
    silence = record(150.0, flags=0x0100)
    assert board.tone_frame_rate([silence, tone, voice, voice]) == 0.5
    assert board.tone_frame_rate([]) is None


def check_configuration(board):
    seen = set()
    for names in board.PARTITIONS.values():
        assert not seen & set(names), "partitions must not share utterances"
        seen |= set(names)
    gates = json.loads((Path(board.__file__).with_name("dvsi_gates.json")).read_text())
    rules = {"max", "min", "decrease", "increase", "no_increase", "no_decrease", "abs_delta_max",
             "all_no_increase", "all_no_decrease"}
    prefixes = ("ours.", "dvsi.", "vs_dvsi.", "pitch.", "encoder.", "tone.", "speech.", "channel.")
    for name, gate_set in gates.items():
        assert gate_set["partition"] in board.PARTITIONS, name
        for check in gate_set["checks"]:
            assert check["rule"] in rules, (name, check)
            assert check["metric"].startswith(prefixes), (name, check)
            modes = check["mode"] if isinstance(check["mode"], list) else [check["mode"]]
            assert set(modes) <= set(board.MODES), (name, check)


def tone(freqs, levels_db, seconds=1.0, seed=3):
    rng = random.Random(seed)
    amplitudes = [32768.0 * 10 ** (level / 20.0) * math.sqrt(2.0) for level in levels_db]
    phases = [rng.uniform(0, 2 * math.pi) for _ in freqs]
    return [sum(a * math.sin(2 * math.pi * f * n / 8000.0 + p) for f, a, p in zip(freqs, amplitudes, phases))
            for n in range(int(8000 * seconds))]


def check_tone_detector(board):
    # A DTMF "1" with 2 dB twist: both components, their balance and the total level.
    found = board.detect_tone(tone((697.0, 1209.0), (-12.0, -10.0))[1000:1320])
    assert found is not None
    assert abs(found["f1"] / 697.0 - 1) < 2e-3 and abs(found["f2"] / 1209.0 - 1) < 2e-3, found
    assert abs(found["balance_db"] - 2.0) < 0.2, found
    total = 10 * math.log10(10 ** (-1.2) + 10 ** (-1.0))
    assert abs(found["level_db"] - total) < 0.2, found
    single = board.detect_tone(tone((1000.0,), (-20.0,))[1000:1320])
    assert single is not None and single["f2"] is None and abs(single["f1"] / 1000.0 - 1) < 2e-3, single
    rng = random.Random(9)
    assert board.detect_tone([rng.gauss(0, 3000.0) for _ in range(320)]) is None
    assert board.detect_tone(harmonic_signal(140.0)[1000:1320]) is None
    assert board.detect_tone([0.0] * 320) is None


def check_tone_comparison(board):
    # DVSI plays the tone 1 dB louder and starts 40 ms later; ours is delayed 37 samples overall.
    dvsi = [0.0] * 1600 + tone((770.0, 1336.0), (-13.0, -11.0), seconds=0.8)
    ours = [0.0] * 37 + [0.0] * 1280 + tone((770.0, 1336.0), (-14.0, -12.0), seconds=0.84)
    metrics = board.compare_tones(ours, dvsi, 37)
    assert abs(metrics["level_abs_error_db"] - 1.0) < 0.1, metrics
    assert metrics["freq_rel_error"] < 1e-3, metrics
    assert metrics["balance_abs_error_db"] < 0.1, metrics
    assert 0.9 < metrics["detector_agreement"] < 0.99, metrics  # the 40 ms early start disagrees
    same = board.compare_tones(dvsi, dvsi, 0)
    assert same["detector_agreement"] == 1.0 and same["level_abs_error_db"] < 1e-6, same


def check_tone_aggregation(board):
    def tone_result(name, level, agreement):
        return {"mode": "r33", "name": name, "partition": "tones",
                "tone": {"level_abs_error_db": level, "freq_rel_error": 0.0, "balance_abs_error_db": None,
                         "detector_agreement": agreement, "steady_tone_windows": 10}}
    summary = board.aggregate([tone_result("dtmf", 0.2, 1.0), tone_result("alltone", 3.0, 0.95)])
    means = summary["r33/tones"]["mean"]
    # The worst tone vector decides: errors take the maximum, agreement the minimum.
    assert means["tone.level_abs_error_db"] == 3.0 and means["tone.detector_agreement"] == 0.95, means
    assert means["tone.balance_abs_error_db"] is None, means

    # Speech PCM identity against the baseline, counted per mode under the tones group.
    before = vector("dam", "development", 7.0, 0.88)
    before["ours"]["pcm_fnv1a"] = "0x00000001"
    after = json.loads(json.dumps(before))
    after["ours"]["pcm_fnv1a"] = "0x00000002"
    counts = board.pcm_changes({"vectors": [after]}, {"vectors": [before]})
    assert counts == {"dstar": 1}, counts
    assert board.pcm_changes({"vectors": [before]}, {"vectors": [before]}) == {"dstar": 0}


def check_review_regressions(board):
    # Pooled pitch medians come from the frames, not from averaging per-vector medians.
    def pitch_result(name, pairs):
        empty = board.summarize_ratios([])
        return {"mode": "dstar", "name": name, "partition": "validation", "ours": {}, "dvsi": {}, "vs_dvsi": {},
                "pitch": {"ours": empty, "dvsi": board.summarize_ratios(pairs), "input": empty},
                "pitch_pairs": {"ours": [], "dvsi": pairs, "input": []}}
    pooled = board.aggregate([pitch_result("a", [(125.0, 1.02)] * 40), pitch_result("b", [(125.0, 0.96)] * 20)])
    entry = pooled["dstar/validation"]
    assert abs(entry["pitch_dvsi_by_f0"]["100_150"]["median"] - 1.02) < 1e-9, entry["pitch_dvsi_by_f0"]
    assert abs(entry["mean"]["pitch.dvsi_max_band_deviation"] - 0.02) < 1e-9, entry["mean"]

    # A wrong tone for 40% of the time is caught, though level and balance match.
    dvsi = tone((697.0, 1209.0), (-13.0, -13.0), seconds=1.0)
    ours = tone((697.0, 1209.0), (-13.0, -13.0), seconds=0.6) + tone((770.0, 1336.0), (-13.0, -13.0), seconds=0.4)
    metrics = board.compare_tones(ours, dvsi, 0)
    assert 0.35 < metrics["wrong_tone_fraction"] < 0.45, metrics
    assert metrics["detector_agreement"] < 0.7, metrics
    # Adjacent single tones are 31.25 Hz apart: 968.75 Hz in place of 937.5 Hz is a wrong tone.
    single = board.compare_tones(tone((937.5,), (-13.0,), seconds=0.6) + tone((968.75,), (-13.0,), seconds=0.4),
                                 tone((937.5,), (-13.0,), seconds=1.0), 0)
    assert 0.35 < single["wrong_tone_fraction"] < 0.45, single
    # DVSI's off-nominal dual tones (two harmonics of one fundamental, up to 3.4% off) still agree.
    near = board.compare_tones(tone((440.0, 620.0), (-13.0, -13.0)), tone((425.0, 640.0), (-13.0, -13.0)), 0)
    assert near["wrong_tone_fraction"] == 0.0 and near["detector_agreement"] == 1.0, near

    # Speech identity needs the baseline's full coverage and real hashes.
    before = vector("dam", "development", 7.0, 0.88)
    before["ours"]["pcm_fnv1a"] = "0x00000001"
    other = vector("clean", "development", 7.0, 0.88)
    other["ours"]["pcm_fnv1a"] = "0x00000002"
    assert board.pcm_changes({"vectors": [before]}, {"vectors": [before, other]}) == {"dstar": 1}
    unhashed = json.loads(json.dumps(before))
    unhashed["ours"]["pcm_fnv1a"] = None
    assert board.pcm_changes({"vectors": [unhashed]}, {"vectors": [unhashed]}) == {"dstar": 1}

    # The low-band error against DVSI is the size of the 0-250 Hz band difference.
    low = vector("dam", "development", 7.0, 0.88)
    low["vs_dvsi"] = {"band_delta_db_0_250": -1.5}
    assert board.flatten(low)["vs_dvsi.lf_abs_error_db"] == 1.5

    # Quartiles without statistics.quantiles (Python 3.7).
    assert board.quartiles([1.0, 2.0, 3.0, 4.0, 5.0]) == (2.0, 4.0)


def tone_bits(tone_id, ad):
    """A TIA-102.BABA-1 Table 10 tone frame's 49 parameter bits."""
    bits = "111111" + format(ad >> 1, "06b") + format(tone_id, "08b") * 4 + str(ad & 1) + "0000"
    assert len(bits) == 49
    return bits


def dstar_tone_bits(index, volume):
    """A D-STAR tone frame's 49 parameter bits: b0 126, the index's top three bits by the selector in
    bits 6..8, its low five at 9, 42, 43, 10 and 11, and the volume at 12..16, 44, 45 and 17."""
    selector = {4: 0, 0: 1, 1: 2, 2: 3, 3: 4, 7: 5, 6: 6, 5: 7}[index >> 5]
    bits = ["1"] * 6 + ["0"] * 43
    for position, value in zip((6, 7, 8), format(selector, "03b")):
        bits[position] = value
    for position, value in zip((9, 42, 43, 10, 11), format(index & 31, "05b")):
        bits[position] = value
    for position, value in zip((12, 13, 14, 15, 16, 44, 45, 17), format(volume, "08b")):
        bits[position] = value
    return "".join(bits)


def check_dstar_encoder_tones(board):
    for index, volume in ((133, 203), (32, 1), (147, 255), (5, 128), (128, 77)):
        assert board.dstar_tone_fields(dstar_tone_bits(index, volume)) == (index, volume), (index, volume)
    assert board.dstar_tone_fields("0" * 49) is None and board.dstar_tone_fields("1" * 88) is None
    voice = {"bits": "0" * 49}

    def frames(spec):
        return [voice if entry is None else {"bits": dstar_tone_bits(*entry)} for entry in spec]

    dvsi = frames([None, (133, 200), (133, 200), (133, 201), None, None])
    ours = frames([None, (133, 201), (133, 200), (133, 201), None, (32, 180)])
    metrics = board.encoder_tone_agreement(ours, dvsi, mode="dstar")
    assert metrics["encoder_agreement"] == 1.0, metrics
    assert abs(metrics["encoder_extra_rate"] - 1 / 6) < 1e-9, metrics
    assert abs(metrics["encoder_level_error_volume"] - 1 / 3) < 1e-9, metrics
    assert "encoder_level_error_ad" not in metrics, metrics


def check_encoder_tones(board):
    assert board.tone_fields(tone_bits(129, 101)) == (129, 101)
    assert board.tone_fields("0" * 49) is None and board.tone_fields("1" * 88) is None
    assert set(board.ENCODER_COMPARE_DELAY) == set(board.MODES)
    voice = {"bits": "0" * 49}

    def frames(spec):
        return [voice if entry is None else {"bits": tone_bits(*entry)} for entry in spec]

    dvsi = frames([None, (129, 100), (129, 100), (129, 100), (129, 101), None, None, None, None, None])
    # One frame later: three of DVSI's four tone frames agree (one sent as another tone), one level a
    # step off, and one tone frame where DVSI sends voice.
    ours = frames([None, None, (129, 100), (129, 101), (130, 100), (129, 101), None, (129, 99), None, None])
    metrics = board.encoder_tone_agreement(ours, dvsi)
    assert metrics["encoder_agreement"] == 0.75, metrics
    assert abs(metrics["encoder_extra_rate"] - 1 / 9) < 1e-9, metrics
    assert abs(metrics["encoder_level_error_ad"] - 1 / 3) < 1e-9, metrics
    assert board.encoder_tone_agreement(ours, frames([None] * 10)) == {"encoder_extra_rate": 0.5}
    # A tone from the vector's start: shifts 0 and -1 agree equally, and -1, where DVSI's encoder sits,
    # leaves no tone frame where DVSI sends voice.
    dvsi = frames([None] + [(129, 100)] * 9)
    ours = frames([(129, 100)] * 10)
    metrics = board.encoder_tone_agreement(ours, dvsi)
    assert metrics["encoder_agreement"] == 1.0 and metrics["encoder_extra_rate"] == 0.0, metrics
    assert metrics["shift"] == -1, metrics
    # A shift fixed from the baseline is used as given: an encoder that sends DVSI's tone frames three
    # frames late no longer agrees, though the free search would line it up again.
    dvsi = frames([None] * 3 + [(129, 100)] * 4 + [None] * 3)
    late = frames([None] * 5 + [(129, 100)] * 4 + [None])
    assert board.encoder_tone_agreement(late, dvsi)["encoder_agreement"] == 1.0
    fixed = board.encoder_tone_agreement(late, dvsi, shift=-1)
    assert fixed["shift"] == -1 and fixed["encoder_agreement"] == 0.25, fixed
    # The encoders a build offers, from its usage line.
    assert board.usage_codecs("Usage: x --codec ambe2400|ambe2450|imbe7200 --in A") == {"ambe2400", "ambe2450", "imbe7200"}
    assert board.usage_codecs("Usage: x --codec ambe2400 --in A") == {"ambe2400"}
    assert board.usage_codecs("no usage") == frozenset()
    # Encoder-row identity: counted only where the baseline recorded hashes.
    def encoded(name, digest, mode="dstar"):
        return {"mode": mode, "name": name, "encoder": {"rows_sha256": digest}}
    base = {"vectors": [encoded("t03", "a"), encoded("t04", "b"), {"mode": "p25", "name": "t03", "encoder": {}}]}
    same = {"vectors": [encoded("t03", "a"), encoded("t04", "b")]}
    changed = {"vectors": [encoded("t03", "a"), encoded("t04", "c")]}
    missing = {"vectors": [encoded("t03", "a")]}
    assert board.encoder_row_changes(same, base) == {"dstar": 0}
    assert board.encoder_row_changes(changed, base) == {"dstar": 1}
    assert board.encoder_row_changes(missing, base) == {"dstar": 1}
    summary = board.aggregate([
        {"mode": "r33", "name": "dtmf", "partition": "tones", "tone": {"encoder_agreement": 0.98}},
        {"mode": "r33", "name": "alltone", "partition": "tones", "tone": {"encoder_agreement": 0.9}},
    ])
    assert summary["r33/tones"]["mean"]["tone.encoder_agreement"] == 0.9, summary
    assert abs(summary["r33/tones"]["mean"]["tone.encoder_agreement_mean"] - 0.94) < 1e-9, summary


def check_recordings(board):
    # Vectors of one source recording (its level, filter and noise variants) form one group, and no group
    # spans two speech partitions, so validation and final stay independent of development.
    assert board.utterance_group("lvl") == board.utterance_group("dam") == board.utterance_group("damcar")
    assert board.utterance_group("ucarf15a") == board.utterance_group("fambm22a") == board.utterance_group("ucarm20a")
    assert board.utterance_group("tambf32b") == board.utterance_group("tambf22b")
    assert board.utterance_group("p01mirs") == board.utterance_group("tia11")
    assert board.utterance_group("t03") == "t03" and board.utterance_group("irstia") == "irstia"
    assert board.utterance_group("dam_e10_hd") == board.utterance_group("dam")
    partitions = {}
    for partition in board.SPEECH_PARTITIONS:
        for name in board.PARTITIONS[partition]:
            group = board.utterance_group(name)
            assert partitions.setdefault(group, partition) == partition, (name, group, partition)


def check_group_bootstrap(board):
    # Four vectors of one recording improve and three others get slightly worse. Per vector the mean is
    # clearly negative, but resampled by recording (four groups) the interval reaches above zero.
    deltas = {"dam": -1.0, "dam10": -1.0, "dam20": -1.0, "lvl": -1.0, "t01": 0.1, "t02": 0.1, "clean": 0.1}
    baseline = {"vectors": [vector(n, "development", 7.0, 0.88) for n in deltas]}
    candidate = {"vectors": [vector(n, "development", 7.0 + d, 0.88) for n, d in deltas.items()]}
    entry = board.compare(candidate, baseline)["dstar/development"]["ours.lsd_db"]
    assert entry["n"] == 7 and entry["groups"] == 4, entry
    assert abs(entry["mean"] + 3.7 / 7) < 1e-9, entry
    low, high = entry["ci95"]
    assert low < 0 < high, entry
    gates = {"partition": "development", "checks": [{"mode": "dstar", "metric": "ours.lsd_db", "rule": "decrease"}]}
    passed, lines = board.check_gates(gates, board.aggregate(candidate["vectors"]), board.compare(candidate, baseline))
    assert not passed, lines


def check_every_vector_rules(board):
    # The mean improves while one vector gets worse: the per-vector rules catch it, within a tolerance.
    names = [f"v{i}" for i in range(4)]
    baseline = {"vectors": [vector(n, "validation", 7.0, 0.88) for n in names]}
    candidate = {"vectors": [vector(n, "validation", 7.0 + (0.2 if i == 0 else -0.3), 0.88 + (-0.01 if i == 1 else 0.01))
                             for i, n in enumerate(names)]}
    comparison = board.compare(candidate, baseline)
    entry = comparison["dstar/validation"]["ours.lsd_db"]
    assert abs(entry["max"] - 0.2) < 1e-9 and abs(entry["min"] + 0.3) < 1e-9, entry
    rising = {"partition": "validation", "checks": [{"mode": "dstar", "metric": "ours.lsd_db", "rule": "all_no_increase"}]}
    assert not board.check_gates(rising, {}, comparison)[0]
    rising["checks"][0]["value"] = 0.25
    assert board.check_gates(rising, {}, comparison)[0]
    falling = {"partition": "validation", "checks": [{"mode": "dstar", "metric": "ours.env_corr", "rule": "all_no_decrease"}]}
    assert not board.check_gates(falling, {}, comparison)[0]
    falling["checks"][0]["value"] = 0.02
    assert board.check_gates(falling, {}, comparison)[0]


def check_frozen_lag_pairing(board):
    # With lags frozen from the baseline, a vector measured at other lags (one whose baseline recorded
    # none, or a lag the run had to search) is not a pair: all its metrics go unpaired.
    def result(name, lags):
        entry = vector(name, "validation", 7.0, 0.88)
        entry["lags"] = lags
        return entry

    baseline = {"vectors": [result("t03", {"ours": 352, "vs_dvsi": 80}), result("t04", {"ours": 344, "vs_dvsi": -8})]}
    same = {"lags_from": "baseline.json", "vectors": [result("t03", {"ours": 352, "vs_dvsi": 80}),
                                                      result("t04", {"ours": 344, "vs_dvsi": -8})]}
    assert board.compare(same, baseline)["dstar/validation"]["ours.lsd_db"]["unpaired"] == []
    moved = {"lags_from": "baseline.json", "vectors": [result("t03", {"ours": 352, "vs_dvsi": 80}),
                                                       result("t04", {"ours": 344, "vs_dvsi": -8, "ref": 16})]}
    comparison = board.compare(moved, baseline)
    entry = comparison["dstar/validation"]["ours.lsd_db"]
    assert entry["unpaired"] == ["t04"] and entry["n"] == 1, entry
    gates = {"partition": "validation", "checks": [{"mode": "dstar", "metric": "ours.lsd_db", "rule": "no_increase"}]}
    passed, lines = board.check_gates(gates, {}, comparison)
    assert not passed and any("unpaired" in line for line in lines), lines
    # Without frozen lags (an old baseline) the lags are free to differ.
    free = dict(moved)
    del free["lags_from"]
    assert board.compare(free, baseline)["dstar/validation"]["ours.lsd_db"]["unpaired"] == []


def check_joint_lags(board):
    # Where two runs' alignment searches disagree, both are rescored at the lag of the run that aligned
    # better (a decoder fixed by the change can align where the old one could not); equal correlations
    # keep the baseline's lag, and the encoder-tone frame shift stays the baseline's.
    def run(entries):
        return {"vectors": [{"mode": "dstar", "name": name, "lags": lags, "alignment": corr}
                            for name, lags, corr in entries]}

    baseline = run([("t03", {"ours": 352, "vs_dvsi": 56}, {"ours": 0.97, "vs_dvsi": 0.985}),
                    ("t04", {"ours": 344, "vs_dvsi": 56}, {"ours": 0.98, "vs_dvsi": 0.99}),
                    ("dtmf", {"vs_dvsi": 40, "encoder_tone_shift": -1}, {"vs_dvsi": 0.957}),
                    ("t01", {"ours": 300}, {"ours": 0.9}),
                    ("t02", {"ours": 300, "enc": 152}, {"ours": 0.9, "enc": None})])
    candidate = run([("t03", {"ours": 352, "vs_dvsi": -8}, {"ours": 0.96, "vs_dvsi": 0.995}),
                     ("t04", {"ours": 336, "vs_dvsi": 64}, {"ours": 0.97, "vs_dvsi": 0.99}),
                     ("dtmf", {"vs_dvsi": -64, "encoder_tone_shift": 0}, {"vs_dvsi": 0.963}),
                     ("t02", {"ours": 300, "enc": 144}, {"ours": 0.9, "enc": 0.98})])
    joint = board.joint_lags(baseline, candidate)
    assert joint[("dstar", "t03")] == {"ours": 352, "vs_dvsi": -8}, joint
    assert joint[("dstar", "t04")] == {"ours": 344, "vs_dvsi": 56}, joint
    assert joint[("dstar", "dtmf")] == {"vs_dvsi": -64, "encoder_tone_shift": -1}, joint
    assert joint[("dstar", "t02")] == {"ours": 300, "enc": 144}, joint  # an unaligned run loses
    assert ("dstar", "t01") not in joint  # scored by one run only
    with tempfile.TemporaryDirectory() as directory:
        first, second = Path(directory) / "a.json", Path(directory) / "b.json"
        first.write_text(json.dumps(baseline))
        second.write_text(json.dumps(candidate))
        assert board.load_lag_source([first]) == {(r["mode"], r["name"]): r["lags"] for r in baseline["vectors"]}
        assert board.load_lag_source([first, second]) == joint


def check_decoder_identity(board):
    with tempfile.TemporaryDirectory() as directory:
        first, second = Path(directory) / "a", Path(directory) / "b"
        for tools in (first, second):
            tools.mkdir()
            (tools / "mbe_quality_eval").write_bytes(b"evaluator")
            (tools / "libmbe-neo.so.2").write_bytes(b"library 1")
        assert board.decoder_identity(first) == board.decoder_identity(second)
        (second / "libmbe-neo.so.2").write_bytes(b"library 2")
        assert board.decoder_identity(first) != board.decoder_identity(second)
        (second / "libmbe-neo.so.2").write_bytes(b"library 1")
        (second / "mbe_quality_eval").write_bytes(b"other evaluator")
        assert board.decoder_identity(first) != board.decoder_identity(second)

    # Encoder deltas are measured only when one decoder decoded both runs' encodings.
    def run_with(identity):
        result = vector("t03", "validation", 7.0, 0.88)
        result["encoder"] = {"speech": {"lsd_db": 5.0}, "alignment_failures": []}
        return {"decoder_identity": identity, "vectors": [result]}

    assert abs(board.compare(run_with("x"), run_with("x"))["dstar/validation"]["encoder.lsd_db"]["mean"]) < 1e-12
    # An unrecorded decoder (a scoreboard from before identities were recorded) is not the same decoder.
    for candidate, baseline in ((run_with("x"), run_with("y")), (run_with("x"), run_with(None)),
                                (run_with(None), run_with(None))):
        try:
            board.compare(candidate, baseline)
        except ValueError as error:
            assert "decoder" in str(error), error
        else:
            raise AssertionError("encoder deltas across decoders must be refused")
    plain = {"decoder_identity": "y", "vectors": [vector("t03", "validation", 7.0, 0.88)]}
    board.compare(plain, dict(plain, decoder_identity="x"))  # no encoder metrics: nothing to refuse


class FakeTools:
    """Stands in for the mbe_quality_* tools: records every command and writes the outputs the scoreboard
    reads. Each comparison's automatic lag differs, and a forced lag is reported as the lag used."""

    def __init__(self):
        self.commands = []
        self.bits = "0" * 49  # every decoded frame's parameter bits

    def __call__(self, command):
        command = [str(part) for part in command]
        self.commands.append(command)
        options, index = {}, 1
        while index < len(command):
            if command[index] == "--soft":
                options["--soft"], index = True, index + 1
            else:
                options[command[index]], index = command[index + 1], index + 2
        tool = Path(command[0]).name
        if tool in ("mbe_quality_reframe", "mbe_quality_encode"):
            Path(options["--out"]).write_text("rows\n")
        elif tool == "mbe_quality_eval":
            if "--out" in options:
                samples = [int(1000 * math.sin(n / 5.0)) for n in range(1600)]
                Path(options["--out"]).write_bytes(struct.pack(f"<{len(samples)}h", *samples))
            if "--params" in options:
                records = [{"f": f, "flags": 0, "errors": 0, "bits": self.bits, "decoded": {"status": 0}}
                           for f in range(10)]
                Path(options["--params"]).write_text("".join(json.dumps(r) + "\n" for r in records))
            if "--json" in options:
                automatic = 8 * len(self.commands)
                metrics = {"alignment_corr": 0.9 + 0.001 * len(self.commands), "alignment_at_limit": False,
                           "auto_lag_samples": automatic,
                           "lag_samples": int(options.get("--lag-samples", automatic)), "level_offset_db": 0.0,
                           "pcm_fnv1a": "0x00000001"}
                Path(options["--json"]).write_text(json.dumps(metrics))
        else:
            raise AssertionError(command)
        return subprocess.CompletedProcess(command, 0, "", "")


def check_scoring_commands(board):
    fake = FakeTools()
    board.run, board.encoder_codecs = fake, lambda *_: frozenset(board.MODES.values())
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        vectors, work, tools, decoder = root / "vectors", root / "work", root / "candidate", root / "decoder"
        for path in (vectors / "dstar", vectors / "r33", work):
            path.mkdir(parents=True)
        silence = struct.pack("<1600h", *([0] * 1600))
        for path in (vectors / "t03.pcm", vectors / "dstar" / "t03.pcm", vectors / "dstar" / "t03.bit"):
            path.write_bytes(silence)

        # With a baseline's lags every comparison is measured at the baseline's lag, and the encoder
        # comparisons decode both encodings (DVSI's bits again, too) with the fixed decoder.
        frozen = {"ours": 352, "dvsi": 256, "vs_dvsi": 80, "ref": 344, "enc": 152, "enc_vs_dvsi": 136}
        result = board.score_vector(tools, vectors, "dstar", "t03", work, True, decode_tools=decoder, lags=frozen)
        assert result["lags"] == frozen, result["lags"]
        by_output = {Path(c[c.index("--json") + 1]).name: c for c in fake.commands if "--json" in c}
        # Each comparison's alignment correlation is recorded beside its lag, for joint_lags().
        assert set(result["alignment"]) == set(frozen), result["alignment"]
        for key in frozen:
            index = fake.commands.index(by_output[f"dstar_t03_{key}.json"]) + 1
            assert abs(result["alignment"][key] - (0.9 + 0.001 * index)) < 1e-12, (key, result["alignment"])
        expected = {"ours": tools, "dvsi": tools, "vs_dvsi": tools, "ref": decoder, "enc": decoder, "enc_vs_dvsi": decoder}
        assert set(by_output) == {f"dstar_t03_{key}.json" for key in expected}, sorted(by_output)
        for key, used in expected.items():
            command = by_output[f"dstar_t03_{key}.json"]
            assert command[command.index("--lag-samples") + 1] == str(frozen[key]), (key, command)
            assert Path(command[0]).parent == used, (key, command)
        encodes = [c for c in fake.commands if Path(c[0]).name == "mbe_quality_encode"]
        assert len(encodes) == 1 and Path(encodes[0][0]).parent == tools, encodes

        # Without them each comparison finds its own lag, which the result records for a later run.
        fake.commands.clear()
        result = board.score_vector(tools, vectors, "dstar", "t03", work, True)
        assert not any("--lag-samples" in command for command in fake.commands), fake.commands
        assert set(result["lags"]) == set(frozen) and result["lags"]["ours"] != frozen["ours"], result["lags"]

        # A soft channel-error vector decodes as soft input, and its alignment with DVSI's decoder comes
        # from the clean vector, at the baseline's lag when there is one.
        for name in ("dam.pcm", "dam.bit", "dam_e10_sd.pcm", "dam_e10_sd.bit"):
            (vectors / "r33" / name).write_bytes(silence)
        fake.commands.clear()
        result = board.score_channel_vector(tools, vectors, "r33", "dam_e10_sd", work, lags={"clean_vs_dvsi": -8})
        assert result["partition"] == "channel" and result["kind"] == "sd", result
        assert result["lags"] == {"clean_vs_dvsi": -8}, result["lags"]
        assert set(result["alignment"]) == {"clean_vs_dvsi"} and 0.9 < result["alignment"]["clean_vs_dvsi"] < 1.0
        decodes = [c for c in fake.commands if "--frames" in c and "dam_e10_sd" in Path(c[c.index("--frames") + 1]).name]
        assert len(decodes) == 1 and "--soft" in decodes[0], fake.commands
        aligned = [c for c in fake.commands if "--json" in c]
        assert len(aligned) == 1 and aligned[0][aligned[0].index("--lag-samples") + 1] == "-8", aligned

        # A tone vector records its lag and the encoder-tone frame shift it chose, and a later run
        # measures at the baseline's.
        for path in (vectors / "dtmf.pcm", vectors / "r33" / "dtmf.pcm", vectors / "r33" / "dtmf.bit"):
            path.write_bytes(silence)
        fake.bits = tone_bits(129, 100)  # all ten frames, DVSI's and ours, carry the same tone
        result = board.score_tone_vector(tools, vectors, "r33", "dtmf", work, True)
        assert result["lags"]["encoder_tone_shift"] == 0 and "vs_dvsi" in result["lags"], result["lags"]
        assert result["tone"]["encoder_agreement"] == 1.0 and "shift" not in result["tone"], result["tone"]
        fake.commands.clear()
        result = board.score_tone_vector(tools, vectors, "r33", "dtmf", work, True, decode_tools=decoder,
                                         lags={"vs_dvsi": 24, "encoder_tone_shift": 2})
        assert result["lags"] == {"vs_dvsi": 24, "encoder_tone_shift": 2}, result["lags"]
        assert set(result["alignment"]) == {"vs_dvsi"}, result["alignment"]
        assert abs(result["tone"]["encoder_agreement"] - 0.8) < 1e-12, result["tone"]  # frames 8 and 9 fall off
        encoded = [c for c in fake.commands if "--frames" in c and "_enc.rows" in c[c.index("--frames") + 1]]
        assert len(encoded) == 1 and Path(encoded[0][0]).parent == decoder, encoded


def check_channel(board):
    rng = random.Random(5)
    clean = [rng.gauss(0, 3000.0) for _ in range(1280)] + [rng.gauss(0, 50.0) for _ in range(320)]
    errored = list(clean)
    for f in (3, 4, 5):  # DVSI mutes frames 3 to 5
        errored[160 * f:160 * (f + 1)] = [0.01 * x for x in clean[160 * f:160 * (f + 1)]]
    loud, dvsi_muted = board.channel_frames(clean, errored)
    assert loud == [True] * 8 + [False] * 2, loud
    assert dvsi_muted == [False] * 3 + [True] * 3 + [False] * 4, dvsi_muted

    # Our mute flags against DVSI's muting, frame for frame on loud frames.
    ours = [False] * 4 + [True] * 4 + [False] * 2
    metrics = board.mute_metrics(loud, dvsi_muted, ours)
    assert abs(metrics["mute_agreement"] - 5 / 8) < 1e-12 and abs(metrics["false_mute_rate"] - 2 / 8) < 1e-12, metrics
    assert abs(metrics["missed_mute_rate"] - 1 / 8) < 1e-12 and metrics["loud_frames"] == 8, metrics
    assert metrics["first_mute_frame"] == 4 and metrics["dvsi_first_mute_frame"] == 3, metrics
    assert board.mute_metrics([False] * 4, [False] * 4, [True] * 4)["mute_agreement"] is None
    # Muting is only seen on loud frames, so our first mute counts from there too.
    assert board.mute_metrics([False, True, True], [False, False, True], [True, False, True])["first_mute_frame"] == 2

    # Per-frame spectral distance: none for the same signal at its lag and level, some where it differs.
    reference = [rng.gauss(0, 2000.0) for _ in range(1600)]
    delayed = [0.0] * 40 + [0.5 * x for x in reference]
    lsd = board.frame_lsd(delayed, reference, 40, 2.0)
    # No 256-sample spectrum is centred in the last frame.
    assert len(lsd) == 10 and lsd[9] is None and max(lsd[:9]) < 1e-6, lsd
    changed = list(delayed)
    changed[40 + 960:40 + 1120] = [rng.gauss(0, 1000.0) for _ in range(160)]
    lsd = board.frame_lsd(changed, reference, 40, 2.0)
    assert lsd[6] > 3.0 and max(lsd[:4]) < 1e-6, lsd

    # Reference frame g holds our samples 160 g + lag onward, so every frame of ours under it must be unmuted.
    assert board.covered_unmuted([0, 0, 0, 1, 0], 80, 5) == [True, True, False, False, False]
    assert board.covered_unmuted([0, 0, 0, 1, 0], -8, 5) == [False, True, True, False, False]

    # LSD on the frames that count in both runs: loud, measured, and with no muted frame (DVSI's, among loud
    # frames, or either run's) at or next to them, since a frame's spectra reach into its neighbours.
    loud = [True] * 18 + [False] * 2
    dvsi_muted = [f in (8, 9) for f in range(20)]

    def run(muted_frame, lsd):
        return {"mode": "r33", "name": "dam_e10_hd", "partition": "channel", "kind": "hd",
                "channel": {"mute_agreement": 0.5}, "lags": {"clean_vs_dvsi": 0},
                "frames": {"loud": loud, "dvsi_muted": dvsi_muted, "muted": [f == muted_frame for f in range(20)],
                           "lsd": lsd}}

    # Frames 1, 3, 7, 10, 13 and 15 neighbour a muted frame: the candidate's 9 dB there must not count.
    near = (1, 3, 7, 8, 9, 10, 13, 14, 15, 18, 19)
    candidate = run(14, [9.0 if f in near else 2.0 for f in range(20)])
    baseline = run(2, [9.0 if f in (2, 8, 9, 14, 18, 19) else 3.0 for f in range(20)])
    comparison = board.compare({"vectors": [candidate]}, {"vectors": [baseline]})
    entry = comparison["r33/channel"]["channel.hd.lsd_common_db"]
    # Common: frames 0, 4, 5, 6, 11, 12, 16 and 17 of the 16 loud frames DVSI leaves unmuted.
    assert abs(entry["mean"] + 1.0) < 1e-9 and abs(entry["coverage"] - 8 / 16) < 1e-12, entry
    own = board.flatten(candidate)
    # The candidate alone: frames 0-6, 11, 12, 16 and 17, two of them (1 and 3) at 9 dB.
    assert abs(own["channel.hd.lsd_db"] - 36 / 11) < 1e-9 and abs(own["channel.hd.lsd_coverage"] - 11 / 16) < 1e-12, own
    assert abs(own["channel.hd.mute_agreement"] - 0.5) < 1e-12, own
    # Frames a run did not measure (a short or missing output) are lost coverage, not a smaller denominator.
    short = run(14, [None if f in (11, 12) else v for f, v in enumerate(candidate["frames"]["lsd"])])
    own = board.flatten(short)
    assert abs(own["channel.hd.lsd_db"] - 32 / 9) < 1e-9 and abs(own["channel.hd.lsd_coverage"] - 9 / 16) < 1e-12, own
    short["frames"]["lsd"] = short["frames"]["lsd"][:4]  # an output that ends after frame 3
    assert abs(board.flatten(short)["channel.hd.lsd_coverage"] - 4 / 16) < 1e-12
    summary = board.aggregate([candidate])
    assert abs(summary["r33/channel"]["mean"]["channel.hd.mute_agreement"] - 0.5) < 1e-12, summary
    # A delta on too small a share of the frames DVSI leaves unmuted fails a gate that sets a floor.
    gate = {"partition": "channel", "checks": [
        {"mode": "r33", "metric": "channel.hd.lsd_common_db", "rule": "no_increase", "min_coverage": 0.7}]}
    assert not board.check_gates(gate, {}, comparison)[0]
    gate["checks"][0]["min_coverage"] = 0.5
    assert board.check_gates(gate, {}, comparison)[0]


def check_unscored_vectors(board):
    # A gate judges its whole partition: every vector of it must be scored, for every mode it gates, in
    # both runs, or the gate set fails.
    gate_set = {"partition": "channel", "checks": [
        {"mode": ["r33", "dstar"], "metric": "channel.hd.mute_agreement", "rule": "all_no_decrease"}]}
    every = {"vectors": [{"mode": mode, "name": name} for mode in ("r33", "dstar")
                         for name in board.PARTITIONS["channel"]]}
    assert board.unscored_vectors(gate_set, every, every) == []
    some = {"vectors": [v for v in every["vectors"] if v["name"] in ("dam_e1_hd", "dam_e10_hd")]}
    missing = board.unscored_vectors(gate_set, some, every)
    assert len(missing) == 12 and ("r33", "dam_e5_sd") in missing and ("dstar", "dam_e2_hd") in missing, missing
    assert board.unscored_vectors(gate_set, every, some) == missing
    summary = {f"{mode}/channel": {"mean": {"channel.hd.mute_rate": 0.1}} for mode in ("r33", "dstar")}
    rate = {"partition": "channel", "checks": [
        {"mode": ["r33", "dstar"], "metric": "channel.hd.mute_rate", "rule": "max", "value": 0.5}]}
    fixed = dict(every, lags_from="baseline.json")
    passed, lines = board.evaluate_gate_sets(["rate"], {"rate": rate}, summary, {}, (fixed, every))
    assert passed and lines[-1] == "gate set rate: PASS", lines
    passed, lines = board.evaluate_gate_sets(["rate"], {"rate": rate}, summary, {}, (fixed, some))
    assert not passed and any("not scored" in line for line in lines), lines
    # A run not scored at fixed lags (its baseline recorded none) cannot pass a gate: its comparisons
    # and the baseline's may sit at different lags.
    passed, lines = board.evaluate_gate_sets(["rate"], {"rate": rate}, summary, {}, (every, every))
    assert not passed and any("lags" in line for line in lines), lines
    try:
        board.evaluate_gate_sets(["other"], {"rate": rate}, summary, {}, (fixed, every))
    except ValueError as error:
        assert "other" in str(error), error
    else:
        raise AssertionError("an undeclared gate set must be refused")


def main():
    board = load_scoreboard()
    check_configuration(board)
    check_recordings(board)
    check_group_bootstrap(board)
    check_every_vector_rules(board)
    check_frozen_lag_pairing(board)
    check_joint_lags(board)
    check_decoder_identity(board)
    check_unscored_vectors(board)
    check_voicing(board)
    check_gates(board)
    check_encoder_metrics(board)
    check_tone_aggregation(board)
    check_encoder_tones(board)
    check_dstar_encoder_tones(board)
    if importlib.util.find_spec("numpy") is None:
        print("NumPy not available: pitch checks skipped")
        return SKIP
    check_pitch(board)
    check_frame_selection(board)
    check_tone_detector(board)
    check_tone_comparison(board)
    check_review_regressions(board)
    check_channel(board)
    check_scoring_commands(board)
    return 0


if __name__ == "__main__":
    sys.exit(main())
