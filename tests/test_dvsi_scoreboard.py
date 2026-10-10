#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the DVSI scoreboard's measurements and gates on synthetic data (no DVSI material)."""

import importlib.util
import json
import math
from pathlib import Path
import random
import struct
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
    candidate, baseline = {"vectors": [failed, clean]}, {"vectors": [failed, clean]}
    passed, lines = board.check_gates(gates, board.aggregate(candidate["vectors"]), board.compare(candidate, baseline))
    assert not passed and any("alignment" in line for line in lines), lines
    # A metric the candidate has but the baseline lacks (or the reverse) is not a pair to skip: the
    # check fails rather than pass on the remaining vectors.
    encoder_gate = {"partition": "validation", "checks": [
        {"mode": "dstar", "metric": "encoder.lsd_db", "rule": "no_increase"}]}
    worse = dict(clean, encoder=dict(clean["encoder"], speech={"lsd_db": 100.0}))
    candidate, baseline = {"vectors": [worse, result]}, {"vectors": [dict(failed, name="v2"), result]}
    passed, lines = board.check_gates(encoder_gate, board.aggregate(candidate["vectors"]),
                                      board.compare(candidate, baseline))
    assert not passed and any("unpaired" in line for line in lines), lines
    # So is a baseline vector the candidate run left out.
    candidate, baseline = {"vectors": [result]}, {"vectors": [clean, result]}
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
    rules = {"max", "min", "decrease", "increase", "no_increase", "no_decrease", "abs_delta_max"}
    prefixes = ("ours.", "dvsi.", "vs_dvsi.", "pitch.", "encoder.", "tone.", "speech.")
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


def main():
    board = load_scoreboard()
    check_configuration(board)
    check_voicing(board)
    check_gates(board)
    check_encoder_metrics(board)
    check_tone_aggregation(board)
    check_encoder_tones(board)
    if importlib.util.find_spec("numpy") is None:
        print("NumPy not available: pitch checks skipped")
        return SKIP
    check_pitch(board)
    check_frame_selection(board)
    check_tone_detector(board)
    check_tone_comparison(board)
    check_review_regressions(board)
    return 0


if __name__ == "__main__":
    sys.exit(main())
