#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the DVSI scoreboard's measurements and gates on synthetic data (no DVSI material)."""

import importlib.util
import json
import math
from pathlib import Path
import random
import sys

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


def check_configuration(board):
    seen = set()
    for names in board.PARTITIONS.values():
        assert not seen & set(names), "partitions must not share utterances"
        seen |= set(names)
    gates = json.loads((Path(board.__file__).with_name("dvsi_gates.json")).read_text())
    rules = {"max", "min", "decrease", "increase", "no_increase", "no_decrease", "abs_delta_max"}
    prefixes = ("ours.", "dvsi.", "vs_dvsi.", "pitch.", "encoder.")
    for name, gate_set in gates.items():
        assert gate_set["partition"] in board.PARTITIONS, name
        for check in gate_set["checks"]:
            assert check["rule"] in rules, (name, check)
            assert check["metric"].startswith(prefixes), (name, check)
            modes = check["mode"] if isinstance(check["mode"], list) else [check["mode"]]
            assert set(modes) <= set(board.MODES), (name, check)


def main():
    board = load_scoreboard()
    check_configuration(board)
    check_voicing(board)
    check_gates(board)
    if importlib.util.find_spec("numpy") is None:
        print("NumPy not available: pitch checks skipped")
        return SKIP
    check_pitch(board)
    check_frame_selection(board)
    return 0


if __name__ == "__main__":
    sys.exit(main())
