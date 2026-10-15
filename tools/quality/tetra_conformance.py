#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Compare TETRA decoder output with a reference, and make bad-frame test inputs.

    tetra_conformance.py compare OURS.pcm REFERENCE.pcm [--gate] [--metrics-only] [--json OUT]
    tetra_conformance.py bfi INPUT.chd OUTPUT.chd [--seed N] [--rate PERCENT]
                             [--burst START:LENGTH ...] [--first N]

compare reads two raw little-endian 16-bit PCM files of equal, nonzero length (240 samples
per 30 ms frame) and reports how many frames match exactly, the first mismatching sample,
the overall and segmental SNR, the level difference, and the frames below the SNR floor. By
default it exits 1 unless the files are identical. With --gate it exits 1 unless the files
meet the acceptance thresholds for this decoder, which is not bit-exact (docs/tetra-acelp.md):

- segmental SNR of at least 30 dB over the frames whose reference level is above -60 dBFS;
- level difference within 0.1 dB;
- no frame above -50 dBFS below 14 dB SNR, except the first 10 frames;
- every all-zero reference frame within 2 LSB in ours.

Frames with zero reference energy are counted separately and never averaged, so no metric
is NaN. --metrics-only always exits 0 after reporting.

bfi copies a .CHD file (138 little-endian 16-bit words per frame: the bad-frame indicator,
then B1..B137, as in EN 300 395-4) and sets the indicator of the chosen frames to 1: a
seeded random PERCENT of all frames, explicit bursts, and/or the first N frames.

The ETSI conformance files and the reference decoder's output stay outside the repository.
"""

import argparse
import json
import math
import os
import random
import struct
import sys
from pathlib import Path

FRAME_SAMPLES = 240
CHD_WORDS = 138
FULL_SCALE = 32768.0
SEG_FLOOR_DB = -60.0
LOW_FLOOR_DB = -50.0
MIN_SEG_SNR_DB = 30.0
MAX_LEVEL_DB = 0.1
MIN_FRAME_SNR_DB = 14.0
SETTLE_FRAMES = 10
SILENT_TOLERANCE = 2
SNR_CAP_DB = 99.0


class ToolError(Exception):
    pass


def read_pcm(path):
    data = Path(path).read_bytes()
    if not data:
        raise ToolError(f"{path}: empty")
    if len(data) % 2:
        raise ToolError(f"{path}: odd byte count {len(data)}")
    return struct.unpack(f"<{len(data) // 2}h", data)


def snr_db(signal, noise):
    if noise == 0:
        return SNR_CAP_DB
    return min(SNR_CAP_DB, 10.0 * math.log10(signal / noise))


def frame_stats(ours, ref):
    signal = sum(int(r) * int(r) for r in ref)
    noise = sum((int(o) - int(r)) ** 2 for o, r in zip(ours, ref))
    return signal, noise


def compare(ours, ref):
    if len(ours) != len(ref):
        raise ToolError(f"lengths differ: {len(ours)} and {len(ref)} samples")
    frames = (len(ref) + FRAME_SAMPLES - 1) // FRAME_SAMPLES
    exact = 0
    seg = []
    low = []
    silent = 0
    silent_bad = []
    total_signal = 0
    total_noise = 0
    ours_energy = 0
    for f in range(frames):
        o = ours[f * FRAME_SAMPLES:(f + 1) * FRAME_SAMPLES]
        r = ref[f * FRAME_SAMPLES:(f + 1) * FRAME_SAMPLES]
        signal, noise = frame_stats(o, r)
        total_signal += signal
        total_noise += noise
        ours_energy += sum(int(v) * int(v) for v in o)
        exact += noise == 0
        if signal == 0:
            silent += 1
            if any(abs(int(v)) > SILENT_TOLERANCE for v in o):
                silent_bad.append(f)
            continue
        level = 10.0 * math.log10(signal / len(r) / (FULL_SCALE * FULL_SCALE))
        if level <= SEG_FLOOR_DB:
            continue
        snr = snr_db(signal, noise)
        seg.append(snr)
        if snr < MIN_FRAME_SNR_DB and level > LOW_FLOOR_DB and f >= SETTLE_FRAMES:
            low.append({"frame": f, "level_dbfs": round(level, 2), "snr_db": round(snr, 2)})
    first = next((i for i, (o, r) in enumerate(zip(ours, ref)) if o != r), None)
    if total_signal and ours_energy:
        level_db = 10.0 * math.log10(ours_energy / total_signal)
    elif total_signal == ours_energy:
        level_db = 0.0
    else:
        level_db = None  # one side is all zero
    return {
        "frames": frames,
        "exact_frames": exact,
        "bit_exact": first is None,
        "first_mismatch": None if first is None else {"frame": first // FRAME_SAMPLES, "sample": first % FRAME_SAMPLES},
        "snr_db": snr_db(total_signal, total_noise) if total_signal else None,
        "segmental_snr_db": sum(seg) / len(seg) if seg else None,
        "segmental_frames": len(seg),
        "level_db": level_db,
        "low_snr_frames": low,
        "silent_reference_frames": silent,
        "silent_frames_not_silent": silent_bad,
    }


def rounded(metrics):
    """A copy for reporting; the gate is always evaluated on the unrounded values."""
    shown = dict(metrics)
    for key, digits in (("snr_db", 2), ("segmental_snr_db", 2), ("level_db", 3)):
        if shown[key] is not None:
            shown[key] = round(shown[key], digits)
    return shown


def gate_failures(metrics):
    failures = []
    seg = metrics["segmental_snr_db"]
    if seg is not None and seg < MIN_SEG_SNR_DB:
        failures.append(f"segmental SNR {seg:.4f} dB < {MIN_SEG_SNR_DB} dB")
    level = metrics["level_db"]
    if level is None or abs(level) > MAX_LEVEL_DB:
        failures.append(f"level difference {'none' if level is None else f'{level:.4f}'} dB beyond {MAX_LEVEL_DB} dB")
    if metrics["low_snr_frames"]:
        failures.append(f"{len(metrics['low_snr_frames'])} frames below {MIN_FRAME_SNR_DB} dB")
    if metrics["silent_frames_not_silent"]:
        failures.append(f"{len(metrics['silent_frames_not_silent'])} silent reference frames not silent")
    return failures


def same_file(a, b):
    a, b = Path(a), Path(b)
    try:
        if a.resolve() == b.resolve():
            return True
    except (OSError, RuntimeError):  # before Python 3.13 a symlink loop raises RuntimeError
        pass
    try:
        return os.path.samefile(a, b)
    except OSError:
        return False


def read_chd(path):
    data = Path(path).read_bytes()
    if not data:
        raise ToolError(f"{path}: empty")
    if len(data) % (2 * CHD_WORDS):
        raise ToolError(f"{path}: {len(data)} bytes is not a whole number of {2 * CHD_WORDS}-byte frames")
    return bytearray(data)


def parse_burst(text):
    try:
        start, length = (int(v, 10) for v in text.split(":", 1))
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"burst {text!r} is not START:LENGTH") from error
    if start < 0 or length < 1:
        raise argparse.ArgumentTypeError(f"burst {text!r} needs START >= 0 and LENGTH >= 1")
    return start, length


def bad_frames(frames, seed, rate, bursts, first):
    rng = random.Random(seed)
    bad = {f for f in range(frames) if rng.random() * 100.0 < rate}
    for start, length in bursts:
        bad.update(range(start, min(frames, start + length)))
    bad.update(range(min(first, frames)))
    return bad


def cmd_compare(args):
    if args.json and (same_file(args.json, args.ours) or same_file(args.json, args.reference)):
        raise ToolError("the JSON report must not be one of the PCM inputs")
    raw = compare(read_pcm(args.ours), read_pcm(args.reference))
    failures = gate_failures(raw)
    metrics = rounded(raw)
    metrics["gate_failures"] = failures
    if args.json:
        Path(args.json).write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    first = metrics["first_mismatch"]
    print(f"frames {metrics['frames']}, exact {metrics['exact_frames']}"
          + ("" if first is None else f", first mismatch at frame {first['frame']} sample {first['sample']}"))
    print(f"SNR {metrics['snr_db']} dB, segmental SNR {metrics['segmental_snr_db']} dB over "
          f"{metrics['segmental_frames']} frames, level {metrics['level_db']} dB")
    print(f"frames below {MIN_FRAME_SNR_DB} dB: {len(metrics['low_snr_frames'])}; silent reference frames: "
          f"{metrics['silent_reference_frames']} ({len(metrics['silent_frames_not_silent'])} not silent in ours)")
    for failure in failures:
        print(f"gate: {failure}")
    if args.metrics_only:
        return 0
    if args.gate:
        return 1 if failures else 0
    return 0 if metrics["bit_exact"] else 1


def cmd_bfi(args):
    if same_file(args.input, args.output):
        raise ToolError("input and output must be different files")
    if not 0.0 <= args.rate <= 100.0:
        raise ToolError(f"rate {args.rate} is not a percentage")
    data = read_chd(args.input)
    frames = len(data) // (2 * CHD_WORDS)
    bad = bad_frames(frames, args.seed, args.rate, args.burst, args.first)
    for f in sorted(bad):
        struct.pack_into("<h", data, f * 2 * CHD_WORDS, 1)
    Path(args.output).write_bytes(bytes(data))
    print(f"{len(bad)} of {frames} frames marked bad")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="TETRA decoder comparison and bad-frame tools.")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("compare", help="compare decoded PCM with a reference")
    p.add_argument("ours", type=Path)
    p.add_argument("reference", type=Path)
    p.add_argument("--gate", action="store_true", help="exit status from the acceptance thresholds")
    p.add_argument("--metrics-only", action="store_true", help="always exit 0 after reporting")
    p.add_argument("--json", type=Path, help="write the metrics as JSON")
    p = sub.add_parser("bfi", help="copy a .CHD file with bad-frame indicators set")
    p.add_argument("input", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--rate", type=float, default=0.0, help="percentage of frames marked bad at random")
    p.add_argument("--burst", type=parse_burst, action="append", default=[], help="START:LENGTH, repeatable")
    p.add_argument("--first", type=int, default=0, help="mark the first N frames bad")
    args = parser.parse_args(argv)
    try:
        return cmd_compare(args) if args.command == "compare" else cmd_bfi(args)
    except (ToolError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
