#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise tools/quality/tetra_conformance.py, and mbe_tetra_decode when its path is given.

    test_tetra_conformance_tools.py [path/to/mbe_tetra_decode]
"""

import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "quality" / "tetra_conformance.py"
FRAME = 240
CHD_WORDS = 138


def run(args, expect):
    result = subprocess.run([sys.executable, "-I", str(TOOL), *map(str, args)], capture_output=True, timeout=60)
    assert result.returncode == expect, (args, result.returncode, result.stdout, result.stderr)
    return result


def write_pcm(path, samples):
    path.write_bytes(struct.pack(f"<{len(samples)}h", *samples))


def write_chd(path, frames):
    path.write_bytes(b"".join(struct.pack(f"<{CHD_WORDS}h", *frame) for frame in frames))


def load_json(path):
    def reject(constant):
        raise AssertionError(f"non-finite value {constant} in {path}")

    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject)


def speech(frames):
    return [int(8000 * (((i * 37) % 200) - 100) / 100) for i in range(frames * FRAME)]


def test_compare(root):
    ref = root / "ref.pcm"
    ours = root / "ours.pcm"
    report = root / "metrics.json"
    samples = speech(20)
    write_pcm(ref, samples)

    write_pcm(ours, samples)
    run(["compare", ours, ref, "--json", report], 0)
    metrics = load_json(report)
    assert metrics["bit_exact"] and metrics["exact_frames"] == 20 and metrics["first_mismatch"] is None
    assert math.isclose(metrics["level_db"], 0.0, abs_tol=1e-9) and metrics["gate_failures"] == []

    # The report never replaces an input, by name or through a hard link.
    before = ref.read_bytes()
    run(["compare", ours, ref, "--json", ref], 2)
    run(["compare", ours, ref, "--json", ours], 2)
    if hasattr(os, "link"):
        alias = root / "ref-link.pcm"
        try:
            os.link(ref, alias)
        except OSError:
            alias = None
        if alias is not None:
            run(["compare", ours, ref, "--json", alias], 2)
            alias.unlink()
    assert ref.read_bytes() == before and ours.read_bytes() == before
    # A report path through a symlink loop is an ordinary error, not a crash.
    loop = root / "loop.json"
    try:
        os.symlink(loop, loop)
    except (OSError, NotImplementedError):
        loop = None
    if loop is not None:
        result = run(["compare", ours, ref, "--json", loop], 2)
        assert b"Traceback" not in result.stderr, result.stderr
        loop.unlink()

    # One sample off: not bit-exact, but well inside the acceptance thresholds.
    changed = list(samples)
    changed[3 * FRAME + 17] += 1
    write_pcm(ours, changed)
    run(["compare", ours, ref, "--json", report], 1)
    run(["compare", ours, ref, "--metrics-only"], 0)
    run(["compare", ours, ref, "--gate", "--json", report], 0)
    metrics = load_json(report)
    assert metrics["first_mismatch"] == {"frame": 3, "sample": 17}, metrics
    assert metrics["exact_frames"] == 19

    # A loud frame decoded badly after the settling frames fails the gate.
    changed = list(samples)
    for n in range(12 * FRAME, 13 * FRAME):
        changed[n] = -changed[n]
    write_pcm(ours, changed)
    run(["compare", ours, ref, "--gate", "--json", report], 1)
    metrics = load_json(report)
    assert [f["frame"] for f in metrics["low_snr_frames"]] == [12], metrics
    # The same damage inside the first 10 frames is tolerated by the frame floor.
    changed = list(samples)
    for n in range(2 * FRAME, 3 * FRAME):
        changed[n] = 0
    write_pcm(ours, changed)
    run(["compare", ours, ref, "--json", report], 1)
    assert load_json(report)["low_snr_frames"] == []

    # A uniform 0.5 dB level error fails the gate.
    write_pcm(ours, [int(round(v * 10 ** (0.5 / 20))) for v in samples])
    run(["compare", ours, ref, "--gate", "--json", report], 1)
    assert abs(load_json(report)["level_db"] - 0.5) < 0.05

    # The gate judges unrounded values: 0.1002 dB fails although it is reported as 0.1.
    write_pcm(ref, [10000] * (20 * FRAME))
    write_pcm(ours, [10116] * (20 * FRAME))
    run(["compare", ours, ref, "--gate", "--json", report], 1)
    metrics = load_json(report)
    assert math.isclose(metrics["level_db"], 0.1, abs_tol=1e-9) and len(metrics["gate_failures"]) == 1, metrics
    write_pcm(ours, [10115] * (20 * FRAME))
    run(["compare", ours, ref, "--gate", "--json", report], 0)
    write_pcm(ref, samples)

    # Silent reference frames: ours may be within 2 LSB.
    quiet = list(samples)
    quiet[5 * FRAME:6 * FRAME] = [0] * FRAME
    write_pcm(ref, quiet)
    noisy = list(quiet)
    noisy[5 * FRAME + 7] = 2
    write_pcm(ours, noisy)
    run(["compare", ours, ref, "--gate", "--json", report], 0)
    metrics = load_json(report)
    assert metrics["silent_reference_frames"] == 1 and metrics["silent_frames_not_silent"] == []
    noisy[5 * FRAME + 7] = -3
    write_pcm(ours, noisy)
    run(["compare", ours, ref, "--gate", "--json", report], 1)
    assert load_json(report)["silent_frames_not_silent"] == [5]

    # All-silent reference against silent and non-silent output: never NaN.
    write_pcm(ref, [0] * (2 * FRAME))
    write_pcm(ours, [0] * (2 * FRAME))
    run(["compare", ours, ref, "--gate", "--json", report], 0)
    metrics = load_json(report)
    assert metrics["snr_db"] is None and metrics["segmental_snr_db"] is None
    assert math.isclose(metrics["level_db"], 0.0, abs_tol=1e-9)
    write_pcm(ours, [5] * (2 * FRAME))
    run(["compare", ours, ref, "--gate", "--json", report], 1)
    assert load_json(report)["level_db"] is None

    # Byte order: 0x0102 is stored as 02 01.
    ref.write_bytes(bytes([2, 1]) * FRAME)
    write_pcm(ours, [0x0102] * FRAME)
    run(["compare", ours, ref], 0)

    # Malformed inputs.
    write_pcm(ref, samples)
    write_pcm(ours, samples[:-1])
    run(["compare", ours, ref], 2)
    ours.write_bytes(b"")
    run(["compare", ours, ref], 2)
    ours.write_bytes(b"\x00\x00\x00")
    run(["compare", ours, ref], 2)
    run(["compare", root / "missing.pcm", ref], 2)


def chd_frames(count):
    frames = []
    for f in range(count):
        frames.append([0] + [((f + i) * 7 % 5) & 1 for i in range(CHD_WORDS - 1)])
    return frames


def bfi_flags(path):
    data = path.read_bytes()
    return [struct.unpack_from("<h", data, f * 2 * CHD_WORDS)[0] for f in range(len(data) // (2 * CHD_WORDS))]


def test_bfi(root):
    src = root / "in.chd"
    out = root / "out.chd"
    again = root / "again.chd"
    write_chd(src, chd_frames(400))
    run(["bfi", src, out, "--seed", 9, "--rate", 10], 0)
    run(["bfi", src, again, "--seed", 9, "--rate", 10], 0)
    assert out.read_bytes() == again.read_bytes()
    flags = bfi_flags(out)
    assert 20 <= sum(flags) <= 60, sum(flags)
    run(["bfi", src, again, "--seed", 10, "--rate", 10], 0)
    assert out.read_bytes() != again.read_bytes()
    # Only indicator words change.
    original = src.read_bytes()
    marked = out.read_bytes()
    for f in range(400):
        base = f * 2 * CHD_WORDS
        assert marked[base + 2:base + 2 * CHD_WORDS] == original[base + 2:base + 2 * CHD_WORDS]

    run(["bfi", src, out, "--burst", "10:5", "--burst", "398:9", "--first", 2], 0)
    flags = bfi_flags(out)
    assert [f for f, v in enumerate(flags) if v] == [0, 1, 10, 11, 12, 13, 14, 398, 399]

    run(["bfi", src, src, "--rate", 5], 2)
    if hasattr(os, "link"):
        link = root / "in-link.chd"
        try:
            os.link(src, link)
        except OSError:
            link = None
        if link is not None:
            run(["bfi", src, link, "--rate", 5], 2)
    partial = root / "partial.chd"
    partial.write_bytes(src.read_bytes()[:-2])
    run(["bfi", partial, out], 2)
    run(["bfi", src, out, "--burst", "5"], 2)
    run(["bfi", src, out, "--rate", 101], 2)


def test_decoder_tool(decoder, root):
    good = root / "good.chd"
    pcm = root / "good.pcm"
    frames = chd_frames(6)
    write_chd(good, frames)
    subprocess.run([decoder, good, pcm], check=True, capture_output=True, timeout=60)
    assert pcm.stat().st_size == 6 * FRAME * 2

    # A good frame's words must be exactly 0 or 1 before narrowing to bits.
    for value in (2, 256, 257, -1, -32768):
        bad = [list(f) for f in frames]
        bad[3][40] = value
        write_chd(root / "bad.chd", bad)
        result = subprocess.run([decoder, root / "bad.chd", root / "bad.pcm"], capture_output=True, timeout=60)
        assert result.returncode != 0, value
        assert not (root / "bad.pcm").exists()
    # A bad frame's payload is not inspected, and any nonzero indicator marks it bad.
    reference = None
    for indicator in (1, 2, -1):
        bad = [list(f) for f in frames]
        bad[3][0] = indicator
        bad[3][40] = 999
        write_chd(root / "bfi.chd", bad)
        subprocess.run([decoder, root / "bfi.chd", root / "bfi.pcm"], check=True, capture_output=True, timeout=60)
        if indicator == 1:
            reference = (root / "bfi.pcm").read_bytes()
        else:
            assert (root / "bfi.pcm").read_bytes() == reference
    # Partial records and an output that is the input are refused.
    (root / "partial.chd").write_bytes(good.read_bytes()[:-1])
    result = subprocess.run([decoder, root / "partial.chd", root / "partial.pcm"], capture_output=True, timeout=60)
    assert result.returncode != 0
    before = good.read_bytes()
    result = subprocess.run([decoder, good, good], capture_output=True, timeout=60)
    assert result.returncode != 0 and good.read_bytes() == before


def main(argv):
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        test_compare(root)
        test_bfi(root)
        if len(argv) > 1:
            test_decoder_tool(argv[1], root)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
