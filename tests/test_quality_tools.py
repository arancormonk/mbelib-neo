#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise quality-tool alignment, frame widths and output permissions."""

import json
import os
from pathlib import Path
import stat
import struct
import subprocess
import sys
import tempfile


def run(executable, *args, success=True):
    result = subprocess.run(
        [executable, *map(str, args)], capture_output=True, timeout=30,
    )
    if success:
        assert result.returncode == 0, result.stderr
    else:
        assert result.returncode != 0, result.stdout
    return result


def private_output(path):
    if os.name == "posix":
        assert stat.S_IMODE(path.stat().st_mode) == 0o600, path


def write_pcm(path, samples):
    path.write_bytes(struct.pack(f"<{len(samples)}h", *samples))


def verify(evaluator, reframer, root):
    reference = root / "reference.raw"
    decoded = root / "decoded.raw"
    report = root / "metrics.json"
    # A changing envelope gives a unique alignment peak. Check both lag signs
    # and identity; all three overlaps contain exactly the reference PCM.
    samples = [(200 + ((i // 160 * 997) % 8000)) * (1 if i % 2 else -1) for i in range(160 * 80)]
    write_pcm(reference, samples)
    for lag in (0, 80, -80):
        write_pcm(decoded, [0] * lag + samples if lag >= 0 else samples[-lag:])
        run(evaluator, "--decoded", decoded, "--ref", reference, "--json", report)
        metrics = json.loads(report.read_text())
        assert metrics["lag_samples"] == lag, metrics
        assert metrics["alignment_corr"] > 0.99, metrics
        assert abs(metrics["lsd_db"]) < 1e-10, metrics
        private_output(report)

    write_pcm(reference, [0] * 640)
    write_pcm(decoded, [0] * 640)
    run(evaluator, "--decoded", decoded, "--ref", reference, "--json", report)
    assert json.loads(report.read_text())["alignment_corr"] is None

    frames = root / "frames.txt"
    wav = root / "decoded.wav"
    # Exercise every accepted Data/Frame width and each neighboring bad width.
    widths = {"imbe7200": (88, 184), "imbe7100": (168,), "ambe2400": (49, 96), "ambe2450": (49, 96)}
    for codec, accepted in widths.items():
        for width in accepted:
            frames.write_text(("0" * width + "\n") * 2)
            run(evaluator, "--codec", codec, "--frames", frames, "--out", wav, "--json", report)
            assert len(wav.read_bytes()) == 44 + 2 * 160 * 2
            private_output(wav)
            for bad in (width - 1, width + 1):
                frames.write_text("0" * bad + "\n")
                result = run(evaluator, "--codec", codec, "--frames", frames, "--out", wav, success=False)
                assert b"length/codec mismatch" in result.stderr

    parameters = root / "parameters.txt"
    reframed = root / "reframed.txt"
    for codec, width in (("imbe7200", 184), ("imbe7100", 168), ("ambe2400", 96), ("ambe2450", 96)):
        parameters.write_text(("0" * (88 if codec.startswith("imbe") else 49) + "\n") * 2)
        run(reframer, "--codec", codec, "--in", parameters, "--out", reframed)
        assert [len(row) for row in reframed.read_text().splitlines()] == [width, width]
        private_output(reframed)
        before = reframed.read_bytes()
        parameters.write_text(parameters.read_text() + "invalid\n")
        run(reframer, "--codec", codec, "--in", parameters, "--out", reframed, success=False)
        assert reframed.read_bytes() == before


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="mbe-quality-tools-") as directory:
        # Children inherit this permissive mask; outputs must still be private.
        old_mask = os.umask(0)
        try:
            verify(sys.argv[1], sys.argv[2], Path(directory))
        finally:
            os.umask(old_mask)
