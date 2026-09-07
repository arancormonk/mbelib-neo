#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Calibrate encoded gain using only the frozen baseline, then cache exact bits.

PCM analysis stays at its original level. OP25 attenuates its log2 spectral-gain
quantizer instead. A -3 dB decoded RMS target leaves headroom before the public
API's int16 clipping. Never recalibrate against the candidate: that would change
both the bitstream and synthesis in the A/B. Metadata binds caches to all inputs.
"""

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def run(command, env=None):
    result = subprocess.run(command, env=env, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"{' '.join(map(str, command))}:\n{result.stderr}")


def main():
    if len(sys.argv) != 7:
        raise ValueError("usage: calibrate_encoder.py encoder evaluator baseline_libdir raw mode frames")
    encoder, evaluator, baseline, raw, mode, frames = sys.argv[1:]
    frames = Path(frames)
    metadata = frames.with_suffix(".calibration.json")
    signature = {
        "raw_sha256": digest(raw),
        "encoder_sha256": digest(encoder),
        "evaluator_sha256": digest(evaluator),
        "baseline_sha256": digest(Path(baseline) / "libmbe-neo.so.2"),
        "calibrator_sha256": digest(__file__),
        "mode": mode,
    }
    if frames.is_file() and metadata.is_file():
        cached = json.loads(metadata.read_text(encoding="utf-8"))
        if cached.get("inputs") == signature and cached.get("frames_sha256") == digest(frames):
            return
    env = dict(os.environ, LD_LIBRARY_PATH=str(Path(baseline).resolve()))
    frames.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".calibrate-", dir=frames.parent) as temp:
        tmp = Path(temp)
        bits, wav, report = tmp / "frames.txt", tmp / "decoded.wav", tmp / "metrics.json"
        gain = math.log2(7)
        for attempt in range(6):
            run([encoder, "--mode", mode, "--in", raw, "--out", str(bits), "--gain-adjust", str(gain)])
            run([evaluator, "--codec", mode, "--frames", str(bits), "--out", str(wav),
                 "--ref", raw, "--json", str(report)], env)
            measured = json.loads(report.read_text(encoding="utf-8"))
            offset = measured["level_offset_db"]
            if not math.isfinite(offset):
                raise ValueError(f"non-finite baseline level for {raw} ({mode})")
            if abs(offset + 3) <= 0.5:
                break
            gain += (offset + 3) / (20 * math.log10(2))
        else:
            raise ValueError(f"baseline gain calibration did not converge: {raw} ({mode}), {offset:.3f} dB")
        if not 0 <= measured["lag_samples"] <= 640:
            raise ValueError(f"baseline alignment out of range: {raw} ({mode}), {measured['lag_samples']}")
        record = {
            "inputs": signature,
            "frames_sha256": digest(bits),
            "gain_adjust_log2": gain,
            "baseline_level_offset_db": offset,
            "baseline_lag_samples": measured["lag_samples"],
            "target_level_offset_db": -3,
        }
        manifest = tmp / "calibration.json"
        manifest.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        os.replace(bits, frames)
        os.replace(manifest, metadata)
        print(f"calibrated {Path(raw).name} {mode}: gain_adjust={gain:.6f}, "
              f"level={offset:.3f} dB, lag={measured['lag_samples']}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(f"calibration: {exc}", file=sys.stderr)
        sys.exit(2)
