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
import subprocess
import sys
import tempfile
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def run(command, env=None):
    result = subprocess.run(command, env=env, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"{' '.join(map(str, command))}:\n{result.stderr}")


def validate_measurement(measured, original_samples):
    if not isinstance(measured, dict) or measured.get("schema_version") != 2:
        raise ValueError("baseline evaluator must report schema_version 2")
    if measured.get("reference_state") != "speech":
        raise ValueError("baseline calibration requires speech; silent controls must bypass calibration")

    def finite_number(key):
        value = measured.get(key)
        if type(value) not in (int, float) or not math.isfinite(value):
            raise ValueError(f"baseline measurement {key} must be non-null and finite")
        return value

    def integer(key):
        value = measured.get(key)
        if type(value) is not int:
            raise ValueError(f"baseline measurement {key} must be an integer")
        return value

    offset = finite_number("level_offset_db")
    finite_number("alignment_corr")
    lag = integer("lag_samples")
    auto_lag = integer("auto_lag_samples")
    if not 0 <= lag <= 640 or auto_lag != lag:
        raise ValueError(f"baseline automatic alignment out of range: lag={lag}, auto_lag={auto_lag}")
    if measured.get("alignment_at_limit") is not False:
        raise ValueError("baseline automatic alignment reached a search boundary")
    if integer("float_nonfinite_samples") != 0:
        raise ValueError("baseline decode contains nonfinite float PCM")
    for key in ("active_frames", "spectral_frames", "join_count"):
        if integer(key) <= 0:
            raise ValueError(f"baseline measurement has no {key} support")
    aligned = integer("aligned_samples")
    if aligned < 256 or integer("decoded_samples") < aligned:
        raise ValueError("baseline measurement has insufficient aligned sample support")
    if (
        integer("reference_samples") != original_samples
        or aligned != original_samples
        or integer("reference_trim_head") != 0
        or integer("reference_trim_tail") != 0
    ):
        raise ValueError("baseline decode does not cover the full original reference; check flush_frames")
    return offset


def main():
    if len(sys.argv) != 8:
        raise ValueError(
            "usage: calibrate_encoder.py encoder evaluator calibration_libdir raw mode frames flush_frames"
        )
    encoder, evaluator, baseline, raw, mode, frames, flush_option = sys.argv[1:]
    flush_frames = int(flush_option)
    if not 0 <= flush_frames <= 50:
        raise ValueError("flush_frames must be an integer within [0, 50]")
    if mode not in ("imbe7200", "ambe2450", "ambe2400"):
        raise ValueError(f"unsupported calibration mode: {mode}")
    original_bytes = Path(raw).stat().st_size
    if original_bytes == 0 or original_bytes % 2:
        raise ValueError("calibration input must contain a nonempty, even number of s16le PCM bytes")
    original_samples = original_bytes // 2
    frames = Path(frames)
    metadata = frames.with_suffix(".calibration.json")
    signature = {
        "raw_sha256": digest(raw),
        "encoder_sha256": digest(encoder),
        "evaluator_sha256": digest(evaluator),
        "baseline_sha256": digest(Path(baseline) / "libmbe-neo.so.2"),
        "calibrator_sha256": digest(__file__),
        "mode": mode,
        "flush_frames": flush_frames,
        "original_samples": original_samples,
    }
    if frames.is_file() and metadata.is_file():
        try:
            cached = json.loads(metadata.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            cached = None
        if (
            isinstance(cached, dict)
            and cached.get("schema_version") == 2
            and cached.get("flush_frames") == flush_frames
            and cached.get("original_samples") == original_samples
            and cached.get("inputs") == signature
            and cached.get("frames_sha256") == digest(frames)
        ):
            return
    env = dict(os.environ, LD_LIBRARY_PATH=str(Path(baseline).resolve()))
    frames.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".calibrate-", dir=frames.parent) as temp:
        tmp = Path(temp)
        bits, wav, report = tmp / "frames.txt", tmp / "decoded.wav", tmp / "metrics.json"
        gain = math.log2(7)
        for attempt in range(6):
            run(
                [
                    encoder,
                    "--mode",
                    mode,
                    "--in",
                    raw,
                    "--out",
                    str(bits),
                    "--gain-adjust",
                    str(gain),
                    "--flush-frames",
                    str(flush_frames),
                ]
            )
            run(
                [
                    evaluator,
                    "--codec",
                    mode,
                    "--frames",
                    str(bits),
                    "--out",
                    str(wav),
                    "--ref",
                    raw,
                    "--json",
                    str(report),
                ],
                env,
            )
            measured = json.loads(report.read_text(encoding="utf-8"))
            offset = validate_measurement(measured, original_samples)
            if abs(offset + 3) <= 0.5:
                break
            gain += (offset + 3) / (20 * math.log10(2))
        else:
            raise ValueError(f"baseline gain calibration did not converge: {raw} ({mode}), {offset:.3f} dB")
        record = {
            "schema_version": 2,
            "flush_frames": flush_frames,
            "original_samples": original_samples,
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
        print(
            f"calibrated {Path(raw).name} {mode}: gain_adjust={gain:.6f}, "
            f"level={offset:.3f} dB, lag={measured['lag_samples']}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(f"calibration: {exc}", file=sys.stderr)
        sys.exit(2)
