#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Exercise quality-tool alignment, frame widths and output permissions."""

import importlib.util
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
        # The top band overlaps 3000_4000 and covers what LSD (to 3687.5 Hz) leaves out.
        assert abs(metrics["band_delta_db_3500_4000"]) < 1e-10, metrics
        assert abs(metrics["band_delta_db_0_250"]) < 1e-10, metrics
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

    # Per-frame model dump: decoded, prediction-history and synthesized records.
    params = root / "params.jsonl"
    for codec, width in (("imbe7200", 88), ("imbe7200", 184), ("ambe2450", 49), ("ambe2400", 49), ("ambe2400", 96)):
        frames.write_text(("0" * width + "\n") * 3)
        run(evaluator, "--codec", codec, "--frames", frames, "--out", wav)
        plain = wav.read_bytes()
        run(evaluator, "--codec", codec, "--frames", frames, "--out", wav, "--params", params)
        assert wav.read_bytes() == plain, "dumping parameters must not change the PCM"
        private_output(params)
        records = [json.loads(line) for line in params.read_text().splitlines()]
        assert [record["f"] for record in records] == [0, 1, 2], records
        for record in records:
            assert set(record) == {"f", "flags", "errors", "decoded", "history", "synth"}, record
            for key in ("history", "synth"):
                model = record[key]
                assert 9 <= model["L"] <= 56 and model["w0"] > 0, record
                assert all(len(model[field]) == model["L"] for field in ("Vl", "Ml", "log2Ml")), record
            decoded = record["decoded"]
            assert decoded["status"] == 0 and decoded["L"] == record["synth"]["L"], record
            assert abs(decoded["w0"] - record["history"]["w0"]) < 1e-9, record
            assert decoded["log2Ml"] == record["history"]["log2Ml"], record
    result = run(evaluator, "--decoded", reference, "--params", params, success=False)
    assert b"--params" in result.stderr
    run(evaluator, "--codec", "ambe2400", "--frames", frames, "--out", wav, "--params", frames, success=False)

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

    # DVSI hard-decision vectors: whole binary frames only; a truncated tail,
    # an unsupported codec or a second input leaves the output untouched.
    vector = root / "vector.bit"
    for codec, size, width in (("imbe7200", 18, 184), ("ambe2450", 9, 96), ("ambe2400", 9, 96)):
        vector.write_bytes(bytes(range(size)) * 3)
        run(reframer, "--codec", codec, "--from-dvsi", vector, "--out", reframed)
        assert [len(row) for row in reframed.read_text().splitlines()] == [width] * 3
        private_output(reframed)
        before = reframed.read_bytes()
        vector.write_bytes(bytes(range(size)) * 3 + b"\x00")
        result = run(reframer, "--codec", codec, "--from-dvsi", vector, "--out", reframed, success=False)
        assert b"truncated DVSI frame" in result.stderr
        assert reframed.read_bytes() == before
    # DVSI's 4-bit soft-decision files (*_sd.bit) have a whole number of hard
    # frames too, so only the explicit rejection keeps them from importing as garbage.
    soft = root / "dam_e1_sd.bit"
    soft.write_bytes(bytes(36) * 2)
    before = reframed.read_bytes()
    result = run(reframer, "--codec", "ambe2400", "--from-dvsi", soft, "--out", reframed, success=False)
    assert b"soft-decision" in result.stderr
    assert reframed.read_bytes() == before
    vector.write_bytes(bytes(168))
    run(reframer, "--codec", "imbe7100", "--from-dvsi", vector, "--out", reframed, success=False)
    run(reframer, "--codec", "ambe2400", "--from-dvsi", vector, "--in", parameters, "--out", reframed, success=False)
    vector.write_bytes(b"")
    run(reframer, "--codec", "ambe2400", "--from-dvsi", vector, "--out", reframed, success=False)


def verify_dvsi_selection():
    """The fetcher takes every supported mode, its no-FEC variant and the transcodes between them."""
    path = Path(__file__).resolve().parents[1] / "tools" / "quality" / "fetch_dvsi_vectors.py"
    spec = importlib.util.spec_from_file_location("fetch_dvsi_vectors", path)
    assert spec is not None and spec.loader is not None, path
    fetcher = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fetcher)
    wanted = [
        "tv-rc/dam.pcm", "tv-rc/cmprc.txt", "tv-rc/dstar/dam.bit", "tv-rc/dstar/dam_e1_sd.pcm",
        "tv-rc/p25/dam.pcm", "tv-rc/p25_nofec/dam.bit", "tv-rc/r33/dam.bit", "tv-rc/r34/dam.pcm",
        "tv-rc/dstar/p25/dam.bit", "tv-rc/r33/p25/alltone.bit", "tv-rc/p25/p25_nofec/sine0_4k.bit",
        "tv-rc/p25_nofec/dstar/dam.bit", "tv-rc/r34/r33/dam.bit",
    ]
    unwanted = [
        "tv-rc/r0/dam.bit", "tv-rc/r62/dam.pcm", "tv-rc/dstar/r0/dam.bit", "tv-rc/r33/r63/dam.bit",
        "tv-rc/r0/p25/dam.bit", "tv-rc/dstar/p25/extra/dam.bit", "tv-rc/readme.txt",
    ]
    def selected(name):
        return any(fetcher.member_matches(name, pattern) for pattern in fetcher.DEFAULT_PATTERNS)
    assert [name for name in wanted if not selected(name)] == []
    assert [name for name in unwanted if selected(name)] == []


if __name__ == "__main__":
    verify_dvsi_selection()
    with tempfile.TemporaryDirectory(prefix="mbe-quality-tools-") as directory:
        # Children inherit this permissive mask; outputs must still be private.
        old_mask = os.umask(0)
        try:
            verify(sys.argv[1], sys.argv[2], Path(directory))
        finally:
            os.umask(old_mask)
