#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""A quality A/B run whose clip metrics cannot be read at the end fails cleanly.

Runner.finalize() re-reads each report's metrics for the table and the per-mode
means. A file that is missing or no longer parses must invalidate its report,
failing the run and keeping the clip out of the means, instead of counting a valid
clip with no numbers (which crashed the means).
"""

import json
import os
import sys
import tempfile
from pathlib import Path


def build_runner(ab, output):
    runner = ab.Runner(None, output)
    manifest = runner.manifest
    manifest["references"] = [{"name": "clip", "path": "refs/clip.wav", "duplicate_of": None}]
    manifest["corpus"] = {"independent_references": ["clip"], "duplicates": []}
    manifest["correctness_tests"] = {"status": "pass", "gates": {}, "tests": []}
    manifest["frame_equivalence"] = [
        {"name": "clip", "variant": variant, "imbe7200_equal": True, "imbe7100_equal": True,
         "frames_equal": True, "invalid_reasons": []}
        for variant in ("baseline", "candidate")
    ]
    for mode in ab.MODES:
        report = {"name": "clip", "mode": mode, "valid": True, "invalid_reasons": []}
        for variant in ("baseline", "candidate"):
            relative = f"metrics/clip.{mode}.{variant}.json"
            path = output / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({metric: 1.0 for metric in ab.SUMMARY_METRICS}), encoding="utf-8")
            report[variant] = {"wav": None, "json": relative}
        manifest["reports"].append(report)
    return runner


def finalize(ab, damage):
    with tempfile.TemporaryDirectory(prefix="mbe_quality_finalize_") as scratch:
        output = Path(scratch)
        runner = build_runner(ab, output)
        damage(output)
        rc = runner.finalize()
        acceptance = json.loads((output / "acceptance.json").read_text(encoding="utf-8"))
        table = json.loads((output / "table.json").read_text(encoding="utf-8"))
        return rc, acceptance, table


def expect_invalid(ab, label, damage, mode, variant, reason):
    rc, acceptance, table = finalize(ab, damage)
    assert rc == 1, (label, rc)
    assert acceptance["correctness"] == "fail", (label, acceptance)
    assert acceptance["evidence"]["supported_full_reference_measurements"] == "fail", (label, acceptance)
    assert any(f"clip/{mode}: {variant}: {reason}" in text for text in acceptance["reasons"]), (label, acceptance)
    rows = {row["mode"]: row for row in table["rows"]}
    assert rows[mode]["valid"] is False and rows[mode][variant] is None, (label, rows[mode])
    assert table["means"][mode]["independent_valid_clips"] == 0, (label, table["means"][mode])
    for other in ab.MODES:
        if other != mode:
            assert table["means"][other]["independent_valid_clips"] == 1, (label, other)
    print(f"PASS {label}: run fails and {mode} leaves the means")


def main():
    sys.path.insert(0, os.path.abspath(sys.argv[1]))
    import run_quality_ab as ab

    rc, acceptance, table = finalize(ab, lambda output: None)
    assert rc == 0 and acceptance["correctness"] == "pass", (rc, acceptance)
    assert all(table["means"][mode]["independent_valid_clips"] == 1 for mode in ab.MODES), table["means"]
    print("PASS control: readable metrics pass")

    expect_invalid(ab, "missing metrics", lambda output: (output / "metrics/clip.imbe7200.baseline.json").unlink(),
                   "imbe7200", "baseline", "metrics unreadable at run end")
    expect_invalid(ab, "corrupt metrics",
                   lambda output: (output / "metrics/clip.ambe2450.candidate.json").write_text("{", encoding="utf-8"),
                   "ambe2450", "candidate", "metrics unreadable at run end")
    return 0


if __name__ == "__main__":
    sys.exit(main())
