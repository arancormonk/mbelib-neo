#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Offline evaluator controls, attributable quality reports and CPU benchmarks."""

import argparse
import itertools
import json
import math
import os
import re
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from quality_support import digest, load_json

ROOT = Path(__file__).resolve().parents[2]


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def fresh_directory(path):
    path = Path(path).resolve()
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"output directory must be new or empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def executable_identity(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"not an executable: {path}")
    return {"path": str(path), "sha256": digest(path)}


def controls(evaluator, out_dir):
    import matplotlib
    import numpy as np
    from quality_audio import numpy_lsd, package_versions, read_pcm, write_wav

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    versions = package_versions()
    source = ROOT / "build/quality/corpus/hts1a.raw"
    if not source.is_file():
        raise ValueError(f"controls require fetched reference: {source}")
    evaluator = executable_identity(evaluator)
    output = fresh_directory(out_dir)
    original = read_pcm(source)
    peak = float(np.max(np.abs(original)))
    if peak <= 0:
        raise ValueError("control reference must contain speech energy")
    reference = np.rint(original * (12000 / peak)).astype("<i2")
    rms = float(np.sqrt(np.mean(reference.astype(float) ** 2)))
    fixtures = {
        "reference": reference,
        "gain2": (reference.astype(np.int32) * 2).astype("<i2"),
        "delay64": np.concatenate((np.zeros(64, dtype="<i2"), reference)),
        "advance64": reference[64:],
        "silence": np.zeros(len(reference), dtype="<i2"),
        "short160": reference[:160],
        "clipped": np.clip(reference.astype(np.int32) * 8, -31128, 31128).astype("<i2"),
        "sine75": np.rint(8000 * np.sin(2 * np.pi * 75 * np.arange(8000) / 8000)).astype("<i2"),
    }
    coefficients = (0, 0.1, 0.25, 0.5, 1)
    alternating = 1 - 2 * ((np.arange(len(reference)) // 160) % 2)
    for index, coefficient in enumerate(coefficients):
        values = np.rint(reference.astype(float) + coefficient * rms * alternating)
        if np.any(np.abs(values) > 32767):
            raise ValueError("specified join control unexpectedly clips")
        fixtures[f"join{index}"] = values.astype("<i2")
    for name, values in fixtures.items():
        values.tofile(output / f"{name}.raw")
    write_wav(output / "identity.wav", reference)
    original_hashes = {path.name: digest(path) for path in output.iterdir() if path.is_file()}
    result = {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "argv": sys.argv,
        "evaluator": evaluator,
        "driver_sha256": digest(__file__),
        "audio_helpers_sha256": digest(Path(__file__).with_name("quality_audio.py")),
        "dependencies": versions,
        "reference": {
            "source": str(source),
            "source_sha256": digest(source),
            "scaled_sha256": digest(output / "reference.raw"),
            "samples": len(reference),
            "scaled_peak": int(np.max(np.abs(reference))),
            "full_reference_rms": rms,
        },
        "cases": [],
        "checks": {},
        "numpy_crosschecks": [],
        "failures": [],
        "notes": [
            "C evaluator metrics are authoritative; NumPy is an independent FFT cross-check.",
            "A continuous 75-Hz sine is not a click solely because boundary_index_db is large.",
            "Join/boundary metrics are diagnostics, not standardized click detectors.",
        ],
    }

    def check(name, condition, detail):
        result["checks"][name] = bool(condition)
        if not condition:
            result["failures"].append(f"{name}: {detail}")

    def invoke(name, decoded, ref, lag, expected_exit=0, extra=None):
        report = output / f"{name}.json"
        command = [evaluator["path"], "--decoded", str(output / decoded), "--json", str(report)]
        if ref:
            command += ["--ref", str(output / ref)]
        if lag is not None:
            command += ["--lag-samples", str(lag)]
        if extra:
            command += extra
        completed = subprocess.run(command, capture_output=True, text=True, timeout=60, check=False)
        (output / f"{name}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (output / f"{name}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
        record = {
            "name": name,
            "command": command,
            "returncode": completed.returncode,
            "expected_exit": expected_exit,
            "decoded": decoded,
            "reference": ref,
            "json": report.name if report.is_file() and completed.returncode == 0 else None,
        }
        result["cases"].append(record)
        check(f"{name}_exit", completed.returncode == expected_exit, completed.stderr.strip())
        if completed.returncode:
            return None
        measured = load_json(report)
        record["metrics"] = measured
        check(f"{name}_schema", measured.get("schema_version") == 2 and measured.get("input_api") == "pcm", measured)
        check(
            f"{name}_float_unobserved",
            all(measured.get(key) is None for key in ("float_nonfinite_samples", "float_peak", "float_clip_samples")),
            measured,
        )
        if ref:
            independent = numpy_lsd(read_pcm(output / ref), read_pcm(output / decoded), measured)
            authoritative = measured.get("lsd_db")
            difference = (
                abs(independent - authoritative) if independent is not None and authoritative is not None else None
            )
            matched = (independent is None and authoritative is None) or (difference is not None and difference <= 1e-6)
            result["numpy_crosschecks"].append(
                {
                    "case": name,
                    "c_lsd_db": authoritative,
                    "numpy_lsd_db": independent,
                    "absolute_difference_db": difference,
                    "status": "pass" if matched else "fail",
                }
            )
            check(f"{name}_independent_lsd", matched, difference)
        return measured

    identity = invoke("identity", "identity.wav", "reference.raw", 0)
    gain = invoke("gain2", "gain2.raw", "reference.raw", 0)
    delayed = invoke("delay64", "delay64.raw", "reference.raw", None)
    advanced = invoke("advance64", "advance64.raw", "reference.raw", None)
    silent = invoke("silence_silence", "silence.raw", "silence.raw", 0)
    muted = invoke("speech_zero", "silence.raw", "reference.raw", 0)
    unreferenced = invoke("no_reference", "reference.raw", None, None)
    clipped = invoke("clipping", "clipped.raw", "reference.raw", 0)
    sine = invoke("continuous_sine75", "sine75.raw", "sine75.raw", 0)
    joins = [invoke(f"join_step_{index}", f"join{index}.raw", "reference.raw", 0) for index in range(len(coefficients))]
    if identity and gain:
        check("identity_lsd_zero", abs(identity["lsd_db"]) <= 1e-6, identity["lsd_db"])
        check("gain_normalized_lsd_zero", abs(gain["lsd_db"]) <= 1e-6, gain["lsd_db"])
        check(
            "exact_double_gain_db", abs(gain["level_offset_db"] - 20 * math.log10(2)) <= 1e-6, gain["level_offset_db"]
        )
        check(
            "gain_join_invariant",
            all(abs(identity[key] - gain[key]) <= 1e-6 for key in ("join_dec_db", "join_ref_db", "join_excess_db")),
            "normalized gain changed join metrics",
        )
    if delayed and advanced:
        check(
            "known_delays_recovered",
            abs(delayed["auto_lag_samples"] - 64) <= 8 and abs(advanced["auto_lag_samples"] + 64) <= 8,
            [delayed["auto_lag_samples"], advanced["auto_lag_samples"]],
        )
    for name, measured, state in (
        ("silence", silent, "silent_reference"),
        ("muted", muted, "silent_decoded"),
        ("unreferenced", unreferenced, "no_reference"),
    ):
        if measured:
            check(
                f"{name}_inapplicable",
                measured["reference_state"] == state
                and all(
                    measured[key] is None
                    for key in ("lsd_db", "level_offset_db", "crest_delta_db", "env_corr", "join_excess_db")
                ),
                measured,
            )
    if silent:
        check("silent_alignment_undefined", silent["alignment_corr"] is None, silent["alignment_corr"])
    if clipped:
        expected_rails = int(np.count_nonzero(np.abs(fixtures["clipped"].astype(np.int32)) >= 31128))
        check(
            "deliberate_rails_counted",
            clipped["pcm_peak"] == 31128
            and clipped["pcm_rail_samples"] == expected_rails
            and clipped["pcm_max_rail_run"] > 0,
            clipped,
        )
    if all(row is not None for row in joins):
        contrasts = [row["join_excess_db"] for row in joins]
        check(
            "injected_join_steps_increase",
            all(after > before for before, after in itertools.pairwise(contrasts)),
            contrasts,
        )
        result["join_controls"] = [
            {"coefficient": coefficient, "join_excess_db": row["join_excess_db"], "lsd_db": row["lsd_db"]}
            for coefficient, row in zip(coefficients, joins)
        ]
        figure, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
        axis.plot(coefficients, contrasts, marker="o")
        axis.set(
            xlabel="injected alternating step / full reference RMS",
            ylabel="join excess (dB)",
            title="Evaluator join diagnostic control (fixed lag 0)",
        )
        figure.savefig(output / "join-controls.png", dpi=120)
        plt.close(figure)
    if sine:
        result["continuous_sine75"] = {
            "frequency_hz": 75,
            "amplitude": 8000,
            "boundary_index_db": sine["boundary_index_db"],
            "join_excess_db": sine["join_excess_db"],
            "interpretation": "continuous sine; no click classification inferred from boundary index",
        }

    protected = output / "short_guard.json"
    protected.write_text("existing destination must survive\n", encoding="utf-8")
    before = digest(protected)
    invoke("short_guard", "short160.raw", "short160.raw", 0, expected_exit=2)
    check("short_input_preserves_destination", digest(protected) == before, "too-short measurement truncated output")
    frames = output / "frames.txt"
    frames.write_text(("0" * 88 + "\n") * 2, encoding="ascii")
    hardlink = output / "reference-hardlink.raw"
    os.link(output / "reference.raw", hardlink)
    (output / "path-component").mkdir()
    collision = output / "collision.wav"
    collision.write_text("preserve both output aliases\n", encoding="utf-8")
    hashes = {path: digest(path) for path in (frames, hardlink, collision, output / "reference.raw")}
    decode_command = [evaluator["path"], "--codec", "imbe7200", "--frames", str(frames)]
    alias_cases = [
        (
            "pcm_input_alias",
            [evaluator["path"], "--decoded", str(output / "reference.raw"), "--json", str(output / "reference.raw")],
        ),
        ("hardlink_alias", [evaluator["path"], "--decoded", str(output / "reference.raw"), "--json", str(hardlink)]),
        ("frame_out_alias", decode_command + ["--out", str(frames)]),
        ("frame_json_alias", decode_command + ["--out", str(output / "unused.wav"), "--json", str(frames)]),
        (
            "reference_out_alias",
            decode_command + ["--out", str(output / "reference.raw"), "--ref", str(output / "reference.raw")],
        ),
        (
            "reference_json_alias",
            decode_command
            + [
                "--out",
                str(output / "unused.wav"),
                "--ref",
                str(output / "reference.raw"),
                "--json",
                str(output / "reference.raw"),
            ],
        ),
        ("output_alias", decode_command + ["--out", str(collision), "--json", str(collision)]),
        (
            "normalized_new_output_alias",
            decode_command + ["--out", str(output / "new.wav"), "--json", str(output / "path-component/../new.wav")],
        ),
    ]
    for name, command in alias_cases:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=60, check=False)
        check(
            name,
            completed.returncode == 2 and all(digest(path) == value for path, value in hashes.items()),
            completed.stderr,
        )
        result["cases"].append(
            {"name": name, "command": command, "returncode": completed.returncode, "stderr": completed.stderr}
        )
    check(
        "aliased_new_output_not_created", not (output / "new.wav").exists(), "aliased nonexistent destination created"
    )
    for key, value in (
        ("--codec", "imbe7200"),
        ("--frames", str(frames)),
        ("--out", str(output / "forbidden.wav")),
        ("--seed", "1"),
    ):
        invoke("forbidden_" + key[2:], "reference.raw", None, None, expected_exit=2, extra=[key, value])
    for name, value in (("lag_low", "-161"), ("lag_high", "801"), ("lag_fraction", "0.5")):
        invoke(name, "reference.raw", "reference.raw", value, expected_exit=2)
    check(
        "all_input_bytes_preserved",
        all(digest(output / name) == value for name, value in original_hashes.items()),
        "measurement rewrote PCM fixture",
    )
    result["status"] = "fail" if result["failures"] else "pass"
    result["files"] = {
        str(path.relative_to(output)): digest(path) for path in sorted(output.rglob("*")) if path.is_file()
    }
    write_json(output / "controls.json", result)
    if result["failures"]:
        raise RuntimeError("evaluator controls failed: " + "; ".join(result["failures"]))
    return {
        "out_dir": str(output),
        "status": result["status"],
        "checks": len(result["checks"]),
        "numpy_crosschecks": len(result["numpy_crosschecks"]),
        "continuous_sine75": result.get("continuous_sine75"),
    }


def build_identity(directory):
    cache = directory / "CMakeCache.txt"
    if not cache.is_file():
        raise ValueError(f"matched benchmark requires CMake build identity: {cache}")
    values = {}
    for line in cache.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("#", "//")) and ":" in line and "=" in line:
            name, value = line.split("=", 1)
            values[name.split(":", 1)[0]] = value
    keys = (
        "CMAKE_BUILD_TYPE",
        "CMAKE_C_COMPILER",
        "CMAKE_C_FLAGS",
        "CMAKE_C_FLAGS_RELEASE",
        "MBELIB_ENABLE_SIMD",
        "MBELIB_ENABLE_FAST_MATH",
        "MBELIB_ENABLE_LTO",
        "MBELIB_ENABLE_HARDENING",
        "NOTONES",
    )
    options = {key: values.get(key) for key in keys}
    if options["CMAKE_BUILD_TYPE"] != "Release":
        raise ValueError("benchmark requires matched Release builds")
    compiler = Path(options["CMAKE_C_COMPILER"]).resolve(strict=True)
    version = subprocess.run([str(compiler), "--version"], capture_output=True, text=True, check=True).stdout
    return {
        "cache": str(cache.resolve()),
        "cache_sha256": digest(cache),
        "options": options,
        "compiler": {"path": str(compiler), "sha256": digest(compiler), "version": version},
    }


def benchmark_identity(executable):
    identity = executable_identity(executable)
    directory = Path(identity["path"]).parent
    library = directory / "libmbe-neo.so.2"
    if not library.is_file():
        raise ValueError(f"benchmark library missing beside executable: {library}")
    env = dict(os.environ, LD_LIBRARY_PATH=str(directory))
    env.pop("LD_TRACE_LOADED_OBJECTS", None)
    trace = subprocess.run(
        [identity["path"]],
        env=dict(env, LD_TRACE_LOADED_OBJECTS="1"),
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    rows = re.findall(r"^\s*libmbe-neo\.so\.2\s+=>\s+(.+?)\s+\(0x[0-9a-fA-F]+\)\s*$", trace.stdout, re.MULTILINE)
    if trace.returncode or len(rows) != 1 or not os.path.samefile(rows[0], library):
        raise ValueError(f"benchmark loader did not select expected library: {trace.stdout}{trace.stderr}")
    identity["library"] = {"path": str(library.resolve()), "sha256": digest(library), "source_revision": None}
    capture_path = ROOT / "build/quality/checkpoint-capture.json"
    if capture_path.is_file():

        def visit(value):
            if isinstance(value, dict):
                if value.get("sha256") == identity["library"]["sha256"] and value.get("source_revision"):
                    identity["library"]["source_revision"] = value["source_revision"]
                for child in value.values():
                    visit(child)
            elif isinstance(value, list):
                for child in value:
                    visit(child)

        visit(load_json(capture_path))
        identity["source_capture_sha256"] = digest(capture_path)
    identity["build"] = build_identity(directory)
    identity["loader_trace"] = trace.stdout
    return identity, env


def benchmark(baseline, candidate, output):
    baseline, baseline_env = benchmark_identity(baseline)
    candidate, candidate_env = benchmark_identity(candidate)
    if baseline["build"]["options"] != candidate["build"]["options"]:
        raise ValueError("benchmark build/compiler/DSP options do not match")
    output = Path(output)
    for path in (baseline["path"], candidate["path"], baseline["library"]["path"], candidate["library"]["path"]):
        if output.resolve() == Path(path).resolve() or (output.exists() and os.path.samefile(output, path)):
            raise ValueError("benchmark output aliases an input executable/library")
    result = {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "argv": sys.argv,
        "driver_sha256": digest(__file__),
        "baseline": baseline,
        "candidate": candidate,
        "units": "reported average CPU seconds per inner benchmark repeat",
        "order": [],
        "interpretation": "Throughput cost only, not algorithmic audio delay.",
        "max_cost_increase_percent": 10,
    }
    samples = {"baseline": [], "candidate": []}
    inner_frames = None
    for index in range(14):
        variant = "baseline" if index % 2 == 0 else "candidate"
        identity, env = (baseline, baseline_env) if variant == "baseline" else (candidate, candidate_env)
        completed = subprocess.run(
            [identity["path"]], env=env, capture_output=True, text=True, timeout=300, check=False
        )
        entry = {
            "index": index,
            "variant": variant,
            "command": [identity["path"]],
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }
        result["order"].append(entry)
        matches = re.findall(r"^avg:\s*([0-9.eE+-]+)\s+s\b", completed.stdout, re.MULTILINE)
        totals = re.findall(r"\(frames=(\d+)\)", completed.stdout)
        repeats = len(re.findall(r"^run\s+\d+:", completed.stdout, re.MULTILINE))
        if completed.returncode or len(matches) != 1 or len(totals) != 1 or not repeats:
            result["status"] = "fail"
            write_json(output, result)
            raise RuntimeError(f"benchmark failed or unrecognized output: {identity['path']}")
        average, total = float(matches[0]), int(totals[0])
        frames = total // repeats
        if (
            not math.isfinite(average)
            or average <= 0
            or total % repeats
            or (inner_frames is not None and inner_frames != frames)
        ):
            result["status"] = "fail"
            write_json(output, result)
            raise ValueError("incomparable benchmark time/frame counts")
        inner_frames = frames
        entry.update(avg_seconds=average, inner_repeats=repeats, frames_per_inner_repeat=frames)
        samples[variant].append(average)
    before, after = statistics.median(samples["baseline"]), statistics.median(samples["candidate"])
    change = 100 * (after / before - 1)
    result.update(
        baseline_median_seconds=before,
        candidate_median_seconds=after,
        cost_increase_percent=change,
        frames_per_inner_repeat=inner_frames,
        status="pass" if after <= before * 1.1 else "fail",
    )
    write_json(output, result)
    return {
        "out": str(output),
        "status": result["status"],
        "baseline_median_seconds": before,
        "candidate_median_seconds": after,
        "cost_increase_percent": change,
        "alternating_process_runs": 14,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    control = commands.add_parser("controls", help="exercise real evaluator PCM controls")
    control.add_argument("--evaluator", required=True)
    control.add_argument("--out-dir", required=True)
    report = commands.add_parser("report", help="validate and analyze an immutable quality run")
    report.add_argument("--run-dir", required=True)
    report.add_argument("--listening-results")
    bench = commands.add_parser("benchmark", help="alternate seven matched Release benchmark runs per variant")
    bench.add_argument("--baseline", required=True)
    bench.add_argument("--candidate", required=True)
    bench.add_argument("--out", required=True)
    args = parser.parse_args()
    if os.environ.get("LD_PRELOAD"):
        raise ValueError("LD_PRELOAD must be unset for controlled analysis")
    if args.command == "controls":
        result = controls(args.evaluator, args.out_dir)
    elif args.command == "report":
        from quality_report import report_run

        result = report_run(args.run_dir, args.listening_results)
    else:
        result = benchmark(args.baseline, args.candidate, args.out)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, ImportError, subprocess.SubprocessError) as error:
        print(f"quality analysis: {error}", file=sys.stderr)
        sys.exit(2)
