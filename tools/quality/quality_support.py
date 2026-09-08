#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Stdlib-only contracts shared by the opt-in quality runner and analysis.

Audio metrics remain implemented by mbe_quality_eval. These helpers validate
recorded support, content provenance, fixed calibration and native regressions.
"""

import hashlib
import json
import math
import os
import re
import statistics
import subprocess
import wave
from pathlib import Path

MODES = ("imbe7200", "imbe7100", "ambe2450", "ambe2400")
CALIBRATION_MODES = ("imbe7200", "ambe2450", "ambe2400")
LEGACY_NAMES = ("hts1a", "hts1", "hts2a", "kristoff", "ve9qrp_10s")
INDEPENDENT_LEGACY = ("hts1", "hts2a", "kristoff", "ve9qrp_10s")
ORIGINAL_REVISION = "5fd3f3dd738b76bc9a8446ff4d7ae99e565e35e0"
CHECKPOINT_REVISION = "8da251d105e50815514b1eee994c2c8471469654"
INPUT_APIS = {
    "imbe7200": "mbe_processImbe7200x4400Framef",
    "imbe7100": "mbe_processImbe7100x4400Framef",
    "ambe2450": "mbe_processAmbe3600x2450Framef",
    "ambe2400": "mbe_processAmbe3600x2400Framef",
}
BAND_METRICS = tuple(f"band_delta_db_{band}" for band in ("0_500", "500_1000", "1000_2000", "2000_3000", "3000_4000"))


def digest(path):
    result = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            result.update(block)
    return result.hexdigest()


def _invalid_constant(value):
    raise ValueError(f"nonfinite JSON constant: {value}")


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json(path):
    with open(path, encoding="utf-8") as stream:
        return json.load(stream, parse_constant=_invalid_constant, object_pairs_hook=_unique_object)


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _sha(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def artifact_path(run_dir, relative):
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError(f"artifact must be a nonempty relative path: {relative!r}")
    if ".." in Path(relative).parts:
        raise ValueError(f"artifact path traversal: {relative!r}")
    root = Path(run_dir).resolve()
    path = (root / relative).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise ValueError(f"missing or outside-run artifact: {relative!r}")
    return path


def load_gain_profile(path):
    profile = load_json(path)
    if (
        not isinstance(profile, dict)
        or type(profile.get("schema_version")) is not int
        or profile["schema_version"] != 1
    ):
        raise ValueError("gain profile requires schema_version 1")
    if not _sha(profile.get("calibration_lib_sha256")):
        raise ValueError("gain profile requires a calibration library SHA256")
    sources = profile.get("source_manifest_sha256")
    required = {f"{name}.{mode}.calibration.json" for name in INDEPENDENT_LEGACY for mode in CALIBRATION_MODES}
    if not isinstance(sources, dict) or set(sources) != required or not all(_sha(value) for value in sources.values()):
        raise ValueError("gain profile requires all twelve independent legacy calibration manifest hashes")
    gains = profile.get("gain_adjust_log2")
    if not isinstance(gains, dict) or set(gains) != set(CALIBRATION_MODES):
        raise ValueError("gain profile requires imbe7200, ambe2450 and ambe2400 gains")
    if not all(_finite(value) and -16 <= value <= 16 for value in gains.values()):
        raise ValueError("gain profile gains must be finite numbers in [-16, 16]")
    return profile


def make_gain_profile(manifest, run_dir):
    """Freeze only original-main native legacy calibrations, never held-out gains."""
    corpus = manifest.get("corpus", {})
    calibration = manifest.get("libraries", {}).get("calibration", {})
    if (
        manifest.get("gain_profile") is not None
        or corpus.get("set") != "legacy"
        or corpus.get("level") != "native"
        or calibration.get("source", {}).get("revision") != ORIGINAL_REVISION
    ):
        return None
    references = {row["name"]: row for row in manifest["references"]}
    if set(references) != set(LEGACY_NAMES):
        return None
    encodings = {(row["name"], row["mode"]): row for row in manifest["encodings"]}
    if len(encodings) != len(manifest["encodings"]):
        raise ValueError("duplicate calibration encoding record")
    if not _sha(calibration.get("sha256")):
        raise ValueError("missing original-main calibration library hash")
    gains = {mode: [] for mode in CALIBRATION_MODES}
    sources = {}
    for name in INDEPENDENT_LEGACY:
        for mode in CALIBRATION_MODES:
            encoding = encodings.get((name, mode))
            if not encoding or not encoding.get("calibration"):
                raise ValueError(f"missing legacy calibration: {name}/{mode}")
            path = artifact_path(run_dir, encoding["calibration"])
            record = load_json(path)
            inputs = record.get("inputs", {})
            gain = record.get("gain_adjust_log2")
            level = record.get("baseline_level_offset_db")
            lag = record.get("baseline_lag_samples")
            if (
                record.get("schema_version") != 2
                or inputs.get("mode") != mode
                or inputs.get("baseline_sha256") != calibration["sha256"]
                or inputs.get("raw_sha256") != references[name]["sha256"]
                or record.get("original_samples") != references[name]["samples"]
                or record.get("flush_frames") != manifest["flush_frames"]
                or record.get("frames_sha256") != encoding["sha256"]
                or not _finite(gain)
                or not -16 <= gain <= 16
                or not _finite(level)
                or abs(level + 3) > 0.5
                or not _integer(lag)
                or not 0 <= lag <= 640
                or record.get("target_level_offset_db") != -3
            ):
                raise ValueError(f"invalid original-main calibration: {path}")
            gains[mode].append(gain)
            sources[f"{name}.{mode}.calibration.json"] = digest(path)
    return {
        "schema_version": 1,
        "calibration_lib_sha256": calibration["sha256"],
        "source_manifest_sha256": sources,
        "gain_adjust_log2": {mode: statistics.median(values) for mode, values in gains.items()},
    }


def measurement_errors(record, expected_samples, expected_api=None):
    """Diagnose unsupported speech measurements without replacing null by zero."""
    errors = []
    if not isinstance(record, dict):
        return ["missing metric object"]
    if type(record.get("schema_version")) is not int or record["schema_version"] != 2:
        errors.append("requires evaluator schema 2")
    if record.get("reference_state") != "speech":
        errors.append(f"inapplicable reference_state={record.get('reference_state')!r}")
    if expected_api is not None and record.get("input_api") != expected_api:
        errors.append(f"wrong public input API: {record.get('input_api')!r}")
    metrics = (
        "alignment_corr",
        "level_offset_db",
        "lsd_db",
        "env_corr",
        "crest_ref_db",
        "crest_dec_db",
        "crest_delta_db",
        "join_ref_db",
        "join_dec_db",
        "join_excess_db",
    ) + BAND_METRICS
    for key in metrics:
        if not _finite(record.get(key)):
            errors.append(f"undefined or nonfinite {key}")
    for key in ("lag_samples", "auto_lag_samples"):
        if not _integer(record.get(key)) or not -160 <= record[key] <= 800:
            errors.append(f"invalid {key}")
    if record.get("alignment_at_limit") is not False:
        errors.append("automatic alignment reaches search boundary or is undefined")
    for key in ("reference_samples", "aligned_samples"):
        if not _integer(record.get(key)) or record[key] != expected_samples:
            errors.append(f"{key} does not cover the original reference")
    if expected_samples < 256:
        errors.append("fewer than 256 aligned reference samples")
    for key in ("reference_trim_head", "reference_trim_tail"):
        if record.get(key) != 0 or not _integer(record.get(key)):
            errors.append(f"reference coverage lost at {key}")
    for key in ("active_frames", "spectral_frames", "join_count"):
        if not _integer(record.get(key)) or record[key] <= 0:
            errors.append(f"empty {key} support")
    decoded = record.get("decoded_samples")
    if not _integer(decoded) or decoded < expected_samples:
        errors.append("insufficient decoded samples")
    elif all(_integer(record.get(key)) for key in ("decoded_trim_head", "decoded_trim_tail")):
        if (
            record["decoded_trim_head"] < 0
            or record["decoded_trim_tail"] < 0
            or record["decoded_trim_head"] + record["decoded_trim_tail"] + expected_samples != decoded
        ):
            errors.append("inconsistent decoded support trims")
    else:
        errors.append("missing decoded support trims")
    pcm_input = record.get("input_api") == "pcm"
    if pcm_input:
        if any(record.get(key) is not None for key in ("float_nonfinite_samples", "float_peak", "float_clip_samples")):
            errors.append("PCM input cannot claim float synthesis observations")
    else:
        if record.get("float_nonfinite_samples") != 0 or not _integer(record.get("float_nonfinite_samples")):
            errors.append("nonfinite float synthesis samples or missing observation")
        if not _finite(record.get("float_peak")) or record["float_peak"] < 0:
            errors.append("missing float peak observation")
    for key in ("pcm_peak", "pcm_rail_samples", "pcm_max_rail_run"):
        if not _finite(record.get(key)) or record[key] < 0:
            errors.append(f"missing or negative {key}")
    if (
        _integer(decoded)
        and all(_integer(record.get(key)) for key in ("pcm_rail_samples", "pcm_max_rail_run"))
        and not 0 <= record["pcm_max_rail_run"] <= record["pcm_rail_samples"] <= decoded
    ):
        errors.append("inconsistent PCM rail counts")
    return errors


def load_manifest(run_dir):
    """Validate immutable input/report snapshots, allowing recorded failed rows."""
    run_dir = Path(run_dir)
    manifest = load_json(run_dir / "manifest.json")
    if (
        not isinstance(manifest, dict)
        or type(manifest.get("schema_version")) is not int
        or manifest["schema_version"] != 2
    ):
        raise ValueError("run manifest requires schema_version 2")
    for key in ("libraries", "tools", "corpus", "files"):
        if not isinstance(manifest.get(key), dict):
            raise TypeError(f"manifest missing {key} object")
    for key in ("references", "encodings", "reports", "frame_equivalence", "argv"):
        if not isinstance(manifest.get(key), list):
            raise TypeError(f"manifest missing {key} list")
    if not _integer(manifest.get("flush_frames")) or not 0 <= manifest["flush_frames"] <= 50:
        raise ValueError("invalid run flush count")
    if not _integer(manifest.get("seed")) or not 0 <= manifest["seed"] <= 0xFFFFFFFF:
        raise ValueError("invalid run seed")
    for role in ("baseline", "candidate", "calibration"):
        library = manifest["libraries"].get(role)
        if not isinstance(library, dict) or not _sha(library.get("sha256")) or "source" not in library:
            raise ValueError(f"missing {role} library identity")
    for relative, expected in manifest["files"].items():
        if not _sha(expected) or digest(artifact_path(run_dir, relative)) != expected:
            raise ValueError(f"run artifact hash mismatch: {relative}")

    def indexed(relative, nullable=False):
        if relative is None and nullable:
            return None
        if relative not in manifest["files"]:
            raise ValueError(f"unindexed run artifact: {relative!r}")
        return artifact_path(run_dir, relative)

    references = {}
    for row in manifest["references"]:
        if row["name"] in references or not _integer(row.get("samples")) or row["samples"] <= 0:
            raise ValueError("duplicate or invalid reference record")
        path = indexed(row["path"])
        if path.stat().st_size != row["samples"] * 2 or manifest["files"][row["path"]] != row["sha256"]:
            raise ValueError(f"reference content/count mismatch: {row['name']}")
        references[row["name"]] = row
    if not references:
        raise ValueError("run has no reference inputs")
    if manifest.get("gain_profile") is not None:
        profile = manifest["gain_profile"]
        path = indexed(profile["path"])
        if digest(path) != profile["sha256"]:
            raise ValueError("fixed gain profile hash mismatch")
        load_gain_profile(path)
    encoding_keys = set()
    for row in manifest["encodings"]:
        key = (row["name"], row["mode"])
        if row["name"] not in references or row["mode"] not in CALIBRATION_MODES or key in encoding_keys:
            raise ValueError(f"unknown/duplicate encoding: {key}")
        encoding_keys.add(key)
        if row.get("path") is None:
            continue  # A failed encoding is represented by invalid report rows.
        bits = indexed(row["path"])
        if digest(bits) != row["sha256"]:
            raise ValueError(f"canonical encoding hash mismatch: {key}")
        if row.get("calibration") is not None:
            calibration_path = indexed(row["calibration"])
            calibration_record = load_json(calibration_path)
            if (
                calibration_record.get("schema_version") != 2
                or calibration_record.get("frames_sha256") != row["sha256"]
                or calibration_record.get("gain_adjust_log2") != row["gain_adjust_log2"]
            ):
                raise ValueError(f"calibration snapshot does not bind its encoded stream: {key}")
    reports = {}
    for row in manifest["reports"]:
        key = (row["name"], row["mode"])
        if row["name"] not in references or row["mode"] not in MODES or key in reports:
            raise ValueError(f"unknown/duplicate report: {key}")
        reports[key] = row
        if type(row.get("valid")) is not bool:
            raise ValueError(f"report missing validity: {key}")
        nullable = not row["valid"]
        frame_path = indexed(row.get("frames"), nullable)
        indexed(row["reference"])
        if frame_path is not None:
            width = {"imbe7200": 184, "imbe7100": 168, "ambe2450": 96, "ambe2400": 96}[row["mode"]]
            count = 0
            with frame_path.open("rb") as stream:
                for line in stream:
                    bits = line.rstrip(b"\r\n")
                    if len(bits) != width or any(bit not in (48, 49) for bit in bits):
                        raise ValueError(f"invalid rectangular frame snapshot: {row['frames']}")
                    count += 1
            expected = (references[row["name"]]["samples"] + 159) // 160 + manifest["flush_frames"]
            if count != expected or row["frame_count"] != count or row["frame_sha256"] != digest(frame_path):
                raise ValueError(f"frame completion/hash mismatch: {key}")
        for variant in ("baseline", "candidate"):
            audio_path = indexed(row[variant].get("wav"), nullable)
            metric_path = indexed(row[variant].get("json"), nullable)
            if audio_path is not None:
                with wave.open(str(audio_path), "rb") as audio:
                    if (audio.getnchannels(), audio.getsampwidth(), audio.getframerate(), audio.getcomptype()) != (
                        1,
                        2,
                        8000,
                        "NONE",
                    ):
                        raise ValueError(f"invalid decoded WAV: {audio_path}")
                    if metric_path is not None and audio.getnframes() != load_json(metric_path).get("decoded_samples"):
                        raise ValueError(f"decoded sample count mismatch: {audio_path}")
        sensitivity = row.get("candidate_auto_sensitivity")
        if sensitivity is not None:
            if sensitivity.get("kind") != "automatic_alignment_sensitivity":
                raise ValueError("unlabelled automatic-alignment sensitivity measurement")
            indexed(sensitivity["json"])
    if set(reports) != {(name, mode) for name in references for mode in MODES}:
        raise ValueError("run does not account for all four frame modes for each reference")
    if manifest.get("identical_libraries") != (
        manifest["libraries"]["baseline"]["sha256"] == manifest["libraries"]["candidate"]["sha256"]
    ):
        raise ValueError("incorrect same-library identity declaration")
    if manifest["identical_libraries"] and manifest.get("identity_check", {}).get("status") == "pass":
        for row in reports.values():
            for kind in ("wav", "json"):
                a, b = indexed(row["baseline"][kind]), indexed(row["candidate"][kind])
                equal = digest(a) == digest(b) if kind == "wav" else load_json(a) == load_json(b)
                if not equal:
                    raise ValueError(f"false identity assertion: {row['name']}/{row['mode']}/{kind}")
    equivalence_keys = set()
    for row in manifest["frame_equivalence"]:
        key = (row["name"], row["variant"])
        if row["name"] not in references or row["variant"] not in ("baseline", "candidate") or key in equivalence_keys:
            raise ValueError(f"unknown/duplicate frame equivalence control: {key}")
        equivalence_keys.add(key)
        canonical = indexed(row.get("canonical_wav"), True)
        canonical_metrics = indexed(row.get("canonical_json"), True)
        counts = [load_json(canonical_metrics).get("frames")] if canonical_metrics is not None else []
        for mode in ("imbe7200", "imbe7100"):
            variant = reports[(row["name"], mode)][row["variant"]]
            audio = indexed(variant.get("wav"), True)
            equal = canonical is not None and audio is not None and digest(canonical) == digest(audio)
            if row.get(f"{mode}_equal") is not equal:
                raise ValueError(f"false recorded {mode} frame equivalence: {key}")
            metrics = indexed(variant.get("json"), True)
            if metrics is not None:
                counts.append(load_json(metrics).get("frames"))
        counts_equal = (
            len(counts) == 3 and all(_integer(count) and count > 0 for count in counts) and len(set(counts)) == 1
        )
        if row.get("frames_equal") is not counts_equal:
            raise ValueError(f"false recorded frame-count equivalence: {key}")
    for duplicate in manifest["corpus"].get("duplicates", []):
        first, parent = references[duplicate["name"]], references[duplicate["duplicate_of"]]
        child = indexed(first["path"]).read_bytes()
        with indexed(parent["path"]).open("rb") as stream:
            prefix = stream.read(len(child))
        if child != prefix or duplicate.get("relation") != "exact_prefix":
            raise ValueError("unproven duplicate-reference relation")
    return manifest


def _cache_value(path, key):
    if not path.is_file():
        return None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(key + ":") and "=" in line:
            return line.split("=", 1)[1]
    return None


def run_correctness_tests(candidate_libdir, out_dir):
    """Bind native regression evidence to matching build artifacts, including ASan.

    Missing tests or a different test source are an unavailable prerequisite,
    never a pass. Static private-DSP tests bind the archive; public frame/noise
    tests additionally verify the actual loader-selected shared library.
    """
    build = Path(candidate_libdir).resolve()
    root = Path(__file__).resolve().parents[2]
    out_dir = Path(out_dir)
    log_dir = out_dir / "correctness"
    log_dir.mkdir(parents=True, exist_ok=True)
    source = _cache_value(build / "CMakeCache.txt", "CMAKE_HOME_DIRECTORY")
    tests = []
    statuses = {}
    if os.environ.get("LD_PRELOAD"):
        raise ValueError("LD_PRELOAD must be unset for correctness evidence")
    env = dict(os.environ, LD_LIBRARY_PATH=str(build))
    env.pop("LD_TRACE_LOADED_OBJECTS", None)
    loader_row = re.compile(r"^\s*(?:(\S+)\s+=>\s+)?(.+?)\s+\(0x[0-9a-fA-F]+\)\s*$")
    for name in ("test_params", "test_frame_paths", "test_noise_determinism"):
        executable = build / name
        record = {"name": name, "command": [str(executable)], "status": "not_established"}
        tests.append(record)
        statuses[name] = "not_established"
        expected_source = root / "tests" / f"{name}.c"
        built_source = Path(source) / "tests" / f"{name}.c" if source else None
        if not executable.is_file() or built_source is None or not built_source.is_file():
            record["reason"] = "matching native regression executable/source unavailable"
            continue
        if digest(built_source) != digest(expected_source):
            record["reason"] = "build test source does not contain this branch's regression contracts"
            continue
        library = build / ("libmbe-neo.a" if name == "test_params" else "libmbe-neo.so.2")
        if not library.is_file() or executable.stat().st_mtime_ns < max(
            library.stat().st_mtime_ns, built_source.stat().st_mtime_ns
        ):
            record["reason"] = "test artifact is older than its current library/source; rebuild required"
            continue
        record.update(
            executable_sha256=digest(executable),
            source_sha256=digest(built_source),
            library=str(library.resolve()),
            library_sha256=digest(library),
            binding="matching-build-static-archive" if name == "test_params" else "loader-samefile",
        )
        trace = subprocess.run(
            [str(executable)],
            env=dict(env, LD_TRACE_LOADED_OBJECTS="1"),
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        trace_name = f"correctness/{name}.loader.txt"
        (out_dir / trace_name).write_text(trace.stdout + trace.stderr, encoding="utf-8")
        record["loader_trace"] = trace_name
        selected = []
        for line in trace.stdout.splitlines():
            match = loader_row.fullmatch(line)
            if match and match.group(1) == "libmbe-neo.so.2":
                selected.append(match.group(2))
        binding_ok = trace.returncode == 0 and (
            (name == "test_params" and not selected)
            or (name != "test_params" and len(selected) == 1 and os.path.samefile(selected[0], library))
        )
        if not binding_ok:
            record["reason"] = "could not verify regression executable library binding"
            continue
        try:
            result = subprocess.run(
                [str(executable)], env=env, capture_output=True, text=True, timeout=240, check=False
            )
            stdout, stderr, returncode = result.stdout, result.stderr, result.returncode
        except subprocess.TimeoutExpired as exc:
            stdout, stderr, returncode = exc.stdout or b"", exc.stderr or b"", None
            stdout = stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout
            stderr = stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr
            stderr += "\nnative correctness check timed out\n"
        record["stdout"] = f"correctness/{name}.stdout.txt"
        record["stderr"] = f"correctness/{name}.stderr.txt"
        (out_dir / record["stdout"]).write_text(stdout, encoding="utf-8")
        (out_dir / record["stderr"]).write_text(stderr, encoding="utf-8")
        record["returncode"] = returncode
        record["status"] = "pass" if returncode == 0 else "fail"
        statuses[name] = record["status"]
    gates = {
        "nonnegative_attenuation": statuses["test_params"],
        "tone_dispatch": statuses["test_frame_paths"],
        "warm_fft_wola": "pass"
        if statuses["test_params"] == statuses["test_noise_determinism"] == "pass"
        else "fail"
        if "fail" in (statuses["test_params"], statuses["test_noise_determinism"])
        else "not_established",
    }
    status = (
        "fail"
        if "fail" in gates.values()
        else "pass"
        if all(value == "pass" for value in gates.values())
        else "not_established"
    )
    return {"status": status, "gates": gates, "tests": tests}
