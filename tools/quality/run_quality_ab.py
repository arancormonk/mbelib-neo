#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Run immutable, content-bound four-frame-mode synthesis comparisons on Linux."""

import argparse
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import wave
from datetime import datetime, timezone
from pathlib import Path

from quality_support import (
    digest,
    load_gain_profile,
    measurement_errors,
    run_correctness_tests,
)

SEED = 0x12345678
FLUSH = 5
MODES = ("imbe7200", "imbe7100", "ambe2450", "ambe2400")
ENCODED_MODES = ("imbe7200", "ambe2450", "ambe2400")
APIS = {
    "imbe7200": "mbe_processImbe7200x4400Framef",
    "imbe7100": "mbe_processImbe7100x4400Framef",
    "ambe2450": "mbe_processAmbe3600x2450Framef",
    "ambe2400": "mbe_processAmbe3600x2400Framef",
}
LEGACY = {"hts1", "hts1a", "hts2a", "kristoff", "ve9qrp_10s"}
SUMMARY_METRICS = ("lsd_db", "env_corr", "crest_delta_db", "join_excess_db", "pcm_rail_samples", "pcm_max_rail_run")


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_json(path):
    def reject(value):
        raise ValueError(f"nonfinite JSON number: {value}")

    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject)


def signature_hash(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def file_identity(path):
    resolved = Path(path).resolve(strict=True)
    return {"requested": str(path), "resolved": str(resolved), "sha256": digest(resolved)}


def pcm_identity(path):
    with wave.open(str(path), "rb") as stream:
        parameters = (
            stream.getnchannels(),
            stream.getsampwidth(),
            stream.getframerate(),
            stream.getnframes(),
            stream.getcomptype(),
        )
        pcm = stream.readframes(stream.getnframes())
    return parameters, hashlib.sha256(pcm).hexdigest()


def frame_count(path, width):
    count = 0
    with path.open("rb") as stream:
        for line in stream:
            bits = line.rstrip(b"\r\n")
            if len(bits) != width or any(bit not in (48, 49) for bit in bits):
                raise ValueError(f"invalid {width}-bit frame at {path}:{count + 1}")
            count += 1
    if not count:
        raise ValueError(f"empty encoded stream: {path}")
    return count


def safe_child(directory, relative):
    if not isinstance(relative, str) or Path(relative).is_absolute():
        raise ValueError(f"invalid corpus-relative path: {relative!r}")
    path = (directory / relative).resolve(strict=True)
    if not path.is_relative_to(directory.resolve()):
        raise ValueError(f"corpus artifact escapes corpus directory: {relative!r}")
    return path


class Runner:
    def __init__(self, args, output):
        self.args = args
        self.output = output
        self.commands = []
        self.profile = None
        self.tools = {}
        self.libdirs = {}
        self.manifest = {
            "schema_version": 2,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "argv": [os.environ.get("MBE_QUALITY_RUNNER_ARGV0", sys.argv[0]), *sys.argv[1:]],
            "seed": SEED,
            "flush_frames": FLUSH,
            "libraries": {},
            "comparison_points": {},
            "tools": {},
            "gain_profile": None,
            "corpus": {},
            "references": [],
            "encodings": [],
            "reports": [],
            "frame_equivalence": [],
            "identical_libraries": False,
            "identity_check": {"status": "not_established", "failures": []},
            "correctness_tests": {"status": "not_established", "gates": {}, "tests": []},
            "failures": [],
            "files": {},
        }

    def relative(self, path):
        return str(path.relative_to(self.output))

    def snapshot(self, source, relative):
        destination = self.output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        before = digest(source)
        shutil.copyfile(source, destination)
        if digest(destination) != before or digest(source) != before:
            raise ValueError(f"input changed while snapshotting: {source}")
        return destination

    def command(self, argv, label, libdir=None):
        argv = [str(value) for value in argv]
        environment = dict(os.environ)
        # Never inherit a request to trace instead of execute, or inherited library search paths.
        environment.pop("LD_TRACE_LOADED_OBJECTS", None)
        environment.pop("LD_LIBRARY_PATH", None)
        if libdir is not None:
            environment["LD_LIBRARY_PATH"] = str(libdir)
        prefix = self.output / "logs" / f"{len(self.commands):05d}-{label}"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        try:
            result = subprocess.run(argv, env=environment, capture_output=True, check=False)
            stdout, stderr, returncode = result.stdout, result.stderr, result.returncode
        except OSError as exc:
            stdout, stderr, returncode = b"", str(exc).encode(), 127
        out_path, err_path = prefix.with_suffix(".stdout"), prefix.with_suffix(".stderr")
        out_path.write_bytes(stdout)
        err_path.write_bytes(stderr)
        self.commands.append(
            {
                "argv": argv,
                "library_directory": str(libdir) if libdir else None,
                "returncode": returncode,
                "stdout": self.relative(out_path),
                "stderr": self.relative(err_path),
            }
        )
        if returncode:
            raise RuntimeError(f"{label} exited {returncode}; see {self.relative(err_path)}")
        return stdout.decode("utf-8", errors="replace")

    def loader_check(self, libdir, role):
        environment = dict(os.environ, LD_LIBRARY_PATH=str(libdir), LD_TRACE_LOADED_OBJECTS="1")
        result = subprocess.run([self.tools["evaluator"]], env=environment, capture_output=True, text=True, check=False)
        log = self.output / "identity" / f"{role}-loader.txt"
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(result.stdout + result.stderr, encoding="utf-8")
        if result.returncode:
            raise ValueError(f"{role} loader tracing failed; see {self.relative(log)}")
        row = re.compile(r"^\s*(?:(\S+)\s+=>\s+)?(.+?)\s+\(0x[0-9a-fA-F]+\)\s*$")
        loaded = [match.groups() for line in result.stdout.splitlines() if (match := row.fullmatch(line))]
        selected = [path for name, path in loaded if name == "libmbe-neo.so.2"]
        expected = libdir / "libmbe-neo.so.2"
        if len(selected) != 1 or not os.path.samefile(selected[0], expected):
            raise ValueError(f"{role} loader must select exactly {expected}, got {selected}")
        extras = [
            path
            for name, path in loaded
            if name != "libmbe-neo.so.2"
            and (Path(name or path).name.startswith("libmbe") or Path(path).name.startswith("libmbe"))
        ]
        if extras:
            raise ValueError(f"{role} loader selected additional libmbe libraries: {extras}")

    def build_identity(self, directory, role):
        cache = directory / "CMakeCache.txt"
        result = {"cache": None, "cache_sha256": None, "options": None, "compiler": None}
        if not cache.is_file():
            return result
        copied = self.snapshot(cache, f"identity/{role}/CMakeCache.txt")
        options = {}
        for line in copied.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith(("//", "#")) or "=" not in line or ":" not in line:
                continue
            key_type, value = line.split("=", 1)
            key = key_type.split(":", 1)[0]
            if key.startswith(("CMAKE_", "MBELIB_")):
                options[key] = value
        result.update(cache=self.relative(copied), cache_sha256=digest(copied), options=options)
        compiler = options.get("CMAKE_C_COMPILER")
        if compiler and Path(compiler).is_file():
            result["compiler"] = file_identity(compiler)
            result["compiler"]["version"] = self.command([compiler, "--version"], f"{role}-compiler")
        for name in (
            "compile_commands.json",
            "CMakeFiles/mbe_shared.dir/flags.make",
            "CMakeFiles/mbe-shared.dir/flags.make",
        ):
            source = directory / name
            if source.is_file():
                copied = self.snapshot(source, f"identity/{role}/{name}")
                result.setdefault("build_files", {})[self.relative(copied)] = digest(copied)
        return result

    def identities(self):
        paths = {
            "encoder": "build/quality/op25_encode",
            "evaluator": self.args.evaluator,
            "reframer": self.args.reframer,
            "calibrator": "tools/quality/calibrate_encoder.py",
            "runner": "tools/quality/run_quality_ab.sh",
            "runner_python": __file__,
            "quality_support": str(Path(__file__).with_name("quality_support.py")),
        }
        for role, path in paths.items():
            self.manifest["tools"][role] = file_identity(path)
            self.tools[role] = self.manifest["tools"][role]["resolved"]
            self.snapshot(path, f"identity/tools/{role}/{Path(path).name}")
        capture_path = Path("build/quality/checkpoint-capture.json")
        captured = {}
        if capture_path.is_file():
            copied = self.snapshot(capture_path, "identity/checkpoint-capture.json")
            captured = read_json(copied).get("libraries", {})
        for role, requested in (
            ("baseline", self.args.baseline),
            ("candidate", self.args.candidate),
            ("calibration", self.args.calibration_libdir or self.args.baseline),
        ):
            directory = Path(requested).resolve(strict=True)
            self.libdirs[role] = directory
            record = file_identity(directory / "libmbe-neo.so.2")
            record["requested"] = requested
            record["source"] = {"revision": None, "status": "unknown"}
            for point in captured.values():
                if point.get("sha256") == record["sha256"] and point.get("source_revision"):
                    record["source"] = {"revision": point["source_revision"], "status": "captured_exact_revision"}
                    break
            record["build"] = self.build_identity(directory, role)
            self.manifest["libraries"][role] = record
            self.loader_check(directory, role)
        for point in ("original", "checkpoint"):
            capture = captured.get(point, {})
            matching = [
                role for role, record in self.manifest["libraries"].items() if record["sha256"] == capture.get("sha256")
            ]
            self.manifest["comparison_points"][point] = {
                "revision": capture.get("source_revision") if matching else None,
                "status": "captured_exact_revision" if matching else "unknown",
                "sha256": capture.get("sha256") if matching else None,
                "library_roles": matching,
            }
        candidate = self.manifest["libraries"]["candidate"]
        self.manifest["comparison_points"]["candidate"] = dict(candidate["source"], sha256=candidate["sha256"])
        self.manifest["identical_libraries"] = self.manifest["libraries"]["baseline"]["sha256"] == candidate["sha256"]

    def references(self):
        corpus = Path(self.args.corpus).resolve(strict=True)
        source_manifest = corpus / "corpus-manifest.json"
        metadata = None
        copied_manifest = None
        entries = {}
        if source_manifest.is_file():
            copied_manifest = self.snapshot(source_manifest, "references/corpus-manifest.json")
            metadata = read_json(copied_manifest)
            if metadata.get("schema_version") != 1 or metadata.get("set") not in ("legacy", "heldout"):
                raise ValueError("corpus manifest must use schema1 and a known corpus set")
            for entry in metadata["references"]:
                if entry["name"] in entries:
                    raise ValueError("duplicate name in corpus manifest")
                entries[entry["name"]] = entry
            for notice in metadata.get("notices", []):
                source = safe_child(corpus, notice["file"])
                if digest(source) != notice["sha256"]:
                    raise ValueError(f"corpus notice hash mismatch: {source}")
                self.snapshot(source, "references/" + notice["file"])
        raw_files = sorted(corpus.glob("*.raw"))
        if not raw_files:
            raise ValueError(f"no *.raw files in {corpus}")
        if metadata is not None and {path.stem for path in raw_files} != set(entries):
            raise ValueError("corpus files differ from exact corpus manifest reference list")
        names = {path.stem for path in raw_files}
        legacy = names == LEGACY
        self.manifest["corpus"] = {
            "set": metadata["set"] if metadata else "legacy" if legacy else "custom",
            "label": metadata["label"] if metadata else "legacy/development" if legacy else "custom",
            "level": metadata["level"] if metadata else "native" if legacy else "unknown",
            "manifest": self.relative(copied_manifest) if copied_manifest else None,
            "duplicates": [],
            "independent_references": [],
            "requested": self.args.corpus,
            "resolved": str(corpus),
        }
        for source in raw_files:
            name = source.stem
            if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
                raise ValueError(f"unsafe reference name: {name!r}")
            copied = self.snapshot(source, f"references/{name}.raw")
            size = copied.stat().st_size
            if size == 0 or size % 2:
                raise ValueError(f"reference must be nonempty s16le PCM: {source}")
            entry = entries.get(name)
            if entry and (
                entry["file"] != source.name or entry["sha256"] != digest(copied) or entry["samples"] != size // 2
            ):
                raise ValueError(f"corpus reference hash/sample mismatch: {source}")
            self.manifest["references"].append(
                {
                    "name": name,
                    "path": self.relative(copied),
                    "sha256": digest(copied),
                    "samples": size // 2,
                    "speaker": entry.get("speaker") if entry else None,
                    "utterance_index": entry.get("utterance_index") if entry else None,
                    "duplicate_of": None,
                    "source": entry or {"original_path": str(source)},
                }
            )
        by_name = {entry["name"]: entry for entry in self.manifest["references"]}
        if "hts1" in by_name and "hts1a" in by_name:
            full = (self.output / by_name["hts1"]["path"]).read_bytes()
            prefix = (self.output / by_name["hts1a"]["path"]).read_bytes()
            if full.startswith(prefix):
                by_name["hts1a"]["duplicate_of"] = "hts1"
                self.manifest["corpus"]["duplicates"].append(
                    {"name": "hts1a", "duplicate_of": "hts1", "relation": "exact_prefix", "samples": len(prefix) // 2}
                )
        self.manifest["corpus"]["independent_references"] = [
            entry["name"] for entry in self.manifest["references"] if entry["duplicate_of"] is None
        ]
        if self.args.gain_profile:
            copied = self.snapshot(self.args.gain_profile, "gain-profile.json")
            self.profile = load_gain_profile(copied)
            self.manifest["gain_profile"] = {"path": self.relative(copied), "sha256": digest(copied)}
        if self.manifest["corpus"]["set"] == "heldout" and self.profile is None:
            raise ValueError("held-out runs require a frozen --gain-profile; held-out calibration is forbidden")

    def encode(self, reference, mode):
        signature = {
            "raw_sha256": reference["sha256"],
            "encoder_sha256": self.manifest["tools"]["encoder"]["sha256"],
            "evaluator_sha256": self.manifest["tools"]["evaluator"]["sha256"],
            "baseline_sha256": self.manifest["libraries"]["calibration"]["sha256"],
            "calibrator_sha256": self.manifest["tools"]["calibrator"]["sha256"],
            "mode": mode,
            "flush_frames": FLUSH,
            "original_samples": reference["samples"],
        }
        if self.profile is not None:
            # Fixed gains produce identical bits across evaluator/library build variants.
            signature = {
                key: signature[key]
                for key in ("raw_sha256", "encoder_sha256", "mode", "flush_frames", "original_samples")
            }
            signature["gain_profile_sha256"] = self.manifest["gain_profile"]["sha256"]
            signature["gain_adjust_log2"] = self.profile["gain_adjust_log2"][mode]
        key = signature_hash(signature)
        cache = Path("build/quality/frame-cache-v3") / key
        cache.mkdir(parents=True, exist_ok=True)
        frames, calibration = cache / "frames.txt", cache / "frames.calibration.json"
        width = {"imbe7200": 88, "ambe2450": 49, "ambe2400": 96}[mode]
        expected_count = (reference["samples"] + 159) // 160 + FLUSH
        with (cache / ".lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            cached = None
            try:
                cached = read_json(calibration)
                if (
                    cached.get("schema_version") != 2
                    or cached.get("inputs") != signature
                    or cached.get("frames_sha256") != digest(frames)
                    or cached.get("original_samples") != reference["samples"]
                    or cached.get("flush_frames") != FLUSH
                ):
                    cached = None
                if cached is not None and frame_count(frames, width) != expected_count:
                    cached = None
            except (OSError, ValueError, TypeError, AttributeError):
                cached = None
            if cached is None:
                # Generate separately, then atomically publish a complete signature-bound pair.
                with tempfile.TemporaryDirectory(prefix=".encode-", dir=cache) as temporary:
                    bits = Path(temporary) / "frames.txt"
                    raw = self.output / reference["path"]
                    if self.profile is None:
                        self.command(
                            [
                                sys.executable,
                                self.tools["calibrator"],
                                self.tools["encoder"],
                                self.tools["evaluator"],
                                self.libdirs["calibration"],
                                raw,
                                mode,
                                bits,
                                FLUSH,
                            ],
                            f"{reference['name']}-{mode}-calibrate",
                        )
                    else:
                        gain = self.profile["gain_adjust_log2"][mode]
                        self.command(
                            [
                                self.tools["encoder"],
                                "--mode",
                                mode,
                                "--in",
                                raw,
                                "--out",
                                bits,
                                "--gain-adjust",
                                repr(gain),
                                "--flush-frames",
                                FLUSH,
                            ],
                            f"{reference['name']}-{mode}-encode",
                        )
                        write_json(
                            bits.with_suffix(".calibration.json"),
                            {
                                "schema_version": 2,
                                "calibration_method": "fixed_profile",
                                "gain_profile_sha256": self.manifest["gain_profile"]["sha256"],
                                "gain_adjust_log2": gain,
                                "original_samples": reference["samples"],
                                "flush_frames": FLUSH,
                                "frames_sha256": digest(bits),
                                "inputs": signature,
                            },
                        )
                    if frame_count(bits, width) != expected_count:
                        raise ValueError("encoder did not emit every input block plus exactly five flush frames")
                    os.replace(bits, frames)
                    os.replace(bits.with_suffix(".calibration.json"), calibration)
            copied_frames = self.snapshot(frames, f"encoded/{reference['name']}.{mode}.txt")
            copied_calibration = self.snapshot(calibration, f"calibration/{reference['name']}.{mode}.calibration.json")
        metadata = read_json(copied_calibration)
        gain = metadata.get("gain_adjust_log2")
        if type(gain) not in (int, float) or not math.isfinite(gain) or not -16 <= gain <= 16:
            raise ValueError(f"invalid gain in {copied_calibration}")
        if metadata.get("inputs") != signature or metadata.get("frames_sha256") != digest(copied_frames):
            raise ValueError(f"calibration signature/content mismatch: {copied_calibration}")
        if self.profile is not None and (
            metadata.get("calibration_method") != "fixed_profile"
            or metadata.get("gain_profile_sha256") != self.manifest["gain_profile"]["sha256"]
            or gain != self.profile["gain_adjust_log2"][mode]
        ):
            raise ValueError(f"fixed gain profile was not preserved: {copied_calibration}")
        record = {
            "name": reference["name"],
            "mode": mode,
            "path": self.relative(copied_frames),
            "sha256": digest(copied_frames),
            "frame_count": frame_count(copied_frames, width),
            "calibration": self.relative(copied_calibration),
            "gain_adjust_log2": gain,
            "cache_signature": signature,
        }
        self.manifest["encodings"].append(record)
        return record

    def evaluate(self, reference, mode, frames, variant, lag=None, canonical=False):
        label = f"{reference['name']}.{'imbe-data' if canonical else mode}.{variant}"
        wav, report = self.output / f"audio/{label}.wav", self.output / f"metrics/{label}.json"
        wav.parent.mkdir(parents=True, exist_ok=True)
        report.parent.mkdir(parents=True, exist_ok=True)
        argv = [
            self.tools["evaluator"],
            "--codec",
            mode,
            "--frames",
            frames,
            "--out",
            wav,
            "--ref",
            self.output / reference["path"],
            "--json",
            report,
            "--seed",
            SEED,
        ]
        if lag is not None:
            argv.extend(["--lag-samples", lag])
        errors = []
        try:
            self.command(argv, label, self.libdirs[variant])
        except RuntimeError as exc:
            errors.append(str(exc))
        measured = None
        if report.is_file():
            try:
                measured = read_json(report)
                api = "mbe_processImbe4400Dataf" if canonical else APIS[mode]
                errors.extend(measurement_errors(measured, reference["samples"], api))
            except (ValueError, TypeError, KeyError) as exc:
                errors.append(f"invalid metric JSON: {exc}")
        else:
            errors.append("evaluator produced no metric JSON")
        if not wav.is_file():
            errors.append("evaluator produced no WAV")
        elif isinstance(measured, dict):
            try:
                parameters, _ = pcm_identity(wav)
                if parameters[:3] != (1, 2, 8000) or parameters[4] != "NONE":
                    errors.append("decoded WAV is not mono 8-kHz signed 16-bit PCM")
                if parameters[3] != measured.get("decoded_samples"):
                    errors.append("decoded WAV sample count differs from evaluator report")
            except (OSError, ValueError, wave.Error) as exc:
                errors.append(f"invalid decoded WAV: {exc}")
        return (
            {
                "wav": self.relative(wav) if wav.is_file() else None,
                "json": self.relative(report) if report.is_file() else None,
            },
            measured,
            errors,
        )

    def compare(self, reference, mode, encoding):
        name = reference["name"]
        entry = {
            "name": name,
            "mode": mode,
            "reference": reference["path"],
            "frames": None,
            "frame_sha256": None,
            "frame_count": None,
            "baseline": {"wav": None, "json": None},
            "candidate": {"wav": None, "json": None},
            "candidate_auto_sensitivity": None,
            "valid": False,
            "invalid_reasons": [],
        }
        self.manifest["reports"].append(entry)
        if encoding is None:
            entry["invalid_reasons"].append("canonical encoding failed; see run failures and command logs")
            return entry
        frames = self.output / f"frames/{name}.{mode}.txt"
        frames.parent.mkdir(parents=True, exist_ok=True)
        try:
            if mode == "ambe2400":
                self.snapshot(self.output / encoding["path"], self.relative(frames))
            else:
                self.command(
                    [self.tools["reframer"], "--codec", mode, "--in", self.output / encoding["path"], "--out", frames],
                    f"{name}-{mode}-reframe",
                )
            count = frame_count(frames, {"imbe7200": 184, "imbe7100": 168, "ambe2450": 96, "ambe2400": 96}[mode])
            if count != encoding["frame_count"]:
                raise ValueError("reframing changed frame count")
            entry.update(frames=self.relative(frames), frame_sha256=digest(frames), frame_count=count)
            baseline, before, errors = self.evaluate(reference, mode, frames, "baseline")
            entry["baseline"] = baseline
            entry["invalid_reasons"].extend("baseline: " + error for error in errors)
            lag = before.get("lag_samples") if isinstance(before, dict) else None
            if type(lag) is not int or not -160 <= lag <= 800:
                entry["invalid_reasons"].append(
                    "baseline lag unavailable; candidate common-lag comparison not possible"
                )
                return entry
            candidate, after, errors = self.evaluate(reference, mode, frames, "candidate", lag)
            entry["candidate"] = candidate
            entry["invalid_reasons"].extend("candidate: " + error for error in errors)
            for variant, measured in (("baseline", before), ("candidate", after)):
                if not isinstance(measured, dict) or measured.get("frames") != count:
                    entry["invalid_reasons"].append(f"{variant}: decoded frame count differs from input")
            if isinstance(after, dict) and after.get("lag_samples") != lag:
                entry["invalid_reasons"].append("candidate did not use baseline lag")
            if (
                isinstance(after, dict)
                and before.get("auto_lag_samples") != after.get("auto_lag_samples")
                and candidate["wav"]
            ):
                sensitivity = self.output / f"metrics/{name}.{mode}.candidate.auto-alignment-sensitivity.json"
                self.command(
                    [
                        self.tools["evaluator"],
                        "--decoded",
                        self.output / candidate["wav"],
                        "--ref",
                        self.output / reference["path"],
                        "--json",
                        sensitivity,
                    ],
                    f"{name}-{mode}-auto-sensitivity",
                    self.libdirs["candidate"],
                )
                entry["candidate_auto_sensitivity"] = {
                    "kind": "automatic_alignment_sensitivity",
                    "json": self.relative(sensitivity),
                }
            if self.manifest["identical_libraries"]:
                identity_errors = self.manifest["identity_check"]["failures"]
                if before != after:
                    identity_errors.append(
                        f"{name}/{mode}: complete metric JSON differs (including null/support fields)"
                    )
                if not baseline["wav"] or not candidate["wav"]:
                    identity_errors.append(f"{name}/{mode}: missing WAV for identity assertion")
                elif digest(self.output / baseline["wav"]) != digest(self.output / candidate["wav"]) or pcm_identity(
                    self.output / baseline["wav"]
                ) != pcm_identity(self.output / candidate["wav"]):
                    identity_errors.append(f"{name}/{mode}: WAV/PCM identity mismatch")
            entry["valid"] = not entry["invalid_reasons"]
        except (OSError, ValueError, RuntimeError, wave.Error) as exc:
            entry["invalid_reasons"].append(str(exc))
        return entry

    def equivalence(self, reference, encoding, reports):
        for variant in ("baseline", "candidate"):
            result = {
                "name": reference["name"],
                "variant": variant,
                "canonical_wav": None,
                "canonical_json": None,
                "imbe7200_equal": False,
                "imbe7100_equal": False,
                "frames_equal": False,
                "invalid_reasons": [],
            }
            self.manifest["frame_equivalence"].append(result)
            try:
                if encoding is None:
                    raise ValueError("canonical IMBE encoding unavailable")
                baseline_path = reports["imbe7200"]["baseline"]["json"]
                if baseline_path is None:
                    raise ValueError("IMBE7200 baseline lag unavailable")
                lag = read_json(self.output / baseline_path)["lag_samples"]
                paths, measured, errors = self.evaluate(
                    reference, "imbe7200", self.output / encoding["path"], variant, lag, canonical=True
                )
                result.update(canonical_wav=paths["wav"], canonical_json=paths["json"])
                result["invalid_reasons"].extend(errors)
                counts = [encoding["frame_count"], measured.get("frames") if isinstance(measured, dict) else None]
                for mode in ("imbe7200", "imbe7100"):
                    comparison = reports[mode][variant]
                    if not paths["wav"] or not comparison["wav"] or not comparison["json"]:
                        raise ValueError(f"missing {mode} output for frame equality")
                    result[f"{mode}_equal"] = pcm_identity(self.output / paths["wav"]) == pcm_identity(
                        self.output / comparison["wav"]
                    )
                    counts.extend([reports[mode]["frame_count"], read_json(self.output / comparison["json"])["frames"]])
                result["frames_equal"] = None not in counts and len(set(counts)) == 1
            except (OSError, ValueError, RuntimeError, TypeError, KeyError, wave.Error) as exc:
                result["invalid_reasons"].append(str(exc))

    def run(self):
        self.identities()
        self.references()
        for reference in self.manifest["references"]:
            encoded = {}
            for mode in ENCODED_MODES:
                try:
                    encoded[mode] = self.encode(reference, mode)
                except (OSError, ValueError, RuntimeError) as exc:
                    self.manifest["failures"].append(f"{reference['name']}/{mode}: {exc}")
                    encoded[mode] = None
            reports = {
                mode: self.compare(reference, mode, encoded["imbe7200" if mode == "imbe7100" else mode])
                for mode in MODES
            }
            self.equivalence(reference, encoded["imbe7200"], reports)
        try:
            self.manifest["correctness_tests"] = run_correctness_tests(self.libdirs["candidate"], self.output)
        except (OSError, ValueError, RuntimeError) as exc:
            self.manifest["correctness_tests"] = {"status": "fail", "gates": {}, "tests": [], "error": str(exc)}

    def finalize(self):
        manifest = self.manifest
        for category in ("libraries", "tools"):
            for role, record in manifest[category].items():
                try:
                    if digest(record["resolved"]) != record["sha256"]:
                        manifest["failures"].append(f"{category}/{role} changed during the run")
                except OSError as exc:
                    manifest["failures"].append(f"{category}/{role} identity unavailable at run end: {exc}")
        expected = len(manifest["references"]) * len(MODES)
        reports_complete = expected > 0 and len(manifest["reports"]) == expected
        valid = reports_complete and all(report["valid"] for report in manifest["reports"])
        equivalence = (
            len(manifest["frame_equivalence"]) == len(manifest["references"]) * 2
            and bool(manifest["frame_equivalence"])
            and all(
                item["imbe7200_equal"]
                and item["imbe7100_equal"]
                and item["frames_equal"]
                and not item["invalid_reasons"]
                for item in manifest["frame_equivalence"]
            )
        )
        identity = manifest["identity_check"]
        if manifest["identical_libraries"]:
            if not valid:
                identity["failures"].append("identity control has incomplete or invalid frame-mode measurements")
            identity["status"] = "fail" if identity["failures"] else "pass"
        failed = bool(manifest["failures"]) or not valid or not equivalence or identity["status"] == "fail"
        tests = manifest["correctness_tests"].get("status", "not_established")
        failed = failed or tests == "fail"
        status = "fail" if failed else "pass" if tests == "pass" else "not_established"
        acceptance = {
            "schema_version": 1,
            "correctness": status,
            "objective_nonregression": "not_established",
            "perceptual": "not_established",
            "evidence": {
                "supported_full_reference_measurements": "pass" if valid else "fail",
                "frame_equivalence": "pass" if equivalence else "fail",
                "same_library_identity": identity["status"],
                "source_correctness_tests": tests,
            },
            "reasons": list(manifest["failures"]),
            "assessment": "Runner correctness evidence only; objective and perceptual assessment requires report.",
        }
        for report in manifest["reports"]:
            acceptance["reasons"].extend(
                f"{report['name']}/{report['mode']}: {reason}" for reason in report["invalid_reasons"]
            )
        acceptance["reasons"].extend(identity["failures"])
        for item in manifest["frame_equivalence"]:
            if not (item["imbe7200_equal"] and item["imbe7100_equal"] and item["frames_equal"]):
                acceptance["reasons"].append(f"{item['name']}/{item['variant']}: canonical/frame equality failed")
            acceptance["reasons"].extend(item["invalid_reasons"])
        rows = []
        groups = {mode: [] for mode in MODES}
        for report in manifest["reports"]:
            row = {
                "name": report["name"],
                "mode": report["mode"],
                "valid": report["valid"],
                "independent": report["name"] in manifest["corpus"].get("independent_references", []),
                "baseline": None,
                "candidate": None,
                "invalid_reasons": report["invalid_reasons"],
            }
            for variant in ("baseline", "candidate"):
                path = report[variant]["json"]
                if path:
                    try:
                        row[variant] = read_json(self.output / path)
                    except (OSError, ValueError):
                        pass
            rows.append(row)
            if row["valid"] and row["independent"]:
                groups[row["mode"]].append(row)
        means = {}
        for mode, values in groups.items():
            means[mode] = {"independent_valid_clips": len(values), "baseline": {}, "candidate": {}}
            for variant in ("baseline", "candidate"):
                for metric in SUMMARY_METRICS:
                    samples = [row[variant].get(metric) for row in values]
                    means[mode][variant][metric] = (
                        sum(samples) / len(samples)
                        if samples and all(type(value) in (int, float) and math.isfinite(value) for value in samples)
                        else None
                    )
        write_json(
            self.output / "table.json",
            {
                "schema_version": 1,
                "rows": rows,
                "means": means,
                "aggregation": "Unweighted valid independent clips; duplicate prefixes excluded.",
            },
        )
        write_json(self.output / "commands.json", self.commands)
        summary = [
            f"Run directory: {self.output}",
            f"Correctness: {status}",
            "Objective nonregression: not_established; perceptual: not_established",
            "Per-clip baseline/candidate LSD and envelope correlation:",
        ]

        def cell(row, metric):
            return "/".join(
                "null"
                if not isinstance(row[variant], dict) or type(row[variant].get(metric)) not in (int, float)
                else f"{row[variant][metric]:.4f}"
                for variant in ("baseline", "candidate")
            )

        for row in rows:
            summary.append(
                f"{row['name']} {row['mode']}: LSD={cell(row, 'lsd_db')} env={cell(row, 'env_corr')} "
                f"{'valid' if row['valid'] else 'INVALID'}"
                f"{' [duplicate; excluded from independent means]' if not row['independent'] else ''}"
            )
        summary.append("Independent valid-clip means (incomplete groups are not acceptance evidence):")
        for mode, group in means.items():
            summary.append(
                f"{mode}: n={group['independent_valid_clips']} "
                f"LSD={group['baseline']['lsd_db']}/{group['candidate']['lsd_db']} "
                f"env={group['baseline']['env_corr']}/{group['candidate']['env_corr']}"
            )
        text = "\n".join(summary) + "\n"
        (self.output / "summary.txt").write_text(text, encoding="utf-8")
        manifest["files"] = {
            self.relative(path): digest(path)
            for path in sorted(self.output.rglob("*"))
            if path.is_file() and self.relative(path) not in ("manifest.json", "acceptance.json")
        }
        write_json(self.output / "manifest.json", manifest)
        write_json(self.output / "acceptance.json", acceptance)
        print(text, end="")
        return 1 if failed else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", help="baseline directory containing libmbe-neo.so.2")
    parser.add_argument("candidate", help="candidate directory containing libmbe-neo.so.2")
    parser.add_argument("corpus", nargs="?", default="build/quality/corpus")
    parser.add_argument("--out-dir")
    parser.add_argument("--evaluator", default="build/dev-debug/mbe_quality_eval")
    parser.add_argument("--reframer", default="build/dev-debug/mbe_quality_reframe")
    parser.add_argument("--calibration-libdir")
    parser.add_argument("--gain-profile")
    args = parser.parse_args()
    if platform.system() != "Linux":
        parser.error("this runner requires Linux (LD_LIBRARY_PATH)")
    if os.environ.get("LD_PRELOAD"):
        parser.error("LD_PRELOAD must be unset for controlled library selection")
    for directory in (args.baseline, args.candidate, args.calibration_libdir or args.baseline):
        path = Path(directory).resolve()
        if ":" in str(path):
            parser.error("library directories must not contain ':' (LD_LIBRARY_PATH separator)")
        if not (path / "libmbe-neo.so.2").is_file() or not os.access(path / "libmbe-neo.so.2", os.R_OK):
            parser.error(f"expected readable libmbe-neo.so.2 in {directory}")
    for executable in ("build/quality/op25_encode", args.evaluator, args.reframer):
        if not Path(executable).is_file() or not os.access(executable, os.X_OK):
            parser.error(f"missing executable prerequisite: {executable}")
    if not Path(args.corpus).is_dir():
        parser.error(f"missing corpus directory: {args.corpus}")
    if args.gain_profile:
        load_gain_profile(args.gain_profile)
    if args.out_dir:
        output = Path(args.out_dir).resolve()
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            parser.error("--out-dir must not exist or must be an empty directory; refusing overwrite")
        output.mkdir(parents=True, exist_ok=True)
    else:
        parent = Path("build/quality/runs")
        parent.mkdir(parents=True, exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix="paired-", dir=parent)).resolve()
    print(f"Run directory: {output}", file=sys.stderr)
    runner = Runner(args, output)
    try:
        runner.run()
    except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
        runner.manifest["failures"].append(str(exc))
        print(f"quality runner: {exc}", file=sys.stderr)
    return runner.finalize()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as exc:
        print(f"quality runner: {exc}", file=sys.stderr)
        sys.exit(2)
