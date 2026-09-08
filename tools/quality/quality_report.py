#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Opt-in, immutable-run spectral and intelligibility assessment.

C evaluator snapshots remain authoritative. All comparisons use the baseline's
recorded alignment and full reference support; none optimizes regenerated phase.
"""

import csv
import json
import math
import platform
import statistics
import warnings
from collections import defaultdict
from pathlib import Path

import quality_audio as audio
import quality_support as support

VARIANTS = ("baseline", "candidate")
AXES = ("correctness", "objective_nonregression", "perceptual")
METRICS = ("lsd_db", "env_corr", "crest_delta_db", "join_excess_db") + support.BAND_METRICS
LIMITATIONS = [
    "STOI/ESTOI are intelligibility diagnostics, not naturalness or MOS estimates.",
    "Boundary index and synthesis-join statistics are diagnostics, not standardized click detectors.",
    "IMBE7100 uses reframed canonical IMBE7200 parameters, not an independent native transmitter capture.",
    "Patched OP25 AMBE2400 is a decoder-targeted fixture, not universal wire-conformance evidence.",
    "Shared IMBE7100/7200 PCM is not an additional independent speech mode or listening sample.",
    "Failed screens are completed negative assessments; held-out evidence must not be used to retune DSP or gains.",
]


def _status(values):
    values = list(values)
    return "fail" if "fail" in values else "pass" if values and all(v == "pass" for v in values) else "not_established"


def _gate(status, reason, **details):
    return {"status": status, "reason": reason, **details}


def _output_path(root, manifest, relative):
    """Never let a derived artifact overwrite an immutable input or escape a run."""
    path = root / relative
    if relative in manifest["files"] or relative == "manifest.json":
        raise ValueError(f"derived output collides with immutable run input: {relative}")
    resolved = path.resolve()
    if not resolved.is_relative_to(root) or any(resolved == (root / name).resolve() for name in manifest["files"]):
        raise ValueError(f"unsafe derived output path: {relative}")
    if path.exists() and path.stat().st_nlink > 1:
        raise ValueError(f"derived output has hard-link aliases: {relative}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _coverage(manifest):
    references = manifest["references"]
    names = {row["name"] for row in references}
    corpus = manifest["corpus"]
    errors = []
    if corpus.get("set") == "legacy" and names != set(support.LEGACY_NAMES):
        errors.append("legacy corpus does not contain the complete five-reference development set")
    if corpus.get("set") == "heldout":
        actual = [(row.get("speaker"), row.get("utterance_index")) for row in references]
        expected = {(speaker, index) for speaker in ("bdl", "rms", "clb", "slt") for index in range(1, 11)}
        if len(actual) != 40 or set(actual) != expected:
            errors.append("held-out corpus requires exactly four speakers times utterances 1 through 10")
        if corpus.get("level") not in ("native", "minus12db"):
            errors.append("held-out corpus level is not identified as native or minus12db")
    return errors


def _independent_names(root, manifest):
    names = {row["name"] for row in manifest["references"]}
    references = {row["name"]: row for row in manifest["references"]}
    excluded = []
    # Do not allow an arbitrary manifest exclusion to hide an unfavorable clip.
    if {"hts1", "hts1a"} <= names:
        child = support.artifact_path(root, references["hts1a"]["path"]).read_bytes()
        with support.artifact_path(root, references["hts1"]["path"]).open("rb") as stream:
            if child == stream.read(len(child)):
                names.remove("hts1a")
                excluded.append(
                    {
                        "name": "hts1a",
                        "duplicate_of": "hts1",
                        "relation": "exact_prefix",
                        "samples": references["hts1a"]["samples"],
                    }
                )
    if set(manifest["corpus"].get("independent_references", names)) != names:
        raise ValueError("independent-reference declaration differs from verified exact hts1a duplication")
    return names, excluded


def _rail(record):
    if record is None:
        return None
    count, total = record.get("pcm_rail_samples"), record.get("decoded_samples")
    return {
        "samples": count,
        "total_samples": total,
        "rate": count / total
        if isinstance(count, (int, float)) and isinstance(total, (int, float)) and total > 0
        else None,
        "max_run": record.get("pcm_max_rail_run"),
        "peak": record.get("pcm_peak"),
        "float_clip_samples": record.get("float_clip_samples"),
        "float_observation": "post-library-limiter, pre-int16",
    }


def _screens(baseline, candidate, mean_abs_lsd, crest_baseline, crest_candidate):
    """The same threshold definitions serve aggregate gates and clip diagnostics."""
    checks = {
        "lsd": (mean_abs_lsd, 0.3, "mean absolute paired LSD change (dB)"),
        "envelope": (baseline["env_corr"] - candidate["env_corr"], 0.01, "envelope correlation decrease"),
        "join": (
            candidate["join_excess_db"] - max(0, baseline["join_excess_db"]),
            0.1,
            "join excess above max(0, baseline) (dB; diagnostic)",
        ),
        "crest_mean_absolute_error": (
            crest_candidate - crest_baseline,
            0.0,
            "mean per-clip absolute crest-error increase (dB)",
        ),
    }
    checks.update(
        {
            key: (abs(candidate[key] - baseline[key]), 0.5, "absolute paired mean band change (dB)")
            for key in support.BAND_METRICS
        }
    )
    return {
        key: {
            "status": "pass" if value <= limit else "fail",
            "value": value,
            "limit": limit,
            "description": description,
        }
        for key, (value, limit, description) in checks.items()
    }


def _clip_exceedances(clip):
    b, c = clip["baseline"], clip["candidate"]
    checks = _screens(b, c, abs(c["lsd_db"] - b["lsd_db"]), abs(b["crest_delta_db"]), abs(c["crest_delta_db"]))
    clip["screens"] = checks
    return [
        f"{key}: {row['value']:.9g} > {row['limit']:.9g} ({row['description']})"
        for key, row in checks.items()
        if row["status"] == "fail"
    ]


def _load_clips(root, manifest, independent):
    import numpy as np
    from pystoi import stoi

    references = {row["name"]: row for row in manifest["references"]}
    pcm_references = {
        name: audio.read_pcm(support.artifact_path(root, row["path"])) for name, row in references.items()
    }
    clips, crosschecks = [], []
    for row in manifest["reports"]:
        ref = references[row["name"]]
        clip = {
            "name": row["name"],
            "mode": row["mode"],
            "speaker": ref.get("speaker"),
            "utterance_index": ref.get("utterance_index"),
            "level": manifest["corpus"].get("level", "unknown"),
            "reference": row["reference"],
            "included_in_aggregate": row["name"] in independent,
            "threshold_exceedances": [],
            "invalid_reasons": list(row.get("invalid_reasons", [])),
            "stoi": {v: None for v in VARIANTS},
            "estoi": {v: None for v in VARIANTS},
            "intelligibility_warnings": [],
            "pcm_equal": None,
            "plots": [],
            "candidate_auto_sensitivity": row.get("candidate_auto_sensitivity"),
        }
        clips.append(clip)
        if row["reference"] != ref["path"]:
            raise ValueError(f"report reference differs from named reference: {row['name']}/{row['mode']}")
        decoded = {}
        for variant in VARIANTS:
            record = row[variant]
            clip[f"{variant}_wav"] = record.get("wav")
            clip[variant] = (
                support.load_json(support.artifact_path(root, record["json"])) if record.get("json") else None
            )
            errors = support.measurement_errors(clip[variant], ref["samples"], support.INPUT_APIS[row["mode"]])
            clip["invalid_reasons"].extend(f"{variant}: {error}" for error in errors)
            if record.get("wav"):
                decoded[variant] = audio.read_pcm(support.artifact_path(root, record["wav"]))
            else:
                clip["invalid_reasons"].append(f"{variant}: missing decoded WAV")
            if clip[variant] is not None and clip[variant].get("frames") != row.get("frame_count"):
                clip["invalid_reasons"].append(f"{variant}: decoded frame count differs from frame snapshot")
        if not row["valid"] and not clip["invalid_reasons"]:
            clip["invalid_reasons"].append("runner recorded an invalid measurement")
        b, c = clip["baseline"], clip["candidate"]
        common_keys = (
            "lag_samples",
            "reference_trim_head",
            "reference_trim_tail",
            "decoded_trim_head",
            "aligned_samples",
            "active_frames",
            "spectral_frames",
            "join_count",
        )
        if b and c and any(b.get(key) != c.get(key) for key in common_keys):
            clip["invalid_reasons"].append(
                "candidate does not use baseline common lag/reference activity/spectral/join support"
            )
        if b and b.get("lag_samples") != b.get("auto_lag_samples"):
            clip["invalid_reasons"].append("baseline primary lag is not its automatic estimate")
        if set(decoded) == set(VARIANTS):
            clip["pcm_equal"] = bool(np.array_equal(decoded["baseline"], decoded["candidate"]))
        clip["rails"] = {v: _rail(clip[v]) for v in VARIANTS}
        increases = []
        if all(clip["rails"][v] is not None for v in VARIANTS):
            for key in ("samples", "rate", "max_run", "float_clip_samples"):
                a, z = (clip["rails"][v][key] for v in VARIANTS)
                if a is not None and z is not None and z > a:
                    increases.append({"metric": key, "baseline": a, "candidate": z, "increase": z - a})
        clip["rail_increases"] = increases
        clip["valid"] = not clip["invalid_reasons"]
        if not clip["valid"]:
            continue
        clip["threshold_exceedances"] = _clip_exceedances(clip)
        reference = pcm_references[row["name"]]
        for variant in VARIANTS:
            r, d = audio.aligned_signals(reference, decoded[variant], b)
            for key, extended in (("stoi", False), ("estoi", True)):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    score = float(stoi(r, d, 8000, extended=extended))
                messages = [str(w.message) for w in caught]
                clip["intelligibility_warnings"].extend(f"{variant}/{key}: {msg}" for msg in messages)
                if not math.isfinite(score) or any("not enough" in msg.lower() for msg in messages):
                    clip["invalid_reasons"].append(
                        f"{variant}/{key}: unsupported full-overlap intelligibility measurement"
                    )
                else:
                    clip[key][variant] = score
            # One paired real clip per actual mode also exposes differing FFT paths.
            if not any(check["mode"] == row["mode"] and check["variant"] == variant for check in crosschecks):
                independent_lsd = audio.numpy_lsd(reference, decoded[variant], clip[variant])
                if independent_lsd is not None and not math.isfinite(independent_lsd):
                    independent_lsd = None
                difference = abs(independent_lsd - clip[variant]["lsd_db"]) if independent_lsd is not None else None
                crosschecks.append(
                    {
                        "name": row["name"],
                        "mode": row["mode"],
                        "variant": variant,
                        "numpy_lsd_db": independent_lsd,
                        "evaluator_lsd_db": clip[variant]["lsd_db"],
                        "absolute_difference_db": difference,
                        "tolerance_db": 1e-6,
                        "status": "pass"
                        if difference is not None and math.isfinite(difference) and difference <= 1e-6
                        else "fail",
                    }
                )
        clip["valid"] = not clip["invalid_reasons"]
    return clips, crosschecks


def _aggregate(clips):
    included = [clip for clip in clips if clip["included_in_aggregate"]]
    valid = [clip for clip in included if clip["valid"]]
    result = {
        "clips": len(included),
        "valid_clips": len(valid),
        "excluded_duplicates": len(clips) - len(included),
        "means": None,
        "stoi": {v: None for v in VARIANTS},
        "estoi": {v: None for v in VARIANTS},
        "screens": {},
        "worst_clips": {},
        "objective_nonregression": "not_established",
        "reasons": [],
    }
    if valid:
        means = {
            variant: {metric: statistics.mean(clip[variant][metric] for clip in valid) for metric in METRICS}
            for variant in VARIANTS
        }
        result["means"] = means
        crest_abs = {v: statistics.mean(abs(clip[v]["crest_delta_db"]) for clip in valid) for v in VARIANTS}
        result["crest_mean_perclip_absolute_error_db"] = crest_abs
        result["crest_absolute_signed_mean_error_db"] = {v: abs(means[v]["crest_delta_db"]) for v in VARIANTS}
        result["crest_absolute_signed_mean_decreased"] = abs(means["candidate"]["crest_delta_db"]) < abs(
            means["baseline"]["crest_delta_db"]
        )
        delta = statistics.mean(abs(clip["candidate"]["lsd_db"] - clip["baseline"]["lsd_db"]) for clip in valid)
        result["mean_absolute_paired_lsd_change_db"] = delta
        result["screens"] = _screens(
            means["baseline"], means["candidate"], delta, crest_abs["baseline"], crest_abs["candidate"]
        )
        for key in ("stoi", "estoi"):
            result[key] = {v: statistics.mean(clip[key][v] for clip in valid) for v in VARIANTS}
            result[key]["candidate_minus_baseline"] = result[key]["candidate"] - result[key]["baseline"]
        for key in result["screens"]:
            worst = max(valid, key=lambda clip: clip["screens"][key]["value"])
            result["worst_clips"][key] = {
                "name": worst["name"],
                "value": worst["screens"][key]["value"],
                "limit": worst["screens"][key]["limit"],
            }
        result["objective_nonregression"] = _status(check["status"] for check in result["screens"].values())
        result["reasons"].extend(
            f"{key}: {check['value']:.9g} > {check['limit']:.9g}"
            for key, check in result["screens"].items()
            if check["status"] == "fail"
        )
    if len(valid) != len(included) or any(not clip["valid"] for clip in clips):
        result["objective_nonregression"] = "fail"
        result["reasons"].append("incomplete or unsupported measurements cannot establish objective non-regression")
    if not included:
        result["reasons"].append("no independent reference clips")
    result["rails"] = {}
    for variant in VARIANTS:
        observations = [clip["rails"][variant] for clip in included if clip["rails"][variant] is not None]
        complete = (
            len(observations) == len(included)
            and bool(observations)
            and all(
                all(row[key] is not None for key in ("samples", "total_samples", "max_run")) for row in observations
            )
        )
        count = sum(row["samples"] for row in observations) if complete else None
        total = sum(row["total_samples"] for row in observations) if complete else None
        result["rails"][variant] = {
            "samples": count,
            "total_samples": total,
            "rate": count / total if total else None,
            "max_run": max(row["max_run"] for row in observations) if complete else None,
        }
    result["rail_increases"] = [
        {"name": clip["name"], "increases": clip["rail_increases"]} for clip in clips if clip["rail_increases"]
    ]
    return result


def _native_gates(root, manifest):
    evidence = manifest.get("correctness_tests", {})
    tests = {row.get("name"): row for row in evidence.get("tests", [])}
    statuses = {}
    for name in ("test_params", "test_frame_paths", "test_noise_determinism"):
        row = tests.get(name, {})
        claimed = row.get("status", "not_established")
        status = "fail" if claimed == "fail" else "not_established"
        logs = all(
            isinstance(row.get(key), str) and row[key] in manifest["files"]
            for key in ("stdout", "stderr", "loader_trace")
        )
        hashes = all(
            isinstance(row.get(key), str) and len(row[key]) == 64
            for key in ("executable_sha256", "source_sha256", "library_sha256")
        )
        binding = (
            row.get("binding") == "matching-build-static-archive"
            if name == "test_params"
            else row.get("binding") == "loader-samefile"
            and row.get("library_sha256") == manifest["libraries"]["candidate"]["sha256"]
        )
        if claimed == "pass" and row.get("returncode") == 0 and logs and hashes and binding:
            status = "pass"
        statuses[name] = status
    requirements = {
        "nonnegative_attenuation": ("test_params",),
        "tone_dispatch": ("test_frame_paths",),
        "warm_fft_wola": ("test_params", "test_noise_determinism"),
    }
    result = {}
    for gate, names in requirements.items():
        claimed = evidence.get("gates", {}).get(gate, "not_established")
        status = _status([claimed] + [statuses[name] for name in names])
        result[gate] = _gate(
            status, "recorded candidate-bound native regressions; no tests executed by report", tests=list(names)
        )
    return result


def _correctness(root, manifest, clips, crosschecks, coverage_errors):
    gates = _native_gates(root, manifest)
    gates["corpus_coverage"] = _gate(
        "fail" if coverage_errors else "pass",
        "; ".join(coverage_errors) or "all declared references cover four actual frame modes",
    )
    expected = {(row["name"], variant) for row in manifest["references"] for variant in VARIANTS}
    actual = {(row["name"], row["variant"]) for row in manifest["frame_equivalence"]}
    equivalence_ok = (
        actual == expected
        and len(manifest["frame_equivalence"]) == len(expected)
        and all(
            row.get("imbe7200_equal") is True and row.get("imbe7100_equal") is True and row.get("frames_equal") is True
            for row in manifest["frame_equivalence"]
        )
    )
    gates["frame_equivalence"] = _gate(
        "pass" if equivalence_ok else "fail",
        "canonical Dataf and both IMBE Framef PCM/count equality, checked from immutable snapshots",
        expected_controls=len(expected),
        recorded_controls=len(manifest["frame_equivalence"]),
        missing_controls=[{"name": name, "variant": variant} for name, variant in sorted(expected - actual)],
    )
    modes = {}
    checkpoint = manifest["libraries"]["baseline"].get("source", {}).get("revision") == support.CHECKPOINT_REVISION
    for mode in support.MODES:
        selected = [clip for clip in clips if clip["mode"] == mode]
        local = dict(gates)
        errors = [{"name": clip["name"], "reasons": clip["invalid_reasons"]} for clip in selected if not clip["valid"]]
        local["supported_measurements"] = _gate(
            "fail" if errors else "pass",
            "finite supported speech, full reference coverage after flushing, common lag and nonfinite float count zero",
            invalid_clips=errors,
        )
        checks = [row for row in crosschecks if row["mode"] == mode]
        local["independent_numpy_lsd"] = _gate(
            _status(row["status"] for row in checks),
            "paired real-clip NumPy versus authoritative evaluator LSD within 1e-6 dB",
            comparisons=checks,
        )
        if manifest["identical_libraries"]:
            failures = [
                clip["name"]
                for clip in selected
                if clip["pcm_equal"] is not True or clip["baseline"] != clip["candidate"]
            ]
            local["same_library_identity"] = _gate(
                "fail" if failures else "pass",
                "exact PCM and all metrics/null/support/frame counts for identical SHA256 libraries",
                failures=failures,
            )
        else:
            local["same_library_identity"] = _gate(
                "pass", "not applicable: distinct library hashes are not evidence of improvement", applicable=False
            )
        if checkpoint:
            clean_keys = (
                "tone_frames",
                "erasure_frames",
                "repeat_frames",
                "mute_frames",
                "c0_errors",
                "protected_errors",
                "c4_errors",
                "total_errors",
            )
            clean = [
                clip
                for clip in selected
                if all(clip[v] is not None and all(clip[v].get(k) == 0 for k in clean_keys) for v in VARIANTS)
            ]
            changed = [clip["name"] for clip in clean if clip["pcm_equal"] is not True]
            local["checkpoint_clean_speech_unchanged"] = _gate(
                "fail" if changed else "pass" if clean else "not_established",
                "the two production fixes must preserve checkpoint PCM on clean fixed-bit speech",
                clean_clips=len(clean),
                changed_clips=changed,
            )
        modes[mode] = {"status": _status(row["status"] for row in local.values()), "gates": local}
    return modes, equivalence_ok


def _plot_clip(root, manifest, clip, index):
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    reference = audio.read_pcm(support.artifact_path(root, clip["reference"]))
    signals = []
    gains = {}
    for variant in VARIANTS:
        decoded = audio.read_pcm(support.artifact_path(root, clip[f"{variant}_wav"]))
        r, d, gain = audio.normalized_signals(reference, decoded, clip["baseline"])
        if not signals:
            signals.append(r)
        signals.append(d)
        gains[variant] = gain
    labels = ("reference", "baseline", "candidate")
    colors = ("black", "#31688e", "#c44e52")
    powers = [audio.stft_power(signal, size=256, hop=80) for signal in signals]
    top = 10 * np.log10(max(float(np.max(powers[0])), 1e-20))
    bottom = top - 80
    fig, axes = plt.subplots(5, 2, figsize=(15, 18), constrained_layout=True)
    title = (
        f"{clip['name']} / {clip['mode']} / {clip['level']} — baseline lag {clip['baseline']['lag_samples']} samples"
    )
    fig.suptitle(title + "\nEqual reference-active RMS; regenerated phase is not optimized")
    t = np.arange(len(signals[0])) / 8000
    for signal, label, color in zip(signals, labels, colors):
        axes[0, 0].plot(t, signal, color=color, label=label, alpha=0.7, linewidth=0.5)
    axes[0, 0].set(title="Common-support waveform", xlabel="Time (s)", ylabel="Matched PCM units")
    axes[0, 0].legend()
    image = None
    for axis, power, label in zip((axes[0, 1], axes[1, 0], axes[1, 1]), powers, labels):
        image = axis.imshow(
            10 * np.log10(np.maximum(power.T, 1e-20)),
            origin="lower",
            aspect="auto",
            extent=(128 / 8000, (128 + 80 * len(power)) / 8000, 0, 4000),
            vmin=bottom,
            vmax=top,
            cmap="magma",
            interpolation="nearest",
        )
        axis.set(title=f"{label}: 256-point Hann / 80-hop", xlabel="Time (s)", ylabel="Frequency (Hz)")
    fig.colorbar(image, ax=[axes[0, 1], axes[1, 0], axes[1, 1]], label="STFT power (dB), one reference-derived scale")
    if len(signals[0]) >= 1024:
        harmonic = [audio.stft_power(signal, size=1024, hop=160) for signal in signals]
        strongest = int(np.argmax(np.sum(harmonic[0], axis=1)))
        frequency = np.fft.rfftfreq(1024, 1 / 8000)
        floor = max(float(np.max(harmonic[0])) * 1e-10, 1e-20)
        for power, label, color in zip(harmonic, labels, colors):
            axes[2, 0].plot(frequency, 10 * np.log10(power[strongest] + floor), label=label, color=color, linewidth=0.8)
            axes[2, 1].plot(
                frequency, 10 * np.log10(np.mean(power, axis=0) + floor), label=label, color=color, linewidth=0.8
            )
        axes[2, 0].set(
            title=f"Harmonic detail: reference-selected frame at {strongest * 0.02:.3f}s",
            xlabel="Frequency (Hz)",
            ylabel="1024-Hann / 160-hop power (dB)",
            xlim=(0, 2000),
        )
        axes[2, 1].set(
            title="Full-overlap long-term 1024-Hann / 160-hop spectra",
            xlabel="Frequency (Hz)",
            ylabel="Mean power (dB)",
        )
        axes[2, 0].legend()
    else:
        for axis in axes[2]:
            axis.text(0.5, 0.5, "Fewer than 1024 aligned samples: harmonic detail unavailable", ha="center", wrap=True)
            axis.set_axis_off()
    band_x = np.arange(len(support.BAND_METRICS))
    for offset, variant, color in ((-0.18, "baseline", colors[1]), (0.18, "candidate", colors[2])):
        axes[3, 0].bar(
            band_x + offset,
            [clip[variant][key] for key in support.BAND_METRICS],
            width=0.36,
            label=variant,
            color=color,
        )
    axes[3, 0].set_xticks(band_x, ["0–500", "500–1000", "1–2k", "2–3k", "3–4k"], rotation=20)
    axes[3, 0].set(
        title="Authoritative active long-term band differences", xlabel="Band (Hz)", ylabel="Decoded / reference (dB)"
    )
    axes[3, 0].axhline(0, color="black", linewidth=0.5)
    axes[3, 0].legend()
    axes[3, 1].set_axis_off()
    details = ["Per-clip threshold diagnostics:"] + (
        clip["threshold_exceedances"] or ["No per-clip screen exceedances"]
    )
    details += [
        f"Rail increases: {json.dumps(clip['rail_increases'])}",
        f"STOI baseline/candidate: {clip['stoi']['baseline']:.6f} / {clip['stoi']['candidate']:.6f}",
        f"ESTOI baseline/candidate: {clip['estoi']['baseline']:.6f} / {clip['estoi']['candidate']:.6f}",
        "Join/boundary metrics are not standardized click detectors.",
    ]
    axes[3, 1].text(0, 1, "\n".join(details), va="top", fontsize=8, wrap=True)
    activity = audio.reference_activity(signals[0])
    offset = clip["baseline"]["decoded_trim_head"]
    joins = [
        j
        for j in range(160 - offset % 160, len(signals[0]) - 4, 160)
        if j > 4 and j // 160 < len(activity) and activity[j // 160]
    ]
    # Inspect largest candidate join step and largest baseline join step. These
    # choose diagnostic views, not an alignment or a pass/fail click threshold.
    chosen = [max(joins, key=lambda j: abs(signals[v][j] - signals[v][j - 1])) for v in (2, 1)] if joins else []
    for axis, center, label in zip(axes[4], chosen, ("largest candidate join step", "largest baseline join step")):
        lo, hi = max(0, center - 80), min(len(signals[0]), center + 80)
        for signal, variant, color in zip(signals, labels, colors):
            axis.plot(np.arange(lo, hi), signal[lo:hi], label=variant, color=color, linewidth=0.8)
        for join in range(160 - offset % 160, hi, 160):
            if join >= lo:
                axis.axvline(join, color="gray", linestyle="--", linewidth=0.7)
        axis.set(
            title=f"160-sample-grid join diagnostic: {label}",
            xlabel="Aligned reference sample",
            ylabel="Matched PCM units",
        )
    if not chosen:
        for axis in axes[4]:
            axis.text(0.5, 0.5, "No supported active joins", ha="center")
    relative = f"analysis/figures/clip-{index:04d}.png"
    try:
        fig.savefig(_output_path(root, manifest, relative), dpi=130)
    finally:
        plt.close(fig)
    return relative, gains


def _write_tables(root, manifest, clips, groups):
    clip_path = "analysis/per-clip.csv"
    fields = ["name", "mode", "speaker", "utterance_index", "level", "valid", "included_in_aggregate"]
    fields += [
        f"{v}_{key}"
        for v in VARIANTS
        for key in METRICS + ("stoi", "estoi", "pcm_rail_samples", "pcm_max_rail_run", "pcm_rail_rate")
    ]
    fields += ["threshold_exceedances", "invalid_reasons", "rail_increases", "plots"]
    with _output_path(root, manifest, clip_path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for clip in clips:
            row = {key: clip[key] for key in fields[:7]}
            for variant in VARIANTS:
                for key in METRICS + ("pcm_rail_samples", "pcm_max_rail_run"):
                    row[f"{variant}_{key}"] = (clip[variant] or {}).get(key)
                for key in ("stoi", "estoi"):
                    row[f"{variant}_{key}"] = clip[key][variant]
                row[f"{variant}_pcm_rail_rate"] = (clip["rails"][variant] or {}).get("rate")
            for key in ("threshold_exceedances", "invalid_reasons", "rail_increases", "plots"):
                row[key] = json.dumps(clip[key], allow_nan=False)
            writer.writerow(row)
    grouped_path = "analysis/grouped.csv"
    fields = ["grouping", "mode", "speaker", "level", "clips", "valid_clips", "objective_nonregression"]
    fields += [
        f"{v}_{key}"
        for v in VARIANTS
        for key in METRICS + ("stoi", "estoi", "crest_mean_absolute_error_db", "crest_absolute_signed_mean_db")
    ]
    fields += ["failed_screens", "worst_clips", "rail_increases"]
    with _output_path(root, manifest, grouped_path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for group in groups:
            row = {key: group[key] for key in fields[:7]}
            for variant in VARIANTS:
                for key in METRICS:
                    row[f"{variant}_{key}"] = (group["means"] or {}).get(variant, {}).get(key)
                for key in ("stoi", "estoi"):
                    row[f"{variant}_{key}"] = group[key][variant]
                row[f"{variant}_crest_mean_absolute_error_db"] = group.get(
                    "crest_mean_perclip_absolute_error_db", {}
                ).get(variant)
                row[f"{variant}_crest_absolute_signed_mean_db"] = group.get(
                    "crest_absolute_signed_mean_error_db", {}
                ).get(variant)
            row["failed_screens"] = json.dumps(group["reasons"])
            row["worst_clips"] = json.dumps(group["worst_clips"])
            row["rail_increases"] = json.dumps(group["rail_increases"])
            writer.writerow(row)
    return [clip_path, grouped_path]


def report_run(run_dir, listening_results=None):
    """Produce a complete positive/negative assessment without modifying inputs.

    Missing offline packages or malformed immutable snapshots are prerequisites,
    not quality results. Failed measurement and screening gates are retained in
    the report and acceptance output rather than aborting a negative assessment.
    """
    versions = audio.package_versions()
    import quality_listening

    root = Path(run_dir).resolve()
    manifest = support.load_manifest(root)
    manifest_sha = support.digest(root / "manifest.json")
    independent, excluded = _independent_names(root, manifest)
    coverage_errors = _coverage(manifest)
    clips, crosschecks = _load_clips(root, manifest, independent)
    correctness, shared_pcm = _correctness(root, manifest, clips, crosschecks, coverage_errors)
    level = manifest["corpus"].get("level", "unknown")
    groups, by_mode = [], {}
    for mode in support.MODES:
        selected = [clip for clip in clips if clip["mode"] == mode]
        aggregate = _aggregate(selected)
        if coverage_errors:
            aggregate["objective_nonregression"] = "fail"
            aggregate["reasons"].extend(coverage_errors)
        by_mode[mode] = aggregate
        groups.append({"grouping": "mode/level", "mode": mode, "speaker": None, "level": level, **aggregate})
        speakers = defaultdict(list)
        for clip in selected:
            speakers[clip["speaker"]].append(clip)
        for speaker, subset in sorted(speakers.items(), key=lambda pair: pair[0] or ""):
            groups.append(
                {
                    "grouping": "speaker/mode/level",
                    "mode": mode,
                    "speaker": speaker,
                    "level": level,
                    **_aggregate(subset),
                }
            )
    generated = []
    plot_records = []
    for mode in support.MODES:
        valid = [clip for clip in clips if clip["mode"] == mode and clip["valid"]]
        flagged = [clip for clip in valid if clip["threshold_exceedances"] or clip["rail_increases"]]
        chosen = flagged or valid[:1]
        for clip in chosen:
            relative, gains = _plot_clip(root, manifest, clip, len(plot_records) + 1)
            clip["plots"].append(relative)
            generated.append(relative)
            plot_records.append(
                {
                    "name": clip["name"],
                    "mode": mode,
                    "path": relative,
                    "active_rms_match_gains": gains,
                    "selection": "flagged" if clip in flagged else "representative",
                }
            )
    generated.extend(_write_tables(root, manifest, clips, groups))
    profile = manifest.get("gain_profile")
    if profile is not None:
        # Read/validate only: fixed-profile bytes are immutable and never serialized again.
        support.load_gain_profile(support.artifact_path(root, profile["path"]))
        gain_profile = {**profile, "method": "unchanged supplied fixed profile"}
    else:
        frozen = support.make_gain_profile(manifest, root)
        gain_profile = None
        if frozen is not None:
            relative = "gain-profile.json"
            _json(_output_path(root, manifest, relative), frozen)
            generated.append(relative)
            gain_profile = {
                "path": relative,
                "sha256": support.digest(root / relative),
                "method": "original-main independent legacy calibration medians",
            }
    prepared = quality_listening.prepare_listening(root, manifest, clips)
    listening = quality_listening.analyze_listening(prepared, listening_results)
    generated.extend(prepared["artifacts"])
    shared_pcm = shared_pcm and prepared["imbe_shared_pcm"]
    dependencies = {
        "imbe7100": {
            "speech_mode": "imbe7200",
            "exact_shared_pcm": shared_pcm,
            "independent_sample": False,
            "reason": "same canonical parameters through separately exercised frame APIs",
        }
    }
    modes = {}
    for mode in support.MODES:
        modes[mode] = {
            "correctness": correctness[mode]["status"],
            "objective_nonregression": by_mode[mode]["objective_nonregression"],
            "perceptual": listening["modes"][mode]["status"],
            "reasons": {
                "correctness": [
                    f"{key}: {gate['reason']}"
                    for key, gate in correctness[mode]["gates"].items()
                    if gate["status"] != "pass"
                ],
                "objective_nonregression": by_mode[mode]["reasons"],
                "perceptual": listening["modes"][mode].get("reason"),
            },
        }
        if mode in dependencies:
            modes[mode]["dependency"] = dependencies[mode]
    acceptance = {
        "schema_version": 2,
        "manifest_sha256": manifest_sha,
        **{axis: _status(row[axis] for row in modes.values()) for axis in AXES},
        "modes": modes,
        "dependencies": dependencies,
        "limitations": LIMITATIONS,
        "level": level,
        "assessment_complete": True,
        "meaning": "A completed assessment may fail or leave human preference not established; no DSP/gain retuning is authorized.",
    }
    source_files = (
        "quality_report.py",
        "quality_audio.py",
        "quality_listening.py",
        "quality_support.py",
        "analyze_quality.py",
    )
    sources = {}
    for name in source_files:
        path = Path(__file__).resolve().with_name(name)
        if path.is_file():
            sources[name] = {"path": str(path), "sha256": support.digest(path)}
    artifacts = {relative: support.digest(support.artifact_path(root, relative)) for relative in sorted(set(generated))}
    analysis = {
        "schema_version": 2,
        "manifest_sha256": manifest_sha,
        "corpus": manifest["corpus"],
        "duplicate_exclusions": excluded,
        "aggregation": "unweighted paired clip means; verified hts1a exact prefix excluded only from independent statistics; levels never pooled",
        "metric_authority": "immutable C evaluator metric snapshots; independent NumPy LSD cross-check only",
        "alignment": "baseline automatic lag used for both variants, full aligned overlap for STOI/ESTOI, no concatenated VAD snippets",
        "provenance": {
            "packages": versions,
            "python": platform.python_version(),
            "platform": platform.platform(),
            "analysis_sources": sources,
            "libraries": manifest["libraries"],
            "tools": manifest["tools"],
            "comparison_points": manifest.get("comparison_points"),
        },
        "clips": clips,
        "groups": groups,
        "modes": by_mode,
        "correctness": correctness,
        "numpy_lsd_crosschecks": crosschecks,
        "plots": plot_records,
        "unavailable_plots": [
            {"name": clip["name"], "mode": clip["mode"], "reasons": clip["invalid_reasons"]}
            for clip in clips
            if not clip["valid"]
        ],
        "gain_profile": gain_profile,
        "listening": {
            "preparation": {key: value for key, value in prepared.items() if key != "items"},
            "assessment": listening,
        },
        "acceptance": acceptance,
        "generated_output_sha256": artifacts,
        "hash_scope": "generated figures/tables/listening/profile; excludes analysis.json and acceptance.json to avoid recursive hashes",
        "limitations": LIMITATIONS,
    }
    if support.digest(root / "manifest.json") != manifest_sha:
        raise ValueError("immutable manifest changed during report generation")
    _json(_output_path(root, manifest, "analysis.json"), analysis)
    _json(_output_path(root, manifest, "acceptance.json"), acceptance)
    return {
        "run_dir": str(root),
        "analysis": "analysis.json",
        "acceptance": "acceptance.json",
        **{axis: acceptance[axis] for axis in AXES},
        "modes": modes,
        "clips": len(clips),
        "figures": len(plot_records),
        "gain_profile": gain_profile,
        "listening_csv": prepared["listener_csv"],
    }
