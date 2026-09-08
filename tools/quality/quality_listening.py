#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Blinded engineering comparisons, not an ITU-compliant MOS/MUSHRA study.

Only submitted human judgments enter acceptance. Numerical bootstrap callers may
supply arrays directly; that helper neither creates nor submits listener rows.
"""

import csv
import json
import random
from pathlib import Path

import numpy as np
from quality_audio import normalized_signals, read_pcm, write_wav
from quality_support import MODES, artifact_path, digest

SEED = 0x12345678
REPLICATES = 10000
SPEAKERS = ("bdl", "rms", "clb", "slt")
CSV_COLUMNS = (
    "listener_id",
    "item_id",
    "preference",
    "naturalness_a",
    "naturalness_b",
    "intelligibility_a",
    "intelligibility_b",
    "artifact_notes",
)
RATING_COLUMNS = CSV_COLUMNS[3:7]


def hierarchical_bootstrap(values, speakers, *, replicates=REPLICATES, seed=SEED):
    """Return crossed listener/clip bootstrap means and percentile intervals.

    ``values`` has shape (listeners, clips[, metrics]); every entry must be
    observed and finite. Listeners and clips are sampled independently with
    replacement; clip sampling preserves each speaker's original item count.
    Multiple metrics share each replicate's draws. ``ci`` is [low, high] for
    a two-dimensional input, or a list of such intervals in metric order.
    """
    scores = np.asarray(values, dtype=np.float64)
    scalar = scores.ndim == 2
    if scalar:
        scores = scores[:, :, None]
    if scores.ndim != 3 or any(size == 0 for size in scores.shape):
        raise ValueError("bootstrap requires nonempty listeners x clips [x metrics]")
    if not np.isfinite(scores).all():
        raise ValueError("bootstrap cannot accept missing or nonfinite judgments")
    labels = list(speakers)
    if len(labels) != scores.shape[1] or any(not isinstance(s, str) or not s for s in labels):
        raise ValueError("bootstrap requires one nonempty speaker label per clip")
    if type(replicates) is not int or replicates <= 0:
        raise ValueError("bootstrap replicate count must be a positive integer")
    strata = [np.array([i for i, label in enumerate(labels) if label == speaker]) for speaker in sorted(set(labels))]
    rng = np.random.default_rng(seed)
    means = np.empty((replicates, scores.shape[2]), dtype=np.float64)
    # Bound the temporary crossed sample array, including for larger panels.
    batch = max(1, min(256, 1000000 // scores.size))
    for start in range(0, replicates, batch):
        count = min(batch, replicates - start)
        listeners = rng.integers(scores.shape[0], size=(count, scores.shape[0]))
        clips = np.concatenate([group[rng.integers(len(group), size=(count, len(group)))] for group in strata], axis=1)
        means[start : start + count] = scores[listeners[:, :, None], clips[:, None, :], :].mean(axis=(1, 2))
    intervals = np.percentile(means, [2.5, 97.5], axis=0).T
    mean = scores.mean(axis=(0, 1))
    return {
        "mean": float(mean[0]) if scalar else mean.tolist(),
        "ci": intervals[0].tolist() if scalar else intervals.tolist(),
        "replicates": replicates,
        "seed": seed,
    }


def _protocol(manifest):
    return {
        "study": "Small blinded engineering comparison; not ITU-compliant MOS/MUSHRA",
        "corpus_set": manifest["corpus"]["set"],
        "level": manifest["corpus"]["level"],
        "speakers": list(SPEAKERS),
        "utterance_indices": [1, 2, 3],
        "balanced_items_per_mode": 12,
        "minimum_complete_listeners": 5,
        "seed": SEED,
        "bootstrap_replicates": REPLICATES,
        "bootstrap": "Independent crossed listener and speaker-stratified clip resampling with replacement",
        "interval": "2.5th and 97.5th percentiles",
        "preference": "Candidate minus baseline: +1 preference, -1 rejection, 0 tie",
        "ratings": "Integer 1..5; higher is better; candidate minus baseline differences",
        "pass_rule": "Preference lower bound > 0 and intelligibility upper bound >= 0",
        "fail_rule": "Preference or intelligibility upper bound < 0",
        "naturalness": "Diagnostic only; no extra acceptance gate",
        "levels": "Each run/level is analyzed separately, never pooled",
        "legacy": "Exceedance items are diagnostic only and excluded from mandatory bootstrap",
        "instructions": "Listen to reference/A/B, then submit one CSV row per listener/item. Use stable human listener IDs; preference A, B, or tie. Keep the private answer key away from listeners.",
        "human_provenance": "The submitting investigator must attest rows are actual human judgments; software cannot authenticate listeners.",
        "presentation": "Common aligned support; equal active reference RMS before one common peak-safe attenuation; s16 mono 8000 Hz",
    }


def _shared_pcm_evidence(run_dir, manifest, clips):
    """Require every named canonical/frame control, not vacuous any/all success."""
    references = manifest.get("references", [])
    names = [ref["name"] for ref in references]
    if not names or len(set(names)) != len(names):
        return False
    expected = {(name, variant) for name in names for variant in ("baseline", "candidate")}
    controls = manifest.get("frame_equivalence", [])
    if len(controls) != len(expected):
        return False
    if {(row.get("name"), row.get("variant")) for row in controls} != expected:
        return False
    indexed = {}
    for clip in clips:
        key = (clip["name"], clip["mode"])
        if key in indexed:
            return False
        indexed[key] = clip
    for row in controls:
        if any(row.get(key) is not True for key in ("imbe7200_equal", "imbe7100_equal", "frames_equal")) or row.get(
            "invalid_reasons"
        ):
            return False
        if any(row.get(key) not in manifest["files"] for key in ("canonical_wav", "canonical_json")):
            return False
        first = indexed.get((row["name"], "imbe7200"))
        second = indexed.get((row["name"], "imbe7100"))
        if not first or not second or not first["valid"] or not second["valid"]:
            return False
        variant = row["variant"]
        if any(not isinstance(clip.get(variant), dict) for clip in (first, second)):
            return False
        support = ("lag_samples", "reference_trim_head", "decoded_trim_head", "aligned_samples", "frames")
        if any(
            first[variant].get(key) is None or first[variant].get(key) != second[variant].get(key) for key in support
        ):
            return False
        try:
            canonical = read_pcm(artifact_path(run_dir, row["canonical_wav"]))
            for clip in (first, second):
                path = clip[f"{variant}_wav"]
                if path not in manifest["files"] or not np.array_equal(
                    canonical, read_pcm(artifact_path(run_dir, path))
                ):
                    return False
        except (OSError, ValueError):
            return False
    return True


def _item_audio(run_dir, clip):
    if not clip["valid"]:
        raise ValueError("clip has invalid supported measurements")
    baseline, candidate = clip["baseline"], clip["candidate"]
    if not isinstance(baseline, dict) or not isinstance(candidate, dict):
        raise TypeError("missing baseline or candidate measurements")
    for key in ("lag_samples", "reference_trim_head", "decoded_trim_head", "aligned_samples"):
        if baseline.get(key) is None or baseline[key] != candidate.get(key):
            raise ValueError(f"baseline/candidate lack common aligned support: {key}")
    reference = read_pcm(artifact_path(run_dir, clip["reference"]))
    r, b, b_gain = normalized_signals(reference, read_pcm(artifact_path(run_dir, clip["baseline_wav"])), baseline)
    other_r, c, c_gain = normalized_signals(
        reference, read_pcm(artifact_path(run_dir, clip["candidate_wav"])), candidate
    )
    if not np.array_equal(r, other_r) or r.shape != b.shape or r.shape != c.shape:
        raise ValueError("reference and both variants must have identical listening support")
    if not all(np.isfinite(signal).all() for signal in (r, b, c)):
        raise ValueError("nonfinite listening PCM")
    peak = max(float(np.max(np.abs(signal))) for signal in (r, b, c))
    if peak <= 0:
        raise ValueError("no audible speech support")
    # Symmetric 32767 headroom also avoids the asymmetric negative int16 rail.
    attenuation = min(1.0, 32767.0 / peak)
    return (
        r * attenuation,
        b * attenuation,
        c * attenuation,
        {"reference": 1.0, "baseline": float(b_gain), "candidate": float(c_gain)},
        attenuation,
    )


def prepare_listening(run_dir, manifest, clips):
    """Write blinded materials; retain missing/invalid items as unavailable."""
    root = Path(run_dir)
    protocol = _protocol(manifest)
    selected = []
    unavailable = []
    corpus_set, level = protocol["corpus_set"], protocol["level"]
    if corpus_set == "heldout":
        for mode in MODES:
            for speaker in SPEAKERS:
                for index in (1, 2, 3):
                    descriptor = {
                        "mode": mode,
                        "speaker": speaker,
                        "utterance_index": index,
                        "level": level,
                        "kind": "balanced_heldout",
                    }
                    matches = [
                        clip
                        for clip in clips
                        if clip["mode"] == mode
                        and clip["speaker"] == speaker
                        and clip["utterance_index"] == index
                        and clip["level"] == level
                    ]
                    if len(matches) != 1:
                        unavailable.append(
                            {
                                **descriptor,
                                "name": None,
                                "reason": "missing clip" if not matches else "ambiguous duplicate clips",
                            }
                        )
                    else:
                        selected.append((matches[0], descriptor))
    elif corpus_set == "legacy":
        for clip in sorted(clips, key=lambda entry: (entry["mode"], entry["name"])):
            if clip["threshold_exceedances"]:
                selected.append(
                    (
                        clip,
                        {
                            "mode": clip["mode"],
                            "speaker": clip["speaker"],
                            "utterance_index": clip["utterance_index"],
                            "level": clip["level"],
                            "kind": "legacy_exceedance",
                        },
                    )
                )
    rng = random.Random(SEED)
    rng.shuffle(selected)
    presentation = root / "listening" / "presentation"
    private = root / "listening" / "private"
    presentation.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    items, artifacts = [], []
    for clip, descriptor in selected:
        item_id = f"{rng.getrandbits(128):032x}"
        a_variant = "candidate" if rng.getrandbits(1) else "baseline"
        b_variant = "baseline" if a_variant == "candidate" else "candidate"
        record = {
            **descriptor,
            "item_id": item_id,
            "name": clip["name"],
            "a_variant": a_variant,
            "b_variant": b_variant,
            "threshold_exceedances": list(clip["threshold_exceedances"]),
        }
        try:
            reference, baseline, candidate, gains, attenuation = _item_audio(root, clip)
        except (OSError, ValueError, TypeError, KeyError) as exc:
            unavailable.append({**record, "reason": str(exc)})
            continue
        signals = {"reference": reference, "baseline": baseline, "candidate": candidate}
        wavs = {}
        for suffix, variant in (("reference", "reference"), ("A", a_variant), ("B", b_variant)):
            path = presentation / f"{item_id}_{suffix}.wav"
            write_wav(path, signals[variant])
            wavs[suffix] = path.relative_to(root).as_posix()
            artifacts.append(wavs[suffix])
        record.update(
            wavs=wavs,
            active_rms_gains=gains,
            common_attenuation=attenuation,
            common_attenuation_db=float(20.0 * np.log10(attenuation)),
            samples=len(reference),
            available=True,
        )
        items.append(record)
    imbe_shared = _shared_pcm_evidence(root, manifest, clips)
    protocol["imbe_dependency"] = (
        "IMBE7100 inherits IMBE7200 only: all canonical/frame PCM and frame-count controls match; no independent sample-size increment"
        if imbe_shared
        else "Exact shared IMBE PCM evidence incomplete; no inherited claim permitted"
    )
    protocol["available_items"] = len(items)
    protocol["unavailable_items"] = len(unavailable)
    answer_key = (private / "answer-key.json").relative_to(root).as_posix()
    listener_csv = "listening/listener.csv"
    key = {
        "schema_version": 1,
        "items": items,
        "unavailable": unavailable,
        "protocol": protocol,
        "imbe_shared_pcm": imbe_shared,
    }
    (root / answer_key).write_text(json.dumps(key, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    # Do not erase real judgments if an investigator has filled the template.
    if not (root / listener_csv).exists():
        with (root / listener_csv).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            writer.writerows({"item_id": item["item_id"]} for item in items)
    artifacts.extend((answer_key, listener_csv))
    return {**key, "answer_key": answer_key, "listener_csv": listener_csv, "artifacts": artifacts}


def _read_judgments(path, items):
    judgments, incomplete, seen = {}, [], set()
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != list(CSV_COLUMNS):
            raise ValueError("listener CSV must have exactly these ordered columns: " + ",".join(CSV_COLUMNS))
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"listener CSV line {line}: wrong number of columns")
            row = {key: value.strip() for key, value in row.items()}
            if not any(row.values()):
                continue
            listener, item_id = row["listener_id"], row["item_id"]
            if not item_id:
                raise ValueError(f"listener CSV line {line}: item_id is required")
            if item_id not in items:
                raise ValueError(f"listener CSV line {line}: unknown or unavailable item {item_id!r}")
            if not listener and not any(row[field] for field in CSV_COLUMNS if field != "item_id"):
                continue  # An untouched template row is not a judgment.
            if not listener:
                raise ValueError(f"listener CSV line {line}: submitted judgments require listener_id")
            pair = (listener, item_id)
            if pair in seen:
                raise ValueError(f"listener CSV line {line}: duplicate listener/item judgment")
            seen.add(pair)
            if row["preference"] and row["preference"] not in ("A", "B", "tie"):
                raise ValueError(f"listener CSV line {line}: preference must be A, B, or tie")
            for field in RATING_COLUMNS:
                if row[field] and row[field] not in ("1", "2", "3", "4", "5"):
                    raise ValueError(f"listener CSV line {line}: {field} must be an integer 1..5")
            if any(not row[field] for field in ("preference", *RATING_COLUMNS)):
                incomplete.append(
                    {
                        "listener_id": listener,
                        "item_id": item_id,
                        "line": line,
                        "reason": "incomplete preference or ratings",
                    }
                )
                continue
            item = items[item_id]
            sign = 1 if item["a_variant"] == "candidate" else -1
            preference = 0 if row["preference"] == "tie" else sign * (1 if row["preference"] == "A" else -1)
            judgments[pair] = (
                preference,
                sign * (int(row["naturalness_a"]) - int(row["naturalness_b"])),
                sign * (int(row["intelligibility_a"]) - int(row["intelligibility_b"])),
            )
    return judgments, incomplete


def _mode_result(items, judgments, mode, level):
    balanced = sorted(
        (
            item
            for item in items
            if item["mode"] == mode and item["kind"] == "balanced_heldout" and item["level"] == level
        ),
        key=lambda item: (item["speaker"], item["utterance_index"]),
    )
    ids = [item["item_id"] for item in balanced]
    listeners = sorted({listener for listener, _ in judgments})
    complete = [listener for listener in listeners if ids and all((listener, item_id) in judgments for item_id in ids)]
    result = {
        "status": "not_established",
        "reason": "Twelve balanced heldout items are unavailable",
        "complete_listeners": 0,
        "items": len(ids),
        "preference_ci": None,
        "naturalness_ci": None,
        "intelligibility_ci": None,
        "level": level,
        "eligible_listener_ids": [],
        "partial_listener_ids": [],
    }
    required = {(speaker, index) for speaker in SPEAKERS for index in (1, 2, 3)}
    if len(balanced) != 12 or {(item["speaker"], item["utterance_index"]) for item in balanced} != required:
        return result
    result["complete_listeners"] = len(complete)
    result["eligible_listener_ids"] = complete
    result["partial_listener_ids"] = [
        listener
        for listener in listeners
        if listener not in complete and any((listener, item_id) in judgments for item_id in ids)
    ]
    if len(complete) < 5:
        result["reason"] = "At least five human listeners must each rate all twelve balanced items"
        return result
    values = np.array([[judgments[(listener, item_id)] for item_id in ids] for listener in complete])
    estimate = hierarchical_bootstrap(values, [item["speaker"] for item in balanced])
    for index, metric in enumerate(("preference", "naturalness", "intelligibility")):
        result[f"{metric}_ci"] = estimate["ci"][index]
        result[f"{metric}_mean"] = estimate["mean"][index]
    preference, intelligibility = result["preference_ci"], result["intelligibility_ci"]
    if preference[1] < 0 or intelligibility[1] < 0:
        result.update(status="fail", reason="Preference or intelligibility interval is wholly below zero")
    elif preference[0] > 0:
        result.update(
            status="pass",
            reason="Preference interval strictly favors candidate; intelligibility interval is not wholly below zero",
        )
    else:
        result["reason"] = "Preference interval does not establish candidate benefit or harm"
    return result


def analyze_listening(prepared, results_path=None):
    """Assess submitted judgments, keeping all absent evidence unestablished."""
    items = prepared["items"]
    indexed = {item["item_id"]: item for item in items}
    if len(indexed) != len(items):
        raise ValueError("prepared listening items contain duplicate IDs")
    for item in items:
        if {item["a_variant"], item["b_variant"]} != {"baseline", "candidate"}:
            raise ValueError("prepared listening answer key has invalid A/B variants")
    judgments, incomplete, results_hash = {}, [], None
    missing_reason = "No human listening results supplied"
    if results_path is not None:
        path = Path(results_path)
        if path.exists():
            results_hash = digest(path)
            judgments, incomplete = _read_judgments(path, indexed)
            missing_reason = "No complete human judgments in supplied results"
        else:
            missing_reason = f"Human listening results file is missing: {path}"
    protocol = dict(prepared["protocol"])
    independent = ("imbe7200", "ambe2450", "ambe2400") if prepared["imbe_shared_pcm"] else MODES
    modes = {mode: _mode_result(items, judgments, mode, protocol["level"]) for mode in independent}
    # Never pool two labels of identical IMBE PCM as independent observations.
    if prepared["imbe_shared_pcm"]:
        modes["imbe7100"] = {
            **modes["imbe7200"],
            "dependency": "Inherited from imbe7200: exact canonical and both frame-path PCM/count controls for every run reference; no extra independent samples",
            "direct_7100_judgments_used": False,
        }
    else:
        modes["imbe7100"]["dependency"] = (
            "Shared PCM evidence incomplete; independently assessed, not inferred from imbe7200"
        )
    states = [modes[mode]["status"] for mode in independent]
    status = "fail" if "fail" in states else "pass" if all(state == "pass" for state in states) else "not_established"
    return {
        "status": status,
        "modes": modes,
        "protocol": protocol,
        "results_sha256": results_hash,
        "incomplete_rows": incomplete,
        "complete_judgment_rows": len(judgments),
        "diagnostic_judgment_rows": sum(indexed[item_id]["kind"] == "legacy_exceedance" for _, item_id in judgments),
        "independent_modes": list(independent),
        "reason": missing_reason
        if not judgments
        else "Per-mode blinded comparison using only complete balanced panels",
    }
