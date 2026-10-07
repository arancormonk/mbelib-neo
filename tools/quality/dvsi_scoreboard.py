#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Score this library against DVSI's AMBE-3000 test vectors.

For every speech vector of the D-STAR, P25 and AMBE+2 (rate 33) sets this decodes
DVSI's bits with mbe_quality_eval and measures, against the original input speech,
both our output and DVSI's own decoded output; it also compares our output with
DVSI's directly. A harmonic-peak pitch estimator checks that our decoded pitch
matches DVSI's and the input's on the same frames. For D-STAR it also encodes the
input with this library's encoder and compares pitch, voicing and spectral error
with DVSI's encoding.

The vectors are fetched locally with fetch_dvsi_vectors.py and must never be
committed or redistributed. This tool writes only aggregate numbers and per-run
work files under --out (git-ignored by default).

Needs NumPy for the pitch measurements.
"""

import argparse
import importlib.util
import json
import math
import random
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

MODES = {"dstar": "ambe2400", "p25": "imbe7200", "r33": "ambe2450"}
# Partitions by original utterance. Development and validation were both examined
# while earlier changes were made, so they are retrospective; "final" vectors are
# scored once per change set and never used for tuning.
PARTITIONS = {
    "development": ("dam", "clean", "fambf22c", "fambm22a", "t01", "t02", "tia11", "tambf22a", "ucarm15a"),
    "validation": ("t03", "t04", "tambf22b", "tambf32b", "tambf32e", "ucarf15a", "ucarm20a"),
    "final": ("irstia", "p01mirs", "mark"),
}
BAND_KEYS = tuple(
    f"band_delta_db_{band}" for band in ("0_500", "500_1000", "1000_2000", "2000_3000", "3000_4000", "3500_4000")
)
SPEECH_KEYS = ("level_offset_db", "lsd_db", "env_corr", "crest_delta_db", "pcm_rail_samples") + BAND_KEYS
F0_BANDS = ((0.0, 100.0), (100.0, 150.0), (150.0, 200.0), (200.0, 300.0), (300.0, 1000.0))
# MBE_PROCESS_FLAG_TONE, _ERASURE, _REPEAT, _MUTE and _SILENCE: frames not synthesized from their own voice model.
FLAG_SKIP = 0x0010 | 0x0020 | 0x0040 | 0x0080 | 0x0100
SAMPLE_RATE = 8000
FRAME = 160


def partition_of(name):
    for partition, names in PARTITIONS.items():
        if name in names:
            return partition
    return None


# --------------------------------------------------------------------------------------------
# Pitch measurement


def harmonic_f0(signal, center, f0_guess, length=400, nfft=16384, fmax=2000.0):
    """Fundamental (Hz) of the harmonic series near f0_guess around sample `center`, or None.

    A harmonic comb over 0.9..1.1 x f0_guess finds the series; each harmonic's peak
    is then located with parabolic interpolation and f0 is the weighted
    least-squares slope through the origin of peak frequency against harmonic number.
    """
    import numpy as np

    start = int(center) - length // 2
    if start < 0 or start + length > len(signal) or not f0_guess or f0_guess <= 0:
        return None
    segment = np.asarray(signal[start:start + length], dtype=np.float64)
    if math.sqrt(float(np.mean(segment * segment))) < 100.0:
        return None
    power = np.abs(np.fft.rfft(segment * np.hanning(length), nfft)) ** 2 + 1e-9
    log_power = np.log(power)
    bins = np.arange(len(power))
    count = max(3, int(fmax / f0_guess))
    harmonics = np.arange(1, count + 1)
    ratios = np.linspace(0.90, 1.10, 401)
    centers = np.outer(ratios * f0_guess, harmonics) * nfft / SAMPLE_RATE
    between = centers + 0.5 * (ratios * f0_guess)[:, None] * nfft / SAMPLE_RATE
    score = np.mean(np.interp(centers, bins, log_power) - np.interp(between, bins, log_power), axis=1)
    best = int(np.argmax(score))
    if score[best] < 2.0:
        return None
    f0_coarse = ratios[best] * f0_guess
    half_width = 0.3 * f0_coarse * nfft / SAMPLE_RATE
    numerator = denominator = 0.0
    for k in harmonics:
        position = k * f0_coarse * nfft / SAMPLE_RATE
        low, high = int(position - half_width), int(position + half_width) + 1
        if low < 1 or high >= len(power) - 1:
            continue
        peak = low + int(np.argmax(power[low:high]))
        a, b, c = log_power[peak - 1], log_power[peak], log_power[peak + 1]
        curvature = a - 2.0 * b + c
        offset = 0.5 * (a - c) / curvature if curvature < 0 else 0.0
        frequency = (peak + offset) * SAMPLE_RATE / nfft
        weight = float(power[peak])
        numerator += weight * k * frequency
        denominator += weight * k * k
    return numerator / denominator if denominator > 0 else None


def model_f0(model):
    return model["w0"] * SAMPLE_RATE / (2.0 * math.pi)


def stable_voiced_frames(records, tolerance=0.012):
    """Indices of voice frames with L >= 9, >= 80% voiced bands and w0 stable over +-2 frames."""
    usable = []
    for index in range(2, len(records) - 2):
        window = records[index - 2:index + 3]
        if any(record["flags"] & FLAG_SKIP or record["decoded"].get("status") != 0 for record in window):
            continue
        model = records[index]["decoded"]
        if "w0" not in model or model["L"] < 9 or sum(model["Vl"]) < 0.8 * model["L"]:
            continue
        w0 = model["w0"]
        if all("w0" in r["decoded"] and abs(r["decoded"]["w0"] / w0 - 1.0) <= tolerance for r in window):
            usable.append(index)
    return usable


def pitch_ratios(records, signal, shift):
    """(nominal f0, measured/nominal) for stable voiced frames; `shift` maps our sample index to `signal`."""
    pairs = []
    for index in stable_voiced_frames(records):
        nominal = model_f0(records[index]["decoded"])
        measured = harmonic_f0(signal, FRAME * index + FRAME // 2 - shift, nominal)
        if measured:
            pairs.append((nominal, measured / nominal))
    return pairs


def summarize_ratios(pairs):
    summary = {}
    for low, high in F0_BANDS:
        values = sorted(ratio for nominal, ratio in pairs if low <= nominal < high)
        key = f"{int(low)}_{int(high)}"
        if values:
            quartiles = statistics.quantiles(values, n=4) if len(values) > 1 else [values[0]] * 3
            summary[key] = {"n": len(values), "median": statistics.median(values), "q1": quartiles[0],
                            "q3": quartiles[2]}
        else:
            summary[key] = {"n": 0, "median": None, "q1": None, "q3": None}
    values = [ratio for _, ratio in pairs]
    summary["all"] = {"n": len(values), "median": statistics.median(values) if values else None}
    return summary


def max_ratio_deviation(summary, minimum=20):
    """Largest |median - 1| over f0 bands with at least `minimum` frames, or None."""
    deviations = [abs(entry["median"] - 1.0) for key, entry in summary.items()
                  if key != "all" and entry["n"] >= minimum]
    return max(deviations) if deviations else None


# --------------------------------------------------------------------------------------------
# Encoder comparisons


def band_voicing(model, bands=4):
    """Majority voicing of each 1 kHz band from a decoded model's per-harmonic flags."""
    f0 = model_f0(model)
    votes = [[0, 0] for _ in range(bands)]
    for harmonic, voiced in enumerate(model["Vl"], start=1):
        band = min(bands - 1, int(harmonic * f0 // 1000.0))
        votes[band][0] += voiced
        votes[band][1] += 1
    return [count and 2 * voiced >= count for voiced, count in votes]


def voicing_agreement(encoded, reference, frame_offset):
    """Per-band agreement between two decoded streams, encoded[f] against reference[f - frame_offset]."""
    agree = [0, 0, 0, 0]
    total = 0
    for index, record in enumerate(encoded):
        other = index - frame_offset
        if not 0 <= other < len(reference):
            continue
        a, b = record["decoded"], reference[other]["decoded"]
        if a.get("status") != 0 or b.get("status") != 0 or "Vl" not in a or "Vl" not in b:
            continue
        if record["flags"] & FLAG_SKIP or reference[other]["flags"] & FLAG_SKIP:
            continue
        for band, (x, y) in enumerate(zip(band_voicing(a), band_voicing(b))):
            agree[band] += x == y
        total += 1
    return [count / total for count in agree] if total else [None] * 4


def octave_errors(records, signal, shift):
    """Fraction of stable voiced frames whose input has a stronger harmonic comb at 2x or 0.5x the coded f0."""
    import numpy as np

    errors = checked = 0
    for index in stable_voiced_frames(records):
        nominal = model_f0(records[index]["decoded"])
        start = FRAME * index + FRAME // 2 - shift - 200
        if start < 0 or start + 400 > len(signal):
            continue
        segment = np.asarray(signal[start:start + 400], dtype=np.float64)
        if math.sqrt(float(np.mean(segment * segment))) < 100.0:
            continue
        log_power = np.log(np.abs(np.fft.rfft(segment * np.hanning(400), 16384)) ** 2 + 1e-9)
        bins = np.arange(len(log_power))

        def comb(f0):
            harmonics = np.arange(1, max(3, int(2000.0 / f0)) + 1)
            on = np.interp(harmonics * f0 * 16384 / SAMPLE_RATE, bins, log_power)
            off = np.interp((harmonics + 0.5) * f0 * 16384 / SAMPLE_RATE, bins, log_power)
            return float(np.mean(on - off))

        scores = {factor: max(comb(nominal * factor * r) for r in (0.97, 1.0, 1.03)) for factor in (0.5, 1.0, 2.0)}
        checked += 1
        errors += max(scores[0.5], scores[2.0]) > scores[1.0] + 1.0
    return errors / checked if checked else None


# --------------------------------------------------------------------------------------------
# Running the tools


def run(command):
    result = subprocess.run([str(part) for part in command], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, command))} failed: {result.stderr.strip()}")
    return result


def read_pcm(path):
    """int16 samples of a raw s16le file or a 16-bit mono WAV."""
    import numpy as np

    data = Path(path).read_bytes()
    if data[:4] == b"RIFF":
        position = 12
        while position + 8 <= len(data):
            chunk, size = data[position:position + 4], int.from_bytes(data[position + 4:position + 8], "little")
            if chunk == b"data":
                return np.frombuffer(data[position + 8:position + 8 + size], dtype="<i2").astype(np.float64)
            position += 8 + size + (size & 1)
        raise RuntimeError(f"{path}: WAV without data chunk")
    return np.frombuffer(data[: len(data) // 2 * 2], dtype="<i2").astype(np.float64)


def load_json(path):
    return json.loads(Path(path).read_text())


def load_records(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line]


def score_vector(tools, vectors, mode, name, work, encoder):
    """All measurements of one vector in one mode."""
    codec = MODES[mode]
    source = vectors / mode / f"{name}.bit"
    base = work / f"{mode}_{name}"
    # The evaluator takes .raw or .wav paths; DVSI's PCM files are raw s16le.
    speech = work / f"{name}_input.raw"
    dvsi_pcm = Path(f"{base}_dvsi.raw")
    shutil.copyfile(vectors / f"{name}.pcm", speech)
    shutil.copyfile(vectors / mode / f"{name}.pcm", dvsi_pcm)
    run([tools / "mbe_quality_reframe", "--codec", codec, "--from-dvsi", source, "--out", f"{base}.frames"])
    run([tools / "mbe_quality_eval", "--codec", codec, "--frames", f"{base}.frames", "--out", f"{base}_ours.wav",
         "--ref", speech, "--json", f"{base}_ours.json", "--params", f"{base}_ours.jsonl"])
    run([tools / "mbe_quality_eval", "--decoded", dvsi_pcm, "--ref", speech, "--json", f"{base}_dvsi.json"])
    run([tools / "mbe_quality_eval", "--decoded", f"{base}_ours.wav", "--ref", dvsi_pcm,
         "--json", f"{base}_vs_dvsi.json"])
    ours, dvsi, versus = (load_json(f"{base}_{tag}.json") for tag in ("ours", "dvsi", "vs_dvsi"))
    records = load_records(f"{base}_ours.jsonl")
    ours_pcm, dvsi_samples, speech_samples = read_pcm(f"{base}_ours.wav"), read_pcm(dvsi_pcm), read_pcm(speech)
    result = {
        "mode": mode, "name": name, "partition": partition_of(name),
        "ours": {key: ours.get(key) for key in SPEECH_KEYS + ("pcm_fnv1a",)},
        "dvsi": {key: dvsi.get(key) for key in SPEECH_KEYS},
        "vs_dvsi": {key: versus.get(key) for key in SPEECH_KEYS},
        "pitch": {
            "ours": summarize_ratios(pitch_ratios(records, ours_pcm, 0)),
            "dvsi": summarize_ratios(pitch_ratios(records, dvsi_samples, versus["lag_samples"])),
            "input": summarize_ratios(pitch_ratios(records, speech_samples, ours["lag_samples"])),
        },
    }
    if encoder and mode == "dstar":
        result["encoder"] = score_encoder(tools, speech, base, records, ours, speech_samples)
    return result


def score_encoder(tools, speech, base, dvsi_records, dvsi_eval, speech_samples):
    """Our D-STAR encoding of the input against DVSI's, both decoded by this library."""
    run([tools / "mbe_quality_encode", "--codec", "ambe2400", "--in", speech, "--out", f"{base}_enc.rows"])
    run([tools / "mbe_quality_eval", "--codec", "ambe2400", "--frames", f"{base}_enc.rows",
         "--out", f"{base}_enc.wav", "--ref", speech, "--json", f"{base}_enc.json",
         "--params", f"{base}_enc.jsonl"])
    encoded = load_json(f"{base}_enc.json")
    records = load_records(f"{base}_enc.jsonl")

    def cents(pairs):
        values = [abs(1200.0 * math.log2(ratio)) for _, ratio in pairs]
        return statistics.median(values) if values else None

    ours_pairs = pitch_ratios(records, speech_samples, encoded["lag_samples"])
    dvsi_pairs = pitch_ratios(dvsi_records, speech_samples, dvsi_eval["lag_samples"])
    frame_offset = round((encoded["lag_samples"] - dvsi_eval["lag_samples"]) / FRAME)
    return {
        "speech": {key: encoded.get(key) for key in SPEECH_KEYS},
        "pitch_cents_median": cents(ours_pairs),
        "dvsi_pitch_cents_median": cents(dvsi_pairs),
        "pitch_input_over_coded": summarize_ratios(ours_pairs),
        "dvsi_pitch_input_over_coded": summarize_ratios(dvsi_pairs),
        "octave_error_rate": octave_errors(records, speech_samples, encoded["lag_samples"]),
        "dvsi_octave_error_rate": octave_errors(dvsi_records, speech_samples, dvsi_eval["lag_samples"]),
        "voicing_agreement": voicing_agreement(records, dvsi_records, frame_offset),
        "frame_offset": frame_offset,
    }


# --------------------------------------------------------------------------------------------
# Aggregation, comparison and gates


def mean(values):
    values = [value for value in values if value is not None]
    return sum(values) / len(values) if values else None


def flatten(result):
    """Scalar metrics of one vector as dotted keys (the namespace the gates use)."""
    flat = {}
    for group in ("ours", "dvsi", "vs_dvsi"):
        for key, value in result[group].items():
            if isinstance(value, (int, float)):
                flat[f"{group}.{key}"] = value
    bands = [result["vs_dvsi"].get(f"band_delta_db_{band}") for band in ("2000_3000", "3000_4000")]
    if all(isinstance(value, (int, float)) for value in bands):
        flat["vs_dvsi.hf_abs_error_db"] = (abs(bands[0]) + abs(bands[1])) / 2.0
    for source, summary in result["pitch"].items():
        flat[f"pitch.{source}_ratio"] = summary["all"]["median"]
    encoder = result.get("encoder")
    if encoder:
        for key, value in encoder["speech"].items():
            if isinstance(value, (int, float)):
                flat[f"encoder.{key}"] = value
        for key in ("pitch_cents_median", "dvsi_pitch_cents_median", "octave_error_rate", "dvsi_octave_error_rate"):
            flat[f"encoder.{key}"] = encoder[key]
        agreement = [value for value in encoder["voicing_agreement"] if value is not None]
        flat["encoder.voicing_agreement_mean"] = mean(agreement)
    return flat


def pooled_pitch(results, source):
    pairs = []
    for result in results:
        for key, entry in result["pitch"][source].items():
            if key != "all" and entry["n"]:
                pairs.append((key, entry["n"], entry["median"]))
    summary = {}
    for low, high in F0_BANDS:
        key = f"{int(low)}_{int(high)}"
        entries = [(n, median) for band, n, median in pairs if band == key]
        count = sum(n for n, _ in entries)
        summary[key] = {"n": count, "median": (sum(n * m for n, m in entries) / count) if count else None}
    return summary


def aggregate(results):
    summary = {}
    for mode in MODES:
        for partition in PARTITIONS:
            group = [result for result in results if result["mode"] == mode and result["partition"] == partition]
            if not group:
                continue
            flats = [flatten(result) for result in group]
            keys = sorted({key for flat in flats for key in flat})
            entry = {"vectors": [result["name"] for result in group],
                     "mean": {key: mean([flat.get(key) for flat in flats]) for key in keys}}
            for source in ("ours", "dvsi", "input"):
                pooled = pooled_pitch(group, source)
                entry[f"pitch_{source}_by_f0"] = pooled
                entry["mean"][f"pitch.{source}_max_band_deviation"] = max_ratio_deviation(pooled)
            summary[f"{mode}/{partition}"] = entry
    return summary


def bootstrap_interval(deltas, replicates=2000, seed=0x5EED):
    if len(deltas) < 2:
        return None
    rng = random.Random(seed)
    means = sorted(sum(rng.choice(deltas) for _ in deltas) / len(deltas) for _ in range(replicates))
    return means[int(0.025 * replicates)], means[int(0.975 * replicates) - 1]


def compare(candidate, baseline):
    """Paired per-vector deltas (candidate - baseline) for every shared metric, by mode and partition."""
    base = {(result["mode"], result["name"]): flatten(result) for result in baseline["vectors"]}
    report = {}
    for result in candidate["vectors"]:
        before = base.get((result["mode"], result["name"]))
        if before is None:
            continue
        after = flatten(result)
        group = report.setdefault(f"{result['mode']}/{result['partition']}", {})
        for key, value in after.items():
            if value is not None and before.get(key) is not None:
                group.setdefault(key, []).append((result["name"], value - before[key]))
    summary = {}
    for group, metrics in report.items():
        summary[group] = {}
        for key, items in metrics.items():
            deltas = [delta for _, delta in items]
            summary[group][key] = {"mean": mean(deltas), "ci95": bootstrap_interval(deltas),
                                   "worst": max(items, key=lambda item: abs(item[1])), "n": len(items)}
    return summary


def check_gates(gates, candidate_summary, comparison):
    """Evaluate a predeclared gate set; returns (passed, lines)."""
    lines = []
    passed = True
    for check in gates["checks"]:
        modes = check["mode"] if isinstance(check["mode"], list) else [check["mode"]]
        for mode in modes:
            group = f"{mode}/{gates['partition']}"
            metric = check["metric"]
            ok, detail = evaluate_check(check, metric, candidate_summary.get(group, {}).get("mean", {}),
                                        comparison.get(group, {}))
            passed &= ok
            lines.append(f"{'PASS' if ok else 'FAIL'} {group} {metric} {check['rule']}: {detail}")
    return passed, lines


def evaluate_check(check, metric, means, deltas):
    rule = check["rule"]
    if rule == "max":
        value = means.get(metric)
        return value is not None and value <= check["value"], f"value={value}"
    if rule == "min":
        value = means.get(metric)
        return value is not None and value >= check["value"], f"value={value}"
    entry = deltas.get(metric)
    if entry is None:
        return False, "no paired data"
    delta, interval = entry["mean"], entry["ci95"]
    detail = f"mean delta={delta:+.4g} ci95={interval} worst={entry['worst']}"
    if rule == "decrease":
        return interval is not None and interval[1] < 0 and delta <= -check.get("value", 0.0), detail
    if rule == "increase":
        return interval is not None and interval[0] > 0 and delta >= check.get("value", 0.0), detail
    if rule == "no_increase":
        return delta <= check.get("value", 0.0), detail
    if rule == "no_decrease":
        return delta >= -check.get("value", 0.0), detail
    if rule == "abs_delta_max":
        return abs(delta) <= check["value"] and abs(entry["worst"][1]) <= check.get("worst", math.inf), detail
    raise ValueError(f"unknown gate rule {rule!r}")


def markdown(summary):
    lines = []
    for group, entry in summary.items():
        means = entry["mean"]
        lines.append(f"### {group} ({len(entry['vectors'])} vectors)\n")
        lines.append("| Metric | Ours vs input | DVSI vs input | Ours vs DVSI |")
        lines.append("|---|---|---|---|")
        for key in SPEECH_KEYS:
            cells = [means.get(f"{source}.{key}") for source in ("ours", "dvsi", "vs_dvsi")]
            lines.append(f"| {key} | " + " | ".join("" if v is None else f"{v:.3f}" for v in cells) + " |")
        lines.append("")
        lines.append("| f0 band (Hz) | Ours / coded | DVSI / coded | Input / coded |")
        lines.append("|---|---|---|---|")
        for low, high in F0_BANDS:
            key = f"{int(low)}_{int(high)}"
            cells = []
            for source in ("ours", "dvsi", "input"):
                pooled = entry[f"pitch_{source}_by_f0"][key]
                cells.append("" if pooled["median"] is None else f"{pooled['median']:.4f} (n={pooled['n']})")
            lines.append(f"| {key.replace('_', '–')} | " + " | ".join(cells) + " |")
        encoder_keys = sorted(key for key in means if key.startswith("encoder."))
        if encoder_keys:
            lines.append("")
            lines.append("| Encoder metric | Mean |")
            lines.append("|---|---|")
            for key in encoder_keys:
                value = means[key]
                lines.append(f"| {key[8:]} | {'' if value is None else f'{value:.4f}'} |")
        lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Score this library against DVSI's AMBE-3000 test vectors.")
    parser.add_argument("--vectors", type=Path, required=True, help="directory holding tv-rc/ from fetch_dvsi_vectors")
    parser.add_argument("--tools", type=Path, required=True, help="build directory with mbe_quality_* tools")
    parser.add_argument("--out", type=Path, required=True, help="work and report directory")
    parser.add_argument("--modes", default=",".join(MODES))
    parser.add_argument("--partitions", default="development,validation",
                        help="comma list; include 'final' only once per change set")
    parser.add_argument("--no-encoder", action="store_true", help="skip the D-STAR encoder comparison")
    parser.add_argument("--compare", type=Path, help="baseline scoreboard.json for paired deltas")
    parser.add_argument("--gates", help="gate set from dvsi_gates.json to evaluate against --compare")
    args = parser.parse_args()

    if importlib.util.find_spec("numpy") is None:
        sys.exit("dvsi_scoreboard.py needs NumPy")
    vectors = args.vectors / "tv-rc" if (args.vectors / "tv-rc").is_dir() else args.vectors
    args.out.mkdir(parents=True, exist_ok=True)
    work = args.out / "work"
    work.mkdir(exist_ok=True)
    results = []
    for mode in args.modes.split(","):
        if mode not in MODES:
            sys.exit(f"unknown mode {mode!r}")
        for partition in args.partitions.split(","):
            for name in PARTITIONS[partition]:
                needed = (vectors / mode / f"{name}.bit", vectors / mode / f"{name}.pcm", vectors / f"{name}.pcm")
                if not all(path.is_file() for path in needed):
                    print(f"skip {mode}/{name}: not fetched", file=sys.stderr)
                    continue
                results.append(score_vector(args.tools, vectors, mode, name, work, not args.no_encoder))
                print(f"scored {mode}/{name}", file=sys.stderr)
    summary = aggregate(results)
    report = {"schema_version": 1, "vectors": results, "summary": summary}
    if args.compare:
        report["comparison"] = compare(report, load_json(args.compare))
    (args.out / "scoreboard.json").write_text(json.dumps(report, indent=1) + "\n")
    (args.out / "scoreboard.md").write_text(markdown(summary) + "\n")
    print(markdown(summary))
    if args.gates:
        if not args.compare:
            sys.exit("--gates needs --compare")
        gates = load_json(Path(__file__).with_name("dvsi_gates.json"))[args.gates]
        passed, lines = check_gates(gates, summary, report["comparison"])
        print("\n".join(lines))
        print(f"gate set {args.gates}: {'PASS' if passed else 'FAIL'}")
        return 0 if passed else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
