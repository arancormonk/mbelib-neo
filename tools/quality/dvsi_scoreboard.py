#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Score this library against DVSI's AMBE-3000 test vectors.

For every speech vector of the D-STAR, P25 and AMBE+2 (rate 33) sets this decodes
DVSI's bits with mbe_quality_eval and measures, against the original input speech,
both our output and DVSI's own decoded output; it also compares our output with
DVSI's directly. A harmonic-peak pitch estimator checks that our decoded pitch
matches DVSI's and the input's on the same frames. It also encodes the input with
this library's encoder for each mode and compares pitch, voicing and spectral error
with DVSI's encoding, and for the rate-33 and D-STAR tone vectors compares the
encoder's tone frames with DVSI's.

The vectors are fetched locally with fetch_dvsi_vectors.py and must never be
committed or redistributed. This tool writes only aggregate numbers and per-run
work files under --out (git-ignored by default).

Needs NumPy for the pitch measurements.
"""

import argparse
import functools
import hashlib
import importlib.util
import json
import math
import random
import re
import shutil
import statistics
import struct
import subprocess
import sys
from pathlib import Path

MODES = {"dstar": "ambe2400", "p25": "imbe7200", "r33": "ambe2450"}
# Partitions by source recording (RECORDINGS): vectors made from one recording sit in one
# partition. Development and validation were both examined while earlier changes were
# made, so they are retrospective; "final" vectors are scored once per change set and
# never used for tuning. Not scored as speech: xfer (not speech), the *_dtx variants and
# the synthetic ptrain*/sine*/zero/grbge.
PARTITIONS = {
    "development": (
        "dam", "dam10", "dam20", "dam40", "dam80", "dam_10ov", "damcar", "lvl",
        "clean", "noisy", "fambf22c", "fambm22a", "t01", "t02", "tia11", "tambf22a", "ucarm15a",
        "ucarf15a", "ucarm20a", "p01mirs",
    ),
    "validation": ("t03", "t04", "tambf22b", "tambf32b", "tambf32e"),
    "final": ("irstia", "mark"),
    # Not speech: scored with the tone detector against DVSI's decoder only.
    "tones": tuple(
        ["alert", "alltone", "cp0", "cp1", "cp2", "cp31", "cpvbad"]
        + ["dtmf", "dtmf4", "dtmf4r", "dtmf8", "dtmf8r", "dtmf15n", "dtmf15p", "dtmf25", "dtmf50", "dtmf60"]
        + ["dtmf100", "dtmf180", "dtmfvbad"]
        + [f"dtone_{i}" for i in range(1, 17)]
        + [f"knox_{i}" for i in range(1, 17)]
    ),
    # dam's bits with channel errors (about 1.1, 2.2, 4.9 and 10.6% bit errors), as hard decisions and
    # as 4-bit soft decisions: our muting and output against DVSI's decoder on the same bits.
    "channel": tuple(f"dam_e{rate}_{kind}" for rate in (1, 2, 5, 10) for kind in ("hd", "sd")),
}
SPEECH_PARTITIONS = ("development", "validation", "final")
# Vectors made from one source recording, found by matching 3 s blocks of the inputs at any lag (envelope
# correlation 0.85 or more and waveform coherence 0.8 or more over 300-3400 Hz): fambm22a, tambf22a,
# ucarm15a, ucarf15a and ucarm20a are one recording through different filters and car noise, tambf22b and
# tambf32b one recording 1680 samples apart, and p01mirs is tia11 through another filter. No other pair of
# speech vectors shares a block. dam's variants (level, overload, noise and level steps in lvl) follow dam's
# envelope at lag 0, and its channel-error vectors are its bits.
RECORDINGS = {
    "dam": ("dam", "dam10", "dam20", "dam40", "dam80", "dam_10ov", "damcar", "lvl") + PARTITIONS["channel"],
    "fambm22a": ("fambm22a", "tambf22a", "ucarm15a", "ucarf15a", "ucarm20a"),
    "tambf22b": ("tambf22b", "tambf32b"),
    "tia11": ("tia11", "p01mirs"),
}
BAND_KEYS = tuple(
    f"band_delta_db_{band}"
    for band in ("0_250", "0_500", "500_1000", "1000_2000", "2000_3000", "3000_4000", "3500_4000")
)
SPEECH_KEYS = ("level_offset_db", "lsd_db", "env_corr", "crest_delta_db", "pcm_rail_samples") + BAND_KEYS
# DVSI's encoders delay the speech more than ours (about 200 samples more in D-STAR and rate 33, 125 in
# P25), past the evaluator's -160-sample lag search, so our decoded encoding is delayed by this much
# before DVSI's is used as its reference.
ENCODER_COMPARE_DELAY = {"dstar": 320, "r33": 320, "p25": 240}
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


def utterance_group(name):
    """The recording a vector was made from (RECORDINGS), else the vector itself."""
    for group, names in RECORDINGS.items():
        if name in names:
            return group
    return name


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


def quartiles(values):
    """(q1, q3) by linear interpolation between order statistics (statistics.quantiles needs 3.8)."""
    ordered = sorted(values)

    def at(fraction):
        position = fraction * (len(ordered) - 1)
        low = int(math.floor(position))
        high = min(low + 1, len(ordered) - 1)
        return ordered[low] + (ordered[high] - ordered[low]) * (position - low)

    return at(0.25), at(0.75)


def summarize_ratios(pairs):
    summary = {}
    for low, high in F0_BANDS:
        values = sorted(ratio for nominal, ratio in pairs if low <= nominal < high)
        key = f"{int(low)}_{int(high)}"
        if values:
            q1, q3 = quartiles(values)
            summary[key] = {"n": len(values), "median": statistics.median(values), "q1": q1, "q3": q3}
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
# Tones


def detect_tone(segment, nfft=8192, floor_db=-60.0):
    """{level_db, f1, f2, balance_db} for a window holding one or two steady sinusoids, else None.

    level_db is the window's RMS in dB relative to 32768. f1 < f2 for a dual tone,
    whose balance_db is the high component's power over the low one's.
    """
    import numpy as np

    x = np.asarray(segment, dtype=np.float64)
    rms = math.sqrt(float(np.mean(x * x))) if len(x) else 0.0
    if rms <= 32768.0 * 10 ** (floor_db / 20.0):
        return None
    power = np.abs(np.fft.rfft(x * np.hanning(len(x)), nfft)) ** 2
    hz = np.arange(len(power)) * SAMPLE_RATE / nfft
    total = float(np.sum(power))

    def peak(spectrum):
        k = int(np.argmax(spectrum))
        if 0 < k < len(spectrum) - 1 and spectrum[k - 1] > 0 and spectrum[k + 1] > 0:
            a, b, c = (math.log(float(spectrum[i])) for i in (k - 1, k, k + 1))
            curvature = a - 2.0 * b + c
            if curvature < 0:
                return (k + 0.5 * (a - c) / curvature) * SAMPLE_RATE / nfft
        return k * SAMPLE_RATE / nfft

    def energy(frequency):
        return float(np.sum(power[np.abs(hz - frequency) <= 40.0]))

    first = peak(power)
    # Past the first component's Hann main lobe (+-50 Hz for 320 samples) and its
    # energy window. Components closer than that (440 + 480 Hz ringback) read as one.
    rest = np.where(np.abs(hz - first) <= 100.0, 0.0, power)
    second = peak(rest)
    e1, e2 = energy(first), energy(second)
    level_db = 20.0 * math.log10(rms / 32768.0)
    if e2 >= 0.01 * total and e1 + e2 >= 0.9 * total:
        (low, e_low), (high, e_high) = sorted(((first, e1), (second, e2)))
        return {"level_db": level_db, "f1": low, "f2": high, "balance_db": 10.0 * math.log10(e_high / e_low)}
    if e1 >= 0.9 * total:
        return {"level_db": level_db, "f1": first, "f2": None, "balance_db": None}
    return None


def same_tone(a, b, tolerance=0.01):
    if a is None or b is None:
        return a is None and b is None
    if (a["f2"] is None) != (b["f2"] is None):
        return False
    pairs = [(a["f1"], b["f1"])] + ([(a["f2"], b["f2"])] if a["f2"] is not None else [])
    return all(abs(x / y - 1.0) <= tolerance for x, y in pairs)


def same_tone_identity(a, b):
    """Whether two detected tones are the same tone ID despite DVSI's approximations.

    Single tones sit 31.25 Hz apart, so they must agree within 10 Hz. DVSI plays a
    dual tone as two harmonics of one fundamental, up to 20 Hz off each nominal
    component, while neighbouring DTMF frequencies are at least 73 Hz apart.
    """
    if (a["f2"] is None) != (b["f2"] is None):
        return False
    if a["f2"] is None:
        return abs(a["f1"] - b["f1"]) <= 10.0
    return abs(a["f1"] - b["f1"]) <= 25.0 and abs(a["f2"] - b["f2"]) <= 25.0


def compare_tones(ours, dvsi, shift, window=320, hop=160):
    """Tone agreement on DVSI's steady windows; `shift` is our delay relative to DVSI's output."""
    starts = range(max(0, shift), len(ours) - window + 1, hop)
    pairs = []
    for start in starts:
        other = start - shift
        if other + window > len(dvsi):
            break
        pairs.append((detect_tone(ours[start:start + window]), detect_tone(dvsi[other:other + window])))
    level, frequency, balance = [], [], []
    agree = checked = wrong = 0
    for index in range(1, len(pairs) - 1):
        ours_tone, dvsi_tone = pairs[index]
        if not (same_tone(pairs[index - 1][1], dvsi_tone) and same_tone(dvsi_tone, pairs[index + 1][1])):
            continue  # DVSI's output changes here: a transition, not a steady window
        checked += 1
        both = ours_tone is not None and dvsi_tone is not None
        agree += (ours_tone is None and dvsi_tone is None) or (both and same_tone_identity(ours_tone, dvsi_tone))
        if not both:
            continue
        wrong += not same_tone_identity(ours_tone, dvsi_tone)
        level.append(abs(ours_tone["level_db"] - dvsi_tone["level_db"]))
        errors = [abs(ours_tone["f1"] / dvsi_tone["f1"] - 1.0)]
        if ours_tone["f2"] is not None and dvsi_tone["f2"] is not None:
            errors.append(abs(ours_tone["f2"] / dvsi_tone["f2"] - 1.0))
            balance.append(abs(ours_tone["balance_db"] - dvsi_tone["balance_db"]))
        elif (ours_tone["f2"] is None) != (dvsi_tone["f2"] is None):
            errors.append(1.0)  # one plays a single tone, the other a dual tone
        frequency.append(max(errors))
    return {
        "level_abs_error_db": statistics.median(level) if level else None,
        "freq_rel_error": statistics.median(frequency) if frequency else None,
        "balance_abs_error_db": statistics.median(balance) if balance else None,
        "detector_agreement": agree / checked if checked else None,
        "wrong_tone_fraction": wrong / len(level) if level else None,
        "steady_tone_windows": len(level),
    }


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


def usage_codecs(usage):
    """The codecs a tool's usage line offers after --codec."""
    match = re.search(r"--codec (\S+)", usage)
    return frozenset(match.group(1).split("|")) if match else frozenset()


@functools.lru_cache(maxsize=None)
def encoder_codecs(tools):
    """Codecs the build's mbe_quality_encode can encode, so an older build still scores what it has."""
    result = subprocess.run([str(tools / "mbe_quality_encode")], capture_output=True, text=True, check=False)
    return usage_codecs(result.stderr + result.stdout)


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


def decoder_identity(tools):
    """SHA-256 over a build's mbe_quality_eval and the mbelib-neo libraries beside it (the decoder it runs):
    libmbe-neo* where GCC and Clang build them, mbe-neo*.dll for MSVC."""
    digest = hashlib.sha256()
    tools = Path(tools)
    evaluator = tools / "mbe_quality_eval"
    if not evaluator.is_file() and (tools / "mbe_quality_eval.exe").is_file():
        evaluator = tools / "mbe_quality_eval.exe"
    libraries = [p for pattern in ("libmbe-neo*", "mbe-neo*.dll") for p in tools.glob(pattern) if p.is_file()]
    for path in [evaluator] + sorted(libraries):
        digest.update(path.name.encode() + b"\0" + hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


def lag_option(lags, key):
    """The evaluator's --lag-samples for comparison `key` when `lags` fixes it; else it searches."""
    return ["--lag-samples", lags[key]] if lags and lags.get(key) is not None else []


def joint_lags(baseline, candidate):
    """Per vector both runs scored, one lag per comparison both measured: their common lag, or where their
    alignment searches disagree, the lag of the run that aligned with the higher correlation. Equal
    correlations, and keys with none recorded (the encoder-tone frame shift), keep the baseline's. A
    change can let a decoder align where the old one could not, as when a gain law that played onsets
    16 dB low pulled the old decoder's alignment off by up to about 100 samples, so neither run's lag is
    right by default."""
    def quality(result, key):
        value = (result.get("alignment") or {}).get(key)
        return value if isinstance(value, (int, float)) else -math.inf

    other = {(result["mode"], result["name"]): result for result in candidate["vectors"]}
    joint = {}
    for result in baseline["vectors"]:
        twin = other.get((result["mode"], result["name"]))
        if twin is None or not result.get("lags") or not twin.get("lags"):
            continue
        lags = {}
        for key in result["lags"].keys() & twin["lags"].keys():
            mine, theirs = result["lags"][key], twin["lags"][key]
            if mine == theirs or quality(result, key) >= quality(twin, key):
                lags[key] = mine
            else:
                lags[key] = theirs
        joint[(result["mode"], result["name"])] = lags
    return joint


def load_lag_source(paths):
    """Lags to score at: one scoreboard.json's own, or two runs' joint_lags() (baseline first)."""
    reports = [load_json(path) for path in paths]
    if len(reports) == 1:
        return {(result["mode"], result["name"]): result["lags"] for result in reports[0]["vectors"]
                if result.get("lags")}
    if len(reports) == 2:
        return joint_lags(*reports)
    raise ValueError("lags come from one scoreboard.json, or two to join")


def alignment_ok(metrics):
    """True when the evaluator found an alignment inside its lag search (not none, and not at a limit)."""
    corr = metrics.get("alignment_corr")
    return isinstance(corr, (int, float)) and math.isfinite(corr) and not metrics.get("alignment_at_limit")


def aligned(metrics, label):
    """The evaluator's metrics, unless it found no usable alignment (a silent decode, or a lag past its search)."""
    if not alignment_ok(metrics):
        raise RuntimeError(f"{label}: no usable alignment (lag {metrics.get('auto_lag_samples')}, "
                           f"correlation {metrics.get('alignment_corr')})")
    return metrics


def write_delayed(path, samples, delay):
    """Raw s16le PCM of `samples` after `delay` zero samples."""
    values = [0] * delay + [int(sample) for sample in samples]
    Path(path).write_bytes(struct.pack(f"<{len(values)}h", *values))


def load_records(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line]


def score_vector(tools, vectors, mode, name, work, encoder, decode_tools=None, lags=None):
    """All measurements of one vector in one mode.

    `lags` (a baseline's recorded lags) fixes each comparison's alignment, so that both runs measure the
    same samples; `decode_tools` decodes the bits of the encoder comparisons (default `tools`).
    """
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
         "--ref", speech, "--json", f"{base}_ours.json", "--params", f"{base}_ours.jsonl"] + lag_option(lags, "ours"))
    run([tools / "mbe_quality_eval", "--decoded", dvsi_pcm, "--ref", speech, "--json", f"{base}_dvsi.json"]
        + lag_option(lags, "dvsi"))
    run([tools / "mbe_quality_eval", "--decoded", f"{base}_ours.wav", "--ref", dvsi_pcm,
         "--json", f"{base}_vs_dvsi.json"] + lag_option(lags, "vs_dvsi"))
    ours, dvsi, versus = (aligned(load_json(f"{base}_{tag}.json"), f"{base.name}_{tag}")
                          for tag in ("ours", "dvsi", "vs_dvsi"))
    records = load_records(f"{base}_ours.jsonl")
    ours_pcm, dvsi_samples, speech_samples = read_pcm(f"{base}_ours.wav"), read_pcm(dvsi_pcm), read_pcm(speech)
    result = {
        "mode": mode, "name": name, "partition": partition_of(name),
        "ours": {key: ours.get(key) for key in SPEECH_KEYS + ("pcm_fnv1a",)},
        "dvsi": {key: dvsi.get(key) for key in SPEECH_KEYS},
        "vs_dvsi": {key: versus.get(key) for key in SPEECH_KEYS},
        "lags": {"ours": ours["lag_samples"], "dvsi": dvsi["lag_samples"], "vs_dvsi": versus["lag_samples"]},
        "alignment": {"ours": ours["alignment_corr"], "dvsi": dvsi["alignment_corr"],
                      "vs_dvsi": versus["alignment_corr"]},
    }
    pairs = {
        "ours": pitch_ratios(records, ours_pcm, 0),
        "dvsi": pitch_ratios(records, dvsi_samples, versus["lag_samples"]),
        "input": pitch_ratios(records, speech_samples, ours["lag_samples"]),
    }
    result["pitch"] = {source: summarize_ratios(values) for source, values in pairs.items()}
    result["pitch_pairs"] = pairs
    if encoder and codec in encoder_codecs(tools):
        result["encoder"] = score_encoder(tools, decode_tools or tools, mode, speech, base, speech_samples,
                                          lags, result["lags"], result["alignment"])
    return result


def tone_frame_rate(records):
    """Fraction of frames decoded as tone or silence frames rather than voice."""
    tone_or_silence = 0x0010 | 0x0100  # MBE_PROCESS_FLAG_TONE, MBE_PROCESS_FLAG_SILENCE
    return sum(1 for record in records if record["flags"] & tone_or_silence) / len(records) if records else None


def score_encoder(tools, decode_tools, mode, speech, base, speech_samples, lags, measured, quality):
    """Our encoding of the input (by `tools`) against DVSI's, both decoded by `decode_tools`.

    One decoder for a baseline's and a candidate's encodings leaves only the encoders' difference. The
    lags this uses (fixed by `lags` when given) are added to `measured`, and the correlation each
    comparison's own alignment search reached to `quality`.
    """
    codec = MODES[mode]
    run([decode_tools / "mbe_quality_eval", "--codec", codec, "--frames", f"{base}.frames", "--out", f"{base}_ref.wav",
         "--ref", speech, "--json", f"{base}_ref.json", "--params", f"{base}_ref.jsonl"] + lag_option(lags, "ref"))
    dvsi_eval = aligned(load_json(f"{base}_ref.json"), f"{base.name}_ref")
    dvsi_records = load_records(f"{base}_ref.jsonl")
    measured["ref"], quality["ref"] = dvsi_eval["lag_samples"], dvsi_eval["alignment_corr"]
    run([tools / "mbe_quality_encode", "--codec", codec, "--in", speech, "--out", f"{base}_enc.rows"])
    run([decode_tools / "mbe_quality_eval", "--codec", codec, "--frames", f"{base}_enc.rows",
         "--out", f"{base}_enc.wav", "--ref", speech, "--json", f"{base}_enc.json",
         "--params", f"{base}_enc.jsonl"] + lag_option(lags, "enc"))
    encoded = load_json(f"{base}_enc.json")
    measured["enc"], quality["enc"] = encoded.get("lag_samples"), encoded.get("alignment_corr")
    if not alignment_ok(encoded):
        # Our encoding cannot be aligned with the input (say, it is mostly silence frames), so nothing
        # measured through that alignment is scored. The count still reaches the summary and the gates.
        print(f"warning: {base.name}_enc: no usable alignment with the input; not scored", file=sys.stderr)
        return {"alignment_failures": [f"{base.name}_enc"],
                "rows_sha256": hashlib.sha256(Path(f"{base}_enc.rows").read_bytes()).hexdigest()}
    # Both encodings through the same decoder: what remains is the encoders' difference, measured with
    # DVSI's encoding as the reference (its active frames, its level).
    delayed = Path(f"{base}_enc_delayed.raw")
    write_delayed(delayed, read_pcm(f"{base}_enc.wav"), ENCODER_COMPARE_DELAY[mode])
    run([decode_tools / "mbe_quality_eval", "--decoded", delayed, "--ref", f"{base}_ref.wav",
         "--json", f"{base}_enc_vs_dvsi.json"] + lag_option(lags, "enc_vs_dvsi"))
    failures = []
    versus = load_json(f"{base}_enc_vs_dvsi.json")
    measured["enc_vs_dvsi"], quality["enc_vs_dvsi"] = versus.get("lag_samples"), versus.get("alignment_corr")
    if not alignment_ok(versus):
        print(f"warning: {base.name}_enc_vs_dvsi: no usable alignment; not scored", file=sys.stderr)
        failures.append(f"{base.name}_enc_vs_dvsi")
        versus = {}
    records = load_records(f"{base}_enc.jsonl")

    def cents(pairs):
        values = [abs(1200.0 * math.log2(ratio)) for _, ratio in pairs]
        return statistics.median(values) if values else None

    ours_pairs = pitch_ratios(records, speech_samples, encoded["lag_samples"])
    dvsi_pairs = pitch_ratios(dvsi_records, speech_samples, dvsi_eval["lag_samples"])
    frame_offset = round((encoded["lag_samples"] - dvsi_eval["lag_samples"]) / FRAME)
    return {
        "rows_sha256": hashlib.sha256(Path(f"{base}_enc.rows").read_bytes()).hexdigest(),
        "speech": {key: encoded.get(key) for key in SPEECH_KEYS},
        "vs_dvsi_bits": {key: versus.get(key) for key in SPEECH_KEYS},
        "rail_excess_samples": max(0, (encoded.get("pcm_rail_samples") or 0) - (dvsi_eval.get("pcm_rail_samples") or 0)),
        "silence_frame_rate": tone_frame_rate(records),
        "dvsi_silence_frame_rate": tone_frame_rate(dvsi_records),
        "alignment_failures": failures,
        "pitch_cents_median": cents(ours_pairs),
        "dvsi_pitch_cents_median": cents(dvsi_pairs),
        "pitch_input_over_coded": summarize_ratios(ours_pairs),
        "dvsi_pitch_input_over_coded": summarize_ratios(dvsi_pairs),
        "octave_error_rate": octave_errors(records, speech_samples, encoded["lag_samples"]),
        "dvsi_octave_error_rate": octave_errors(dvsi_records, speech_samples, dvsi_eval["lag_samples"]),
        "voicing_agreement": voicing_agreement(records, dvsi_records, frame_offset),
        "frame_offset": frame_offset,
    }


def tone_fields(bits):
    """(tone index, AD) of an AMBE+2 tone frame's parameter bits (TIA-102.BABA-1 Table 10), else None."""
    if len(bits) != 49 or bits[:6] != "111111":
        return None
    return int(bits[12:20], 2), (int(bits[6:12], 2) << 1) | int(bits[44])


DSTAR_TONE_HIGH_BITS = (4, 0, 1, 2, 3, 7, 6, 5)  # the index's top three bits, by the selector in bits 6..8
DSTAR_TONE_VOLUME_BITS = (12, 13, 14, 15, 16, 44, 45, 17)  # most significant first


def dstar_tone_fields(bits):
    """(tone index, volume) of a D-STAR tone frame's parameter bits (b0 126 or 127), else None."""
    if len(bits) != 49 or bits[:6] != "111111":
        return None
    low = int(bits[9] + bits[42] + bits[43] + bits[10] + bits[11], 2)
    return (DSTAR_TONE_HIGH_BITS[int(bits[6:9], 2)] << 5) | low, int("".join(bits[i] for i in DSTAR_TONE_VOLUME_BITS), 2)


# DVSI's encoders run about a frame behind ours (on every speech vector our frame f lines up with DVSI's
# f + 1), so where several frame shifts match DVSI's tone frames equally well, the tone comparison takes
# the one nearest this. A tone that starts with the vector otherwise shows a tone frame where DVSI sends
# voice at the vector's first frame.
ENCODER_TONE_SHIFT = -1
# Per mode with tone frames from our encoder: how to read one, and the name of its level-step error.
ENCODER_TONE_FIELDS = {
    "r33": (tone_fields, "encoder_level_error_ad"),
    "dstar": (dstar_tone_fields, "encoder_level_error_volume"),
}


def encoder_tone_agreement(encoded, reference, max_shift=3, mode="r33", shift=None):
    """Our encoder's tone frames against DVSI's, at the frame shift where the most agree (of equals, the
    one nearest ENCODER_TONE_SHIFT), or at `shift` when given (a baseline's).

    Returns the fraction of DVSI's tone frames we send as the same tone, the fraction of all compared
    frames we send as a tone where DVSI sends voice, the mean absolute level difference, in the
    mode's level steps (AD, or D-STAR's volume), where the tones agree, and the shift used. When DVSI
    sends no tone frame only the second is defined.
    """
    fields, level_key = ENCODER_TONE_FIELDS[mode]
    ours = [fields(record.get("bits", "")) for record in encoded]
    theirs = [fields(record.get("bits", "")) for record in reference]
    tone_frames = [index for index, value in enumerate(theirs) if value is not None]
    if not tone_frames:
        compared = min(len(ours), len(theirs))
        return {"encoder_extra_rate": (sum(1 for value in ours[:compared] if value) / compared) if compared else None}

    def mine(index, shift):
        return ours[index + shift] if 0 <= index + shift < len(ours) else None

    def agreeing(shift):
        return [index for index in tone_frames if mine(index, shift) and mine(index, shift)[0] == theirs[index][0]]

    if shift is None:
        shift = max(range(-max_shift, max_shift + 1),
                    key=lambda candidate: (len(agreeing(candidate)), -abs(candidate - ENCODER_TONE_SHIFT)))
    compared = [index for index in range(len(theirs)) if 0 <= index + shift < len(ours)]
    extra = [index for index in compared if theirs[index] is None and mine(index, shift)]
    same = agreeing(shift)
    return {
        "encoder_agreement": len(same) / len(tone_frames),
        "encoder_extra_rate": len(extra) / len(compared) if compared else 0.0,
        level_key: mean([abs(mine(index, shift)[1] - theirs[index][1]) for index in same]),
        "shift": shift,
    }


def score_tone_vector(tools, vectors, mode, name, work, encoder=False, decode_tools=None, lags=None):
    """Our decode of DVSI's bits for a tone vector against DVSI's decoded output, and for rate 33 and
    D-STAR our encoding of the input against DVSI's (decoded by `decode_tools`, default `tools`)."""
    codec = MODES[mode]
    base = work / f"{mode}_{name}"
    dvsi_pcm = Path(f"{base}_dvsi.raw")
    shutil.copyfile(vectors / mode / f"{name}.pcm", dvsi_pcm)
    run([tools / "mbe_quality_reframe", "--codec", codec, "--from-dvsi", vectors / mode / f"{name}.bit",
         "--out", f"{base}.frames"])
    run([tools / "mbe_quality_eval", "--codec", codec, "--frames", f"{base}.frames", "--out", f"{base}_ours.wav",
         "--params", f"{base}_ours.jsonl"])
    run([tools / "mbe_quality_eval", "--decoded", f"{base}_ours.wav", "--ref", dvsi_pcm,
         "--json", f"{base}_vs_dvsi.json"] + lag_option(lags, "vs_dvsi"))
    versus = aligned(load_json(f"{base}_vs_dvsi.json"), f"{base.name}_vs_dvsi")
    records = load_records(f"{base}_ours.jsonl")
    metrics = compare_tones(read_pcm(f"{base}_ours.wav"), read_pcm(dvsi_pcm), versus.get("lag_samples") or 0)
    metrics["tone_frames"] = sum(1 for record in records if record["flags"] & 0x0010)
    metrics["frames"] = len(records)
    metrics["alignment_corr"] = versus.get("alignment_corr")
    used = {"vs_dvsi": versus["lag_samples"]}
    if encoder and mode in ENCODER_TONE_FIELDS and codec in encoder_codecs(tools):
        speech = work / f"{name}_input.raw"
        shutil.copyfile(vectors / f"{name}.pcm", speech)
        run([tools / "mbe_quality_encode", "--codec", codec, "--in", speech, "--out", f"{base}_enc.rows"])
        run([(decode_tools or tools) / "mbe_quality_eval", "--codec", codec, "--frames", f"{base}_enc.rows",
             "--out", f"{base}_enc.wav", "--params", f"{base}_enc.jsonl"])
        agreement = encoder_tone_agreement(load_records(f"{base}_enc.jsonl"), records, mode=mode,
                                           shift=(lags or {}).get("encoder_tone_shift"))
        if "shift" in agreement:
            used["encoder_tone_shift"] = agreement.pop("shift")
        metrics.update(agreement)
    return {"mode": mode, "name": name, "partition": "tones", "tone": metrics, "lags": used,
            "alignment": {"vs_dvsi": versus["alignment_corr"]}}


# --------------------------------------------------------------------------------------------
# Channel errors

LOUD_FRAME_ENERGY = FRAME * 300.0 ** 2  # a frame of DVSI's clean decode at RMS 300 or more
MUTED_ENERGY_RATIO = 1e-3  # DVSI's decode of the errored bits 30 dB or more below its clean decode
MUTE_FLAG, REPEAT_FLAG = 0x0080, 0x0040  # MBE_PROCESS_FLAG_MUTE, MBE_PROCESS_FLAG_REPEAT


def channel_frames(clean, errored):
    """Per frame of DVSI's output, (loud, muted): its clean decode is loud, its decode of the errored bits
    is muted (meaningful on loud frames only)."""
    import numpy as np

    clean, errored = np.asarray(clean, dtype=np.float64), np.asarray(errored, dtype=np.float64)
    loud, muted = [], []
    for f in range(min(len(clean), len(errored)) // FRAME):
        energy = float(np.sum(clean[f * FRAME:(f + 1) * FRAME] ** 2))
        loud.append(energy >= LOUD_FRAME_ENERGY)
        muted.append(float(np.sum(errored[f * FRAME:(f + 1) * FRAME] ** 2)) < MUTED_ENERGY_RATIO * max(energy, 1.0))
    return loud, muted


def mute_metrics(loud, dvsi_muted, muted):
    """Our mute decisions against DVSI's on loud frames, frame for frame (frame f of both outputs comes from
    the same bits)."""
    frames = [f for f in range(min(len(loud), len(dvsi_muted), len(muted))) if loud[f]]

    def share(test):
        return sum(1 for f in frames if test(f)) / len(frames) if frames else None

    return {
        "loud_frames": len(frames),
        "mute_agreement": share(lambda f: bool(muted[f]) == bool(dvsi_muted[f])),
        "false_mute_rate": share(lambda f: muted[f] and not dvsi_muted[f]),
        "missed_mute_rate": share(lambda f: dvsi_muted[f] and not muted[f]),
        "mute_rate": share(lambda f: muted[f]),
        "dvsi_mute_rate": share(lambda f: dvsi_muted[f]),
        "first_mute_frame": next((f for f in frames if muted[f]), None),
        "dvsi_first_mute_frame": next((f for f in frames if dvsi_muted[f]), None),
    }


def frame_lsd(decoded, reference, lag, gain):
    """Per frame of `reference`, the mean log-spectral distance (dB) of `decoded` from it.

    Reference sample n is compared with decoded sample n + lag scaled by `gain`, on mbe_quality_eval's
    spectra: 256-sample Hann windows every 80 samples, bins 2 to 118, floored 60 dB below the reference's
    strongest bin, each counted in the frame holding its centre. None where no spectrum is centred.
    """
    import numpy as np

    d = np.asarray(decoded, dtype=np.float64) * gain
    r = np.asarray(reference, dtype=np.float64)
    count = len(r) // FRAME
    positions = [p for p in range(0, len(r) - 255, 80) if p + lag >= 0 and p + lag + 256 <= len(d)]
    if not positions:
        return [None] * count
    window = np.hanning(256)
    ref = np.abs(np.fft.rfft(np.stack([r[p:p + 256] for p in positions]) * window, axis=1)) ** 2
    dec = np.abs(np.fft.rfft(np.stack([d[p + lag:p + lag + 256] for p in positions]) * window, axis=1)) ** 2
    floor = max(float(ref.max()) * 1e-6, 1e-20)
    delta = 10.0 * np.log10(np.maximum(ref[:, 2:119] + floor, 1e-20) / np.maximum(dec[:, 2:119] + floor, 1e-20))
    sums, counts = [0.0] * count, [0] * count
    for position, distance in zip(positions, np.sqrt(np.mean(delta * delta, axis=1))):
        frame = (position + 128) // FRAME
        if frame < count:
            sums[frame] += float(distance)
            counts[frame] += 1
    return [total / n if n else None for total, n in zip(sums, counts)]


def covered_unmuted(muted, lag, count):
    """Per reference frame g, whether all our samples under it (from 160 g + lag) come from frames we did
    not mute."""
    cover = []
    for g in range(count):
        first, last = (FRAME * g + lag) // FRAME, (FRAME * g + FRAME - 1 + lag) // FRAME
        cover.append(first >= 0 and last < len(muted) and not any(muted[first:last + 1]))
    return cover


def channel_lsd(*runs):
    """(mean LSD of each run, coverage) on the frames that count in every run of one channel vector: loud
    in DVSI's clean decode, not muted by DVSI, measured by every run, and with no muted frame at or
    next to them (a frame's spectra reach into its neighbours): neither a loud frame DVSI mutes nor one
    under a frame a run muted. Coverage is their share of all the loud frames DVSI leaves unmuted, so
    frames a run did not measure count as lost; LSD is None without any such frame."""
    frames = runs[0]["frames"]
    count = len(frames["loud"])
    eligible = [frames["loud"][g] and not frames["dvsi_muted"][g] for g in range(count)]
    clear = [not (frames["loud"][g] and frames["dvsi_muted"][g]) for g in range(count)]
    for run in runs:
        cover = covered_unmuted(run["frames"]["muted"], run["lags"]["clean_vs_dvsi"], count)
        clear = [ok and covered for ok, covered in zip(clear, cover)]

    def measured(run, g):
        values = run["frames"]["lsd"]
        return g < len(values) and values[g] is not None

    chosen = [g for g in range(count) if eligible[g] and all(clear[max(0, g - 1):g + 2])
              and all(measured(run, g) for run in runs)]
    coverage = len(chosen) / sum(eligible) if any(eligible) else None
    if not chosen:
        return [None] * len(runs), coverage
    return [sum(run["frames"]["lsd"][g] for g in chosen) / len(chosen) for run in runs], coverage


@functools.lru_cache(maxsize=None)
def channel_reference(tools, vectors, mode, work, lag):
    """(lag, gain, alignment correlation) aligning and level-matching our decodes of dam's bits with
    DVSI's, from the clean bits; `lag` fixes the lag when not None."""
    codec = MODES[mode]
    base = work / f"{mode}_channel_dam"
    dvsi_pcm = Path(f"{base}_dvsi.raw")
    shutil.copyfile(vectors / mode / "dam.pcm", dvsi_pcm)
    run([tools / "mbe_quality_reframe", "--codec", codec, "--from-dvsi", vectors / mode / "dam.bit",
         "--out", f"{base}.frames"])
    run([tools / "mbe_quality_eval", "--codec", codec, "--frames", f"{base}.frames", "--out", f"{base}_ours.wav"])
    run([tools / "mbe_quality_eval", "--decoded", f"{base}_ours.wav", "--ref", dvsi_pcm,
         "--json", f"{base}_vs_dvsi.json"] + lag_option({"clean": lag}, "clean"))
    versus = aligned(load_json(f"{base}_vs_dvsi.json"), f"{base.name}_vs_dvsi")
    return versus["lag_samples"], 10.0 ** (-versus["level_offset_db"] / 20.0), versus["alignment_corr"]


def score_channel_vector(tools, vectors, mode, name, work, lags=None):
    """Our decode of one of dam's channel-error vectors against DVSI's decode of the same bits."""
    codec = MODES[mode]
    kind = name.rsplit("_", 1)[1]
    lag, gain, corr = channel_reference(tools, vectors, mode, work, (lags or {}).get("clean_vs_dvsi"))
    base = work / f"{mode}_{name}"
    run([tools / "mbe_quality_reframe", "--codec", codec, "--from-dvsi", vectors / mode / f"{name}.bit",
         "--out", f"{base}.frames"])
    run([tools / "mbe_quality_eval", "--codec", codec, "--frames", f"{base}.frames", "--out", f"{base}_ours.wav",
         "--params", f"{base}_ours.jsonl"] + (["--soft"] if kind == "sd" else []))
    records = load_records(f"{base}_ours.jsonl")
    dvsi_errored = read_pcm(vectors / mode / f"{name}.pcm")
    loud, dvsi_muted = channel_frames(read_pcm(vectors / mode / "dam.pcm"), dvsi_errored)
    muted = [bool(record["flags"] & MUTE_FLAG) for record in records]
    metrics = mute_metrics(loud, dvsi_muted, muted)
    metrics["repeat_rate"] = sum(1 for r in records if r["flags"] & REPEAT_FLAG) / len(records) if records else None
    metrics["errors_per_frame"] = mean([record["errors"] for record in records])
    return {"mode": mode, "name": name, "partition": "channel", "kind": kind, "channel": metrics,
            "lags": {"clean_vs_dvsi": lag}, "alignment": {"clean_vs_dvsi": corr},
            "frames": {"loud": loud, "dvsi_muted": dvsi_muted, "muted": muted,
                       "lsd": frame_lsd(read_pcm(f"{base}_ours.wav"), dvsi_errored, lag, gain)}}


# --------------------------------------------------------------------------------------------
# Aggregation, comparison and gates


def mean(values):
    values = [value for value in values if value is not None]
    return sum(values) / len(values) if values else None


def flatten(result):
    """Scalar metrics of one vector as dotted keys (the namespace the gates use)."""
    flat = {}
    if "tone" in result:
        return {f"tone.{key}": value for key, value in result["tone"].items()}
    if "channel" in result:
        prefix = f"channel.{result['kind']}."
        flat = {prefix + key: value for key, value in result["channel"].items()
                if isinstance(value, (int, float)) and not isinstance(value, bool)}
        (lsd,), coverage = channel_lsd(result)
        flat[prefix + "lsd_db"], flat[prefix + "lsd_coverage"] = lsd, coverage
        return flat
    for group in ("ours", "dvsi", "vs_dvsi"):
        for key, value in result[group].items():
            if isinstance(value, (int, float)):
                flat[f"{group}.{key}"] = value
    low = result["vs_dvsi"].get("band_delta_db_0_250")
    if isinstance(low, (int, float)):
        flat["vs_dvsi.lf_abs_error_db"] = abs(low)
    bands = [result["vs_dvsi"].get(f"band_delta_db_{band}") for band in ("2000_3000", "3000_4000")]
    if all(isinstance(value, (int, float)) for value in bands):
        flat["vs_dvsi.hf_abs_error_db"] = (abs(bands[0]) + abs(bands[1])) / 2.0
    for source, summary in result["pitch"].items():
        flat[f"pitch.{source}_ratio"] = summary["all"]["median"]
    encoder = result.get("encoder")
    if encoder:
        flat["encoder.alignment_failures"] = len(encoder.get("alignment_failures", []))
        for key, value in encoder.get("speech", {}).items():
            if isinstance(value, (int, float)):
                flat[f"encoder.{key}"] = value
        for key, value in encoder.get("vs_dvsi_bits", {}).items():
            if isinstance(value, (int, float)):
                flat[f"encoder.vs_dvsi_bits.{key}"] = value
        level = encoder.get("vs_dvsi_bits", {}).get("level_offset_db")
        if isinstance(level, (int, float)):
            flat["encoder.level_gap_db"] = abs(level)
        for key in ("pitch_cents_median", "dvsi_pitch_cents_median", "octave_error_rate", "dvsi_octave_error_rate",
                    "rail_excess_samples", "silence_frame_rate", "dvsi_silence_frame_rate"):
            if key in encoder:
                flat[f"encoder.{key}"] = encoder[key]
        agreement = [value for value in encoder.get("voicing_agreement", []) if value is not None]
        flat["encoder.voicing_agreement_mean"] = mean(agreement)
    return flat


def pooled_pitch(results, source):
    """Per-f0-band summary of every vector's frames pooled together."""
    pairs = [tuple(pair) for result in results for pair in result.get("pitch_pairs", {}).get(source, [])]
    return {key: entry for key, entry in summarize_ratios(pairs).items() if key != "all"}


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
            if partition == "channel":
                entry["rows"] = {result["name"]: flat for result, flat in zip(group, flats)}
                summary[f"{mode}/{partition}"] = entry
                continue
            if partition == "tones":
                # The worst tone vector decides: errors take the maximum, agreement the minimum. The
                # encoder's tone metrics also get their mean over vectors, since the vectors of
                # deliberately malformed tones are ones DVSI itself is inconsistent on.
                for key in keys:
                    values = [flat[key] for flat in flats if flat.get(key) is not None]
                    pick = min if key.endswith("agreement") else max
                    entry["mean"][key] = pick(values) if values else None
                    if key.startswith("tone.encoder_"):
                        entry["mean"][f"{key}_mean"] = mean(values)
                summary[f"{mode}/{partition}"] = entry
                continue
            for source in ("ours", "dvsi", "input"):
                pooled = pooled_pitch(group, source)
                entry[f"pitch_{source}_by_f0"] = pooled
                entry["mean"][f"pitch.{source}_max_band_deviation"] = max_ratio_deviation(pooled)
            summary[f"{mode}/{partition}"] = entry
    return summary


def pcm_changes(candidate, baseline):
    """Per mode, the baseline's speech vectors whose decoded PCM differs in the candidate.

    A vector missing from the candidate, or without a hash on either side, counts
    as changed, so incomplete runs cannot pass an identity gate.
    """
    after = {(r["mode"], r["name"]): r["ours"].get("pcm_fnv1a") for r in candidate["vectors"] if "ours" in r}
    counts = {}
    for result in baseline["vectors"]:
        if "ours" not in result:
            continue
        key = (result["mode"], result["name"])
        before, now = result["ours"].get("pcm_fnv1a"), after.get(key)
        counts[result["mode"]] = counts.get(result["mode"], 0) + (before is None or now is None or before != now)
    return counts


def encoder_row_changes(candidate, baseline):
    """Per mode, the baseline's encoded vectors whose encoder bits differ in the candidate.

    Only modes whose baseline recorded encoder-row hashes are counted; there a vector the candidate
    did not encode, or encoded without a hash, counts as changed.
    """
    after = {(r["mode"], r["name"]): (r.get("encoder") or {}).get("rows_sha256") for r in candidate["vectors"]}
    counts = {}
    for result in baseline["vectors"]:
        before = (result.get("encoder") or {}).get("rows_sha256")
        if before is None:
            continue
        now = after.get((result["mode"], result["name"]))
        counts[result["mode"]] = counts.get(result["mode"], 0) + (now is None or now != before)
    return counts


def add_encoder_row_changes(summary, candidate, baseline):
    """Expose encoder_row_changes() as encoder.rows_changed_vectors in every speech group of each mode."""
    for mode, count in encoder_row_changes(candidate, baseline).items():
        for group in [key for key in summary if key.startswith(f"{mode}/") and not key.endswith("/tones")]:
            summary[group]["mean"]["encoder.rows_changed_vectors"] = count


def add_pcm_changes(summary, candidate, baseline):
    """Expose pcm_changes() as speech.pcm_changed_vectors in every group of each mode."""
    for mode, count in pcm_changes(candidate, baseline).items():
        groups = [key for key in summary if key.startswith(f"{mode}/")] or [f"{mode}/tones"]
        for group in set(groups) | {f"{mode}/tones"}:
            summary.setdefault(group, {"vectors": [], "mean": {}})["mean"]["speech.pcm_changed_vectors"] = count


def bootstrap_interval(groups, replicates=2000, seed=0x5EED):
    """95% interval of the mean per-vector delta, resampling whole groups of deltas (the vectors of one
    recording, utterance_group()), since those do not vary independently."""
    if len(groups) < 2:
        return None
    rng = random.Random(seed)
    means = []
    for _ in range(replicates):
        sample = [delta for _ in groups for delta in rng.choice(groups)]
        means.append(sum(sample) / len(sample))
    means.sort()
    return means[int(0.025 * replicates)], means[int(0.975 * replicates) - 1]


def check_encoder_decoders(candidate, baseline):
    """Refuse encoder deltas unless one decoder (--decode-tools) decoded both runs' encodings."""
    if not any("encoder" in result for result in candidate["vectors"]):
        return
    if not any("encoder" in result for result in baseline["vectors"]):
        return
    before, after = baseline.get("decoder_identity"), candidate.get("decoder_identity")
    if before is None or before != after:
        raise ValueError(f"the encoder comparisons were decoded by different decoders (baseline {before}, "
                         f"this run {after}); score both runs with the same --decode-tools, or use --no-encoder")


def compare(candidate, baseline):
    """Paired per-vector deltas (candidate - baseline) for every shared metric, by mode and partition."""
    check_encoder_decoders(candidate, baseline)
    raw = {(result["mode"], result["name"]): result for result in baseline["vectors"]}
    base = {key: flatten(result) for key, result in raw.items()}
    report, unpaired, coverage = {}, {}, {}
    for result in candidate["vectors"]:
        before = base.get((result["mode"], result["name"]), {})
        after = flatten(result)
        name = f"{result['mode']}/{result['partition']}"
        group = report.setdefault(name, {})
        twin = raw.get((result["mode"], result["name"]))
        if candidate.get("lags_from") and twin is not None and result.get("lags") != twin.get("lags"):
            # Measured at other lags than the baseline's: not a pair.
            for key in sorted(set(after) | set(before)):
                unpaired.setdefault(name, {}).setdefault(key, []).append(result["name"])
            continue
        if "channel" in result and twin is not None and "channel" in twin:
            # LSD against DVSI on the frames that count in both runs, since each run mutes its own frames.
            key = f"channel.{result['kind']}.lsd_common_db"
            (value, previous), share = channel_lsd(result, twin)
            after[key], before = value, dict(before, **{key: previous})
            coverage.setdefault(name, {}).setdefault(key, []).append(share)
        for key in sorted(set(after) | set(before)):
            value, previous = after.get(key), before.get(key)
            if value is not None and previous is not None:
                group.setdefault(key, []).append((result["name"], value - previous))
            elif value is not None or previous is not None:
                # Measured in one run only (an alignment failure, a vector missing from the baseline).
                unpaired.setdefault(name, {}).setdefault(key, []).append(result["name"])
    seen = {(result["mode"], result["name"]) for result in candidate["vectors"]}
    for result in baseline["vectors"]:
        if (result["mode"], result["name"]) in seen:
            continue
        # A baseline vector the candidate run left out: every metric it has goes unpaired.
        name = f"{result['mode']}/{result['partition']}"
        report.setdefault(name, {})
        for key, value in sorted(flatten(result).items()):
            if value is not None:
                unpaired.setdefault(name, {}).setdefault(key, []).append(result["name"])
    summary = {}
    for group, metrics in report.items():
        summary[group] = {}
        for key in sorted(set(metrics) | set(unpaired.get(group, {}))):
            items = metrics.get(key, [])
            deltas = [delta for _, delta in items]
            recordings = {}
            for name, delta in items:
                recordings.setdefault(utterance_group(name), []).append(delta)
            summary[group][key] = {"mean": mean(deltas), "ci95": bootstrap_interval(list(recordings.values())),
                                   "worst": max(items, key=lambda item: abs(item[1])) if items else None,
                                   "min": min(deltas) if deltas else None, "max": max(deltas) if deltas else None,
                                   "n": len(items), "groups": len(recordings),
                                   "unpaired": unpaired.get(group, {}).get(key, [])}
            shares = coverage.get(group, {}).get(key)
            if shares is not None:
                summary[group][key]["coverage"] = None if None in shares else min(shares)
    return summary


def unscored_vectors(gates, *runs):
    """(mode, name) of the gate set's partition, for every mode it gates, that one of `runs` did not score:
    a gate measured on part of its partition says nothing about the rest."""
    modes = {mode for check in gates["checks"]
             for mode in (check["mode"] if isinstance(check["mode"], list) else [check["mode"]])}
    missing = set()
    for run in runs:
        scored = {(result["mode"], result["name"]) for result in run["vectors"]}
        missing |= {(mode, name) for mode in modes for name in PARTITIONS[gates["partition"]]
                    if (mode, name) not in scored}
    return sorted(missing)


def evaluate_gate_sets(names, declared, candidate_summary, comparison, runs):
    """Every named gate set of `declared`: its checks, every vector of its partition scored in each of
    `runs` (the candidate and the baseline), and the candidate scored at fixed lags, so that compare()
    paired only vectors measured at the baseline's lags. Returns (passed, lines)."""
    passed, lines = True, []
    for name in names:
        if name not in declared:
            raise ValueError(f"unknown gate set {name!r}")
        ok, set_lines = check_gates(declared[name], candidate_summary, comparison)
        missing = unscored_vectors(declared[name], *runs)
        if missing:
            ok = False
            set_lines.append(f"FAIL {declared[name]['partition']}: not scored in both runs: {missing}")
        if not runs[0].get("lags_from"):
            ok = False
            set_lines.append(f"FAIL {declared[name]['partition']}: not scored at fixed lags (--lags-from)")
        lines += set_lines + [f"gate set {name}: {'PASS' if ok else 'FAIL'}"]
        passed &= ok
    return passed, lines


def check_gates(gates, candidate_summary, comparison):
    """Evaluate a predeclared gate set; returns (passed, lines)."""
    lines = []
    passed = True
    # A partition with unscorable vectors cannot pass: their metrics are missing, not measured.
    gated = {mode for check in gates["checks"]
             for mode in (check["mode"] if isinstance(check["mode"], list) else [check["mode"]])}
    for mode in sorted(gated):
        group = f"{mode}/{gates['partition']}"
        for key, value in sorted(candidate_summary.get(group, {}).get("mean", {}).items()):
            if key.endswith("alignment_failures") and value:
                passed = False
                lines.append(f"FAIL {group} {key}: vectors that could not be aligned")
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
    if entry is None or entry["mean"] is None:
        return False, "no paired data"
    if entry.get("unpaired"):
        return False, f"unpaired vectors {entry['unpaired']}"
    delta, interval = entry["mean"], entry["ci95"]
    detail = f"mean delta={delta:+.4g} ci95={interval} worst={entry['worst']}"
    if "min_coverage" in check:
        # A delta on too few common frames (say, both runs mute most of a vector) does not count.
        share = entry.get("coverage")
        if share is None or share < check["min_coverage"]:
            return False, f"{detail} coverage={share} below {check['min_coverage']}"
        detail += f" coverage={share:.3f}"
    if rule == "all_no_increase":
        return entry["max"] <= check.get("value", 0.0), f"{detail} max={entry['max']:+.4g}"
    if rule == "all_no_decrease":
        return entry["min"] >= -check.get("value", 0.0), f"{detail} min={entry['min']:+.4g}"
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
        if group.endswith("/channel"):
            columns = ("mute_agreement", "false_mute_rate", "missed_mute_rate", "mute_rate", "dvsi_mute_rate",
                       "first_mute_frame", "dvsi_first_mute_frame", "repeat_rate", "lsd_db", "lsd_coverage")
            lines.append("| Vector | " + " | ".join(columns) + " |")
            lines.append("|---" * (len(columns) + 1) + "|")
            for name, flat in entry.get("rows", {}).items():
                prefix = f"channel.{name.rsplit('_', 1)[1]}."
                cells = [flat.get(prefix + column) for column in columns]
                lines.append(f"| {name} | " + " | ".join("" if v is None else f"{v:.4g}" for v in cells) + " |")
            lines.append("")
            continue
        if group.endswith("/tones"):
            lines.append("| Metric | Worst vector |")
            lines.append("|---|---|")
            for key in sorted(means):
                value = means[key]
                lines.append(f"| {key} | {'' if value is None else f'{value:.4f}'} |")
            lines.append("")
            continue
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
    from reference_provenance import retain_dataset

    parser = argparse.ArgumentParser(description="Score this library against DVSI's AMBE-3000 test vectors.")
    parser.add_argument("--vectors", type=Path, required=True, help="directory holding tv-rc/ from fetch_dvsi_vectors")
    parser.add_argument("--tools", type=Path, required=True, help="build directory with mbe_quality_* tools")
    parser.add_argument("--out", type=Path, required=True, help="work and report directory")
    parser.add_argument("--modes", default=",".join(MODES))
    parser.add_argument("--partitions", default="development,validation",
                        help="comma list; include 'final' only once per change set")
    parser.add_argument("--no-encoder", action="store_true", help="skip the encoder comparisons")
    parser.add_argument("--decode-tools", type=Path,
                        help="build whose decoder decodes both encodings in the encoder comparisons (default "
                             "--tools); give a baseline and a candidate run the same one")
    parser.add_argument("--compare", type=Path,
                        help="baseline scoreboard.json for paired deltas; its recorded lags fix this run's "
                             "alignments unless --lags-from is given")
    parser.add_argument("--lags-from",
                        help="scoreboard.json whose lags fix this run's alignments, or two (baseline,candidate; "
                             "comma-separated) to score at their joint lags")
    parser.add_argument("--gates", help="comma list of gate sets from dvsi_gates.json to evaluate against --compare")
    args = parser.parse_args()

    if importlib.util.find_spec("numpy") is None:
        sys.exit("dvsi_scoreboard.py needs NumPy")
    vectors = args.vectors / "tv-rc" if (args.vectors / "tv-rc").is_dir() else args.vectors
    args.out.mkdir(parents=True, exist_ok=True)
    acquisition_manifest = args.vectors / "dvsi-manifest.json"
    if not acquisition_manifest.exists() and args.vectors.name == "tv-rc":
        acquisition_manifest = args.vectors.parent / "dvsi-manifest.json"
    reference_data = retain_dataset(acquisition_manifest, args.out)
    work = args.out / "work"
    work.mkdir(exist_ok=True)
    baseline = load_json(args.compare) if args.compare else None
    # Both runs measure the same samples: each comparison at one lag for both.
    sources = [Path(path) for path in args.lags_from.split(",")] if args.lags_from else (
        [args.compare] if args.compare else [])
    try:
        frozen = load_lag_source(sources) if sources else {}
    except ValueError as error:
        sys.exit(str(error))
    if sources and not frozen:
        print("warning: no lags recorded there, so each comparison searches its own", file=sys.stderr)
    encoder = not args.no_encoder
    decode_tools = args.decode_tools or args.tools
    results = []
    for mode in args.modes.split(","):
        if mode not in MODES:
            sys.exit(f"unknown mode {mode!r}")
        for partition in args.partitions.split(","):
            for name in PARTITIONS[partition]:
                needed = [vectors / mode / f"{name}.bit", vectors / mode / f"{name}.pcm"]
                # A channel-error vector is measured against DVSI's decode of the clean bits.
                needed += ([vectors / mode / "dam.bit", vectors / mode / "dam.pcm"] if partition == "channel"
                           else [vectors / f"{name}.pcm"])
                if not all(path.is_file() for path in needed):
                    print(f"skip {mode}/{name}: not fetched", file=sys.stderr)
                    continue
                lags = frozen.get((mode, name))
                if partition == "channel":
                    results.append(score_channel_vector(args.tools, vectors, mode, name, work, lags))
                elif partition == "tones":
                    results.append(score_tone_vector(args.tools, vectors, mode, name, work, encoder, decode_tools, lags))
                else:
                    results.append(score_vector(args.tools, vectors, mode, name, work, encoder, decode_tools, lags))
                print(f"scored {mode}/{name}", file=sys.stderr)
    summary = aggregate(results)
    report = {"schema_version": 2, "tools_identity": decoder_identity(args.tools),
              "reference_data": reference_data,
              "decoder_identity": decoder_identity(decode_tools) if encoder else None,
              "lags_from": ",".join(map(str, sources)) if frozen else None, "vectors": results, "summary": summary}
    if baseline:
        try:
            report["comparison"] = compare(report, baseline)
        except ValueError as error:
            sys.exit(str(error))
        add_pcm_changes(summary, report, baseline)
        add_encoder_row_changes(summary, report, baseline)
    (args.out / "scoreboard.json").write_text(json.dumps(report, indent=1) + "\n")
    (args.out / "scoreboard.md").write_text(markdown(summary) + "\n")
    print(markdown(summary))
    if args.gates:
        if not args.compare:
            sys.exit("--gates needs --compare")
        declared = load_json(Path(__file__).with_name("dvsi_gates.json"))
        try:
            passed, lines = evaluate_gate_sets(args.gates.split(","), declared, summary, report["comparison"],
                                               (report, baseline))
        except ValueError as error:
            sys.exit(str(error))
        print("\n".join(lines))
        return 0 if passed else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
