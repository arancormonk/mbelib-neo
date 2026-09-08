#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Offline audio views and independent NumPy cross-checks, not replacement metrics."""

import importlib
import importlib.metadata
import math
import wave
from pathlib import Path

import numpy as np


def package_versions():
    try:
        importlib.import_module("numpy")
        importlib.import_module("scipy")
        importlib.import_module("matplotlib")
        importlib.import_module("pystoi")
        return {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib", "pystoi")}
    except (ImportError, importlib.metadata.PackageNotFoundError) as error:
        raise RuntimeError(f"offline analysis prerequisite unavailable: {error}") from error


def read_pcm(path):
    path = Path(path)
    if path.suffix == ".raw":
        content = path.read_bytes()
    elif path.suffix == ".wav":
        with wave.open(str(path), "rb") as stream:
            if (stream.getnchannels(), stream.getsampwidth(), stream.getframerate(), stream.getcomptype()) != (
                1,
                2,
                8000,
                "NONE",
            ):
                raise ValueError(f"expected mono 8-kHz signed-16 PCM: {path}")
            content = stream.readframes(stream.getnframes())
            if len(content) != stream.getnframes() * 2:
                raise ValueError(f"truncated PCM WAV: {path}")
    else:
        raise ValueError(f"PCM input extension must be .raw or .wav: {path}")
    if not content or len(content) % 2:
        raise ValueError(f"empty or incomplete signed-16 PCM: {path}")
    return np.frombuffer(content, dtype="<i2").astype(np.float64)


def aligned_signals(reference, decoded, metrics):
    r0, d0, count = (metrics[key] for key in ("reference_trim_head", "decoded_trim_head", "aligned_samples"))
    if any(type(value) is not int or value < 0 for value in (r0, d0, count)):
        raise ValueError("invalid recorded alignment support")
    if count == 0 or r0 + count > len(reference) or d0 + count > len(decoded):
        raise ValueError("recorded alignment exceeds signal support")
    return reference[r0 : r0 + count], decoded[d0 : d0 + count]


def reference_activity(reference):
    count = len(reference) // 160
    if not count:
        return np.zeros(0, dtype=bool)
    frames = np.asarray(reference[: count * 160], dtype=np.float64).reshape(count, 160)
    energy = np.sum(frames * frames, axis=1)
    return (energy > 0) & (energy >= float(np.max(energy)) * 1e-4)


def active_rms(reference, signal):
    if len(reference) != len(signal):
        raise ValueError("reference and signal support must match")
    active = reference_activity(reference)
    count = int(np.count_nonzero(active))
    if not count:
        raise ValueError("active RMS has no reference support")
    frames = np.asarray(signal[: len(active) * 160], dtype=np.float64).reshape(len(active), 160)
    return math.sqrt(float(np.sum(frames[active] ** 2)) / (count * 160))


def normalized_signals(reference, decoded, metrics):
    r, d = aligned_signals(reference, decoded, metrics)
    reference_rms, decoded_rms = active_rms(r, r), active_rms(r, d)
    if reference_rms <= 0 or decoded_rms <= 0:
        raise ValueError("speech normalization requires positive reference and decoded energy")
    gain = reference_rms / decoded_rms
    return r, d * gain if gain != 1 else d, gain


def stft_power(signal, size=256, hop=80):
    if type(size) is not int or type(hop) is not int or size < 2 or hop <= 0:
        raise ValueError("invalid STFT size/hop")
    if len(signal) < size:
        return np.empty((0, size // 2 + 1), dtype=np.float64)
    windows = np.lib.stride_tricks.sliding_window_view(np.asarray(signal, dtype=np.float64), size)[::hop]
    spectrum = np.fft.rfft(windows * np.hanning(size), axis=1)
    return spectrum.real**2 + spectrum.imag**2


def numpy_lsd(reference, decoded, metrics):
    if (
        metrics.get("reference_state") != "speech"
        or not metrics.get("active_frames")
        or not metrics.get("spectral_frames")
    ):
        return None
    r, d, _ = normalized_signals(reference, decoded, metrics)
    reference_power, decoded_power = stft_power(r), stft_power(d)
    if not len(reference_power):
        return None
    epsilon = max(float(np.max(reference_power)) * 1e-6, 1e-20)
    active = reference_activity(r)
    centers = (np.arange(len(reference_power)) * 80 + 128) // 160
    supported = centers < len(active)
    selected = np.zeros(len(centers), dtype=bool)
    selected[supported] = active[centers[supported]]
    if int(np.count_nonzero(selected)) != metrics["spectral_frames"]:
        raise ValueError("independent STFT support differs from the evaluator")
    if not np.any(selected):
        return None
    ratio = (reference_power[selected, 2:119] + epsilon) / (decoded_power[selected, 2:119] + epsilon)
    difference = 10 * np.log10(ratio)
    return float(np.mean(np.sqrt(np.mean(difference * difference, axis=1))))


def write_wav(path, samples):
    values = np.asarray(samples, dtype=np.float64)
    if values.ndim != 1 or not len(values) or not np.all(np.isfinite(values)):
        raise ValueError("WAV output requires a nonempty finite mono signal")
    rounded = np.rint(values)
    if np.any(rounded < -32768) or np.any(rounded > 32767):
        raise ValueError("WAV output would clip; apply a documented common attenuation first")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with wave.open(str(path), "wb") as stream:
        stream.setparams((1, 2, 8000, 0, "NONE", "not compressed"))
        stream.writeframes(rounded.astype("<i2").tobytes())
