#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Prepare pinned development speech or untouched, independently held-out speech."""

import argparse
import hashlib
import http.client
import json
import shutil
import struct
import subprocess
import sys
import tempfile
import wave
from pathlib import Path
from urllib.parse import urlsplit

QUALITY = Path("build/quality")
CODEC2_REVISION = "310777b1c6f1af0bc7c72f5b32f80f6fd9136962"
CODEC2_URL = "https://github.com/drowe67/codec2.git"
LEGACY = {"hts1a": 48000, "hts2a": 48000, "hts1": 96000, "kristoff": 80000, "ve9qrp_10s": 160000}
SPEAKERS = ("bdl", "rms", "clb", "slt")
ARCTIC_ROOT = "http://festvox.org/cmu_arctic/cmu_arctic"


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            result.update(block)
    return result.hexdigest()


def run(argv):
    result = subprocess.run([str(part) for part in argv], capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(f"command failed: {argv!r}\n{result.stderr.decode('utf-8', errors='replace')}")
    return result.stdout


def write_manifest(directory, manifest):
    (directory / "corpus-manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def publish(pairs, staging):
    """Publish only after every asset is ready; restore previous sets on failure."""
    for _, destination in pairs:
        if destination.is_symlink() or (destination.exists() and not destination.is_dir()):
            raise ValueError(f"corpus destination must be a real directory: {destination}")
    installed = []
    backups = []
    try:
        for index, (prepared, destination) in enumerate(pairs):
            if destination.exists():
                backup = staging / f"previous-{index}"
                destination.rename(backup)
                backups.append((backup, destination))
            prepared.rename(destination)
            installed.append((prepared, destination))
    except BaseException:
        for prepared, destination in reversed(installed):
            destination.rename(prepared)
        for backup, destination in reversed(backups):
            backup.rename(destination)
        raise


def prepare_legacy(staging):
    source = QUALITY / "codec2-src"
    if source.exists():
        actual = run(["git", "-C", source, "rev-parse", "HEAD"]).decode().strip()
        if actual != CODEC2_REVISION:
            raise ValueError(f"existing {source} is not pinned Codec2 {CODEC2_REVISION}; preserving it unchanged")
        checkout = source
    else:
        checkout = staging / "codec2-src"
        run(["tools/fetch-pinned-git.sh", CODEC2_URL, CODEC2_REVISION, checkout])
    prepared = staging / "legacy"
    prepared.mkdir()
    references = []
    # Read the pinned Git objects, not potentially edited checkout files.
    for name, expected_size in LEGACY.items():
        content = run(["git", "-C", checkout, "show", f"{CODEC2_REVISION}:raw/{name}.raw"])
        if len(content) != expected_size:
            raise ValueError(
                f"unexpected pinned asset size: {name}.raw ({len(content)} bytes; expected {expected_size})"
            )
        output = prepared / f"{name}.raw"
        output.write_bytes(content)
        sha256 = digest(output)
        source_asset = checkout / "raw" / f"{name}.raw"
        if source_asset.read_bytes() != content:
            raise ValueError(f"checkout asset differs from pinned Git object: {source_asset}")
        references.append(
            {
                "name": name,
                "file": output.name,
                "samples": len(content) // 2,
                "sha256": sha256,
                "speaker": None,
                "utterance_index": None,
                "source_url": f"https://raw.githubusercontent.com/drowe67/codec2/{CODEC2_REVISION}/raw/{name}.raw",
                "source_sha256": sha256,
                "source_file": str(source / "raw" / output.name),
                "source_revision": CODEC2_REVISION,
                "duplicate_of": "hts1" if name == "hts1a" else None,
            }
        )
    prefix = (prepared / "hts1a.raw").read_bytes()
    if (prepared / "hts1.raw").read_bytes()[: len(prefix)] != prefix:
        raise ValueError("pinned hts1a is not the expected exact 24000-sample prefix of hts1")
    duplicates = [{"name": "hts1a", "duplicate_of": "hts1", "relation": "exact_prefix", "samples": len(prefix) // 2}]
    write_manifest(
        prepared,
        {
            "schema_version": 1,
            "set": "legacy",
            "label": "legacy/development",
            "level": "native",
            "references": references,
            "notices": [],
            "duplicates": duplicates,
        },
    )
    if checkout != source:
        checkout.rename(source)
    publish([(prepared, QUALITY / "corpus")], staging)
    return [QUALITY / "corpus"]


def validate_wav(path):
    try:
        with wave.open(str(path), "rb") as audio:
            if (audio.getnchannels(), audio.getframerate(), audio.getsampwidth(), audio.getcomptype()) != (
                1,
                16000,
                2,
                "NONE",
            ):
                raise ValueError(f"expected mono 16-kHz signed-16 PCM WAV: {path}")
            frames = audio.getnframes()
            if frames <= 0 or len(audio.readframes(frames)) != frames * 2:
                raise ValueError(f"empty or truncated PCM WAV: {path}")
    except (wave.Error, EOFError) as error:
        raise ValueError(f"invalid WAV {path}: {error}") from error


def validate_notice(path):
    content = path.read_bytes()
    if not content.strip() or b"<html" in content.lower() or b"<!doctype html" in content.lower():
        raise ValueError(f"missing or invalid COPYING notice: {path}")
    content.decode("utf-8")


def fetch_asset(url, destination, validator):
    # Cache entries are individually atomic and validated even on reuse. Invalid
    # cache entries fail explicitly rather than quietly substituting other audio.
    parsed = urlsplit(url)
    if parsed.scheme != "http" or parsed.netloc != "festvox.org" or parsed.query or parsed.fragment:
        raise ValueError(f"expected a fixed Festvox HTTP corpus URL: {url}")
    if destination.exists():
        validator(destination)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".download-", dir=destination.parent) as temporary:
        staged = Path(temporary) / destination.name
        connection = http.client.HTTPConnection("festvox.org", timeout=60)
        try:
            connection.request("GET", parsed.path)
            with connection.getresponse() as response:
                if response.status != 200:
                    raise RuntimeError(f"corpus download failed: {url} (HTTP {response.status})")
                with staged.open("wb") as output:
                    shutil.copyfileobj(response, output)
        finally:
            connection.close()
        validator(staged)
        staged.rename(destination)


def attenuate(native, output):
    scale = 10 ** (-12 / 20)
    source = native.read_bytes()
    if not source or len(source) % 2:
        raise ValueError(f"conversion did not produce nonempty s16le PCM: {native}")
    converted = bytearray(len(source))
    for index, (sample,) in enumerate(struct.iter_unpack("<h", source)):
        # Python round, like np.rint, rounds exact ties to the even integer.
        rounded = round(sample * scale)
        saturated = max(-32768, min(32767, rounded))
        if saturated != rounded:
            raise ValueError(f"-12 dB attenuation unexpectedly clips: {native} sample {index}")
        struct.pack_into("<h", converted, index * 2, saturated)
    output.write_bytes(converted)


def prepare_heldout(staging):
    ffmpeg_version = run(["ffmpeg", "-version"]).decode("utf-8").strip()
    assets = []
    notices = []
    # Fetch and validate the entire exact source set before beginning conversion.
    for speaker in SPEAKERS:
        root = f"{ARCTIC_ROOT}/cmu_us_{speaker}_arctic"
        source_directory = QUALITY / "arctic-src" / speaker
        notice = source_directory / "COPYING"
        notice_url = f"{root}/COPYING"
        fetch_asset(notice_url, notice, validate_notice)
        notices.append(
            {
                "speaker": speaker,
                "file": f"notices/{speaker}.COPYING",
                "sha256": digest(notice),
                "source_url": notice_url,
            }
        )
        for index in range(1, 11):
            filename = f"arctic_a{index:04d}.wav"
            source = source_directory / "wav" / filename
            url = f"{root}/wav/{filename}"
            fetch_asset(url, source, validate_wav)
            assets.append((speaker, index, source, url))
    native = staging / "heldout"
    reduced = staging / "heldout-minus12db"
    for directory in (native, reduced):
        (directory / "notices").mkdir(parents=True)
        for notice in notices:
            shutil.copyfile(QUALITY / "arctic-src" / notice["speaker"] / "COPYING", directory / notice["file"])
    native_references = []
    reduced_references = []
    for speaker, index, source, url in assets:
        name = f"{speaker}_arctic_a{index:04d}"
        output = native / f"{name}.raw"
        argv = [
            "ffmpeg",
            "-nostdin",
            "-v",
            "error",
            "-i",
            str(source),
            "-map",
            "0:a:0",
            "-ac",
            "1",
            "-ar",
            "8000",
            "-c:a",
            "pcm_s16le",
            "-f",
            "s16le",
            str(output),
        ]
        run(argv)
        attenuate(output, reduced / output.name)
        reference = {
            "name": name,
            "file": output.name,
            "samples": output.stat().st_size // 2,
            "sha256": digest(output),
            "speaker": speaker,
            "utterance_index": index,
            "source_url": url,
            "source_sha256": digest(source),
            "source_file": str(source),
            "conversion": {"argv": argv, "ffmpeg_version": ffmpeg_version},
        }
        native_references.append(reference)
        reduced_references.append(
            {
                **reference,
                "sha256": digest(reduced / output.name),
                "attenuation_db": -12,
                "native_sha256": reference["sha256"],
                "attenuation": {
                    "formula": "round(sample * 10**(-12/20))",
                    "rounding": "ties_even",
                    "clipped_samples": 0,
                },
            }
        )
    for directory, references, level in (
        (native, native_references, "native"),
        (reduced, reduced_references, "minus12db"),
    ):
        write_manifest(
            directory,
            {
                "schema_version": 1,
                "set": "heldout",
                "label": "heldout",
                "level": level,
                "references": references,
                "notices": notices,
                "duplicates": [],
            },
        )
    destinations = [QUALITY / "corpus-heldout", QUALITY / "corpus-heldout-minus12db"]
    publish(list(zip((native, reduced), destinations)), staging)
    return destinations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--set", choices=("legacy", "heldout"), default="legacy")
    args = parser.parse_args()
    if not Path("tools/fetch-pinned-git.sh").is_file():
        parser.error("run from the repository root")
    QUALITY.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".corpus-stage-", dir=QUALITY) as temporary:
        staging = Path(temporary)
        directories = prepare_legacy(staging) if args.set == "legacy" else prepare_heldout(staging)
    for directory in directories:
        print(directory / "corpus-manifest.json")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, http.client.HTTPException) as error:
        print(f"corpus preparation failed: {error}", file=sys.stderr)
        sys.exit(1)
