#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Require a complete version-specific release set and checksum every payload."""

import argparse
import hashlib
from pathlib import Path
import re


def needs_source(version):
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise ValueError("Expected a numeric X.Y.Z release version")
    return tuple(map(int, version.split("."))) >= (2, 6, 0)


def collect(directory, version):
    prefix = f"mbelib-neo-{version}-"
    expected = {"linux", "macos", "windows"}
    if needs_source(version):
        expected.add("source")
    packages = sorted([*directory.glob("*.tar.gz"), *directory.glob("*.zip")])
    platforms = set()
    for package in packages:
        match = re.fullmatch(re.escape(prefix) + r"(linux|macos|windows)-(x64|arm64)\.(tar\.gz|zip)", package.name)
        if match:
            platform = match[1]
            if (platform == "windows") != (match[3] == "zip"):
                raise ValueError("Incorrect platform archive format")
        elif package.name == prefix + "source.tar.gz" and "source" in expected:
            platform = "source"
        else:
            raise ValueError(f"Unsupported release package: {package.name}")
        if platform in platforms:
            raise ValueError(f"Duplicate platform: {platform}")
        platforms.add(platform)
    if platforms != expected:
        raise ValueError(f"Expected release packages for {sorted(expected)}, found {sorted(platforms)}")
    sboms = sorted(directory.glob("*.spdx.json"))
    if {s.name for s in sboms} != {p.name + ".spdx.json" for p in packages}:
        raise ValueError("Expected exactly one named SBOM for every archive")
    if any(p.is_symlink() or not p.is_file() for p in [*packages, *sboms]):
        raise ValueError("Release assets must be ordinary files")
    return packages, sboms


def prepare(directory, version):
    packages, sboms = collect(directory, version)
    lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in sorted(packages + sboms)]
    (directory / "SHA256SUMS.txt").write_text("".join(lines), encoding="utf-8")
    (directory / "RELEASE-NOTES.md").write_text(release_notes(version), encoding="utf-8")
    return "".join(lines)


def release_notes(version):
    if not needs_source(version):
        return "Historical package format: three platform archives.\n"
    base = "https://github.com/arancormonk/mbelib-neo"
    source = f"mbelib-neo-{version}-source.tar.gz"
    return (f"Complete corresponding source: [{source}]({base}/releases/download/v{version}/{source}).\n\n"
            f"Read the [component inventory and unresolved licensing questions]({base}/blob/v{version}/legal/status.md) "
            f"and [patent research]({base}/blob/v{version}/legal/patents.md). "
            "Source delivery and conformance do not establish legal clearance.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--directory", type=Path, default=Path("dist"))
    args = parser.parse_args()
    if not args.prefix.startswith("mbelib-neo-"):
        parser.error("Expected mbelib-neo-X.Y.Z prefix")
    print(prepare(args.directory, args.prefix[len("mbelib-neo-"):]), end="")


if __name__ == "__main__":
    main()
