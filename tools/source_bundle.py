#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Bundle an already validated release tag and its pinned dependency source."""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile

from pffft_metadata import read_pins
from update_pffft import download_archive

ROOT = Path(__file__).resolve().parents[1]


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args])


def check_path(name, inventory):
    path = PurePosixPath(name)
    if (path.is_absolute() or ".." in path.parts or not path.parts
            or path.parts[0] not in inventory["tracked_roots"]
            or any(part in inventory["forbidden_parts"] for part in path.parts)
            or path.suffix.lower() in inventory["forbidden_suffixes"]):
        raise ValueError(f"Source inventory excludes {name}")


def tracked_source(root, commit):
    # Read Git objects, never dirty/untracked working files or Git metadata.
    payload = git(root, "archive", "--format=tar", commit)
    files = {}
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:") as archive:
        for entry in archive:
            if entry.isdir():
                continue
            if not entry.isfile():
                raise ValueError(f"Unsupported tracked entry: {entry.name}")
            files[entry.name] = (archive.extractfile(entry).read(), entry.mode & 0o777)
    inventory = json.loads(files["legal/source-inventory.json"][0])
    for name in files:
        check_path(name, inventory)
    missing = set(inventory["required"]) - files.keys()
    if missing:
        raise ValueError(f"Missing required source/notice: {sorted(missing)}")
    # git archive obeys export-ignore/export-subst; reject any lost or altered
    # tracked blobs, since corresponding source must be the exact tag tree.
    tree = git(root, "ls-tree", "-rz", commit).split(b"\0")
    tracked = {}
    for row in filter(None, tree):
        metadata, name = row.split(b"\t", 1)
        mode, kind, oid = metadata.split()
        if kind != b"blob" or mode not in (b"100644", b"100755"):
            raise ValueError(f"Unsupported source tree object: {name!r}")
        tracked[name.decode("utf-8")] = (oid.decode("ascii"), int(mode, 8) & 0o777)
    if files.keys() != tracked.keys():
        raise ValueError("Archive must contain every tracked source file")
    for name, (data, _) in files.items():
        # Git blob identity only; SHA256 verifies every delivered source file.
        # nosemgrep: python.lang.security.insecure-hash-algorithms.insecure-hash-algorithm-sha1
        blob = hashlib.sha1(f"blob {len(data)}\0".encode() + data, usedforsecurity=False).hexdigest()
        if blob != tracked[name][0]:
            raise ValueError(f"Git archive changed tracked source: {name}")
        # Git archive's tar.umask is user configurable; use the tree mode so
        # two builders do not emit different permissions for the same tag.
        files[name] = (data, tracked[name][1])
    return files


def deterministic_archive(path, prefix, files, timestamp):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        with temporary.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for name, (data, mode) in sorted(files.items()):
                    entry = tarfile.TarInfo(f"{prefix}/{name}")
                    entry.size, entry.mode, entry.mtime = len(data), mode, timestamp
                    entry.uid = entry.gid = 0
                    entry.uname = entry.gname = ""
                    archive.addfile(entry, io.BytesIO(data))
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def create(root, tag, expected_commit, destination, dependency=None):
    if not re.fullmatch(r"v\d+\.\d+\.\d+", tag) or not re.fullmatch(r"[0-9a-f]{40}", expected_commit):
        raise ValueError("Expected a validated vX.Y.Z tag and full commit SHA")
    commit = git(root, "rev-parse", f"refs/tags/{tag}^{{commit}}").decode().strip()
    if commit != expected_commit:
        raise ValueError("Release tag differs from validated commit")
    files = tracked_source(root, commit)
    version = tag[1:]
    if not re.search(r"project\(mbelib-neo\s+VERSION\s+" + re.escape(version) + r"\s",
                     files["CMakeLists.txt"][0].decode()):
        raise ValueError("Release tag does not match source version")
    pins = read_pins_text(files["tools/ci-dependency-pins.env"][0])
    data = dependency.read_bytes() if dependency else download_archive(pins["PFFFT_SHA"])
    if hashlib.sha256(data).hexdigest() != pins["PFFFT_ARCHIVE_SHA256"]:
        raise ValueError("PFFFT archive SHA256 differs from release pin")
    files[f"third_party/pffft-{pins['PFFFT_SHA']}.tar.gz"] = (data, 0o644)
    timestamp = int(git(root, "show", "-s", "--format=%ct", commit))
    identity = {
        "SPDX-License-Identifier": "GPL-2.0-or-later", "schema_version": 1,
        "tag": tag, "commit": commit, "version": version, "source_date_epoch": timestamp,
        "pffft": pins, "files": [
            {"path": name, "sha256": hashlib.sha256(content).hexdigest(), "mode": mode}
            for name, (content, mode) in sorted(files.items())
        ],
    }
    files["SOURCE-IDENTITY.json"] = ((json.dumps(identity, indent=2) + "\n").encode(), 0o644)
    prefix = f"mbelib-neo-{version}-source"
    output = destination / f"{prefix}.tar.gz"
    deterministic_archive(output, prefix, files, timestamp)
    return output


def read_pins_text(data):
    # Same pin parser used by dependency tooling, without a working-tree input.
    return read_pins(text=data.decode("utf-8"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--output", type=Path, default=Path("dist"))
    parser.add_argument("--pffft-archive", type=Path)
    args = parser.parse_args()
    print(create(ROOT, args.tag, args.commit, args.output, args.pffft_archive))


if __name__ == "__main__":
    main()
