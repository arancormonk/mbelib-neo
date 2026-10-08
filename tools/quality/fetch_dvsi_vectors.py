#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
"""Fetch selected DVSI reference test-vector files for local measurement only.

DVSI publishes AMBE-3000 test vectors (input speech, DVSI-encoded bit streams and
DVSI-decoded speech) on https://www.dvsinc.com/dlapps/appsoft.shtml. The site's
terms grant no license to copy, modify, transfer or mirror its materials, so the
files are fetched into the git-ignored build tree and must never be committed,
redistributed or placed in published listening material. Only URLs, hashes and
aggregate measurements belong in the repository.

The archive is about 1.9 GB; this tool reads its central directory with HTTP
range requests (via curl) and downloads only the requested members, verifying each CRC-32
and recording SHA-256 digests in a manifest.
"""

import argparse
import fnmatch
import hashlib
import json
import struct
import subprocess
import sys
import zlib
from pathlib import Path, PurePosixPath

ARCHIVE_URL = "https://www.dvsinc.com/get-usb/tv-rc.zip"
DEFAULT_OUT = Path("build/quality/dvsi")
# The modes this library decodes: D-STAR, P25 full rate with and without FEC,
# and AMBE+2 rate 33 (with FEC) and rate 34 (without).
SUPPORTED_MODES = ("dstar", "p25", "p25_nofec", "r33", "r34")
# Speech inputs, each mode's encodings, and DVSI's rate conversions between
# those modes (tv-rc/<from>/<to>/, made with "-rc -rd <from> -re <to>").
DEFAULT_PATTERNS = (
    ("tv-rc/*.pcm", "tv-rc/cmprc.txt")
    + tuple(f"tv-rc/{mode}/*" for mode in SUPPORTED_MODES)
    + tuple(f"tv-rc/{src}/{dst}/*" for src in SUPPORTED_MODES for dst in SUPPORTED_MODES if src != dst)
)
EOCD_SIGNATURE = 0x06054B50
ZIP64_LOCATOR_SIGNATURE = 0x07064B50
ZIP64_EOCD_SIGNATURE = 0x06064B50
CENTRAL_SIGNATURE = 0x02014B50
LOCAL_SIGNATURE = 0x04034B50
USER_AGENT = "mbelib-neo-quality/1"


def split_response(raw):
    """Return the final header block and the body, skipping interim (1xx) and proxy CONNECT blocks."""
    head, separator, body = raw.partition(b"\r\n\r\n")
    while separator and body.startswith(b"HTTP/"):
        status = head.split(b"\r\n", 1)[0].split(b" ")
        code = status[1] if len(status) > 1 else b""
        if not (code.startswith(b"1") or b"connection established" in head.split(b"\r\n", 1)[0].lower()):
            break
        head, separator, body = body.partition(b"\r\n\r\n")
    if not separator:
        raise RuntimeError(f"{ARCHIVE_URL}: malformed response")
    return head, body


def fetch_range(start, length):
    """Read bytes [start, start+length) of the fixed archive URL over HTTPS only."""
    result = subprocess.run(
        [
            "curl", "--fail", "--silent", "--show-error", "--proto", "=https", "--max-time", "300",
            "--suppress-connect-headers", "--user-agent", USER_AGENT, "--range", f"{start}-{start + length - 1}",
            "--dump-header", "-", "--output", "-", ARCHIVE_URL,
        ],
        capture_output=True, check=False,
    )
    if result.returncode:
        raise RuntimeError(f"{ARCHIVE_URL}: curl failed: {result.stderr.decode('utf-8', errors='replace')}")
    head, data = split_response(result.stdout)
    if b" 206" not in head.split(b"\r\n", 1)[0]:
        raise RuntimeError(f"{ARCHIVE_URL}: server ignored the range request")
    total = b""
    for line in head.split(b"\r\n"):
        if line.lower().startswith(b"content-range:"):
            total = line.rpartition(b"/")[2].strip()
    if len(data) != length:
        raise RuntimeError(f"{ARCHIVE_URL}: short range read at {start} ({len(data)} of {length} bytes)")
    return data, int(total) if total.isdigit() else None


def archive_size():
    _, total = fetch_range(0, 1)
    if total is None:
        raise RuntimeError(f"{ARCHIVE_URL}: missing Content-Range total")
    return total


def find_central_directory(size):
    tail_length = min(size, 65536 + 22)
    tail, _ = fetch_range(size - tail_length, tail_length)
    position = tail.rfind(struct.pack("<I", EOCD_SIGNATURE))
    if position < 0:
        raise RuntimeError("end of central directory not found")
    (_, _, _, _, count, cd_size, cd_offset, _) = struct.unpack_from("<IHHHHIIH", tail, position)
    if 0xFFFFFFFF in (cd_size, cd_offset) or count == 0xFFFF:
        locator = position - 20
        if locator < 0 or struct.unpack_from("<I", tail, locator)[0] != ZIP64_LOCATOR_SIGNATURE:
            raise RuntimeError("ZIP64 locator not found")
        zip64_offset = struct.unpack_from("<Q", tail, locator + 8)[0]
        record, _ = fetch_range(zip64_offset, 56)
        if struct.unpack_from("<I", record, 0)[0] != ZIP64_EOCD_SIGNATURE:
            raise RuntimeError("ZIP64 end of central directory not found")
        count, cd_size, cd_offset = struct.unpack_from("<QQQ", record, 32)
    return cd_offset, cd_size, count


def zip64_values(extra, size, compressed, offset):
    position = 0
    while position + 4 <= len(extra):
        header, length = struct.unpack_from("<HH", extra, position)
        body = extra[position + 4:position + 4 + length]
        if header == 0x0001:
            fields = list(struct.unpack_from(f"<{len(body) // 8}Q", body))
            if size == 0xFFFFFFFF:
                size = fields.pop(0)
            if compressed == 0xFFFFFFFF:
                compressed = fields.pop(0)
            if offset == 0xFFFFFFFF:
                offset = fields.pop(0)
        position += 4 + length
    return size, compressed, offset


def parse_central_directory(blob, count):
    entries = []
    position = 0
    for _ in range(count):
        if struct.unpack_from("<I", blob, position)[0] != CENTRAL_SIGNATURE:
            raise RuntimeError("corrupt central directory")
        fields = struct.unpack_from("<IHHHHHHIIIHHHHHII", blob, position)
        method, crc, compressed, size = fields[4], fields[7], fields[8], fields[9]
        name_length, extra_length, comment_length = fields[10], fields[11], fields[12]
        offset = fields[16]
        name = blob[position + 46:position + 46 + name_length].decode("cp437")
        extra = blob[position + 46 + name_length:position + 46 + name_length + extra_length]
        size, compressed, offset = zip64_values(extra, size, compressed, offset)
        entries.append(
            {"name": name, "method": method, "crc32": crc, "compressed": compressed, "size": size, "offset": offset}
        )
        position += 46 + name_length + extra_length + comment_length
    return entries


def member_matches(name, pattern):
    """Glob per path segment, so '*' never crosses a directory separator."""
    parts = name.split("/")
    globs = pattern.split("/")
    return len(parts) == len(globs) and all(fnmatch.fnmatchcase(p, g) for p, g in zip(parts, globs))


def safe_relative(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"unsafe archive member name: {name!r}")
    return Path(*path.parts)


def extract(entry, destination):
    header, _ = fetch_range(entry["offset"], 30)
    if struct.unpack_from("<I", header, 0)[0] != LOCAL_SIGNATURE:
        raise RuntimeError(f"bad local header for {entry['name']}")
    name_length, extra_length = struct.unpack_from("<HH", header, 26)
    start = entry["offset"] + 30 + name_length + extra_length
    payload = fetch_range(start, entry["compressed"])[0] if entry["compressed"] else b""
    if entry["method"] == 0:
        data = payload
    elif entry["method"] == 8:
        data = zlib.decompress(payload, -15)
    else:
        raise RuntimeError(f"{entry['name']}: unsupported compression method {entry['method']}")
    if len(data) != entry["size"] or (zlib.crc32(data) & 0xFFFFFFFF) != entry["crc32"]:
        raise RuntimeError(f"{entry['name']}: size or CRC-32 mismatch")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    temporary.write_bytes(data)
    temporary.replace(destination)
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Fetch selected DVSI test-vector files for local measurement.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--include", action="append", help="archive member glob (repeatable)")
    parser.add_argument("--list", action="store_true", help="list matching members without downloading")
    args = parser.parse_args()
    patterns = tuple(args.include or DEFAULT_PATTERNS)

    size = archive_size()
    cd_offset, cd_size, count = find_central_directory(size)
    directory, _ = fetch_range(cd_offset, cd_size)
    entries = [
        entry
        for entry in parse_central_directory(directory, count)
        if not entry["name"].endswith("/") and any(member_matches(entry["name"], p) for p in patterns)
    ]
    if args.list:
        for entry in entries:
            print(f"{entry['size']:>10} {entry['name']}")
        print(f"{len(entries)} members, {sum(e['size'] for e in entries)} bytes", file=sys.stderr)
        return 0

    manifest = {
        "schema_version": 1,
        "url": ARCHIVE_URL,
        "archive_bytes": size,
        "patterns": list(patterns),
        "notice": "DVSI materials: local measurement only; do not commit or redistribute.",
        "members": [],
    }
    for entry in entries:
        relative = safe_relative(entry["name"])
        digest = extract(entry, args.out / relative)
        manifest["members"].append(
            {"name": entry["name"], "bytes": entry["size"], "crc32": f"{entry['crc32']:08x}", "sha256": digest}
        )
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "dvsi-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"fetched {len(entries)} members into {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
