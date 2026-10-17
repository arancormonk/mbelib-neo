#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Describe release contents and embedded components using the maintained inventory."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
import zipfile

from legal_inventory import load, owners, validate
from pffft_metadata import read_pins

ROOT = Path(__file__).resolve().parents[1]
AGGREGATE = "SPDXRef-Package-mbelib-neo"


def spdx_sha1(content):
    # SPDX 2.3 package verification codes require SHA1. Security checks and
    # release attestations use the separate SHA256 checksums.
    # nosemgrep: python.lang.security.insecure-hash-algorithms.insecure-hash-algorithm-sha1
    return hashlib.sha1(content, usedforsecurity=False).hexdigest()


def package_files(path):
    files = {}
    if path.name.endswith(".zip"):
        with zipfile.ZipFile(path) as archive:
            entries = [(e.filename, archive.read(e)) for e in archive.infolist() if not e.is_dir()]
    else:
        with tarfile.open(path, "r:gz") as archive:
            entries = []
            for entry in archive:
                if entry.isdir():
                    continue
                if entry.issym():
                    # Installed Unix libraries use versioned relative symlinks.
                    target = PurePosixPath(entry.linkname)
                    if target.is_absolute() or ".." in target.parts:
                        raise ValueError("Unsafe package symlink")
                    entries.append((entry.name, entry.linkname.encode()))
                elif entry.isfile():
                    entries.append((entry.name, archive.extractfile(entry).read()))
                else:
                    raise ValueError("Unsupported release member")
    roots = set()
    for name, data in entries:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or len(path.parts) < 2:
            raise ValueError(f"Unsafe package member: {name}")
        roots.add(path.parts[0])
        relative = str(PurePosixPath(*path.parts[1:]))
        if relative in files:
            raise ValueError(f"Duplicate package member: {name}")
        files[relative] = data
    if len(roots) != 1:
        raise ValueError("Expected one package root")
    return files


def validate_package(files, kind, root=ROOT):
    inventory = load(root)
    if kind == "source":
        identity = json.loads(files["SOURCE-IDENTITY.json"])
        records = identity["files"]
        paths = [r["path"] for r in records]
        if len(set(paths)) != len(paths) or set(paths) != files.keys() - {"SOURCE-IDENTITY.json"}:
            raise ValueError("Source identity does not cover exact contents")
        for record in records:
            if hashlib.sha256(files[record["path"]]).hexdigest() != record["sha256"]:
                raise ValueError("Source identity hash mismatch")
        inventory = json.loads(files["legal/components.json"])
        names = [name for name in files if not name.startswith("third_party/") and name != "SOURCE-IDENTITY.json"]
        validate(inventory, files.__getitem__, names)
        staging = json.loads(files["legal/source-inventory.json"])
        from source_bundle import check_path
        for name in names:
            check_path(name, staging)
        if set(staging["required"]) - files.keys():
            raise ValueError("Required source notice or build input is absent")
        pins = read_pins(text=files["tools/ci-dependency-pins.env"].decode())
        dependency = f"third_party/pffft-{pins['PFFFT_SHA']}.tar.gz"
        if set(name for name in files if name.startswith("third_party/")) != {dependency}:
            raise ValueError("Unexpected bundled dependency")
        if hashlib.sha256(files[dependency]).hexdigest() != pins["PFFFT_ARCHIVE_SHA256"]:
            raise ValueError("PFFFT archive hash mismatch")
        if identity["pffft"] != pins:
            raise ValueError("Source identity dependency mismatch")
        return inventory, pins, identity, files.__getitem__
    # Installed notices are validated by contents, independent of lib/lib64 layout.
    for c in inventory["components"]:
        if "binary" in c["distribution"]:
            for notice in c["notices"]:
                expected = (root / notice).read_bytes().replace(b"\r\n", b"\n")
                if not any(content.replace(b"\r\n", b"\n") == expected for content in files.values()):
                    raise ValueError(f"Binary package lacks retained notice {notice}")
    retained = (root / "LICENSES/LicenseRef-PFFFT.txt").read_bytes().replace(b"\r\n", b"\n")
    marker = b"PFFFT: https://github.com/arancormonk/pffft\n"
    expected_notices = retained[retained.index(marker):]
    if not any(name.endswith("/pffft/PFFFT-NOTICES.txt") and content.replace(b"\r\n", b"\n") == expected_notices
               for name, content in files.items()):
        raise ValueError("Missing compiled PFFFT notices")
    forbidden = {".wav", ".pcm", ".raw", ".mp3", ".pdf", ".zip"}
    if any(PurePosixPath(n).suffix.lower() in forbidden or ".git" in PurePosixPath(n).parts
           or not n.startswith(("bin/", "lib/", "lib64/", "include/", "share/")) for n in files):
        raise ValueError("Unexpected binary package contents")
    return inventory, read_pins(root / "tools/ci-dependency-pins.env"), {}, lambda n: (root / n).read_bytes()


def create(path, kind, version, commit, epoch, root=ROOT):
    files = package_files(path)
    inventory, pins, identity, read = validate_package(files, kind, root)
    if inventory["release"] != version:
        raise ValueError("Component inventory has not been reviewed for this release version")
    if identity and (identity["version"] != version or identity["commit"] != commit
                     or identity["tag"] != f"v{version}"):
        raise ValueError("Source identity differs from validated release")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    document = {
        "spdxVersion": "SPDX-2.3", "dataLicense": "CC0-1.0", "SPDXID": "SPDXRef-DOCUMENT",
        "name": path.name, "documentNamespace": f"https://spdx.org/spdxdocs/mbelib-neo-{digest}",
        "creationInfo": {"creators": ["Tool: mbelib-neo-release-sbom"],
                         "created": datetime.fromtimestamp(epoch, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
        "packages": [{"name": "mbelib-neo", "SPDXID": AGGREGATE, "versionInfo": version,
                      "downloadLocation": f"https://github.com/arancormonk/mbelib-neo/releases/download/v{version}/{path.name}",
                      "filesAnalyzed": True, "licenseDeclared": "NOASSERTION", "licenseConcluded": "NOASSERTION",
                      "copyrightText": "NOASSERTION", "sourceInfo": f"Release commit {commit}; see legal/status.md",
                      "checksums": [{"algorithm": "SHA256", "checksumValue": digest}]}],
        "relationships": [{"spdxElementId": "SPDXRef-DOCUMENT", "relationshipType": "DESCRIBES",
                           "relatedSpdxElement": AGGREGATE}], "files": [], "hasExtractedLicensingInfos": [],
    }
    sha1s = []
    for index, (name, content) in enumerate(sorted(files.items())):
        sha1 = spdx_sha1(content)
        sha1s.append(sha1)
        file_id = f"SPDXRef-File-{index}"
        document["files"].append({"fileName": "./" + name, "SPDXID": file_id,
                                  "checksums": [{"algorithm": "SHA1", "checksumValue": sha1},
                                                {"algorithm": "SHA256", "checksumValue": hashlib.sha256(content).hexdigest()}],
                                  "licenseConcluded": "NOASSERTION", "licenseInfoInFiles": ["NOASSERTION"],
                                  "copyrightText": "NOASSERTION"})
        document["relationships"].append({"spdxElementId": AGGREGATE, "relationshipType": "CONTAINS", "relatedSpdxElement": file_id})
    document["packages"][0]["packageVerificationCode"] = {
        "packageVerificationCodeValue": spdx_sha1("".join(sorted(sha1s)).encode())}
    add_components(document, inventory, pins, kind, commit, read)
    return document


def add_components(document, inventory, pins, kind, commit, read):
    for c in inventory["components"]:
        if kind not in c["distribution"]:
            continue
        package_id = "SPDXRef-Component-" + c["id"]
        package = {"name": c["name"], "SPDXID": package_id, "filesAnalyzed": False,
                   "versionInfo": commit if c["revision"] == "release-commit" else c["revision"],
                   "downloadLocation": c["origin"] if c["origin"].startswith("https://") else "NOASSERTION",
                   "licenseDeclared": c["license_declared"], "licenseConcluded": c["license_concluded"],
                   "copyrightText": "NOASSERTION", "licenseComments": c["assessment"],
                   "sourceInfo": "; ".join(c["sources"])}
        if c["id"] == "pffft":
            package.update(versionInfo=pins["PFFFT_SHA"],
                           downloadLocation=f"https://codeload.github.com/arancormonk/pffft/tar.gz/{pins['PFFFT_SHA']}",
                           checksums=[{"algorithm": "SHA256", "checksumValue": pins["PFFFT_ARCHIVE_SHA256"]}])
        document["packages"].append(package)
        relation = "STATIC_LINK" if c["id"] == "pffft" and kind == "binary" else "CONTAINS"
        document["relationships"].append({"spdxElementId": AGGREGATE, "relationshipType": relation,
                                          "relatedSpdxElement": package_id})
        license_id = c["license_declared"]
        if license_id.startswith("LicenseRef-"):
            notice = next(n for n in c["notices"] if license_id in n)
            document["hasExtractedLicensingInfos"].append({"licenseId": license_id, "name": c["name"],
                                                          "extractedText": read(notice).decode("utf-8")})
    if kind == "source":
        # Preserve overlapping first-party/ISC provenance without asserting
        # an automatic copyright conclusion for each file.
        for f in document["files"]:
            for c in owners(f["fileName"][2:], inventory):
                document["relationships"].append({"spdxElementId": "SPDXRef-Component-" + c["id"],
                                                  "relationshipType": "CONTAINS", "relatedSpdxElement": f["SPDXID"]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--kind", choices=("source", "binary"), required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--epoch", type=int, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"\d+\.\d+\.\d+", args.version) or not re.fullmatch(r"[0-9a-f]{40}", args.commit):
        parser.error("Expected release version and full commit SHA")
    result = create(args.archive, args.kind, args.version, args.commit, args.epoch)
    Path(str(args.archive) + ".spdx.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
