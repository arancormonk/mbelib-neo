#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Preserve reference-data evidence without inferring redistribution permission."""

import argparse
from datetime import date, datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            result.update(block)
    return result.hexdigest()


def acquisition(source_url, terms_url, *, retrieved_at=None, status="unknown", notices=()):
    if retrieved_at is not None:
        date.fromisoformat(retrieved_at)
    if status not in ("allowed", "unknown", "restricted"):
        raise ValueError("Unknown redistribution status")
    return {"source_url": source_url, "terms_url": terms_url, "retrieved_at": retrieved_at,
            "recorded_at": datetime.now(timezone.utc).isoformat(), "redistribution": status,
            "notices": list(notices), "intended_use": "local measurement; no permission inferred"}


def distribution(manifest):
    references = []
    for reference in manifest.get("references", []):
        source = reference.get("source", reference)
        evidence = source.get("acquisition", {})
        references.append({"name": reference["name"], "acquisition": evidence,
                           "redistribution": evidence.get("redistribution", "unknown"),
                           "permission_source": evidence.get("permission_source")})
    allowed = bool(references) and all(r["redistribution"] == "allowed" and r["permission_source"] for r in references)
    return {"schema_version": 1, "redistribution": "allowed" if allowed else "local-only",
            "references": references, "notices": manifest.get("corpus", {}).get("notices", []),
            "assessment": "Recorded evidence only; unknown permissions remain unknown."}


def safe_path(root, relative):
    path = PurePosixPath(relative)
    target = root / relative
    if path.is_absolute() or ".." in path.parts or not target.resolve().is_relative_to(root.resolve()):
        raise ValueError("Unsafe reference notice path")
    return target


def carry_notices(root, manifest, presentation):
    record = distribution(manifest)
    # Older immutable runs already retained notices beside the corpus manifest.
    # Read that evidence without changing the run or upgrading unknown rights.
    corpus_manifest = manifest.get("corpus", {}).get("manifest")
    if not record["notices"] and corpus_manifest:
        copied = safe_path(root, corpus_manifest)
        metadata = json.loads(copied.read_text(encoding="utf-8"))
        record["notices"] = [{**n, "path": str(PurePosixPath(corpus_manifest).parent / n["file"])}
                             for n in metadata.get("notices", [])]
    retained = []
    for index, notice in enumerate(record["notices"]):
        source = safe_path(root, notice["path"])
        if digest(source) != notice["sha256"]:
            raise ValueError("Reference notice hash mismatch")
        relative = f"notices/{index}-{source.name}"
        target = safe_path(presentation, relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        retained.append({**notice, "path": relative})
    record["notices"] = retained
    presentation.mkdir(parents=True, exist_ok=True)
    (presentation / "DISTRIBUTION.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return ["DISTRIBUTION.json", *(n["path"] for n in retained)]


def check_public(directory):
    record = json.loads((directory / "DISTRIBUTION.json").read_text(encoding="utf-8"))
    refs = record.get("references", [])
    if (record.get("redistribution") != "allowed" or not refs
            or any(r.get("redistribution") != "allowed" or not r.get("permission_source") for r in refs)):
        raise ValueError("Public distribution lacks explicit permission evidence for every reference")
    actual = {n["sha256"] for n in record.get("notices", [])}
    required = {n["sha256"] for r in refs for n in r.get("acquisition", {}).get("notices", [])}
    if not required.issubset(actual):
        raise ValueError("Required corpus notices are absent")
    for notice in record.get("notices", []):
        if digest(safe_path(directory, notice["path"])) != notice["sha256"]:
            raise ValueError("Reference notice hash mismatch")


def retain_dataset(manifest_path, output):
    """Carry acquisition evidence and notices into a local measurement report."""
    if not manifest_path.is_file():
        return {"redistribution": "unknown", "acquisition": {}, "notices": []}
    source = json.loads(manifest_path.read_text(encoding="utf-8"))
    evidence = source.get("acquisition", {})
    retained = []
    for index, notice in enumerate(evidence.get("notices", [])):
        original = safe_path(manifest_path.parent, notice["file"])
        if digest(original) != notice["sha256"]:
            raise ValueError("Reference notice hash mismatch")
        relative = f"notices/{index}-{original.name}"
        target = safe_path(output, relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(original.read_bytes())
        retained.append({**notice, "file": relative})
    return {"redistribution": evidence.get("redistribution", "unknown"), "acquisition": evidence,
            "manifest_sha256": digest(manifest_path), "notices": retained}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-public", type=Path)
    parser.add_argument("--file", type=Path)
    parser.add_argument("--source-url")
    parser.add_argument("--terms-url")
    parser.add_argument("--retrieved-at", help="actual YYYY-MM-DD retrieval date, otherwise unknown")
    parser.add_argument("--notice", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check_public:
        check_public(args.check_public)
        return
    if not all((args.file, args.source_url, args.terms_url, args.output)):
        parser.error("recording requires --file, --source-url, --terms-url and --output")
    if args.output.resolve() in {p.resolve() for p in [args.file, *args.notice]}:
        parser.error("output must not overwrite reference data or notices")
    notices = [{"file": str(p), "sha256": digest(p)} for p in args.notice]
    record = acquisition(args.source_url, args.terms_url, retrieved_at=args.retrieved_at, notices=notices)
    record.update(file=str(args.file), sha256=digest(args.file))
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
