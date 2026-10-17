#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Validate component evidence and render the inventory used by release SBOMs."""

import argparse
from datetime import date
import fnmatch
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FIELDS = {"id", "name", "paths", "distribution", "origin", "revision", "license_declared",
          "license_concluded", "notices", "copyright_permission", "gpl_compatibility",
          "patent_status", "assessment", "sources", "review_date", "notice_sha256"}


def load(root=ROOT):
    return json.loads((root / "legal/components.json").read_text(encoding="utf-8"))


def owners(name, inventory):
    return [component for component in inventory["components"]
            if any(fnmatch.fnmatchcase(name, pattern) for pattern in component["paths"])
            and not any(fnmatch.fnmatchcase(name, pattern) for pattern in component.get("exclude", []))]


def validate(inventory, read, names):
    if inventory.get("schema_version") != 1 or inventory.get("aggregate_license_concluded") != "NOASSERTION":
        raise ValueError("Expected reviewed inventory schema and conservative aggregate conclusion")
    ids = set()
    for component in inventory["components"]:
        if FIELDS - component.keys() or component["id"] in ids:
            raise ValueError("Incomplete or duplicate component record")
        ids.add(component["id"])
        date.fromisoformat(component["review_date"])
        for notice in component["notices"]:
            content = read(notice)
            if not content.strip():
                raise ValueError(f"Empty notice: {notice}")
            if hashlib.sha256(content.replace(b"\r\n", b"\n")).hexdigest() != component["notice_sha256"].get(notice):
                raise ValueError(f"Retained notice hash mismatch: {notice}")
        for source in component["sources"]:
            if "://" not in source:
                read(source.split("#", 1)[0])
        if component["distribution"] and component["copyright_permission"] == "unresolved":
            exception = component.get("exception", "")
            if not exception or exception not in read("legal/status.md").decode():
                raise ValueError("Distributed unresolved component lacks a disclosed exception")
        declared = component["license_declared"]
        if declared.startswith("LicenseRef-") and not any(declared in n for n in component["notices"]):
            raise ValueError(f"Custom license lacks extracted text: {declared}")
    required = {"first-party", "isc-origins", "pffft", "etsi-data", "op25-tooling",
                "dvsi-reference", "etsi-reference", "speech-corpora", "license-texts"}
    if not required.issubset(ids):
        raise ValueError("Component inventory is incomplete")
    for name in names:
        if not owners(name, inventory):
            raise ValueError(f"Uncovered source: {name}")
        if name.endswith((".c", ".h", ".cc", ".patch")):
            first_line = read(name).decode("utf-8").splitlines()[0]
            licenses = {c["license_declared"] for c in owners(name, inventory)}
            identifier = re.search(r"SPDX-License-Identifier:\s*(.*?)\s*(?:\*/)?$", first_line)
            if identifier is None or identifier[1] not in licenses:
                raise ValueError(f"Unrecorded source license: {name}")
    gpl3 = read("LICENSES/GPL-3.0-or-later.txt").decode()
    if "Version 3, 29 June 2007" not in gpl3 or "END OF TERMS AND CONDITIONS" not in gpl3:
        raise ValueError("Complete GPLv3 text is required")
    for name in ("LICENSES/ISC.txt", "docs/upstream-mbelib-attribution.md"):
        if "Permission to use, copy, modify" not in read(name).decode():
            raise ValueError("ISC permission notice is missing")


def table(inventory):
    lines = ["<!-- SPDX-License-Identifier: GPL-2.0-or-later -->", "# Component inventory", "",
             "Generated from `components.json` by `tools/legal_inventory.py --write-table`.", "",
             "| Component | Declared license | Assessed license | Distributed in | Permission / GPL compatibility |",
             "| --- | --- | --- | --- | --- |"]
    for c in inventory["components"]:
        lines.append(f"| {c['id']} | {c['license_declared']} | {c['license_concluded']} | "
                     f"{', '.join(c['distribution']) or 'external only'} | "
                     f"{c['copyright_permission']} / {c['gpl_compatibility']} |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-table", action="store_true")
    args = parser.parse_args()
    inventory = load()
    # Include new working files for local preflight; release bundling uses Git objects.
    names = subprocess.check_output(["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
                                    cwd=ROOT).decode().split("\0")
    validate(inventory, lambda name: (ROOT / name).read_bytes(), filter(None, names))
    output = ROOT / "legal/component-table.md"
    expected = table(inventory)
    if args.write_table:
        output.write_text(expected, encoding="utf-8")
    elif not output.is_file() or output.read_text(encoding="utf-8") != expected:
        raise ValueError("Component table is stale; run tools/legal_inventory.py --write-table")
    print("Component evidence and notices checked; this is not legal clearance.")


if __name__ == "__main__":
    main()
