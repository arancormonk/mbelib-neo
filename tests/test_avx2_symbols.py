#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check that the AVX2 translation units export only *_avx2 symbols.

A baseline caller must never reach code compiled for AVX2 except through the
run-time dispatch. A non-static (or inline, COMDAT) definition in an AVX2 unit
could otherwise be the copy the linker keeps for every caller. Reads the static
library's members with nm. With --expect-code, every AVX2 unit must define at
least one *_avx2 function (the configuration dispatches); without it, the
units may be empty. Exits 77 (skip) for an LTO configuration, whose objects
hold compiler IR rather than final symbols.
"""

import argparse
import subprocess
import sys

AVX2_UNITS = (
    "mbe_voiced_avx2.c",
    "mbe_tone_fit_avx2.c",
    "ecc_soft_search_avx2.c",
    "tetra_acelp_avx2.c",
)
# Global symbol types: text, data, read-only, BSS, common, weak, unique.
GLOBAL_TYPES = set("TDRBCGSVWiu")


def unit_of(member):
    for unit in AVX2_UNITS:
        if member.startswith(unit):
            return unit
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nm")
    parser.add_argument("archive")
    parser.add_argument("--expect-code", default="0")
    parser.add_argument("--lto", default="0")
    args = parser.parse_args()
    if args.lto == "1":
        print("LTO configuration: object files hold compiler IR; skipping")
        return 77
    output = subprocess.run([args.nm, "-A", "--defined-only", args.archive], capture_output=True, text=True,
                            check=True).stdout
    defined = {unit: [] for unit in AVX2_UNITS}
    seen = set()
    for line in output.splitlines():
        # archive:member:  value type name
        head, _, rest = line.rpartition(":")
        member = head.rsplit(":", 1)[-1] if ":" in head else head
        unit = unit_of(member)
        if unit is None:
            continue
        seen.add(unit)
        fields = rest.split()
        if len(fields) >= 2 and fields[-2] in GLOBAL_TYPES:
            defined[unit].append((fields[-2], fields[-1]))
    missing = [unit for unit in AVX2_UNITS if unit not in seen and args.expect_code == "1"]
    failures = []
    for unit in AVX2_UNITS:
        for kind, name in defined[unit]:
            if not name.endswith("_avx2") or kind not in "TDRB":
                failures.append(f"{unit}: {kind} {name}")
        if args.expect_code == "1" and not any(kind == "T" for kind, _ in defined[unit]):
            failures.append(f"{unit}: defines no *_avx2 function in a dispatching configuration")
    if missing:
        failures.append(f"members not found in {args.archive}: {', '.join(missing)}")
    for unit in AVX2_UNITS:
        print(f"{unit}: {' '.join(name for _, name in defined[unit]) or '(empty)'}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
