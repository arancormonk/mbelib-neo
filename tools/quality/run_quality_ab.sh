#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
# Linux only: swap libmbe-neo.so.2 using LD_LIBRARY_PATH, with one evaluator.
# Run from the repository root. Requires Bash, coreutils, and Python 3 (tooling
# only), plus the built encoder and quality evaluator. Gain is calibrated only
# against the baseline. Frame caches are bound to corpus, encoder, evaluator,
# baseline library, and calibration code hashes; both variants use identical bits.
set -euo pipefail

fail() {
    echo "$*" >&2
    exit 2
}

[[ $# -ge 2 && $# -le 3 ]] || fail \
    "usage: $0 <baseline_libdir> <candidate_libdir> [corpus_dir=build/quality/corpus]"
[[ $(uname -s) == Linux ]] || fail "this runner requires Linux (LD_LIBRARY_PATH)"
command -v python3 >/dev/null || fail "Python 3 is required for JSON/table rendering"
encoder=build/quality/op25_encode
evaluator=build/dev-debug/mbe_quality_eval
[[ -x "$encoder" && -x "$evaluator" ]] || fail \
    "run from repository root after building $encoder and $evaluator"

for dir in "$1" "$2"; do
    [[ -d "$dir" && -r "$dir/libmbe-neo.so.2" && -f "$dir/libmbe-neo.so.2" ]] || fail \
        "expected a library directory containing readable libmbe-neo.so.2: $dir"
done
baseline=$(realpath -e -- "$1")
candidate=$(realpath -e -- "$2")
[[ "$baseline" != *:* && "$candidate" != *:* ]] || fail \
    "library directories must not contain ':' (LD_LIBRARY_PATH separator)"
corpus=${3:-build/quality/corpus}
[[ -d "$corpus" && -r "$corpus" && -x "$corpus" ]] || fail "unreadable corpus directory: $corpus"
corpus=$(realpath -e -- "$corpus")
shopt -s nullglob
raw_files=("$corpus"/*.raw)
[[ ${#raw_files[@]} -gt 0 ]] || fail "no *.raw files in corpus: $corpus"
for raw in "${raw_files[@]}"; do
    [[ -f "$raw" && -r "$raw" && -s "$raw" ]] || fail "empty or unreadable raw file: $raw"
done

mkdir -p build/quality/frames build/quality/out
reports=()
for raw in "${raw_files[@]}"; do
    name=${raw##*/}
    name=${name%.raw}
    for mode in imbe7200 ambe2450 ambe2400; do
        frames="build/quality/frames/$name.$mode.txt"
        python3 tools/quality/calibrate_encoder.py \
            "$encoder" "$evaluator" "$baseline" "$raw" "$mode" "$frames"
        prefix="build/quality/out/$name.$mode"
        for variant in baseline candidate; do
            libdir=$baseline
            [[ "$variant" != candidate ]] || libdir=$candidate
            LD_LIBRARY_PATH="$libdir" "$evaluator" --codec "$mode" \
                --frames "$frames" --out "$prefix.$variant.wav" \
                --ref "$raw" --json "$prefix.$variant.json" >&2
        done
        reports+=("$name" "$mode" "$prefix.baseline.json" "$prefix.candidate.json")
    done
done

python3 - "${reports[@]}" <<'PY'
import json
import math
import sys

metrics = [
    ("lsd_db", "lsd_db"),
    ("env_corr", "env_corr"),
    ("crest_delta_db", "crest_delta_db"),
    ("boundary_index_db", "boundary_index_db"),
    ("join_excess_db", "join_excess_db"),
    ("band_delta_db_0_500", "band_0_500"),
    ("band_delta_db_500_1000", "band_500_1000"),
    ("band_delta_db_1000_2000", "band_1000_2000"),
    ("band_delta_db_2000_3000", "band_2000_3000"),
    ("band_delta_db_3000_4000", "band_3000_4000"),
]
rows = []
groups = {mode: [] for mode in ("imbe7200", "ambe2450", "ambe2400")}
for offset in range(1, len(sys.argv), 4):
    name, mode, baseline_path, candidate_path = sys.argv[offset:offset + 4]
    pair = []
    for path in (baseline_path, candidate_path):
        with open(path, encoding="utf-8") as stream:
            record = json.load(stream)
        for key, _ in metrics:
            value = record.get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise SystemExit(f"missing or non-finite {key} in {path}")
        if not isinstance(record.get("pcm_fnv1a"), str):
            raise SystemExit(f"missing pcm_fnv1a string in {path}")
        pair.append(record)
    groups[mode].append(pair)
    rows.append([name, mode] + [
        f"{pair[0][key]:.3f}/{pair[1][key]:.3f}" for key, _ in metrics
    ] + [f"{pair[0]['pcm_fnv1a']}/{pair[1]['pcm_fnv1a']}"])

mean_rows = []
for mode, pairs in groups.items():
    mean_rows.append(["MEAN", mode] + [
        "/".join(f"{sum(pair[side][key] for pair in pairs) / len(pairs):.3f}"
                 for side in (0, 1))
        for key, _ in metrics
    ] + ["-"])

headers = ["file", "mode"] + [label for _, label in metrics] + ["pcm_fnv1a"]
widths = [max(len(row[column]) for row in [headers] + rows + mean_rows)
          for column in range(len(headers))]

def render(row):
    return " | ".join(value.ljust(width) for value, width in zip(row, widths))

print("Each cell: baseline/candidate. Band columns are band_delta_db_<range>.")
print("Per-mode means are unweighted across corpus files; hashes are not averaged.")
print(render(headers))
print("-+-".join("-" * width for width in widths))
for row in rows:
    print(render(row))
print("-+-".join("-" * width for width in widths))
for row in mean_rows:
    print(render(row))
PY
