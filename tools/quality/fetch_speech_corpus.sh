#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
# Run from the repository root. Requires Bash, git, coreutils, and network access.
set -euo pipefail

if [[ $# -ne 0 || ! -x tools/fetch-pinned-git.sh ]]; then
    echo "usage: tools/quality/fetch_speech_corpus.sh (from repository root)" >&2
    exit 2
fi

source_dir=build/quality/codec2-src
corpus_dir=build/quality/corpus
files=(hts1a.raw hts2a.raw hts1.raw kristoff.raw ve9qrp_10s.raw)
sizes=(48000 48000 96000 80000 160000)

tools/fetch-pinned-git.sh https://github.com/drowe67/codec2.git \
    310777b1c6f1af0bc7c72f5b32f80f6fd9136962 "$source_dir"

# Check all pinned assets before replacing any existing corpus file.
for i in "${!files[@]}"; do
    path="$source_dir/raw/${files[i]}"
    if [[ ! -f "$path" || ! -r "$path" ]]; then
        echo "missing corpus asset: $path" >&2
        exit 1
    fi
    actual=$(stat -c %s -- "$path")
    if [[ "$actual" -ne "${sizes[i]}" ]]; then
        echo "unexpected corpus size: $path ($actual bytes; expected ${sizes[i]})" >&2
        exit 1
    fi
done

mkdir -p "$corpus_dir"
tmp=
trap 'if [[ -n "$tmp" ]]; then rm -f -- "$tmp"; fi' EXIT
for file in "${files[@]}"; do
    tmp=$(mktemp "$corpus_dir/.corpus.XXXXXX")
    cp -- "$source_dir/raw/$file" "$tmp"
    mv -f -- "$tmp" "$corpus_dir/$file"
    tmp=
    printf '%s\n' "$corpus_dir/$file"
done
