#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
set -euo pipefail

# Run from the repository root. OP25 is GPL-3.0-or-later; its bundled mbelib
# stays in this standalone executable and is never linked with mbelib-neo.
# Requires Bash, git, patch, coreutils, gcc/g++ (C++17), and network access.
# Each invocation refetches the pinned checkout and rebuilds the encoder.
if [[ ! -f tools/quality/op25_encode_driver.cc || ! -x tools/fetch-pinned-git.sh ]]; then
  echo "Run $0 from the repository root." >&2
  exit 2
fi

tools/fetch-pinned-git.sh https://github.com/boatbod/op25.git \
  71abcd0ead32f86f51615ea6cc8a6a4dba4c949a build/quality/op25-src

lib=build/quality/op25-src/op25/gr-op25_repeater/lib
objdir=build/quality/op25-obj
mkdir -p "$objdir"
# Adapt disposable copies only; preserve the pinned upstream checkout verbatim.
# This is a mbelib-neo quality bridge, not an interoperable D-STAR transmitter.
cp "$lib/ambe_encoder.cc" "$lib/p25p2_vf.cc" "$objdir/"
patch --batch --fuzz=0 -d "$objdir" -p0 < tools/quality/op25_mbelib_gain.patch
objects=()
for unit in ambe mbelib; do
  object="$objdir/$unit.o"
  gcc -O2 -I"$lib" -c "$lib/$unit.c" -o "$object"
  objects+=("$object")
done

# The complete imbe_vocoder/CMakeLists.txt manifest at the pinned revision.
vocoder_units=(
  aux_sub basicop2 ch_decode ch_encode dc_rmv decode dsp_sub encode
  imbe_vocoder math_sub pe_lpf pitch_est pitch_ref qnt_sub rand_gen
  sa_decode sa_encode sa_enh tbls uv_synt v_synt v_uv_det
)
sources=()
for unit in "${vocoder_units[@]}"; do
  sources+=("$lib/imbe_vocoder/$unit.cc")
done
sources+=("$objdir/ambe_encoder.cc" "$objdir/p25p2_vf.cc" "$lib/rs.cc"
  tools/quality/op25_encode_driver.cc)
for source in "${sources[@]}"; do
  unit=${source##*/}
  object="$objdir/${unit%.cc}.o"
  g++ -std=c++17 -O2 -I"$lib" -I"$lib/imbe_vocoder" -c "$source" -o "$object"
  objects+=("$object")
done

g++ "${objects[@]}" -lm -o build/quality/op25_encode
printf 'Built %s\n' build/quality/op25_encode
