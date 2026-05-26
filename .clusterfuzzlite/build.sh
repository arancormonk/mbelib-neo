#!/bin/bash -eu

PROJECT_DIR="$SRC/mbelib-neo"
OBJ_DIR="$WORK/obj"
GEN_INCLUDE_DIR="$WORK/include"

mkdir -p "$OBJ_DIR" "$GEN_INCLUDE_DIR/mbelib-neo" "$OUT"
sed 's/@PROJECT_VERSION@/0.0.0-fuzz/' \
  "$PROJECT_DIR/include/mbelib-neo/version.h.in" \
  > "$GEN_INCLUDE_DIR/mbelib-neo/version.h"

read -r -a FUZZ_CFLAGS <<< "${CFLAGS:-}"
read -r -a FUZZ_CXXFLAGS <<< "${CXXFLAGS:-}"
read -r -a FUZZING_ENGINE <<< "${LIB_FUZZING_ENGINE:-}"

COMMON_CFLAGS=(
  -std=gnu99
  -DPFFFT_SIMD_DISABLE=1
  -I"$PROJECT_DIR/include"
  -I"$GEN_INCLUDE_DIR"
  -I"$PROJECT_DIR"
  -I"$PROJECT_DIR/src/internal"
  -I"$PROJECT_DIR/src/external/pffft"
)

PROJECT_WARNING_FLAGS=(
  -Wall
  -Wextra
  -Wpedantic
  -Werror
)

sources=(
  src/ambe/ambe3600x2400.c
  src/ambe/ambe3600x2450.c
  src/ambe/ambe_common.c
  src/ecc/ecc.c
  src/ecc/ecc_const.c
  src/imbe/imbe7100x4400.c
  src/imbe/imbe7200x4400.c
  src/core/mbelib.c
  src/core/mbe_adaptive.c
  src/core/mbe_unvoiced_fft.c
  src/external/pffft/pffft.c
  src/external/pffft/fftpack.c
)

objects=()
for source in "${sources[@]}"; do
  object="$OBJ_DIR/${source//\//_}.o"
  "$CC" "${FUZZ_CFLAGS[@]}" "${COMMON_CFLAGS[@]}" "${PROJECT_WARNING_FLAGS[@]}" -c "$PROJECT_DIR/$source" -o "$object"
  objects+=("$object")
done

build_fuzzer() {
  local source=$1
  local name
  name=$(basename "$source" .cc)

  "$CXX" "${FUZZ_CXXFLAGS[@]}" "${PROJECT_WARNING_FLAGS[@]}" -std=c++17 \
    -I"$PROJECT_DIR/include" \
    -I"$GEN_INCLUDE_DIR" \
    "$PROJECT_DIR/$source" \
    "${objects[@]}" \
    "${FUZZING_ENGINE[@]}" \
    -lm \
    -o "$OUT/$name"
}

fuzzers=(
  fuzz/fuzz_process_frame.cc
  fuzz/fuzz_frame_decode.cc
)

for fuzzer in "${fuzzers[@]}"; do
  build_fuzzer "$fuzzer"
done
