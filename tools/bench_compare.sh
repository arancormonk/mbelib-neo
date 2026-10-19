#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
set -euo pipefail

# Compare synthesis, encoding and PCM conversion with matching build options:
# scalar and SIMD portable builds (MBELIB_ENABLE_NATIVE=OFF; on x86-64 the SIMD
# build runs its AVX2 kernel copies on AVX2 CPUs), and a SIMD build tuned for
# this host (the MBELIB_ENABLE_NATIVE default).

root_dir=$(cd "$(dirname "$0")/.." && pwd)
scalar_build="$root_dir/build/bench-compare-scalar"
simd_build="$root_dir/build/bench-compare-simd"
native_build="$root_dir/build/bench-compare-native"
benchmarks=(bench_synth bench_unvoiced bench_convert bench_encode bench_replay bench_tetra)

configure_and_build() {
  local label="$1" build_dir="$2"
  shift 2
  echo "== Configure $label benchmark build =="
  cmake --preset dev-release -B "$build_dir" -DMBELIB_BUILD_BENCHMARKS=ON "$@" > /dev/null
  echo "== Build $label benchmark =="
  cmake --build "$build_dir" -j --target "${benchmarks[@]}" > /dev/null
}

configure_and_build scalar "$scalar_build" -DMBELIB_ENABLE_SIMD=OFF -DMBELIB_ENABLE_NATIVE=OFF
configure_and_build "SIMD (portable)" "$simd_build" -DMBELIB_ENABLE_SIMD=ON -DMBELIB_ENABLE_NATIVE=OFF
configure_and_build "SIMD (host-tuned)" "$native_build" -DMBELIB_ENABLE_SIMD=ON -DMBELIB_ENABLE_NATIVE=ON

run_benchmarks() {
  local build_dir="$1"
  local workload voicing
  for workload in steady gradual jump; do
    for voicing in mixed voiced; do
      "$build_dir/bench_synth" "$workload" "$voicing"
    done
  done
  "$build_dir/bench_unvoiced" 36 20000 7
  "$build_dir/bench_convert" 10000000
  "$build_dir/bench_encode"
}

echo
echo "== Running scalar benchmark =="
run_benchmarks "$scalar_build"

echo
echo "== Running SIMD (portable) benchmark =="
run_benchmarks "$simd_build"

echo
echo "== Running SIMD (host-tuned) benchmark =="
run_benchmarks "$native_build"

echo
echo "Replay speech and frames with bench_replay, or TETRA CHD files with bench_tetra, in these same trees; see docs/testing.md."
