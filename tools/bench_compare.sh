#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
set -euo pipefail

# Compare synthesis, encoding and PCM conversion with matching build options.

root_dir=$(cd "$(dirname "$0")/.." && pwd)
scalar_build="$root_dir/build/bench-compare-scalar"
simd_build="$root_dir/build/bench-compare-simd"
benchmarks=(bench_synth bench_unvoiced bench_convert bench_encode bench_replay)

echo "== Configure scalar benchmark build =="
cmake --preset dev-release \
  -B "$scalar_build" \
  -DMBELIB_ENABLE_SIMD=OFF \
  -DMBELIB_BUILD_BENCHMARKS=ON > /dev/null

echo "== Build scalar benchmark =="
cmake --build "$scalar_build" -j --target "${benchmarks[@]}" > /dev/null

echo "== Configure SIMD benchmark build =="
cmake --preset dev-release \
  -B "$simd_build" \
  -DMBELIB_ENABLE_SIMD=ON \
  -DMBELIB_BUILD_BENCHMARKS=ON > /dev/null

echo "== Build SIMD benchmark =="
cmake --build "$simd_build" -j --target "${benchmarks[@]}" > /dev/null

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
echo "== Running SIMD benchmark =="
run_benchmarks "$simd_build"

echo
echo "Replay real speech and hard/soft frames with bench_replay in these same trees; see docs/testing.md."
