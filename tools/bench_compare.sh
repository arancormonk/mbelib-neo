#!/usr/bin/env bash
set -euo pipefail

# Build and run the synthesis benchmark in explicit scalar and SIMD modes for a quick comparison.

root_dir=$(cd "$(dirname "$0")/.." && pwd)
scalar_build="$root_dir/build/bench-compare-scalar"
simd_build="$root_dir/build/bench-compare-simd"
scalar_bench="$scalar_build/bench_synth"
simd_bench="$simd_build/bench_synth"

echo "== Configure scalar benchmark build =="
cmake --preset dev-release \
  -B "$scalar_build" \
  -DMBELIB_ENABLE_SIMD=OFF \
  -DMBELIB_BUILD_BENCHMARKS=ON >/dev/null

echo "== Build scalar benchmark =="
cmake --build "$scalar_build" -j --target bench_synth >/dev/null || cmake --build "$scalar_build" -j >/dev/null

echo "== Configure SIMD benchmark build =="
cmake --preset dev-release \
  -B "$simd_build" \
  -DMBELIB_ENABLE_SIMD=ON \
  -DMBELIB_BUILD_BENCHMARKS=ON >/dev/null

echo "== Build SIMD benchmark =="
cmake --build "$simd_build" -j --target bench_synth >/dev/null || cmake --build "$simd_build" -j >/dev/null

echo
echo "== Running scalar benchmark =="
"$scalar_bench"

echo
echo "== Running SIMD benchmark =="
"$simd_bench"
