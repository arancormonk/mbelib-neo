# mbelib-neo

Performance‑enhanced IMBE/AMBE vocoder primitives with a modern CMake build, installable headers, and pkg-config/CMake package integration.

Project homepage: https://github.com/arancormonk/mbelib-neo

[![CI](https://github.com/arancormonk/mbelib-neo/actions/workflows/ci.yml/badge.svg)](https://github.com/arancormonk/mbelib-neo/actions/workflows/ci.yml)

## Patent Notice

```
This source code is provided for educational purposes only. It
is a written description of how certain voice encoding/decoding
algorithms could be implemented. Executable objects compiled or
derived from this package may be covered by one or more patents.
Readers are strongly advised to check for any patent restrictions
or licensing requirements before compiling or using this source code.
```

This notice is advisory and does not modify the license. See `LICENSE` for terms.

## License

- mbelib-neo is distributed under GPL-2.0-or-later (`LICENSE`). All compiled binaries and source distributions fall under that license.
- Portions originate from the original ISC-licensed mbelib; the retained ISC notice is in `COPYRIGHT` for attribution and clarity. Those contributions are redistributed under GPL-2.0-or-later as permitted by the ISC terms.
- Source files carry SPDX identifiers reflecting their license.

## Overview

- A performance‑enhanced fork of [lwvmobile/mbelib](https://github.com/lwvmobile/mbelib), which is a fork of [szechyjs/mbelib](https://github.com/szechyjs/mbelib)
- Supports IMBE 7200x4400 (P25 Phase 1), IMBE 7100x4400 (ProVoice), AMBE (D‑STAR), and AMBE+2 (DMR, NXDN, P25 Phase 2, dPMR, etc.).
- Stable public API in `#include <mbelib-neo/mbelib.h>` with version macro `MBELIB_VERSION`.
- Ships as shared and static libraries: `libmbe-neo.{so|dylib}` on Unix-like platforms, and `mbe-neo.dll` + import lib / `mbe-neo-static.lib` on Windows.
- Installable CMake package (`mbe_neo::mbe_shared` / `mbe_neo::mbe_static`) and pkg-config file (`libmbe-neo`).

## Build From Source

Requirements

- C compiler with C99 support.
- CMake ≥ 3.20.
- On non‑Windows platforms, links against `-lm` automatically.

Using CMake presets (recommended)

```
# From the repository root:

# Debug build with tests/examples
cmake --preset dev-debug
cmake --build --preset dev-debug -j
ctest --preset dev-debug -V

# Release build (SIMD + fast-math + LTO)
cmake --preset dev-release
cmake --build --preset dev-release -j

# Optional sanitizers (Linux/Clang/GCC)
cmake --preset asan-ubsan-debug
cmake --build --preset asan-ubsan-debug -j
ctest --preset asan-ubsan-debug -V

# Disable tone synthesis (AMBE tones)
cmake --preset notones-debug
cmake --build --preset notones-debug -j
ctest --preset notones-debug -V
```

Notes

- Presets create out-of-source builds under `build/<preset>/`. Run the above from the repo root.
- Alternatively, you can run from anywhere using `cmake -S <repo> -B <builddir>`.

Manual configure/build

```
# From the repository root (single-config generators like Ninja/Make):
cmake -S . -B build/manual-debug -DCMAKE_BUILD_TYPE=Debug -DMBELIB_BUILD_TESTS=ON -DMBELIB_BUILD_EXAMPLES=ON
cmake --build build/manual-debug -j
ctest --test-dir build/manual-debug -V

# From anywhere (equivalent):
cmake -S <repo-path> -B <build-dir> -DCMAKE_BUILD_TYPE=Debug -DMBELIB_BUILD_TESTS=ON -DMBELIB_BUILD_EXAMPLES=ON
cmake --build <build-dir> -j
ctest --test-dir <build-dir> -V
```

## Install / Uninstall

```
# Single-config generators (Unix Makefiles/Ninja):
cmake --install build/dev-release

# Multi-config generators (Visual Studio/Xcode):
cmake --install build/dev-release --config Release

# Uninstall from the same build directory
cmake --build build/dev-release --target uninstall
```

## Configuration Options

- `-DNOTONES=ON` — Disable AMBE/AMBE+2 tone synthesis (adds `-DDISABLE_AMBE_TONES`).
- `-DMBELIB_BUILD_TESTS=ON` — Build CTest executables (default ON).
- `-DMBELIB_BUILD_EXAMPLES=ON` — Build `examples/` (default ON).
- `-DMBELIB_ENABLE_WARNINGS=ON` — Enable common warnings (default ON).
- `-DMBELIB_WARNINGS_AS_ERRORS=ON` — Treat warnings as errors.
- `-DMBELIB_ENABLE_ASAN=ON` — Enable AddressSanitizer in Debug builds.
- `-DMBELIB_ENABLE_UBSAN=ON` — Enable UndefinedBehaviorSanitizer in Debug builds.
- `-DMBE_ENABLE_DEBUG_LOGS=ON` — Verbose debug logging in codec sources.
- `-DMBELIB_BUILD_DOCS=ON` — Add `docs` target (requires Doxygen).
- `-DMBELIB_ENABLE_FAST_MATH=ON` — Enable fast-math (`-ffast-math`/`/fp:fast`) on library targets.
- `-DMBELIB_ENABLE_LTO=ON` — Enable IPO/LTO in Release builds when supported.
- `-DMBELIB_ENABLE_SIMD=ON` — Enable SIMD-accelerated routines (SSE2 on x86/x86_64, NEON on ARM64) in hot paths (including float→int16 conversion and unvoiced FFT synthesis loops). Falls back to portable scalar code when unavailable.
- `-DMBELIB_STRICT_ORDER=ON` — Compatibility flag that currently only defines `MBELIB_STRICT_ORDER` at compile time (no behavior change in this tree).
- Note: the `dev-release` preset enables SIMD, fast-math, and LTO by default when supported.
- `-DMBELIB_BUILD_BENCHMARKS=ON` — Build optional local micro‑benchmarks (not run in CI): `bench_synth`, `bench_unvoiced`, and `bench_convert`.

## Using The Library

Headers and linking

- Public header: `#include <mbelib-neo/mbelib.h>`
- pkg-config: `pkg-config --cflags --libs libmbe-neo`
- CMake package:
  ```cmake
  find_package(mbe-neo CONFIG REQUIRED)
  target_link_libraries(your_target PRIVATE mbe_neo::mbe_shared) # or mbe_neo::mbe_static
  ```

Minimal example

```c
#include <stdio.h>
#include <mbelib-neo/mbelib.h>

int main(void) {
  char ver[32] = {0};
  mbe_printVersion(ver);
  printf("mbelib version: %s\n", ver);
  // Or: puts(mbe_versionString());
  return 0;
}
```

You can also build and run the bundled example:

```
cmake --build build/dev-debug --target example_print_version
./build/dev-debug/example_print_version
```

### Frame/Data API Quick Reference

| Codec path | Frame API input | Data API input | Scratch/output bits |
| --- | --- | --- | --- |
| AMBE 3600x2400 | `char ambe_fr[4][24]` | `char ambe_d[49]` | `ambe_d[49]` |
| AMBE 3600x2450 | `char ambe_fr[4][24]` | `char ambe_d[49]` | `ambe_d[49]` |
| IMBE 7200x4400 | `char imbe_fr[8][23]` | `char imbe_d[88]` | `imbe_d[88]` |
| IMBE 7100x4400 | `char imbe_fr[7][24]` | n/a (frame-only public path) | `imbe_d[88]` (converted to 7200 layout) |

Use `mbe_process*Frame*` when you have interleaved codec frames, and `mbe_process*Data*` when you already have decoded parameter bits.

### Stateful Decode Workflow

- Keep one `mbe_parms` state triplet per audio stream/thread: `cur_mp`, `prev_mp`, and `prev_mp_enhanced`.
- Initialize once before decoding with `mbe_initMbeParms(&cur_mp, &prev_mp, &prev_mp_enhanced)`.
- Prefer `mbe_process*Frame*` APIs when you have raw vocoder frames. These paths run demod/ECC and write `errs`/`errs2`.
- Use `mbe_process*Data*` APIs only when you already have unpacked parameter bits:
  - `errs` and `errs2` are caller-provided inputs for these parameter-only paths.
  - `mbe_processAmbe2450Data*` and `mbe_processImbe4400Data*` currently ignore `*errs`, but require a valid pointer for API compatibility.
- `err_str` is a compact status trace: `'='` repeated `*errs2` times, optional suffix markers (`E`, `T`, `R`, `M` depending on codec/path), then a trailing NUL.
- Size `err_str` for at least `(*errs2 + 3)` bytes so there is room for up to two suffix markers plus NUL (a conservative `char err_str[128]` works well for public entry points).
- `mbe_printVersion(char *str)` uses a legacy fixed write width of 32 bytes. Pass a buffer of at least 32 bytes, or prefer `mbe_versionString()` when possible.

### Audio Sample Scaling (Float vs int16)

- The library synthesizes 8 kHz audio in 160-sample frames.
- The `short` entry points write 16-bit PCM with soft clipping (~95% full-scale).
- The `*f` entry points return mbelib’s historical float scale (not normalized `[-1, +1]`). `mbe_floattoshort()` applies the same `* 7.0` scaling and clipping used by the `short` APIs.
- To feed a normalized float pipeline, scale each float sample by `(7.0f / 32768.0f)` (range is approximately `[-0.95, +0.95]` after soft clipping).

## Windows (MSVC) Quickstart

- Open the "x64 Native Tools Command Prompt for VS".
- From the repository root:

```
cmake --preset dev-release
cmake --build --preset dev-release --config Release -j
ctest --preset dev-release -C Release
cmake --install build/dev-release --config Release
```

Consuming with CMake on Windows (link static to avoid DLL path issues):

```cmake
cmake_minimum_required(VERSION 3.20)
project(consumer C)
find_package(mbe-neo CONFIG REQUIRED)
add_executable(consumer consumer.c)
target_link_libraries(consumer PRIVATE mbe_neo::mbe_static) # or mbe_neo::mbe_shared
```

## pkg-config on Windows (MSYS2/MinGW)

- Install MSYS2 and open a MinGW64 shell (e.g., UCRT64 or MINGW64).
- Install toolchain and pkg-config:

```
pacman -S --needed mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake mingw-w64-x86_64-pkg-config
```

- Configure, build, and install using the MinGW generator (example to the default /mingw64 prefix):

```
cmake -G "MinGW Makefiles" -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build --prefix /mingw64
```

- Compile a consumer using pkg-config:

```
cat > consumer.c << 'EOF'
#include <stdio.h>
#include <mbelib-neo/mbelib.h>
int main(void) { char ver[32]={0}; mbe_printVersion(ver); printf("%s\n", ver); }
EOF
cc consumer.c $(pkg-config --cflags --libs libmbe-neo) -o consumer.exe
```

Notes

- If you install to a custom prefix, set `PKG_CONFIG_PATH="<prefix>/lib/pkgconfig:$PKG_CONFIG_PATH"`.
- For shared builds, ensure `<prefix>/bin` is on `PATH` (so the DLL is found) or copy the DLL next to your `.exe`.
- To prefer static linking via pkg-config, you can try `pkg-config --static --libs libmbe-neo` (requires static libs available) or use the CMake package and link `mbe_neo::mbe_static`.

## Audio Quality Improvements

mbelib-neo implements JMBE-compatible audio synthesis algorithms for improved audio quality:

- **FFT-based unvoiced synthesis** (Algorithms #117-126): Uses 256-point FFT with band-level scaling instead of the legacy oscillator bank approach. Provides cleaner, more natural unvoiced sounds with proper spectral shaping.

- **Weighted Overlap-Add (WOLA)** (Algorithm #126): Smooth frame transitions for unvoiced synthesis using a 211-element trapezoidal synthesis window. Eliminates audible discontinuities between frames.

- **Adaptive smoothing** (Algorithms #111-116): Error-rate-based parameter smoothing that gracefully handles corrupted frames. Includes local energy tracking, adaptive voicing thresholds, and amplitude scaling.

- **Voiced phase/amplitude interpolation** (Algorithms #134-138): Smooth interpolation of pitch and amplitude for low-frequency harmonics during stable pitch periods. Reduces "buzzy" artifacts in voiced speech.

- **Codec-specific frame repeat/muting parity**: Matches JMBE behavior where IMBE mutes on max repeats or error-rate threshold, while AMBE muting is repeat-driven in the synth path. IMBE prolonged repeat headroom resets to a default model state.

- **AMBE tone/erasure fallback parity**: AMBE 3600x2450 tone classification uses JMBE-style verification plus BER gating (`errorCount < 6`). Unverified/high-BER tone candidates fall back to erasure/repeat behavior, and invalid tone IDs reuse prior voice model until repeat max is reached.

- **LCG noise generator with buffer overlap**: JMBE-compatible Linear Congruential Generator for deterministic noise, with 96-sample overlap for smooth continuity between frames.

- **Comfort-noise parity model**: Muted-frame noise follows JMBE’s low-level uniform white-noise model (`0.003` gain semantics) using Java `Random`-compatible per-thread RNG behavior.

These improvements bring mbelib-neo's audio quality closer to the reference JMBE (Java Multi-Band Excitation) implementation.

## API Notes

- Public API is prefixed `mbe_` and declared in `mbelib.h`.
- Float and 16‑bit PCM variants are provided (e.g., `mbe_processAmbe3600x2400Framef` and `mbe_processAmbe3600x2400Frame`).
- Version macro `MBELIB_VERSION` is defined in the generated header `mbelib-neo/version.h` and returned by `mbe_printVersion`.
- `mbe_versionString()` returns a const pointer to the version string.
- Unvoiced synthesis uses a JMBE-style LCG state carried in `mbe_parms` (`noiseSeed`/`noiseOverlap`), with per-thread cold-start seeding via `mbe_setThreadRngSeed(uint32_t)`.
- ABI: shared library `SOVERSION` follows the project major version; minor updates aim to remain ABI compatible.

### Determinism & RNG

- `mbe_setThreadRngSeed(seed)` seeds both thread-local comfort-noise RNG state and the next unvoiced LCG cold start. For reproducible output, set it per thread before synthesis.
- `mbe_setThreadRngSeed(0)` is accepted and remapped internally to a non-zero seed. Use an explicit non-zero seed when exact reproducibility matters across builds.
- Unvoiced noise progression after cold start is driven by the per-frame LCG state in `mbe_parms`, so deterministic playback requires carrying frame state forward consistently.
- Runtime synthesis helpers (RNG state and FFT plan) are thread-local; do not share mutable decode state (`mbe_parms`) across threads unless you synchronize externally.
- AMBE/IMBE frame handling intentionally differs for JMBE parity: IMBE uses error-rate muting and repeat-headroom reset behavior, while AMBE tone/erasure/repeat transitions follow AMBE-specific JMBE gating rules.
- Enabling `MBELIB_ENABLE_SIMD=ON` selects vectorized math on supported CPUs. This can change
  floating‑point rounding at the bit level. Tests enforce exactness for int16 on x86 in Debug, and
  sanity bounds elsewhere.

## Benchmarks (Optional)

Build micro-benchmarks:

```
cmake --preset dev-release-simd -DMBELIB_BUILD_BENCHMARKS=ON
cmake --build --preset dev-release-simd -j --target bench_synth bench_unvoiced bench_convert
```

Run benchmark executables:

```
./build/dev-release-simd/bench_synth
./build/dev-release-simd/bench_unvoiced
./build/dev-release-simd/bench_convert
```

Quick scalar-vs-SIMD comparison helper:

```
tools/bench_compare.sh
```

## Tests and Examples

- Run tests with `ctest --preset dev-debug -V` (or `ctest -V` from the build directory).
- Included tests: `test_api` (version/headers), `test_ecc` (Golay and Hamming), `test_noise_determinism` (unvoiced RNG/frame-state determinism), `test_params` (parameter/synthesis behavior), `test_golden_pcm` (golden hash regression checks).
- Example: `examples/print_version.c` shows linking and header usage.
- Golden hash helper: `gen_golden` (available when `MBELIB_BUILD_TESTS=ON`, default) prints current FNV-1a reference values for synthesis/conversion regression workflows (`cmake --build build/dev-debug --target gen_golden && ./build/dev-debug/gen_golden`).

## Documentation

Optional Doxygen documentation can be generated:

```
cmake --preset dev-debug -DMBELIB_BUILD_DOCS=ON
cmake --build --preset dev-debug --target docs
# Output in build/docs/html
```

The generated site focuses on the public API (`include/`) and bundled examples.

## Project Layout

- Public headers: `include/mbelib-neo/`
- Sources: `src/core/`, `src/ecc/`, `src/ambe/`, `src/imbe/`
- Internal headers: `src/internal/`
- External libraries: `src/external/pffft/` (bundled PFFFT + FFTPACK sources, BSD-like license)
- Tests: `tests/` • Examples: `examples/` • Benchmarks: `bench/`
- Developer tooling and CI parity scripts: `tools/`

## Contributing

- Follow `.clang-format` (LLVM style, 4‑space indent, 120 cols). You can run `tools/format.sh`.
- Static analysis scripts (CI-aligned):
  - `tools/clang_tidy.sh` (supports `--strict` and targeted TUs).
  - `tools/cppcheck.sh` (supports `--strict` and targeted TUs).
  - `tools/iwyu.sh` (include hygiene via include-what-you-use; excludes `src/external`).
  - `tools/gcc_fanalyzer.sh` (GCC `-fanalyzer` diagnostics; excludes `src/external`).
  - `tools/scan_build.sh` (Clang Static Analyzer via `scan-build`; excludes `src/external`).
  - `tools/semgrep.sh` (additional SAST rules; use `--strict` to fail on findings).
- Git hooks: `tools/install-git-hooks.sh` enables auto-format on commit and CI-style pre-push checks.
- Optional full scan-build pre-push/preflight pass: set `MBE_HOOK_RUN_SCAN_BUILD=1`.
- Manual preflight runner: `tools/preflight_ci.sh` runs the same pre-push checks without pushing.
- Prefer keeping internal symbols `static` and declarations in headers where shared.
- Before sending changes: build locally, run `ctest -V`, and ensure examples still link.
