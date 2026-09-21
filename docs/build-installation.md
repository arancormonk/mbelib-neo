# Build and Installation Policy

mbelib-neo uses CMake as its primary build and installation system.

## Standard Build Variables

CMake honors standard compiler and linker selection mechanisms. Users may pass
toolchain choices through environment variables or CMake cache variables,
including:

- `CC`
- `CFLAGS`
- `LDFLAGS`
- `CMAKE_C_COMPILER`
- `CMAKE_C_FLAGS`
- `CMAKE_EXE_LINKER_FLAGS`
- `CMAKE_SHARED_LINKER_FLAGS`

The project may add warning, sanitizer, SIMD, LTO, hardening, or platform flags
based on explicit CMake options, but should not discard user-supplied compiler
or linker settings.

The example binaries are built locally and are not installed. They accept no path arguments and support pipes: `dstar_encode < input.wav > output.dstar`
and `dstar_decode < output.dstar > output.wav`. Diagnostics go to stderr. The decoder buffers
PCM in memory before writing the WAV header and data to stdout.

## SIMD Target Selection

`MBELIB_ENABLE_SIMD=ON` uses SSE2 or NEON routines when the compiler target
provides those instructions, with scalar routines otherwise. The bundled PFFFT
library uses its own supported backends (including SSE and AltiVec); CMake
explicitly selects scalar FFT code when none is available, even with SIMD ON.
Warnings-as-errors remains enabled by default in either case.

On 32-bit ARM, the option preserves the toolchain's FPU selection. A VFP-only
target builds with scalar FFT code and a configure status message explaining
the fallback. To enable NEON on a compatible CPU, pass
`-DCMAKE_C_FLAGS=-mfpu=neon` along with `-DMBELIB_ENABLE_SIMD=ON`, preserving
any other compiler flags your toolchain needs. Targets that already enable NEON,
including AArch64, use it without additional flags. CMake does not change the
ARM floating-point ABI.

The existing 32-bit x86 behavior is different: enabling the option requests
SSE2 code generation and requires an SSE2-capable CPU. Leave the option OFF for
baseline i386 portability. SIMD selection is based on the compiler target,
not runtime CPU detection; deploy to CPUs that support the chosen instructions.
Measure scalar and SIMD performance on your own core with a representative
workload before choosing a configuration.

## Debug Information

The install rules do not strip binaries. If users request debug information
through their compiler flags or build type, the build and installation process
should preserve it.

## Non-Recursive Build Structure

The project uses one top-level CMake build graph with explicit source lists. It
does not rely on recursive make invocations to build cross-dependent
subdirectories.

## Repeatable Builds

To repeat a build from the same source tree, toolchain, options, and environment:

```sh
rm -rf build/repeatable
cmake -S . -B build/repeatable -DCMAKE_BUILD_TYPE=Release
cmake --build build/repeatable -j
ctest --test-dir build/repeatable --output-on-failure
```

For release packages, use the tagged release workflow or reproduce its documented
steps from `.github/workflows/ci.yml`.

Known limits:

- Floating-point optimization and SIMD choices can affect generated code and
  output across architectures.
- Binary identity across different compilers, operating systems, or CMake
  generators is not guaranteed.
- Release artifacts should be compared within the same source, toolchain,
  platform, build options, and packaging environment.

## Installation Conventions

Install with CMake:

```sh
cmake --install build/dev-release --prefix /usr/local
```

On POSIX systems, staged packaging can use `DESTDIR`:

```sh
DESTDIR="$PWD/pkgroot" cmake --install build/dev-release --prefix /usr
```

Uninstall from the same build directory:

```sh
cmake --build build/dev-release --target uninstall
```

## Developer Setup

For local development:

```sh
cmake --preset dev-debug
cmake --build --preset dev-debug -j
ctest --preset dev-debug --output-on-failure
```

Optional local hooks:

```sh
tools/install-git-hooks.sh
```
