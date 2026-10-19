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

The project may add warning, sanitizer, SIMD, LTO, hardening, host-tuning or
platform flags based on CMake options, but should not discard user-supplied
compiler or linker settings. A target ISA in `CMAKE_C_FLAGS` disables the host
tuning (see [host tuning](#host-tuning-mbelib_enable_native)).

The example binaries are built locally and are not installed. They accept no path arguments and support pipes: `dstar_encode < input.wav > output.dstar`
and `dstar_decode < output.dstar > output.wav`. Diagnostics go to stderr. The decoder buffers
PCM in memory before writing the WAV header and data to stdout.

## SIMD Target Selection

`MBELIB_ENABLE_SIMD=ON` uses SSE2/AVX2 or NEON routines when the compiler target
provides those instructions, with scalar routines otherwise. The option defaults
to ON for Release, RelWithDebInfo, and MinSizeRel builds of 64-bit targets, and
for 64-bit targets under a multi-config generator (Ninja Multi-Config, Visual
Studio, Xcode), which ignores `CMAKE_BUILD_TYPE`; it defaults to OFF for
single-config Debug builds and 32-bit targets. The downloaded PFFFT library uses its own supported
backends (including SSE and AltiVec); CMake explicitly selects scalar FFT code
when none is available, even with SIMD ON. Warnings-as-errors remains enabled by
default in either case.

### Run-time AVX2/FMA dispatch

A SIMD build for baseline x86-64 still runs on any SSE2 CPU and uses AVX2/FMA
where the CPU has it. CMake compiles a second copy of the kernels whose output
depends on the vector width (the voiced oscillators and the tone projection,
the voiced copy with IEEE semantics like its baseline) and of the exact integer
kernels that gain from AVX2 (the soft-decision ECC search and the TETRA
adaptive codebook) with `-mavx2 -mfma`, or `/arch:AVX2` for MSVC and clang-cl,
and the library picks one per call from a CPU probe cached on first use. The
copies are built for every configuration without host tuning; a configuration
whose compiler target already has AVX2 (also one set only in a
`CMAKE_C_FLAGS_<CONFIG>`) leaves them unused, as the preprocessor decides per
translation unit. The
probe requires the whole x86-64-v3 set (AVX2, FMA, BMI1/2, LZCNT, MOVBE, F16C)
and OS support for the YMM state. A build whose compiler target already has
AVX2 (`-march=native` on an AVX2 host, `-march=x86-64-v3`, `/arch:AVX2`) calls
its kernels directly with no check. The rest of the library, including the
encoders' codebook searches, PCM conversion and PFFFT, keeps the target's
instructions. 32-bit x86 builds, universal macOS builds (several
`CMAKE_OSX_ARCHITECTURES`) and non-x86 targets do not dispatch; MinGW GCC
dispatches only when its assembler accepts `-muse-unaligned-vector-move`
(binutils 2.38 or later), since GCC does not align the Windows x64 stack for
32-byte spills (MinGW and Cygwin). GCC's LTO link drops `-Wa` options when its
inputs disagree on them, so those AVX2 copies are always compiled without LTO
(`-fno-lto`), however LTO was enabled (`MBELIB_ENABLE_LTO`,
`CMAKE_INTERPROCEDURAL_OPTIMIZATION` or `-flto` in `CMAKE_C_FLAGS`), and keep
the option in the shipped objects, static library included. The copies are built as the object library
`mbe-avx2-kernels` inside both libraries; its target options carry everything
that depends on the configuration, since the Xcode generator cannot generate
configuration-dependent per-source options (configure refuses such source
properties). The configure log reports the decision (`AVX2 dispatch: ...`).

In Debug builds the AVX2 copy produces exactly the output of a Debug build for
`-march=x86-64-v3` (both round every multiply and add separately), so the
golden hashes cover both paths. In optimized builds the AVX2 copy contracts
multiply-adds into FMAs as a v3 build does.

### Host tuning (`MBELIB_ENABLE_NATIVE`)

`MBELIB_ENABLE_NATIVE` (default ON) tunes the libraries' Release,
RelWithDebInfo and MinSizeRel configurations for the CPU of the build host:
`-march=native` on x86 and `-mcpu=native` on ARM with GCC and Clang, and with
MSVC or clang-cl on an x86/x64 host the highest of `/arch:AVX512`, `/arch:AVX2`
and `/arch:AVX` that a configure-time CPUID/XGETBV probe finds the CPU and OS
support (cached for the build directory and compiler) and that the compiler
accepts; an ignored or rejected option (an older toolset without
`/arch:AVX512`) falls back to the next lower level. Debug configurations are never tuned,
so their golden hashes, sanitizers and coverage keep the baseline target. The
flags stay private to the library targets (and the tests that select their
kernels); they never reach the exported CMake targets or `libmbe-neo.pc`. In a
tuned configuration the AVX2 copies are not built; an AVX2 host calls the AVX2
kernels directly.

CMake tunes nothing, and says why, when it would be unsafe or would override
your choice: when cross-compiling or `CMAKE_C_COMPILER_TARGET` is set, for
Android, iOS or Emscripten targets (also an Apple device SDK selected only by
`CMAKE_OSX_SYSROOT`), when `CMAKE_OSX_ARCHITECTURES` names several
architectures or another than the host's, when the compiler command
(`CC="gcc -march=x86-64"`), `CMAKE_C_FLAGS` or a release-like
`CMAKE_C_FLAGS_<CONFIG>` already contains `-march=`, `-mcpu=`, `/arch:` or
`-arch` (distribution flags such as makepkg's `-march=x86-64` included, quoted
or not), when MSVC's `CL` or `_CL_` environment variable does, or when the
compiler rejects the option. MinGW and Cygwin GCC are tuned only together with
the assembler's `-muse-unaligned-vector-move`, for the same stack-alignment
reason as the dispatch; since that option reaches every unit and would not
survive an LTO link, such a tuned build compiles every unit without LTO,
whatever enabled it.

A binary built with the default runs only on CPUs with the build host's
instruction set. **Packagers and CI that ship binaries pass
`-DMBELIB_ENABLE_NATIVE=OFF`**, which keeps the compiler's default target
(x86-64 v1/SSE2 on x86-64) and, on x86-64, the run-time AVX2 dispatch. To
choose another target, pass it in `CMAKE_C_FLAGS`, which also turns the tuning
off.

### 32-bit targets

On 32-bit ARM, the option preserves the toolchain's FPU selection. A VFP-only
target builds with scalar FFT code and a configure status message explaining
the fallback. To enable NEON on a compatible CPU, pass
`-DCMAKE_C_FLAGS=-mfpu=neon` along with `-DMBELIB_ENABLE_SIMD=ON`, preserving
any other compiler flags your toolchain needs. Targets that already enable NEON,
including AArch64, use it without additional flags. CMake does not change the
ARM floating-point ABI.

The existing 32-bit x86 behavior is different: enabling the option requests
SSE2 code generation and requires an SSE2-capable CPU. Leave the option OFF for
baseline i386 portability. Measure scalar and SIMD performance on your own core
with a representative workload before choosing a configuration.

## Debug Information

The install rules do not strip binaries. If users request debug information
through their compiler flags or build type, the build and installation process
should preserve it.

## Non-Recursive Build Structure

PFFFT is fetched as a SHA-256-verified commit archive, or supplied with
`FETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT` for offline builds. Its top-level CMake
project is not configured; the build remains C99-only. See
[dependency management](dependencies.md) for cache, update and rollback steps.

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
