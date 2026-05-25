# Architecture

mbelib-neo is a C library for IMBE/AMBE vocoder primitives. It is organized so
public consumers use installed headers while codec internals remain private to
the library.

## Main Actors

- Library consumer: application code that links `mbe_neo::mbe_shared`,
  `mbe_neo::mbe_static`, or `libmbe-neo` through pkg-config.
- Public API: declarations under `include/mbelib-neo/`.
- Codec implementation: project-owned C sources under `src/`.
- Vendored FFT implementation: PFFFT/FFTPACK sources under
  `src/external/pffft/`.
- Build and packaging system: CMake, CMake presets, install rules, CMake package
  export, and pkg-config metadata.
- Test and validation system: CTest executables under `tests/`, fuzz targets
  under `fuzz/`, and CI workflows under `.github/workflows/`.

## Source Layout

- `include/mbelib-neo/`: installed public headers.
- `src/core/`: shared synthesis, parameter, adaptive, and unvoiced FFT logic.
- `src/ecc/`: error-correction helpers and constants.
- `src/ambe/`: AMBE frame and data processing.
- `src/imbe/`: IMBE frame and data processing.
- `src/internal/`: private headers shared by implementation units.
- `src/external/pffft/`: vendored third-party FFT sources.
- `tests/`: automated unit and regression tests.
- `examples/`: small consumer examples.
- `bench/`: optional local micro-benchmarks.
- `tools/`: local and CI-aligned quality scripts.

## Public Interface

The supported external software interface is the C API declared in
`include/mbelib-neo/mbelib.h`, plus the generated version header
`include/mbelib-neo/version.h` after configure/install.

Consumers normally use:

- `#include <mbelib-neo/mbelib.h>`
- CMake target `mbe_neo::mbe_shared` or `mbe_neo::mbe_static`
- pkg-config package `libmbe-neo`

The library writes 8 kHz audio frames and keeps stream state in caller-provided
`mbe_parms` objects. Processing APIs are reentrant when each stream has its own
state.

## Data Flow

1. The caller provides hard-decision frames, soft-decision frames, or unpacked
   parameter bits.
2. Codec-specific decode logic validates frame layout, applies error correction
   where applicable, and produces parameter bits or process status.
3. Synthesis logic updates caller-owned stream state and emits float or int16
   PCM samples.
4. Optional status helpers format process results for diagnostics.

## Trust Boundaries

Frame and data buffers supplied by callers are untrusted input. Public functions
must validate assumptions about dimensions, pointers, and bounded writes before
using those buffers.

Build scripts, CI workflow inputs, release credentials, and publication keys are
separate trust boundaries. CI workflows use least-privilege permissions and
dedicated analysis jobs to reduce the chance that untrusted changes can reach
release credentials.

## Build Architecture

CMake builds shared and static libraries from explicit source lists. The project
does not use source globbing for library sources, which keeps build membership
reviewable.

Install rules export:

- library artifacts
- public headers
- generated version header
- CMake package metadata
- pkg-config metadata
- license files

## Dependency Architecture

The main compiled third-party code is vendored PFFFT/FFTPACK. CI and analysis
tools are managed through pinned workflow actions, hashed Python requirement
files, and documented supply-chain guardrails.
