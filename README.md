<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# mbelib-neo

Performance‑enhanced IMBE/AMBE vocoder primitives and a TETRA ACELP speech decoder, with a modern CMake build, installable headers, and pkg-config/CMake package integration.

Project homepage: https://github.com/arancormonk/mbelib-neo

[![CI](https://github.com/arancormonk/mbelib-neo/actions/workflows/ci.yml/badge.svg)](https://github.com/arancormonk/mbelib-neo/actions/workflows/ci.yml)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/12975/badge)](https://www.bestpractices.dev/projects/12975)
[![OpenSSF Baseline](https://www.bestpractices.dev/projects/12975/baseline)](https://www.bestpractices.dev/projects/12975)

## Patent Notice

Use and distribution of the vocoders may implicate patents. The project has
not established third-party patent clearance; educational or research use is
not asserted as an exemption. AMBE+2 and TETRA remain enabled with the open
questions recorded in [the dated patent review](legal/patents.md).

## License

The project's contributions are GPL-2.0-or-later; original ISC mbelib notices
and pffft-neo's GPL and inherited PFFFT notices are retained. The source-only OP25 quality patch
is GPL-3.0-or-later. TETRA codebook data retains its ETSI attribution, while
redistribution permission and GPL compatibility remain unresolved. No blanket
license conclusion is asserted for the combined distribution.

Read the [component inventory](legal/component-table.md) and
[licensing status and retained-component exception](legal/status.md). Starting
with 2.6.0, releases include a **complete corresponding-source archive** with
the pinned FFT source, build inputs, tests and notices alongside the binaries.
See [offline source builds](docs/source-build.md) and
[release verification](docs/release-verification.md).

## Overview

- A performance‑enhanced fork of [lwvmobile/mbelib](https://github.com/lwvmobile/mbelib), which is a fork of [szechyjs/mbelib](https://github.com/szechyjs/mbelib)
- Supports IMBE 7200x4400 (P25 Phase 1), IMBE 7100x4400 (ProVoice), AMBE (D‑STAR), and AMBE+2 (DMR, NXDN, P25 Phase 2, dPMR, etc.).
- Decodes TETRA full-rate speech (ETSI EN 300 395-2 ACELP, 137 bits per 30 ms frame) with bad-frame concealment. The decoder is integer fixed-point and passes the EN 300 395-4 decoder test bit for bit (see [TETRA decoding](#tetra-decoding)).
- Public API in `#include <mbelib-neo/mbelib.h>` with version macro `MBELIB_VERSION`.
- v2.0.0 is an ABI-breaking release that adds strict bit-array validation, soft-decision ECC/frame decode, staged frame-decode helpers, and result-structured synthesis APIs; use older tags for the v1 API.
- Ships as shared and static libraries: `libmbe-neo.{so|dylib}` on Unix-like platforms, and `mbe-neo.dll` + import lib / `mbe-neo-static.lib` on Windows.
- Installable CMake package (`mbe_neo::mbe_shared` / `mbe_neo::mbe_static`) and pkg-config file (`libmbe-neo`).

## Build From Source

Requirements

- C compiler with C99 support.
- CMake ≥ 3.20.
- On non‑Windows platforms, links against `-lm` automatically.
- The first configure downloads a pinned [pffft-neo](https://github.com/arancormonk/pffft-neo) revision and verifies its SHA-256. No submodule or C++ compiler is needed.

Release source bundles automatically use their hash-verified FFT archive.
For a fresh offline build from a Git checkout, supply an extracted copy of the pinned
pffft-neo revision (see [dependency updates](docs/dependencies.md)):

```sh
cmake --preset dev-debug -DFETCHCONTENT_SOURCE_DIR_MBELIB_PFFFT=/absolute/path/to/pffft-neo
```

An already populated build directory can be reconfigured offline with
`-DFETCHCONTENT_FULLY_DISCONNECTED=ON`. pffft-neo is compiled into both libraries;
consumers need no separate FFT library.

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

## Tagged Releases

GitHub Actions builds release package artifacts when you push a tag matching `vX.Y.Z`.
The tag must match the project version declared in `CMakeLists.txt` and verify as an annotated tag signed by one of the trusted release keys, or the release workflow fails before packaging.

On a matching tag, CI rebuilds `dev-release`, reruns the release preset tests, packages the installed library tree for Linux, macOS, and Windows, generates SPDX SBOMs, creates GitHub artifact attestations, and creates or overwrites the matching GitHub Release assets.
Release verification instructions are in `docs/release-verification.md`.

Typical release flow:

```
# Update project(mbelib-neo VERSION X.Y.Z ...) first, commit, and push.

git tag -s vX.Y.Z -m "mbelib-neo X.Y.Z"
git push origin vX.Y.Z
```

## Configuration Options

- `-DNOTONES=ON` — Disable AMBE/AMBE+2 tone synthesis (adds `-DDISABLE_AMBE_TONES`).
- `-DMBELIB_BUILD_TESTS=ON` — Build CTest executables (default ON).
- `-DMBELIB_BUILD_EXAMPLES=ON` — Build `examples/` (default ON).
- `-DMBELIB_ENABLE_WARNINGS=ON` — Enable common warnings (default ON).
- `-DMBELIB_WARNINGS_AS_ERRORS=ON` — Treat warnings as errors (default ON).
- `-DMBELIB_ENABLE_ASAN=ON` — Enable AddressSanitizer in the Debug configuration.
- `-DMBELIB_ENABLE_UBSAN=ON` — Enable UndefinedBehaviorSanitizer in the Debug configuration.
- `-DMBE_ENABLE_DEBUG_LOGS=ON` — Verbose debug logging in codec sources.
- `-DMBELIB_BUILD_DOCS=ON` — Add `docs` target (requires Doxygen).
- `-DMBELIB_ENABLE_FAST_MATH=ON` — Enable fast-math (`-ffast-math`/`/fp:fast`) on library targets.
  PCM conversion is always compiled with IEEE semantics so NaN/Inf handling is preserved.
- `-DMBELIB_ENABLE_LTO=ON` — Enable IPO/LTO for the Release and RelWithDebInfo configurations when supported (single- and multi-config generators).
- `-DMBELIB_ENABLE_HARDENING=ON` — Enable supported Release-like compiler/linker hardening (default ON).
- `-DMBELIB_ENABLE_SIMD=ON` — Enable SSE2/AVX2/NEON routines when the compiler target provides them, otherwise use scalar routines. Defaults to ON for Release, RelWithDebInfo, and MinSizeRel builds of 64-bit targets (and for multi-config generators without a build type), and OFF otherwise. On baseline x86-64 the width-dependent kernels also get an AVX2/FMA copy that is chosen at run time on CPUs with AVX2. The bundled FFT also falls back to scalar when its compiler target has no supported SIMD backend. On 32-bit ARM, pass `-DCMAKE_C_FLAGS=-mfpu=neon` to enable NEON; the option itself preserves the toolchain's FPU selection. Measure performance on your own core. On 32-bit x86, this option compiles the library for SSE2 and therefore requires an SSE2-capable CPU; leave it OFF for baseline i386 portability. See [SIMD target selection](docs/build-installation.md#simd-target-selection).
- `-DMBELIB_ENABLE_NATIVE=ON` — Tune the libraries' release-like configurations for the build host's CPU (`-march=native`/`-mcpu=native`, or a probed `/arch:` level for MSVC and clang-cl) when that is safe: not when cross-compiling or `CMAKE_C_COMPILER_TARGET` is set, for Android/iOS/Emscripten, for foreign or several `CMAKE_OSX_ARCHITECTURES`, or when the compiler command, `CMAKE_C_FLAGS` or MSVC's `CL`/`_CL_` already select a target ISA (default ON). Debug is never tuned, and the flags never reach the exported targets or pkg-config. **Pass `-DMBELIB_ENABLE_NATIVE=OFF` to build portable binaries for distribution.** See [host tuning](docs/build-installation.md#host-tuning-mbelib_enable_native).
- `-DMBELIB_TEST_REQUIRE_AVX2=ON` — Fail, rather than skip, `test_golden_pcm_avx2` when the build or CPU cannot run the AVX2 kernels (for CI on AVX2 runners).
- Note: the `dev-release` preset enables SIMD, fast-math, and LTO by default when supported, and keeps the native default, so its binaries are tuned for the build host.
- `-DMBELIB_BUILD_BENCHMARKS=ON` — Build optional local micro‑benchmarks (not run in CI): `bench_synth`, `bench_unvoiced`, `bench_convert`, and `bench_encode` (D-STAR, AMBE+2 and IMBE encoders).
- `-DMBELIB_BUILD_TOOLS=ON` — Build the opt-in `mbe_quality_eval` public-API decoder and spectral analyzer, and the `mbe_tetra_decode` TETRA conformance-file decoder. See the [quality evaluation workflow](docs/testing.md#speech-quality-evaluation) and [TETRA conformance](docs/testing.md#tetra-conformance).

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
  printf("mbelib version: %s\n", mbe_versionString());
  return 0;
}
```

You can also build and run the bundled example:

```
cmake --build build/dev-debug --target example_print_version
./build/dev-debug/example_print_version
```

### Frame/Data API Quick Reference

| Codec path | Hard frame input | Soft frame input | Data API input | Scratch/output bits |
| --- | --- | --- | --- | --- |
| AMBE 3600x2400 | `char ambe_fr[4][24]` | `mbe_soft_bit ambe_fr[4][24]` | `char ambe_d[49]` | `ambe_d[49]` |
| AMBE 3600x2450 | `char ambe_fr[4][24]` | `mbe_soft_bit ambe_fr[4][24]` | `char ambe_d[49]` | `ambe_d[49]` |
| IMBE 7200x4400 | `char imbe_fr[8][23]` | `mbe_soft_bit imbe_fr[8][23]` | `char imbe_d[88]` | `imbe_d[88]` |
| IMBE 7100x4400 | `char imbe_fr[7][24]` | `mbe_soft_bit imbe_fr[7][24]` | `char imbe_d[88]` | `imbe_d[88]` (converted to 7200 layout) |

Use hard-frame `mbe_process*Frame*` APIs when you want one-call decode+synthesis from interleaved codec frames. These APIs accept const frame arrays, run demod/ECC, populate `mbe_process_result`, and write the scratch/output parameter bits.

Use `mbe_decode*Frame()` / `mbe_decode*SoftFrame()` when you need staged decode. These helpers accept const frame input, fill hard parameter bits, optionally populate `mbe_process_result`, and return the corrected error total without synthesizing audio. Use `mbe_process*Data*` for later synthesis when you want to carry that result context forward.

Use `mbe_process*Data*` when you already have unpacked parameter bits.

IMBE 7100x4400 frame decoders convert their `imbe_d[88]` output to the 7200x4400/IMBE 4400 layout; synthesize converted data with the IMBE 4400 data APIs, passing along the decode's `mbe_process_result` so its C0/C4 context and `MBE_PROCESS_FLAG_PROVOICE` apply (ProVoice keeps TIA-102.BABA eq 46's nominal pitch, where P25 plays a fit to the pitch DVSI's P25 decoder plays on speech, and mutes with library comfort noise).

- `mbe_ambe2400EncoderAlloc()` creates a caller-owned encoder context and FFT plan; `mbe_ambe2400EncoderReset(enc)` restarts its analysis state, and `mbe_ambe2400EncoderFree(enc)` releases it.
- `mbe_encodeAmbe2400Parms(enc, samples, ambe_d, cur_mp, prev_mp)` encodes 160 float PCM samples into 49 AMBE 2400 parameter bits. It is bit-compatible with this library's `mbe_decodeAmbe2400Parms()`/`mbe_processAmbe3600x2400*()` path and follows the D-STAR AMBE bit layout (interleave, scrambler and Golay parity cross-checked against the MMDVM tables). Its analysis follows the method of TIA-102.BABA chapter 5, and the spectral reconstruction it targets matches DVSI's AMBE-3000 D-STAR test vectors; on-air interoperability has not been verified.
- `mbe_encodeAmbe2400ParmsShort(enc, samples, ambe_d, cur_mp, prev_mp)` encodes 160 signed 16-bit PCM samples into 49 AMBE 2400 parameter bits.
- `mbe_encodeAmbe3600x2400Frame()` adds FEC and interleaving to 49 parameter bits, producing `char ambe_fr[4][24]`.
- `mbe_encodeDStarDVData()` packs a frame into nine D-STAR DV data bytes in air order, LSB first, without the sync word.
- `mbe_decodeDStarDVData()` unpacks nine D-STAR DV data bytes into `char ambe_fr[4][24]`.
- `mbe_ambe2450EncoderAlloc()`, `mbe_ambe2450EncoderReset()` and `mbe_ambe2450EncoderFree()` manage an AMBE+2 3600x2450 (DMR, NXDN, YSF, P25 Phase 2) encoder context. `mbe_encodeAmbe2450Parms(enc, samples, ambe_d, cur_mp, prev_mp)` and `mbe_encodeAmbe2450ParmsShort()` encode 160 float or 16-bit PCM samples into 49 parameter bits, quantized as TIA-102.BABA-1 clause 4 describes, and return `MBE_AMBE2450_FRAME_VOICE`, or `MBE_AMBE2450_FRAME_TONE` for a DTMF, KNOX, call-progress or single tone sent as a §7.2 tone frame. A call-progress tone is sent once it has filled three consecutive 20 ms analysis spans (an 80 ms tone always does; DVSI's encoder starts one frame earlier) and for up to two frames after it ends. A DTMF, KNOX or single tone is sent for the frames whose analysis span it fills at least 45%; where one changes directly to another, the frame whose span holds both carries the newer tone if it is a DTMF or KNOX tone filling the span's newest half, else the older tone if it fills the oldest half, as DVSI's encoder sends them. `mbe_encodeAmbe3600x2450Frame()` adds the Golay codes and the modulation of C1, producing `char ambe_fr[4][24]`.
- `mbe_imbe4400EncoderAlloc()`, `mbe_imbe4400EncoderReset()` and `mbe_imbe4400EncoderFree()` manage an IMBE 7200x4400 (P25 Phase 1) encoder context. `mbe_encodeImbe4400Parms(enc, samples, imbe_d, cur_mp, prev_mp)` and `mbe_encodeImbe4400ParmsShort()` encode 160 float or 16-bit PCM samples into 88 parameter bits, quantized as TIA-102.BABA chapter 6 describes, with the alternating synchronization bit of §6.5 in `imbe_d[87]`. `mbe_encodeImbe7200x4400Frame()` adds the Golay and Hamming codes and the modulation of vectors 1–6, producing `char imbe_fr[8][23]`.
- All three encoders share one speech analysis, and the two AMBE encoders one quantizer. Measured on DVSI's validation vectors against DVSI's own encoding of the same input, both decoded by this library, the AMBE+2 encoder matches its level within 0.1 dB (median) and the D-STAR and IMBE encoders within 0.2 dB, all three sit closer to the input speech (log-spectral distance 4.6 against 4.9 dB for D-STAR, 4.7 against 4.8 dB for AMBE+2, 4.5 against 4.7 dB for IMBE) and they make the same voicing decision for 82–89% of each 1 kHz band. At pitches below about 90 Hz, where the voicing test cannot tell noise from a voice, a 500 Hz band holding little more than the background noise the encoders track (steady mains hum among it) is not voiced ([DVSI parity](docs/dvsi-parity.md#encoder-voicing-over-background-noise)). A voice frame with no voiced band is sent at one fixed pitch code, b0 123 (D-STAR, 56 harmonics), 119 (AMBE+2, 56) or 127 (IMBE, 37), with its amplitudes measured on that code's harmonic grid, rather than at whatever pitch the tracker found in the noise: DVSI's D-STAR and P25 decoders read the pitch code of such a frame, this library's decoders read the D-STAR and IMBE codes as letting the voicing of the frame before last through the frame's first 10 ms in the bands whose 500 Hz column, or the column 250 Hz higher, it covers, and DVSI's are inferred to read them so from other codes outside their unvoiced ranges (a ThumbDV plays D-STAR b0 120 to 123 with no voiced band as low-pitch voice, but whether b0 123 or 127 lets the voicing of the frame before last through the frame's first 10 ms has not been measured on hardware, and no P25 hardware measurement exists) ([DVSI parity](docs/dvsi-parity.md#encoder-pitch-code-of-unvoiced-frames)). The AMBE+2 and D-STAR encoders choose a frame's six spectral vectors together: each starts as the codebook row nearest its own target, as TIA-102.BABA-1 4.3.1–4.3.2 describe, and two passes then refine them against the error of the amplitudes the decoder will reconstruct, weighted toward the strong harmonics; the bits stay valid for any decoder ([DVSI parity](docs/dvsi-parity.md#joint-choice-of-the-ambe-spectral-vectors)). On DVSI's rate-33 and D-STAR tone vectors the AMBE+2 and D-STAR encoders send DVSI's tone index for 98% of DVSI's tone frames, and no tone frame for any speech vector. The frame encoders reproduce DVSI's P25 and rate-33 frames bit for bit from their parameter bits. On-air interoperability has not been verified.

### Encoder Workflow

- Encoder state: use one `mbe_ambe2400_encoder` context per stream, with any number of contexts per thread; concurrent use of the same context requires external synchronization. Initialize with `mbe_ambe2400EncoderAlloc()` plus `mbe_initMbeParms()`, advance prediction with `mbe_moveMbeParms(cur_mp, prev_mp)` between frames, and restart with `mbe_ambe2400EncoderReset()` plus `mbe_initMbeParms()`. Encoding never allocates and does not modify `prev_mp`. For a voice frame, `cur_mp` equals the `cur_mp` of a bare `mbe_decodeAmbe2400Parms()` chain on the emitted bits, and the `prev_mp` the `mbe_processAmbe2400*` path leaves (that path keeps the prediction model in `prev_mp` and enhances its `cur_mp` in place), in `w0`, `L`, `gamma`, the frame mark (`PHIl[0]`), and `Vl`, `log2Ml` and `Ml` over harmonics 1..L; `K` is not among them: the encoder writes the number of voiced harmonics there, and the decoder leaves `K` as the caller has it. The analysis delay is about 10 ms; feed one final zero frame to flush the tail (a call-progress tone that ends with the input is held for two frames, so encode two to send them all). As with DVSI's encoder, quiet and silent input is coded as low-level voice and the decoded level follows the input level; there is no AGC. DTMF, call-progress and single tones are sent as D-STAR tone frames (b0 126), decided and timed as the AMBE+2 encoder decides its tone frames, since DVSI's D-STAR and rate-33 encoders send tone frames on the same frames; KNOX tones, which D-STAR cannot carry, are coded as voice. A tone frame leaves the decoder's prediction history unchanged, so `cur_mp` is then a copy of `prev_mp`.
- The AMBE+2 and IMBE encoders use the same calling sequence and timing, with `mbe_ambe2450_encoder` or `mbe_imbe4400_encoder` contexts; `cur_mp` must not be `prev_mp`. For a voice frame their `cur_mp` is, in every field (`K` included), what `mbe_decodeAmbe2450Parms()` or `mbe_decodeImbe4400Parms()` leaves in it decoding the emitted bits given `prev_mp`. Every IMBE frame and every AMBE+2 frame other than a tone frame is a voice frame. A tone frame leaves the decoder's prediction history unchanged, so `cur_mp` is then a copy of `prev_mp` and the usual `mbe_moveMbeParms(cur_mp, prev_mp)` keeps the history; to hear what was sent, run the emitted bits through `mbe_processAmbe2450Data()` with its own decoder state.

### Stateful Decode Workflow

- Keep one `mbe_parms` state triplet per audio stream/thread: `cur_mp`, `prev_mp`, and `prev_mp_enhanced`. The three must be distinct objects.
- In AMBE 3600x2450, `prev_mp` is the last valid voice frame (the prediction history) and `prev_mp_enhanced` is the last synthesized frame (through a mute or a tone, the model the decoder kept, which was not heard). When you drive `mbe_decodeAmbe2450Parms()` yourself, copy `cur_mp` into `prev_mp` only when `mbe_classifyAmbe2450Frame()` returns `MBE_AMBE2450_FRAME_VOICE`. `mbe_decodeAmbe2450Parms()` returns `MBE_AMBE2450_FRAME_VOICE` for silence frames too, as in 2.1, since both decode a model. A voice frame of V/UV rows 18 to 31 comes with its row and pitch code in `PHIl[0]` (below, **AMBE 3600x2450 frame types**). Since 2.6.13 it returns a frame of V/UV row 17 with every band unvoiced and every amplitude at the unvoiced scale, as TIA-102.BABA-1 Annex B lists the row (2.6.7 to 2.6.12 returned its 500 Hz column 0 voiced); the column voiced in the first half of the frame is carried in the frame mark (`PHIl[0]`), which only the AMBE 3600x2450 `mbe_process*` functions read, so `mbe_synthesizeSpeechf()` plays the decoded model as decoded, all unvoiced, unless its adaptive smoothing, which the error counts in the state drive, voices some of its harmonics ([DVSI parity](docs/dvsi-parity.md#rate-33-rows-13-15-and-17)).
- Initialize once before decoding with `mbe_initMbeParms(&cur_mp, &prev_mp, &prev_mp_enhanced)`.
- Prefer `mbe_process*Frame*` APIs when you have raw hard-decision vocoder frames and do not need staged parameter handling.
- For soft-decision input, fill `mbe_soft_bit.bit` with the hard decision and `mbe_soft_bit.reliability` with confidence (`0` is erasure-like, `255` is highly reliable), then call `mbe_decode*SoftFrame()` or `mbe_process*SoftFrame*()`.
- Helper constructors are available for common inputs: `mbe_softBitsFromHard()` assigns a fixed reliability and validates hard bits, while `mbe_softBitsFromLlr()` maps positive signed LLRs to bit 1 with magnitude clamped to `0..255`.
- Processing APIs report frame state through `mbe_process_result`: `c0_errors`, `protected_errors`, IMBE-specific `c4_errors`, `total_errors`, and flags such as `MBE_PROCESS_FLAG_SOFT_INPUT`, `MBE_PROCESS_FLAG_C0_VALID`, `MBE_PROCESS_FLAG_C4_VALID`, `MBE_PROCESS_FLAG_TONE`, `MBE_PROCESS_FLAG_ERASURE`, `MBE_PROCESS_FLAG_REPEAT`, `MBE_PROCESS_FLAG_MUTE`, `MBE_PROCESS_FLAG_SILENCE` (AMBE 3600x2450 silence frames), and `MBE_PROCESS_FLAG_PROVOICE` (context set by the IMBE 7100x4400 decoders).
- Use `mbe_formatProcessResult()` when you need a compact status string from a result. It writes `'='` repeated `total_errors` times, then any `E`, `T`, `R`, and `M` flags in that order, truncated to the supplied buffer size. The silence flag is not rendered.
- Hard input bits must be exactly `0` or `1`; soft `mbe_soft_bit.bit` values must also be exactly `0` or `1`. Invalid pointers/counters return `MBE_STATUS_INVALID_ARGUMENT`; invalid bit values return `MBE_STATUS_INVALID_BITS`. Of the caller's `mbe_parms` state, two reads are guarded: a sample of a saved unvoiced overlap (`previousUw`) that is not a finite number plays as 0, and so does a sample of the `mbe_synthesizeSpeechf()` noise overlap outside the range its noise generator writes (+0 to 53124 by its bit pattern) (defects inherited from 2.6.12 and earlier, found in review). `mbe_synthesizeSpeechf()` does not validate the rest of a caller's state: a previous model with a NaN amplitude (`Ml`) or an infinite phase (`PSIl`) still gives NaN output, as it did before.
- For parameter-only processing, seed `mbe_process_result.total_errors` and any valid C0/C4 context before calling `mbe_process*Data*`.
- Use `mbe_versionString()` for the static NUL-terminated library version string.

### TETRA Decoding

```c
mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();      /* one context per stream */
uint8_t bits[MBE_TETRA_FRAME_BITS];                    /* B1..B137, one 0/1 byte each */
int16_t pcm[MBE_TETRA_FRAME_SAMPLES];                  /* 240 samples, 30 ms at 8 kHz */
/* ... per received speech frame ... */
int status = mbe_tetraDecodeFrame(dec, crc_failed ? NULL : bits, crc_failed, pcm);
/* ... at the end of the call ... */
mbe_tetraDecoderFree(dec);
```

- Input is the 137 type-1 bits of one speech frame in the order of EN 300 395-2 clause 4.2.2.7, tables 3 and 4 (B1 is the most significant bit of the first LSP index), as a TETRA channel decoder delivers them. Each byte must be 0 or 1.
- Pass a nonzero `bad_frame` when the class-2 CRC failed or the frame was stolen for signalling; the frame is then concealed as clause 4.2.3.2 describes (previous LSPs, the last subframe-4 pitch, energies lowered 1.5 dB per subframe, previous innovation indices) and `bits` is not read. Bad frames before the first good frame since allocation or `mbe_tetraDecoderReset()` give silence; as in ETSI's reference decoder, the output after them can start with a burst.
- Output is plain 16-bit PCM at the level of ETSI's reference decoder; there is no float variant and no mbelib float scale.
- Returns 0, `MBE_STATUS_INVALID_ARGUMENT` (NULL `dec` or `pcm`, NULL `bits` on a good frame) or `MBE_STATUS_INVALID_BITS`. A rejected frame leaves the context and `pcm` untouched, so it can be retried as a bad frame.
- The context holds all decoder state and is not shared: use one per stream, any number per thread, and lock externally for concurrent use of one context. Decoding never allocates. The arithmetic is integer only, so the output is identical on every platform and build configuration.
- Conformance: the decoder passes the EN 300 395-4 decoder test: TE1_46.CHD decodes to TE2_128.PCM byte for byte. Its output also matches ETSI's reference decoder on every other stream compared: TE3, the tested loss patterns (including bad frames before the first good frame) and about 2.1 million frames of synthetic stress streams. `docs/tetra-acelp.md` records how the details the standard's text leaves open were resolved, and the measurements.

### Audio Sample Scaling (Float vs int16)

This section covers the IMBE/AMBE decoders. The TETRA decoder writes 240-sample int16 frames with plain saturation and has no float variant ([TETRA decoding](#tetra-decoding)).

- The IMBE/AMBE decoders synthesize 8 kHz audio in 160-sample frames.
- The `short` entry points write 16-bit PCM with soft clipping (~95% full-scale).
- The `*f` entry points return mbelib’s historical float scale (not normalized `[-1, +1]`). Every `short` API scales by 7 and clips at 95% of int16 full scale (NaN becomes 0). The ProVoice `short` APIs, `mbe_synthesizeSpeech()` and `mbe_floattoshort()` then truncate toward zero, so `mbe_floattoshort()` of a ProVoice float output equals the ProVoice `short` output. The D-STAR, AMBE 3600x2450 and P25 `short` APIs round to nearest instead, halves away from zero, as DVSI's D-STAR, AMBE 3600x2450 and P25 decoders round: to get their 16-bit samples from the float APIs, round `7 * x` to nearest (e.g. `roundf()`) after the same clip. DVSI's output separates rounding to nearest from truncation, not one tie rule from another; halves away from zero is this library's choice. The D-STAR and AMBE 3600x2450 mute noise lies under one 16-bit step in the float output, so a float caller who truncates gets silence where a muted frame plays the noise alone ([DVSI parity](docs/dvsi-parity.md#16-bit-output)).
- To feed a normalized float pipeline, scale each float sample by `(7.0f / 32768.0f)` (range is approximately `[-0.95, +0.95]` after soft clipping).
- Decoded speech comes out at the level of the original speech, the same level DVSI's AMBE-3000 decoders produce on identical bits. Releases up to 2.2.x played speech about 17 dB louder and soft-clipped 5–7% of the samples of normally recorded speech; adjust any downstream gain that compensated for that.

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
int main(void) { printf("%s\n", mbe_versionString()); }
EOF
cc consumer.c $(pkg-config --cflags --libs libmbe-neo) -o consumer.exe
```

Notes

- If you install to a custom prefix, set `PKG_CONFIG_PATH="<prefix>/lib/pkgconfig:$PKG_CONFIG_PATH"`.
- For shared builds, ensure `<prefix>/bin` is on `PATH` (so the DLL is found) or copy the DLL next to your `.exe`.
- To prefer static linking via pkg-config, you can try `pkg-config --static --libs libmbe-neo` (requires static libs available) or use the CMake package and link `mbe_neo::mbe_static`.

## Audio Quality Improvements

mbelib-neo combines regenerated MBE voiced phase with TIA-102.BABA smoothing and FFT-based unvoiced synthesis. The codec process functions render speech as DVSI's decoders do, in 10 ms subframes (see **DVSI parity** below); the public `mbe_synthesizeSpeechf()` renders it as TIA-102.BABA describes:

- **FFT-based unvoiced synthesis** (Algorithms #117-126): Uses 256-point FFT with band-level scaling instead of the legacy oscillator bank approach. Provides cleaner, more natural unvoiced sounds with proper spectral shaping. The codec paths use the same transform on DVSI's subframes.

- **Weighted Overlap-Add (WOLA)** (Algorithm #126): In `mbe_synthesizeSpeechf()`, unvoiced and voiced synthesis share the symmetric spec-exact trapezoidal window, with complementary weights across each 160-sample frame transition.

- **Adaptive smoothing** (Algorithms #111-116): Error-rate-based parameter smoothing that gracefully handles corrupted frames. Includes local energy tracking, adaptive voicing thresholds, and amplitude scaling.

- **Regenerated voiced phase**: The odd log-magnitude kernel from [US 5,701,390, Eqs. 7–9](https://patents.google.com/patent/US5701390A/en) derives harmonic phase from a frame's log spectral magnitudes. The codec process functions apply it, as DVSI's decoders do and as the patent's text describes, to the frame's decoded amplitudes (`log2Ml`, before spectral enhancement, the high-band shelf and any scaling; [DVSI parity](docs/dvsi-parity.md#phase-regeneration)); `mbe_synthesizeSpeechf()` applies it to the enhanced, smoothed amplitudes it renders. Neither changes the public API or parameter layout. Fully voiced output does not depend on the unvoiced noise.

- **Codec-specific frame repeat/muting**: IMBE mutes on max repeats or the error-rate threshold, and IMBE prolonged repeat headroom resets to a default model state. D-STAR mutes on max repeats and, for hard-decision input, while its error level (0.99 of the last level plus the frame's corrected bits, parity included) is above 342, which matches the muting of DVSI's AMBE-3000 under sustained channel errors on hardware and on DVSI's vectors (see [DVSI parity](docs/dvsi-parity.md#d-star-muting-under-channel-errors)). In every codec, a repeated frame keeps its own error accounting and continues the unvoiced-noise sequence; only the model is repeated. The exception is a rate-33 lost frame while a tone sounds: it repeats the tone, and DVSI's noise generator stands still through it, as through a tone frame.

- **AMBE 3600x2450 frame types per TIA-102.BABA-1** (P25 Half-Rate Vocoder Addendum):
  - **Silence frames** (b0 124/125) use ω₀ = 2π/32 and all bands unvoiced (§4.1). As on DVSI's decoder they play on 15 bands, where §4.1 gives L = 14, with their amplitudes a further 13.8 dB below the voice-frame law, and the previous frame's model plays on through their first 10 ms ([DVSI parity](docs/dvsi-parity.md#rate-33-silence-frames)). They are synthesized but never used for prediction, so gain and log-magnitude history stay with the last voice frame (§4.3, §4.4.1 eq. 26, §4.4.3 eq. 43).
  - **Error counts** (§5.5): every bit of the received C0 and C1 codewords, parity bits included, that differs from the corrected codeword, plus C0's extended parity bit when it disagrees with the corrected codeword (4 marks a detected, uncorrectable C0). Soft decisions count each such bit by its reliability.
  - **Repeats** (§5.6): erasure (b0 120–123), or C0 errors ≥ 4. DVSI's AMBE-3000 does not apply eq. 57–58's second condition (C0 ≥ 2 with ≥ 6 total), and neither does this library. These criteria apply to every frame, before tone classification. A repeat replays the last frame played unchanged (no re-enhancement or adaptive smoothing) and leaves the history untouched; after a tone it replays the tone, and after a mute the model the decoder kept, rising from silence (below).
  - **Departures followed from DVSI's AMBE-3000**: the pitch law log2 f0 = −4.322605 − 0.0220427 b0 in place of the Annex A table (12 harmonics, not 11, at b0 17), and u3's bits 8–5 carrying b4(2), b4(1), b4(0), b3(0) where the v1.0.5 text of Table 8 prints b3(0) first; DVSI's decoder and encoder both use them ([DVSI parity](docs/dvsi-parity.md)).
  - **V/UV rows 13, 15 and 17**: each row's own voicing, Annex B's, is the decoded voicing at the frame's end (the adaptive smoothing that channel errors bring can voice harmonics of any frame afterwards, a row-17 frame's included, whose end is then voiced), and the first half of the frame voices row 13's 500 Hz columns 0 and 1, no band for row 15 (whatever the frame before) and column 0 for row 17, as on DVSI's decoder; TIA-102.BABA-1 holds a frame's voicing for the whole frame and defines no model at mid-frame, which the library's synthesis in subframes needs. A lost frame that repeats such a frame plays the repeated model's own voicing through. The encoder never sends these rows, so its bits do not change; the parameter decoder's row 17 is under **Stateful Decode Workflow** above ([DVSI parity](docs/dvsi-parity.md#rate-33-rows-13-15-and-17)).
  - **V/UV rows 18 to 31**: Annex B decodes them as no band voiced, and the library plays them on DVSI's band grid (`w0 = 2.81 / L`) in rate 33's form of the pulsed state, by a rule of its own fitted to DVSI's decodes of its own encodings (what DVSI's decoder computes is not known, and the fit leaves a residual): the row names the pulsed 500 Hz channels of the frame's own subframe and of its mid-frame subframe, whose bands play one pulse instead of noise, and the pitch code sets the pulse's time (rows 24 and 25, which pulse both subframes, code two times on a grid of 11); each pulsed model plays its own pulse, through its own bands only. A frame whose adaptive smoothing voices a harmonic plays as before; a lost frame that repeats such a frame plays the repeated frame's own subframe again in both halves (observed on six lost frames, repeats of rows 19, 26 and 27; for the other rows the rule is an extension). The parameter decoder carries the row and code in `PHIl[0]` in place of a frame mark (the float 1024 + 128 (b1 − 18) + b0; a repeat's form adds 2048), which only the AMBE 3600x2450 `mbe_process*` functions read, and only on the band-grid model of the word's code, its pitch and harmonic count (a word on any other model plays as no word; on that model a caller's word is the history of such a frame, whatever the model's voicing in a lost frame's repeat, which gives it its held form or clears it, though a model with a voiced band plays no pulse; a word of another code of the same harmonic count passes too, codes 40 and 41 sharing a grid model); `mbe_synthesizeSpeechf()` and the other codecs treat it as no mark. The encoder never sends these rows, so its bits do not change ([DVSI parity](docs/dvsi-parity.md#rate-33-the-pulsed-rows)).
  - **Tones** (§7, §7.3): a frame whose first six bits of u0 equal 63 is a tone frame, whatever its b0 would read as (a single tone can read as silence). An invalid tone index is an erasure (flagged `ERASURE` and `TONE`), and so is a tone frame whose redundant tone fields disagree (a local consistency check; the spec leaves this to the decoder). Tone ID 255, a zero-amplitude tone in §7.3 and Table 9, is an erasure too, as on DVSI's decoder ([DVSI parity](docs/dvsi-parity.md#rate-33-lost-frames-and-tones)). Each tone component plays at the §7.2 level, 0.711 dB per AD step from +3.17 dBm0 at AD 127 (at the reference DVSI's decoder uses), and dual tones carry both components at that level. The frequencies are those DVSI's decoder plays, not Table 9's nominal ones: harmonics l1 and l2 of a fundamental on a grid of 2^-19 cycle per sample, Annex J's but for IDs 145 and 147, whose fundamentals are DVSI's own (DTMF "1" plays 693.909 + 1214.340 Hz, where Table 9's nominal 697 + 1209 Hz is what the library played before 2.6.13; DTMF tones lie within 0.48% of nominal, KNOX tones 1.55%, call progress 3.41%), and a single tone by its own rule, harmonic ID / 16 + 1 of a fundamental under 500 Hz, not Annex J's harmonic number (within 0.1 Hz of nominal). DVSI's decoder starts a tone as a fresh pulse train only after a frame with no voiced band, and continues a voiced frame's train into the tone; the library starts every tone that follows a frame that is not a tone frame fresh, at the phase DVSI's decoder starts one at after a frame with no voiced band ([DVSI parity](docs/dvsi-parity.md#tones-as-harmonics-of-a-fundamental)). A held tone in a caller's state counts only with a tone ID the decoder plays. The public `mbe_synthesizeTonef()` plays the same frequencies and continues the phases its state carries: its output differs from earlier releases', its signature does not.
  - **Muting** (§5.7): when the error rate exceeds 0.096, or instead of the 4th consecutive repeat. Where §5.7 outputs uniform noise in [−5, 5], the frame plays DVSI's mute noise (see **Mute noise** below), and the first muted frame after a heard one plays the last frame once more, fading out over its first 80 samples, under that noise; a tone left sounding fades out over the same samples. Underneath, the decoder's model goes on: each muted valid frame's own model becomes the one the next frame starts from ([DVSI parity](docs/dvsi-parity.md#rate-33-mutes)).
  - **Tone transitions**: as on DVSI's decoder, a tone joins the frames around it over 80 samples. It fades in over the first half of its first frame while the last synthesized frame, played once more, fades out (after a mute, from silence), and it goes on through the first half of the frame after its last, fading out on top of that frame. A tone that changes while it sounds changes at once. A lost frame while a tone sounds repeats the tone at its level, up to the repeat limit, whose mute fades it out. TIA-102.BABA-1 does not say how a tone frame joins the frames around it ([DVSI parity](docs/dvsi-parity.md#rate-33-tone-transitions), [lost frames and tones](docs/dvsi-parity.md#rate-33-lost-frames-and-tones)).
  - **Frames after a mute or a tone**: the first frame heard again is synthesized from the model the decoder kept through the mute or the tone (the last speech model before a tone), rising from silence over its first 10 ms with a new pulse train, as on DVSI's decoder after a mute; after a tone, the tone's tail fades out on top. After a tone DVSI's decoder continues the tone's pulse train into a voiced frame, where the library starts that frame from silence with a new train. A silence frame there rises to its own model. The spec doesn't define the synthesis state across these frames ([DVSI parity](docs/dvsi-parity.md#rate-33-frames-after-a-mute-or-a-tone)).

- **DVSI parity**: on DVSI's AMBE-3000 test vectors, D-STAR's gain is absolute (2 Dg[b2] − 0.09, not rate 33's differential law, which played onsets after silence up to 16 dB low) and D-STAR, P25 and rate 33 render DVSI's low band: a third-order Butterworth high-pass at 75 Hz in magnitude and phase, and an adaptive low cut whose corner follows the pitch of the stream's fully voiced frames (from 110 Hz, at most 200 Hz), under which voiced harmonics are scaled by (f / corner)^2 and the unvoiced bins removed; ProVoice and `mbe_synthesizeSpeechf()` keep the second-order filter at 100 Hz the codec paths played before. With it the waveform correlation of our decode of DVSI's bits with DVSI's decode rises by 0.0010 / 0.0016 / 0.0006 (D-STAR / P25 / rate 33, development vectors) and 0.0015 / 0.0018 / 0.0011 (validation), partition means, the log-spectral distance from DVSI's decode falls by 0.050 to 0.064 dB, and the distance from the input speech rises by 0.022 to 0.030 dB, the cost of removing fundamentals the input has as DVSI's decoders do ([DVSI parity](docs/dvsi-parity.md#dvsis-low-band), [measured](docs/dvsi-parity.md#november-2026-dvsis-low-band)). P25 plays a frame whose decoded magnitudes are low lower still, as DVSI's P25 decoder plays the quiet frames of its own vectors: every magnitude it renders times 2^clamp(1.84 (x + 1.73), -2, 0), x the mean of the frame's decoded log2 magnitudes, a fit to DVSI's output read through its rounding (not on a tone model, not in ProVoice; the prediction history and the encoders keep the decoded values), and P25's 16-bit entry points round to nearest as DVSI's do (ProVoice's keep truncating); with it P25's waveform correlation with DVSI's decode rises from 0.989085 to 0.990694 on the development vectors (the two near-silent recordings, `dam80` 0.968 to 0.994 and `dam40` 0.982 to 0.988) and does not fall on validation, and the log-spectral distance from DVSI's decode falls by 0.22 dB and that from the input by 0.07 dB on development (partition means; one declared check, no vector's distance from DVSI's decode rising by more than 0.002 dB, fails on `noisy`, +0.0024 dB, by the rounding; [DVSI parity](docs/dvsi-parity.md#p25-quiet-frames), [measured](docs/dvsi-parity.md#november-2026-p25s-quiet-frames)). P25 plays each pitch code at a fit to the pitch DVSI's decoder plays on speech, TIA-102.BABA eq 46's moved by at most 0.071% ([DVSI parity](docs/dvsi-parity.md#p25-pitch)), and a tone sent as P25 voice frames as DVSI's decoder plays it (below); ProVoice keeps the nominal pitch. Where a voiced model's pitch period about doubles or halves from one model to the next, the pitch-pulse train is first taken an octave down or up, as on DVSI's decoders ([DVSI parity](docs/dvsi-parity.md#the-pulse-train-across-an-octave)). After a frame with no voiced band, a P25 frame starts with the frame only the voiced harmonics in, or within 250 Hz below, the 500 Hz columns its voicing covers, and the others at mid-frame, a rule fitted to DVSI's P25 decodes (it agrees with 1,838 of the 1,843 onsets labelled decisively) ([DVSI parity](docs/dvsi-parity.md#p25-the-bands-that-start-with-the-frame)). Rate 33's V/UV rows 13, 15 and 17 voice the first half of the frame as DVSI's decoder does (columns 0 and 1, none, and column 0), each row's own voicing, Annex B's, holding at the frame's end as decoded ([DVSI parity](docs/dvsi-parity.md#rate-33-rows-13-15-and-17)). Each of the two marks acts only in its own codec's path; P25 reads the voicing it renders, after adaptive smoothing, and gives a voiced model it renders a mark of its own, which replaces any mark a caller's state put there. With the two, the waveform correlation of our decode of DVSI's bits with DVSI's decode rises from 0.9687 to 0.9849 (P25, development vectors) and 0.9855 to 0.9896 (validation), and from 0.9737 to 0.9806 and 0.9788 to 0.9834 (rate 33), partition means ([DVSI parity](docs/dvsi-parity.md#november-2026-mid-frame-voicing)). With rate 33's mid-frame pitch wave (below) its correlation rises further, from 0.9806 to 0.9893 (development), 0.9834 to 0.9921 (validation) and 0.9783 to 0.9860 (final, held out), partition means ([DVSI parity](docs/dvsi-parity.md#november-2026-rate-33s-mid-frame-pitch-wave)). In D-STAR and P25 the mid-frame model voices each band by one rule, fitted to per-harmonic labels of DVSI's decodes (which of two test decodes is nearer DVSI's in each harmonic's band), with which a DVstick 30's step curves come closer: a band in the same state in both frames keeps it, and one that changes is voiced by 500 Hz columns read in both frames; with it the correlation rises from 0.9907 to 0.9916 (D-STAR, development) and 0.9820 to 0.9935 (validation), and from 0.9849 to 0.9867 and 0.9896 to 0.9908 (P25), partition means ([DVSI parity](docs/dvsi-parity.md#d-star-and-p25-mid-frame-voicing-by-columns), [measured](docs/dvsi-parity.md#november-2026-mid-frame-voicing-by-columns)). In rate 33 the pitch-pulse train runs on from a frame of V/UV row 17 rendered with no voiced band into a frame voiced at mid-frame after it (a frame with a voiced band of its own, or another row-17 frame), where since the mid-frame voicing change a new train had started (where channel errors' adaptive smoothing voices the row-17 frame's end, it carried before too); with that, rate 33's correlation rises from 0.9893 to 0.9924 (development) and stays at 0.9921 (validation, where one decode changes in two frames), partition means ([DVSI parity](docs/dvsi-parity.md#november-2026-the-pulse-train-through-a-row-17-frame)). In D-STAR and P25 a frame with no voiced band at one of the 56 pitch codes DVSI's encoders keep for such frames plays DVSI's pulsed state, a single pulse in place of noise in the 500 Hz channels its code names (below; [DVSI parity](docs/dvsi-parity.md#d-star-and-p25-the-pulsed-state)). Rate 33's V/UV rows 18 to 31 play rate 33's form of it, each pulsed model its own pulse at a time the pitch code gives; with that, rate 33's whole-file correlation with DVSI's decode rises from 0.992438 to 0.993226 (development) and 0.992119 to 0.993591 (validation), every vector rising, at a cost of 0.008 and 0.016 dB of log-spectral distance from the input ([DVSI parity](docs/dvsi-parity.md#rate-33-the-pulsed-rows), [measured](docs/dvsi-parity.md#november-2026-rate-33s-pulsed-rows)). The encoders quantize with the same laws. What was measured, and what still differs, is in [`docs/dvsi-parity.md`](docs/dvsi-parity.md).

- **P25 tones sent as voice frames**: TIA-102.BABA has no tone frames, and DVSI's encoder sends a DTMF, KNOX, call-progress or single tone as voice frames at the pitch code of the tone's fundamental (Annex J's f0 for a dual tone, but DVSI's own for IDs 145 and 147; 31.25 ID / l for a single tone, l Annex J's harmonic number, where Annex J prints a multiplier for two ranges, `10.417 ID` and `5.2803 ID`) with a voicing word that is the same in every steady frame of a tone (12,555 steady frames compared; 112 frames at tone boundaries carry another pair). A P25 voice frame whose pitch code and voicing word are a DTMF, KNOX or single tone's, and whose tone harmonics dominate the others (two limits measured on DVSI's decodes), plays as DVSI's decoder plays it: the tone's harmonics at one level at eq 46's fundamental cut to a 2^-19 grid, the other harmonics silent, the frame's own model at mid-frame; call-progress tones, which DVSI's decoder does not play as equal-amplitude tone models, and ProVoice frames play as decoded. The decoded frame stays the prediction history, and the tone it carries is left in its otherwise unused `Vl[0]` (0 after any other frame decoded; a status error writes nothing; a word in a caller's state counts only if the decode's own test, run again on the model it sits in, gives that word) ([DVSI parity](docs/dvsi-parity.md#p25-tones-sent-as-voice-frames)).

- **Synthesis in 10 ms subframes** (every `mbe_process*` function): DVSI's decoders synthesize a frame as two subframes, and so do the codec paths here ([DVSI parity](docs/dvsi-parity.md#synthesis-in-subframes)).
  - A frame's own model stands at its end, and a model between the previous frame's and this one's at mid-frame (pitch the geometric mean, magnitudes met in the log domain; voicing, in D-STAR and P25, the state both frames agree on and otherwise set by 500 Hz columns read in both frames, and in rate 33 this frame's but where the frame marks below say otherwise). Between two voiced rate-33 frames of different pitch the mid-frame pitch is the geometric mean times 1 − 3.0e-4 times 1 + (h(Pa) + h(Pb)) / 2, where h is a wave in each frame's pitch period P that repeats every half octave, the factor and the wave fitted on DVSI's decodes of the development vectors; what in DVSI's decoder makes the wave is not known ([DVSI parity](docs/dvsi-parity.md#rate-33-mid-frame-pitch)). Each subframe crossfades linearly from one model to the next.
  - Voiced bands play as a train of pitch pulses, each sounding its model's harmonics for one pitch period and, fading, 16 samples into the next pulse's; a subframe with no voiced band at either end starts the train afresh. Voiced harmonics render 0.55 dB louder.
  - **DVSI's low band** (D-STAR, P25 and rate 33, not ProVoice): every voiced and pulsed harmonic and every unvoiced bin takes the gain and phase of a third-order Butterworth high-pass at 75 Hz. A tracker, 110 Hz at a stream's start, moves toward the fundamental of each frame rendered with every band voiced (`T + 0.043 (f0 - T)`; not a P25 tone model or a silence frame) and stands still through every other frame; the corner a frame applies is the smaller of it and 200 Hz, a voiced harmonic under it is scaled by (f / corner)^2 and the unvoiced bins under it are zero. The tracker and the corner travel in the stream's state with the pulse train; a value no frame stores there reads as a fresh stream's. The decoded models keep their amplitudes: the cut belongs to the rendering ([DVSI parity](docs/dvsi-parity.md#dvsis-low-band)).
  - Unvoiced bands play noise from DVSI's 16-bit generator (`x <- 173 x + 13849 mod 65536` from 0, the same sequence in every stream) under a square-root triangular window two subframes long. `mbe_setThreadRngSeed()` does not reach it.
  - **The pulsed state** (D-STAR and P25, not ProVoice or rate 33): a frame with no voiced band at one of the 56 pitch codes DVSI's encoders keep for such frames (P25 b0 128..183, D-STAR b0 64..119) plays the pattern its code names, as DVSI's decoders play it on the streams compared (DVSI's own encodings, on which P25's counter never falls below 0): for the frame's own subframe and its mid-frame subframe, the unvoiced bands in the pulsed 500 Hz channels (the patterns are Table 2 of US 7,970,606; D-STAR's codes carry the same ones) play no noise but the response to one pulse, 24 samples after the model's point, a model emitting only when the one before it did not. The code comes from the model as rendered, after adaptive smoothing (its pitch compared by bit pattern with the decoder's own). D-STAR reads it whenever the model qualifies; P25 only while a validity counter, as US 7,970,606 describes DVSI's P25 decoder, allows it (frames with no voiced band at another P25 code's pitch lower it, those at the table's codes raise it, and a model at no code's pitch, such as the one the repeat limit installs, leaves it), a guard against a stream from an encoder that does not use the codes, taken from the document (the counter, and that a code it does not read takes no part in the mid-frame pitch, are this library's reading of the published description, not measurements): on DVSI's own encodings it never stops a code. Lost frames leave the counter where the model they replay and the repeat limit's reset model leave it, so after a run of them the first code is read only if the counter stood at 0 or more before it. This library's encoders send none of the codes. With it the waveform correlation of our decode of DVSI's bits with DVSI's decode rises from 0.9916 to 0.9926 (D-STAR, development vectors), 0.9935 to 0.9947 (validation) and 0.9856 to 0.9893 (final, held out and scored once, after the reviews), and from 0.9867 to 0.9875, 0.9908 to 0.9920 and 0.9869 to 0.9908 (P25), partition means, every vector rising; the log-spectral distance from DVSI's decode falls by 0.16 to 0.19 dB, and the distance from the input speech rises by 0.005 to 0.026 dB, the cost of matching DVSI here (on the library's decodes PESQ against the input falls by 0.02, to within 0.01 of the score of DVSI's own decode; no equivalence was tested) ([DVSI parity](docs/dvsi-parity.md#d-star-and-p25-the-pulsed-state), [measured](docs/dvsi-parity.md#november-2026-the-pulsed-state)).
  - Where a voiced run starts or ends inside a frame follows what DVSI's decoders read from fields the specifications leave free: the pitch code of a P25 or D-STAR frame with no voiced band ([DVSI parity](docs/dvsi-parity.md#voiced-onsets-and-offsets)) and rate 33's V/UV rows 13, 15 and 17 ([DVSI parity](docs/dvsi-parity.md#rate-33-rows-13-15-and-17)). In D-STAR and P25 the mid-frame voicing also follows a rule of 500 Hz columns read in the two frames, which this library fitted to labels of DVSI's decodes; it is not something DVSI's decoders are known to read ([DVSI parity](docs/dvsi-parity.md#d-star-and-p25-mid-frame-voicing-by-columns)).
  - On DVSI's speech vectors our decode of DVSI's bits then follows DVSI's waveform sample for sample: correlation 0.98 to 0.99 (the correlation of the two waveforms sample against sample at lag 0) per codec and partition (means over each partition's vectors, with the mid-frame columns change and the row-17 train change; D-STAR / P25 / rate 33 development 0.9916 / 0.9867 / 0.9924, validation 0.9935 / 0.9908 / 0.9921, final 0.9856 / 0.9869 / 0.9860, the final partition's decodes unchanged by the row-17 train change, in every mode; with the pulsed state D-STAR / P25 development 0.9926 / 0.9875, validation 0.9947 / 0.9920, final 0.9893 / 0.9908, rate 33 unchanged; with rate 33's pulsed rows rate 33 development 0.9932, validation 0.9936, final 0.9896, from 0.9924 / 0.9921 / 0.9860, D-STAR and P25 unchanged; with DVSI's low band D-STAR / P25 / rate 33 development 0.9936 / 0.9891 / 0.9938, validation 0.9961 / 0.9938 / 0.9947, final 0.9905 / 0.9919 / 0.9904, from 0.9893 / 0.9908 / 0.9896; with P25's quiet frames P25 development 0.9907 and validation 0.9938, the final partition not yet scored) on the development, validation and final partitions, the final one, `irstia` and `mark`, held out and scored once for each of the mid-frame voicing change, the wave, the columns change, the row-17 train change, the pulsed state, rate 33's pulsed rows and DVSI's low band, each after its review (for the row-17 train change it shows no harm and cannot show the gain, nothing in it changing; [DVSI parity](docs/dvsi-parity.md#november-2026-the-pulse-train-through-a-row-17-frame); for the pulsed state the gain declared in advance is there on both vectors in both codecs, [DVSI parity](docs/dvsi-parity.md#the-final-partition); for rate 33's pulsed rows on both vectors, [DVSI parity](docs/dvsi-parity.md#november-2026-rate-33s-pulsed-rows); for DVSI's low band on both vectors in all three codecs, the least rate 33's on `irstia`, +0.000046, [DVSI parity](docs/dvsi-parity.md#november-2026-dvsis-low-band)) ([DVSI parity](docs/dvsi-parity.md#november-2026-mid-frame-voicing-by-columns)), and 0.99 or better on the median frame with no voiced band ([DVSI parity](docs/dvsi-parity.md#synthesis-in-subframes)).
  - `mbe_synthesizeSpeechf()` keeps TIA-102.BABA eq 117-141, with voiced and unvoiced harmonics 0.3 dB louder and 0.4 dB quieter and its own noise generator.

- **D-STAR tones**: single tones 5–122 as AMBE+2, DTMF 128–143 (128 + 4 × column + row) and call-progress tones 144–147, at the level the 8-bit tone volume gives on DVSI's decoder, volume 0 included. As on DVSI's decoder, a tone crossfades with the speech over the first 80 samples of a frame, and a tone that goes on at another volume moves to it over those samples. In the first tone frame, the last frame repeated fades out under the tone, if that frame was heard. In the first voice frame after it, the tone goes on and fades out under the voice ([DVSI parity](docs/dvsi-parity.md#d-star-tone-transitions)). Each tone plays the frequencies DVSI's decoder plays for the AMBE+2 tone of the same nominal pair (above), starting as AMBE+2's do, and `mbe_synthesizeTonefdstar()` plays them too, continuing the phases its state carries. DVSI's decoder starts a tone as a fresh pulse train only after a frame with no voiced band and continues a voiced frame's train into the tone; the library starts every tone after a frame that is not a tone frame fresh. After a tone DVSI's decoder continues the tone's train into a voiced frame; in the library, when the frame before the tone was heard (muted neither by the error level nor by the repeat limit), that frame is played once more under the tone's first frame, and a D-STAR voice frame after the tone goes on from the speech pulse train that repeat left; rate 33 starts that frame from silence with a new train.

- **D-STAR lost frames**: DVSI's AMBE-3000 plays b0 125, a b0 126 tone frame whose index carries no tone, and a b0 127 frame (among them the zero-volume index-128 frame that earlier versions of this library's encoder sent for silence) exactly as it plays a lost frame, and so does this library: it repeats the previous frame (flagged `ERASURE` and `REPEAT`, plus `TONE` for b0 126/127), and the 4th in a row mutes. The all-ones DV frame (`FF` x 9) is one of those b0 127 frames: DVSI's host software writes it into bit files for the frames its DTX simulation does not send, and only that software's decoding side plays it as comfort noise, by flagging it to the chip. A frame whose C0 word has four errors, which the extended Golay code detects, plays as a lost frame too (flagged `C0_UNCORRECTABLE` by the frame decoders), while the next frame still predicts from its decoded amplitudes; corrected errors repeat no frame, and a tone frame with corrections plays, as on DVSI's AMBE-3000 ([DVSI parity](docs/dvsi-parity.md#d-star-channel-errors)). A played tone mutes with the error level as voice does, and lost frames right after a tone follow DVSI's AMBE-3000: after a single tone the tone's tail fades out and silence follows, while a DTMF or call-progress tone goes on until the repeat limit's mute. The voice after that mute starts from the prediction history about 10 dB down, as on DVSI's AMBE-3000, rather than fading in from silence. b0 120–123 are voice at the pitch law. b1's four bits voice the 1 kHz bands from 0 Hz up, each voiced band reaching 250 Hz into the band below, as on DVSI's AMBE-3000 ([DVSI parity](docs/dvsi-parity.md#d-star-voicing-bands)). b0 124, the null frame among them, plays with every band unvoiced (b1 ignored) as on DVSI's decoder, its amplitudes decoded and predicted as voice; the library uses the AMBE+2 silence model (ω₀ = 2π/32, L = 14), the best fit to DVSI's output of the models measured, which the hardware supports over the 9-harmonic alternative ([DVSI parity](docs/dvsi-parity.md#d-star-b0-124)). The encoder sends b0 0–123 for voice, as DVSI's encoder does. Measurements: [DVSI parity](docs/dvsi-parity.md#d-star-lost-frames).

- **LCG noise generator with buffer overlap**: in `mbe_synthesizeSpeechf()`, the TIA-102.BABA §11.2 linear congruential generator for deterministic noise, with 96-sample overlap for smooth continuity between frames. The codec paths draw from DVSI's generator instead and keep their pitch-pulse state in that buffer.

- **Mute noise**: P25 IMBE mutes with the spec level, uniform in [−5, 5] at the int16 output (±5/7 in the float output, which is int16 / 7; TIA-102.BABA §7.8). D-STAR and AMBE 3600x2450 mute at the level of DVSI's muted output in both modes, not TIA-102.BABA-1 §5.7's [−5, 5]: uniform in [−0.75, 0.75] on the int16 scale, which their `short` APIs round to 0 on two thirds of the samples and ±1 on the rest (a caller of the float APIs who truncates gets silence). D-STAR's repeat limit (the 4th lost or errored frame in a row) fades the last frame out under that noise and keeps the prediction history, as DVSI's decoder does, and an AMBE 3600x2450 mute fades the last frame out under it the same way. ProVoice (not a TIA-102 codec; the 7100x4400 frame APIs, or the IMBE 4400 data API given a result carrying `MBE_PROCESS_FLAG_PROVOICE`) and direct `mbe_synthesizeSpeechf()` callers use the library comfort-noise level (`0.003` gain semantics, also used by `mbe_synthesizeComfortNoisef()`). All of them draw from a Java `Random`-compatible per-thread RNG.

- **Nyquist guard**: the shared synthesizer skips voiced harmonics at or above Nyquist (`l * w0 >= pi`) instead of rendering them as aliases. Decoded models never contain such harmonics; the guard protects caller-supplied models.

The opt-in [encode/decode quality pipeline](docs/testing.md#speech-quality-evaluation) compares a frozen baseline and candidate using identical encoded frames. It reports spectral, crest-factor, envelope, and reference-matched frame-join measurements; these are diagnostics, not a substitute for listening or a standardized perceptual score.


## API Notes

- Public API is prefixed `mbe_` and declared in `mbelib.h`.
- Float and 16‑bit PCM variants are provided (e.g., `mbe_processAmbe3600x2400Framef` and `mbe_processAmbe3600x2400Frame`).
- Version macro `MBELIB_VERSION` is defined in the generated header `mbelib-neo/version.h`.
- `mbe_versionString()` returns a const pointer to the version string.
- Processing APIs return an `int` error total and use `mbe_process_result` for decode/synthesis status.
- Soft-decision ECC helpers are public for Golay(23,12), Hamming(15,11), and the IMBE 7100x4400 Hamming mapping.
- Unvoiced synthesis keeps its generator state in `mbe_parms` (`noiseSeed`): DVSI's 16-bit generator on the codec paths, which every stream starts from the same state, and the TIA-102.BABA §11.2 LCG in `mbe_synthesizeSpeechf()`, with per-thread cold-start seeding via `mbe_setThreadRngSeed(uint32_t)`. Treat `mbe_parms` as opaque stream state between calls: the otherwise unused element 0 of `Vl`, `PHIl` and `PSIl`, and `noiseOverlap`, carry the codec paths' own state (in `noiseOverlap` the pulse train, the pulsed state and DVSI's low cut's tracker). A stream moved between a codec path and `mbe_synthesizeSpeechf()` keeps its model and noise position; the unvoiced overlap of the rendering it leaves is not played by the other, so it fades in.
- ABI: v2.0.0 breaks ABI from v1. The shared library `SOVERSION` follows the project major version; minor updates within a major version aim to remain ABI compatible.

### Determinism & RNG

- `mbe_setThreadRngSeed(seed)` seeds both thread-local comfort-noise RNG state and the next unvoiced LCG cold start of `mbe_synthesizeSpeechf()`. The comfort-noise RNG also draws every mute noise of the `mbe_process*` functions (P25's TIA-102.BABA §7.8 noise, DVSI's AMBE mute level of D-STAR and rate 33, ProVoice's comfort noise), so a muted frame's samples depend on that thread's RNG state, not on the stream's state alone. For reproducible output of either, set it per thread before synthesis. The unvoiced noise of the `mbe_process*` functions takes no seed: every stream draws the same sequence from a fresh state, as DVSI's decoders do.
- `mbe_setThreadRngSeed(0)` is accepted and remapped internally to a non-zero seed. Use an explicit non-zero seed when exact reproducibility matters across builds.
- Unvoiced noise progression after cold start is driven by the per-frame generator state in `mbe_parms`, so deterministic playback requires carrying frame state forward consistently, and, for streams with muted frames, the same thread RNG state (`mbe_setThreadRngSeed()`) before them, since mute noise is drawn from it.
- Runtime synthesis helpers (RNG state and FFT plan) are thread-local; do not share mutable decode state (`mbe_parms`) across threads unless you synchronize externally.
- Frame handling differs per codec: IMBE uses error-rate muting and a repeat-headroom reset (with the TIA-102.BABA mute-noise level for P25, library comfort noise for ProVoice), D-STAR uses repeat-driven muting plus an error-level mute for hard-decision input, and AMBE 3600x2450 follows TIA-102.BABA-1 (see Audio Quality Improvements).
- Enabling `MBELIB_ENABLE_SIMD=ON` selects vectorized math when the compiler target supports it, with
  scalar fallback otherwise. This can change
  floating‑point rounding at the bit level. Tests enforce exactness for int16 on x86 in Debug, and
  sanity bounds elsewhere.

## Benchmarks (Optional)

Build micro-benchmarks:

```
cmake --preset dev-release -DMBELIB_BUILD_BENCHMARKS=ON
cmake --build --preset dev-release -j --target bench_synth bench_unvoiced bench_convert bench_encode
```

Run benchmark executables:

```
./build/dev-release/bench_synth
./build/dev-release/bench_unvoiced
./build/dev-release/bench_convert
./build/dev-release/bench_encode
```

Scalar-vs-SIMD comparison for all three encoders, voiced/unvoiced synthesis and PCM conversion (builds ad hoc comparison trees under `build/bench-compare-scalar/` and `build/bench-compare-simd/`):

```
tools/bench_compare.sh
```

## Tests and Examples

- Run tests with `ctest --preset dev-debug -V` (or `ctest -V` from the build directory).
- Included tests: `test_api` (version/header/result helpers), `test_ecc` (hard and soft Golay/Hamming), `test_noise_determinism` (unvoiced RNG/frame-state determinism), `test_unvoiced_subframes` and `test_subframe_synthesis` (the codec paths' synthesis in 10 ms subframes: DVSI's noise generator, the mid-frame model, the pitch-pulse train, the frame marks the decoders set and the 500 Hz columns a P25 onset starts with), `test_params` (parameter/synthesis behavior, soft frame decode, and v2 wrappers), `test_floattoshort_parity` (exact float-to-int16 conversion parity), `test_pcm_rounding` (the AMBE paths' round-to-nearest int16 conversion in every compiled variant against an independent reference), `test_ambe2450_frame_types` (AMBE 3600x2450 silence, repeat, erasure, tone and mute handling per TIA-102.BABA-1 and DVSI's decoder: silence frames, mutes, the frame after a mute or a tone, lost frames around tones, the mid-frame voicing of V/UV rows 13, 15 and 17), `test_ambe2400_frame_types` (D-STAR lost-frame, null-frame, tone and all-ones frame classes), `test_tones` (AMBE+2 and D-STAR tone levels and the frequencies DVSI's decoders play, against the DVSI-calibrated laws), `test_p25_tones` (P25 tones sent as voice frames: the decoder's tone test, the tone model and its rendering), `test_golden_pcm` (golden hash regression checks, including voice-only AMBE 2450 and D-STAR sequences), `test_tetra_decoder` (TETRA decoding stages, concealment, API contract and stream hashes), `test_fixed_point` (the TETRA decoder's saturating operators), `test_tetra_shared_api` (TETRA entry points through the shared library), `test_tetra_tables` (table generators and provenance), `test_tetra_conformance_tools` (the TETRA comparison tool).
- Example: `examples/print_version.c` shows linking and header usage.
- WAV round trip: `./build/dev-debug/dstar_encode < input.wav > output.dstar`, then `./build/dev-debug/dstar_decode < output.dstar > output.wav`. These example binaries are not installed. They use binary stdin/stdout, accept no path arguments, and send diagnostics to stderr. The private `.dstar` container is not a D-STAR air-interface stream.
- When Python 3.7 or newer is available, CTest also checks the D-STAR examples through stdin/stdout pipes, including malformed WAVs, exact silence, and the final flush frame.
- Golden hash helper: `gen_golden` (available when `MBELIB_BUILD_TESTS=ON`, default) prints current FNV-1a reference values for synthesis/conversion regression workflows (`cmake --build build/dev-debug --target gen_golden && ./build/dev-debug/gen_golden`).

## Documentation

Optional Doxygen documentation can be generated:

```
cmake --preset dev-debug -DMBELIB_BUILD_DOCS=ON
cmake --build --preset dev-debug --target docs
# Output in build/docs/html
```

The generated site focuses on the public API (`include/`) and bundled examples.

Project documentation and policies:

- Contribution process, DCO, coding standards, review requirements, and test policy: `CONTRIBUTING.md`
- Defect reporting: `docs/issue-reporting.md`
- Architecture and trust boundaries: `docs/architecture.md`
- Build and installation policy: `docs/build-installation.md`
- Security requirements and design: `docs/security-requirements.md`
- Dependency management and monitoring: `docs/dependencies.md`
- OpenSSF OSPS Baseline evidence: `docs/openssf-baseline.md`
- Testing, CI, coverage target, and dynamic analysis: `docs/testing.md`
- Parity with DVSI's AMBE-3000 decoders and encoders (measurements, departures from the specification, known differences): `docs/dvsi-parity.md`
- Release verification and signing: `docs/release-verification.md`
- Supply-chain guardrails: `docs/supply-chain-guardrails.md`
- Code-quality guardrails: `docs/code-quality-guardrails.md`

## Project Layout

- Public headers: `include/mbelib-neo/`
- Sources: `src/core/`, `src/ecc/`, `src/ambe/`, `src/imbe/`, `src/tetra/`
- Internal headers: `src/internal/`
- External libraries: [pffft-neo](https://github.com/arancormonk/pffft-neo), pinned by `cmake/Pffft.cmake` (`GPL-2.0-or-later AND LicenseRef-PFFFT`)
- Tests: `tests/` • Examples: `examples/` • Benchmarks: `bench/`
- Developer tooling and CI parity scripts: `tools/`

## Contributing

- Read `CONTRIBUTING.md` before opening a pull request.
- Follow `.clang-format` (LLVM style, 4‑space indent, 120 cols). You can run `tools/format.sh`.
- Static analysis scripts (CI-aligned; CI and pre-push/preflight use strict mode where applicable):
  - `tools/clang_tidy.sh` (promotes broad bugprone/performance/portability findings; targeted TUs supported).
  - `tools/cppcheck.sh` (strict uses exhaustive value-flow, inconclusive findings, and targeted TUs).
  - `tools/iwyu.sh` (strict include hygiene via include-what-you-use, including project headers; excludes `src/external`).
  - `tools/gcc_fanalyzer.sh` (strict GCC `-fanalyzer` plus additional high-signal GCC warnings; excludes `src/external`).
  - `tools/scan_build.sh` (strict Clang Static Analyzer with extra security/portability checkers; excludes `src/external`).
  - `tools/semgrep.sh` (strict SAST with default, C, and security-audit rule packs).
  - `tools/lizard.sh` (strict complexity/length/argument budgets; excludes `src/external`).
- Git hooks: `tools/install-git-hooks.sh` enables auto-format on commit and CI-style pre-push checks.
- Manual preflight runner: `tools/preflight_ci.sh` runs the local CI quality gates without pushing.
- Prefer keeping internal symbols `static` and declarations in headers where shared.
- Before sending changes: build locally, run `ctest -V`, and ensure examples still link.
