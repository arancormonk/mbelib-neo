<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Testing Policy

mbelib-neo uses automated tests, static analysis, sanitizers, and fuzzing to
reduce regression and security risk.

## Test Suites

The CTest suite includes:

- API/version/result helper checks
- ECC tests for hard and soft Golay/Hamming paths
- AMBE 2400 encoder round trips, exact `log2Ml` parity, pitch endpoints, independent contexts and reset replay in one thread, plus behaviour on synthetic speech: noise of any colour and level stays unvoiced, steady vowels from 70 to 310 Hz are voiced and on pitch, harmonics below 2 kHz with noise above voice only the lower bits, and 240 Hz, missing-fundamental, strong-second-harmonic, glide and onset cases show no octave errors; the decoded level follows the input down 40 dB and silent or -70 dBFS input is coded as voice frames that decode near silence; the 16 DTMF digits, the call-progress tones and single tones at several phases are sent as tone frames in DVSI's layout with the index the decoder plays at their frequencies, decode within 0.5 dB of the input level and leave the prediction history unchanged, KNOX tones stay voice, voice-tone-voice prediction stays with the decoder's, a malformed history is rejected on a tone frame without touching state, and call-progress tones are timed as below (`test_ambe2400_encoder`)
- AMBE+2 2450 encoder: frame FEC round trips of all 49 bits and single-error correction; decoded frames requantized to the same model (gain, prediction, block DCTs, PRBA and higher-order quantizers, packing); the TIA-102.BABA-1 eq 4 voicing quantizer against an independent evaluation; DTMF, KNOX, call-progress and single tones at several phases sent as tone frames the decoder plays at the input level, the DTMF frequency and twist limits, voice-tone-voice prediction, and a malformed history rejected on a tone frame without touching state; call-progress timing, checked for the D-STAR encoder too: a tone is sent once it has filled three spans and for two frames after, 45 ms bursts at any alignment stay voice, a frame rejected during that confirmation leaves the stream unchanged, a 20 ms dropout is bridged and a DTMF digit that follows is sent at once; tone edges, also checked for the D-STAR encoder: 60 ms DTMF and single tones are sent for three or four frames wherever they start, and a DTMF digit that changes directly to another leaves no voice frame between them wherever the change falls, also at a level just above the detector's floor; plus the checks shared with IMBE below (`test_ambe2450_encoder`)
- IMBE 7200x4400 encoder: the (15,11) Hamming encoder over every data word, frame FEC round trips of all 88 bits and single-error correction, the alternating synchronization bit, and decoded frames of every harmonic count requantized to identical bits (`test_imbe7200_encoder`)
- Shared by both (`tests/encoder_test_support.h`): `cur_mp` bit-identical to the decoder chain with `prev_mp` untouched, rejected samples leaving the stream unchanged, prediction boundary taps, 16-bit/float parity, independent contexts and reset replay, noise staying unvoiced, steady vowels from 70 to 310 Hz voiced and on pitch, no octave errors (240 Hz, missing fundamental, strong 2nd harmonic, glide, and 133, 200 and 390 Hz harmonic series), the decoded level following the input down 40 dB, quiet input coded as voice, and a one-frame flush
- MBE speech analysis numerics: Kaiser windows, the window transform, harmonic fits and magnitudes of off-grid harmonics, noise statistics without magnitude ripple, the FFT autocorrelation against the direct sum, and the voicing thresholds (`test_speech_analysis`)
- Encoder context, FFT and autocorrelation plan allocation failures (each of the eleven, for the D-STAR, AMBE+2 and IMBE encoders) clean up fully; after successful allocation, encoding (tone frames included) and reset allocate nothing (`test_ambe2400_encoder_oom`, `test_ambe2450_encoder_oom`, `test_imbe7200_encoder_oom`, GNU link wrapping when LTO is disabled)
- noise determinism and frame-state determinism checks
- parameter and synthesis behavior checks, including an `L * w0 < pi` bound for every model the decoders emit
- AMBE 3600x2450 frame-type handling per TIA-102.BABA-1: frame classification, silence history freeze, repeats, erasures, tones (including unverified tone frames), muting, and recovery after mutes and tones (`test_ambe2450_frame_types`)
- D-STAR frame classes outside plain voice: b0 125, tone frames without a tone and b0 127 frames repeating like error repeats through the hard, soft and parameter-bit APIs and sharing the repeat limit, the null frame and b0 120..124 as voice, zero-volume tones, tones that fail the error check playing as lost frames, tones muted with the error rate, the repeat-limit mute (fade-out, DVSI-level noise, prediction kept) and the error-rate mute level, and the all-ones DV frame as a lost frame for runs of 1, 3 and 20, and the 80-sample crossfades into and out of a tone against a stream with a lost frame in its place (`test_ambe2400_frame_types`)
- float-to-int16 conversion parity checks
- golden PCM hash regression checks, by kernel ISA: `test_golden_pcm` runs the CPU probe's choice,
  `test_golden_pcm_sse2` and `test_golden_pcm_avx2` force either copy of a run-time AVX2 dispatch build
  and skip (exit 77) where the build or CPU cannot run it; `MBELIB_TEST_REQUIRE_AVX2=ON` makes the AVX2
  skip a failure
- SIMD architecture detection checks, including where the AVX2 dispatch applies
- run-time AVX2 dispatch (`test_cpu_dispatch`): the CPU probe against the compiler's CPU-feature
  runtime, the ISA selection rules, voiced oscillators that differ between SSE2 and AVX2 only within
  tolerance, tone fits within 1e-9, and soft ECC decoders and the TETRA adaptive codebook bit-identical
  across ISAs; `test_simd_numerics`, `test_tetra_simd` and `test_tetra_decoder` repeat their oracles and
  stream hashes under every ISA the build and CPU can select
- the AVX2 translation units define only `*_avx2` externals, so no AVX2 code can be linked into a
  baseline caller (`test_avx2_symbols`, Linux, non-LTO configurations)
- the host-tuning rules in CMake script mode, then configure-only checks of per-configuration
  settings (`test_build_configuration`; the configure checks need Ninja and skip when the forwarded
  compiler command selects a target ISA): every skip rule of `cmake/MbeNativeArch.cmake` (including
  a target ISA in the compiler command, quoted or tab-separated flags and MSVC's `CL`/`_CL_`) and the
  MinGW/Cygwin assembler requirement; Ninja Multi-Config, with or without `CMAKE_BUILD_TYPE` and
  from a source path with spaces, gets SIMD in both configurations (64-bit targets) and LTO only in
  Release; Release gets the native flag (when the host allows it) and Debug never does; a user
  `-march` or `MBELIB_ENABLE_NATIVE=OFF` turns it off; and a portable Release keeps its AVX2 copies
  when only Debug targets AVX2; one build directory reconfigured for `-m32` or (Clang) AArch64 probes
  the architecture again; a universal `CMAKE_OSX_ARCHITECTURES` off Apple keeps the dispatch; and the
  parent's `CMAKE_OSX_ARCHITECTURES` is forwarded. Exported build settings (`CFLAGS`, `CL`, `CMAKE_GENERATOR`, ...)
  are removed from the nested configures' environment
- quality-tool filesystem helper checks (`test_quality_fs`, built with `MBELIB_BUILD_TOOLS=ON`)
- quality-tool CLI alignment, accepted/rejected frame widths, and private output
  permissions (`test_quality_tools`, with tools enabled and Python 3 available)
- DVSI test-vector import: deinterleave bijection, TIA-102.BABA 7.5 code spacing,
  fixture round trips and preserved channel bits (`test_quality_dvsi`)
- TETRA decoder: bit unpacking (tables 3 and 4), every pitch index, pulse positions and absent
  pulses, sharpening, the 16-bit wrap of an overflowing codevector, LSP lookup, boundary widening
  and the rejection of crossed vectors,
  interpolation, LSP to A(z) against a double-precision expansion and against coefficients read
  back from ETSI's reference decoder, the synthesis filter, the
  adaptive codebook, log2/pow2, the gains against values ETSI's reference decoder produced on the
  test's own streams (each ruling out another plausible computation) and at the extremes
  (including a code gain past 16 bits), the floor of the code-vector energy, the 27/25 energy
  limits, concealment state (energy decay and floor, held LSPs, lag and indices, silence and the
  reference decoder's LSP state before the first good frame, ended by any good frame), the
  API contract (NULL arguments, rejected frames leaving state and output untouched, ignored
  bad-frame payloads, reset, interleaved contexts), and stream hashes asserted on every build
  (`test_tetra_decoder`); the saturating fixed-point operators at their edges
  (`test_fixed_point`); independent ordered arithmetic oracles for SIMD and scalar
  synthesis, fractional pitch, pulse shaping, excitation, energy and LSP interpolation
  (`test_tetra_simd`), including saturation, wraps, short-lag dependencies, unaligned
  buffers and tails; the TETRA entry points through the shared library
  (`test_tetra_shared_api`); TETRA argument and bit validation (`test_input_validation`); the
  table generators and the committed tables' shape and provenance (`test_tetra_tables`); the
  TETRA comparison tool, and with tools enabled `mbe_tetra_decode` (`test_tetra_conformance_tools`,
  `test_tetra_decode_tool`). `fuzz_tetra_decode` fuzzes the decoder through the public API.

PCM conversion is always compiled with IEEE semantics, including fast-math builds, so NaN/Inf handling is preserved.

Test assertions remain enabled in every configuration. MSVC force-includes
`tests/enable_assertions.h` to undefine `NDEBUG` before reading test sources,
avoiding the D9025 warning from combining `/DNDEBUG` and `/UNDEBUG`.
`test_assertions` checks that the assertion expression actually executes.

Run the default test suite with:

```sh
cmake --preset dev-debug
cmake --build --preset dev-debug -j
ctest --preset dev-debug --output-on-failure
```

## Continuous Integration

GitHub Actions runs tests and quality checks on pull requests and pushes to the
primary branch. Required checks include cross-platform builds, sanitizer tests,
static analysis, CodeQL (C/C++, workflows and Python), PR fuzzing, repository
security guardrails, workflow linting, dependency review, secret scanning, OSV
scanning, and install/consume checks.

Fast-math changes what the optimizer may assume, so the shipped configuration
gets jobs of its own. On Linux the whole suite, quality tools included, also
runs in `dev-release` builds (fast-math, LTO, SIMD) with x86-64-v3 code
generation, once with GCC and once with a current Clang from the Arch
toolchain. Clang folds bit tests on floats it computed (an exponent-mask
NaN/Inf test once passed overflowed encoder models on macOS this way) without
a diagnostic, so only a test that feeds the guard a NaN or infinity catches
that. With `-Werror`, its `-Wnan-infinity-disabled` does turn `isnan()`,
`isinf()` and NaN or infinity constants in fast-math code into a build error;
the job first probes that the compiler still enforces it. A third job builds and links everything for x86-64-v4 with
GCC and LTO and warnings as errors, since GCC reports some diagnostics only for
AVX-512 targets at LTO link time; the runner may lack AVX-512, so that build is
never run. To test for NaN or infinity in fast-math code, compare the
magnitude of the bit pattern against a bound, or compile the translation unit
with IEEE semantics as `mbe_pcm_convert.c` is.

Consumers do not build a preset. dsd-neo's CI, distribution packages and the
vcpkg port configure a plain `CMAKE_BUILD_TYPE=Release`, whose defaults are
SIMD on, fast-math off, warnings as errors and host tuning
(`MBELIB_ENABLE_NATIVE`) on; packagers that ship binaries pass
`-DMBELIB_ENABLE_NATIVE=OFF`. The IEEE-only SIMD batches
(`!defined(__FAST_MATH__)`, such as the AMBE codebook search) compile only
there, so a fourth Linux job builds and tests exactly that configuration with
the runner's GCC and fails if those defaults change. v2.5.0's codebook search
failed to build under it (GCC 13 `-Waggressive-loop-optimizations` on an `int`
row offset) while every preset job passed. The same job also builds and tests
the portable variant (`-DMBELIB_ENABLE_NATIVE=OFF`), which keeps baseline x86-64
code generation and dispatches the AVX2 kernels at run time.

It then configures Ninja Multi-Config without a build type, with LTO, and runs
the suite for Release and Debug, so a multi-config Release keeps SIMD, LTO and
the host tuning.

Release assets and the install-and-consume jobs build `dev-release` with
`-DMBELIB_ENABLE_NATIVE=OFF`, so the shipped libraries need only x86-64 v1
(SSE2) or the ARM64 baseline and use AVX2 through the run-time dispatch; the
install-and-consume jobs also run that build's suite (on Linux and Windows
with `MBELIB_TEST_REQUIRE_AVX2=ON`). The `dev-release` legs of the build matrix
keep the native default and also require the AVX2 golden test on Linux and
Windows. A Debug job builds a portable SIMD tree and checks the strict golden
hashes for both dispatch paths (`test_golden_pcm_sse2`, `test_golden_pcm_avx2`
with `MBELIB_TEST_REQUIRE_AVX2=ON`), then builds a Debug `-march=x86-64-v3`
tree and requires its `gen_golden` output to equal the dispatched
`gen_golden avx2`.

Every required check runs on every pull request. A required check that a path
filter keeps from starting never reports, so the pull request waits on it
forever. A job that only matters for some paths therefore runs on every pull
request and skips its work when `tools/ci_changed_files.sh` finds nothing for
it, as PR fuzzing does. A job that `if:` skips still reports success under its
name, so a required job that only makes sense on a pull request, such as
dependency review, lives in `guardrails-pr`, a workflow that only pull requests
start, with no event condition.

## Regression Test Requirement

At least 50% of bugs fixed in the last six months should include regression
tests. A pull request that fixes a bug should add a regression test unless:

- the behavior cannot be reproduced reliably in automation
- the fix is entirely documentation or packaging metadata
- a better guardrail exists, such as a static-analysis rule or workflow check

When no regression test is added for a bug fix, the pull request must explain
why.

## Major Functionality Test Requirement

Major new functionality must add or update automated tests. Major functionality
includes:

- public API changes
- codec or DSP behavior changes
- external input handling changes
- dependency changes that affect compiled code
- security-sensitive workflow or release changes
- installation and packaging behavior changes

## Coverage Target

The project target is at least 80% statement coverage for project-owned source
files when measured with a FLOSS C coverage tool. Vendored code under
`src/external/` is excluded from project-owned coverage accounting.

Coverage evidence should be generated from an instrumented Debug build and
recorded in a pull request, issue, or CI artifact before claiming coverage-based
badge criteria.

## Dynamic Analysis

For memory-safety-sensitive C changes, run sanitizer tests:

```sh
cmake --preset asan-ubsan-debug
cmake --build --preset asan-ubsan-debug -j
ctest --preset asan-ubsan-debug --output-on-failure
```

Frame-processing paths are also covered by ClusterFuzzLite PR fuzzing with
AddressSanitizer and UndefinedBehaviorSanitizer, including fixed-size hard/soft
frame decode paths and parameter synthesis paths. The fuzzers run on pull
requests that change the library, its headers, the fuzz harnesses or the fuzz
build (`fuzz_targets` in `tools/ci_changed_files.sh`).

## SIMD and replay benchmarks

Build the opt-in benchmarks with the library configuration being measured:

```sh
cmake --preset dev-release -DMBELIB_BUILD_BENCHMARKS=ON
cmake --build --preset dev-release -j
build/dev-release/bench_synth
build/dev-release/bench_synth steady mixed
build/dev-release/bench_synth gradual voiced
build/dev-release/bench_synth jump voiced
build/dev-release/bench_replay encode ambe2400 3 7 < speech.raw
build/dev-release/bench_replay hard ambe2450 3 7 < ambe2450.frames
build/dev-release/bench_replay soft imbe7200 3 7 mixed < imbe7200.frames
build/dev-release/bench_tetra 3 7 < /path/to/TE1_46.CHD
```

`bench_synth [legacy|steady|gradual|jump] [mixed|voiced]` defaults to legacy
mixed, the original workload, unchanged: w0 jumps between 0.09 and 0.11 at a
fixed L = 40 (harmonics above Nyquist are skipped by the synthesizer), with no
warm-up and state carried across runs. The quality A/B runner times
`bench_synth` without arguments against older builds, so that default must stay
the same work. In the other workloads every model keeps its harmonics below
Nyquist: steady and gradual pitch exercise interpolated low harmonics, jumps
exercise windowed synthesis, and each of ten measured runs starts with the same
state after one untimed warm-up. Every workload prints the `avg:` line the A/B
runner reads.

`bench_replay encode|hard|soft ambe2400|ambe2450|imbe7200|imbe7100
[repeats=1] [runs=7] [uniform|mixed]` loads standard input into memory, then
allocates the encoder and warms up before timing. Encoding takes mono 8 kHz s16le PCM, zero pads
a partial last frame and rejects partial samples; ProVoice has no encoder.
Decoding takes one rectangular, row-major bit row per line: 96 bits for AMBE,
184 for P25 or 168 for ProVoice, as written by `mbe_quality_reframe` (parameter
rows are not accepted). Soft replay uses confidence 255 or a reproducible mix
of all 256 confidence values, including zero. Use channel-error fixtures to
measure noisy decoding. Inputs are capped at 65,536 frames.

The replay benchmark reports median CPU microseconds per 20 ms frame, excludes
loading, allocation and stream resets, and checks a checksum against its warm-up
on every run. Repeats concatenate the input within a run; each run resets the
encoder, parameter triplet and RNG. Checksumming is included in the timing.
Pin runs to one CPU when possible and alternate matched baseline/candidate
builds. Compare identical compiler, fast-math, LTO and `MBELIB_ENABLE_NATIVE`
settings: a default Release build is tuned for the build host, so pass
`-DMBELIB_ENABLE_NATIVE=OFF` to measure the portable build. `-march=x86-64-v3`
is a separate compile target; portable SIMD builds still require only SSE2 on
x86 (running the dispatched kernels as AVX2 on AVX2 CPUs) and use NEON only
when enabled by the existing target selection. Native ARM
measurements are required before claiming NEON performance improvements.

`bench_tetra [repeats=10] [runs=7]` applies the same in-memory replay method to
TETRA's 30 ms frames. Standard input is CHD format: 138 little-endian 16-bit
words per frame (bad-frame indicator followed by 137 bits). It rejects empty,
truncated or invalid good-frame input, caps input at 65,536 frames, and ignores
bad-frame payloads as the public decoder does. The reported median excludes
loading, allocation, resets and one warm-up run. Repeats concatenate the input
without resetting; every measured run starts from reset and must reproduce the
warm-up checksum. One PCM sample per frame contributes to the timed checksum;
use `mbe_tetra_decode` and byte comparisons to check complete PCM identity.
Include speech, losses and overloaded synthetic streams in comparisons. See
[TETRA SIMD measurements](tetra-acelp.md#simd-optimization) for the measured
speed and bit-exactness checks.

The soft ECC codebooks are immutable 32-bit codewords (32 KiB in all) in the
original enumeration order, with the codeword bit of each data bit. Regenerate
them with Python 3.10+ using `python3 tools/gen_ecc_codebooks.py >
src/internal/mbe_ecc_codebooks.h`. Each decode builds small cost tables on the
stack (data bits by enumeration index, parity bits by value) and starts from
the cost of the codeword carrying the hard decision's data bits; it skips only
codewords that cost more, so the result is exactly that of a full scan, ties
included. SIMD builds screen candidates in groups of four (SSE2) or eight (AVX2),
then resolve surviving candidates in enumeration order with the same exact
integer scores and tie rules. AVX2 gathers pairs of 16-bit parity costs without
reading beyond the cost table. Other targets use the scalar search. The search
has no lazy initialization or shared mutable state; a portable x86-64 build
only reads the run-time dispatch's cached ISA choice, which every thread
computes identically. `test_ecc` checks every clean codeword, then compares
1024 pseudo-random inputs per code against an independent oracle that scores
every codeword, with zero, full, binary, single-bit, small tied, mostly
confident and random confidences.

`test_simd_numerics` compares tone fits to direct trigonometric projection and
an independent solver (relative tolerance 1e-9), and, where the target has
vector voiced kernels, voiced output to the scalar formulas (absolute
tolerances 3e-5 for unit-scale windowed components and 2e-4 for
interpolation). It exercises unaligned inputs and scalar tails. Targets
without vector kernels compile the historical voiced loops in `mbelib.c`
unchanged, so their PCM stays bit-identical, under x87 excess precision and
fast-math too; SIMD recurrence batching permits bounded PCM rounding
differences and has separate Debug goldens. Each SIMD interpolation lane
seeds from the phase's own cosine and sine, rotated by its offset in double,
so caller-supplied phases far beyond the decoder's range, up to near the
largest float, follow the exact phase. `mbe_voiced.c` is compiled with IEEE
semantics even in fast-math builds, as `mbe_pcm_convert.c` is, since a
fast-math math library need not reduce huge cos/sin arguments exactly. Phase and RNG updates stay outside
these kernels. Tone refinement
trials remain sequential and double precision; ARM32 uses scalar tone fitting.
Without fast-math, the pitch low-pass filter assigns independent output samples
to SIMD lanes; its 21 taps still accumulate in order. Fast-math builds retain
the compiler-vectorized filter/energy loop. IEEE AMBE codebook searches assign
complete candidate distances to lanes and resolve equal minima by original row
index; fast-math builds retain the compiler-vectorized codebook loop.
Neither kernel reduces partial floating-point sums across lanes. Remaining
analysis reductions retain their existing evaluation order.

`mbe_floattoshort` classifies NaNs in integer SIMD lanes before arithmetic,
converts them to zero, and clips infinities to the signed limit. This avoids a
separate scalar scan and keeps exceptional frames vectorized. Its tests include
arbitrary IEEE bit patterns, signaling NaNs, unaligned inputs/outputs and frame
boundary canaries. Hard/soft bit validation uses compiler-vectorized integer OR
reductions in SIMD builds, with tests for every invalid byte value and tail
length. All three frame encoders use the same hard-bit validator. Reliability
bytes are excluded from validation.

Unvoiced noise generation advances four independent integer LCG lanes (eight
on AVX2), preserving the original sample order and stream state exactly.
`test_simd_numerics` compares every one of the 53,125 states with a scalar
recurrence across three consecutive frames, including the overlap and next
seed. This oracle also runs on ARM32 and AArch64 under QEMU.

On AVX2 compile targets, float kernels process eight lanes, including voiced
oscillators, unvoiced windowing/overlap and PCM conversion. Portable x86-64
builds still require only SSE2; their voiced oscillators, tone fits, soft ECC
search and TETRA adaptive codebook also have an AVX2/FMA copy that runs on
CPUs with the x86-64-v3 set (see
[run-time dispatch](build-installation.md#run-time-avx2fma-dispatch)); the
other kernels keep SSE2 widths. Debug AVX2 has its own float and AMBE sequence
PCM goldens because wider oscillator recurrences change rounding, and the
dispatched AVX2 copy reproduces a Debug `-march=x86-64-v3` build exactly.
Regenerate them in a Debug SIMD build with `gen_golden sse2` and
`gen_golden avx2` (`auto`, the default, prints the probe's choice), and check
`gen_golden avx2` against `gen_golden` from a Debug build with
`-DCMAKE_C_FLAGS=-march=x86-64-v3`.

`tools/bench_compare.sh` runs the encoder, PCM conversion, unvoiced synthesis
and all six steady/gradual/jump mixed/voiced workloads in matching Release
builds with explicit SIMD disabled, enabled (both portable, so the SIMD build
dispatches AVX2 on AVX2 CPUs) and enabled with host tuning. Compilers can still automatically
vectorize C loops when explicit SIMD is disabled. Use `bench_replay` in those
same trees for real speech, all four decoders, tone frames and channel errors;
synthetic synthesis timings alone do not establish end-to-end gains. It also
builds `bench_tetra`; supply local CHD inputs to compare TETRA in the same trees.

The October 2026 SIMD audit compared against `a962d3d` on a Ryzen 7 5800HS,
with GCC 16.2.1, Release, fast-math and LTO enabled on both sides. These are
medians of three alternating comparisons pinned to one CPU; lower is better.
Replay uses the DVSI `dam` speech input. SSE2 and AVX2 are separate compile
targets, each compared with the same target before the changes.

| Workload | Target | Before | After | Time reduction |
| --- | --- | ---: | ---: | ---: |
| PCM conversion, 10 million frames | SSE2 | 0.6210 s | 0.2845 s | 54% |
| PCM conversion, 10 million frames | AVX2 | 0.2800 s | 0.2033 s | 27% |
| Mixed synthesis, steady pitch, 2,000 frames | AVX2 | 9.751 ms | 8.531 ms | 13% |
| Mixed synthesis, pitch jumps, 2,000 frames | AVX2 | 9.612 ms | 7.761 ms | 19% |
| D-STAR hard-frame replay | AVX2 | 3.726 us/frame | 3.466 us/frame | 7% |
| AMBE+2 hard-frame replay | AVX2 | 4.025 us/frame | 3.768 us/frame | 6% |
| IMBE hard-frame replay | AVX2 | 4.101 us/frame | 3.889 us/frame | 5% |

Explicit codebook batching reduced the complete AMBE residual quantizer's time
by 9–17% with GCC and IEEE math. Its fast-math version measured slower than the
compiler-vectorized loop, so explicit batching is enabled only with IEEE math.
Fast-math encoder timings stayed near their previous values. With Clang 23.1.1,
SSE2 and neither fast-math nor LTO, the pitch-filter changes reduced whole-frame
encoding time by about 7% for the AMBE codecs and 10% for IMBE. FFTs, tone
fitting and analysis reductions already had SIMD or effective compiler
vectorization; changes there require their own measurements. ARM numerical
checks run under QEMU establish correctness, not native ARM speedups.

For the same audit, all three encoders produced identical parameter-bit rows
on all 79 cached DVSI inputs (27 speech and 52 tone vectors): 256,707 frames
per build pair. This held for GCC SSE2/AVX2 with fast-math and LTO, GCC SSE2
with IEEE math, Clang AVX2/IEEE SSE2 configurations, and Android NDK Clang 19
IEEE builds for ARM32 and AArch64 under QEMU. Compare each compiler and build
configuration against its own baseline; cross-compiler identity is not the
acceptance criterion.

Decoder comparisons used every cached DVSI frame stream, a fixed seed and the
baseline's alignment against both the original speech and DVSI's decoded PCM.
SSE2 output was byte-identical. Across GCC and Clang AVX2 builds, each changed
int16 sample differed by at most one unit, with no added clipping or nonfinite
samples. The largest absolute change in log-spectral distance was below
0.00018 dB against the input and 0.00007 dB against DVSI; envelope correlation
changed by less than 0.000008. Preserve these per-codec comparisons when
changing SIMD kernels: passing synthetic numerical tests alone does not
establish speech quality. Keep DVSI vectors and generated audio local.

The subsequent noise-generation and frame-validation changes were also
compared directly with the earlier PR head, `c204ce4`. Across GCC SSE2/AVX2
with fast-math and LTO, and Clang 23.1.1 SSE2 with IEEE math, all 79 DVSI inputs
produced identical encoder parameter bits, packed frames and decoded PCM:
770,121 compared frames and 123,105,600 PCM samples, with no nonfinite output.
The established DVSI quality scores therefore remain unchanged by these two
optimizations. Native ARM performance remains a separate measurement.

For this follow-up, five alternating comparisons pinned to one CPU on the
same Ryzen/GCC setup measured the following additional reductions relative
to `c204ce4`. The first three codecs replay the DVSI `dam` stream through
`bench_replay hard <codec> 10 7`; ProVoice uses a converted frame fixture.

| Hard-frame replay | Target | Before | After | Time reduction |
| --- | --- | ---: | ---: | ---: |
| D-STAR | SSE2 | 3.954 us/frame | 3.791 us/frame | 4.1% |
| AMBE+2 | SSE2 | 4.305 us/frame | 4.145 us/frame | 3.7% |
| IMBE | SSE2 | 4.394 us/frame | 4.238 us/frame | 3.6% |
| ProVoice | SSE2 | 4.721 us/frame | 4.568 us/frame | 3.2% |
| D-STAR | AVX2 | 3.493 us/frame | 3.140 us/frame | 10.1% |
| AMBE+2 | AVX2 | 3.802 us/frame | 3.456 us/frame | 9.1% |
| IMBE | AVX2 | 3.911 us/frame | 3.576 us/frame | 8.6% |
| ProVoice | AVX2 | 4.109 us/frame | 3.778 us/frame | 8.1% |

`bench_unvoiced 36 20000 7` took 7.2% less time on SSE2 and 20.8% less on
AVX2. A separate frame-packing microbenchmark cycled 256 deterministic input
bit arrays through two million calls to each `mbe_encode*Frame` API: AMBE
packing took 17–19% less time (about 14–15 ns/frame), and IMBE took 5–9% less
(about 12–25 ns/frame). These packing figures exclude PCM analysis and
parameter quantization; whole-frame decoding above establishes the gain
from batching the noise generator. With Clang 23.1.1 SSE2 and IEEE math,
the same replay comparisons took 6.6–7.9% less time.

### PR #130 synthesis and ECC optimization

The follow-up to the DVSI parity changes compares with PR #130 at `d2293ba`.
It batches the spectral-envelope phase filter across independent harmonics
using SSE2, AVX2 or NEON, keeping each harmonic's ten taps in order. Synthesis
computes the high-pass phase only for voiced oscillators, including harmonics
fading between different model lengths. On SIMD targets the gain loop is
separate from the phase calls, allowing compiler vectorization. The scalar
path retains its combined expression, including x87 rounding. The required
phases still use the double-precision `atan2` law.

AMBE+2's counted C1 path re-encodes the data from its first Golay decode to
count parity errors, removing the second decode and, for soft input, the
second codebook search. It saves the received word before packing output,
preserving callers that put the 49 output bits at the start of the input
frame. Encoder analysis, quantization, FFTs and oscillator kernels retain
the implementations from the earlier SIMD work.

These measurements use a Ryzen 9 5950X, GCC 16.2.1, Release with fast-math
and LTO on both sides. Each figure is the median of five alternating
baseline/candidate comparisons pinned to CPU 6. Each decode invocation uses
`bench_replay hard|soft <codec> 10 5` on DVSI's `dam` stream; ProVoice uses
the same speech encoded to IMBE parameter rows and reframed as ProVoice.
Input loading, allocation and warm-up are outside the timer. Lower is better.

| Complete frame decode | Target | Before | After | Time reduction |
| --- | --- | ---: | ---: | ---: |
| D-STAR, hard | SSE2 | 4.541 us | 3.739 us | 17.7% |
| AMBE+2, hard | SSE2 | 4.991 us | 4.130 us | 17.3% |
| IMBE, hard | SSE2 | 4.957 us | 4.155 us | 16.2% |
| ProVoice, hard | SSE2 | 5.314 us | 4.532 us | 14.7% |
| AMBE+2, soft | SSE2 | 7.502 us | 5.819 us | 22.4% |
| D-STAR, hard | AVX2 | 3.191 us | 2.964 us | 7.1% |
| AMBE+2, hard | AVX2 | 3.519 us | 3.283 us | 6.7% |
| IMBE, hard | AVX2 | 3.574 us | 3.349 us | 6.3% |
| ProVoice, hard | AVX2 | 3.761 us | 3.555 us | 5.5% |
| AMBE+2, soft | AVX2 | 5.948 us | 4.947 us | 16.8% |

The other soft decoders take 8.5–12.8% less time on SSE2 and 2.5–4.7% less
on AVX2. Scalar Release with IEEE math takes 8.1–9.3% less time for hard
decoding. Encoder medians on SSE2 and AVX2 stay within 0.2% of baseline.
With Clang 23.1.1, hard decoding takes 14.5–18.0% less time on SSE2/IEEE
and 17.3–21.5% less on AVX2 with fast-math and LTO. Longer Clang AVX2
encoder runs put AMBE medians within 0.5% of baseline, inside the observed
timing spread, and IMBE 1.8% faster.
Speedups depend on the compiler, target and speech; QEMU correctness checks
do not establish native ARM performance.

For each of six host build pairs, all three encoders produce identical
parameter-bit rows on all 97 cached DVSI inputs: 296,772 encoded frames.
The configurations are GCC scalar/IEEE, GCC SSE2/IEEE, GCC SSE2 and AVX2
with fast-math and LTO, Clang SSE2/IEEE, and Clang AVX2 with fast-math and
LTO. All 321 cached frame streams, including tones, synthetic inputs and
hard/soft channel-error vectors, contribute 55,073,120 compared PCM samples
per pair. Five configurations produce identical PCM. GCC SSE2 with
fast-math changes 298 samples, each by one int16 unit, as the compiler
vectorizes the gain calculation. There are no nonfinite samples, additional
clips, or changed tone, erasure, repeat, silence, mute or error counts.

For the changed SSE2 streams, scoring both outputs at the baseline's lag
changes log-spectral distance by at most 0.000004 dB against the input and
0.000002 dB against DVSI; envelope correlation changes by less than
0.0000004. Tone detector agreement, wrong-tone fraction and steady-tone
window counts are unchanged. The parity PR's existing tone-gate result
remains as recorded in [DVSI parity](dvsi-parity.md#gates). The converted
ProVoice fixture also produces identical PCM in all six host configurations.

ARM32 and AArch64 NEON each compare the same 321 streams under QEMU.
ARM32 PCM is identical. AArch64 changes 6,045 of the 55,073,120 samples,
each by one int16 unit. Its fixed-lag log-spectral distance changes by less
than 0.00004 dB against both the input and DVSI, and envelope correlation
by less than 0.000003. Tone detector agreement, wrong-tone fraction and
steady-tone window counts are unchanged; there are no additional PCM rail hits.

The phase oracle covers every supported model length (1–56), mirrored and
extrapolated envelope edges, zero magnitudes and both voicing classes. ECC
checks cover 512 arbitrary received words with hard and soft decisions,
including zero-reliability ties, and 32 in-place layouts. Scalar, SSE2 and
AVX2 Debug golden hashes are unchanged. Checks include the full Debug and
Release suites, scalar and AVX2 ASan/UBSan, NOTONES, x87 Debug, the six
Release configurations above, ARM32 and AArch64 NEON numerical/API tests
under QEMU with Android NDK Clang 18, and `tools/quality_preflight.sh`.
ARM32 NDK builds disable hardening because that compiler accepts the
stack-clash flag in CMake's probe but rejects it when compiling for ARM32.
ARM32 QEMU runs use a PID namespace to stay within Bionic's 16-bit PID limit.
Gitleaks history and worktree-diff scans are also clean; the history scan
mounts the worktree's common Git directory as well as the worktree.
DVSI inputs, audio, comparison reports and benchmark logs stay local.

## TETRA Conformance

EN 300 395-4 tests a TETRA speech decoder bit for bit: TE1_46.CHD (138 little-endian 16-bit
words per frame: the bad-frame indicator, then B1..B137) decoded must equal TE2_128.PCM. The
files are in `orignl.zip` inside ETSI's en_30039504v010301p0.zip. Keep them, and any build of
ETSI's reference decoder and its output, outside the checkout; never commit them.

```sh
cmake --preset dev-debug -DMBELIB_BUILD_TOOLS=ON
cmake --build --preset dev-debug -j --target mbe_tetra_decode
./build/dev-debug/mbe_tetra_decode TE1_46.CHD /tmp/ours.pcm
python3 -I tools/quality/tetra_conformance.py compare /tmp/ours.pcm TE2_128.PCM
```

`compare` exits 0 only for identical files, as TE1_46.CHD now gives; `--gate` instead takes the
exit status from tolerance thresholds, for runs that may differ (a candidate change, say);
`--json` writes the metrics. For concealment,
`tetra_conformance.py bfi TE1_46.CHD /tmp/lossy.chd --seed 1 --rate 3` (also `--burst START:LENGTH`
and `--first N`) marks frames bad; decode the copy with both decoders and compare. Measured
results are in [docs/tetra-acelp.md](tetra-acelp.md#measured-results).

## Speech Quality Evaluation

Developer-only tooling; no new library runtime dependency. CI builds the two
quality tools and tests their filesystem helper; full quality runs remain opt-in.
The C99 `mbe_quality_eval` uses the installed public API and an independent 256-point FFT.
`mbe_quality_reframe` is a statically linked private fixture generator. A/B runs
additionally require Bash, Python 3, git, gcc/g++ (C++17), `patch`, and
coreutils. Fetching held-out speech requires network access and `ffmpeg`.
Offline reports require NumPy, SciPy, matplotlib, and pystoi; install these in a
local virtual environment rather than making them library dependencies.

Quality-tool outputs are created with owner read/write permissions on POSIX
(`0600`, further restricted by the caller's umask). Reopening an output truncates
its contents while retaining its existing permissions. Windows uses CRT
read/write permissions together with the directory's inherited ACL.

### Windows

With `MBELIB_BUILD_TOOLS=ON`, `mbe_quality_eval` and `mbe_quality_reframe`
build with MSVC or MinGW for local single measurements and fixture generation.
The evaluator always links the shared library; static linking would embed one
implementation and prevent measuring a different library at runtime.

The A/B runner, calibration, and OP25 encoder build remain Linux-only. The
comparison workflow uses `LD_TRACE_LOADED_OBJECTS` to verify which shared
library was actually loaded; Windows DLL placement does not reproduce that
provenance check. Windows evaluator paths must be ordinary ASCII drive paths;
non-ASCII names, UNC, `\\?\` and `\\.\` namespaces, alternate data streams and
leaves ending in a dot or space are rejected.

### Reproducible comparison workflow

```sh
cmake --preset dev-debug -DMBELIB_BUILD_TOOLS=ON
cmake --build --preset dev-debug --parallel
tools/quality/build_op25_encoder.sh
tools/quality/fetch_speech_corpus.sh --set legacy
python3 -m venv build/quality/venv
build/quality/venv/bin/python -m pip install numpy scipy matplotlib pystoi

# Preserve source revisions in external worktrees; do not reuse a different revision.
git worktree add --detach /tmp/mbelib-quality-original 5fd3f3dd738b76bc9a8446ff4d7ae99e565e35e0
git worktree add --detach /tmp/mbelib-quality-checkpoint 8da251d105e50815514b1eee994c2c8471469654
cmake -S /tmp/mbelib-quality-original -B build/quality/original-debug \
  -DCMAKE_BUILD_TYPE=Debug -DMBELIB_BUILD_TESTS=OFF -DMBELIB_BUILD_EXAMPLES=OFF
cmake --build build/quality/original-debug --parallel
cmake -S /tmp/mbelib-quality-checkpoint -B build/quality/checkpoint-debug \
  -DCMAKE_BUILD_TYPE=Debug -DMBELIB_BUILD_TESTS=OFF -DMBELIB_BUILD_EXAMPLES=OFF
cmake --build build/quality/checkpoint-debug --parallel

tools/quality/run_quality_ab.sh build/quality/original-debug build/quality/original-debug \
  --out-dir build/quality/runs/identity-debug
tools/quality/run_quality_ab.sh build/quality/original-debug build/dev-debug \
  --out-dir build/quality/runs/original-to-final-debug
tools/quality/run_quality_ab.sh build/quality/checkpoint-debug build/dev-debug \
  --calibration-libdir build/quality/original-debug \
  --out-dir build/quality/runs/checkpoint-to-final-debug
build/quality/venv/bin/python tools/quality/analyze_quality.py report \
  --run-dir build/quality/runs/original-to-final-debug
```

The runner accepts `baseline_libdir candidate_libdir [corpus_dir]`, followed by
`--out-dir`, `--evaluator`, `--reframer`, `--calibration-libdir`, and
`--gain-profile`. Without `--out-dir`, it creates a unique directory under
`build/quality/runs/`; an explicit nonempty output directory is rejected.
Historical `build/quality/out/` and `build/quality/frames/` are not overwritten.
Current cache entries live under `build/quality/frame-cache-v3/`.

Both libraries must have compatible installed signatures/state layouts and
matched compiler/DSP options. Loader tracing checks the actual selected library,
including inode identity, and rejects `LD_PRELOAD` and unexpected extra libmbe
libraries. A sanitizer evaluator must use its instrumented library for both
operands, never loader-swap an uninstrumented baseline into that process.

Each schema-2 `manifest.json` binds library/tool/build identities, seed, flushing,
references, encoded/reframed bits, calibration records, WAVs, metric JSON, command
logs, and tables. Inputs and reports are copied into the run directory and indexed
by SHA256; analysis does not depend on the mutable frame cache. Source revision
labels require a captured matching binary hash; a checkout or directory name
alone does not prove binary provenance. Unknown source identities remain explicit.
`build/quality/checkpoint-capture.json`, when present, supplies these SHA-bound
comparison points. Same-library runs require exact PCM, frame-count, and complete
metric/null/support equality. Original-only identity controls can pass that gate
without establishing the later corrected-tone/adaptive regression gates.

### Encoder provenance, frame paths, and calibration

- Codec2 revision: `310777b1c6f1af0bc7c72f5b32f80f6fd9136962`
  (<https://github.com/drowe67/codec2>). The named **legacy/development** set is
  `hts1a`, `hts2a`, `hts1`, `kristoff`, and `ve9qrp_10s`. `hts1a` is exactly the
  first 24,000 samples of `hts1`; both remain available, but the duplicate is
  excluded from independent aggregate statistics when its parent is present.
- OP25 revision: `71abcd0ead32f86f51615ea6cc8a6a4dba4c949a`
  (<https://github.com/boatbod/op25>). Its GPL-3.0-or-later encoder and bundled
  mbelib remain in a separate executable. Never mix their incompatible
  `mbe_parms` with mbelib-neo. `op25_mbelib_gain.patch` changes disposable build
  copies only. Its AMBE2400 bridge is a decoder-targeted fixture, not an
  interoperable D-STAR transmitter or independent wire-conformance proof.
- IMBE7200 is encoded/calibrated once, then its canonical 88-bit parameters are
  reframed into both 184-bit IMBE7200 and 168-bit IMBE7100 rectangular frames.
  AMBE2450's 49-bit parameters become 96-bit frames; OP25's existing 96-bit
  AMBE2400 output is retained. Reports exercise all four actual `Framef` APIs.
  Canonical IMBE `Dataf` is an exact-PCM equality control, not a fifth independent
  codec. Reframed IMBE7100 is not a native transmitter/capture validation.
- `--gain-adjust` is log2 encoder attenuation, default `log2(7)`. Per-clip
  calibration uses only the calibration library, targeting decoded active RMS
  −3 ± 0.5 dB relative to reference, with at most six passes and lag 0…640.
  Schema-2 speech support, finite alignment/level, full reference coverage,
  and zero nonfinite float samples are mandatory. Invalid calibration fails;
  silent controls bypass calibration rather than selecting a fallback gain.
- `--flush-frames N` accepts integers 0…50, default 5. After the final nonempty
  input block, the encoder submits exactly N zero PCM blocks through its normal
  output path; ordinary partial-block padding remains. Empty input fails.
  Original sample counts and flushing are separately recorded and cache-bound.
  Five blocks cover alignment without discarding the reference tail; this is
  encoder fixture completion, not added decoder latency.

### Evaluator interface and measurements

```sh
build/quality/op25_encode --mode imbe7200 \
  --in build/quality/corpus/hts1a.raw --out build/quality/example.data.txt \
  --gain-adjust 3.4 --flush-frames 5
build/dev-debug/mbe_quality_reframe --codec imbe7100 \
  --in build/quality/example.data.txt --out build/quality/example.frames.txt
build/dev-debug/mbe_quality_eval --codec imbe7100 \
  --frames build/quality/example.frames.txt --out build/quality/example.wav \
  --ref build/quality/corpus/hts1a.raw --seed 0x12345678 \
  --json build/quality/example.json
build/dev-debug/mbe_quality_eval --decoded build/quality/example.wav \
  --ref build/quality/corpus/hts1a.raw --lag-samples 456 \
  --json build/quality/example.pcm.json
```

Use baseline calibration or a frozen profile for comparisons, not the example
gain. Frame rows are literal ASCII 0/1, including rectangular padding:
IMBE7200 accepts 88-bit Dataf or 184-bit Framef; IMBE7100 only 168-bit Framef;
either AMBE mode accepts 49-bit Dataf or 96-bit Framef. One input file uses one
selected public API, reported as `input_api`. The reframer accepts parameter
rows only and validates everything before changing its destination.

`--decoded input.raw|input.wav` is measurement-only: it rejects `--codec`,
`--frames`, `--out`, and `--seed`, and never rewrites PCM. Raw/WAV inputs must
be mono 8-kHz s16le; WAV chunks are parsed rather than assuming a fixed header.
Output WAVs preserve pre-normalization PCM. Invalid input, aliases (including
hard links and normalized not-yet-created output paths), or fewer than 256
aligned samples fail with exit 2 **before output creation**. A negative codec
status or nonfinite float PCM exits 3 with the frame index.

Schema-2 metrics preserve the original supported formulas:

- `lag_samples` is the lag actually used. Optional `--lag-samples` accepts
  −160…800. Automatic alignment always runs independently: `auto_lag_samples`,
  nullable `alignment_corr`, and `alignment_at_limit` describe its estimate.
  It correlates log-RMS envelopes from centered 160-sample windows at 8-sample
  hops. Correlations within `1e-12` are tied, favoring the smallest absolute lag.
  Positive lag means decoded speech is delayed. The runner measures the
  candidate at the baseline lag, preserving common reference support/activity
  and the synthesis-join grid; differing automatic lags generate a separately
  labelled auto-alignment sensitivity report, not a claimed decoder-delay change.
- `reference_samples`, `decoded_samples`, `aligned_samples`, both signals'
  `*_trim_head`/`*_trim_tail`, `active_frames`, `spectral_frames`, and `join_count`
  expose support. Active complete 160-sample reference frames have positive
  energy within 40 dB of maximum. Boundary alignment or incomplete reference
  coverage invalidates a comparison.
- `reference_state` is exactly `no_reference`, `silent_reference`,
  `silent_decoded`, or `speech`. Speech requires positive active energy in both
  signals. Undefined reference metrics/correlation serialize as JSON `null`,
  never favorable zero LSD or zero crest error. Pure silence succeeds as an
  explicitly inapplicable speech-quality measurement.
- `level_offset_db` is decoded/reference active RMS before normalization.
  Reference-based measurements then use active-RMS-matched int16-scale signals.
  `lsd_db` is mean active-frame log-spectral distance: symmetric 256-point Hann,
  80-sample hop, bins 2…118, reference-derived power floor
  `max(1e-6 * max(reference power), 1e-20)`.
- `env_corr` correlates log 40-sample RMS envelopes, floored 60 dB below each
  maximum. `crest_ref_db`, `crest_dec_db`, and `crest_delta_db` are mean 20-ms
  peak/RMS values and decoded-minus-reference difference. `band_delta_db_*`
  measures active STFT power ratios in 0–500, 500–1000, 1000–2000, 2000–3000,
  and 3000–4000 Hz, plus overlapping 0–250 Hz and 3500–4000 Hz bands; the
  top one covers what `lsd_db` (up to 3687.5 Hz) leaves out. Without reference, decoded absolute
  `band_db_*` remains available where supported.
- `join_dec_db`/`join_ref_db` compare squared across-join derivatives with the
  eight neighboring derivatives on the same reference mask and actual retained
  synthesis grid; `join_excess_db` is their difference. `boundary_index_db`
  compares unshifted boundary-region and mid-frame derivative energy. These are
  **diagnostics, not standardized click detectors**: a continuous 75-Hz sine
  produces boundary index +16.61 dB and zero reference-matched join excess.
- `float_nonfinite_samples`, `float_peak`, and `float_clip_samples` inspect
  **post-library-limiter, pre-int16** PCM, using bitwise finite classification
  even under fast-math and the existing `(32767 * 0.95) / 7` threshold. They do
  not reveal inaccessible pre-limiter peaks. PCM-input mode reports them null.
  `pcm_peak`, `pcm_rail_samples` (`abs(sample) >= 31128`), and
  `pcm_max_rail_run` are measured before normalization. Tone/erasure/repeat/mute/
  silence frame counts and C0/protected/C4/total error sums retain process-result
  context; clips with any of them are not clean speech for the checkpoint gate.

### Held-out speech and fixed operating points

```sh
tools/quality/fetch_speech_corpus.sh --set heldout
tools/quality/run_quality_ab.sh build/quality/original-debug build/dev-debug \
  build/quality/corpus-heldout \
  --gain-profile build/quality/runs/original-to-final-debug/gain-profile.json \
  --out-dir build/quality/runs/heldout-debug
build/quality/venv/bin/python tools/quality/analyze_quality.py report \
  --run-dir build/quality/runs/heldout-debug
```

The fetcher preserves `arctic_a0001.wav`…`arctic_a0010.wav` for each of
`bdl`, `rms`, `clb`, and `slt` from
<http://festvox.org/cmu_arctic/cmu_arctic/>: 40 US-English recordings from two
male and two female speakers. Original mono16-kHz/s16 WAVs and each `COPYING`
remain under `build/quality/arctic-src/`; manifests bind source/converted hashes,
conversion argv and ffmpeg version. There is no loudness normalization.
`build/quality/corpus-heldout-minus12db/` contains
`round(sample * 10**(-12/20))` copies with preserved notices and native hashes.
Missing/invalid assets fail preparation; development speech is not substituted.

For a complete original-main native legacy calibration run, `report` emits
`gain-profile.json`: schema 1, calibration library SHA256, twelve source-manifest
hashes, and the per-mode median `gain_adjust_log2` across
`hts1`, `hts2a`, `kristoff`, `ve9qrp_10s`. Fixed-profile runs validate and preserve
these exact bytes, bypass per-clip calibration, and use the same gains across
native/−12 dB levels and Debug/Release. IMBE7100 inherits IMBE7200. Held-out
failures are findings, not permission to retune gains or phase constants.

Repeat the command with the attenuated corpus and a distinct output directory.
Release uses `--evaluator build/dev-release/mbe_quality_eval` and
`--reframer build/dev-release/mbe_quality_reframe`, matched Release library
operands, and the same profile; encoded-frame hashes must match Debug.

### Controls, acceptance, and listening

```sh
build/quality/venv/bin/python tools/quality/analyze_quality.py controls \
  --evaluator build/dev-debug/mbe_quality_eval --out-dir build/quality/control-run
build/quality/venv/bin/python tools/quality/analyze_quality.py benchmark \
  --baseline build/quality/original-release/bench_synth \
  --candidate build/dev-release/bench_synth --out build/quality/benchmark-original-final.json
```

Controls exercise the real PCM-input evaluator: identity, exact ×2 gain, ±64
sample delays, silence/silence, speech/zero, too-short input, aliases, clipping,
continuous 75-Hz sine, and increasing alternating frame-boundary steps. NumPy
independently cross-checks LSD on controls and real speech. The benchmark
alternates seven process runs per variant, records stdout/build/binary/library
identities, and compares median reported `avg` CPU seconds; the cost gate is
≤10% increase, not an algorithmic-delay measurement.

`report` validates immutable snapshots and writes `analysis.json`,
`acceptance.json`, per-clip/grouped CSV tables, shared-support waveform and
reference-scaled spectrogram figures, harmonic detail, band differences, and
zoomed joins for every flagged clip. STOI/ESTOI use the complete aligned overlap,
not concatenated active snippets, and are intelligibility diagnostics—not MOS
or naturalness estimates.

Acceptance separates `correctness`, `objective_nonregression`, and `perceptual`,
each `pass`, `fail`, or `not_established`, for every actual mode:

- Correctness requires supported/full-coverage measurements, no nonfinite PCM,
  exact canonical/frame equivalence, same-library identity where applicable,
  and candidate-bound native regression evidence for corrected tones,
  nonnegative attenuation, and warmed FFT/WOLA. Checkpoint comparisons also
  require unchanged clean-speech PCM.
- Objective screens use independent per-mode clip aggregates:
  mean absolute candidate-minus-baseline LSD change ≤0.3 dB; mean envelope
  correlation decrease ≤0.01; absolute mean per-band change ≤0.5 dB; candidate
  mean join excess ≤`max(0, baseline mean) + 0.1` dB; and **mean per-clip absolute
  crest error must not increase**. Absolute signed-mean crest error is reported
  separately. Every per-clip exceedance, worst clip, rail rate/run and increase
  remains visible even when aggregate screens pass.
- Listening material is reference/A/B with shared timing, equal active RMS and
  one common extra attenuation if needed to avoid clipping. Opaque presentation
  filenames and deterministic seed `0x12345678` randomize labels. The answer key
  is under `listening/private/`; `listening/listener.csv` contains item IDs and
  blank `listener_id,preference,naturalness_a,naturalness_b,intelligibility_a,
  intelligibility_b,artifact_notes` fields. Do not expose the key to listeners.
- Each held-out run prepares twelve balanced items per codec (each speaker's
  utterances 1,2,3); legacy per-clip exceedances are separate diagnostic items.
  At least five actual listeners must each rate all twelve items before a mode
  claim. Preferences are `A`, `B`, or `tie`; ratings are integers 1…5.
  Run `report --run-dir DIR --listening-results FILE` after collection.
  Ten thousand seeded bootstrap replicates independently resample listeners
  and clips, with clip resampling stratified by speaker. Preference CI must be
  strictly positive and intelligibility CI not wholly negative for a pass;
  wholly negative preference/intelligibility means fail, otherwise preference
  is not established. Shared IMBE7100/7200 PCM does not inflate independent
  sample size. This is a small blinded engineering comparison, not an
  ITU-compliant MOS/MUSHRA study. Playback alone is not a judgment.

### DVSI reference vectors

DVSI publishes AMBE-3000 test vectors on
<https://www.dvsinc.com/dlapps/appsoft.shtml>: input speech, DVSI-encoded bit
streams and DVSI-decoded speech for D-STAR, P25 and AMBE+2 (rate 33), plus
channel-error variants. They are the closest available stand-in for real
hardware, and they let decoder output be compared with DVSI's own decoder on
identical bits.

The [site terms](https://www.dvsinc.com/terms_of_use.shtml) grant no license to
copy, modify, transfer or mirror its materials. Local measurement is the tool's
intended use; it does not establish permission. Never commit the vectors,
attach them to issues or pull requests, or use DVSI audio in published
listening material. Only URLs, SHA-256 digests and aggregate numbers belong in
this repository.

```sh
# about 223 MB of the 1.9 GB archive; a directory outside the checkout survives
# rebuilds and can be shared by several checkouts
python3 tools/quality/fetch_dvsi_vectors.py --out ../dvsi-vectors
build/dev-debug/mbe_quality_reframe --codec ambe2400 \
  --from-dvsi ../dvsi-vectors/tv-rc/dstar/dam.bit --out dam.frames
```

- The fetcher needs `curl`. It reads only the needed members of
  `get-usb/tv-rc.zip` through HTTPS range requests. It verifies each CRC-32 and writes
  `dvsi-manifest.json` with SHA-256 digests, retrieval date, source/terms URLs,
  notice hashes and restricted redistribution metadata into the output directory
  (default `../dvsi-vectors`). Reports retain this evidence and the terms notice.
  See [reference-data handling](../legal/reference-data.md) for corpus notices,
  listening-pack checks and recording user-acquired ETSI inputs.
- `tv-rc/*.pcm` are the inputs (s16le, 8 kHz). Each mode directory holds
  `X.bit` (DVSI encoder) and `X.pcm` (DVSI decoder). The fetcher takes the modes
  this library decodes: `dstar`, `p25`, `p25_nofec`, `r33` and `r34`.
  `cmprc.txt` records the rate words: D-STAR uses
  `0x0130 0x0763 0x4000 0x0000 0x0000 0x0048`.
- `tv-rc/<from>/<to>/` (only `alltone`, `dam` and `sine0_4k`) hold DVSI's rate
  conversions between those modes, made with `-rc -rd <from> -re <to>`: DVSI's
  decoded model of the `<from>` bits, re-quantized into `<to>`. They are lossy
  probes of DVSI's model, not an exact readout of it.
- `--from-dvsi` takes 9-byte D-STAR DV data (LSB first), 18-byte P25 frames, or
  9-byte rate-33 frames (MSB first, DMR interleave). It writes Framef rows for
  `mbe_quality_eval`. Bits are only deinterleaved, so the `_eN` channel-error
  variants still exercise FEC and repeat handling. A truncated final frame is
  rejected.
- The `_eN_hd` and `_eN_sd` variants carry the same channel errors, at bit
  error rates of about 1.1, 2.2, 4.9 and 10.6% for `e1`, `e2`, `e5` and `e10`.
  `_sd.bit` (DVSI's `-sdbits 4`) holds one 4-bit value per bit of the `_hd.bit`
  file, two per byte, high nibble first, in that file's MSB-first bit order;
  values 8 to 15 decide 1, so the `_hd` file is exactly their hard decisions.
  `--from-dvsi` recognizes the `_sd.bit` name and writes one lowercase hex digit
  per frame cell; `mbe_quality_eval --soft` decodes those rows through the
  `*SoftFramef` APIs, with each value's distance from 7.5 as the reliability
  (16 to 240).
- Copy the inputs and DVSI's decoded files to `.raw` names before using them as
  `--ref` or `--decoded` operands. The evaluator selects the format by
  extension.

### DVSI scoreboard

`tools/quality/dvsi_scoreboard.py` runs every speech vector of the D-STAR, P25
and rate-33 sets through the tools above and writes `scoreboard.json` and
`scoreboard.md` (needs NumPy). A change is judged by scoring a baseline build
and the candidate, each aligning freely, then both again at their joint lags,
the candidate against the baseline:

```sh
SB="python3 tools/quality/dvsi_scoreboard.py --vectors ../dvsi-vectors --decode-tools build/tools-debug"
Q=build/quality/scoreboard
$SB --tools build/baseline --out $Q/baseline-free
$SB --tools build/tools-debug --out $Q/candidate-free
JOINT=$Q/baseline-free/scoreboard.json,$Q/candidate-free/scoreboard.json
$SB --tools build/baseline --out $Q/baseline --lags-from $JOINT
$SB --tools build/tools-debug --out $Q/candidate --lags-from $JOINT \
  --compare $Q/baseline/scoreboard.json --gates pitch-law
```

- **Speech:** ours vs the input, DVSI's decoder vs the input, and ours vs
  DVSI's decoder, with the `mbe_quality_eval` metrics above (absolute level and
  RMS-matched spectra).
- **Pitch:** a harmonic-peak f0 estimator measures our output, DVSI's output
  and the input on stable voiced frames and reports each against the pitch our
  decoder coded, by f0 band. A decoder that plays sharp or flat shows up as a
  ratio away from 1.
- **Encoder:** our encoding of the input (D-STAR, IMBE and AMBE+2 encoders)
  against DVSI's encoding of it, both decoded by the `--decode-tools` build
  (default `--tools`): pitch error
  against the input, octave errors, per-1 kHz-band voicing agreement, the speech
  metrics against the input, the same metrics of our decoded encoding against
  DVSI's (`encoder.vs_dvsi_bits.*`, with `encoder.level_gap_db` the unsigned
  level difference), rail samples beyond DVSI's encoding and the fraction of
  tone or silence frames. DVSI's encoding is the reference there (its active
  frames and level); DVSI's encoders delay the speech about 200 samples more
  than ours (125 for P25), past the evaluator's -160 lag limit, so ours is
  first delayed by 320 samples (240 for P25). It cannot show how DVSI hardware
  decodes our bits. Each encoding's rows are hashed (`rows_sha256`); with
  `--compare` against a baseline that recorded them, `encoder.rows_changed_vectors`
  counts the vectors whose encoder bits changed. A build whose
  `mbe_quality_encode` does not offer a codec skips that encoder, so this
  script can score a baseline built from an older checkout. Give the baseline
  and the candidate run the same `--decode-tools`, so that only the encoders
  differ (where a change moves both, decode with the candidate's decoder: it
  stands for the decoders that will receive the bits). The run records that
  decoder's identity (`decoder_identity`, a SHA-256 of its `mbe_quality_eval`
  and the library beside it), and `--compare` refuses encoder deltas between
  runs decoded by different or unrecorded decoders.
- **Encoder tones (rate 33 and D-STAR):** our encoding of each tone vector
  against DVSI's bits, frame by frame at the best frame shift (of equally good
  shifts, the one where DVSI's encoder runs a frame behind ours, as on every
  speech vector): the share of
  DVSI's tone frames we send with the same tone index
  (`tone.encoder_agreement`), the share of frames we send as a tone where DVSI
  sends voice (`tone.encoder_extra_rate`) and the mean level difference where
  both agree, in AD steps for rate 33 (`tone.encoder_level_error_ad`) and in
  steps of D-STAR's 8-bit volume (`tone.encoder_level_error_volume`).
  The tones group reports the worst vector of each and, for these, also the
  mean over vectors (`_mean`), since DVSI itself is inconsistent on its
  deliberately malformed tone vectors (`cpvbad`, `dtmfvbad`).
- Every comparison must find an alignment inside the evaluator's lag search.
  A decoder comparison without one (a silent decode, or a lag at a limit)
  fails the run. An encoding without one (one that is mostly silence frames,
  say) is left unscored, counted in `encoder.alignment_failures`, and fails
  any gate on its partition.
- **Channel errors** (`--partitions channel`): `dam`'s `_eN_hd` and `_eN_sd`
  vectors, decoded as hard and as soft decisions, against DVSI's decode of the
  same bits. A loud frame (RMS 300 or more in DVSI's decode of the clean bits)
  is muted by DVSI when its decode of the errored bits is 30 dB or more below
  that; ours is muted when the frame carries `MBE_PROCESS_FLAG_MUTE`. The
  section reports, per vector, how often our mute decisions agree with DVSI's
  on loud frames, the false and missed mute rates, the first loud frame each
  mutes, our repeat rate, and the log-spectral distance from DVSI's output on
  the loud frames with no muted frame at or next to them (a frame's spectra
  reach into its neighbours), with our output aligned and level-matched by the
  clean bits' comparison. Each run mutes its own frames, so `--compare`
  measures that distance for both runs on the frames that count in both
  (`channel.hd.lsd_common_db`) and reports their share of all the loud frames
  DVSI plays (`coverage`), so frames a run did not measure count as lost. Soft
  vectors need a baseline build whose tools read them (`--from-dvsi` of an
  `_sd.bit` file and `mbe_quality_eval --soft`): for an older library, build it
  with this checkout's `tools/quality` sources.
- **Partitions** are by source recording. Several vectors are one recording
  through different filters, noise or levels (found by matching 3 s blocks of
  the inputs at any lag): `dam` and its variants (`lvl` among them),
  `fambm22a`, `tambf22a`, `ucarm15a`, `ucarf15a` and `ucarm20a`, `tambf22b` and
  `tambf32b`, and `tia11` and `p01mirs`. Each such group sits in one partition,
  and the bootstrap intervals resample whole groups. `development` and
  `validation` were both examined during earlier changes; `final` (`irstia`
  and `mark`) is scored once per change set (`--partitions final`). Not scored
  as speech: `xfer`, the DTX and channel-error variants and the synthetic
  pulse trains and sweeps. The `tones` partition covers the DTMF, single,
  KNOX, alert and call-progress vectors.
- **Tones** are compared with one acoustic detector run on both outputs, on
  windows where DVSI's output is steady: level, frequencies, balance, how
  often both agree on the tone playing, and the fraction of windows playing a
  different tone. Single tones must agree within 10 Hz (they are 31.25 Hz
  apart); each dual-tone component within 25 Hz, since DVSI plays dual tones
  up to 20 Hz off nominal and neighbouring DTMF frequencies are 73 Hz or more
  apart.
- `--compare` adds paired per-vector deltas with bootstrap 95% intervals and the
  worst vector, and counts the baseline's speech vectors whose PCM changed
  (a vector missing from the run counts as changed). Both runs are measured at
  one lag per comparison, so a change that moves the evaluator's alignment
  search does not also move the samples measured: `--lags-from` takes the lags
  a scoreboard recorded (`lags` in each vector's entry, with the frame shift of
  the encoder-tone comparison), or, given two, their joint lags. Where the two
  free runs' searches disagree, the joint lag is that of the run whose
  alignment reached the higher correlation (`alignment`), since a change can
  let a decoder align where the old one could not: main's D-STAR gain law
  played onsets up to 16 dB low and pulled its alignment with DVSI's output to
  40-96 samples, where the corrected decoder aligns at -16 to 0 (P25 and rate
  33 at -16 to +16 in both builds). Without `--lags-from`, `--compare` uses
  the baseline's lags. A vector
  whose recorded lags still differ from the baseline's is not paired.
  `--gates` evaluates gate sets (a comma list) from
  `tools/quality/dvsi_gates.json`, declared before the change it judges, and
  exits non-zero when a gate fails. A gate set fails unless the run was scored
  at fixed lags and both runs scored every vector of its partition in every
  mode it gates, and a paired check fails when a metric was measured for a
  vector in only one of the two runs.
  `all_no_increase` and `all_no_decrease` hold every vector to the limit, not
  just the mean, and a check's `min_coverage` fails it when its common frames
  are too few.

`mbe_quality_eval --params frames.jsonl` writes one JSON record per frame:
`bits` (the frame's parameter bits after FEC), `decoded` (the frame's own
model from the public `mbe_decode*Parms` decoder,
with its status), `history` (the prediction state after the frame) and `synth`
(the enhanced model synthesized). Each holds `w0` (radians per sample), `L` and
per-band `Vl`, `Ml` (linear amplitude as synthesized; for AMBE unvoiced bands
it includes the unvoiced scale factor) and `log2Ml` (the codec's log2 spectral
amplitude before that factor), plus `gamma`. `synth` also holds `PHIl` and
`PSIl`, the model's harmonic phases for the frame. Synthesis renders each voiced
harmonic with `PHIl` plus the output high-pass's phase at its frequency
(`mbe_low_band_phase()`: 180 degrees toward DC, 90 at 100 Hz, under 9 above
1 kHz), so phase readbacks can be checked against the decoder itself. Dumping
does not change the PCM.

### Encoder evaluation

`mbe_quality_encode` (built with `MBELIB_BUILD_TOOLS=ON`) encodes 8 kHz s16le
speech with the D-STAR (`--codec ambe2400`), AMBE+2 (`ambe2450`) or IMBE
(`imbe7200`) encoder and writes the 49- or 88-bit rows `mbe_quality_eval`
decodes, so encoder output is scored against the input with the same metrics
as decoder output. Comparing with DVSI's own encoding of the same input:

```sh
build/dev-debug/mbe_quality_encode --codec ambe2400 --in dam.raw --out dam.rows
build/dev-debug/mbe_quality_eval --codec ambe2400 --frames dam.rows --out ours.wav --ref dam.raw --json ours.json
build/dev-debug/mbe_quality_reframe --codec ambe2400 --from-dvsi ../dvsi-vectors/tv-rc/dstar/dam.bit --out dvsi.frames
build/dev-debug/mbe_quality_eval --codec ambe2400 --frames dvsi.frames --out dvsi.wav --ref dam.raw --json dvsi.json
```

Both bit streams then go through the same decoder, which itself matches DVSI's
output on DVSI's bits, so the band, LSD and envelope figures compare the
encoders. One zero flush frame (`--flush-frames`, default 1) covers the
encoders' analysis delay.

### Attribution and limitations

Shared phase regeneration and corrected WOLA were already present at checkpoint
`8da251d`; interpolation, enhancement, adaptive smoothing and most excitation
state handling predate this quality-completion work. The additional production
corrections are limited to clamping the **applied** adaptive attenuation to
nonnegative values while preserving signed threshold recurrence, and dispatching
accepted AMBE2400 tone IDs 5/6 consistently with IDs 7…122 while setting `TONE`.
Tone mapping, NOTONES behavior, limiter, codec concealment policy, public headers,
state layout and exports are unchanged. Regressions cover positive attenuation,
fade/recovery, tones and rejection thresholds, warm replay/interleaving, all
160 WOLA reconstruction samples, all 256 IMBE fundamentals, and corrected ECC.

Historical three-mode figures under the old tooling were absolute **signed-mean**
crest errors 1.225/1.941/2.086 → 0.301/0.205/0.317 dB, not mean per-clip absolute
errors. That 27-second corpus included the duplicated `hts1a`; historic tooling
hashes do not bind every report to the current evaluator/calibrator. Preserve
those artifacts as historical input, not regenerated four-mode evidence.
Spectral/envelope improvements alone never establish listener preference.
Without the specified actual human judgments, perceptual acceptance remains
`not_established`; an objective failure is a completed negative assessment,
not grounds for claiming “all modes improved.”

The completed fixed-profile assessment is preserved under
`build/quality/runs/{original-to-final-debug,checkpoint-to-final-debug,heldout-debug,heldout-minus12db-debug,heldout-release,heldout-minus12db-release}/`.
Each comparison contains its own immutable manifest and generated acceptance,
analysis, tables, figures, and listening material. `build/quality/evidence-index.json`
indexes the delivered runs, controls, compatibility proof, and benchmark results.

On the 40 held-out references (118.032 seconds), Debug and matched Release agree
on the following aggregate decisions; each level uses the same frozen gains:

| Actual frame mode | Correctness | Native objective | −12 dB objective | Perceptual |
| --- | --- | --- | --- | --- |
| IMBE7200 | pass | fail | fail | not established |
| IMBE7100 | pass | fail | fail | not established; shared IMBE PCM |
| AMBE2450 | pass | fail | fail | not established |
| AMBE2400 | pass | pass | pass | not established |

Debug IMBE mean per-clip absolute crest error increases 1.072→1.145 dB at native
level and 1.007→1.023 dB at −12 dB. AMBE2450 mean join excess increases
0.332→0.564 dB and 0.043→0.303 dB respectively, exceeding its +0.1 dB allowance.
AMBE2400 passes the aggregate screens but retains 17 native and 15 attenuated
per-clip exceptions. Mean STOI/ESTOI decrease slightly in all three independent
speech modes; they do not establish perceptual harm or preference by themselves.
There are no candidate rail-count increases in these held-out runs.
The two additional production fixes preserve all legacy clean checkpoint PCM;
the held-out screen failures concern the pre-existing branch synthesis changes,
not evidence for retuning either correction. The result does **not** establish
“all modes improved.”

After the automated assessment, the maintainer completed all 48 native-level
blinded presentations as an informal single-listener preference session:
14 baseline preferences, 14 candidate preferences, and 20 ties. These are
presentation counts, not 48 independent pairs: the two IMBE modes repeat the
same PCM comparisons. Only the first two items received full numerical ratings;
the remaining scores were left blank. The strong second-played preference
(A: 2, B: 26, ties: 20) is a playback-order concern, not proof of its cause.
The formal multi-listener perceptual result remains `not_established`.
The maintainer explicitly chose to retain the candidate despite the mixed
objective and personal listening results; that decision does not alter the
recorded measurement gates.
