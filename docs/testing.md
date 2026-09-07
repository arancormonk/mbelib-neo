# Testing Policy

mbelib-neo uses automated tests, static analysis, sanitizers, and fuzzing to
reduce regression and security risk.

## Test Suites

The CTest suite includes:

- API/version/result helper checks
- ECC tests for hard and soft Golay/Hamming paths
- noise determinism and frame-state determinism checks
- parameter and synthesis behavior checks
- float-to-int16 conversion parity checks
- golden PCM hash regression checks
- SIMD architecture detection checks

Run the default test suite with:

```sh
cmake --preset dev-debug
cmake --build --preset dev-debug -j
ctest --preset dev-debug --output-on-failure
```

## Continuous Integration

GitHub Actions runs tests and quality checks on pull requests and pushes to the
primary branch. Required checks include cross-platform builds, sanitizer tests,
static analysis, repository security guardrails, workflow linting, dependency
review, secret scanning, OSV scanning, and install/consume checks.

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
AddressSanitizer, including fixed-size hard/soft frame decode paths and
parameter synthesis paths.

## Speech Quality Evaluation

Developer-only tooling; not a CTest registration or a CI dependency. The C99
`mbe_quality_eval` executable uses only the installed public API and an independent
256-point FFT. The Linux A/B scripts additionally require Bash, Python 3, git,
gcc/g++ (C++17), `patch`, coreutils, and network access for pinned source fetches.
An audio player such as `aplay` is needed for listening.

```sh
cmake --preset dev-debug -DMBELIB_BUILD_TOOLS=ON
cmake --build --preset dev-debug --parallel
tools/quality/fetch_speech_corpus.sh
tools/quality/build_op25_encoder.sh

# Use an external worktree to avoid recursively copying build outputs.
git worktree add /tmp/mbelib-quality-baseline 5fd3f3dd738b76bc9a8446ff4d7ae99e565e35e0
cmake -S /tmp/mbelib-quality-baseline -B build/quality/baseline \
  -DCMAKE_BUILD_TYPE=Debug -DMBELIB_BUILD_TESTS=OFF -DMBELIB_BUILD_EXAMPLES=OFF
cmake --build build/quality/baseline --parallel

# Harness identity check: both columns must have identical PCM hashes.
tools/quality/run_quality_ab.sh build/quality/baseline build/quality/baseline \
  | tee build/quality/baseline.txt
# Candidate comparison: same executable, switched libmbe-neo.so.2 search path.
tools/quality/run_quality_ab.sh build/quality/baseline build/dev-debug \
  | tee build/quality/candidate.txt
```

Use separate baseline/candidate library directories with the same build flags.
The runner writes per-file WAVs and flat JSON metrics to `build/quality/out/`,
prints baseline/candidate columns and unweighted per-mode means, and caches
encoded bits under `build/quality/frames/`. A third positional argument selects
another directory containing 8 kHz mono s16le `*.raw` references.

### Encoder provenance and calibration

- Codec2 corpus commit: `310777b1c6f1af0bc7c72f5b32f80f6fd9136962`
  from <https://github.com/drowe67/codec2>. The five checked-size references are
  `hts1a`, `hts2a`, `hts1`, `kristoff`, and `ve9qrp_10s`.
- OP25 encoder commit: `71abcd0ead32f86f51615ea6cc8a6a4dba4c949a`
  from <https://github.com/boatbod/op25>. Its GPL-3.0-or-later encoder and bundled
  mbelib are compiled into a separate process, never linked into mbelib-neo.
  Their incompatible `mbe_parms` layouts must not be mixed.
- `op25_mbelib_gain.patch` adapts disposable compilation copies, leaving the
  fetched checkout unchanged. D-STAR gain quantization predicts from previous
  gamma, and b8 quantizes/transmits the even codebook indices reconstructed by
  this decoder. These are **decoder-targeted quality fixtures, not an
  interoperable D-STAR transmitter or a conformance corpus**.
- Passing OP25's unadjusted gain through mbelib-neo's historical PCM conversion
  produces substantial clipping. `calibrate_encoder.py` adjusts the encoder's
  log2 gain quantizer, leaving input speech analysis unscaled. It targets
  baseline decoded active RMS at −3 dB relative to reference (±0.5 dB tolerance),
  with up to six calibration passes. Alignment must be in 0–640 samples.
  The candidate never participates in calibration.
- Calibration manifests bind cached frames to hashes of the reference, encoder,
  evaluator, baseline library, and calibration implementation. Both library
  variants decode exactly the same cached bits. Changing any calibration input
  invalidates its cache; failed calibration is an error, not a fallback.

IMBE 7200, AMBE 2450, and AMBE 2400 have encoder modes. No open IMBE 7100 encoder
is used: that decoder converts parameters into the shared IMBE 4400 synthesis
path, covered here by IMBE 7200, while existing frame-path tests cover conversion.

### Evaluator interface and measurements

```sh
build/quality/op25_encode --mode imbe7200 \
  --in build/quality/corpus/hts1a.raw --out build/quality/example.frames \
  --gain-adjust 3.4
build/dev-debug/mbe_quality_eval --codec imbe7200 \
  --frames build/quality/example.frames --out build/quality/example.wav \
  --ref build/quality/corpus/hts1a.raw --seed 0x12345678 \
  --json build/quality/example.json
```

`--gain-adjust` is log2 attenuation; the encoder's default is `log2(7)`.
For comparisons, use the runner's baseline calibration rather than guessing it.
Raw encoder inputs have an even byte count; the last partial frame is zero-padded.
Frame files contain one row-major literal binary string per line:
88 bits for IMBE 4400 Dataf, 49 bits for AMBE 2400/2450 Dataf, or 96 bits for the
AMBE 3600x2400 Framef API (including zero rectangular padding).
The evaluator rejects invalid digits, lengths, or codec combinations with exit 2;
a negative process status exits 3 with the frame index. It requires at least
256 aligned samples to report spectral metrics.

References may be headerless `.raw` s16le or RIFF/WAVE PCM `.wav`, mono, 8000 Hz,
16-bit. WAV chunks are walked rather than assuming a 44-byte input header.
Output WAVs are 16-bit PCM with a 44-byte header. `pcm_fnv1a` hashes the
little-endian output sample bytes, before measurement gain normalization.

- Alignment maximizes correlation of log-RMS envelopes: centered 160-sample
  windows, 8-sample hop, lag search −160…800 samples. Reference and decoded are
  trimmed to common overlap; active 160-sample reference frames are within 40 dB
  of maximum energy. Positive lag means decoded speech is delayed.
- `level_offset_db` is decoded/reference active RMS before normalization.
  Reference-based measurements then use RMS-matched signals in int16 units.
- `lsd_db`: mean active-frame log-spectral distance, 256-point symmetric Hann,
  hop 80, bins 2…118, common power floor `1e-6 * max(reference power)`.
- `env_corr`: correlation of log 40-sample RMS envelopes, floored 60 dB below
  each signal's maximum.
- `crest_ref_db`, `crest_dec_db`, `crest_delta_db`: mean 20 ms peak/RMS values
  and decoded-minus-reference difference.
- `band_delta_db_*`: active STFT power ratios in 0–500, 500–1000, 1000–2000,
  2000–3000, and 3000–4000 Hz bands.
- `boundary_index_db`: unshifted decoded derivative energy near each frame
  boundary divided by that near the frame center. Keep as a diagnostic only:
  a perfectly continuous 75 Hz sine yields +16.61 dB, so the index is not a
  standalone click detector and lower is not necessarily better.
- `join_dec_db` and `join_ref_db`: at the actual synthesis joins retained after
  alignment, compare summed squared across-join derivatives with the mean
  squared derivatives at the eight neighboring positions (four on each side).
  Both use the same aligned active-reference mask.
  `join_excess_db = join_dec_db - join_ref_db`; positive values indicate stronger
  localized join derivatives than the reference. Inspect per-file values too:
  this statistic is not a standardized perceptual measure.

Without `--ref`, crest, absolute mean STFT band powers (`band_db_*`, unnormalized
int16 units), the original boundary index, and decoded join contrast are reported.
Silent inputs have no active frames; empty-aggregate metrics are zero except
absolute band levels, which use the numerical power floor.

### Quality acceptance and interpretation

For the five-file corpus, compare unweighted per-mode means:

- Absolute `crest_delta_db` decreases.
- `join_excess_db(candidate) <= max(0, join_excess_db(baseline)) + 0.1 dB`.
- `lsd_db` changes by at most 0.3 dB; `env_corr` decreases by at most 0.01.
- Every mean band delta changes by at most 0.5 dB.
- Matched Release `bench_synth` cost increases by at most 10%.

The reference-matched join gate replaces the original monotonic boundary-index
gate; both metrics remain visible. Investigation controls showed zero join excess
for identity/gain scaling and increasing excess for injected frame-boundary steps
in all five recordings. Per-file join excess can still be positive when a mode's
mean is negative; neither an aggregate pass nor successful playback proves
subjective quality.

The regenerated-phase evaluation at λ=0.44, D=19, γ=0.72 passed these gates.
Across modes, absolute mean crest error changed from 1.225/1.941/2.086 dB to
0.301/0.205/0.317 dB (IMBE/AMBE2450/AMBE2400), envelope correlation increased,
LSD changed by less than 0.015 dB, and mean band deltas changed by less than
0.079 dB. The original boundary-index increases were 1.780/1.257/0.360 dB:
they are not hidden by the revised gate. Seven alternating matched Release
benchmark runs gave a 1.04% median-mean cost increase on the development x86-64
machine. These are corpus/build-specific results, not universal quality claims.

Regression coverage also checks RNG-independent fully voiced synthesis after
the third noise-buffer advance, flat/rising/falling spectral phase behavior,
and AMBE prediction at the 56-harmonic interpolation endpoint. The endpoint
regression fails under UBSan before the bounds fix; the real 15-file/mode corpus
passes ASan/UBSan afterward. Public header bytes and parameter layout are unchanged.

Listen to the reference and both variants before a perceptual acceptance:

```sh
aplay -f S16_LE -r 8000 -c 1 build/quality/corpus/hts1a.raw
aplay build/quality/out/hts1a.imbe7200.baseline.wav
aplay build/quality/out/hts1a.imbe7200.candidate.wav
```

Repeat for both AMBE modes and the other speakers. Audio playback was exercised
during development; a human listening preference was not established.
