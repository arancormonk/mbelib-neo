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

Developer-only tooling; no new library runtime or default-CI dependency. The C99
`mbe_quality_eval` uses the installed public API and an independent 256-point FFT.
`mbe_quality_reframe` is a statically linked private fixture generator. Linux A/B
runs additionally require Bash, Python 3, git, gcc/g++ (C++17), `patch`, and
coreutils. Fetching held-out speech requires network access and `ffmpeg`.
Offline reports require NumPy, SciPy, matplotlib, and pystoi; install these in a
local virtual environment rather than making them library dependencies.

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
  hops. Positive lag means decoded speech is delayed. The runner measures the
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
  and 3000–4000 Hz. Without reference, decoded absolute `band_db_*` remains
  available where supported.
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
  `pcm_max_rail_run` are measured before normalization. Tone/erasure/repeat/mute
  frame counts and C0/protected/C4/total error sums retain process-result context.

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
