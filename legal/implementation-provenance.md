<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Implementation provenance review

The 2.6.0 review began from commit
`8cc203314419ae38377bfc469041c3b2d534d9c5` on 2026-10-17. Historical source,
integration changes and exposure records were examined by the integrating
Codex agent and retained in a local sibling audit directory. They are excluded
from release artifacts. This engineering procedure is not a legal certification
or a finding that historical releases were cleared.

An independent Codex implementation agent received no inherited conversation
or prior implementation source. Its separate workspace contained the
TIA-102.BABA and BABA-1 specifications, public/private interface declarations,
project portability macros, pinned pffft-neo, and mathematical/state contracts
prepared by the reviewer. The material hashes, authorship record and
equation-derived fixtures are retained locally. No prior vocoder implementation
was supplied. The integrating agent reviewed the new implementations separately.

| Scope | Basis and treatment |
| --- | --- |
| Adaptive smoothing | New implementation of TIA-102.BABA §§8–9, equations 111–116. |
| Unvoiced synthesis | New implementation of §11.2, equations 117–126: LCG, windowed FFT, band normalization, inverse FFT and weighted overlap-add. |
| Noise state, initialization, limiting | Independently implemented library contracts described below. |
| Tone oscillator and redundancy check | Independent implementation of oscillator/state contracts and the local check over BABA-1 §7.2 Table 10 fields. |
| Frame repetition | TIA-102.BABA §7.7 equations 97–104 and BABA-1 §§5.6–5.7; new default-model reset. Existing per-frame accounting and generator continuity remain. |
| Spectral enhancement and voiced synthesis | Retained specification-supported equations 105–110 and 134–138, compatible ISC-origin code and subsequent project optimizations. |
| Regenerated phase and calibration | Retained project implementation using US 5,701,390 equations 7–9 (on the codec paths over each frame's decoded log2 magnitudes, in the public `mbe_synthesizeSpeechf()` over the enhanced ones), with independently measured low-band, shelf and output gains. Patent status is separate. |
| Tone frequencies and parameter decoding | Retained specification tables, ISC-origin code and project measurements; no public signature or structure changes. Since 2.6.14 the rate-33 and D-STAR decoders play each dual tone at TIA-102.BABA-1 Annex J's harmonics l1 and l2 of a fundamental on DVSI's grid of 2^-19 cycle per sample (`src/internal/mbe_tone.h`): 34 of the 36 dual-tone fundamentals are Annex J's printed f0 rounded to that grid, and the two that differ (indices 145 and 147) were measured on DVSI's rate-33 decodes of its tone vectors (D-STAR's 20 dual tones give the same values). Single tones follow a rule of their own (harmonic ID / 16 + 1, not Annex J's harmonic number), measured on the same decodes. Table 9's nominal frequencies stay for the encoders, tone detection and index validity. |
| P25 tones sent as voice frames | Since 2.6.14 (`src/imbe/imbe7200x4400.c`): the tone table's harmonic numbers are Annex J's; its pitch codes are `floor(4 pi / w0 - 39)` of the fundamental Annex J gives, taken for a dual tone as the printed f0 (with the two fundamentals measured on DVSI's decodes for indices 145 and 147) and for a single tone as 31.25 ID / l, the quotient Annex J prints rounded (built from the printed values, nine codes would differ: IDs 32, 65, 66, 68, 69, 71, 73, 74 and 76); so built they equal the codes in DVSI's bit streams of its tone vectors for all 154 tones; its voicing words were read from those bit streams; the two dominance limits, the level constant and the 2^-19 pitch grid were measured on DVSI's decodes of its tone and speech vectors. The test follows the description of US 7,970,606 ([patents.md](patents.md)). |

The specification's equations support algorithms; they do not establish the
provenance of every historical implementation. The identified implementation
portions were replaced, rather than having only comments or names changed.

Library policies preserved for existing streams include a zero first unvoiced
buffer, the 96-sample state representation, a 48-bit comfort-noise recurrence,
the pitch-index-134 initial model, repeat overflow reset, 0.003 comfort gain
for applicable codecs, and a hard limiter at 95% of full scale. These are not
claimed as TIA requirements. The 48-bit recurrence is also documented in
[Java Random's API](https://docs.oracle.com/en/java/javase/21/docs/api/java.base/java/util/Random.html).
The P25 IMBE mute-noise level remains that of TIA-102.BABA §7.8. Since 2.6.9 the
AMBE 3600x2450 path mutes at the level of DVSI's decoder, as D-STAR does, in
place of BABA-1 §5.7's ([DVSI parity](../docs/dvsi-parity.md#rate-33-mutes)).

The integration review retained the independently written algorithms and
adapted includes, index types and private helper names to repository conventions.
The shared FFT workspace reduces the encoder's allocation sites from thirteen
to seven; failure injection covers every remaining site and verifies cleanup,
retry and allocation-free frame processing. Synthetic equation tests compare
the transform with a direct DFT, autocorrelation with a direct sum, smoothing
with double-precision equations, and noise with its documented recurrences.

## Regression evidence

Compared with the baseline above, all 237 DVSI decoder streams in the
development, validation, final and tone partitions produce byte-identical PCM
in scalar Debug. All 185 encoder streams are bit-identical: 79 D-STAR, 79
AMBE+2 and 27 IMBE, including the applicable tone vectors. All four existing
DVSI parity gate sets pass without changing their thresholds.

The independent unvoiced band-power reduction uses one scalar accumulation
order across targets instead of vector horizontal reductions. This preserves
the specification's normalization and changes SIMD rounding. SSE2 Debug
comparison across the same 237 streams found 165 samples differing by one
int16 unit among 41,035,200 samples (RMS difference 0.00201 int16 units).
Scalar golden hashes and all int16 golden hashes remain unchanged. Four
AMBE float golden hashes were regenerated only after this comparison; SSE2
and AVX2 values agree on the host and Ubuntu 24.04's glibc 2.39.

Local validation passed Debug, optimized GCC and Clang (including x86-64-v3
with fast-math/LTO), Clang ASan/UBSan, NOTONES, Ubuntu 24.04 Debug, and ARM32
default-FPU and NEON under QEMU. Public declarations and structure definitions are unchanged;
all 112 exported symbol names match the baseline and `SOVERSION` remains 2.
TETRA's 5,340 TE1 frames match TE2 PCM byte for byte. The complete-source
rehearsal builds, tests, installs and compiles both shared and static consumers
in a fresh directory with networking disabled. Source and binary SBOMs validate
against the SPDX 2.3.1 JSON schema. Broad local quality and supply-chain
guardrails passed. Native macOS/Windows execution and signed release
attestations remain the release workflow's checks; no release was published.

Historical disclosure and original-source backfill drafts for the eleven
affected published tags remain local. They are not clearance findings or
verified replacements for the original release artifacts. TETRA data and
patent questions remain governed by [status.md](status.md) and
[patents.md](patents.md), independently of these engineering results.
