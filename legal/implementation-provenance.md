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
project portability macros, pinned PFFFT, and mathematical/state contracts
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
| Regenerated phase and calibration | Retained project implementation using US 5,701,390 equations 7–9, with independently measured low-band, shelf and output gains. Patent status is separate. |
| Tone frequencies and parameter decoding | Retained specification tables, ISC-origin code and project measurements; no public signature or structure changes. |

The specification's equations support algorithms; they do not establish the
provenance of every historical implementation. The identified implementation
portions were replaced, rather than having only comments or names changed.

Library policies preserved for existing streams include a zero first unvoiced
buffer, the 96-sample state representation, a 48-bit comfort-noise recurrence,
the pitch-index-134 initial model, repeat overflow reset, 0.003 comfort gain
for applicable codecs, and a hard limiter at 95% of full scale. These are not
claimed as TIA requirements. The 48-bit recurrence is also documented in
[Java Random's API](https://docs.oracle.com/en/java/javase/21/docs/api/java.base/java/util/Random.html).
P25 mute-noise levels remain those of TIA-102.BABA §7.8 and BABA-1 §5.7.

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
