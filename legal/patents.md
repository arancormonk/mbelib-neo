<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Patent research record

Reviewed 2026-10-17. These are engineering observations, not claim construction,
an exhaustive search, or a legal clearance opinion. Copyright licensing does
not resolve patent rights. AMBE+2 and TETRA remain available with open questions.
No educational or local-measurement exemption is asserted.

## US 8,359,197

[The published patent and database record](https://patents.google.com/patent/US8359197B2/en)
show “Active” and a projected date of 2028-05-20. Database status is not a
verified legal conclusion. Check the official register, maintenance history,
term adjustments, family members and relevant jurisdictions before relying on it.

| Claims | Feature to review | Implementation observation |
| --- | --- | --- |
| 1; 2–41 | First FEC word carries a subset of pitch, voicing and gain quantizer bits | `ambe3600x2450_enc.c` quantizes these parameters; `ambe_encoder.c` protects selected bits in C0 and packs the remainder elsewhere. |
| 42; 43–59 | Decode that partial parameter word, reconstruct parameters and speech | `ambe_common.c` corrects C0; `ambe3600x2450.c` reconstructs pitch, voicing, gain and spectral information for shared synthesis. |
| 60; 61–71 | Tone/voice selection and reconstruction | `ambe3600x2450.c` tests the tone marker and reconstructs tone or voice. |
| 72; 73–87 | Codebook index dependent on frame length | Public codec entry points select fixed layouts; whether that satisfies the claimed dependency is unresolved. |

This feature map flags review targets. It does not assert that any complete
claim is met, valid, enforceable or infringed. No patent license was established.

## Separate watchlist

| Record | Reason for review | Status of project assessment |
| --- | --- | --- |
| [US 8,595,002](https://patents.google.com/patent/US8595002B2/en) | Related half-rate family; parameter/FEC coding | Open; independently verify family and legal status. |
| [US 5,701,390](https://patents.google.com/patent/US5701390A/en) | Regenerated phase equations used in synthesis | Historical algorithm reference; no clearance inferred from age. |
| [US 2025/0336405](https://patents.google.com/patent/US20250336405A1/en) | Voicing smoothing | Open; compare claimed method with speech analysis and adaptive smoothing. |
| ETSI TETRA declarations ISLD-200101-005 and ISLD-190001-073 | Full-rate ACELP, including codebook and gain processing | Prior search is recorded in [tetra-acelp.md](../docs/tetra-acelp.md#ipr-search); completeness and legal conclusions remain open. |

The watchlist is not exhaustive. Revisit these records when synthesis, tone
detection, quantization, FEC or TETRA processing changes and before each release.
Keep research dates, implementation revision and evidence separate from legal
conclusions. An anticipated expiration date passing never automatically changes
an issue to resolved. The baseline revision is recorded in
[implementation-provenance.md](implementation-provenance.md).
