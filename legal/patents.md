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
| [US 5,701,390](https://patents.google.com/patent/US5701390A/en) | Regenerated phase equations used in synthesis: the codec paths apply its edge-detection kernel to each frame's decoded spectral magnitudes (`log2Ml`), as its text describes the computation; the public `mbe_synthesizeSpeechf()` applies it to the enhanced, smoothed magnitudes it renders | Historical algorithm reference; no clearance inferred from age. |
| [US 8,200,497](https://patents.google.com/patent/US8200497B2/en) | Synthesis in 10 ms subframes (square-root triangular windows; an impulse sequence convolved with a pitch-dependent impulse response), which the codec paths render. Since 2.6.13 `imbe7200x4400.c` also applies the droop of its K = 256 synthesis table, read with linear interpolation (`sinc^2(pi l / 256)` at harmonic l), to the harmonics of a P25 tone model. The pulse train's octave step (`mbe_subframe.c`) was measured on DVSI's decodes and is not taken from this or another document. Read again for the mid-frame voicing change: its eqs 6-7 give the pulse phase per subframe from the counts of voiced harmonics of the two models that bound it, and reinitialise the pulse positions when both have none, which is the model the library already had; it says nothing about which bands a mid-frame model voices. The P25 onset columns and the rate-33 rows' mid-frame voicing (`mbe_subframe.c`) were measured on DVSI's decodes and are not taken from this or another document | Algorithm reference. The database record showed “Expired - Fee Related” on 2026-11-09; not verified against the official register, and no clearance inferred. |
| [US 7,970,606](https://patents.google.com/patent/US7970606B2/en) | Decoder-side tone identification: its description checks a frame's reconstructed parameters against its tone table (TIA-102.BABA-1's Annex J) by fundamental, voicing and the dominance of the tone's harmonics over the other spectral magnitudes, then attenuates the others. Since 2.6.13 `imbe7200x4400.c` does this for a P25 tone sent as voice frames (`imbe_tone_word()`, with limits measured on DVSI's decodes). Read again for the mid-frame voicing change: its encoder decides voicing in eight 500 Hz channels; the decoder's 500 Hz onset columns were measured on DVSI's decodes, not taken from it | Algorithm reference. The database record showed “Expired - Lifetime”, adjusted expiration 2025-09-08, on 2026-11-09; not verified against the official register, and no clearance inferred. |
| [US 6,912,495](https://patents.google.com/patent/US6912495B2/en) | Read for the mid-frame voicing change: it lists band edges at 375, 875, ... 3375 Hz. In DVSI's P25 bit streams of the development and validation vectors, the centres of the bands at voicing transitions (the first voiced band of a run above the first band and the first unvoiced band after a run, the last band of a frame left out, as it can be longer) cluster just above 350 to 375 Hz modulo 500 (none at 225 to 275 Hz for f0 100 to 140 Hz), which is consistent with an encoder that maps such channel decisions onto the bands. The library's P25 onset columns, whose edges lie at 505 + 500 k Hz, were measured on DVSI's decodes and are not taken from this or another document | Algorithm reference. The database record showed it expired on 2023-11-28 (the date of that reading is not recorded); not verified against the official register, and no clearance inferred. |
| [US 2025/0336405](https://patents.google.com/patent/US20250336405A1/en) | Voicing smoothing | Open; compare claimed method with speech analysis and adaptive smoothing. |
| ETSI TETRA declarations ISLD-200101-005 and ISLD-190001-073 | Full-rate ACELP, including codebook and gain processing | Prior search is recorded in [tetra-acelp.md](../docs/tetra-acelp.md#ipr-search); completeness and legal conclusions remain open. |

The watchlist is not exhaustive. Revisit these records when synthesis, tone
detection, quantization, FEC or TETRA processing changes and before each release.
Keep research dates, implementation revision and evidence separate from legal
conclusions. An anticipated expiration date passing never automatically changes
an issue to resolved. The baseline revision is recorded in
[implementation-provenance.md](implementation-provenance.md).
