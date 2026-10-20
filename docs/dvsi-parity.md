<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Parity with DVSI's AMBE-3000

DVSI's AMBE-3000 vocoder chip sits in most D-STAR, P25 Phase 2, DMR, NXDN and dPMR radios, so
its decoders set what these codecs sound like and its encoders set the bits other radios receive.
DVSI publishes test vectors for the chip ([DVSI reference vectors](testing.md#dvsi-reference-vectors)):
input speech, DVSI's encoding of it and DVSI's decoding of those bits, for D-STAR, P25 (IMBE) and
AMBE+2 rate 33, with channel-error and rate-conversion variants. This page records what was
changed to bring this library's decoders and encoders closer to them, where DVSI departs from the
TIA-102 text, and what still differs.

The rule followed: the specification first (TIA-102.BABA, and TIA-102.BABA-1 for rate 33);
where it is silent and DVSI's output and the input speech disagree, DVSI. Where DVSI departs from
the specification text, the choice made is recorded in [Departures from the
specification](#departures-from-the-specification).

The vectors are measured locally and never committed or redistributed; only aggregate numbers
appear here.

## Method

- **Same bits, two decoders.** Our decode of DVSI's bits is compared with DVSI's decode of the
  same bits, so any difference is the decoder's.
- **Readback of DVSI's model.** On steady voiced frames, DVSI's harmonic amplitudes and phases
  are read from its output by fitting the harmonics of the coded fundamental at frame
  boundaries. The readback was validated on this library's own output first: it recovers our
  rendered phases within a median 0.12-0.23 degrees on pulse trains (0.9-1.6 on speech), and our
  amplitudes with a ratio of 1.000 (interquartile range +-0.4% on pulse trains, +-3% on speech).
- **Rate conversions as probes.** DVSI's rate-converted streams re-quantize DVSI's decoded model
  of one mode into another, so they show in the target mode's terms (pitch code, gain index,
  vector indices) what DVSI's decoder made of the source bits.
- **Channel errors.** `dam`'s bits at 1.1, 2.2, 4.9 and 10.6% bit errors, as hard decisions and
  as 4-bit soft decisions, show when DVSI repeats and mutes frames.
- **Scoreboard.** `tools/quality/dvsi_scoreboard.py` measures every change on all of it
  ([DVSI scoreboard](testing.md#dvsi-scoreboard)). Its development and validation partitions
  were examined while these changes were made, so their results are retrospective; the final
  partition (`irstia`, `mark`) was scored only after the changes were fixed (twice; see [Scoring
  history](#scoring-history)).

## Changes

### D-STAR gain

D-STAR's gain is absolute, `gamma = 2 Dg[b2] - 0.09`, with no term from the previous frame.
mbelib carried over rate 33's differential law (`Dg[b2] + 0.5 gamma(-1)`, TIA-102.BABA-1 eq 26),
whose steady state is the same `2 Dg`: levels matched on average, but every onset after silence
played low, by up to 2.6 in log2 amplitude (16 dB). In DVSI's re-quantizations of its decoded
D-STAR frames into rate 33, the gain implied for fully voiced frames sits a median 0.037 log2 from
the absolute law per index, against 0.35 for the differential one; the offset -0.09 is where the
rate-33 gain indices agree best (80.5% exact). The encoder quantizes the gain against the same
levels, so that DVSI's decoders play our bits at the intended level.

### D-STAR muting under channel errors

DVSI's decoder mutes D-STAR under sustained channel errors near 10% and never at 5% or below.
With hard decisions this library now mutes a frame while the error rate it
already keeps (`0.95 eps_R + 0.001064 eps_T`, with D-STAR's own error counts) is above 0.042.
On DVSI's vectors that agrees with DVSI's mute decisions on all loud frames at 1-5% and on 90% at
10%. It is an approximation: at 10% it starts muting 0.7 s into the stream, where DVSI starts at
3.3 s, since DVSI's rule integrates the errors more slowly; reproducing it exactly needs per-stream
state that `mbe_parms` cannot gain within a minor release. Soft decisions do not mute on it: at
10% DVSI mutes 53% of the loud frames of the soft stream against 88% of the hard one, and no rule
tried reproduced that. A muted frame is synthesized as usual (model, phase, overlap and noise
state advance) and played as mute noise at the level of DVSI's: DVSI's muted D-STAR output is 0 on
67% of the samples and +-1 on the rest (`dam_e10_hd`, RMS 0.60), and ours is uniform noise in
[-1.5, 1.5] on the int16 scale, which the int16 conversion truncates to the same mix
(every muted frame of `dam_e10_hd` plays at RMS 0.5-0.7, where the library's comfort noise played
56).

Tone frames follow the same mute: a frame that reads as a playable tone is muted while the error
rate is above the threshold, and one that fails the error check (2 C0 corrections, or 3 or more
in all) is a lost frame, repeated with its error count kept. On `dam_e10_hd` a channel error read
as a tone at frame 367 used to play at RMS 8497 inside a stretch DVSI keeps muted (0.5), and a
rejected tone at frame 477 reset the decoder, clearing the error rate, so our output unmuted for
the frames after it (RMS 370-562) where DVSI's stays muted (0.6). With both fixed, our mute
decisions on `dam_e10_hd`'s loud frames agree with DVSI's on 89.6% (86.4% before; missed mutes
4.1% to 1.0%, false mutes unchanged); no other vector's mute decisions change. Soft decisions do
not mute on the error rate, so on `dam_e10_sd` the rejected tone at frame 477 now repeats and the
frames after it follow the speech (RMS 2965, 2599, 1979; DVSI's decode of the clean bits 3151,
2677, 2362) where the reset played 55, 345, 494; DVSI mutes there.

### D-STAR lost frames

DVSI's test vectors never carry b0 124 or 125 outside the channel-error vectors, so how DVSI's
D-STAR decoder plays them, and tone frames with no tone, was measured on hardware for
[#106](https://github.com/arancormonk/mbelib-neo/issues/106): one ThumbDV (AMBE-3000R, firmware
V120.E100) decoding three passages of `dstar/dam.bit` and `dstar/clean.bit`, with 1, 3 or 5
crafted frames inserted after eight voice frames and 16 voice frames after them. Every capture
was repeated, and every repeat was sample-identical.

- **As a lost frame** (whole stream, voice after included, sample-identical to the same number
  of frames sent with the chip's LOST_FRAME flag): b0 125 (b2 0, 4, 32, 63); a b0 126 tone frame
  with index 0, 4, 123, 127, 148 or 255; the b0 127 frame with index 128 and volume 0 that this
  library's encoder used to send for silence. A lost frame repeats the last frame at level for
  three frames; the 4th fades out and the 5th is near silent.
- **Not as a lost frame**: the null frame (`9E 8D 32 88 26 1A 3F 61 E8`, b0 124), b0 124 with
  b1 15, b0 120..123, and b0 126 tone frames with a playable index at any volume, 0 included.

This library now decodes b0 125, b0 126 tone frames without a tone and b0 127 frames as lost
frames: they repeat the last frame exactly as a frame with too many errors does and count toward
the same repeat limit. It takes every b0 127 frame as lost because the one measured plays
differently from the same fields under b0 126. On the measured streams our decoder plays each of
them exactly as its own repeat for errors, for runs of 1 and 3; its repeats are close to DVSI's
lost frames (`dam` run of three, frame RMS: ours 4372, 4153, 4352; DVSI 4421, 4931, 4928).

Every b0 127 frame in DVSI's `dstar/*_dtx.bit` vectors is the all-ones DV frame (`FF` x 9; b0 127,
tone index 191, volume 183 after two corrections), which DVSI's decoder plays near silent (frame RMS
about 0.6) rather than as a lost frame. Those vectors were decoded with `-dtx 10000`, while the
hardware above ran with DTX off, so the all-ones frame keeps its earlier handling here (comfort
noise and a decoder reset) until it is measured with DTX off. Our encoder sent b0 124 and 125 for the
lowest pitches and for quiet frames (16 and 156 of `dam`'s 1276 frames), and DVSI's decoder plays
b0 125 as a lost frame; it now sends b0 0..123, as DVSI's encoder does. A played tone frame now
also ends a run of repeats, so a repeat after it starts a new count. On the validation vectors, with
both encodings decoded by the earlier decoder (which plays b0 125 as voice, so both are scored on
the same frames), the change moves the encoder's log-spectral distance from the input by +0.015 dB
and its envelope correlation by -0.0014, and leaves pitch, octave errors and voicing agreement no
worse (`dstar-lost-frames` gate set); DVSI's own speech bits decode unchanged.

The repeat limit follows DVSI too. On the hardware's runs of five lost frames DVSI plays the 4th as
the last frame fading out (`dam` 2316), then near silence, and keeps the prediction history: the
voice after the run follows the voice before it. The library used to play comfort noise (57) and
reset the model, so the voice after it started from scratch. It now synthesizes the 4th frame as
the last frame fading into a silent copy of itself, under the mute noise, keeps the prediction
history, and stays muted until a frame plays; while the error rate also mutes, the fade is not
played, since the speech it would fade out was never heard. From the 3rd lost frame to the 8th
voice frame after the run, the mean of `|20 log10((ours + 1) / (DVSI + 1))|` over the frame RMS
values falls from 8.2 to 3.2 dB (`dam`), 3.4 to 2.1 dB and 3.0 to 0.8 dB (two `clean` passages);
resetting the model after the fade instead was further from DVSI on all three.

Still different from DVSI on these streams: a tone frame cuts off the previous frame's tail where
DVSI lets it fade out, and the voice after the tone plays at once where DVSI's fades in. Adding
the fade-out of the previous frame under the tone (and fading the voice after it in) cuts the
mean frame-level error on zero-volume tones from about 8 to 2.5 dB, but DVSI's fade-out there is
about 5 dB quieter than this library's synthesis window gives, its fade-in is faster, and the change
moves tones with a nonzero volume (the only kind DVSI's encoder sends) further from DVSI (1.8-2.1
to 1.9-2.5 dB); it waits for waveforms of DVSI's tone transitions (#147).

### D-STAR high-band shelf

The shelf above the spectral enhancement is +2.5 dB, not +3.0 dB, which balances D-STAR's
3-4 kHz band against DVSI's once the gain is right.

### Rate 33: field assignment of u3

DVSI's rate-33 decoder and encoder place `b4(2)`, `b4(1)`, `b4(0)` and then `b3(0)` in u3's bits 8
to 5 (`ambe_d[40..43]`); the v1.0.5 text of TIA-102.BABA-1 Table 8 prints `b3(0)` first. DVSI's
rate-34 frames (no FEC) carry exactly the 49-bit vector this library deinterleaves from a rate-33
frame, so the difference lies in the field assignment, not the interleave. With DVSI's
assignment, DVSI's own re-quantizations of a decoded model pick the PRBA vectors this library's
quantizer picks (84-91% exact for b3 and b4, against 8-34% with the printed one), and our decode
of DVSI's bits comes closer to DVSI's output (log-spectral distance on development 3.92 to 3.64
dB). One table now serves the decoder's reads and the encoder's packing.

### Rate 33: pitch

DVSI's decoder and encoder use a smooth law, `log2 f0 = -4.322605 - 0.0220427 b0`. TIA-102.BABA-1
Annex A tabulates f0 up to 0.2% higher for some codes (b0 70-76) and has two oversized steps
(30-31 and 76-77). DVSI's decoded harmonics sit on the law (fit residual 1.2e-4 in log2 f0 over
81 codes), and its encoder's re-quantizations fall on the law's decision points. The harmonic
count follows the rule Annex A's own L column follows (every harmonic at or below 0.9254 of
Nyquist), which on this law differs from the column only at b0 17: 12, not 11. DVSI's output
carries the twelfth harmonic there (closer to an L = 12 decode in 9 of 11 frames). The encoder
quantizes to the nearest code in log2.

### Rate 33: error counts and repeats

The decoder counts C0 and C1 errors as clause 5.5 does: every bit of the received codeword,
parity bits included, that differs from the corrected one, plus C0's extended parity bit when it
disagrees with the corrected codeword, so that 4 means a detected, uncorrectable pattern. The
Golay decoders count only the data bits they correct, which undercounted. With the clause 5.5
counts, clause 5.7's mute (`eps_R > 0.096`) reproduces DVSI's muting on its 10% hard-decision
vector (agreement on loud frames 0.62 to 0.925, with the same first muted frame). Soft decisions
count each disagreeing bit by its reliability, in units of a fully reliable bit, so a fully
confident soft stream counts as its hard twin does and a flip the channel marked unreliable counts
little; their muting stays as close to DVSI's soft-decision muting as before (0.969 at 10%).

A frame repeats when C0 shows four errors. DVSI's decoder does not apply eq 57-58's second
condition (`eps_0 >= 2` and `eps_T >= 6`): its output repeats where `eps_0` is 4 and not where only
that condition holds, and dropping it brings the frames we play closer to both DVSI's output and
the input speech.

P25 and D-STAR keep the data-bit counts: on DVSI's vectors full counts make P25 repeat where DVSI's
output shows none, and make D-STAR repeat at 5% where DVSI does not.

### Output high-pass phase

DVSI's decoders, in every mode, roll off the band below 200 Hz like a second-order high-pass
near 100 Hz, and this library already rendered the magnitude of a second-order Butterworth
high-pass at 100 Hz for it. It now renders that filter's phase too, on every voiced harmonic of
both frames being interpolated, at render time, so the model, its regenerated phase and the
stored state keep the decoded values. On the validation vectors the phase DVSI renders differs
from ours by a median 8.9, 6.6 and 7.7 degrees (D-STAR, P25, rate 33), 4-5 below 300 Hz, where it
differed by 12.1, 10.0 and 11.5 (12-14 below 300 Hz): per harmonic of steady voiced frames, once
each frame's time shift and common phase are removed. The filter's corner and shape come from the
magnitude already measured; nothing in it was fitted to the phase.

### Voiced and unvoiced level

With the changes above in place, DVSI's decoders still play voiced frames 0.32-0.38 dB louder
and unvoiced frames 0.32-0.40 dB quieter than this synthesis, in all three modes (per-frame
medians on the development vectors, identical bits). Voiced harmonics are now rendered 0.3 dB
louder and unvoiced ones 0.4 dB quieter, applied with the high-pass, so the model and the
adaptive smoothing see the decoded amplitudes.

## Departures from the specification

| Item | Specification text | DVSI's AMBE-3000 | This library |
|---|---|---|---|
| Rate-33 u3 bits 8-5 | Table 8: b3(0), b4(2), b4(1), b4(0) | b4(2), b4(1), b4(0), b3(0) | DVSI |
| Rate-33 pitch | Annex A table | smooth log2 law | DVSI |
| Rate-33 harmonic count, b0 17 | Annex A: 11 | 12 | DVSI |
| Rate-33 repeat | eq 57-58: `eps_0 >= 4`, or `eps_0 >= 2` and `eps_T >= 6` | `eps_0 >= 4` | DVSI |
| P25 output under channel errors | eq 112-116 | attenuates progressively | specification |

The P25 attenuation is not adopted: no rule tried (eq 112-116 as printed, and instantaneous,
clamped or leaky variants of it) explains it, and under channel errors our output stays closer to
the input speech than DVSI's (log-spectral distance from the input at 10% bit errors 10.4 against
12.1 dB).

## Known differences

- D-STAR muting starts earlier than DVSI's at sustained 10% bit errors, and soft-decision D-STAR
  streams are not muted on the error rate (see above).
- After a fully unvoiced gap DVSI repositions the voiced pulse train in a way no simple rule
  reproduced; within voiced runs, including pitch steps, our pulse positions follow DVSI's.
- DVSI hard-clips at full scale; this library clips at 0.95 of full scale. Clipping at full
  scale moved no metric, so it is not adopted.
- Small pitch differences remain: P25 codes up to 0.07% from DVSI's (rounding-like), D-STAR codes
  within 0.03%.
- DVSI's rate-33 encoder sends V/UV indices 17-31 on some voice frames, where clause 4.2
  searches 0-16 and keeps the lowest index on ties. Those rows are all unvoiced, as row 16 is, so
  the bits differ and the decoded voicing does not.

## Measured results

The baseline is main's library (this branch's first commit, `ce62a89`, which adds only the tools that
read DVSI's soft-decision vectors); the candidate is this branch. Both builds' encodings are decoded by
the candidate's decoder, and both runs are scored at their joint lags
([below](#scoring-history)). Figures are means over the vectors of each partition (development 20,
validation 5, final 2); LSD is log-spectral distance, "vs DVSI" compares with DVSI's decode of the same
bits.

### Decoding DVSI's bits

| Codec | Partition | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| D-STAR | development | 4.48 → 3.70 | 0.882 → 0.952 | -0.59 → +0.05 | 7.02 → 6.58 |
| D-STAR | validation | 4.05 → 3.19 | 0.867 → 0.948 | -0.49 → +0.02 | 5.69 → 4.99 |
| D-STAR | final | 5.54 → 4.45 | 0.815 → 0.922 | -0.84 → +0.17 | 6.78 → 5.76 |
| P25 | development | 3.85 → 3.80 | 0.955 → 0.955 | -0.25 → +0.03 | 6.49 → 6.51 |
| P25 | validation | 3.22 → 3.18 | 0.956 → 0.956 | -0.32 → -0.03 | 4.91 → 4.90 |
| P25 | final | 4.18 → 4.15 | 0.930 → 0.929 | -0.13 → +0.12 | 5.65 → 5.64 |
| Rate 33 | development | 3.94 → 3.62 | 0.947 → 0.953 | -0.42 → -0.24 | 6.62 → 6.51 |
| Rate 33 | validation | 3.43 → 3.06 | 0.943 → 0.952 | -0.19 → -0.07 | 5.20 → 5.03 |
| Rate 33 | final | 4.63 → 4.20 | 0.921 → 0.932 | -0.14 → +0.07 | 6.02 → 5.80 |

The 3-4 kHz band now sits within 0.25 dB of DVSI's in every mode and partition (main: up to +1.1 dB
for D-STAR).

### Channel errors

`dam` at 10.6% bit errors (at 1.1-4.9% neither build mutes, and DVSI mutes only P25, on at most 0.5%
of the loud frames):

| Codec | Decisions | Mute agreement on loud frames | DVSI mutes | We mute | First muted frame (ours / DVSI) |
|---|---|---|---|---|---|
| D-STAR | hard | 0.119 → 0.935 | 88% | 0% → 86% | none → 34 / 165 |
| D-STAR | soft | 0.474 → 0.474 | 53% | 0% → 0% | none / 325 |
| P25 | hard | 0.986 → 0.986 | 1.4% | 0% | none / 62 |
| P25 | soft | 0.985 → 0.985 | 1.5% | 0% | none / 76 |
| Rate 33 | hard | 0.618 → 0.925 | 38% | 0% → 46% | none → 68 / 68 |
| Rate 33 | soft | 0.969 → 0.969 | 3.1% | 0% | none / 109 |

These are the figures measured for that change; the current D-STAR hard-decision agreement is in
[D-STAR muting under channel errors](#d-star-muting-under-channel-errors).

On the loud frames that neither run nor DVSI mutes, our output comes closer to DVSI's at every error
rate: by 1.17 dB of LSD for D-STAR (hard decisions; 1.08 soft), 0.36 dB for rate 33 (0.34) and 0.02 dB
for P25 (0.02).

### Encoders

Our encoding of the input against DVSI's encoding of it, both decoded by the candidate's decoder (the
one that reads DVSI's bits as DVSI does), so the baseline column is how main's encoders' bits play on a
decoder that follows DVSI:

| Codec | Partition | LSD vs DVSI's encoding (dB) | Level gap to DVSI's encoding (dB) | LSD vs input (dB) |
|---|---|---|---|---|
| D-STAR | validation | 6.09 → 4.76 | 2.79 → 0.32 | 6.09 → 4.79 |
| D-STAR | final | 7.27 → 5.98 | 2.11 → 0.54 | 6.78 → 5.61 |
| Rate 33 | validation | 4.81 → 4.61 | 0.37 → 0.13 | 5.00 → 4.83 |
| Rate 33 | final | 6.03 → 5.75 | 0.41 → 0.10 | 5.85 → 5.62 |

The IMBE encoder's bits are unchanged on every vector.

### Gates

Declared in `tools/quality/dvsi_gates.json` before the final scoring: per codec and partition, LSD
against DVSI not worse by more than 0.01 dB and envelope correlation not lower by more than 0.002, the
IMBE encoder's bits unchanged, the D-STAR and rate-33 encoders' LSD against DVSI's encoding not worse
by more than 0.01 dB; on every tone vector, detector agreement and encoder tone agreement not lower and
wrong-tone windows not more frequent; at every channel error rate, mute agreement not lower, and
hard-decision LSD on common frames not worse by more than 0.01 dB on at least a tenth of the frames.

`dvsi-parity-development`, `-validation`, `-final` and `-channel` pass. `dvsi-parity-tones` fails:
detector agreement drops on six of the 104 tone vectors (D-STAR's `knox_12`; rate 33's `dtmf60`,
`cpvbad`, `dtmfvbad`, `dtmf` and `dtmf180`), by one to eleven windows each (at most 1.7% of a
vector's steady windows), and rises on two. Every window that lost agreement covers only frames
DVSI's encoder coded as voice next to its tones (D-STAR carries KNOX tones as voice; the partial-tone
frames of the DTMF vectors and the malformed call-progress tones of `cpvbad` come out as harmonics
of one fundamental), where the detector's test, two components holding at least 90% of the window's
energy, flips. In those windows our output is closer to DVSI's than main's on `dtmf60` (LSD
13.4 → 11.5 dB), `dtmf` (10.9 → 7.7) and `dtmfvbad` (10.4 → 10.1), and slightly further on `cpvbad`
(6.5 → 6.7), `knox_12` (2.5 → 2.6, all three at the threshold) and `dtmf180` (20.9 → 21.1, at
-50 dB). Tone frames play as before (level error on the DTMF vectors within 0.002 dB of main's), no
window plays a wrong tone, and the encoders' tone frames are identical. On D-STAR's KNOX vectors,
whose tones are voice frames, the level error against DVSI's output changes by -0.25 to +0.13 dB per
vector (-0.02 dB on average over the 11 with steady tone windows).

### Tone-gate investigation

A fresh comparison of `ce62a89` with `58cbfb3` (including the SIMD optimization), using GCC 16.2.1,
SSE2, Release, fast-math and LTO, reproduces the same six failures on all 104 tone vectors at joint
lags. Counts of windows agreeing with DVSI's detector result are:

| Codec / vector | Baseline | Candidate | Windows checked |
|---|---:|---:|---:|
| D-STAR / `knox_12` | 74 | 73 | 76 |
| Rate 33 / `cpvbad` | 354 | 350 | 358 |
| Rate 33 / `dtmf` | 4466 | 4465 | 4470 |
| Rate 33 / `dtmf60` | 630 | 619 | 639 |
| Rate 33 / `dtmf180` | 4474 | 4473 | 4477 |
| Rate 33 / `dtmfvbad` | 1813 | 1811 | 1817 |

The detector accepts any one or two sufficiently pure sinusoids; it does not require a valid DTMF or
call-progress frequency pair. Its agreement denominator also includes windows where it detects no
tone in DVSI's output and neither neighbour. On the five rate-33 vectors, all 21 newly disagreeing
windows are additional detections in our speech-coded output where DVSI's output is not sufficiently
sinusoidal for this detector. Two other windows regain agreement, leaving 19 fewer matches across
those five vectors. Every newly disagreeing window contains only voice frames.

The one D-STAR loss is a missed detection in `knox_12` at 2.16 s. The two strongest components hold
90.0414% of the baseline window's spectral energy and 89.9669% of the candidate's, crossing the
detector's 90% threshold; DVSI's window holds 91.0610%. The component frequencies remain close:
739.2 and 1414.8 Hz in the candidate against 738.0 and 1412.5 Hz in DVSI's output.

Diagnostic builds reverse one change at a time, keeping the same input bits and joint lags. Restoring
the old rate-33 u3 field assignment alone recovers the baseline agreement counts on all five affected
rate-33 vectors and removes 19 of the 21 individual added detections. Two windows need combined
reversions: restoring the old field assignment, pitch table, output phase and voicing gains together
removes all 21. Restoring D-STAR's old differential gain alone recovers its `knox_12` count. Reversing
the output high-pass phase, voicing gains, rate-33 pitch table or D-STAR shelf individually does not
recover those six counts. These experiments locate the changes responsible. With DVSI parity as the
priority, the corrected field assignment, absolute gain law and synthesis changes are retained.

Across the 52,400 decoded frames in the 104 vectors, frame flags and extracted parameter bits are
unchanged. All 28,473 explicit tone frames (13,212 D-STAR and 15,261 rate 33) produce byte-identical
int16 PCM: 4,555,680 samples. Wrong-tone detections remain zero in both runs, and encoder tone
agreement is unchanged. The remaining discrepancy belongs to speech synthesis on voice-coded tone
regions, whose spectral comparisons above have mixed results. The detector gate remains failing with
its original threshold; improving DVSI waveform parity in those regions is the remaining work.

### Scoring history

The final partition was scored twice, with the library unchanged in between. The first scoring fixed
every comparison at the lag the baseline's own alignment search had found, as the gates' protocol then
required. That failed the D-STAR checks on every partition (envelope correlation 0.02-0.04 lower,
final LSD 0.10 dB worse) because main's D-STAR decoder had aligned with DVSI's output at 40-96
samples (correlation 0.970-0.993): its gain law played onsets up to 16 dB low, which pulled the
envelope correlation off. The corrected decoder aligns at -16 to 0 samples, as P25 (-16 to +16) and
rate 33 (-16 to -8) do in both builds, with higher correlation on every speech vector
(0.980-0.997), and main's own decode
comes closer to DVSI's at those lags too (development LSD 4.48 dB against 4.97 at its own). Both
runs are now scored at joint lags: each comparison at the lag of whichever run aligned with the
higher correlation
([DVSI scoreboard](testing.md#dvsi-scoreboard)). Nothing else changed between the two scorings, and the
tone gate failed in both, on windows of the same kind.
