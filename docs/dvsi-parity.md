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
  were examined while these changes were made, so their results are retrospective; for the
  October 2026 changes the final partition (`irstia`, `mark`) was scored only after they were
  fixed (twice; see [Scoring history](#scoring-history)). For the November 2026 synthesis changes
  ([Unvoiced synthesis](#unvoiced-synthesis), [Synthesis in subframes](#synthesis-in-subframes),
  [Voiced onsets and offsets](#voiced-onsets-and-offsets)) no partition was held out: the study
  that found them read `irstia` and `mark` along with the others, so every figure for them is
  retrospective. The same holds for the 2.6.9 rate-33 frame-handling and 16-bit output changes
  ([November 2026: rate-33 frame handling](#november-2026-rate-33-frame-handling)). For the three
  changes of 2.6.12 ([Encoder pitch code of unvoiced
  frames](#encoder-pitch-code-of-unvoiced-frames), [Joint choice of the AMBE spectral
  vectors](#joint-choice-of-the-ambe-spectral-vectors) and the band grid of [Rate 33: V/UV
  vectors](#rate-33-vuv-vectors)) the final partition was held out again: each was measured and
  gated on the development, validation, tone and channel partitions only, and `irstia` and `mark`
  were scored once, for the three together after all were fixed, against a gate set declared
  before ([November 2026: encoders and rate-33 rows 18 to
  31](#november-2026-encoders-and-rate-33-rows-18-to-31)). The three decoder changes of 2.6.13
  ([Phase regeneration](#phase-regeneration), [P25 pitch](#p25-pitch) and [Rate 33: mid-frame
  pitch](#rate-33-mid-frame-pitch)) held the final partition out too: the study that found them
  read the development and validation vectors, DVSI's synthetic pulse-train vectors and the tone
  vectors, never `irstia` or `mark`; each change was measured and gated on the development,
  validation, tone and channel partitions, and the final partition was scored once for the three
  together, after the review's fixes, against the gate set `decoder-phase-pitch-final`, declared
  before, which it passed ([November 2026: phase regeneration and
  pitch](#november-2026-phase-regeneration-and-pitch)). The three decoder changes of the tone model
  change set (`c736edf` to `c2ba220`; [The pulse train across an
  octave](#the-pulse-train-across-an-octave), [Tones as harmonics of a
  fundamental](#tones-as-harmonics-of-a-fundamental) and [P25 tones sent as voice
  frames](#p25-tones-sent-as-voice-frames)) hold it out as well: the study that found them read
  DVSI's tone vectors and the development and validation vectors, never `irstia` or `mark`; the
  change set was measured and gated on the development, validation, tone and channel partitions,
  and the final partition was scored once for it, after the review, against the gate set
  `tone-model-final`, declared before any build was scored, which it passed ([November 2026: tones
  and the octave step](#november-2026-tones-and-the-octave-step)). The two decoder changes of the
  mid-frame voicing change set (from `6acba79`; [P25: the bands that start with the
  frame](#p25-the-bands-that-start-with-the-frame) and
  [Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)) hold it out too: the study that found
  them read the development and validation vectors and DVSI's tone vectors, never `irstia` or
  `mark`; the change set was measured and gated on the development, validation, tone and channel
  partitions, and the final partition was scored once for it, after the review, against the gate
  set `mid-voicing-final`, declared before any build was scored, which it passed ([November 2026:
  mid-frame voicing](#november-2026-mid-frame-voicing)). Rate 33's mid-frame pitch wave
  ([Rate 33: mid-frame pitch](#rate-33-mid-frame-pitch)) holds it out too: the study that found it
  read the development and validation vectors and, for the checks of the change, DVSI's rate-33
  tone vectors, never `irstia` or `mark`, and its constants were fitted on the development vectors
  alone; the study's scratch build had been measured on the development and validation vectors
  before the gate sets were declared, in the branch's first commit, `8942924`, before any build of
  the branch was scored, so those partitions are not out-of-sample for it. The change was measured
  and gated on the development, validation, tone and channel partitions, and the final partition
  was scored once for it, after the review, against the gate set `r33-pitch-wave-final`, which it
  passed ([November 2026: rate 33's mid-frame pitch wave](#november-2026-rate-33s-mid-frame-pitch-wave)).
  The mid-frame columns change ([D-STAR and P25: mid-frame voicing by
  columns](#d-star-and-p25-mid-frame-voicing-by-columns)) holds it out too: the study that found it
  read the development and validation vectors (its per-harmonic labels), DVSI's tone, channel-error
  and DTX vectors and the DVstick 30 curves of #157, never `irstia` or `mark`, and the rule was found
  on the labels of the development and validation vectors together; the study's scratch build had
  been measured on both before the gate sets were declared, in the branch's first commit,
  `3593e8c`, before any build of the branch was scored, so neither partition is out of sample for
  it. The change was measured and gated on the development, validation, tone and channel
  partitions, and the final partition was scored once for it, after the review and its two rounds
  of fixes, against the gate set `mid-columns-final`, which it passed ([November 2026: mid-frame
  voicing by columns](#november-2026-mid-frame-voicing-by-columns)). The row-17 train change
  ([Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)) holds it out too: the study that
  found it read the development and validation vectors and DVSI's rate-33 tone vectors, never
  `irstia` or `mark`; its scratch build had been measured on the development and validation
  vectors before the gate sets were declared, in the branch's first commit, `6a32fb8`, before any
  build of the branch was scored, so neither partition is out of sample for it. The change was
  measured and gated on the development, validation, tone and channel partitions, and the final
  partition was scored once for it, after the review and its fixes, against the gate set
  `r33-row17-train-final`, which it passed; no decode of the two held-out vectors changes, so that
  partition shows no harm and cannot show the gain ([November 2026: the pulse train through a
  row-17 frame](#november-2026-the-pulse-train-through-a-row-17-frame)). Rate 33's pulsed rows ([Rate 33: the pulsed
  rows](#rate-33-the-pulsed-rows)) hold it out too: the study that found them read the development
  and validation vectors and DVSI's rate-33 tone, channel-error (hard and soft decisions) and DTX
  vectors, never `irstia` or `mark`; its scratch build had been measured on all of those, with the
  scoreboard's own measures for rate 33, before the gate sets were declared, in the branch's first
  commit, `5af6391`, before any build of the branch was scored, so none of those partitions is out of
  sample for it. The change was measured and gated on the development, validation, tone and channel
  partitions; the final partition was scored once for it, on 2026-11-27, after the two reviews and
  their fixes, against the gate set `r33-pulsed-rows-final` (the build of `1502dc0` against the
  baseline's): passed, 12 of 12 checks ([November 2026: rate 33's pulsed
  rows](#november-2026-rate-33s-pulsed-rows)). DVSI's low band ([DVSI's low band](#dvsis-low-band))
  holds it out too: the study that found it read the development and validation vectors and DVSI's
  synthetic pulse-train vectors, never `irstia` or `mark`; its scratch build had been measured on the
  development, validation, tone, channel and pulse-train vectors, with the scoreboard's own measures
  too, before the gate sets were declared, in the branch's first commit, `4677a41`, before any build of
  the branch was scored, so none of those partitions is out of sample for it. The change was measured
  and gated on the development, validation, tone and channel partitions; the final partition is to be
  scored once for it, after the review, against `low-band-final` ([November 2026: DVSI's low
  band](#november-2026-dvsis-low-band)). The four rules of DVSI's P25 decoder ([P25
  pitch](#p25-pitch), [Phase regeneration](#phase-regeneration), [D-STAR and P25: mid-frame voicing
  by columns](#d-star-and-p25-mid-frame-voicing-by-columns)) hold it out too: the study that found
  them read the development and validation vectors and DVSI's synthetic pulse-train vectors, never
  `irstia` or `mark`; its scratch build had been measured on the development, validation, tone and
  channel vectors, with the scoreboard's own measures too, before the gate sets were declared, in the
  branch's first commit, `f1be05d`, before any build of the branch was scored, so none of those
  partitions is out of sample for them. The change set was measured and gated on the development,
  validation, tone and channel partitions; the final partition is to be scored once for it, after the
  review, against `p25-decoder-final` ([November 2026: four rules of DVSI's P25
  decoder](#november-2026-four-rules-of-dvsis-p25-decoder)).
- **PESQ and ESTOI.** Where a section gives them, PESQ is ITU-T P.862 narrow-band, from the Python
  package pesq 0.0.4, on the decoded signal and the input from their starts; ESTOI is from pystoi
  0.4.1, at the lag within 60 samples of the envelope alignment, in steps of 10, that maximises it.
  Neither is part of this repository's tooling. A change is given as the mean paired change over the
  vectors ± its standard error, with the number of vectors it improves. Two sets are used: the 25
  development and validation vectors' input speech, and a second set of 70 recordings from the Open
  Speech Repository, CMU ARCTIC and LibriSpeech test-clean, resampled to 8 kHz, which DVSI's vectors
  do not contain. The 25 vectors come from 12 recordings (eight are versions of `dam`, five of one
  other recording, and two pairs; [DVSI scoreboard](testing.md#dvsi-scoreboard)), and the vectors of
  one recording do not vary independently, so their standard error is clustered by recording: for
  per-vector changes d_i with mean m over n vectors in G recordings, sqrt(G / (G - 1) * sum over
  recordings of (sum over the recording's vectors of (d_i - m))²) / n. Each of the 70 recordings is
  its own cluster, so there it is the plain standard error over the recordings; 30 of them (the CMU
  ARCTIC ones) are by two speakers and the 20 LibriSpeech ones by four, so it is an error over
  recordings, not speakers. Where a section says a ± is over the vectors, the per-vector scores were
  not kept and it is the plain standard error over the 25 vectors. An encoder change is scored with
  both encodings decoded by one decoder, named in its section. [Rate 33: V/UV
  vectors](#rate-33-vuv-vectors), a decoder change, scores its decode of DVSI's bits on the 25
  vectors, with ESTOI at the scoreboard's alignment of the baseline's decode with the input instead
  of the search above.
- **Hardware.** Where a contributor measured a real AMBE-3000 (a DVMEGA DVstick 30 or a ThumbDV)
  on constructed frames, those measurements decide between rules DVSI's vectors leave open; each
  section names the issue that carries the data. The rate-33 frame-handling sections also cite a third
  party's published measurements of an AMBE-3000R on constructed frames, from the documents of the
  [mbevocoder](https://github.com/srcejon/mbevocoder) project, by document and section; only those
  published measurements are used.

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
With hard decisions this library mutes a frame while an error level is above 342. The level is
`0.99 * level + errors`, where `errors` is the frame's corrected bit errors with the parity bits of
C0 and C1 included (the count rate 33 keeps), held in element 0 of `PSIl`, which is indexed from 1.
A constant count takes the level to 100 times that count, so 3 corrected bits a frame never mute,
and 4, 5 and 6 mute after 193, 115 and 84 frames. The rule was fitted by the reporter of
[#155](https://github.com/arancormonk/mbelib-neo/issues/155) to a DVMEGA DVstick 30 decoding 500-frame
`dam` streams with 3 to 7 corrected bits a frame, as data bits, parity bits and both, in C0, C1 and
both: the output drops 4 dB below the clean stream at frame 201, 122 and 91 for 4, 5 and 6 bits
whatever the split (this rule: 200, 122, 91), never for 3 bits, and never for 3 + 3 bits on every 2nd
or 3rd frame (nor here). Over those 41 streams the mean frame-level difference from the unit falls
from 10.30 to 2.36 dB. The rule before this one, the error rate `0.95 eps_R + 0.001064 eps_T` above
0.042 on corrected data bits only, muted from 2 data bits a frame, never on parity bits, and 0.7 s into
a 10% channel.

DVSI's own channel-error vectors, which played no part in fitting it, agree: on `dam_e10_hd` our mute
decisions now match DVSI's on all 603 loud frames (540 before), and the first loud frame muted is frame
165 in both decodes (34 before; the mute itself starts at frame 136, on quiet frames); at 1-5% neither
mutes. The log-spectral distance from DVSI's decode on the loud frames neither mutes is unchanged where
both runs measure it, and now covers all 72 loud frames DVSI plays at 10% (14 before); `dam_e1_hd` moves
from 3.097 to 3.109 dB, through the adaptive smoothing that now also sees parity errors, and `dam_e2_sd`
and `dam_e5_sd` from 3.133 to 3.121 and 3.533 to 3.522. Still different: streams where C1 holds four or
more errors that it miscorrects drop 20 dB on the unit earlier than here (the reporter's `c1par4`, `c1par6`
and C0 3 + C1 4 streams), and the unit's output departs slightly from the clean stream well before the
mute (from frame 22 or 16 at 5 or 6 bits), which no rule here reproduces.

Soft decisions do not mute on it: at
10% DVSI mutes 53% of the loud frames of the soft stream against 88% of the hard one. Counting the bits
the soft decoder corrects mutes it as the hard one (agreement 0.65, false mutes 35%); its
reliability-weighted costs reach the threshold of 342 never, and a threshold fitted to that one vector
(182) would agree on 90% of its loud frames, from frame 295 (DVSI 325), on a reliability scale that
each caller sets, so neither is adopted. A muted frame is synthesized as usual (model, phase, overlap and noise
state advance) and played as mute noise at the level of DVSI's: DVSI's muted D-STAR output is 0 on
67% of the samples and +-1 on the rest (`dam_e10_hd`, RMS 0.60), and ours is uniform noise in
[-0.75, 0.75] on the int16 scale, which the 16-bit output rounds to the same mix ([16-bit
output](#16-bit-output); until 2.6.8, noise in [-1.5, 1.5] truncated)
(every muted frame of `dam_e10_hd` plays at RMS 0.5-0.7, where the library's comfort noise played
56).

Tone frames follow the same mute: a frame that reads as a playable tone is muted while the error
rate is above the threshold. One that failed an error check (2 C0 corrections, or 3 or more in
all) used to be a lost frame, repeated with its error count kept; a tone frame with correctable
errors now plays ([D-STAR channel errors](#d-star-channel-errors)). On `dam_e10_hd` a channel error read
as a tone at frame 367 used to play at RMS 8497 inside a stretch DVSI keeps muted (0.5), and a
rejected tone at frame 477 reset the decoder, clearing the error rate, so our output unmuted for
the frames after it (RMS 370-562) where DVSI's stays muted (0.6). With both fixed, our mute
decisions on `dam_e10_hd`'s loud frames agree with DVSI's on 89.6% (86.4% before; missed mutes
4.1% to 1.0%, false mutes unchanged); no other vector's mute decisions change. Soft decisions do
not mute on the error rate, so on `dam_e10_sd` the rejected tone at frame 477 then repeated and the
frames after it followed the speech (RMS 2965, 2599, 1979; DVSI's decode of the clean bits 3151,
2677, 2362) where the reset played 55, 345, 494; DVSI mutes there. Since tone frames with
correctable errors play, that frame plays as a tone (RMS 1923) and the frames after it again follow
the speech (2236, 2028, 1782).

### D-STAR channel errors

A DVMEGA DVstick 30 (AMBE-3000, firmware V121.E100), which decodes all 104 of DVSI's D-STAR
hard-decision vectors sample for sample as DVSI's `.pcm` files, decoded crafted streams with exact
C0 and C1 bit errors ([#155](https://github.com/arancormonk/mbelib-neo/issues/155)). Each stream
was captured twice from a fresh start, and the two captures were identical.

- **Uncorrectable C0.** A frame with four errors in C0 (data bits, parity bits, the overall parity
  bit among them, or a mix) plays sample for sample as LOST_FRAME, for single frames and bursts of 3
  to 6, voice and tone frames alike. The overall parity bit tells these from three errors. The
  library decoded them as Golay (23,12) does: it miscorrected the 23 bits to a wrong codeword, or,
  with the overall parity bit among the four, played the frame as if it were clean.
- **Corrected errors.** Voice with 1+3, 2+2, 3+1 or 3+3 corrections in C0 and C1 plays as the
  clean frame (18 of 20 frames sample for sample, the other two within 0.02 dB). Tone frames
  with 0+3, 2+0, 1+2 or 3+3 corrections play exactly as the clean tone. The library repeated a
  voice frame with more than 3 corrections, and played a tone frame with 2 C0 corrections, or 3 in
  all, as lost. Four errors in C1, which Golay (23,12) miscorrects, do not make a lost frame.
- **After an uncorrectable frame**, the chip's output differs from what it plays after the same
  frames sent as LOST_FRAME, by up to 7 dB.

The frame decoders now flag a C0 word with four errors (`MBE_PROCESS_FLAG_C0_UNCORRECTABLE`,
carried to the parameter-bit API like the C0 context; the soft decoder checks its hard decisions).
Such a frame plays as a lost frame and counts toward the repeat limit. A voice frame decoded from it
still becomes the prediction history: the next frame predicts from its decoded amplitudes, and the
amplitudes a repeat plays stay those of the last good frame. The error counts, and so the mute, were
left as they were (the mute now counts parity bits too; see [D-STAR muting under channel
errors](#d-star-muting-under-channel-errors)). A frame is no longer repeated for its count of corrections, and a tone frame
with corrections plays. The hardware streams are hard decisions; soft decisions follow the same
rules, which DVSI's soft-decision vectors below favour over the earlier ones, and fully reliable
soft bits keep decoding exactly as hard ones.

Measured against the hardware, as frame levels (mean `|dB|` from the chip over the test frames and
those after them):

| Streams | Before | After |
|---|---|---|
| Four-error C0 bursts of 3-6 frames (`g8`) | 3.17 | 0.89 |
| Four-error C0 patterns, 1 or 3 frames (`f2`) | 1.40 | 0.79 |
| Tone frames with corrections (`g6`) | 7.47 | 0.54 |
| Voice with C0 + C1 corrections (`g2`) | 1.94 | 1.09 |
| C0 + C1 parity errors on every frame (`g4`) | 17.35 | 1.22 |

How each stream compares with its control, identical or not, now matches the chip's verdict in
every comparison of these groups but three: the two `g2` frames the chip changes by 0.02 dB, and
`g4`'s 3 + 3 parity errors on every frame, which change the chip's output by up to 0.05 dB and
leave ours as it was (parity errors are not all counted here). Predicting from the bad frame's
decode is what makes the frames after it differ from the LOST_FRAME version, as on the chip: over
the 400 frames after the 25 bad-frame/LOST pairs, the level difference between the two versions
follows the chip's (RMS 1.73 dB) with a correlation of 0.962 and an RMS error of 0.49 dB; repeating
with the history untouched gives 0.888 and 0.81 dB and leaves 18 of the 25 pairs sample-identical,
where the chip's never are. This fits
the measurements; it is not shown to be the chip's internal rule.

DVSI's channel-error vectors contain such frames: 1 in `dam_e2_hd`, 32 in `dam_e5_hd` and 209 in
`dam_e10_hd`. On the 32 of `dam_e5_hd` the log-spectral distance from DVSI's decode falls from 4.79
to 3.19 dB (level 3.88 to 1.86 dB), and on the third frame after them from 3.62 to 3.13 dB.
On the loud frames neither decoder mutes, the distance falls from 3.143 to 3.133 dB (`dam_e2_hd`),
3.700 to 3.521 (`dam_e5_hd`) and 6.003 to 4.816 (`dam_e10_hd`, on the 14 of the 72 loud frames
DVSI plays that qualify), and with soft decisions from 3.143 to 3.133 (`dam_e2_sd`), 3.553 to
3.533 (`dam_e5_sd`) and 6.250 to 5.610 (`dam_e10_sd`). Mute decisions agree with DVSI's as often
as before, and on one more loud frame of `dam_e10_sd`; `dam_e1` and every speech and tone vector
are unchanged, and so are the D-STAR encoder's bits for all 97 of DVSI's inputs.

Of five tone frames with four C1 errors in a row, none plays as a lost frame, and their volume
changes ramp over 80 samples as this library now plays them ([D-STAR tone
transitions](#d-star-tone-transitions)).

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
- **Not as a lost frame** (b0 124 plays as a silence frame, see [D-STAR b0 124](#d-star-b0-124)):
  the null frame (`9E 8D 32 88 26 1A 3F 61 E8`, b0 124), b0 124 with
  b1 15, b0 120..123, and b0 126 tone frames with a playable index at any volume, 0 included.

This library now decodes b0 125, b0 126 tone frames without a tone and b0 127 frames as lost
frames: they repeat the last frame exactly as a frame with too many errors does and count toward
the same repeat limit. It takes every b0 127 frame as lost because the one measured plays
differently from the same fields under b0 126. On the measured streams our decoder plays each of
them exactly as its own repeat for errors, for runs of 1 and 3; its repeats are close to DVSI's
lost frames (`dam` run of three, frame RMS: ours 4372, 4153, 4352; DVSI 4421, 4931, 4928).

The all-ones DV frame (`FF` x 9; b0 127, tone index 191, volume 183 after two corrections) plays
as a lost frame too. On the same ThumbDV, runs of 1, 3 and 20 of it after voice are sample-identical
to as many LOST_FRAME frames on all three passages, voice after included, with the encoder's DTX bit
(ECMODE `0x0800`, an encoder option) clear and set alike. Runs of 1 and 3 of the same fields sent
with clean FEC, under b0 127 and under b0 126 (index 191 has no tone), are sample-identical to
LOST_FRAME too (DTX bit clear). The frame is not something DVSI's
encoder sends on the channel: DVSI's USB-3000 host software, simulating DTX (`-enc -dtx N`), writes
all ones into the bit file for each silence frame it does not "transmit", and when decoding
(`-dec -dtx N`) it hands each all-ones frame to the chip with the CNI flag set, so the chip plays
comfort noise instead of decoding it (USB-3000 manual, `-dtx`; DVSI's `dstar/*_dtx` vectors were
made that way). DVSI's `_dtx` decodes therefore show the CNI flag: the first all-ones frame of a run
fades the last frame out and the rest play at RMS 0.57. A decoder given only the bits, as this
library is, gets the chip's lost frame. The library used to play the all-ones frame as comfort noise
(RMS about 56) and reset the decoder, which neither does. It now decodes it as the other b0 127 frames.
On the `dam` passage with 20 of them, frame RMS from the first all-ones frame: DVSI 4421, 4931, 4928,
2316, 12, then 0.5-0.6, and 2787, 2779 on the voice after; ours 4372, 4153, 4352, 2483, 0.5, then
0.5-0.6, and 1680, 3163 (mean `|20 log10((ours + 1) / (DVSI + 1))|` over the run and the 16 voice
frames 1.4 dB, 18.7 dB before). Against DVSI's `_dtx` decodes, where the first frames of a run
differ by design (a repeat against a fade), the mean over the all-ones frames falls from 29-31 dB to
1.0 dB (`dam_dtx`), 0.7 dB (`clean_dtx`) and 4.8 dB (`noisy_dtx`, whose repeats of the background
noise play at level where CNI fades out) and on the first voice frame after a run from 14-32 dB to
2.8-3.2 dB.

Our encoder sent b0 124 and 125 for the lowest pitches and for quiet frames (16 and 156 of `dam`'s 1276 frames), and DVSI's decoder plays
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
history, and stays muted until a frame plays; while the error level also mutes, the fade is not
played, since the speech it would fade out was never heard. From the 3rd lost frame to the 8th
voice frame after the run, the mean of `|20 log10((ours + 1) / (DVSI + 1))|` over the frame RMS
values falls from 8.2 to 3.2 dB (`dam`), 3.4 to 2.1 dB and 3.0 to 0.8 dB (two `clean` passages);
resetting the model after the fade instead was further from DVSI on all three.

The voice after the mute does not fade in from silence on DVSI's AMBE-3000 either
([#156](https://github.com/arancormonk/mbelib-neo/issues/156)). On a DVMEGA DVstick 30, voice after
a CNI run of five or more plays at the level of voice after as many lost frames (within 0.1 dB on
six pairs of streams; levels only, the samples were not compared), so DVSI's `_dtx` decodes stand
in for that frame's shape: over their 65 onsets after five or more all-ones frames, ours was
digital zero for its first 48 samples, where DVSI's, averaged in 16-sample blocks, starts about
15 dB below its later level in the first block and rises over about 110 samples. On the DVstick 30
the first voice frame played up to 12 dB louder than ours after loud voice, and also after a single
tone, whose lost frames are silent; the only state that can carry that energy there is the
prediction history, the voice from before the tone. The library now starts that frame from the
prediction history 10.5 dB down (gain 0.3): the last synthesized frame, which the mute left silent,
becomes the history's model at that gain (its amplitudes built from log2Ml as the decoder builds
them, then the spectral enhancement), with that model's own unvoiced overlap. Only that frame
changes; the state it commits and every frame after it are as before. A history held from an
uncorrectable C0 word, whose pitch and voicing are not kept, and the initial state keep the fade-in
from silence. A caller that keeps one state as `cur_mp` and `prev_mp_enhanced`, outside the
library's contract of three distinct objects, has no last synthesized frame apart from the frame
itself, and there the frame plays and commits its own decoded model ([Phase
regeneration](#phase-regeneration), "D-STAR held histories", says what is checked for that
arrangement). The gain is a fit to frame levels; how the chip does this is not known. The mean
`|dB|` of the first voice frame after a run of four or more:

| Data | Fade-in from silence | From the history |
| --- | --- | --- |
| DVMEGA DVstick 30 (#155, #156; 29 runs, uncorrectable C0 words included) | 6.54 | 2.39 |
| ThumbDV (#106; runs of 5 and 20 on three passages) | 3.64 | 1.85 |
| DVSI's `_dtx` decodes (110 runs), frame level | 3.19 | 2.92 |
| DVSI's `_dtx` decodes, log-spectral distance (the 57 onsets above -60 dBFS) | 3.52 | 3.42 |

Each source alone favours another gain (0.2 on the DVstick 30, 0.3 on the ThumbDV, 0.5 or more on
the `_dtx` decodes), and every gain from 0.2 to 0.5 beats the fade from silence on all three. Still
different: after the DVstick 30's one run of 30, our first voice frame plays 12 dB louder than the
unit's (26 dB quieter before); from a fresh start and after a run of 10 it plays 3.3 dB louder, as
before; and after uncorrectable C0 words, with the fade from silence, it plays 3.5 dB louder and
8.3 dB quieter on the two bursts that differ. No other DVSI vector changes, soft decisions included.

Right after a tone the lost frame repeats the tone, which ends
([#156](https://github.com/arancormonk/mbelib-neo/issues/156)). On the DVstick 30 above, after one or
five tone frames (index 32, volume 0, 180 or 255) the first lost frame plays the tone's tail fading out
over its first 80 samples (-30.7 dBFS after a tone at -23.0), the next two are digital silence, and
the repeat limit's 4th and 5th follow at the mute level. The all-ones frame sent as a normal frame
plays exactly as LOST_FRAME there, so a decoder given only the bits gets this. The library repeated
the speech from before the tone at its level (-13 to -15 dBFS on those streams). It now plays the
tone's tail fading out as when voice follows a tone, over a silent repeat, and the frames repeated after
it are silent; the prediction history stays the one from before the tone. Over every frame of those
streams (voice, the tone, the run, voice) the mean `|dB|` from the chip falls from 9.19 to 1.18 with
runs of five and from 4.22 to 0.74 with runs of two. At volume 255 the chip's second lost frame still plays at -91.5 dBFS; ours is
silent.

That holds for the single tone measured. After a dual-frequency tone, DTMF (index 128 to 143) or
call progress (144 to 147), the same unit keeps the tone playing through the first three lost frames,
within 0.5 dB of its level though not as exact repeats, and the 4th is 6 to 8 dB lower (#156, five
tone frames of index 128, 143, 144 or 147 at volume 180, and 128 at volume 255, then runs of 1, 3 and
4 all-ones frames). The library followed the single tone there, a fade and then silence. It now lets a
dual tone go on, from its phases, through the lost frames until the repeat limit, and fades it out
over the first 80 samples of the 4th under the mute noise; single tones keep the fade and silence.
Over those 15 streams (the run and the five voice frames after it) the mean frame-level difference
from the unit falls from 19.13 to 1.16 dB, the 4th lost frame from 72.4 to 1.7 dB (4.4 dB at index
144) and the first voice frame after the run from 15.1 to 5.9 dB; after a run of four, that voice
frame played 12 dB low, as the voice after any repeat-limit mute faded in from silence, and it now
plays 3.2 dB loud, started from the prediction history (above). The reporter's patch
was the starting point. Only frame levels were measured, so the fade's shape and the continuing phases
are the library's choice. A later round measured single tones 5, 50, 100 and 122 and DTMF 135 the same
way (five tone frames at volume 180, then four all-ones frames): the single tones end as index 32
does, the tail and then silence, and DTMF 135 goes on as the other dual tones do. The library follows
all five: its first lost frame is within 0.6 dB of the unit's (2.2 dB quieter after index 5), the
next two are digital silence on both (except the unit's second after index 5, at -77 dBFS), the 4th
is at the mute level (-99.1 dBFS on the unit, -95.8 here), and DTMF 135's four are within 0.2 dB.
Other volumes, runs of 2 or more than 4 and other indices were not measured. No DVSI vector changes.

### D-STAR b0 124

b0 124, the null frame among them, plays on DVSI's decoder as a frame with every band unvoiced, its
spectral amplitudes decoded and predicted as a voice frame's. This library decodes it with the AMBE+2
silence model of TIA-102.BABA-1 4.1 (w0 = 2 pi / 32, L = 14), the model that fits DVSI's output best
of those tried below; the hardware rules out the 9-harmonic alternative and supports 14 harmonics
(see "Checked on the hardware" below). The ThumbDV of #106 played ten b0 124 frames from a fresh start
(`clean` frame 110's spectral bits) at b1 15 and 0 and gains b2 12, 34 and 59, each repeated
sample-identically:

- **b1 is ignored**: b1 15 and 0 play sample for sample alike, the voice after the run included.
- **b2 is the gain**: from the second frame of the run on, b2 12 to 34 adds 19.7 dB and 34 to 59
  22.9 dB on every frame (17.2 and 22.8 dB on the first, which starts from silence).
- **The level climbs** through the run (b2 59: frame RMS 1314, 1778, 2529, ... 7022 on frames
  2-9), as the prediction from the fresh start's history builds up.
- **No steady pitch**: at most 2 of the 10 frames look periodic.

This library decoded b0 124 as voice at the pitch law's 64.4 Hz with 56 harmonics and b1's voicing.
It now decodes it with the silence model above and commits it to the prediction history, as before.
Mean `|20 log10((ours + 1) / (DVSI + 1))|` over frame RMS:

| Streams | b0 124 frames | Next 6 frames |
|---|---|---|
| Fresh start: 6 runs of 10 (their first frame left out, see below), 8 streams after | 3.00 to 1.87 dB | 2.96 to 1.62 dB |
| After voice, 1 or 3 frames (12 streams on six anchors, round 3) | 2.25 to 2.09 dB | 2.04 to 1.40 dB |
| Null frame x3 after voice (3 streams) | 4.40 to 4.09 dB | 2.19 to 1.94 dB |

Over the fresh runs' frames 2-10, the error in the share of power in each 1 kHz band falls from 2.57
to 0.35 dB (mean of the six runs). An all-unvoiced model with the 9 harmonics of the highest pitch
code matches the frame RMS a little better (1.11 dB on the fresh runs' b0 124 frames) but misses the
bands by 2.84 dB, and the hardware rules it out (below).
Leaving the silence frame out of the prediction history, as AMBE+2 does, is further from DVSI on all
three sets (6.82, 2.42 and 4.46 dB on the b0 124 frames), and the 80-sample crossfades of tones help
the frames after a fresh run but not the b0 124 frames after voice (3.23 dB).

Still different from DVSI:
- **The first frame from a fresh start**: an all-unvoiced first frame plays silent here (the noise
  generator starts from silence), where DVSI's plays RMS 394 at b2 59 (29 at b2 34, 4 at b2 12).
- **The first voice frame after b0 124** plays louder than DVSI's inside speech (`clean` anchor 110,
  three frames: 6815 against 2760). After a fresh run, the part of that frame that depends on the
  run's gain (the outgoing noise; b2 59 minus b2 12) last reaches 10 at sample 87 of the frame
  (counting from 0) on DVSI's decoder and at sample 104 here, the end of the unvoiced overlap.
- **The second of three null frames** plays at RMS 4 against DVSI's 123 (`dam`).

DVSI's vectors carry b0 124 only once, at frame 1102 of `dam_e5_hd` (from channel errors): its 112
changed frames differ from DVSI's decode by 1.63 dB against 1.61 before, and every other vector
decodes bit for bit as before.

**Checked on the hardware.** A further ThumbDV run (#106 round 5: the same fresh-start streams at
b1 0, b2 34, each field changed one at a time, every stream repeated sample-identically):

- **Not 9 harmonics**: b5 carries no coefficient at L = 9 and does at L = 14 (and at other counts
  above 9). On b0 124, b5 0 and 15 play differently (up to 14.3 dB per frame, the voice after the
  run included; 11.2 dB within the run; b5 5 and 10 differ from both); on b0 0 (L = 9) they play
  sample for sample alike.
- **Not pitch code 0**: b0 124 and b0 0 with the same bits play differently (up to 16.6 dB, the
  voice after the run included; 16.3 dB within the run).
- **Only b0 124**: b0 120..123 with the same bits play as low-pitch voice, each different from b0
  124, and leave the voice after the run louder (+28 to +29 dB on its first frame against +23.9).
- **Every spectral field is decoded**: b3, b4, b6, b7 and b8 each change the output when set to 0
  or to its maximum.
- **A 250 Hz, 14-harmonic frame**: b0 34 (250.3 Hz, L 14) with b1 0 and the same bits leaves the
  voice after the run nearly sample-identical to b0 124's from its second frame on (no sample more
  than 32 apart at b5 0, 1 at b5 15), so the state b0 124 leaves is very close to an ordinary
  14-harmonic unvoiced frame's: the strongest support for 14 harmonics. The
  run itself plays up to 3.6 dB (b5 0) and 4.6 dB (b5 15) apart, so DVSI synthesizes the b0 124 frame
  itself a little differently.

This library reproduces each of those relations: b5 changes b0 124 and not b0 0, b0 124 differs
from b0 0, and b0 120..123 stay voice. Over the 24 streams of that run, the error on frames 2-10 of
the run falls from 2.20 to 1.33 dB and on the six frames after it from 3.46 to 2.21 dB (b0 124 with
b5 0: 2.15 to 1.60 dB on the run; b5 15: 5.12 to 2.08 dB). Two differences remain: b0 124 and b0 34
leave the voice after the run at the same level here but not with the same waveform (DVSI's nearly
match), and on the run with b3 at its maximum, which plays near full scale, DVSI's level falls
faster after the run (616 against 14369 on its fourth frame); the earlier pitch-law decode misses
there too.

### D-STAR voicing bands

b1's four bits voice the 1 kHz bands from 0 Hz up. On a DVMEGA DVstick 30
([#157](https://github.com/arancormonk/mbelib-neo/issues/157)), each voiced band also reaches about
250 Hz into the band below. Two sweeps measured it: b0 0 to 60 in steps of 5, and b0 0, 10, 20, 30,
50, 70, 90, 100, 110 and 119 (18 distinct pitches), each with every b1 value, as the harmonicity of each 500 Hz band of
the output (voiced at 16 dB or more, unvoiced at 8 dB or less). The bands' upper edges sit at 1, 2
and 3 kHz, as this library had them; their lower edges sit 250 Hz lower. A free fit of each band's
run of voiced harmonics at each pitch finds those edges: an overlap in hertz fits 66 of the 72
pitch-and-band cells, one counted in harmonics 60, and the old edges 53.

The library now voices a harmonic when the bit of its band, or of the band 250 Hz (half a 500 Hz
column) above it, is set (`mbe_ambe2400_harmonic_voiced()`, used by the decoder and by the encoder's
decoder-equivalent model). Against the unit's clear band calls the mismatches fall from 186 to 56 of
2847 (84 to 30 of 1238 in the first sweep, 102 to 26 of 1609 in the second), and 250 Hz is the best
overlap in each sweep alone (240 and 260 Hz give 68 and 64; no upward overlap helps). Frame levels on
those streams move closer to the unit's (mean `|dB|` 0.608 to 0.547 and 0.769 to 0.690); the unit's
other #155 and #156 streams play unchanged, since their b1 values never voice a band above an
unvoiced one, which is the only case the overlap changes. On DVSI's vectors the D-STAR log-spectral
distance from DVSI falls by 0.009 dB (development), 0.006 dB (validation), 0.010 dB (final) and
0.022 dB (hard-decision channel errors). The encoder's b1 decisions are unchanged (its bits are
identical on all of DVSI's inputs). One tone-detector window of `dtmfvbad` flips: a strong 874 Hz
harmonic, now voiced at the bottom of an overlap, reads as a dual tone there (agreement 0.9909 to
0.9902).

Still different: 36 of the 56 are boundary bands that read as partly voiced here (between 8 and
16 dB) where the unit's are clear; 20 are opposite calls. At b0 30 and 55 the unit's
2000-3000 Hz band reaches lower still (1500-2000 Hz voiced with b1 2), and at b0 55 it leaves
2500-3000 Hz unvoiced, as if that band sat at 1500-2500 Hz; at b0 100 the unit never voices 3500 to
4000 Hz. No rule tried accounts for those.

### D-STAR tone transitions

DVSI's decoder crossfades into and out of a tone over the first 80 samples of a frame (#147). The
same ThumbDV played tones after voice on the `dam` and `clean` passages: indexes 32, 128 and 144 at
volumes 0, 30, 90, 180 and 255, one frame or five, plus zero-volume tones from a fresh start. Each
capture was repeated, and every repeat was sample-identical.

- **The speech under the tone** is the last frame repeated, as a lost frame plays it, faded out. On
  all three passages, the zero-volume tone frame divided by the lost frame at the same point gives,
  per 8 samples, 0.96, 0.92, 0.78, 0.67, 0.57, 0.46, 0.40, 0.29, 0.12, 0.07, 0.02, 0.01 and then 0
  (the second `clean` passage; the first agrees within 0.07, `dam` over its first 32 samples): a
  near-linear fade to zero at sample 80. A tone stream minus the zero-volume stream scales with the
  volume (the curves at volumes 90 and 180 are identical), so the speech under the tone does not
  depend on the volume.
- **The tone** fades in over the same 80 samples (index 32 on `clean`: 0.10 of its steady level at
  sample 8, 0.49 at 40, 0.98 at 80) and, once the voice returns, goes on under it and fades out over
  that frame's first 80 samples (0.97 at sample 4, 0.06 at 76). The voice fades in over those 80
  samples too: per 16 samples, the median ratio to the voice with no tone is 0.08, 0.22, 0.54, 0.40,
  0.94, and scattered around 1 after that. A tone that follows a tone plays on without a crossfade.
  (Which frequencies a tone plays and the phase it starts at: [Tones as harmonics of a
  fundamental](#tones-as-harmonics-of-a-fundamental).)

This library used to cut from the speech to the tone and from the tone to the voice at the frame
boundary. It now plays the first tone frame as `g(n) tone(n) + (1 - g(n)) repeat(n)` and the first
voice frame after a tone as `g(n) voice(n) + (1 - g(n)) tone(n)`, with `g(n) = n / 80` for n < 80
and 1 after it, `repeat` the last frame repeated as a lost frame plays it, and `voice` the frame as
it would play after that repeat. A tone after a muted frame has no speech under it. A tone muted by
the error level still holds its state, so the next tone plays on without a crossfade, but it leaves no
tail. Each frame's state records whether the error level muted it, so both hold when the next frame
comes from the other kind of decisions (soft decisions do not mute on it). Over the tone frames and the four voice frames after them, the mean
`|20 log10((ours + 1) / (DVSI + 1))|` over frame RMS falls on 76 hardware streams:

| Streams | Before | After |
|---|---|---|
| Tone, volume 30-255, one frame (27) | 4.15 dB | 0.76 dB |
| Tone, volume 30-255, five frames (24) | 2.80 dB | 0.61 dB |
| Tone, volume 0, one frame (7) | 12.38 dB | 0.86 dB |
| Tone, volume 0, five frames (18) | 4.01 dB | 1.50 dB |

On DVSI's 52 D-STAR tone vectors (tones as DVSI's encoder sends them, decoded by DVSI's decoder),
compared frame for frame as the hardware streams are, the same error over frames above RMS 30 falls
from 1.73 to 0.18 dB on the 1219 frames where a tone starts, from 7.99 to 0.60 dB on the 1190 where
one ends, and from 1.27 to 0.74 dB over all of them (about 20,250). Both decoders put out one frame
per frame of bits; the scoreboard's lags for these vectors (-72 to +8 samples) come from correlating
periodic tones, and at those lags the figures are 1.78 to 0.41, 8.10 to 1.29 and 1.28 to 0.79 dB.
Of the 2409 frames where a tone starts or ends, 2173 come more than 1 dB closer to DVSI's and 25
move more than 1 dB further (1891 and 89 at the scoreboard's lags, 57 of those in `dtmf`, whose
tails match at lag 0). 18 of the 25 are in `dtmfvbad`, where tones of one or two frames alternate
with voice frames that carry the tones' edges, and the tail adds to that voice in a phase DVSI's
does not share. Steady tone level, frequency and balance are unchanged (mean level error 0.096
against 0.094 dB), and DVSI's speech and channel-error vectors decode bit for bit as before. One
steady window of `cp1` (of 46) now reads as a tone where DVSI's has none, so the
`dvsi-parity-tones` gate fails on it: there, DVSI's encoder coded the tone's onset as voice, and
that voice frame fading out under the first tone frame partly cancels the tone's own fade-in (frame
RMS 4757 against DVSI's 6552; the cut played 6879).

The crossfade is specific to tones. Moving the speech synthesis window `Ws` to an 80-sample crossfade
at the start of the frame moves D-STAR speech away from DVSI's decode (mean development
log-spectral distance from DVSI's 3.70 to 3.86 dB). At the repeat limit, fading the repeat out over
the same 80 samples plays `dam`'s faded frame at 1552 against DVSI's 2316, where the current fade
gives 2483, so the repeat limit keeps its fade.

A tone that goes on at another volume, the same index in consecutive tone frames, moves to the new
level over the first 80 samples too, and the output is clipped after that. On a DVMEGA DVstick 30
(#155), five tone frames of index 32 with volumes 100, 164, 37, 245 and 117 play -45.4, -30.4,
-36.3, -3.3 and -8.2 dBFS, where switching volume at once played -28.8, -74.5, -1.7 and -45.5 on the
last four; with the ramp they play -30.4, -36.5, -3.2 and -8.3 (mean difference over the five
15.95 to 0.30 dB). Ramping after the clip instead plays the 245-to-117 frame at -9.4. On DVSI's
tone vectors, the 942 frames where a tone's volume changes play 0.250 dB from DVSI's decode in
frame level before and 0.218 dB after. Twelve of the 14 vectors with such frames come closer;
`dtmf50` and `dtmf100`, three such frames between them, move from 0.099 to 0.103 dB, and the steady
tones stay within the `dstar-tone-volume-fade` gate sets. A new tone index starts at
once: crossfading those as well moved DVSI's sweeps further from its decode (`sine0_4k` 0.08 to
0.33 dB on its 117 index changes).

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
quantizes to the nearest code in log2, except on a frame with no voiced band, which it sends at b0
119 ([Encoder pitch code of unvoiced frames](#encoder-pitch-code-of-unvoiced-frames)).

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

### Rate 33: V/UV vectors

DVSI's rate-33 decoder voices a harmonic when its own 500 Hz column of the V/UV vector (Annex B, eq
5) is voiced or when the column 250 Hz above it is: each voiced column reaches half a column (250
Hz) into an unvoiced one below it, as D-STAR's 1 kHz bands reach 250 Hz into theirs ([D-STAR voicing
bands](#d-star-voicing-bands)). The decoder, and the encoder's eqs 4 and 8, now read the voicing
that way. From 2.6.7 to 2.6.12, and in this series up to the mid-frame voicing change (`adecb29`),
the decoder also read row 17 with its lowest column voiced, where Annex B lists rows 17 to 31 all
unvoiced like row 16; the encoder searches rows 0 to 16 (clause 4.2), so row 17 only changed how
DVSI's bits play. On DVSI's speech vectors that reading brought our decode of frames coded 17 0.5 dB
of log-spectral distance closer to DVSI's (27 frames), and of voiced frames generally 0.03 dB. Since
the mid-frame voicing change row 17's own voicing is Annex B's again, all unvoiced, and as decoded
its lowest column is voiced in the first half of the frame only; the adaptive smoothing that
channel errors bring can voice harmonics of any frame afterwards, a row-17 frame's included, and
its end is then voiced ([Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)).

As measured for 2.6.7, rows 1, 3, 13, 15 and 17 played some columns neither voiced nor unvoiced
on DVSI's decoder; voiced was the nearer binary value (log-spectral distance against DVSI's decode
on frames coded 13 or 15: 6.87 dB voiced, 7.56 dB unvoiced), so the table kept Annex B's voiced
columns there. For rows 13, 15 and 17 what played so is a first half of the frame voiced otherwise
than its end ([Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)); rows 1 and 3 are not
resolved, and the table keeps Annex B's voiced columns for them ([Known
differences](#known-differences)).

**Rows 18 to 31: the band grid.** Annex B lists rows 18 to 31 all unvoiced, like row 16, and clause
4.1 gives every voice frame the fundamental of its pitch code. The published AMBE-3000R measurements
of the [mbevocoder](https://github.com/srcejon/mbevocoder) project (`docs/decoder-quality.md`,
"Unvoiced rows 18-31: the band grid") report that the chip plays these rows on a band grid whose L
bands always span 3.58 kHz, `w0 = 2.81 / L` radians per sample, the pitch code setting only L. The
decoder now does so: a voice frame of rows 18 to 31 takes that w0, with L from b0 as before, and its
unvoiced factor `0.2046 / sqrt(w0)` from it. The grid's w0 is the model's pitch like any other and
enters the history; the prediction (4.4) depends on L alone, so for the same history the next voice
frame's predicted log2 amplitudes do not change. Rows 0 to 17, silence and tone frames decode to the
same parameters as before, and lost frames are classified and repeated by the same rule, so a lost
frame after one of rows 18 to 31 repeats that frame's changed model. A voice or silence frame
after one of these rows plays differently too: its synthesis starts from the changed model (the
mid-frame model between the two and the unvoiced overlap carry it). The encoder still never sends
these rows, so its bits do not change.

DVSI's vectors bear the grid out. On the plain voice frames of the 25 development and validation
vectors (DVSI's bits, decoded by both), the log-spectral distance from DVSI's decode falls on every
group of these rows: by 0.18 dB on rows 18, 19 and 21 to 23 (1869 frames), 0.20 on row 20 (634),
0.05 on rows 24 and 25 (875), 0.10 on rows 26, 27, 30 and 31 (1591) and 0.13 on rows 28 and 29
(1101); the waveform correlation of the plain rows rises by 0.004. Taking the unvoiced factor from
the coded w0 instead cannot be told apart (-0.0086 against -0.0081 dB on the scoreboard). The
constant was checked once, on the development vectors alone (4568 frames of rows 18 to 31): the
change in mean log-spectral distance from DVSI's decode on those frames, from 5.94 dB at the coded
pitch, is

| `L w0` | 2.75 | 2.78 | 2.81 | 2.84 | 2.87 |
|---|---|---|---|---|---|
| Change (dB) | -0.088 | -0.136 | -0.133 | -0.073 | +0.031 |

Paired against 2.81 frame for frame, 2.78 is -0.002 ± 0.010 dB (standard error), so the published
value stands.

The same document's other statement, that row 17 plays as row 16, does not hold on DVSI's vectors:
played all unvoiced, row-17 frames lose 0.064 of waveform correlation with DVSI's decode and gain
0.35 dB of log-spectral distance from it, and PESQ against the input falls by 0.0022 ± 0.0011 (worse
on 17 of the 25 vectors; this ± is over the vectors, [Method](#method)). Row 17 kept its lowest
column voiced. That measurement compares two readings that each hold for the whole frame; on DVSI's
decoder row 17 voices its lowest column in the frame's first half and no band at its end, which
neither reading is, and the decoder reads it so since the mid-frame voicing change ([Rate 33: rows
13, 15 and 17](#rate-33-rows-13-15-and-17)).

On the scoreboard (joint lags, against the commit before the change), rate 33 decoding DVSI's bits,
mean paired change with its 95% interval:

| Partition | LSD vs DVSI (dB) | Waveform corr. vs DVSI | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| development | -0.0081 [-0.0130, -0.0053] | +0.00001 | +0.00004 | -0.0001 | -0.0014 [-0.0024, -0.0006] |
| validation | -0.0079 [-0.0101, -0.0043] | +0.00001 | +0.00002 | -0.0001 | -0.0021 [-0.0030, -0.0007] |

PESQ against the input ([Method](#method)) rises by 0.0035 ± 0.0016 (better on 19 of the 25
vectors) and ESTOI by 0.0002 ± 0.0001 (17 of 25); DVSI's own decode of the same bits scores 0.0217 ±
0.0046 above this library's before the change. Under channel errors (`dam` at 1.1 to 10.6%) the
mute decisions do not change and the output on the frames neither DVSI nor this library mutes comes
0.003 dB (hard decisions) and 0.005 dB (soft) closer to DVSI's. The decodes of nine tone vectors change,
but no tone metric moves. The gate sets `r33-row-grid-development`,
`-validation`, `-tones` and `-channel` were declared before scoring and pass: D-STAR and P25 decodes
and every encoder's bits are unchanged. The final partition was not scored for this change alone
([November 2026: encoders and rate-33 rows 18 to
31](#november-2026-encoders-and-rate-33-rows-18-to-31)). Beyond the grid they play rate 33's form of
the pulsed state: the row names pulsed 500 Hz channels of each subframe, whose bands play one pulse
in place of noise, and the pitch code the pulse's time ([Rate 33: the pulsed
rows](#rate-33-the-pulsed-rows)).

### Rate 33 tone transitions

DVSI's rate-33 decoder joins a tone to the frames around it over 80 samples, as its D-STAR decoder
does ([D-STAR tone transitions](#d-star-tone-transitions)); TIA-102.BABA-1 says nothing on it. The
tone fades in linearly over the first 80 samples of its first frame while the last frame's speech,
played once more, fades out; and it goes on from its phases through the first 80 samples of the
frame after its last, fading out on top of that frame, which rises from silence over the same
samples ([Rate 33 frames after a mute or a tone](#rate-33-frames-after-a-mute-or-a-tone)). A tone
that changes while it sounds changes at once. The library cut the tone in and out at the frame
edges. Which frequencies a tone plays and the phase it starts at are in [Tones as harmonics of a
fundamental](#tones-as-harmonics-of-a-fundamental).

When the crossfades came in (2.6.8), over DVSI's 52 tone vectors, in 20-sample blocks relative to
the tone's level, DVSI's medians were 0.15, 0.37, 0.63, 0.87 going in (1102 starts) and 0.90, 0.65,
0.39, 0.16 going out (1099 ends), and ours 0.14, 0.38, 0.63, 0.88 and 0.89, 0.65, 0.40, 0.15. Frame
for frame, the level of a tone's first frame came from 1.83 to 0.14 dB of DVSI's and that of the
frame after a tone from 22.1 to 0.58 dB; where a tone starts on loud speech the first 40 samples
came from 9.7 to 0.9 dB (123 frames). Since the frame after a tone starts from silence, the voice
frames right after a tone lie 0.37 dB from DVSI's level, and a tone's first frame is unchanged
([Rate 33 frames after a mute or a tone](#rate-33-frames-after-a-mute-or-a-tone)). Declared gate
that failed then (`ambe2450-tone-fade-tones`): on `dtmfvbad` the tone detector agrees with DVSI's
output on 98.96% of its windows, from 99.93%, while the frames of that vector with a 40-sample block
more than 6 dB from DVSI's fall from 459 to 256; the other 51 vectors keep their agreement, and no
speech vector changes. With the rate-33 frame-handling changes of 2.6.9, `dtmfvbad` agrees on
99.45% of its windows.

The tone a frame leaves sounding, its ID and AD, is kept in the otherwise unused element 0 of `Vl`
of `prev_mp_enhanced`, with the request that the next frame synthesized start from silence;
`prev_mp_enhanced` itself keeps the last speech model before the tone. A mute while a tone sounds
fades the tone out over 80 samples under the mute noise ([Rate 33 mutes](#rate-33-mutes)), lost
frames repeat the tone ([Rate 33 lost frames and tones](#rate-33-lost-frames-and-tones)), and a
tone after a mute fades in from silence.

### Rate 33 silence frames

TIA-102.BABA-1 4.1 decodes a silence frame (b0 124 or 125; the encoder sends 124) at
w0 = 2 pi / 32 with L = 14 and every band unvoiced, and the opening of 4.3 decodes its spectral
amplitudes as a voice frame's. It says nothing on how the frame joins the one before. This library
decoded it so and rendered it as any frame with no voiced band: the mid-frame model between the two
frames' (their log-mean magnitudes, the silence frame's voicing), so the voicing of the frame before
ended where the silence frame started.

DVSI's decoder plays the frame on 15 bands at that pitch, its amplitudes a further factor 0.2046
(13.8 dB) below the voice-frame law, predicted from the last voice frame without advancing the
prediction state. The published AMBE-3000R measurements ([Method](#method)) show this
(mbevocoder's `docs/frontend.md`, "VAD and DTX", with a standard deviation of 0.4 dB over 18
combinations of preceding frame and shape codes; `docs/decoder-quality.md`, "Silence frames").
DVSI's rate-33 DTX vectors show how the frame joins the one before: DVSI's decode of them is matched
best when the previous frame's whole model (pitch, harmonic count, voicing, amplitudes, phases)
plays on through the silence frame's first 10 ms and the second 10 ms crossfade to the silence
model.

The library now plays it so:

- Both silence codes decode on 15 bands (`MBE_AMBE2450_SILENCE_L`); b0 125, which DVSI's encoder
  does not send and none of DVSI's vectors holds, is given the same model. The residual is
  reconstructed on 15 bands and predicted from the history at 15 bands' taps by the voice frame's
  law; the factor 0.2046 scales `Ml` alone, and a silence frame still never becomes the prediction
  history. `mbe_decodeAmbe2450Parms()` returns that model.
- The decoder marks a silence frame (`MBE_FRAME_MARK_PREVIOUS_HELD`, in element 0 of `PHIl` like the
  marks of [Voiced onsets and offsets](#voiced-onsets-and-offsets)). The subframe synthesis then
  takes the previous frame's own model, unchanged, as the mid-frame model, and lets the pulse train
  run on from it, so a voiced frame right after one silence frame continues the train and one
  after two starts it afresh.
- D-STAR's b0 124 keeps 14 bands ([D-STAR b0 124](#d-star-b0-124)).

The three DTX vectors (`clean_dtx`, `dam_dtx`, `noisy_dtx`) hold 321 silence frames, each the first
frame of a pause. The pause's other frames are the all-ones frames DVSI's host software writes for
the frames its DTX simulation does not send ([D-STAR lost frames](#d-star-lost-frames); in rate 33
the all-ones frame decodes, after error correction, to a tone frame with index 180, which carries no
tone, so it is a lost frame, as it was before these changes). DVSI's `_dtx` decodes play comfort noise
through them, at the host software's request, where this library plays lost frames ([Rate 33 lost
frames and tones](#rate-33-lost-frames-and-tones)), so only the silence frames themselves are
compared. Our decode against DVSI's decode of the same bits:

| Silence frames | Main | Now |
|---|---|---|
| Frame correlation, median, after a voiced frame (53) | 0.631 | 0.979 |
| Frame correlation, median, after a frame with no voiced band (84) | 0.852 | 0.994 |
| Second half's level against DVSI's, median, after a voiced frame | +7.15 dB | +0.05 dB |
| Second half's level against DVSI's, median, after a frame with no voiced band | +9.20 dB | +0.21 dB |
| `noisy_dtx`, 750-2500 Hz, per 50 Hz bin around the frame's end, mean (range) | +12.75 dB (+9.9 to +14.2) | +0.60 dB (-0.13 to +1.21) |

Per 250 Hz band (centred on 250, 500, ... 3750 Hz), the second half of 318 silence frames lay
-1.0, -2.5, +6.9, +10.1, +9.6, +8.9, +9.9, +8.8, +7.5, +7.5, +7.0, +3.6, -2.5, -7.8 and -5.9 dB
from DVSI's, and lies -0.6, -0.8, +0.4, +0.9, +1.0, +0.4, +0.5, +0.5, +0.3, -0.0, -0.2, -0.2, -0.5,
-0.8 and -0.6 dB from it now.

The mid-frame rule was chosen on these vectors (correlations after a voiced and after an unvoiced
frame): the level and band count alone, with the ordinary mid-frame model, 0.747 and 0.958; the
voicing before lasting to mid-frame, as for `MBE_FRAME_MARK_PITCH_ABSENT`, 0.867 and 0.958; the
previous frame's whole model at mid-frame 0.979 and 0.994 (0.964 after a voiced frame with P25's
505 Hz condition added).

DVSI's vectors do not settle the band count: a single silence frame's second half still plays
mostly the model before it, and a 14-band build differs from the 15-band one by 0.4 dB in the top
band on these frames. The 15 bands rest on the AMBE-3000R measurements. None of DVSI's 27 speech
vectors, 52 tone vectors or channel vectors up to 4.9% bit errors holds a silence frame, and their
rate-33 decodes are byte-identical before and after this change; bit errors make some in
`dam_e10_hd` and `dam_e10_sd`. A run of consecutive silence frames is not verified (DVSI's vectors
hold none), and a silence frame right after a mute or a tone is in [Rate 33 frames after a mute or
a tone](#rate-33-frames-after-a-mute-or-a-tone).

### Rate 33 mutes

TIA-102.BABA-1 5.7 mutes a frame when the error rate is above 0.096, or in place of the 4th
consecutive repeat. Its recommended method computes the repeat's update equations, bypasses the
synthesis and sets the output to noise uniform over [-5, 5]; it says nothing on how a mute begins
or ends. This library did that, and also zeroed the amplitudes of the last synthesized frame, so
the frames after a mute faded in from nothing.

DVSI's decoder plays its muted frames quieter. On `dam_e10_hd` (`dam` at 10.6% bit errors, hard
decisions), over the 185 frames inside muted runs that both decoders mute, DVSI's output has a
median RMS of 0.59 on the 16-bit scale: 0 on 58.6% of the samples, +-1 on 18.4% and 18.5%, +-2 on
2.0% each. Main's played 2.46, spread evenly over -4..+4 (0 on 20.4%). And the first muted frame
after one that was heard plays the last frame once more, fading out over its first 80 samples. At
the 13 mute onsets that fall on the same frame in both decodes, the level per 20-sample block over
the previous frame's RMS (medians):

| Decoder | Blocks 1-8 of the muted frame |
|---|---|
| DVSI | 0.73 0.61 0.31 0.13 0.01 0.00 0.00 0.00 |
| Main | 0.00 in every block |
| Now | 0.84 0.58 0.24 0.12 0.00 0.00 0.00 0.00 |

The AMBE-3000R measurements agree: of four lost frames in a row, the fourth fades
(mbevocoder's `docs/decoder-quality.md`, "AMBE+2 erasures and unusable tone frames").

The library now plays a mute as DVSI's decoder does:

- **The noise** is D-STAR's, uniform in [-0.75, 0.75] on the 16-bit scale for both AMBE paths
  (`MBE_AMBE_MUTE_NOISE_AMPLITUDE`), which their 16-bit output rounds to 0 on two thirds of the
  samples and +-1 on the rest ([16-bit output](#16-bit-output)). Ours now: median RMS 0.58, 0 on
  66.5% of the samples, +-1 on 16.7% and 16.8%. P25 keeps [-5, 5].
- **The fade.** The first muted frame after a heard one plays a repeat of the last synthesized frame
  toward a silent copy of itself, as D-STAR's repeat limit does ([D-STAR lost
  frames](#d-star-lost-frames)), so the subframe synthesis fades it out over the first 80 samples,
  under the noise. The silent copy carries no frame mark, so a silence frame before the mute does
  not hold its model to mid-frame and move the fade into the second subframe. A tone left sounding
  goes on from its phases instead, fading out linearly over 80 samples under noise that does not
  fade. Later muted frames, and a mute before anything was played, are the noise alone. Over the
  first 60 samples of the 13 onsets our output correlates with DVSI's at a median 0.94 (main
  -0.08). Two candidates were measured at the six onsets on a valid frame: the last frame repeated
  and fading out correlates 0.97 with DVSI's first 60 samples, the muted frame's own synthesis under
  a falling gain 0.94. The first is used.
- **What a mute leaves.** A mute silences DVSI's output, not its decoder: the model goes on
  underneath. Through a mute, each valid voice or silence frame's own model becomes the last
  synthesized model as an unheard synthesis of it leaves it (enhanced and smoothed), and an invalid
  frame keeps the model it would have repeated; only the unvoiced overlap is cleared. DVSI's noise
  generator runs one frame per muted frame. A valid voice frame still becomes the prediction history
  when the error rate mutes it. The first frame heard again starts from that kept model, and the
  measurements behind this are in [Rate 33 frames after a mute or a
  tone](#rate-33-frames-after-a-mute-or-a-tone).

The mute decisions and the `MUTE` flag are unchanged. Still different: DVSI's muted samples include
+-2 (2.0% each); ours are 0 and +-1 only.

### Rate 33 frames after a mute or a tone

TIA-102.BABA-1 does not say what the frame after a mute or a tone starts from. This library zeroed
the last synthesized frame's amplitudes at a mute or a tone, so the mid-frame model of the next
frame was silent: the frame played nothing for 80 samples and faded in over its second half, its
pulse train running on from before, and after a tone the voice was also multiplied by the
crossfade gain.

DVSI's decoder synthesizes the first frame heard again from the model it kept ([Rate 33
mutes](#rate-33-mutes)) in the ordinary way, with two differences: the kept model itself is not
played, so the first subframe rises from silence to the mid-frame model with no unvoiced overlap
from before, and the pulse train starts anew at the frame's start. The mid-frame model is the
ordinary one between the kept model and the frame's own, frame marks included. After a tone, the
tone's tail fades out on top over the same 80 samples ([Rate 33 tone
transitions](#rate-33-tone-transitions)). The AMBE-3000R measurements agree: after lost frames and
their mute, the next voice frame plays at its level (mbevocoder's `docs/decoder-quality.md`,
"AMBE+2 erasures and unusable tone frames").

The library now renders that frame so. The subframe synthesis takes a per-call option, separate
from the rendering family (`MBE_SYNTHESIS_FROM_SILENCE`): the previous frame's own model plays
nothing and counts as having no voiced band, the pulse train the previous state carries is dropped
and a new one starts at the frame's start at phase 0, and no unvoiced overlap comes from before.
The new train's first pulse sounds at the frame's start when the mid-frame model has a voiced band,
and otherwise where the first subframe with a voiced band begins: a V/UV row 15 frame, whose first
half voices no band, launches its first pulse at sample 80 ([Rate 33: rows 13, 15 and
17](#rate-33-rows-13-15-and-17)). DVSI's noise generator runs as for any frame. The rate-33 path
uses it for the first valid voice or
silence frame after a mute or a tone, with a tone's tail added over its first 80 samples and the sum
limited, and for a lost frame below the repeat limit after a mute, which repeats the kept model from
silence; either one uses up the request.

One case rests on no measurement. A silence frame holds the previous model to mid-frame ([Rate 33
silence frames](#rate-33-silence-frames)), and that hold continues what was heard; the kept model
was not heard. From silence, a silence frame takes its own model at mid-frame instead: it rises from
silence to it over the first subframe and plays it over the second. Holding the kept model would
play speech from before the mute or the tone in the middle of a pause. No DVSI vector holds a
silence frame right after a mute or a tone where DVSI's decoder decodes the frames before it (in the
DTX vectors those frames are the host software's comfort noise), so this is a design choice.

On `dam_e10_hd`, at the 14 frames where DVSI's decode and ours leave a mute together, the rules
tried compare as follows (medians):

| Frame after a mute | Correlation with DVSI's output: first half / second half / next frame | First half's level against DVSI's |
|---|---|---|
| Main: the mid-frame model silent | 0.00 / 0.04 / 0.02 | -54.2 dB |
| The frame's own model at mid-frame | 0.84 / 0.81 / 0.63 | +2.2 dB |
| The model from before the mute kept through it | 0.70 / 0.76 / -0.01 | -4.5 dB |
| The muted frames' own models kept, the ordinary mid-frame model, a new pulse train (now) | 0.95 / 0.95 / 0.95 | +0.95 dB |

Per 20-sample block over the frame's second-half RMS, the medians are 0.14 0.27 0.46 0.56 0.99 1.04
0.75 1.17 on DVSI's decoder, 0.00 0.00 0.00 0.00 0.22 0.70 1.11 1.44 on main and 0.17 0.31 0.53 0.56
1.03 0.96 0.81 1.02 now.

After a tone: seven of DVSI's 52 rate-33 tone vectors (`alert`, `alltone`, `dtmf`, `dtmf25`,
`dtmf60`, `dtmf180`, `dtmfvbad`) have a voice frame right after a tone frame, 1206 such frames in
all. Their mean absolute level error against DVSI's decode:

| Voice frames right after a tone | First half (main → now) | Second half | Whole frame |
|---|---|---|---|
| All 1206 | 0.35 → 0.40 dB | 2.31 → 1.61 dB | 0.58 → 0.37 dB |
| 205 with a voiced band | | 1.48 → 1.06 dB | 0.66 → 0.57 dB |
| 184 of V/UV row 16 | | 3.98 → 1.34 dB | 1.46 → 0.43 dB |
| 817 of V/UV rows 17-31 | | 2.15 → 1.81 dB | 0.36 → 0.30 dB |

The second half's signed median moves from -0.41 dB to +0.03 dB against DVSI's. With the frame's
own model at mid-frame instead, the second half plays 4.4 dB loud (median), with a mean absolute
error of 4.42 dB. A tone's first frame is unchanged (1208 frames: 0.38 / 0.21 / 0.14 dB).

Still different after a tone: the pulse train of a voiced frame does not line up with DVSI's. Over
the 205 voiced frames the next frame correlates with DVSI's waveform at a median 0.26, on main and
now alike, where frames in running speech reach 0.9 and more. DVSI's decoder plays its tones through
the voiced synthesizer, and a voiced frame after a tone continues the tone's pulse train there
([Tones as harmonics of a fundamental](#tones-as-harmonics-of-a-fundamental)); the library plays
tones beside its synthesizer, so it does not. DVSI's tone vectors do hold the frame after a tone
with a loud model kept through it: 478 tone starts follow a frame whose DVSI output is above RMS 300
(`alltone` 115, `dtmf` 131, `dtmfvbad` 230, `alert` 1, `dtmf25` 1), and of the voice frames right
after a tone, 477 have such a loud frame as the model kept through the tone, so the level figures
above measure that case too. But in every one of those 478 the previous tone frame lies 2 to 35
frames back: the loud frames are the frames between tones, where the tone signal itself is coded as
voice frames. What the vectors do not hold is a tone inside speech: speech, a tone, then different
speech. For that case this library takes the mid-frame model between the last speech model before
the tone and the frame's own, as for any other.

### Rate 33 lost frames and tones

TIA-102.BABA-1 7.3 sets the amplitude of a tone frame with ID 255 to zero whatever its AD (Table 9
lists 255 as "Zero Amplitude"), and makes a frame with an invalid tone index an erasure, for which
5.6 sets the model equal to the previous frame's decoded model. Neither clause says what a repeat
plays when the frame before was a tone. This library played ID 255 as silence (fading out a tone
left sounding); after a tone, a lost frame played the tone's tail fading out, and the frames
repeated after it were silent.

DVSI's decoder does not play ID 255 as a tone. The AMBE-3000R measurements (mbevocoder's
`docs/decoder-quality.md`, "AMBE+2 erasures and unusable tone frames"; `docs/tone-frames.md`) show
a tone frame whose index it cannot play, 200 and the silent tone 255 alike, playing as a lost frame:
the last frame played repeats three times, the fourth fades, silence follows, and the next voice
frame plays at its level. The frame repeated is the last one played, so after a tone it is the
tone, three times at full level, the fourth faded. DVSI's vectors cannot show either rule: none
holds a tone frame with an unusable index or a lost frame next to a tone as DVSI's decoder plays
it. The all-ones frames of the rate-33 DTX vectors show neither: demodulated and corrected (no C0
error, three corrected in all), each is a tone frame with index 180, which carries no tone, so a lost
frame (`TONE`, `ERASURE` and `REPEAT`) before these changes as after them; only 49 all-ones
parameter bits would read as ID 255. DVSI's decodes there come from its host software's
comfort-noise request, not from the chip decoding those bits, as with D-STAR's all-ones frame
([D-STAR lost frames](#d-star-lost-frames)).

The library now follows those measurements:

- **ID 255 is a lost frame.** `mbe_classifyAmbe2450Frame()` and `mbe_decodeAmbe2450Parms()` return
  `MBE_AMBE2450_FRAME_ERASURE` for it, and the process functions flag it `ERASURE`, `REPEAT` and
  `TONE`, as any tone frame with an unusable index. The repeat criteria still come first: a frame
  they reject (four C0 errors, or, on the parameter-bit API without C0 context, the 2.1 heuristic)
  is flagged `REPEAT` alone, whatever it holds.
- **Lost frames repeat a tone.** While a tone sounds, a lost frame below the repeat limit plays the
  tone again from its phases at its level, its repeat count running; the 4th mutes, fading the tone
  out over 80 samples under the mute noise ([Rate 33 mutes](#rate-33-mutes)), and the voice after
  it starts from silence ([Rate 33 frames after a mute or a
  tone](#rate-33-frames-after-a-mute-or-a-tone)). Lost frames inside a tone are not heard: a valid
  tone after them goes on without a fade and ends the repeat run. A repeated tone carries the flags
  of its own frame (`REPEAT`, with `ERASURE` and `TONE` as the frame holds them), and the 4th adds
  `MUTE`. DVSI's noise generator stands still through a repeated tone, as through a tone frame, and
  runs through the mute that ends the run. The prediction history and the kept speech model stay
  those from before the tone.

On the DTX vectors the all-ones frames were lost frames on main too, and their per-frame flags are
unchanged. What changes in this library's own decode is what those lost frames play: the repeated
silence model is the 15-band one 13.8 dB lower ([Rate 33 silence
frames](#rate-33-silence-frames)), the 4th lost frame fades it out instead of cutting to the noise,
and the mute noise is DVSI's level ([Rate 33 mutes](#rate-33-mutes)). On `dam_dtx`, the RMS at 16
bits of the silence frame's three repeats, the 4th lost frame and the two after it was 18.7, 18.5,
17.8, then 2.35, 2.43, 2.37 on main, and is 3.8, 3.9, 3.7, then 1.74, 0.59, 0.54 now. None of that
has to do with tone ID 255. DVSI's decodes play comfort noise there, so those frames are not
compared. Playing ID 255 as a lost frame rests on the AMBE-3000R measurements alone and is recorded
as a departure from Table 9.

### 16-bit output

DVSI's D-STAR and rate-33 decoders round their output to the nearest 16-bit step; its P25 decoder
does not. TIA-102.BABA and TIA-102.BABA-1 say nothing on how the synthesized speech becomes 16-bit
samples. This library truncated toward zero in every `short` entry point (`mbe_floattoshort()`).

The share of our samples identical to DVSI's decode of the same bits, truncating and rounding to
nearest:

| Codec | `dam80` | `dam` | Mean over the 27 speech vectors |
|---|---|---|---|
| D-STAR | 0.384 → 0.546 | 0.204 → 0.281 | 0.140 → 0.177 |
| Rate 33 | 0.667 → 0.684 | 0.214 → 0.284 | 0.174 → 0.190 |

and the share of zero samples:

| Codec, vector | DVSI | Truncating | Rounding |
|---|---|---|---|
| D-STAR, `dam80` | 17.2% | 30.2% | 15.4% |
| D-STAR, `dam` | 7.0% | 12.8% | 6.5% |
| Rate 33, `dam80` | 61.3% | 66.9% | 61.5% |
| Rate 33, `dam` | 14.1% | 22.6% | 13.4% |
| P25, `dam80` | 44.4% | 24.4% | 12.3% |

The D-STAR and rate-33 `short` entry points now round to nearest, halves away from zero, after the
same 7x scale, the same limit at 0.95 of full scale and the same NaN handling as before. DVSI's
vectors tell rounding to nearest from truncation; they cannot tell one rule for halves from
another, and halves away from zero is this library's choice. Every compiled variant (scalar, SSE2,
AVX2, NEON) gives the same samples. P25, ProVoice, `mbe_synthesizeSpeech()`, the comfort noise and
`mbe_floattoshort()` keep truncating: with rounding, P25's identical share falls (`dam80` 0.264 to
0.149, `dam` 0.153 to 0.086, on a trial build), and DVSI's P25 output holds more zeros than even
truncation gives.

Measured outside the repository's tools, narrow-band PESQ (ITU-T P.862) against
the input speech rises with the rounding on all 27 rate-33 decodes (mean 2.937 to 2.958, DVSI's
2.981; `dam80` 2.548 to 2.822, DVSI's 2.847), and D-STAR's mean stays where it was (2.947 to 2.945,
DVSI's 2.931), its largest falls on `dam80` and `dam40` moving to DVSI's own scores.

The mute noise of both AMBE paths halves to uniform in [-0.75, 0.75] on the 16-bit scale, so that
rounding gives the mix truncation gave D-STAR's [-1.5, 1.5]: 0 on two thirds of the samples and +-1
on a sixth each. DVSI's muted D-STAR output on its vectors is 0 on 66.9% of the samples and +-1 on
16.5% and 16.6%.

The `short` output of the two AMBE paths is no longer `mbe_floattoshort()` of their float output. A
caller of their float entry points who wants the same samples scales by 7, limits at 0.95 of full
scale and rounds to nearest (`roundf()`, for example). The mute noise lies under one 16-bit step in
the float output, so a caller who truncates gets silence where a muted frame plays the noise alone.
`mbe_quality_eval` now takes its 16-bit output from the codec's own `short` entry points
([Evaluator interface and measurements](testing.md#evaluator-interface-and-measurements)), so the
scoreboard measures what those entry points return.

### DVSI's low band

DVSI's decoders play the band under 300 Hz in two parts, in D-STAR, P25 and rate 33 alike: a fixed
high-pass, a third-order Butterworth at 75 Hz, in magnitude and in phase; and an adaptive low cut
whose corner follows the pitch of the stream's fully voiced frames, under which a voiced harmonic is
scaled down and the unvoiced bins are removed. D-STAR's, P25's and rate 33's renderings play both
(`mbe_low_band.h`, `mbe_subframe.c`, `mbe_unvoiced_fft.c`); ProVoice, for which there is no DVSI
material, and the public `mbe_synthesizeSpeechf()` keep the second-order filter at 100 Hz that all
of them played before and no cut. What DVSI's decoders compute is not known: the filter and the
tracker below are the simplest forms found that fit their output, and every constant is a fit to
DVSI's decodes.

**Before.** From #107 the library rendered the magnitude of a second-order Butterworth high-pass at
100 Hz on every voiced harmonic and unvoiced bin, a fit to DVSI's decodes on average (DVSI's
decoders "roll off the band below 200 Hz like a second-order high-pass near 100 Hz"), and later its
phase too, at render time; on the validation vectors that brought the phase DVSI renders within a
median 8.9, 6.6 and 7.7 degrees (D-STAR, P25, rate 33) of ours, 4-5 below 300 Hz, from 12.1, 10.0
and 11.5 (12-14 below 300 Hz), per harmonic of steady voiced frames once each frame's time shift and
common phase were removed (the figures of that change, not measured again). The AMBE tones started
at that filter's phase. Under 140 Hz the average hid two populations: frames in which DVSI's lowest
harmonic lies a decibel or two above that curve, and frames in which it lies up to 12 dB below it
with no change of phase, the same frames differing between recordings of the same speech.

**The measure.** On a voiced frame with harmonics 1 and 2 voiced, a pitch within 4% of the frame
before's, the frame before voiced at harmonic 1, a fundamental under 300 Hz and DVSI's decode over 60
rms in a 320-sample Hann window centred on the frame's middle (samples 160 f - 80 to 160 f + 240):
the complex ratio of DVSI's decode to ours summed over the bins of harmonic 1 (half a harmonic
spacing either side, a 4096-point transform), in dB, less the difference of the two fixed filters at
that frequency (third order at 75 Hz less second order at 100 Hz): the *deficit* of DVSI's
fundamental. A frame the decoder flags as a repeat, mute, tone, erasure or silence frame gives no
reading. With the library's build of the branch point (`9a101f8`) as ours, on the 20 development
and 5 validation vectors and DVSI's 13 synthetic pulse-train vectors (`ptrain0` to `ptrain12`, which
the scoreboard scores in no partition): 81,342 readings, 26,765 P25 (17,139 of speech, 9,626 of the
pulse trains), 28,842 rate 33 (19,181, 9,661) and 25,735 D-STAR (16,453, 9,282); 8,515, 8,738 and
7,698 of them under -1 dB. With no model their rms is 2.593 dB (2.802 P25, 2.482 rate 33, 2.486
D-STAR; speech 2.948, pulse trains 1.757). Where DVSI's fundamental is cut it is a scaled copy of
ours: on the speech frames with a deficit under -1 dB (7,468 P25, 7,808 rate 33, 6,757 D-STAR) the
coherence of DVSI's harmonic-1 band with ours, over the band's bins, has quartiles of 0.9962 to
0.9998 (P25), 0.9959 to 0.9997 (rate 33) and 0.9976 to 0.9998 (D-STAR).

**It is state, not the frame.** `dam` frames 866 to 900 (a phrase of a male talker, f0 131 falling
to 97 Hz, after a sentence of a female one), P25: deficits of -6.7 to -10.1 dB on every frame with a
reading (-7 to -9 on all but the first). DVSI's harmonic 1 over ours (the raw ratio, by the study's
`fundts10b` window, not less the filters' difference) over frames 866 to 873 is the same at every
level of the same speech and with a quiet background (`dam` -6.5 to -7.6 dB, `dam20` -6.4 to -7.9,
`dam_10ov` -6.7 to -7.7, `damcar` -5.9 to -7.7; `dam40` and `dam80`, whose first frames lie near
the quantisation floor, -5.9 to -7.6 and -5.2 to -7.5 from frame 868), while in `dam10`, the same
speech over noise, the same frames read +0.1 to +1.7 dB, two of them (870, 871) fully voiced in both
with harmonic 1's decoded log2 amplitude within 0.6 (8.7 and 8.6 against 9.0 and 9.2). Harmonics 2
to 4 agree within 0.2 dB in both over frames 867 to 873, and harmonic 2's phase, which the phase
kernel forms from harmonic 1's decoded amplitude among others, within 0.04 rad. The frames are the
same; the rendering differs, by what the stream played before: `dam`'s female
sentence (frames 718 to 799, P25) holds 47 fully voiced frames of its 67 voiced ones, `dam10`'s none
of 62.

**It is a corner.** Within one talker the corner a `(f / fc)^2` law implies, `fc = f0 10^(-deficit
/ 40)` where the deficit is under -0.4 dB, stands still while his pitch moves: `t01`'s voiced runs
from frame 440 to 570 (P25, f0 from 204 down to 109 Hz) imply 135 to 148 Hz wherever their
fundamental lies under it; `tia11`'s first talker 102 to 177 Hz and mostly 115 to 135, its second
(from frame 1535) 143 to 213 Hz; `t04`'s runs 86 to 124 Hz.

**It follows fully voiced frames only.** The pulse trains (rate 33, DVSI's harmonic 1 over ours,
every 27th frame): `ptrain0` (period 60, 615 of 616 voiced frames fully voiced) +0.5 to +1.7 dB (mostly
+0.9 to +1.3) from its second reading on, which is the two fixed filters' difference at 133 Hz,
+1.06 dB; `ptrain11` (a sweep of
the period from 40 to 118, 1,996 of 1,997 fully voiced) rising from 0 to +3 to +4 dB, the same
difference; `ptrain5` (period 103, none of its 616 voiced frames fully voiced) -5.5 dB, then -2.0 to
-3.6 for six hundred frames: a deficit of 6.14 dB on average over its readings, where a corner of
110 Hz that never moves gives 5.98; `ptrain12`, steps of 5 in the period of about 75 frames each: at period 30 every frame is
fully voiced (the tracker climbs to 250 Hz, the corner stands at the ceiling, no deficit at 267
Hz); at 35 and 40 none or one is (the tracker holds; the fundamentals, 229 and 200 Hz, are not under
the corner); at 45 and 50 every one is again and the tracker comes down (median deficits -1.14 and
-0.30 dB; in P25 41 of the 76 frames at period 50); at 55, 60 and 65 none or one is and the corner stands where the period-50 frames left it
(median deficits -1.48, -3.17 and -4.34 dB). In speech `t01` holds 98 fully voiced frames of 1,588
voiced (P25), which is why its corner hardly moves.

**The release**, frame by frame (not used to fit the rate): `ptrain12`'s 76 frames at period 45
(177.9 Hz, all fully voiced) after the tracker stood above the ceiling (249.9 Hz at the segment's
start by the rule), deficits every sixth frame, rate 33, measured -2.1 -2.2 -2.0 -2.2 -2.1 -1.6 -1.2
-0.9 -0.7 -0.5 -0.4 -0.4 -0.2, by the rule -2.0 -2.0 -2.0 -2.0 -2.0 -1.8 -1.4 -1.1 -0.8 -0.6 -0.5
-0.4 -0.3, and with the ceiling clipping the tracker itself -1.9 -1.5 -1.2 -0.9 -0.7 -0.5 -0.4 -0.3
-0.2 -0.2 -0.1 -0.1 -0.1 (no plateau: not what DVSI plays); P25 (77 frames, the tracker at 253.7 Hz)
measured -1.9 -2.2 -2.2 -2.2 -2.2 -1.7 -1.3 -1.0 -0.7 -0.5 -0.4 -0.3 -0.2, the rule -2.1 (five
times) -1.9 -1.5 -1.2 -0.9 -0.7 -0.5 -0.4 -0.3; D-STAR, whose encoding brings the tracker to 204.1 Hz
only, measured -2.0 -2.1 -1.5 -1.0 -0.8 -0.7 -0.5 -0.2 -0.4 -0.5 -0.1 -0.1 -0.3, the rule -2.1 -1.9
-1.5 -1.2 -0.9 -0.7 -0.5 -0.4 -0.3 -0.2 -0.2 -0.1 -0.1. (The plateau at -2.2 where `40 log10(177.9 /
200)` is -2.0 allows a ceiling a hertz or two above 200.)

**The rule** the library renders (constants in `mbe_low_band.h`):

1. **The tracker.** T, 110 Hz at the start of a stream. A frame whose model, as the synthesis renders
   it (after adaptive smoothing), has every band 1..L voiced, a positive normal pitch and no tone
   model (`MBE_FRAME_MARK_TONE`) moves it toward that frame's fundamental f0 = w0 8000 / (2 pi):
   `T <- T + 0.043 (f0 - T)`, formed in double from the stored floats and rounded to float once,
   kept within 50 to 500 Hz (every pitch code of the three codecs decodes to a fundamental inside
   that range, so only a caller's pitch reaches a bound). Every other frame leaves it.
2. **The corner** a frame applies is `min(T, 200 Hz)` of the tracker the frame before left: the
   ceiling clips the corner, not the tracker, and a frame's own step shows from the next frame on.
   The mid-frame model and the frame's own model take that corner; the last frame's own model takes
   the corner its frame applied, so that a model at a frame's end plays at one corner on both sides
   of it.
3. **Voiced harmonics.** A voiced harmonic of a model at frequency f under the model's corner is
   scaled by `(f / corner)^2`, its phase untouched; one at or above the corner is not scaled at all
   (`mbe_low_cut_gain()`). f is formed once, w0 8000 / (2 pi) rounded to float times l rounded to
   float, and the fixed filter and the comparison with the corner both read that float. The bands of
   the pulsed state (`BAND_PULSED`) are scaled alike (below).
4. **Unvoiced bins.** In both subframes' transforms the bins whose centre frequency, `k 8000 / 256`
   Hz, lies under the subframe's corner (the mid-frame model's, the frame's own's) are zero, after
   the band they belong to has been scaled against its own noise power: the bins above keep the
   level they had. The noise generator runs as before.
5. **The fixed filter.** `x^3 / sqrt(1 + x^6)` and the argument of `s^3 / ((s + 1)(s^2 + s + 1))` at
   `s = j x`, `x = f / 75 Hz`, in double, rounded to float once (`mbe_low_band_dvsi_gain()`,
   `mbe_low_band_dvsi_phase()`), on every voiced and pulsed harmonic and, as a second table of the
   thread-local transform plan, on every unvoiced bin; the second-order table stays for ProVoice.
6. **What it does not touch.** The decoded models, their regenerated phases, `log2Ml`, the stored
   amplitudes and everything the encoders read keep the decoded values: the cut and the filter are
   properties of the rendering.

**The constants.** Fitted on the readings above (the study fitted them on its scratch build's;
re-measured here on the library's), the rms of the measured deficit less `min(0, 40 log10(f0 /
corner))` over all 81,342: with the library's constants 0.6148 dB (P25 0.558, rate 33 0.653, D-STAR
0.627; speech 0.675, pulse trains 0.483), 0.6145 with a ceiling of 198 Hz. With a ceiling of 198: a
rate of 0.040 / 0.042 / 0.044 / 0.046 / 1/24 0.6182 / 0.6146 / 0.6154 / 0.6197 / 0.6149; a ceiling of
196 or 202 Hz 0.6157 or 0.6163; a first value of 100 / 105 / 115 / 120 Hz 0.6955 / 0.6408 / 0.6182 /
0.6492; two rates in log frequency (0.045 up, 0.035 down) 0.6220; one rate in log frequency at best
0.6568; the corner read after the frame's own step instead of before 0.6141 (no difference: the
library reads it before); the ceiling clipping the tracker 0.6372; a Butterworth magnitude
`(f / fc)^2 / sqrt(1 + (f / fc)^4)` in place of the sharp knee 2.227; the power 1.8 / 1.9 / 2.0 /
2.1 / 2.2: 0.6551 / 0.6223 / 0.6145 / 0.6328 / 0.6750; frames counted as fully voiced from a voiced
share of 0.9 / 0.95 / 0.98 / 1: 0.6967 / 0.6167 / 0.6145 / 0.6145; a floor on the corner of 75 / 90
/ 100 / 110 Hz: 0.6247 / 0.8042 / 1.0754 / 1.4377 (none needed). By the size of the predicted
deficit the model is unbiased: predicted 0, measured -0.02 dB (45,977 readings, rms 0.39, the
measure's own noise); predicted -0.35 (10,972), measured -0.18; -1.97 (10,541), -1.97; -4.53 (9,805),
-4.55; -7.24 (3,135), -7.07; -10.66 (904), -10.28. D-STAR's pulse trains, 0.727, carry level
differences of their own that the cut does not explain (`ptrain1` +0.91 dB on average at no
predicted deficit, `ptrain11` +0.80 against a predicted -0.21, `ptrain5` 0.88 dB beyond its
predicted -5.90): they stand with the library's cut
too (below). The rate was fitted on speech; 0.042 to 0.044 fit alike, 0.043 is the library's.

**Harmonic 2** where the model puts 2 f0 under the corner (51 P25, 48 rate-33, 26 D-STAR readings):
predicted -2.09 / -2.31 / -1.56 dB, measured -2.47 / -2.02 / -1.76, correlation 0.90 / 0.83 / 0.92:
the cut is a function of frequency, not of the harmonic's number.

**Unvoiced bins.** Frames with no voiced band at the default pitch code (P25 248.1 Hz, whose first band
starts at 124 Hz; D-STAR 250.4 Hz) between two such frames, DVSI's decode over 40 rms: the transfer
function from our decode to DVSI's by cross-spectrum over the frames (480-sample Hann, 1024-point
transform; the two coherent because the noise is DVSI's own generator), grouped by the corner the
library's rule gives the frame (from the readings' series). With the branch point's build as ours:
P25, a corner under 135 Hz (2,216 frames), +0.2 to +1.8 dB from 102 to 211 Hz; at 198 Hz or more
(1,116 frames), -28 to -31 dB from 102 to 148 Hz (coherence 0.05 to 0.06), -14.5 at 164, -6.9 at
180, -3.8 at 195, -1.1 at 211, 0 from 227 up; D-STAR (2,537 and 1,046 frames) +0.1 to +1.6 dB, and
-26 to -38 dB from 102 to 148 Hz, -15.8 at 164, -7.0 at 180, -3.5 at 195, -2.0 at 211. With the
library's build, the bins under the corner zero: at a high corner P25 -1.1 to +1.0 dB from 102 to 211
Hz (coherence 0.60 to 0.88), D-STAR -1.3 to +3.5 (+3.5 at 164 Hz, +3.3 at 180; coherence 0.58 to
0.92); at a low corner P25 -0.4 to +1.7, D-STAR -0.1 to +0.6. The cut on the unvoiced bins removes
the -15 to -38 dB; what is left near the corner, where DVSI keeps more than the library's whole bins
let through, is not modelled (a corner of 0.9 or 0.8 of the tracker's, or a scale `(f / fc)^4`, left
-4 to -18 dB in the study's scratch build).

**The fixed filter.** Re-fitted on the library's own decodes, the cut in, on the frames whose
fundamental the library's tracker puts at or above the corner (the steady voiced frames of the
measure above with harmonics 1 to 3 voiced, f0 under 300 Hz), each harmonic under 300 Hz (9,933 P25,
11,717 rate-33, 10,028 D-STAR readings): DVSI's level and phase over ours against what a Butterworth
high-pass of order n at fc would give in place of the third order at 75 Hz, each with a free
constant (a level offset, a phase offset), least squares. The phase picks the third order: at its
best, P25 73.0 Hz (within 1% of the least squared error from 72.0 to 74.0), rate 33 75.5 (74.5 to
76.0) and D-STAR 75.5 (75.0 to 76.0), rms 0.0555, 0.0475 and 0.0324 rad; the second order fits worse
at any corner (0.0589, 0.0536, 0.0411, near 100 Hz), the fourth order a little worse at about 60 Hz
(0.0558, 0.0479, 0.0329), where its level would want 80 Hz. The level, flatter and noisier, puts the
third order's corner at 78.0, 78.0 and 76.5 Hz (74.0 to 81.0 within 1%). One corner for the three
codecs, 75 Hz. By whole-file correlation with DVSI's decode at lag 0 (the 20 development / 5
validation vectors, mean change against 75 Hz, the corner alone changed in library copies): P25 65 Hz
-0.003483 / -0.004521, 70 -0.001324 / -0.001618, 80 +0.000463 / +0.000296, 85 +0.000052 / -0.000740;
D-STAR -0.002129 / -0.003281, -0.000623 / -0.000926, -0.000291 / -0.000536, -0.001508 / -0.002530;
rate 33 -0.001133 / -0.002425, -0.000137 / -0.000510, -0.000750 / -0.000927, -0.002396 / -0.003290. P25's
correlation alone would take 80 Hz: DVSI's harmonics there lead ours by about 0.04 rad that is not
the filter's (P25's phase offset above), and a higher corner hides part of it. The correlation answers
to the filter's phase tail over the whole band, not to the low harmonics, and is not used to set
the corner.

**On the library.** With the low band in, DVSI's harmonic 1 over ours on the same readings has an
rms of 0.600 dB (2.593 with the branch point's build as ours and no model): P25 0.546 (speech 0.637,
pulse trains 0.324), rate 33 0.637 (0.757, 0.274), D-STAR 0.612 (0.532, 0.731); by the predicted
deficit, 0 (45,977 readings) a mean of -0.017 dB and an rms of 0.390, the measure's own level, and
from 0.645 at a predicted -1 to 0 dB to 1.947 under -9 dB (904 readings), the means within 0.16 dB.
The vectors left with the most: P25 `noisy` 1.200 dB (a mean of -0.50), `ucarf15a` 1.087, `ucarm15a`
0.913, `t04` 0.859; rate 33 `ucarm20a` 1.199, `fambm22a` 1.123, `ucarm15a` 1.106; D-STAR `ptrain5`
1.049 (-0.90), `ptrain11` 1.035 (+1.02), `ptrain1` 0.993 (+0.91). By fundamental, P25's harmonic 1
over ours on the steady frames (17,139; the bands 60 to 300 Hz of `corner10w`): -0.03, -0.48,
-0.20, -0.23, -0.12, -0.06, -0.01, -0.05 dB (it was -2.79 to +0.10) with a phase of +0.075, +0.022,
+0.026, +0.019, +0.025, +0.026, +0.048, +0.048 rad (+0.492 to +0.039); its rms 0.637 dB (2.253),
beside harmonic 2's 0.516 (0.566) and harmonic 3's 0.641 (0.640). Rate 33: harmonic 1 rms 0.757
(1.950), D-STAR 0.532 (1.867).

**Which frames move the tracker, path by path** (each the library's choice where DVSI's vectors do
not isolate the case; `test_low_cut_*_path` in `test_subframe_synthesis.c`): the rule reads the model
the synthesis renders or advances, so a frame moves it exactly when it renders or advances a fully
voiced model.

- A lost frame's repeat renders the repeated model: it counts as that model does (D-STAR, P25, rate
  33), and so does D-STAR's repeat-limit mute and rate 33's first muted frame after a heard one, whose
  fades render a silent copy of the last model with its voicing kept. P25's 5th lost frame installs
  the repeat limit's reset model, with no voiced band: it holds.
- A frame the error rate mutes: P25's is advanced (`mbe_advance_subframes()`), D-STAR's synthesized
  before its output is replaced, rate 33's synthesized unheard: each counts as its model does. Rate
  33's muted frames after the first, which keep the last model and render nothing, hold.
- Rate 33's silence frames (no voiced band) hold. A frame that starts from silence (after a rate-33
  mute or tone) renders its own model and counts as it does; the tracker goes on from the state the
  mute or the tone left (it is no part of the pulse train, which is dropped).
- D-STAR's first frame after the repeat limit's mute renders from the resumed history, which stands
  in the last frame and moves nothing; the frame counts by its own model.
- Tones: an AMBE tone frame renders no speech and leaves the state (D-STAR's tone's first frame renders
  the repeat of the last frame under it into the stream's state, and so counts as that model; rate
  33's renders it in a scratch state, and holds). A P25 tone played as its tone model holds (its mark;
  it is never fully voiced either).
- Frames of the pulsed state, D-STAR's and P25's codes and rate 33's rows 18 to 31, have no voiced
  band: they hold.

**The pulsed bands.** Whether DVSI's decoders cut the bands of the pulsed state cannot be read from a
pulse's deficit as above; on the frames where cutting them changes our decode (the pulsed bands under
the corner: 1,911 D-STAR, 2,793 P25 and 2,015 rate-33 frames of the 25 vectors), the residual energy
of our decode against DVSI's at lag 0, pooled over those frames, over DVSI's energy there: cut 0.05200
/ 0.03568 / 0.17433, uncut 0.05353 / 0.03687 / 0.17775, the cut version nearer DVSI's on 1,284 / 1,700 /
1,322 frames and the uncut on 593 / 1,026 / 668 (two library copies alike but for that). The library
cuts them as voiced bands.

**The state.** One stream, one tracker: it travels in the stream's subframe state beside the pulsed
state's slots, `noiseOverlap[38]` (`MBE_LOW_CUT_TRACK_SLOT`, the tracker the frame leaves) and
`noiseOverlap[39]` (`MBE_LOW_CUT_CORNER_SLOT`, the corner the frame applied), written by rendering
and by `mbe_advance_subframes()` alike and carried wherever the subframe state is (`mbe_copy_generator_state()`,
`mbe_subframe_state_adopt()`, which clears a state the other rendering wrote, the repeats, mutes, tones
and starts from silence above). Each is read by its bit pattern, only in a state with the subframe
rendering's mark, and only within what a frame stores (the tracker 50 to 500 Hz, the corner 50 to 200);
any other value, 0, a negative or subnormal value, an infinity or a NaN, a state without the mark (a
fresh one, or one the TIA rendering wrote, whose noise samples sit there), reads as a fresh stream's:
the tracker at 110 Hz, the last corner the frame's own. ProVoice's rendering writes 0 in both and reads
neither. A caller's value chooses at most a corner from 50 to 200 Hz, which scales voiced harmonics
under 200 Hz and zeroes unvoiced bins under 200 Hz and nothing else: no index, loop bound or division
depends on it beyond the ratio `f / corner` of a corner checked to be 50 Hz or more and the bins
compared with it as floats.

**What is not known.** What DVSI's decoder computes (an adaptive high-pass of some kind; the rule is a
fit); the tracker's exact law (one rate in frequency and two in log frequency fit within 0.01 dB of
each other on these readings, one rate in log frequency 0.04 dB worse) and whether it steps once a frame or once a subframe; whether the ceiling
is exactly 200 Hz; what the tracker does on repeated, muted and silence frames and through tones (no
vector isolates them: the library's choices above); what DVSI plays near the corner in the unvoiced
bands (above); whether DVSI's decoders cut the pulsed bands as voiced ones (the decodes come nearer with
it, by little); ProVoice; and what an AMBE-3000 plays on streams from other encoders (the rule reads the
decoded model only, so it applies to them as written). Under channel errors the adaptive smoothing can
voice every band of a frame: whether DVSI's tracker follows the decoded or the rendered voicing could
hardly be told apart on DVSI's channel-error vectors ([November 2026: DVSI's low band](#november-2026-dvsis-low-band)).

**Tones.** The AMBE tones start at the third-order filter's phase (`tone_start_phase()`,
`mbelib.c`). On DVSI's rate-33 and D-STAR tone vectors, a DTMF or KNOX tone (IDs 128 to 163) started
after a frame that is not a tone frame and has no voiced band, its phase extrapolated to the first tone
frame's start from a least-squares fit of the run's steady frames (runs of five frames or more; 68
components, 1,540 onsets, both codecs pooled per component): the measured phase less the prediction
`-16 F + p(F) / (2 pi)` is +0.00178 cycle on average over the components with the second-order
filter's p (rms 0.00285) and +0.00035 with the third-order one's (rms 0.00210); under 800 Hz +0.00298
against +0.00104. The spread within a component is 0.00008 cycle (median standard deviation). The tone
oscillators take no magnitude of either filter and no cut, as before.

### Voiced and unvoiced level

With the changes above in place, DVSI's decoders still play voiced frames 0.32-0.38 dB louder
and unvoiced frames 0.32-0.40 dB quieter than this synthesis, in all three modes (per-frame
medians on the development vectors, identical bits). Voiced harmonics are now rendered 0.3 dB
louder and unvoiced ones 0.4 dB quieter, applied with the high-pass, so the model and the
adaptive smoothing see the decoded amplitudes. These two gains now belong to the public
`mbe_synthesizeSpeechf()` alone: the codec paths render in subframes
([Synthesis in subframes](#synthesis-in-subframes)), with a voiced gain of +0.55 dB and no
unvoiced gain.

### Unvoiced synthesis

DVSI's decoders, in every mode, draw the noise of the unvoiced bands from one 16-bit generator
and synthesize it a 10 ms subframe at a time. TIA-102.BABA 11.2 (eq 117-126) has another
generator and one 256-sample window per frame.

- **The generator** is `x <- 173 x + 13849 (mod 65536)` from 0, read as a signed 16-bit sample;
  output sample t of a stream takes the state after t + 1 steps. The unvoiced parts of DVSI's
  decodes repeat after 65536 samples, and they are the same sequence at the same position in all
  79 speech vectors of the three modes. The generator runs one step per output sample through
  voice, silence, repeated and muted frames, and stands still through tone frames: on `dtmfvbad`
  it falls 160 samples behind for every tone frame, in D-STAR and in rate 33.
- **A subframe** is 160 noise samples under a square-root triangular window two subframes long, a
  256-point transform, each band scaled to its amplitude against that band's own noise power
  between edges rounded to the nearest bin, the inverse transform, the same window again, and
  overlap-add at a hop of 80 samples. A frame's own model stands at its end; the subframe at
  mid-frame takes a model between the previous frame's and this one's
  ([Synthesis in subframes](#synthesis-in-subframes)).
- **The level** is the specification's scale carried to the shorter window, without the 0.4 dB
  unvoiced render gain; the output high-pass acts on every bin, in magnitude and in phase (in
  D-STAR, P25 and rate 33 since the low-band change DVSI's third-order filter at 75 Hz, and the bins
  under the adaptive low cut's corner are zero: [DVSI's low band](#dvsis-low-band)).

The codec paths render this; the public `mbe_synthesizeSpeechf()` keeps eq 117-126 and its own
generator. `mbe_setThreadRngSeed()` no longer reaches the codec paths' unvoiced noise, since every
DVSI stream starts from the same state. The generator's state stays in `noiseSeed` and the second
half of a frame's own subframe in `previousUw`, so `mbe_parms` is unchanged. A mute moves the
generator one frame on and a tone frame leaves it where it was. D-STAR keeps the overlap of the
frame before the tone with it, so the noise goes on after the tone where it stopped; rate 33
clears that overlap, and the frame after a tone starts from silence
([Rate 33 frames after a mute or a tone](#rate-33-frames-after-a-mute-or-a-tone)). A rate-33 lost
frame that repeats a tone leaves the generator where it was too
([Rate 33 lost frames and tones](#rate-33-lost-frames-and-tones)). The two renderings keep different things in `previousUw` and `noiseOverlap`: a
state this one wrote is marked in the last element of `noiseOverlap`, and neither plays the
overlap the other left, so a stream moved from one to the other fades in.

On DVSI's vectors, at joint lags against main, the log-spectral distance from DVSI's decode
(development / validation / final):

| Codec | Before (dB) | After (dB) |
|---|---|---|
| D-STAR | 3.21 / 2.72 / 3.89 | 2.05 / 1.94 / 2.68 |
| P25 | 3.28 / 2.73 / 3.56 | 2.26 / 1.96 / 2.41 |
| Rate 33 | 2.97 / 2.54 / 3.59 | 2.02 / 1.86 / 2.47 |

Every vector comes closer. A frame with no voiced band between two such frames now correlates
with DVSI's waveform, sample for sample, at a median 0.999 for D-STAR, 0.998 for P25 and 0.992 for
rate 33, where the specification's rendering did not correlate at all. Under channel errors the
frames neither decoder mutes come 0.28-0.57 dB closer and the mute decisions are unchanged.
Declared gate that fails (`unvoiced-subframes-tones`): rate 33's `dtmf180` loses 0.3% of its
detector windows.

### Synthesis in subframes

DVSI's decoders synthesize a frame as two 10 ms subframes, the voiced part too. TIA-102.BABA
(eq 127-141) renders a frame in one piece: harmonics 1-7 at steady pitch interpolate in amplitude
and phase, the others crossfade under the window `wS`, and the phase integrates the pitch through
every frame.

- **Two models per frame.** A frame's own model stands at its end. A model between the previous
  frame's and this one's stands at mid-frame: its pitch the geometric mean of the two (in rate 33,
  between two voiced frames of different pitch, within -0.111% and +0.059% of it: [Rate 33: mid-frame
  pitch](#rate-33-mid-frame-pitch)), its harmonic count that pitch's, each band's voicing in D-STAR
  and P25, since the mid-frame columns change, the same state as in both frames where the two agree
  and otherwise set by 500 Hz columns read in both frames ([D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns);
  before it this frame's at the band's frequency, except where a frame mark said otherwise: [Voiced
  onsets and offsets](#voiced-onsets-and-offsets) and, since the mid-frame voicing change, a P25
  frame after a frame with no voiced band, [P25: the bands that start with the
  frame](#p25-the-bands-that-start-with-the-frame)), and in rate 33 this frame's at the band's
  frequency, except for V/UV rows 13, 15 and 17 ([Rate 33: rows 13, 15 and
  17](#rate-33-rows-13-15-and-17)), its
  magnitudes the two frames' at that frequency, met in the log domain. Unvoiced magnitudes enter
  that mean without their factor `0.2046 / sqrt(w0)` and take the mid-frame pitch's on the way
  out. With one frame voiced and the other not, the mid-frame pitch is the voiced frame's.
- **Linear crossfades.** Each subframe plays the model at its start with a falling weight and
  the model at its end with a rising one.
- **Pitch pulses.** The voiced part is a train of pulses. Each pulse plays a model's voiced
  harmonics, their phases the regenerated ones ([Phase regeneration](#phase-regeneration)),
  peaking 16 samples after the pulse, through a window that rises over 16 samples (a raised
  cosine), stays whole until one pitch period after
  the pulse and falls over the 16 samples after that. A response so lasts a period and 16
  samples and overlaps the next by 16, where the two windows add up to one: the responses of a
  steady train add up to the harmonics. The pulse rate moves linearly from one model's
  pitch to the next's across a subframe, after the train is taken an octave down or up where
  the new period is within 20% of twice or of half the old ([The pulse train across an
  octave](#the-pulse-train-across-an-octave)); into or out of voicing it keeps the voiced
  model's. A subframe with no voiced band at either end starts the train afresh. From one frame to
  the next the train runs on from a frame with a voiced band, from one marked as lacking a pitch or
  holding the model before it to mid-frame, and, in the rate-33 rendering, from a V/UV row-17 frame
  ([Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)); after any other it starts afresh.

DVSI describes a synthesis of this shape in US 8,200,497: square-root triangular windows over
10 ms subframes, and an impulse sequence convolved with a pitch-dependent impulse response. The
rules and constants here were measured on DVSI's decoders, from its vectors and, for D-STAR, on
the hardware below.

The codec paths render this, with a voiced gain of +0.55 dB: fully voiced frames then play
within 0.02 dB of DVSI's (frame medians on the development vectors, all three modes; 0.09-0.12 dB
low at the +0.45 dB of the rendering before). The pulse phase and the pulses still sounding
travel in `noiseOverlap`, which these paths no longer use for noise (in D-STAR and P25 also whether
the frame's mid-frame model emitted the pulse of the pulsed state, [D-STAR and P25: the pulsed
state](#d-star-and-p25-the-pulsed-state), and in rate 33 the row and pitch code of a frame whose own
subframe is pulsed, [Rate 33: the pulsed rows](#rate-33-the-pulsed-rows), and in all three the
adaptive low cut's tracker and the corner the frame applied, [DVSI's low band](#dvsis-low-band)), and
`PSIl` and `PHIl` hold the harmonic phases at the frame's end as before. The public `mbe_synthesizeSpeechf()` keeps
the specification's rendering.

This replaces two approximations of the same behaviour. The fourth-root amplitude law (a
harmonic's fourth root linear over the frame, [#157](https://github.com/arancormonk/mbelib-neo/issues/157))
is close to two linear crossfades through a geometric mean; the phase rules (the pitch moving 16
samples after mid-frame, a restart at `l w0 (N - 16)` on onsets) are what the pulse train gives
at the frame's end.

**On hardware.** A DVMEGA DVstick 30 played 179 constructed D-STAR streams for #157, each a step
of one or two fields inside a held frame. Mean error against the unit's curves, per 16-sample
block of the step frame, main and now:

| Steps | Measure | Main | Now |
|---|---|---|---|
| 30 gain steps, 7 to 42 dB | level change as a fraction of the settled change | 0.022 | 0.009 |
| 9 from voiced to unvoiced | difference from the unstepped stream, relative | 0.340 | 0.086 |
| 9 from unvoiced to voiced | same | 0.443 | 0.095 |
| 3 pitch steps | same | 0.231 | 0.123 |
| 27 double steps | same | 0.406 | 0.341 |
| 18 steps of one band's voicing | same | 0.294 | 0.342 |

The single-band steps come out worse; the mid-frame model of a D-STAR frame here takes this frame's
voicing for every band. (An earlier reading of these curves said when the unit's lowest and higher
bands start; the curves are the full-band rms of the difference in 16-sample blocks, so they show
how a step as a whole comes out, not which 500 Hz band starts when.)

The mid-frame columns change ([D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns)), fitted to labels of DVSI's
decodes, has its sequel here: by the rule a band's lower column and its upper column can take
different states at mid-frame (at an onset of the lowest band alone its upper column starts late,
its lower one with the frame), and the mean error of the single-band steps' curves falls. The same
steps and method, the head before the change (`c3c6d42`) and with it, mean errors over the steps of
each kind, to five decimals:

| Steps | Head | With the columns |
|---|---|---|
| 18 steps of one band's voicing (relative) | 0.34565 | 0.18882 |
| 9 from unvoiced to voiced (relative) | 0.03372 | 0.03374 |
| 9 from voiced to unvoiced (relative) | 0.06052 | 0.06052 |
| 30 gain steps (level / relative) | 0.00821 / 0.00626 | 0.00821 / 0.00626 |
| 3 pitch steps (relative) | 0.09050 | 0.09050 |
| double steps (15 with a level curve / 24 with a relative curve) | 0.20285 / 0.25641 | 0.20285 / 0.25641 |

The onset steps' mean moves in the fifth decimal (0.033721 → 0.033742; 0.034 in both to three
decimals); the offset, gain, pitch and double steps' means are equal to nine decimals in the two
builds (the decodes of some of those streams differ, their curves' means do not). 0.189 is the mean
of the 18 steps' errors, not a bound on each. The 18 as six means of three positions each (b1
written in bits; the stream's own b1 is 1111): 0100 → 1100 (the lowest band turns voiced) 0.376 →
0.065; 1111 → 1011 (band 1 turns unvoiced) 0.556 → 0.121; 1111 → 1101 (band 2 turns unvoiced)
0.383 → 0.212; 1100 → 0100 0.147 → 0.129; 1101 → 1111 (band 2 turns voiced) 0.391 → 0.384;
1011 → 1111 0.222 → 0.222 (to three decimals, 0.22197 → 0.22192; the decodes differ inside the step
frame, in 10, 2 and 12 samples at positions 16, 25 and 34). Two of the 18 get worse: 1101 → 1111 at
position 34, 0.426 → 0.478, and 1011 → 1111 at position 16, in the fifth decimal (0.13696 →
0.13697). Variants on the 18 band steps (scratch builds): bands turning unvoiced only 0.242, bands
turning voiced only 0.293, the frame before's state for every band 0.332, steady bands eroded too
0.433; D-STAR read on P25's grid 0.189, but then the 9 onset steps go from 0.034 to 0.088 and the 9
offset steps from 0.061 to 0.138, where with D-STAR's own columns they stay at 0.034 and 0.061 to
three decimals. The library's decodes of the 99 streams the comparison decodes equal those of the
scratch build the rule was measured on, before and after the review's fixes.

On DVSI's vectors, at joint lags against the commit before (the unvoiced synthesis above),
development / validation / final:

| Codec | Waveform corr. vs DVSI, sample for sample | LSD vs DVSI (dB) |
|---|---|---|
| D-STAR | 0.69 / 0.67 / 0.55 → 0.95 / 0.93 / 0.93 | 2.03 / 1.94 / 2.60 → 1.55 / 1.37 / 1.81 |
| P25 | 0.85 / 0.86 / 0.87 → 0.92 / 0.94 / 0.93 | 2.24 / 1.94 / 2.41 → 1.76 / 1.40 / 1.73 |
| Rate 33 | 0.91 / 0.93 / 0.91 → 0.94 / 0.95 / 0.95 | 1.99 / 1.81 / 2.47 → 1.47 / 1.26 / 1.73 |

Every vector comes closer to DVSI's decode and 0.06-0.11 dB closer to the input speech; under
channel errors the frames neither decoder mutes come 0.22-0.52 dB closer with the mute decisions
unchanged. The `voiced-subframes` gate sets were declared before scoring. The first run, at the
old voiced gain, failed their level check in every mode (the level against DVSI's moved by -0.14
to -0.19 dB, allowed 0.1); the gain was then set on development frame medians and the sets scored
again. They now pass except one check: rate 33's validation vectors gain 0.020 in waveform
correlation where the set asks for 0.03.

By frame class, main and now: the correlation of a frame with DVSI's waveform (median / mean
over eleven development vectors, frames above -56 dBFS in DVSI's output). In the first three
classes the frame and both its neighbours are of the kind named.

| Codec | Build | No voiced band | Every band voiced | Some bands voiced | First voiced frame | First unvoiced frame |
|---|---|---|---|---|---|---|
| D-STAR | main | 0.001 / -0.001 | 0.928 / 0.604 | 0.919 / 0.612 | 0.477 / 0.409 | 0.246 / 0.297 |
| D-STAR | now | 0.999 / 0.946 | 0.966 / 0.949 | 0.964 / 0.944 | 0.955 / 0.921 | 0.956 / 0.903 |
| P25 | main | 0.000 / 0.000 | 0.935 / 0.790 | 0.940 / 0.782 | 0.642 / 0.523 | 0.415 / 0.400 |
| P25 | now | 0.998 / 0.947 | 0.955 / 0.890 | 0.965 / 0.922 | 0.953 / 0.913 | 0.945 / 0.871 |
| Rate 33 | main | 0.004 / 0.004 | 0.959 / 0.866 | 0.938 / 0.849 | 0.297 / 0.330 | 0.330 / 0.344 |
| Rate 33 | now | 0.992 / 0.930 | 0.969 / 0.931 | 0.969 / 0.945 | 0.976 / 0.948 | 0.967 / 0.882 |

The first voiced frame after an unvoiced one played a median 1.8-2.1 dB louder than DVSI's on
main and plays 0.10-0.14 dB louder now; the first unvoiced frame after a voiced one 1.2-1.4 dB
louder on main and within 0.14 dB now. The means of the unvoiced frames sit under their medians
because of the frames whose pitch code or V/UV row tells DVSI's decoder more than is modelled
([Known differences](#known-differences)).

### Voiced onsets and offsets

Where a voiced run starts and ends inside a frame is not always in the model. DVSI's decoders
read it from fields the specifications leave free, and the decoders here now mark a frame with
what those fields say (one mark per frame, in the unused element 0 of `PHIl`):

- **Unvoiced frames have pitch codes of their own.** DVSI's encoders give a frame with no voiced
  band either one default pitch code (P25 b0 25, D-STAR b0 34: 71% of the unvoiced frames of
  the speech vectors in both) or a code from a range they use for nothing else (P25 b0 128..183,
  D-STAR b0 64..119). Rate 33 has its V/UV rows 16 to 31 for this. Each of the 56 also names a
  pattern of pulsed 500 Hz channels for the frame itself, which DVSI's decoders play in place of
  noise there ([D-STAR and P25: the pulsed state](#d-star-and-p25-the-pulsed-state)).
- **Outside those ranges the voicing before lasts to mid-frame.** A P25 or D-STAR frame with no
  voiced band at any other pitch code takes the previous frame's voicing for its mid-frame
  model when that voicing reaches 505 Hz, the pulse train runs on at the old pitch through the
  whole frame, and a voiced frame after it continues the train. On the DVstick 30 (#157) a frame
  whose b1 alone is set to 0 plays exactly as its voiced neighbour would through sample 80, at
  b0 20 and at the pitch code of the frame it was made from, and departs from the frame's first
  sample at b0 100. On DVSI's vectors the default-code frames after a voiced frame correlate
  0.84-0.96 with DVSI's waveform, from 0.58-0.88, and voiced runs after one such frame are
  aligned with DVSI's in 98 of 105 cases for P25 (67 of 98 before) and 83 of 91 for D-STAR (47 of
  85). The 505 Hz comes from P25: three voiced harmonics reaching 502.6 Hz do not hold and three
  reaching 507.9 Hz do. Rate 33's row 16 ends the voicing where the frame starts. Since the
  mid-frame columns change the hold is by column, in D-STAR and P25: the bands of the frame before's
  covered 500 Hz columns (and those within 250 Hz below one) hold through the first 10 ms and the
  others do not, the pulse train runs on through the frame when one holds and starts afresh when
  none does, and the 505 Hz test is
  gone (P25's columns start at 505 Hz, so a frame voiced only below it still holds nothing; [D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns)).
  The measurements of this item stand as made then.
- **D-STAR b0 64..71, 88..103 and 112..119 carry their pitch** into a voiced frame after them:
  its first subframe runs at the geometric mean of the two pitches (253 onsets of 253 on DVSI's
  vectors), and at the voiced frame's own pitch after any other unvoiced frame (925 of 925).
- **Late onsets.** A rate-33 frame coded with V/UV row 15 starts its voiced bands in its second
  subframe; before the mid-frame voicing change the library did so only after an unvoiced frame, and
  with it does so whatever the frame before, in the rate-33 rendering only, beside rows 13 and 17,
  whose first halves have voicings of their own ([Rate 33: rows 13, 15 and
  17](#rate-33-rows-13-15-and-17)). Before the mid-frame
  voicing change the same mark started the voiced bands of a P25 frame whose voicing does not reach
  505 Hz in its second subframe after an unvoiced one (all 40 such onsets on DVSI's vectors); with
  it that case is part of the column rule of [P25: the bands that start with the
  frame](#p25-the-bands-that-start-with-the-frame), and P25 no longer sets the mark.

With these, voiced runs after two or more unvoiced frames start where DVSI's do in 95% of cases for
P25 (791 of 831 on ten vectors), 96% for D-STAR (812 of 846) and 96% for rate 33 (775 of 810).
Measured again for the octave step on all 25 development and validation vectors (D-STAR
voiced runs of three voiced frames or more after a frame with no voiced band, loud; out of line
where the two frames after the first correlate with DVSI's under 0.8), 25 of 1,945 D-STAR runs were
out of line before that step and 4 of 1,945 are with it: most of the misses after a pitch-carrying
code were a step of about an octave into the second subframe ([The pulse train across an
octave](#the-pulse-train-across-an-octave)). P25 onsets that DVSI's decoder starts in the frame's
second subframe and the library started with the frame were open before the mid-frame voicing
change; its
column rule and the onsets that rule still does not fit are in [P25: the bands that start with the
frame](#p25-the-bands-that-start-with-the-frame). Since the mid-frame columns change D-STAR's onsets
are read in columns too, its own (a voiced band starts at mid-frame where neither its column nor the
column 250 Hz higher has both its own sample and the next voiced), and P25's onset rule is one case of a rule for every band of both codecs ([D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns)). See also [Known differences](#known-differences)
for what the unvoiced pitch codes further tell DVSI's decoder, and [Rate 33: the pulsed
rows](#rate-33-the-pulsed-rows) for rate 33's rows 18 to 31.

### Phase regeneration

DVSI's decoders regenerate the phase of a voiced harmonic from the frame's decoded log2 spectral
amplitudes, before spectral enhancement, the high-band shelf and any scaling, with the kernel this
library already used (US 5,701,390 eqs 7-9: an odd 1/m edge-detection kernel over the log
magnitudes). The codec paths applied it to the enhanced amplitudes they render.

**Where the difference was.** Over the 25 development and validation vectors, our 16-bit decode of
DVSI's bits differed from DVSI's decode of the same bits, sample for sample at lag 0, by a residual
energy of 0.115 (D-STAR), 0.155 (P25) and 0.110 (rate 33) of DVSI's (correlation 0.943, 0.923 and
0.945). 96-97% of that residual lies in frames that have a voiced band and follow a frame that has
one, which hold 97-98% of the signal energy. Within voiced runs that line up with DVSI's (run
correlation 0.8 or more) the ratio was 0.09 to 0.15 (a pitch change between the frames under 1%:
0.093, 0.111 and 0.090 for D-STAR, P25 and rate 33; over 10%: 0.106, 0.150 and 0.096): it rose
little with the pitch change in D-STAR and rate 33 and to 0.15 in P25 at large changes, and it
did not depend on the voicing mix (every band voiced in both frames: 0.088, 0.109 and 0.088) or the
pitch period. The difference was in the harmonics of every voiced frame, not in the transitions.

**What differed.** On steady, fully voiced frames of the development vectors (every band voiced,
the pitch within 1% over four frames, DVSI's frame RMS 300 or more: 224 D-STAR, 201 P25 and 238
rate-33 frames), fitting one complex gain per harmonic of DVSI's output against ours, over a Hann
window on the frame, takes the windowed residual from a median 0.075, 0.100 and 0.058 of the signal
to 0.0007, 0.0007 and 0.0006. So within a frame each harmonic of ours was a scaled, rotated copy of
DVSI's; what the fit took out was each harmonic's amplitude and phase: the amplitudes by a standard
deviation of 0.4-0.6 dB in D-STAR and P25, the phases by 0.27, 0.31 and 0.21 rad rms (energy
weighted). One gain over one frame does not show the frequencies exact, and [P25 pitch](#p25-pitch)
finds P25's off by up to 0.07%.

**Which magnitudes the phases come from.** Regressing that phase difference, for D-STAR, on two
free 19-tap kernels, one over the enhanced log2 magnitudes the phases were regenerated from and
one over the decoded ones, gives minus the library's kernel on the first and the library's kernel
on the second: odd taps 1, 3, 5 and 7 of -0.438, -0.135, -0.080 and -0.073 over the enhanced
magnitudes and +0.441, +0.135, +0.081 and +0.074 over the decoded ones, where the library's kernel
has 0.440, 0.147, 0.088 and 0.063, and even taps within 0.015 of zero. After a constant the phase
difference left falls from 0.211 to 0.049 rad rms. DVSI's decoder applies the same kernel to the
decoded amplitudes; P25 and rate 33 give the same signs and sizes, more noisily. The patent words
it so: the computation "is applied to the decoded spectral magnitudes for the current frame".

**The rule.** The codec paths regenerate a voiced harmonic's phase from the frame's `log2Ml`, for
every band whatever its voicing, with the kernel unchanged (its floor, tail, mirror and taps;
`mbe_regenerate_phases_log2()`). `log2Ml` holds the decoded amplitudes in every codec, in the
synthesized state too, a P25 tone frame's included ([P25 tones sent as voice
frames](#p25-tones-sent-as-voice-frames)). The mid-frame model ([Synthesis in
subframes](#synthesis-in-subframes)) takes the mean of the two frames' `log2Ml`, each read at the
band's frequency as the magnitudes are, as stored: where the magnitudes' interpolation takes the
unvoiced factor `0.2046 / sqrt(w0)` out of an unvoiced band, `log2Ml` carries no such factor. Taking
it out or putting it in measured no better (correlation with DVSI's decode, development means at lag
0, D-STAR / P25 / rate 33):

| Phases from | Correlation |
|---|---|
| The enhanced amplitudes (main) | 0.9461 / 0.9182 / 0.9428 |
| `log2Ml` as stored | 0.9844 / 0.9386 / 0.9659 |
| `log2Ml`, the unvoiced factor taken out of unvoiced bands | 0.9832 / 0.9376 / 0.9630 |
| `log2Ml`, the unvoiced factor put into unvoiced bands | 0.9829 / 0.9358 / 0.9658 |

A `log2Ml` from a caller's state is read on its bit pattern: not a number counts as the kernel's
floor, and finite values are kept within ±1000. The public `mbe_synthesizeSpeechf()` keeps the
kernel over the enhanced, smoothed magnitudes it renders.

**The kernel's input domain in P25** (the P25 decoder change set). An AMBE decoder's `log2Ml` of an
unvoiced band is already the band's magnitude with the unvoiced factor taken out: TIA-102.BABA-1's
log-magnitude domain adds `0.5 log2(w0) + 2.289` to an unvoiced band before it is coded (2.289 is
`-log2(0.2046)`), and its decoder builds the unvoiced band as `0.2046 / sqrt(w0)` times
`2^log2Ml`. P25's decoded `log2Ml` of an unvoiced band is the IMBE magnitude itself. DVSI's P25
decoder evidently regenerates P25's phases from the AMBE codecs' domain: for a band the frame
renders as unvoiced, the kernel's input is `log2Ml - log2(0.2046 / sqrt(w0))`, 2.289 + 0.5
log2(w0) more, at the frame's own pitch. Read as decoded, the kernel answers a step at every
voicing boundary that is not there, and every voiced harmonic near one is turned by the kernel's
response to that step. On the library before the change (with the grid and the floored mid-frame
parabola in: [P25 pitch](#p25-pitch)), the phase of DVSI's voiced harmonics against ours (steady
voiced frames in line with DVSI's decode, from a run's third frame, DVSI's decode coherent with ours
at 0.8 or more in the harmonic's band; energy-weighted means; 160,888 harmonics in mixed frames,
199,373 in fully voiced ones):

| Voiced share of the frame | under 0.3 | 0.3-0.5 | 0.5-0.7 | 0.7-0.85 | 0.85-0.999 | fully voiced |
|---|---|---|---|---|---|---|
| before, rad | +0.101 | +0.085 | +0.068 | +0.047 | +0.038 | +0.009 |
| with the change, rad | +0.012 | +0.006 | +0.004 | +0.004 | +0.005 | +0.005 |

and in mixed frames by the distance in harmonics to the first unvoiced band above, before → with
it: 1: +0.309 → +0.002, 2: +0.139 → +0.004, 3: +0.125 → +0.008, 4: +0.104 → +0.003, 5: +0.101 →
+0.005, 6: +0.074 → +0.010, 8: +0.067 → +0.003, 10: +0.063 → +0.001, 13: +0.056 → +0.004, 16:
+0.042 → +0.006, 20: +0.017 → +0.006, and where the only unvoiced bands lie below, -0.092 → -0.032:
the kernel's response to a step on the unvoiced bands, tap for tap (its taps are 1/m at odd m to 19,
the sums from distance d 2.13, 1.13, 1.13, 0.80, 0.80, 0.60 ... and nothing from 20). The step's
size, from a regression of the phase on 0.44 times that response and 0.44 times the response times
log2(w0), before the change: `c = 1.894 + 0.413 log2(w0)` (per pitch range: 0.18, 0.23, 0.43, 0.57
and 0.81 at w0 about 0.047, 0.062, 0.082, 0.109 and 0.177), where the unvoiced factor's
`-log2(0.2046 / sqrt(w0))` is `2.289 + 0.5 log2(w0)`; the regression leaves 0.248 of the phase's
mean square, against 0.328 with a constant step. With the change it finds `0.078 + 0.020 log2(w0)`
(0.995 of the mean square left: nothing). The same regression on rate 33 and D-STAR finds no step
(the study's scratch build; their decoded magnitudes are already in that domain). The decoded
voicing decides, as the synthesis renders it (after adaptive smoothing, which can voice a band: a
band it voices is rendered voiced, its `log2Ml` the decoded magnitude, in the voiced domain, and is
read as decoded); a band with no state and a tone model are read as decoded (a tone model is told by
the structure the decoder gives one, its mark, its word in `Vl[0]`, its two harmonics alone voiced
and every other band of no magnitude, not by its mark alone, since the previous model is caller
state). The kernel's floor and its tail past L take the shifted values: the tail follows band L, and
an unvoiced band under the floor can rise above it. The mid-frame envelope reads each frame's
shifted inputs, at the frame's own pitch, before the two meet at the band's frequency. On the
waveform (whole-file correlation at lag 0, mean over the 20 development / 5 validation vectors),
against the library with every rule in (0.994695 / 0.996755): the step left out -0.002170 /
-0.001278 (every vector down); `K + 0.5 log2(w0)` with K 2.0 -0.000278 / -0.000199, 2.2 -0.000025 /
-0.000027, 1.8 -0.000801 / -0.000526, 2.5 -0.000155 / -0.000060; the regression's `1.905 + 0.449
log2(w0)` (the study's) -0.000175 / -0.000126; half the step -0.000542 / -0.000338; the step on the
bands the own model's column reading leaves unvoiced instead of the decoded ones -0.000276 /
-0.000160 (every vector down). The step alone on the baseline: +0.002747 / +0.001540, every vector
up. The earlier measurement of the table above ("the unvoiced factor taken out of unvoiced bands",
-0.0010 against `log2Ml` as stored on P25's development vectors) was made with the earlier pitch
law's error in place, which hid it.

**D-STAR held histories.** A D-STAR frame with an uncorrectable C0 word plays as a lost frame
while the voice frame decoded from it becomes the prediction history ([D-STAR channel
errors](#d-star-channel-errors)). In the first version of the phase change (before any release), the
state it synthesized carried the history's `log2Ml` beside the repeated model's amplitudes, so the
repeat took its phases from another frame (main regenerated them from the `Ml` it rendered). It now
renders the repeated model with that model's own `log2Ml` and leaves them in `prev_mp_enhanced`, and
so do the repeats after it and the speech under a tone's first frame; the prediction history in
`prev_mp` is unchanged, and the model resumed after the repeat limit's mute carries the history's
`log2Ml`.

A caller may keep one state as `cur_mp` and `prev_mp_enhanced` of the D-STAR decoder. That is
outside the library's contract, which wants three distinct objects ([Stateful Decode
Workflow](../README.md#stateful-decode-workflow)); what follows is what is checked for it, and
nothing else is. There the phase change broke a held history from the second repeat on (found in
review, before any release): the state left as the last synthesized frame, with the played model's
`log2Ml`, was then given the history's back, and the next lost frame's decode overwrote it before it
was read, so such a frame no longer played as a lost frame. The repeated model's `log2Ml` are now
read before the frame is decoded, and every commit to `prev_mp` writes the prediction history
explicitly. Over every sequence of four voice frames followed by zero to six voice, lost,
uncorrectable-C0, single-tone and dual-tone frames (19,531 sequences, none with channel errors),
what the decoders predict from (`log2Ml`, L and w0) is the same with one state as with three after
every frame, and so are the rest of the committed model (K, the voicing, gamma, the frame mark) and
its error state; the amplitudes `Ml`, which no decoder predicts from, differ after a tone that
follows the repeat limit's mute, and the phases and synthesis state stored with them (`PSIl` and
`PHIl` above element 0, `localEnergy`, the overlaps) are, with one state, the synthesized frame's.
With three states the PCM and `prev_mp` are identical to the commit before over the same sequences.
With channel errors the voicing and amplitudes of the model a tone's first frame commits can differ
between the arrangements, as on main: with one state the speech under the tone's first frame leaves
the state as the repeat it synthesized, enhanced and adaptively smoothed with that frame's error
state, and that is committed, where three states commit the last frame's model as that frame
smoothed it (a voice frame at b0 90, b1 0, b2 37 with three C0 and three other corrections, then a
single tone with two and three: `prev_mp->Vl[30]` is 1 with one state and 0 with three, L 34 in
both). Over 244,444 sequences of four voice frames followed by up to five such frames (up to four
where the voice frames have no voiced band), the first four all clean or all with three C0 and three
other corrections and each later frame either, `log2Ml`, L and w0 never differed; the voicing
differed in 8,256 of the 22,222 with unvoiced voice frames. Neither value is established against
DVSI's decoder.

Two behaviours change with one state only, and both differed from the three-state behaviour on main
already. The first voice frame after the repeat limit's mute no longer resumes from the history at
gain 0.3 ([D-STAR lost frames](#d-star-lost-frames)). That wrote the resumed model into the frame
being decoded, which then played the old model and committed its L, w0 and voicing as the prediction
history; the frame now plays and commits its own decoded model, as every other voice frame does with
one state. And the repeats after the repeat limit's mute stay muted: each loaded the repeated model
over the silent last frame the mute left, and the mute faded that speech out again on every frame
(from a fresh state, a voice frame and six lost frames: the 5th and 6th lost frames peaked at 119
and 160 in float PCM, flagged `MUTE`, where three states play the mute noise, 0.107). Now, as with
three states, the 4th repeat fades the last frame out and every repeat after it, until a voice or
tone frame, plays the mute noise alone.

**Result.** On the same steady frames the median frame residual falls from 0.0735 to 0.0040
(D-STAR) and from 0.0597 to 0.0067 (rate 33), and from 0.1065 to 0.0420 for P25, whose remaining
difference is its pitch ([P25 pitch](#p25-pitch)); after one offset and one time shift per frame,
0.005 rad of phase is left per harmonic in all three. By frame class (the table of [Synthesis in
subframes](#synthesis-in-subframes), recomputed over all 20 development vectors: median / mean
correlation of a frame with DVSI's waveform, frames above -56 dBFS in DVSI's output, in the first
three classes the frame and both its neighbours of the kind named; P25 also after [P25
pitch](#p25-pitch), rate 33 also after [Rate 33: mid-frame pitch](#rate-33-mid-frame-pitch)):

| Codec | Build | No voiced band | Every band voiced | Some bands voiced | First voiced frame | First unvoiced frame |
|---|---|---|---|---|---|---|
| D-STAR | main | 0.998 / 0.950 | 0.967 / 0.949 | 0.963 / 0.942 | 0.956 / 0.922 | 0.958 / 0.909 |
| D-STAR | phase | 0.998 / 0.950 | 0.997 / 0.982 | 0.996 / 0.981 | 0.986 / 0.948 | 0.979 / 0.926 |
| P25 | main | 0.998 / 0.951 | 0.957 / 0.905 | 0.963 / 0.920 | 0.954 / 0.916 | 0.946 / 0.881 |
| P25 | phase | 0.998 / 0.951 | 0.990 / 0.927 | 0.987 / 0.941 | 0.978 / 0.935 | 0.960 / 0.891 |
| P25 | P25 pitch | 0.998 / 0.951 | 0.998 / 0.973 | 0.989 / 0.953 | 0.977 / 0.936 | 0.963 / 0.901 |
| Rate 33 | main | 0.992 / 0.934 | 0.970 / 0.934 | 0.968 / 0.942 | 0.976 / 0.950 | 0.966 / 0.886 |
| Rate 33 | phase | 0.992 / 0.934 | 0.992 / 0.955 | 0.993 / 0.970 | 0.993 / 0.966 | 0.981 / 0.895 |
| Rate 33 | mid-frame pitch | 0.992 / 0.934 | 0.993 / 0.965 | 0.994 / 0.972 | 0.993 / 0.966 | 0.982 / 0.898 |

The classes hold, D-STAR / P25 / rate 33, 7511 / 7224 / 6462, 6443 / 4697 / 6477, 4958 / 7072 /
5568, 1689 / 1701 / 1483 and 1496 / 1609 / 1356 frames.

**On hardware.** The DVMEGA DVstick 30 step curves of
[#157](https://github.com/arancormonk/mbelib-neo/issues/157) (the table of [Synthesis in
subframes](#synthesis-in-subframes); mean error against the unit's curves per 16-sample block of
the step frame), main and now:

| Steps | Main | Now |
|---|---|---|
| 30 gain steps (level) | 0.009 | 0.008 |
| 9 from voiced to unvoiced | 0.086 | 0.061 |
| 9 from unvoiced to voiced | 0.094 | 0.034 |
| 3 pitch steps | 0.123 | 0.125 |
| 24 double steps (their level curves 0.203 and 0.203) | 0.341 | 0.341 |
| 18 steps of one band's voicing | 0.342 | 0.346 |

Both columns are measured with the method of the earlier table, Main on main at `7fb5ca6` and Now
on this branch, so a last digit can differ from that table's, measured at 2.6.8. The double-step
row counts 24 of the package's 30 double steps: for the six that step there and back (three gain
steps up and down, three from voiced to unvoiced and back) our decode equals the unstepped
stream's again from the twelfth frame after the step in both builds, so the relative curve, which
is normalised by the settled difference, is undefined. The packages of #155 to #158 hold no audio
(hashes, per-frame levels, curves), so they cannot check harmonic phases directly.

On the scoreboard (joint lags, against main), decoding DVSI's bits, partition means:

| Codec | Partition | Waveform corr. vs DVSI, sample for sample | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|---|
| D-STAR | development | 0.9461 → 0.9844 | 1.461 → 1.357 | 0.9877 → 0.9881 | +0.010 → +0.025 | 6.428 → 6.426 |
| D-STAR | validation | 0.9299 → 0.9612 | 1.372 → 1.248 | 0.9846 → 0.9849 | -0.018 → -0.006 | 4.854 → 4.847 |
| P25 | development | 0.9182 → 0.9386 | 1.759 → 1.668 | 0.9886 → 0.9887 | +0.032 → +0.046 | 6.320 → 6.318 |
| P25 | validation | 0.9422 → 0.9619 | 1.402 → 1.291 | 0.9871 → 0.9872 | -0.025 → -0.014 | 4.734 → 4.724 |
| Rate 33 | development | 0.9428 → 0.9659 | 1.432 → 1.295 | 0.9916 → 0.9919 | +0.000 → +0.014 | 6.263 → 6.271 |
| Rate 33 | validation | 0.9524 → 0.9722 | 1.249 → 1.115 | 0.9880 → 0.9883 | -0.053 → -0.041 | 4.836 → 4.833 |

The waveform correlation rises on every vector but one, P25 `t01` (hum under the speech, long
voiced runs at one pitch code: 0.9476 to 0.9389, and 0.9956 once [P25 pitch](#p25-pitch) follows);
every other P25 vector rises by 0.0072 or more. The distance from DVSI's decode falls on every
vector in every mode, by 0.10 dB (D-STAR), 0.09 (P25) and 0.14 (rate 33) on the development
vectors, and the partition means of the distance from the input move by less than 0.01 dB. The gate
sets `phase-source-development`, `-validation`, `-tones` and `-channel` were declared before
scoring. `-validation` and `-channel` pass; `-development` fails on one check, that the waveform
correlation fall by no more than 0.002 on any vector (P25 `t01`, -0.0087), and passes the others;
`-tones` fails on one check, that the level error of the steady tones rise by no more than 0.05 dB
on any D-STAR or rate-33 vector (D-STAR `knox_12`, +0.070 dB), and passes the others.

- **Tones.** The decodes of 25 D-STAR, all 52 P25 and 10 rate-33 tone vectors change, in the voice
  frames they carry (P25 codes every tone as voice, D-STAR its KNOX tones). No D-STAR or rate-33
  window plays a wrong tone. On D-STAR the detector agrees with DVSI's output more often on `knox_9` and `knox_12`
  and as often on the others, and the level error of the steady tones moves by -0.013 to +0.070 dB
  on a vector (KNOX tones, which D-STAR codes as voice: their one or two strong harmonics sit in long
  runs at one pitch code, where the correlation follows the pitch). On rate 33 the detector loses
  0.07% of `dtmfvbad`'s windows and the level error moves by at most 0.0003 dB. The encoders' tone
  frames are unchanged on every vector.
- **Channel errors.** Every decode of the 8 channel-error vectors changes in every mode; the mute
  decisions do not, at any error rate, hard or soft. On the loud frames neither decoder mutes the
  distance from DVSI's decode falls by 0.092 dB for D-STAR (hard decisions; 0.087 soft), 0.037 dB
  for P25 (0.025) and 0.118 dB for rate 33 (0.113).
- **Encoders.** Every encoder's bits are unchanged.

PESQ and ESTOI ([Method](#method)) of the decode of DVSI's bits, main against the change (mean
paired change ± standard error, vectors improved, of 25):

| Codec | PESQ vs input | ESTOI vs input | PESQ vs DVSI's decode |
|---|---|---|---|
| D-STAR | 2.9476 → 2.9468: -0.0008 ± 0.0020 (11/25) | 0.7643 → 0.7647: +0.0004 ± 0.0001 (22/25) | 4.2696 → 4.2729: +0.0034 ± 0.0004 (19/25) |
| P25 | 3.0399 → 3.0376: -0.0023 ± 0.0033 (13/25) | 0.7868 → 0.7870: +0.0002 ± 0.0002 (21/25) | 4.2699 → 4.2739: +0.0040 ± 0.0012 (21/25) |
| Rate 33 | 2.9615 → 2.9556: -0.0059 ± 0.0030 (10/25) | 0.7654 → 0.7652: -0.0002 ± 0.0005 (17/25) | 4.2327 → 4.2378: +0.0051 ± 0.0009 (24/25) |

The decode comes closer to DVSI's by PESQ, as by the distance from DVSI's decode on every vector and
the waveform correlation on every vector but P25's `t01` (above). Against the input it hardly moves:
D-STAR's ESTOI rises on 22 of 25 vectors, by +0.0004 ± 0.0001 (four standard errors), P25's and rate
33's ESTOI move by no more than a standard error, and PESQ against the input moves by less than a
standard error for D-STAR and P25 and falls by two for rate 33. The change is judged by the waveform
against DVSI's decoder, the project's rule where the specification is silent ([Departures from the
specification](#departures-from-the-specification)), not by a gain against the input.

### P25 pitch

DVSI's P25 decoder plays pitch code b0 at TIA-102.BABA eq 46's fundamental cut to a grid of 2^-19
cycle per sample, `f0 = floor(2^21 / (2 b0 + 79)) / 2^19` (eq 46's `4 pi / (b0 + 39.5)` radians per
sample is `4 / (2 b0 + 79)` cycles), from 0 to 2^-19 cycle, a relative 0 to 2.26e-4 (b0 201),
under eq 46's: the value the library already played for a tone its encoder sent as voice frames
([P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)), found there first. Between two
models that both have a voiced band, its mid-frame pitch is not the geometric mean of the two but a
parabola for that square root, floored to the same grid: with `fa fb = y 2^e` the product of the two
fundamentals in cycles per sample, y in [1, 2),

    q(y) = sqrt(1.5) + (sqrt(2) - 1) (y - 1.5) + 2 (sqrt(2) + 1 - 2 sqrt(1.5)) (y - 1.5)^2,

the parabola through sqrt(y) at y = 1, 1.5 and 2, the mid-frame fundamental is
`floor(2^19 sqrt(fa fb) q(y) / sqrt(y)) / 2^19`. Its relative error `q(y) / sqrt(y) - 1` is 0 at
both ends of the octave and at y = 1.5, about +7.7e-4 near y = 1.2 and -1.2e-3 near y = 1.75: the
half-octave wave the library had fitted into each code's own pitch, at twice that law's size and in
the place it measures. It applies at an unchanged code as at a step; at an onset or an offset (one
of the two models with no voiced band) the voiced model's own pitch stands. Nothing in the rule is
fitted: the grid is the tone path's, the nodes are the interpolant's. What the decoder computes is
not known beyond its output (a three-point parabola for a square root, floored to a 2^-19 grid, is
the simplest thing that gives the mid-frame pitch measured; the integer the decoder holds may differ
from the floor of the exact parabola by a unit at some pairs, which the decodes cannot show).
L, K and the bit layout stay on eq 46's nominal pitch, computed in float as the decoder always has.
ProVoice, which reaches the same data path, has no DVSI reference: it keeps eq 46's nominal pitch
and the geometric mean at mid-frame, a decision, not a finding. The public
`mbe_decodeImbe4400Parms()` stores the grid pitch.

**What could not tell them apart.** The library's earlier law (below) was fitted to the drift of
DVSI's pulse train over runs of four frames at one pitch code. Over such a run the pulses move at a
mean of the frame's own pitch and its mid-frame pitch, so the drift cannot tell the two apart: an
offset fitted there is half the mid-frame pitch's and half the own pitch's. What tells them apart
is the lag of DVSI's decode against ours read at each frame boundary, where both play the frame's
own model (from the phases of the voiced harmonics: [Rate 33: mid-frame pitch](#rate-33-mid-frame-pitch)
has the method), and its increment over the frame, which the mid-frame pitch alone sets (a
mid-frame pitch off by a relative d moves the train by 80 d samples over the frame). And one rule's
error hid the next: the phase kernel's reading of P25's unvoiced bands ([Phase
regeneration](#phase-regeneration)) left an error in mixed frames and at onsets that biased every
drift and lag measured before it was found, so the figures below are the library's with every rule
before it in, re-measured on the change set's library (its scratch-switch copy, outside the
repository, decoding byte for byte as the library with each rule on and as the baseline with all
four off).

**The grid.** The per-code drift (DVSI's fundamental against eq 46's over runs at one code: 49
rows of pitch code and corpus with 8 frame pairs or more, from DVSI's speech and pulse-train
vectors; a measurement of DVSI's decoder, which no change of the library's moves), after taking out
a candidate quantization's relative error and fitting a constant and two harmonics over the half
octave to the rest, weighted by each row's standard error (floored at 0.0005%): 0.00338% rms left
for f0 floored to 2^-19, 0.00447% for 2^-20 floored, 0.00452% for 2^-21 rounded, 0.00453% for no
quantization and for every grid of the period (2^-6 to 2^-10 of a sample, floored or rounded) or of
w0 tried. With the four rules in, on the whole-file waveform correlation with DVSI's decode at lag 0
(mean over the 20 development / 5 validation vectors, against the library's 0.994695 / 0.996755):
the own pitch rounded to 2^-19 -0.002065 / -0.001849, no grid (eq 46 as it stands) -0.001948 /
-0.001838, floored to 2^-20 -0.000567 / -0.000560, floored to 2^-18 -0.001646 / -0.001570; the own
pitch a tenth of a grid step higher -0.000037 / -0.000017, a tenth lower -0.000119 / -0.000157, a
quarter step lower -0.000343 / -0.000329. Every vector falls with each but the tenth-step ones.

**The mid-frame pitch.** With the grid and the geometric mean at mid-frame (and the low band, the
pulsed state and every other rule of the branch below), over the 25 vectors, 19,627 frames whose two
boundaries both read a lag, 16,584 with nine harmonics or more at both and an increment under 0.4
samples: the median absolute increment is 0.0630 samples, 0.0623 at an unchanged code (2,121 frames
of the development vectors) and 0.0629 at a step of one or two codes: the deviation is there whether
the code changes or not. With the phase kernel's reading of unvoiced bands in as well (the geometric
mean kept) it is 0.0522 (0.0571 at an unchanged code), and as a function of `u = frac(log2(Pa Pb))`
(the two periods in samples) one curve on the development and validation vectors. A parabola for
the square root over the octave of the fundamentals' product, three coefficients fitted to those
increments of the development vectors, is `1.224750 + 0.414032 (y - 1.5) - 0.071392 (y - 1.5)^2`;
the interpolant at 1, 1.5 and 2 is `1.224745 + 0.414214 (y - 1.5) - 0.070552 (y - 1.5)^2`. The
increments' standard deviation it leaves is 0.03338 samples on development and 0.04390 on
validation, against 0.03337 and 0.04389 for two harmonics in u; the same parabola over the product
of the two periods leaves more (0.03363, 0.04423). With the rule in (the parabola floored), the
median absolute increment falls to 0.0217 samples before the phase change and 0.0117 with it,
0.0116 with all four rules:

| Library (development vectors) | Median absolute increment, samples: unchanged code / step of 1-2 / 3-6 / 7 or more codes |
|---|---|
| baseline (`f725a66`: the old law, the geometric mean) | 0.0156 / 0.0195 / 0.0277 / 0.0336 |
| the grid, the geometric mean | 0.0623 / 0.0629 / 0.0609 / 0.0629 |
| the grid, the floored parabola | 0.0144 / 0.0184 / 0.0264 / 0.0329 |
| ... and the phase kernel's reading of unvoiced bands | 0.0086 / 0.0101 / 0.0134 / 0.0195 |
| ... and the own models by columns (the library) | 0.0086 / 0.0101 / 0.0132 / 0.0185 |

(2,121 to 2,211 frames at an unchanged code, 5,477 to 5,757, 4,043 to 4,244 and 1,453 to 1,534 at
the steps, by build; on the validation vectors the library's medians are 0.0094 / 0.0106 / 0.0148 /
0.0204.) Before the phase
change, frames at a pitch step kept more of the lag's change than frames at one code, increasingly
with the step (the study's scratch build found a fitted parabola drifting from the exact one as the
step grew); with it, the excess of a step over an unchanged code falls to between 0.4 and 0.55 of
what it was (0.0040 → 0.0015 samples at a step of one or two codes, 0.0185 → 0.0099 at seven or
more), and what is left is not explained. The floor against the other
readings, with all four rules in: the mid-frame value rounded to the grid -0.000513 / -0.000675 of
correlation (19 of 20 and 5 of 5 vectors down), left unfloored -0.000512 / -0.000608 (19 and 5
down). The parabola at an onset too (this frame's pitch squared) -0.002205 / -0.001308 (all down), at
an offset -0.000054 / -0.000012: the rule's condition, two voiced models, stands. Each rule alone on
the baseline: the grid with the old law's mid-frame pitch (the geometric mean) -0.023848 / -0.021519
(every vector down: the old law carried the mid-frame deviation in the own pitch), the parabola with
the old law's own pitch -0.033616 / -0.025837; the two together +0.002570 / +0.001383, every vector
up.

**On the pulse-train vectors** (`ptrain0` to `ptrain12`: long runs at one code, where the earlier
law was fitted), the median absolute increment is 0.0163 samples with the library, against 0.0185
with the baseline and 0.0990 with the grid and the geometric mean (0.0179, 0.0176 and 0.0923 at an
unchanged code). The whole-file correlation at lag 0 of those vectors, which over 620 to 2,000
frames at one code mostly shows how far the train has drifted by the time each run ends, rises on 8
and falls on 5 (the study's `parts8.py`: `ptrain1` 0.445 → 0.973, `ptrain4` 0.126 → 0.827, `ptrain7`
0.947 → 0.224, `ptrain12` 0.847 → 0.190): those swings are the phase at a run's end, not the
pitch's error per frame, which falls.

**The floor's edge.** The library forms the mid-frame pitch in double from the stored floats
(`imbe_mid_pitch()`, `src/core/mbe_subframe.c`): their product is exact in double, divided by
(2 pi)^2 for the fundamentals' product, its mantissa by `frexp()`, the square roots from libm in
double, the floor in double, and 2 pi times the whole number of steps over 2^19 rounded to float
once. Over every ordered pair of the 208 codes (equal codes included) the unfloored number of steps
lies at least 6.3e-6 of a step from a whole number (b0 8 and 100, 12866.0000063; next b0 14 and 20,
2.4e-5), about 4.9e-10 of itself, so no build's double arithmetic (x87's extended one, fast-math's
reordering) floors it otherwise; exact arithmetic (decimal, 60 digits) and the library's double
arithmetic floor every pair alike; no product of two codes' fundamentals is an exact power of two
(where q's two ends meet); and 2 pi times any whole number of steps from 3,000 to 27,000 lies at
least 7.6e-13 of itself from a float's rounding edge. At equal codes the parabola keeps the code's
own pitch only where the mantissa of its square lies near 1, 1.5 or 2 (the nearest whole number of
steps there is 0.0038 away, b0 84). A product or a result past the float range, and a mean that
floors to no step, give no pitch, which the caller rejects as every pitch that is not positive and
normal.

**The own pitch's arithmetic.** `imbe_code_pitch()` stores, for P25, `2 pi floor(2^21 / (2 b0 +
79)) / 2^19` formed in double from the integers and rounded to float once (`imbe_grid_fundamental()`):
over the 208 codes that product lies at least 4.9e-11 of itself from a float's rounding edge, so
every build stores the same float. The decode, the tone check (`imbe4400_take_tone()`, which asks a
tone word's model for the pitch of its code by bit pattern) and the pulsed state's code tests
(`mbe_imbe4400_pulsed_code()`, `mbe_imbe4400_pitch_code()`) share that one function, and the tone
model keeps the pitch of its code: one definition. Every decoded pitch is its own code's for the
code tests, the floats next to it on either side are no code's, and eq 46's nominal pitch of any
code (ProVoice's, the repeat limit's reset model and the initial state at b0 134) is no code's, so
the validity counter leaves it (`test_p25_code_tests_on_the_grid`). The quotient `4 pi / w0 - 39.5`
the code tests round lies from 0 to 0.054 above its code (b0 201; 0.046 at most over b0 128..183).

**Decisions read at the pitch played.** The harmonic at or next above 505 Hz (`ceil(505 Hz / f0)`)
is the nominal pitch's for every code (the nearest harmonic to 505 Hz is b0 119's 5th, 0.062% below
it; 0.053% at the nominal pitch). Of the 208 x 14 pairs of code and column edge (505 + 500 k Hz, k 0
to 6, and the same edges 250 Hz lower), the harmonic at or next above the edge at the grid pitch
differs from the nominal pitch's at 2 (the old law's: 14), and every edge lies at least 5.97e-5 of a
harmonic from a whole number (b0 152 at 2005 Hz; 2.86e-4 under the old law, 1.55e-4 at ProVoice's
nominal pitch), so float rounding moves none (`test_p25_onset_edges`). The band edges' transform
bins, the mid-frame harmonic count and the harmonic where the high-band shelf starts follow the
pitch played, as before. The IMBE encoder quantizes the analysis pitch against eq 45 and derives L
and K from eq 46's nominal pitch (`imbe_enc_pitch()`, `imbe_enc_set_pitch()`): no decision of it
reads the decoded `w0`, which moves only in the `cur_mp` it returns; its bits are unchanged on all
95 inputs.

**The earlier law, as history.** From the P25 pitch change (`470374c`) to the branch below, the
decoder stored eq 46's nominal float times `1 + g`, `g = a0 + a1 cos t + b1 sin t + a2 cos 2t + b2
sin 2t`, t = 4 pi log2 of the nominal period (a0 = -1.980564e-4, a1 = 1.263798e-4, b1 = 4.835086e-4,
a2 = 2.352386e-5, b2 = 7.623621e-5; -7.06e-4 at b0 127 to +3.38e-4 at b0 96), the product in double
and rounded to float once: a least-squares fit to that per-code drift (69 rows of pitch code and
corpus, 50 from DVSI's pulse-train vectors and 19 from speech, 5,401 frame pairs on 61 codes, each
row's weight capped at 40 pairs; 0.0064% rms left with those weights). It took the aligned voiced
runs' residual from growing with the run's length (0.023 of the signal over a run's first four
frames to 0.21 past the twentieth) to flat (0.019 to 0.028), and, on the scoreboard, P25's waveform
correlation with DVSI's decode from 0.9386 to 0.9672 (development) and 0.9619 to 0.9855 (validation).
What it fitted was the mean of a frame's own pitch and its mid-frame pitch over runs at one code:
the grid and half the mid-frame parabola's deviation. The tone frames, whose own model stands at
mid-frame and so never take a mid-frame pitch, drifted under it by up to 145 mrad a frame on a
harmonic, which led to the grid ([P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)).

**What is not known.** What the routines are in DVSI's decoder beyond their output; why frames at
pitch steps keep a little more of the lag's change than frames at one code; and what an AMBE-3000
does with pitch codes no DVSI vector holds in three consecutive frames (the grid covers them by
rule).

### P25 tones sent as voice frames

TIA-102.BABA has no tone frames: every P25 voice frame plays its decoded harmonics. DVSI's P25
encoder sends a tone it detects as voice frames built from TIA-102.BABA-1 Annex J: the pitch code
`b0 = floor(4 pi / w0 - 39)` of the tone's fundamental, the tone's harmonics l1 and l2 large, and
a voicing word b1 that is the same in every steady frame of a tone (154 tones, 89 distinct pairs of b0 and b1; index 128: b0 164,
`000011110100`; index 129: b0 53, `1110000`). The fundamental is, for a dual tone, Annex J's
printed f0, with DVSI's own 150.80 and 86.47 Hz for indices 145 and 147 ([Tones as harmonics of a
fundamental](#tones-as-harmonics-of-a-fundamental)); for a single tone it is `31.25 ID / l`, l
Annex J's harmonic number (1 for IDs 5..12, 2 for 13..25, 3 for 26..38, 4 for 39..51, 5 for 52..64,
6 for 65..76, 7 for 77..89, 8 for 90..102, 9 for 103..115, 10 for 116..122; these are the table's
l1 = l2, not the rule `ID / 16 + 1` by which DVSI's rate-33 and D-STAR decoders play a single
tone). For two ranges Annex J prints a multiplier in place of the quotient: `10.417 ID` for IDs
26..38, which is 31.25 / 3 rounded, and `5.2803 ID` for 65..76, where 31.25 / 6 is 5.2083. Built
from the printed values, the codes of nine tones would come out otherwise (IDs 32, 65, 66, 68, 69,
71, 73, 74 and 76). Built so, the codes equal those observed in
DVSI's bit streams for all 154 tones. The voicing word is the same in every steady frame of a tone
(all 12,555 steady frames compared across the codecs agree); 112 frames at tone boundaries carry
another pair. The word is not the bands that hold l1 and l2 under the decoder's band mapping: for
index 128 the band of harmonic 12 is unvoiced.

**The test.** US 7,970,606 describes how a decoder recognises such frames (FIG. 7, step 720): the
reconstructed parameters are checked against its Table 1, which is Annex J's, and a frame is a tone
frame if its fundamental is about an entry's, the voicing bands of the entry's harmonics are voiced
and those harmonics are "dominant over the other spectral magnitudes"; the other harmonics are then
attenuated ("20 dB attenuation is typical"). The patent gives no numbers for the test. Measured on
DVSI's P25 decodes of its tone vectors and of the 25 development and validation speech vectors:

- **Which tones.** The frames selected, by a local script: frames of the P25 tone vectors other
  than a file's first and last, at a dual tone's pitch code, whose decoded magnitudes have a norm
  (the square root of the sum of their squares) of at least 60, whose second half in DVSI's decode
  is a clean pair of components at the tone's two harmonics of the grid pitch (two components leave
  under 5% of its energy), with both tone harmonics decoded voiced and their decoded magnitudes
  differing by more than 0.2 in natural-log ratio. DVSI's decoder plays the two components at equal
  amplitude (the absolute natural-log ratio of the two fitted amplitudes, each over its synthesis
  droop, under 0.04) in 5,092 of 5,094 DTMF frames and 1,344 of 1,344 KNOX frames: those are tone
  frames. Of the 50 call-progress frames (160..163) none meets that test, and 25 lie within 0.1 in
  natural-log ratio of the decoded amplitude ratio: these frames are not played as equal-amplitude
  tone models, and the library leaves them as decoded (they are not in the test). Without the
  condition that both tone harmonics be decoded voiced the counts are 5,482 of 5,484, 1,381 of
  1,381 and 0 of 76 (51 of the 76 within 0.1 of the decoded ratio).
- **The voicing word is needed.** Without it nothing separates a tone from speech: 8 speech frames
  have a dual-tone row's harmonics as their two largest magnitudes (up to 84% of the energy), and
  DVSI's decoder plays all 8 as speech.
- **Dual tones.** With b0 and b1 those of a dual tone there are 13,653 loud frames of the tone
  vectors, 13,557 of them clean pairs, and 48 pairs of a frame and a row in the 25 speech vectors.
  The ratio r of the rms of the other decoded magnitudes to the rms of the tone's is at most 0.296
  on those tone frames (99.9% under 0.242) and at least 0.482 on those speech frames. The library's
  limit is 0.375 (`IMBE_TONE_SHARE_DUAL`). A tone's first frame passes with as little as a third of
  the frame's energy in the tone's harmonics.
- **Single tones.** Their pairs of b0 and b1 also hold speech: 292 frames of the 25 speech vectors,
  high voices with a dominant fundamental. Of the 142 speech frames that pass a single-tone row at
  0.375, the isolated ones, ordered by r: at 0.0159, 0.0206, 0.0210, 0.0213 and 0.0234 the tone path
  is closer to DVSI's output than the speech path, or as close; from 0.0236 up (about 65 frames) the
  speech path is closer in nearly every case (one exception, at 0.0289): no clean split, and the
  limit lies between the two groups. The library's limit is 0.0235 (`IMBE_TONE_SHARE_SINGLE`);
  98.9% of the 918 loud single-tone frames of the tone vectors lie under it (median 0.003).

The library's table (`imbe_tone_rows` in `src/imbe/imbe7200x4400.c`) holds the 116 distinct rows
of b0, b1, l1 and l2 of the DTMF, KNOX and single tones as DVSI's encoder sends them: the harmonics
are Annex J's, the pitch codes are built from its fundamentals as above, and the voicing words were
read from DVSI's bit streams of its tone vectors. Left out are the call-progress tones and three single tones (120..122, at b0 2 and 3)
whose harmonic 10 lies past the frame's L of 9. Single tones on one code and word share a row's pair
at several harmonics; the rows of one pair are all dual or all single tones, and of the rows a frame
passes the one with the smallest r wins. The test reads the decoded magnitudes, before spectral
enhancement. It runs for P25 only: ProVoice reaches the same data path but has no DVSI reference.

**What a tone frame plays**, measured on DVSI's P25 decodes of the tone vectors:

- **Pitch:** eq 46's fundamental cut to the 2^-19 grid, `floor(2^21 / (2 b0 + 79)) / 2^19` cycle
  per sample, not the speech pitch of [P25 pitch](#p25-pitch). On the 14 codes with a long run in
  the KNOX vectors, DVSI's relative offset from eq 46 against the grid's, in units of 1e-5. DVSI's
  column is the mean of the two tone lines' relative offsets from their harmonic of eq 46's
  fundamental, each line's frequency where a Hann-windowed transform of the run's interior peaks,
  over the code's long run (codes 97 and 196 have two each, both counted):

  | b0 | DVSI | Grid | b0 | DVSI | Grid |
  |---|---|---|---|---|---|
  | 56 | -7.77 | -7.77 | 133 | -11.49 | -11.54 |
  | 67 | -7.90 | -7.96 | 146 | -12.35 | -12.40 |
  | 97 | -11.40 | -11.40 | 151 | -5.93 | -6.10 |
  | 112 | -4.26 | -4.24 | 175 | -9.61 | -9.54 |
  | 117 | -2.42 | -2.48 | 195 | -12.09 | -12.06 |
  | 118 | -9.33 | -9.39 | 196 | -12.41 | -12.40 |
  | 128 | -2.53 | -2.48 | 197 | -16.39 | -16.36 |

  The two lines of a tone differ by at most 0.23 (code 151), and the mean lies within 0.18 of the
  grid's offset (code 151; within 0.08 for the other 13).

- **Amplitude:** the tone's harmonics at one amplitude. On DVSI's decoder it is a median 2.0974
  times the rms of their decoded magnitudes on the 3,574 steady frames of the DTMF and KNOX tones
  (2.0975 with the 143 call-progress frames, which the tone path leaves out: 3,717 frames, robust
  spread 1.8%) and 2.108 to 2.119 times for the three single tones of `alert`, where the library
  renders a voiced harmonic at 2.1308 times its magnitude; each harmonic then times
  `sinc^2(pi l / 256)`, the droop of the synthesis table ([Tones as harmonics of a
  fundamental](#tones-as-harmonics-of-a-fundamental); the upper component over the lower is a median
  0.9999 of the sinc^2 ratio). The library takes 0.988 times the rms (`IMBE_TONE_LEVEL`): the steady
  tone frames of `knox_1` to `knox_16` and `dtone_1` to `dtone_16` then play 0.012 dB above DVSI's
  (0.9866 would be exact there; at 0.98 they would play -0.058 dB, at 1.00 +0.117 dB), and those
  of the DTMF, `alert`, `alltone` and call-progress vectors 0.022 dB (dual tones, 2,087 frames)
  and 0.031 dB (single tones, 411 frames) below it.
- **The other harmonics:** two components leave 0.0009 of the frame where the tone holds 60-70% of
  the decoded energy, so the others are at least 27 dB down; the measurement bounds them, it does
  not show them at 0. The library makes them unvoiced at 0, its choice, as clause 8 of
  TIA-102.BABA-1 builds a tone's voice parameters.
- **Timing:** the frame's own model stands at mid-frame. The amplitude moves linearly from the last
  frame's to this one's over the first 80 samples and holds over the second 80; a first tone frame
  plays at half level, then full, and the frame after the last tone frame plays the tone fading out
  over 80 samples.
- **Phases:** the regenerated ones ([Phase regeneration](#phase-regeneration)). The library's kernel
  over the frame's decoded log2 magnitudes, on a pulse train started at the first tone frame's first
  sample, against DVSI's phase from frame to frame: over the 64 components of `knox_1` to `knox_16`
  and `dtone_1` to `dtone_16`, the median of each component's standard deviation of the frame-to-frame
  phase error is 0.0009 cycle (0.018 cycle with zero phases), and the error's mean, an offset that
  stays, is 0.00 to 0.07 cycle by component (median magnitude 0.015), largest for harmonics 4 and 7
  at short periods. A model with zero phases, as rate 33's tone frames have, does not fit:
  its correlations with DVSI's waveform run from -0.5 to +0.85 by vector, at levels within 0.1 dB.
  At mid-frame the frame's own phases follow DVSI's waveform at a median correlation of 0.9924 on
  the 3,840 steady frames of those 32 vectors, the mean of the two frames' `log2Ml` at 0.9903.
- **The noise generator** runs on through a P25 tone frame as through any voice frame, unlike the
  AMBE tone frames: the unvoiced frames two and three after a tone run of `dtmf` and `alltone`
  follow DVSI's waveform at a median correlation of 0.80 and 0.99 (154 and 7 frames) so, and at
  -0.01 and -0.10 with the generator held (`dtmf60` and `dtmf180` hold no such frames).

**How the library carries it.** The decode (`imbe_tone_word()`, P25's path only) leaves the tone it
recognises in the otherwise unused `Vl[0]` of `cur_mp`: a flag with b0 in bits 16-23, l2 in bits
8-15 and l1 in bits 0-7, and 0 on any other frame. The decoded model and its frame mark are left as
they are, so the prediction history is the decoded frame; the public `mbe_decodeImbe4400Parms()`
sets the word too, so the IMBE encoder's state carries it, and the encoder's bits are unchanged. The
synthesis (`imbe4400_synthesize_frame()`) builds the tone model after the spectral enhancement: the
grid pitch (in double, rounded to float once), the tone's harmonics voiced at their level, every
other harmonic unvoiced at 0, L, K and `log2Ml` kept (the phases are regenerated from them), and the
mark `MBE_FRAME_MARK_TONE`, which makes the mid-frame model the frame's own
(`mbe_subframe_mid_model()`). A repeat loads the history, `Vl[0]` included, and plays the tone model
again; the repeat limit's reset clears it. A word counts only if the decode's own test, run again
on the model it sits in, gives exactly that word: the model at the pitch its code decodes to (w0 by
its bits, L and K) and `imbe_tone_word()` on the model's voicing and decoded magnitudes, before the
enhancement, equal to the word. A decode's own word passes on its frame and on every repeat of it;
any other word a caller's state holds plays as no tone (code 1's word in a code-2 tone frame's
history, of the same L, K and harmonics, included). Only the synthesis of a frame that plays the tone model
marks it `MBE_FRAME_MARK_TONE`: no decoder sets the mark, so the data path clears one it finds on
any other frame, P25 and ProVoice alike. The public decoder leaves 0 in `Vl[0]` on every frame it
decodes that is not a tone, special and invalid frames included; a status error
(`MBE_STATUS_INVALID_ARGUMENT`, `MBE_STATUS_INVALID_BITS`) is no frame and writes nothing, so
`cur_mp` keeps what it held. It does not write the word of `prev_mp` (it writes the prediction taps
there, as before). ProVoice's states keep 0 in `Vl[0]`: a state moved between the P25 and ProVoice
entry points in mid-stream is outside the stream contract, and a P25 repeat after a ProVoice frame
plays the repeated model as speech. No result flag marks the frame: it stays a voice frame.

**Since the P25 decoder change set** every P25 code plays the grid pitch ([P25 pitch](#p25-pitch)),
the tone's pitch included: the tone model keeps the pitch of its code, which the decode stored and
the tone check (`imbe4400_take_tone()`) compares by bit pattern with the code's, and the tone path
has no computation of its own any more. Every tone the P25 path plays passes its own check and plays
at the same `w0` as before, by bit pattern (`test_p25_tones`). Tone frames never drifted with the grid
because a tone model stands at mid-frame itself (no mid-frame pitch enters it), which is also why
the grid was found on them first: under the earlier law their drift was the own pitch's alone. All
52 tone decodes change (the frames around and between the tones), and on the eight DTMF vectors
whose tones are separated by short runs of voiced frames (`dtmf4` to `dtmf100`, not `dtmf60` and
`dtmf180`) the whole-file correlation with DVSI's decode at lag 0 falls, from 0.976-0.979 to
0.958-0.966: the voice frame after a tone model takes the floored parabola between the tone's pitch
and its own at mid-frame, which moves the pulse train's phase into the next tone run (the change is
rule 2's alone: with the grid and the geometric mean `dtmf4` is at 0.9777). A geometric mean
next to a tone model (a scratch form) brings six of them back to 0.977-0.978 and takes `dtmf50` and
`dtmf100` down to 0.920, `alert` and `dtmf` lower: no rule found there, and the library keeps the
parabola between any two voiced models. Elsewhere the tones' correlation rises or holds (`cp1`
0.406 → 0.646, `cpvbad` 0.976 → 0.999, `alert` 0.567 → 0.600; `cp31` falls, 0.793 → 0.708, a
call-progress vector P25 plays as decoded speech).

**The first implementation** (`d245efd`) used one limit, 0.375, for single and dual tones. P25
speech then moved away from DVSI's decode: its waveform correlation with DVSI's fell by 0.0044 on
the development vectors (16 of 20 lower, `noisy` by 0.058) and by 0.0072 on the validation ones,
all of it from the 142 speech frames that passed single-tone rows and the frames after them in the
same voiced stretch (1,104 changed frames, their residual 0.0017 to 0.0140 of the signal energy),
and its development and validation gate sets failed ([November 2026: tones and the octave
step](#november-2026-tones-and-the-octave-step)). The single-tone limit above, set on the speech
frames DVSI's decoder plays as speech or as tones, replaced it (`b58d91d`).

**Before this change.** The library (`fac0111`) played every P25 frame as speech.
Over the 52 tone vectors decoded as
P25, frame for frame at lag 0 (a frame is loud when DVSI's frame RMS is above 30; 17,334 steady
frames, loud with every neighbour loud: 17,286 between two loud frames and 48 a file's first or last
frame, with one neighbour only; 791 first loud frames of a tone; 743 first quiet frames after one),
against DVSI's decode (level: mean absolute difference, dB):

| Build | Steady frames: level | Steady frames: correlation with DVSI's waveform, median | First frame of a tone: level | Last frame: level | Frame after: level |
|---|---|---|---|---|---|
| main (`7fb5ca6`, 2.6.12) | 0.555 | 0.847 | 4.27 | 0.46 | 4.82 |
| phase change (`468e582`) | 0.555 | 0.880 | 4.27 | 0.44 | 4.81 |
| P25 pitch (`470374c`) | 0.560 | 0.656 | 4.27 | 0.45 | 4.82 |
| tone model | 0.279 | 0.986 | 0.37 | 0.21 | 3.60 |

(The mean of the steady frames' correlation: 0.582, 0.590, 0.408 and 0.940; the mean over the
vectors of the whole-file correlation: 0.576, 0.550, 0.127 and 0.918.) The row "tone model" is this
change set, measured on the same frames with the same calculation: the like-for-like comparison with
the
rows above. The measurement is the one made at 2.6.8, which
main reproduces. On 2.6.7, before the subframe synthesis, it gave 1.3 dB for the first frame of a
tone (its faster onset is 2 dB too loud on speech), 0.62 dB and a median correlation of 0.76 for the
steady frames, 1.24 dB for a tone's last frame and 6.7 dB for the frame after it. DVSI's decoder
brings such a tone in at its full level over the first 80 samples of its first frame (0.13, 0.40,
0.61, 0.90 of the tone's level in 20-sample blocks, 65 onsets of `dtmf`), and plays the unvoiced
frame before it, whose amplitudes already hold part of the tone, much quieter than the library did
before this change; there
the first frame followed the mid-frame model like any voiced onset (0.01, 0.03, 0.04, 0.04, then
0.18, 0.47, 0.71, 0.95). The measurement of [November 2026: tones and the octave
step](#november-2026-tones-and-the-octave-step) selects other frames (those the tone test passes,
between two such frames, with DVSI's frame above a level); on that selection the median correlation
rises from 0.388 to 0.986, which is not the comparison with the rows above.

### P25: the bands that start with the frame

After a frame with no voiced band, DVSI's P25 decoder does not start every voiced band of the next
frame with that frame. A rule fits what it does: read the frame's voicing in fixed 500 Hz columns,
start the voiced harmonics that lie in a column the voicing covers, or within 250 Hz below one,
with the frame, and start the other voiced harmonics at mid-frame. Of the 1,843 onsets labelled
decisively below, its whole-frame form agrees with 1,838 (70 of the 72 late ones, 1,768 of the
1,771 early ones), and per harmonic it agrees on 0.9565 of the harmonics selected; how DVSI's
decoder comes to it is not known. In TIA-102.BABA a frame's V/UV decisions hold for
the whole frame, which it renders in one piece; the library's synthesis in two subframes, DVSI's
([Synthesis in subframes](#synthesis-in-subframes)), needs a voicing for the mid-frame model, which
the text does not define. The rule below is fitted to DVSI's decodes, and it is listed with the
other departures that follow DVSI ([Departures from the
specification](#departures-from-the-specification), "Where voicing starts and ends inside a
frame").

**Since the mid-frame columns change.** The onset rule of this section is one case of the rule of
[D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns): at a voiced frame after a frame with no voiced band that rule voices a
band at mid-frame exactly where this one does, but for one exception (an isolated voiced last band
of a frame of more than 36 harmonics reads as unvoiced in the columns). The rule reads both frames
in every situation, the hold through a frame lacking a pitch is by column, and
`MBE_FRAME_MARK_ONSET_COLUMNS` no longer chooses bands (below). The measurements of this section stand as they were
made, on the builds they name.

**What was open.** Before this change (up to 2.6.12, and in this series up to `adecb29`; the
measurements below were made on `cc272ef`, the tone model change before its review's fixes, whose
later commits change no decode of a stream a decoder produces) the mid-frame model took this frame's
voicing, so a voiced frame after a frame with no voiced band started its voiced bands, and its pulse
train, with the frame. One exception was known for P25: a frame voiced only below 505 Hz started
them at mid-frame (`MBE_FRAME_MARK_LATE_ONSET`, set in `imbe_decode_voicing()`; [Voiced onsets and
offsets](#voiced-onsets-and-offsets)). The P25 voiced runs that DVSI's decoder starts half a frame
later than the library were an open item of [Known differences](#known-differences), and P25's
voiced runs out of line with DVSI's held the largest single part of the residual left in the three
codecs (the accounting is in [November 2026: mid-frame voicing](#november-2026-mid-frame-voicing)).

**Which onsets start late.** Two test builds outside the repository, of `cc272ef` with environment
switches in `src/core/mbe_subframe.c` (with the switches off their P25 decodes of the 25 vectors
equal that code's byte for byte), decode every onset late (after a frame with no voiced band the
mid-frame model has no voiced band) and none late. An onset is a voiced frame after a frame with no
voiced band, followed by two voiced frames; the two frames after it together (320 samples) have a
standard deviation of at least 300 in DVSI's decode; the first three and the last four frames of a
file are left out. DVSI's decode "follows" a build when their correlation over the 400 samples from
the onset frame's mid-point is above 0.8 and more than 0.1 above the other build's. Over the 1,940
P25 onsets of the 25 development and validation vectors:

| P25, 1,940 onsets | Count | The library at `cc272ef` |
|---|---|---|
| DVSI's decode follows the early start | 1,771 | early in 1,771 |
| DVSI's decode follows the late start | 72 | late in 43 (voiced only below 505 Hz), early in 29 |
| Both fit (within 0.1 of each other, both above 0.8) | 83 | |
| Neither fits | 14 | |

1,898 onsets were in line (above 0.8) at `cc272ef`, and 1,926 would be with the better start for
each. The onsets are selected by the level of the two frames after the onset frame, not by the
onset frame's own (622 of the 1,940 onset frames have a standard deviation under 300 in DVSI's
decode). Within this selection the start goes with the voicing bits and the pitch, and no
dependence on the frame's level was seen:

- A frame voiced in its first band only (harmonics 1 to 3): of the 41 such onsets at b0 56 and
  above, all are labelled late; of the 108 at b0 55 and below, 99 are labelled early and 9 fit both
  starts; none of either group is labelled the other way. At the pitch the library plays ([P25
  pitch](#p25-pitch)) the third harmonic lies at 502.8 Hz for b0 56 and at 508.1 Hz for b0 55
  (502.6 and 507.9 Hz at eq 46's nominal pitch).
- One voicing word starts early or late by pitch: `1000001100` early at b0 86 (f0 127.5 Hz) and late
  at b0 92 (121.7 Hz); `1011` followed by unvoiced bands early at b0 79 and late at b0 92, 99 and
  102.

The other codecs have no such onsets left. Of rate 33's 1,847 onsets labelled so, DVSI's decode
follows the late start in 1,461, all of them V/UV row 15, and the early one in 293 (rows 0 to 14);
both fit 83 and neither 10, and the library already started every labelled one as DVSI's decoder
does: in rate 33 the encoder signals the late start by the row. Of D-STAR's 1,945, DVSI's decode
follows the early start in 1,845 and the late one in none (both 96, neither 4). The rule below is
P25's own: its bit stream cannot signal the late start, and DVSI's decoder infers it.

**The rule.** The columns are fixed: [0, 505), [505, 1005), [1005, 1505), [1505, 2005), [2005,
2505), [2505, 3005) and [3005, 3505) Hz, and nothing above 3505 Hz. A column is covered when the
frame's voicing is voiced at both its ends, the voicing read at harmonic 1 for the end at 0 Hz and
at the harmonic at or next above the edge for the others (unvoiced when that harmonic is past L). At
a voiced frame after a frame with no voiced band, a voiced harmonic at frequency f sounds in the
mid-frame model only if the column that holds f, or the column that holds f + 250 Hz, is covered.
Every other voiced harmonic of the frame is unvoiced in the mid-frame model: it starts in the
frame's second half. If no harmonic sounds at mid-frame, the pulse train starts at mid-frame, as a
late onset's does. The rule contains the late onset it replaces: a frame voiced only below 505 Hz
covers no column.

**The evidence by onset.** A whole-frame form of the rule (an onset late exactly when no column is
covered, every voiced band early otherwise) finds 70 of the 72 late onsets and calls 3 of the 1,771
early ones late, where the library at `cc272ef` found 43 of the 72 and called none of the 1,771
late. (The script that counted these took column 0 as covered when its upper end alone is voiced;
with harmonic 1 checked as well, the counts are the same, 70 and 3.) Five onsets of the development
and validation vectors do not fit it, and they are not
explained: early on DVSI's decoder with no column covered, `clean` frame 477 (word `110000001100`,
b0 181), `t04` frame 208 (`100000001100`, b0 171) and `ucarm15a` frame 769 (`100011000110`, b0 164),
all at f0 72 to 79 Hz with two-band runs above 1.8 kHz; late with a column covered, `p01mirs` frame
2101 (`11110000000`, b0 106) and `ucarf15a` frame 2022 (`100000000001`, b0 153, whose twelfth band
has 11 harmonics). The constants of the whole-frame form, varied one at a time (late onsets found of
72 / early onsets called late of 1,771):

| First edge (Hz) / spacing (Hz) | 505 / 500 | 500 / 500 | 510 / 500 | 505 / 480 | 505 / 520 | 505 / 450 | 505 / 550 |
|---|---|---|---|---|---|---|---|
| Found / called late | 70 / 3 | 68 / 3 | 70 / 4 | 67 / 3 | 70 / 4 | 60 / 2 | 68 / 6 |

The first edge is also bracketed by the frames voiced in their first band only (above), at 502.8
to 508.1 Hz. The library's 505 Hz, 500 Hz and 250 Hz are its choices within what these
measurements allow: the onset labels bracket the first edge, the per-harmonic agreement below is a
little higher with a first edge of 500 Hz (0.9578) than of 505 Hz (0.9565), and reaches from 150 to
500 Hz score alike.

**The evidence by harmonic.** In the 1,771 onsets DVSI's decoder starts early, each voiced
harmonic's band of the first half (a Hann window over samples -40 to 120 of the frame, the band half
a harmonic spacing either side of the harmonic) is compared with the decodes of the two test builds;
the library's unvoiced part is DVSI's generator ([Unvoiced synthesis](#unvoiced-synthesis)), so both
are deterministic and differ only where the voicing does. A harmonic is kept when the residual
energies of the two builds' decodes against DVSI's decode, in its band, differ by more than 20% of
the band's energy. Over the 25,660 voiced harmonics so kept, the share
on which a reading agrees with the build that fits better (a harmonic that sounds with the one that
starts none late, a silent one with the one that starts all late) is 0.9565 with the rule, 0.9494
without the 250 Hz reach, and 0.8441 with every voiced harmonic sounding (the library at `cc272ef`).
By this measure the reach of 250 Hz scores alike anywhere from 150 to 500 Hz (0.9547 to 0.9565,
250 Hz best); a first edge of 500 Hz scores 0.9578, 505 Hz 0.9565, 495 Hz 0.9559 and 510 Hz 0.9538,
and a spacing of 495 Hz 0.9558, 500 Hz 0.9565 and 505 Hz 0.9482. An edge of 500 Hz would make the
frames voiced in the first band only at b0 56, whose third harmonic lies at 502.8 Hz, cover column
0 and start early, where DVSI's decoder starts them late; 505 Hz is kept (above). Above 3 kHz (2,641
of
these harmonics) the agreement is 0.878 with nothing covered above 3505 Hz, as the rule has it; it
would be 0.646 if a column above 3505 Hz counted when its lower end is voiced, and 0.70 if the
voicing were taken as voiced at every edge past L.

The profile below the top of a voiced run is measured on another selection: voiced runs of six
harmonics or more that end below L, the harmonics where the two builds' residual energies differ by
more than 10% of the band's energy, in 50 Hz bins that hold at least 20 of them, reported up to
1,100 Hz below the run's top. There the share of harmonics for which the decode that sounds them
fits better rises with their distance below the run's top: 0.06 to 0.25 within 250 Hz, about 0.5 at
250 to 400 Hz, 0.9 from 450 Hz and 0.99 from 700 to 1,100 Hz. It shows which decode fits better,
not that the top column of every run is silent on DVSI's decoder.

**Where it does not apply.**

- Bands that turn voiced inside a voiced run: applying the rule to every harmonic that turns voiced,
  not only at onsets, is worse than applying it at onsets only (development residual 0.0312 of the
  signal against 0.0309; correlation with DVSI's decode lower on 17 of the 20 development vectors
  and higher on 3). The rule of [D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns)
  differs: it reads the column samples of both frames, and inside a run a band that turns voiced
  starts with the frame where `M` holds at its column or at the column 250 Hz higher.
- Tone frames: a frame played as a P25 tone model keeps its own model at mid-frame ([P25 tones sent
  as voice frames](#p25-tones-sent-as-voice-frames)). On the decodes of the tone vectors the five
  tone measures the gates check keep their figures, while `tone.alignment_corr` moves on nine P25
  vectors ([November 2026: mid-frame voicing](#november-2026-mid-frame-voicing)).
- D-STAR and rate 33: no onset of D-STAR follows the late start, and rate 33's follow the row
  (above). (Labelled per harmonic, D-STAR's onsets fit a reading in its own columns, the top column
  of each voiced region labelled late more often than not: [D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns).) Rate 33's rows give parts of their frames a mid-frame voicing of their own ([Rate 33:
  rows 13, 15 and 17](#rate-33-rows-13-15-and-17)).

**Which pitch the columns are read at.** With eq 46's nominal pitch and with the pitch the library
played when this rule was fitted (the earlier law, up to 0.07% off it) the harmonic at or next above
an edge differed at 8 of the 208 × 7 pairs of pitch code and column edge, and at 6 of the 208 × 7
pairs with the edges 250 Hz lower; the edges are known to about 5 Hz only, so the measurements
cannot tell the two pitches apart. The library reads the columns at the pitch it plays: since the
P25 decoder change set eq 46 on the 2^-19 grid ([P25 pitch](#p25-pitch)), at which the harmonic at
or next above an edge differs from the nominal pitch's at 2 of the 208 × 14 pairs. No pitch code puts
an edge on a harmonic: the nearest lie 5.97e-5 of a harmonic away at the pitch P25 plays (b0 152 at
2005 Hz; 2.86e-4 under the earlier law, b0 3 at 1505 Hz) and 1.55e-4 at ProVoice's nominal pitch (b0
36 at 2755 Hz), so no build's arithmetic (x87 excess precision, fast-math's reordering) moves a
harmonic across an edge (`test_p25_onset_edges`). Since the change set the onset rule's "voiced
harmonics" are those P25's own model plays, by the columns ([D-STAR and P25: mid-frame voicing by
columns](#d-star-and-p25-mid-frame-voicing-by-columns)): a harmonic the word leaves unvoiced in a
covered column, or within 250 Hz below one, sounds with the frame too (at b0 86 with the word
1000001100, harmonic 18 joins 19 to 23). The labels of this section were not re-measured on it.

**How the library carries it.** The decoder marks every P25 and ProVoice voice frame that has a
voiced band `MBE_FRAME_MARK_ONSET_COLUMNS` (5.0f in the otherwise unused `PHIl[0]`, set in
`imbe_decode_voicing()`), and P25 no longer sets `MBE_FRAME_MARK_LATE_ONSET`. Before this change a
voiced P25 frame carried no mark, except one voiced only below 505 Hz, which carried the late
onset; the IMBE encoder's `cur_mp`, which the public decoder fills, carries the decoder's mark too,
no encoder reads it, and the encoder's bits do not change.

P25 renders with a rendering of its own (`MBE_VOICED_SYNTHESIS_DVSI_IMBE`), ProVoice, since the
pulsed state, with `MBE_VOICED_SYNTHESIS_DVSI_PROVOICE`, which shares P25's column rule and its
marking and leaves out only the pulsed state; until
the mid-frame columns change only that rendering used the mark to choose bands; since then none
does, P25's reading every voiced model but a tone model in its columns whatever its mark (it still reads the mark to
clear it from a model with no voiced band, below). D-STAR's and rate 33's play it as no mark, and
P25's plays rate 33's 6 and 7 as none; other marks in a caller's state still act as marks
([D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns), "Marks in a caller's state"). That rendering
also marks the model it renders, after the adaptive smoothing of the synthesis core and before
rendering, in `mbe_synthesize_subframes()` and in `mbe_advance_subframes()` alike: smoothing, which
channel errors bring, can voice harmonics of a frame decoded with none (in `dam_e2_sd`, frame 300
decodes with no voiced band at b0 25 and is marked `MBE_FRAME_MARK_PITCH_ABSENT`, and smoothing
voices its harmonic 12 alone), and the rule was measured on the model as synthesized. So in P25's
rendering the mark of the model rendered is the path's own:

- A model rendered with a voiced band that is not a tone model carries
  `MBE_FRAME_MARK_ONSET_COLUMNS` whatever mark it was decoded or repeated with. A mark in a
  caller's state on such a model is overwritten and has no effect, `MBE_FRAME_MARK_LATE_ONSET` and
  `MBE_FRAME_MARK_PREVIOUS_HELD` included (at the branch point, `adecb29`, a late onset there
  started the voiced bands at mid-frame after a frame with no voiced band, as P25's decoder itself
  then wrote it for a frame voiced only below 505 Hz, and a held previous model stood at mid-frame;
  `test_p25_rendering_marks_the_model_it_renders`).
- A tone model keeps `MBE_FRAME_MARK_TONE`: its own model stands at mid-frame.
- A model with no voiced band loses `MBE_FRAME_MARK_ONSET_COLUMNS` and keeps any other mark. Of
  those, `MBE_FRAME_MARK_PITCH_ABSENT` (the voicing before lasting to mid-frame),
  `MBE_FRAME_MARK_PREVIOUS_HELD` (the model before standing at mid-frame) and
  `MBE_FRAME_MARK_PITCH_CARRIES` (its pitch entering the next voiced frame's mid-frame pitch) still
  act there; a late onset and rate 33's 6 and 7 change nothing, there being no voiced band to
  start late and no row to read.

The prediction history (`prev_mp`) keeps the decoded model and its decoded mark. A frame played as
a P25 tone model carries `MBE_FRAME_MARK_TONE` in the model rendered and
`MBE_FRAME_MARK_ONSET_COLUMNS` in its decoded history. D-STAR's rendering keeps every mark's meaning
as at the branch point.

The mid-frame model (`mbe_subframe_mid_model()` in `src/core/mbe_subframe.c`) applied the rule only
to a frame so marked after a frame with no voiced band, until the mid-frame columns change; after a
voiced frame the mark changed nothing. The edges' positions among the frame's harmonics, `edge * pi / (4000 * w0)`, are formed in
double from the stored (played) `w0`, and they are compared as bit patterns before `ceil()` makes
an index of them (a caller's state can hold any pitch). ProVoice, which shares the data path and has
no DVSI reference, takes the same rule at its nominal pitch. An onset frame that the error rate
mutes, which the synthesis core advances without rendering (`mbe_advance_subframes()`), leaves the
pulse train and the noise generator where the same frame heard leaves them. It leaves no unvoiced
overlap, where a heard frame with unvoiced bands leaves the second half of its own subframe, and
that loss can change the frame after, over its first 80 samples, where the overlap is added. On the
test's input (an onset with unvoiced bands, and the same error history after it in both streams)
the frame after differs in its first 80 samples only (`test_p25_muted_onset_advances`); an onset
with no unvoiced band has nothing to lose, and where the error history after it differs (the
muted stream's error rate staying over the limit) the frame after can differ in its second half
too. A `MBE_FRAME_MARK_TONE` on a frame that plays no tone
model (a caller's state; decoded streams do not reach this) gives way in the history to
`MBE_FRAME_MARK_ONSET_COLUMNS` when the model has a voiced band and to no mark otherwise, and the
model rendered is marked by its own voicing as above. A decoded mark that a caller's state
overwrote is not recovered: a frame with no voiced band decoded as `MBE_FRAME_MARK_PITCH_ABSENT`
and then marked as a tone carries no mark in its history afterwards, and in the model rendered no
mark when it is rendered with no voiced band, or `MBE_FRAME_MARK_ONSET_COLUMNS` when smoothing
voices a harmonic of it (the three states then read 5, 0 and 5).

**What it gives.** Measured with a local script on the scratch builds (the correlation of our decode
with DVSI's, sample for sample over the whole file, and the residual energy over the signal's, means
over the vectors; a vector counts as up or down when its correlation changes by more than 1e-6, here
and in [Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)), P25 decoding DVSI's bits:

| P25 | Correlation | Residual | Vectors up / down |
|---|---|---|---|
| development, `cc272ef` | 0.9687 | 0.0633 | |
| development, whole-frame form | 0.9848 | 0.0311 | 16 / 0 |
| development, the rule | 0.9848 | 0.0309 | 20 / 0 (largest rise `lvl` +0.0879) |
| validation, `cc272ef` | 0.9855 | 0.0289 | |
| validation, the rule | 0.9896 | 0.0208 | 5 / 0 (`t04` +0.0122) |

(The scoreboard gives 0.9849 for the rule on development, [November 2026: mid-frame
voicing](#november-2026-mid-frame-voicing): it correlates only the frames where DVSI's decode lies
within 40 dB of its loudest frame, where the script takes the whole file.) On the 1,863 onset frames
where the rule and
the whole-frame form differ, the rule is closer to DVSI's decode in 1,668, and the residual of those
frames falls from 0.0928 to 0.0733 of their energy (their first half from 0.2262 to 0.1874); the
250 Hz reach is closer in 143 of the 166 onset frames it changes. The library's decode equals the
scratch build's with the rule on all 25 development and validation vectors (16-bit samples), before
and after the review's fixes (the rendering of its own and the marking after smoothing, above),
which change none of these decodes nor any of DVSI's tone vectors. What is left of P25's onsets is
in
[Known differences](#known-differences).

### Rate 33: rows 13, 15 and 17

TIA-102.BABA-1 Annex B gives each V/UV row one voicing, which holds for the whole frame, and the
specification defines no mid-frame model; the library's synthesis in two subframes needs a voicing
at mid-frame, which the text does not define. The readings of rows 13, 15 and 17 below are DVSI's
decoder's, listed with the departures that follow DVSI ([Departures from the
specification](#departures-from-the-specification), "Rate-33 V/UV vectors"), and row 17's own
voicing, at the frame's end, is Annex B's again. Annex B lists three rows with the pattern 1110 0000
(10, 13 and 15, column 0 first), three
all voiced (0, 1 and 3) and rows 17 to 31 all unvoiced, like row 16. On DVSI's decoder rows 13, 15
and 17 give each half of the frame its own voicing: the row's own pattern holds at the frame's end,
and the first half voices row 13's columns 0 and 1, none for row 15 and column 0 for row 17. Before
this change the library decoded each row's own pattern as printed, with two exceptions it had
measured: row 15 marked `MBE_FRAME_MARK_LATE_ONSET`, which voiced no band at mid-frame after a frame
with no voiced band ([Voiced onsets and offsets](#voiced-onsets-and-offsets)), and row 17 with
column 0 voiced in the whole frame ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)).

**Where the alignment was lost.** Rate 33's onsets were right (1,838 of 1,847 in line; [P25: the
bands that start with the frame](#p25-the-bands-that-start-with-the-frame) for the labels). Among
its voiced runs under a correlation of 0.8 with DVSI's decode, the ones looked at, picked from the
runs with the largest residual (examples, not a census), started in line and lost it abruptly at
the same sequences of rows: row 14 then row 15 (`fambm22a` frames 2364-2365 and 2561-2562, `t02`
1095-1096, `noisy` 623-624), row 17 then row 13 (`t02` 1618-1619) or 15 (`damcar` 359-360), and
rows 17, 14 and 13 (`fambf22c` 1474-1476). Those after a row-17 frame that the readings below left
out of line are in line with DVSI's decode once the pulse train is carried through the row-17 frame
(**The pulse train through a row-17 frame**, below).

**Each row**, measured on the same vectors and builds of `cc272ef` with switches of the same kind (a
mid-frame column mask per V/UV row, and the row's own columns; the masks go through the decoder's
own rule, a harmonic voiced when its 500 Hz column or the one 250 Hz above it is voiced), on the
row's loud frames: the frames coded with the row whose DVSI decode has a standard deviation of at
least 100 over the frame and the frame after it together, a file's first two frames and its last
frame left out. The residual is the energy of the difference from DVSI's decode over the signal's,
in the frame's first half, over the frame, and over the frame after it. The rows were read one
after another, so each figure names its baseline:

- **Row 15**: no voiced band at mid-frame, whatever the frame before; its own voicing 1110 0000. On
  its 1,619 loud frames the residual (first half / frame / frame after) is 0.2317 / 0.1253 / 0.0354
  as the library decoded it at `cc272ef` and 0.2306 / 0.1242 / 0.0343 so; with column 0, or columns
  0 and 1, voiced at mid-frame instead, 0.37 / 1.19 / 2.03 or worse. Whole files: development 0.9737
  → 0.9782 (12 vectors up, 3 down, the largest fall 0.0005, the largest rise `fambm22a` +0.0401),
  validation 0.9788 → 0.9833 (5 up). After a voiced frame the pulse train is not restarted: it runs
  through the unvoiced mid-frame at the old pitch and on at the new one, as the library's train
  already does when only one end of a subframe has a voiced band ([Synthesis in
  subframes](#synthesis-in-subframes)). Holding the previous frame's whole model to mid-frame
  instead is worse (development 0.9600).
- **Row 13**: columns 0 and 1 voiced at mid-frame; its own voicing 1110 0000. On its 638 loud frames
  the residual (first half / frame) is 0.0438 / 0.0378 at `cc272ef`; with row 15 already read as
  above, 0.0423 / 0.0352 as the library decoded row 13 (columns 0 to 2 at mid-frame), 0.0221 /
  0.0195 with columns 0 and 1, 0.1427 / 0.1062 with column 0 only and 0.5925 / 0.4598 with none.
- **Row 17**: column 0 voiced at mid-frame and no voiced band at the frame's end as decoded; its own
  voicing is Annex B's. On its 115 loud frames the residual (first half / frame / frame after) is
  0.0572 / 0.0712 / 0.6921 at `cc272ef`; with rows 15 and 13 already read as above, 0.0564 / 0.0703
  / 0.6057 as the library decoded row 17 (column 0 voiced in both halves), 0.0564 / 0.0635 / 0.1656
  so, and 0.1055 / 0.1884 / 0.1656 with nothing voiced at all (Annex B's row played as row 16, which
  the library had measured worse than its own reading: [Rate 33: V/UV
  vectors](#rate-33-vuv-vectors)); columns 0-1, 0-2, 0-3 or all voiced at mid-frame are a little
  worse than column 0 (0.0593 to 0.0631 in the first half). Whether a row-17 frame carries its pitch
  into the next voiced frame's first half cannot be told: with the mid-frame pitch between a row-17
  frame and a voiced frame taken as the geometric mean of the two instead of the voiced frame's, the
  development mean is 0.98064 against 0.98062, and 18 development vectors rise and 2 fall against
  `cc272ef`, against 19 and 1. The library keeps the voiced frame's pitch. (Both were measured
  with the voiced frame after a row-17 frame starting a new pulse train, before the row-17 train
  change below.)
- **Rows 1 and 3**: nothing resolved. With rows 15 and 13 already read as above, row 1 (482 frames):
  first half 0.03298 with every column voiced at mid-frame and in the frame's own voicing, as the
  library decodes it (0.0341 at `cc272ef`), 0.03257 with the mid-frame top column cleared, 0.03294
  with the frame's own top column cleared and every column voiced at mid-frame, and 0.0349 to 0.65
  with fewer columns at mid-frame; row 3 (416 frames): 0.02274 with every column voiced in both
  (0.0232 at `cc272ef`), 0.02406 with the mid-frame top column cleared, 0.02276 with the frame's own
  top column cleared and every column voiced at mid-frame. The library keeps Annex B's rows.

The three readings together are in the speech table of [November 2026: mid-frame
voicing](#november-2026-mid-frame-voicing); with them the residual energy against DVSI's decode,
mean over the vectors, falls from 0.0527 to 0.0388 of the signal on development and from 0.0423 to
0.0331 on validation. D-STAR sends four voicing bits, one per 1 kHz band, and no row table, so it
has nothing of the kind.

**How the library carries it.** The rate-33 decoder (`ambe2450_decode_vuv()`) marks a voice frame of
row 17 `MBE_FRAME_MARK_MID_COLUMN_0` (6.0f) and one of row 13 `MBE_FRAME_MARK_MID_COLUMNS_0_1`
(7.0f); row 15 keeps `MBE_FRAME_MARK_LATE_ONSET`, which in the rate-33 rendering now means no band
voiced at mid-frame, whatever the frame before (`starts_late()` in `src/core/mbe_subframe.c`; P25 no
longer uses the mark). The mid-frame model (`mbe_subframe_mid_model()`) takes each band's voicing
from the row the mark names, through `mbe_ambe2450_harmonic_voiced()` at the mid-frame model's own
pitch, whatever the two frames' voicing; a band with no state keeps none. Only the rate-33 rendering
(`MBE_VOICED_SYNTHESIS_DVSI_AMBE2450`) reads these readings. In D-STAR's rendering marks 6 and 7
play
as no mark and the late onset keeps its earlier meaning, no band voiced at mid-frame after a frame
with no voiced band only, so a mark in a caller's state changes no D-STAR stream; in P25's, which
marks the model it renders, a caller's mark on a voiced model is replaced by the onset mark and
has no effect, and on a model with no voiced band 6, 7 and the late onset change nothing ([P25:
the bands that start with the frame](#p25-the-bands-that-start-with-the-frame)) (what was compared
is in [November 2026: mid-frame voicing](#november-2026-mid-frame-voicing)). The
V/UV table's row 17 (`AmbeVuv[17]` in `src/internal/ambe3600x2450_const.h`) is all unvoiced again,
as in Annex B, so the frame's `Vl` and the unvoiced amplitude scale follow it, and so does the
public parameter decoder: `mbe_decodeAmbe2450Parms()` returns a row-17 frame with no voiced band and
every amplitude at the unvoiced scale `0.2046 / sqrt(w0)`, where 2.6.7 to 2.6.12 (and this series up
to `adecb29`) returned column 0 voiced. The frame's mark carries its mid-frame voicing, and only the
rate-33 process functions read it: the public `mbe_synthesizeSpeechf()` has no mid-frame model and
plays such a model as decoded, all unvoiced, unless its adaptive smoothing, which the error counts
in
the state drive, voices some of its harmonics. "No voiced band at the frame's end" is the frame as
decoded: under channel errors the adaptive smoothing can voice harmonics of any frame afterwards, a
row-17 frame's too (in `dam_e10_sd`, frame 500, row 17, plays its harmonic 4 voiced at its end, and
the voiced frame after it continues that pulse train). Such a frame's end is voiced, and its
mid-frame pitch is the one between two voiced frames when the frame before is voiced (in the
rate-33 rendering that of [Rate 33: mid-frame pitch](#rate-33-mid-frame-pitch)).

A lost-frame repeat and the speech under a tone's first frame go through `ambe2450_load_repeat()`,
which drops the marks of rows 13, 15 and 17: the mid-frame voicing belongs to a decoded frame,
between the frame before it and its own model, and the repeated model plays its own voicing through.
The fade-out of a mute already cleared the mark. A decoded frame under a mute (unheard) and the
first frame from silence keep their row's mid-frame voicing. From silence (the first frame after a
mute or a tone), a row-17 frame voices column 0 at the mid-frame pitch, which, while the frame has
no voiced band after smoothing, is the kept model's when the kept model was voiced and the geometric
mean of the two pitches when it was not (a row-17 frame that smoothing voices takes the pitch
between two voiced frames, or its own after an unvoiced kept model), and the new pulse
train launches its first pulse at the frame's start for rows 13 and 17 and at sample 80 for row
15, whose first half voices no band (`test_rows_13_15_17_after_mute_or_tone`); since the row-17
train change a frame voiced in its first half after such a row-17 frame takes on the train the
row-17 frame started (`test_train_through_row17`). That is a design
choice, not examined against DVSI's decoder: `dam_e10_hd` has six audible row-17 frames right after
a mute of the library's, two of them with a voiced kept model (frames 639 and 1172), and what DVSI's
decoder plays there was not measured. TIA-102.BABA-1 says nothing on either, as it knows no
mid-frame model.

**The pulse train through a row-17 frame.** Since the mid-frame voicing change a row-17 frame's own
model as decoded has no voiced band and its mid-frame model voices column 0. Without channel errors
the pulse train ran through such a frame: over its first half between the frame before and the
voiced mid-frame model, over its second out of the mid-frame model at that model's pitch (a
subframe with one voiced end keeps the train; [Synthesis in subframes](#synthesis-in-subframes)),
and the frame stored the phase it reached and the pulses still sounding. The frame after it then
dropped that train. The library takes the train a stream carries into a frame only from a frame
whose rendered model has a voiced band, or that is marked `MBE_FRAME_MARK_PITCH_ABSENT` or
`MBE_FRAME_MARK_PREVIOUS_HELD`, through whose first half the voicing before may have lasted; a
row-17 frame rendered with no voiced band was neither, so a frame voiced at mid-frame after it
started a new train at its first sample. Before the mid-frame voicing change the library voiced
row 17's own column 0, and the train ran on. Since the row-17 train change (`8cbef45`) the rate-33
rendering also takes the train from a frame marked `MBE_FRAME_MARK_MID_COLUMN_0` whose rendered
model has no voiced band (`train_runs_on()` in `src/core/mbe_subframe.c`; one whose end the
adaptive smoothing of channel errors voices carried the train before as after): the phase and the
pulses carried are those the row-17 frame stored, rendering and `mbe_advance_subframes()` load them
alike, and whether they sound depends on the frame after: they do when its mid-frame model is voiced
(a frame with a voiced band of its own at mid-frame, or another row-17 frame, which voices column 0
there), and a frame whose first subframe has no voiced band at either end starts afresh as before.

The cases were searched for on the development and validation vectors, on the scratch builds of
the study (whose rate-33 decodes of the 25 vectors, without and with the change, are those of the
library's builds of `e446950` and of the change, file for file): a frame decoded as a voice frame
(no tone, erasure, repeat, mute or silence flag) of row 17 whose rendered model has no voiced band,
with a voiced band in the frame before and in the frame after as the baseline renders them, frames
2 to n - 5 of each file. There are 12, all on development vectors; validation holds none. For each,
the voiced run that starts after it (its first 12 frames at most) is correlated with DVSI's decode
under the baseline (`e446950`, whose rate-33 decodes are `c3c6d42`'s) and under the change. For
runs of three frames or more, DVSI's lag against the baseline's decode is measured on the run's
second frame: DVSI's frame against the baseline's same frame, each under a Hann window of 160
samples, 1,024-point transforms, the lag that maximises the cross-correlation limited to the
cross-spectrum between 0.05 and 2.2 rad/sample, searched over half a period either way in steps of
0.05 sample (0.0007 to 0.0016 cycle at these pitches), in cycles of the frame's pitch as the
baseline renders it; the "coherence" is that maximum over the sum of the cross-spectrum's
magnitudes. Beside it stands the library's own train phase at the row-17 frame's end, which is what
a train that runs on brings into the next frame (the baseline's fresh start brings 0):

| Vector, row-17 frame | Next row | Run (frames) | Correlation, baseline → change | DVSI against the baseline's fresh start (cycle) | A train running on brings (cycle) |
|---|---|---|---|---|---|
| `clean` 1416 | 13 | 24 | -0.011 → +0.982 | -0.202 | -0.194 |
| `fambf22c` 1212 | 13 | 13 | -0.293 → +0.980 | +0.377 | +0.383 |
| `t02` 1618 | 13 | 27 | -0.354 → +0.990 | +0.265 | +0.275 |
| `fambf22c` 1474 | 14 | 13 | +0.768 → +0.997 | -0.000 on the run's second frame (coherence 0.96); +0.043, +0.039, +0.040 and +0.040 on its third to sixth | +0.037 |
| `ucarm20a` 609 | 12 | 5 | +0.060 → +0.979 | -0.002 on the run's second frame (coherence 0.79); -0.102 (0.64), -0.132 (0.99) and -0.127 (0.99) on its third to fifth | -0.125 |
| `t02` 1417 | 14 | 1 | +0.941 → +0.997 | run too short | -0.217 |
| `ucarm15a` 609 | 14 | 1 | +0.935 → +0.943 | run too short | -0.058 |
| `dam40` 19 | 15 | 8 | +0.996 → +0.996 | -0.001 | (+0.032) |
| `damcar` 359 | 15 | 20 | +0.995 → +0.995 | -0.001 | (+0.399) |
| `t01` 1120 | 15 | 13 | +0.998 → +0.998 | -0.003 | (-0.081) |
| `tia11` 2130 | 15 | 1 | +0.705 → +0.705 | run too short | (+0.261) |
| `tambf22a` 998 | 15 | 3 | +0.965 → +0.965 | -0.000 | (-0.211) |

The run after the frame comes nearer DVSI's decode by more than 0.02 of correlation in 6 of the 12
cases and farther in none. On the three runs that start with a row-13 frame after the row-17 frame
(`clean` 1416, `fambf22c` 1212, `t02` 1618) DVSI's lag on the run's second frame lies within 0.011
cycle of the phase a train running on brings (observed minus predicted -0.007854, -0.006153 and
-0.010282 cycle), and within 0.021 cycle on the second to sixth frames. On `fambf22c` 1474 DVSI's
lag on the run's second frame is that of a fresh start (-0.000037) and on its third to sixth frames
within 0.007 cycle of the running train's phase (+0.006367, +0.002810, +0.003869, +0.002920); on
`ucarm20a` 609 it is that of a fresh start on the run's second frame (-0.001799, at a coherence of
0.785), and on its fourth and fifth frames, at a coherence of 0.99, within 0.008 cycle of the
running train's (-0.007825 and -0.002084). The five frames followed by row 15 do not change, and are
not expected to: a row-15 frame voices no band in its first half either, so the subframe from the
row-17 frame's end to its mid-frame has no voiced band at either end and the train starts afresh at
mid-frame in both decodes; DVSI's lag on the four runs long enough to measure lies within 0.003
cycle of a fresh start (the bracketed phases do not apply).

The same measure of DVSI's decode against the change's decode, on the second to sixth frames of the
five long runs that change (the second to fifth of `ucarm20a`'s): a mean of -0.009310 cycle
(`clean` 1416; -0.015478 to -0.005516), -0.007554 (`fambf22c` 1212; -0.011379 to -0.004260),
-0.006825 (`t02` 1618; -0.008068 to -0.004488), +0.001588 (`fambf22c` 1474; -0.000037 to
+0.002894) and -0.003759 (`ucarm20a` 609; -0.008761 to -0.000099). The measure's own spread, on
the third to seventh frames of every voiced run of 8 frames or more (as the baseline renders them)
of the 25 vectors, with DVSI's frame rms 200 or more and a coherence over 0.9 (7,179 frames): a
median of -0.000015 cycle, 10th and 90th percentiles -0.002633 and +0.002584, 1st and 99th -0.008487
and +0.007295. The three runs that start with row 13 so keep a small lag of one sign, -0.007 to
-0.009 cycle in the mean, near the 1st percentile of that spread; it is not explained.

The train's rate through the row-17 frame and the frame after it is the library's rule: the frame
before's pitch through the whole row-17 frame (its mid-frame pitch is the voiced frame's), and the
next frame's own pitch through that frame's first half. Other rates, each as the change it would
make to the phase at the next frame's end against that rule, in cycles, for `clean` 1416,
`fambf22c` 1212, `t02` 1618, `fambf22c` 1474 and `ucarm20a` 609 (periods of the frame before / the
row-17 frame as decoded / the frame after: 42.95 / 42.95 / 35.22, 46.36 / 40.41 / 36.31, 35.76 /
35.76 / 32.13, 82.86 / 114.20 / 77.95 and 41.66 / 44.29 / 44.97 samples):

| Rate | `clean` 1416 | `fambf22c` 1212 | `t02` 1618 | `fambf22c` 1474 | `ucarm20a` 609 |
|---|---|---|---|---|---|
| the row-17 frame's second half gliding from the pitch before to its own decoded pitch | 0 | +0.127177 | 0 | -0.132502 | -0.056920 |
| its second half at its own decoded pitch | 0 | +0.254353 | 0 | -0.265003 | -0.113840 |
| the next frame's first half gliding from the row-17 frame's pitch to its own | -0.204611 | -0.111741 | -0.126269 | -0.162926 | +0.013695 |
| both frames as between two voiced frames (the mid-frame pitches with the wave, linear glides; what 2.6.7 to 2.6.12 did) | -0.419711 | +0.021955 | -0.255490 | -0.618585 | -0.089516 |
| the two mid-frame pitches as between two voiced frames, the intervals otherwise as now | -0.322651 | +0.009779 | -0.193832 | -0.484737 | -0.069436 |
| what remains with the library's rule (DVSI against the change, above) | -0.009310 | -0.007554 | -0.006825 | +0.001588 | -0.003759 |

Each alternative lies 0.05 cycle or more from the observation in at least three of the five runs
(the first two cannot be told from the library's rule on `clean` 1416 and `t02` 1618, where the
row-17 frame carries the pitch code of the frame before). What this establishes, and no more: of the
two decodes, the one that carries the train is nearer DVSI's; DVSI's lag agrees with the library's
train through the row-17 frame within 0.016 cycle on the second to sixth frames of all five long
changed runs; a train whose rate followed the row-17 frame's own decoded pitch, or mid-frame pitches
as between two voiced frames, would not. What DVSI's decoder holds internally is not known. Two of
the twelve are the losses of alignment that nothing tried had explained until then (`clean`
1416-1417 and `t02` 1618-1619, each a row-17 frame followed by a row-13 frame with a pitch step):
with the train carried through the row-17 frame they are in line with DVSI's decode (the runs'
correlations above), at the pitch steps as they were.

Where else the change acts. In every rate-33 bit file of DVSI's vectors at the top of its mode's
directory (106 files, `irstia` and `mark` left out), the frames whose decode changes form 20 runs,
each starting right after a row-17 voice frame whose rendered model has no voiced band: 12 where
the frame after has a voiced band of its own (`clean` 1417, `clean_dtx` 1417, `dam_e10_sd` 116, 396
and 527, `fambf22c` 1213 and 1475, `t02` 1418 and 1619, `tambf22b` 1642, `ucarm15a` 610, `ucarm20a`
610) and 8 where it is another row-17 frame (`dam_e10_hd` 387, 413 and 861; `dam_e10_sd` 387, 413,
573 and 861; `t02` 2184). `t02` 2183-2184 is the one such pair without channel errors (not among
the twelve, whose selection asks for a voiced band in the frame after): frame 2184 is quiet (DVSI's
rms 22.2), and its correlation with DVSI's decode goes from 0.5714 to 0.7603 (first half 0.6274 to
0.6793, second half 0.5730 to 0.9952), its residual over the signal from 1.4740 to 0.9991; one frame,
recorded as measured. The 13 changed frames of the channel-error decodes (`dam_e10_hd` 387, 413 and
861; `dam_e10_sd` 116, 117, 387, 396, 397, 413, 527, 528, 573 and 861) against DVSI's decode of the
same bits: each frame's correlation rises on 6 and falls on 7, and ten of them lie under an rms of 4
in DVSI's decode; at 10% channel errors the two decodes are far apart sample for sample, and these
frames decide nothing. Where the adaptive smoothing of channel errors voices a row-17 frame's end,
the train carried before the change as after: `dam_e10_sd` frame 500 (row 17, decoded with no voiced
band, rendered with harmonic 4 voiced) changes no sample after it.

What else holds. A frame whose first subframe has no voiced band at either end after a row-17 frame
rendered with none (a row-15, unvoiced or silence frame, or a lost-frame repeat of the row-17 frame,
which drops the mark through `ambe2450_load_repeat()`) starts afresh as before, since the interval it
starts with has nothing to carry the train. A row-17 frame from silence, after a mute or a tone,
starts a new train at its first sample, and a frame voiced at mid-frame after it now takes that
train on. No DVSI vector shows a row-17 frame right after a mute or a tone followed by a frame with a
voiced band of its own; `dam_e10_hd` 859-861 shows one followed by another row-17 frame, voiced at
mid-frame (a muted frame, an audible row-17 frame from silence, another row-17 frame, ten of whose
samples change), under 10% channel errors, where it decides nothing; the library makes no statement
there beyond what the rule gives. Under a mute a row-17 frame is not heard, and the frame that ends
the mute starts from silence and drops the train the unheard frame stored, as before.

A mark in a caller's state. With `MBE_FRAME_MARK_MID_COLUMN_0` on the last frame played
(`prev_mp_enhanced`) the rate-33 rendering loads the train that frame stored, as the two marks
already honoured there do: the samples and the returned states are those of the same state under
`MBE_FRAME_MARK_PITCH_ABSENT`. In D-STAR's, P25's and ProVoice's paths, in any of the three slots,
and in rate 33's other slots, the samples are those of the same state with no mark in that slot, and
the returned states are alike in everything but the mark word; the mark word is state like any
other and travels where a path copies the state that holds it (6 in `prev_mp` before a lost frame:
D-STAR's, P25's and ProVoice's repeat loads that history, so all three returned states hold 6 where
the stream with no mark holds 0; rate 33's keeps it in its history). A mark put in place of a real
one also removes the real one's effect: 6 over the `MBE_FRAME_MARK_PITCH_ABSENT` of D-STAR's or P25's
last frame played plays as no mark does, not as that frame with its own mark
(`test_row17_mark_in_codec_states`, `test_marks_act_in_their_own_rendering`). The rendering
honours the mark from the caller's state as written, without a check of where it came from; what it
decides is bounded (whether the train that state stored is loaded, its count, times and phase
range-checked as `load_pulse_train()` reads them). Never dropping a carried train at all (a scratch
variant: after a subframe with no voiced band at either end the stored phase is 0 anyway) gives the
same rate-33 decodes of the 25 vectors as honouring the row-17 mark alone; the library honours the
mark alone, so that a caller's state with no voiced band and none of these marks carries no train,
as before. What the change does to the speech vectors is in [November 2026: the pulse train through
a row-17 frame](#november-2026-the-pulse-train-through-a-row-17-frame).

The library's rate-33 encoder never sends rows 1, 3, 13, 15 or 17 (it searches rows 0 to 16, and
ties keep the lower index), so its bits, and the decode of its own bit streams, do not change.

### Rate 33: mid-frame pitch

When both frames have a voiced band and their pitches differ, DVSI's rate-33 decoder plays the
mid-frame model ([Synthesis in subframes](#synthesis-in-subframes)) at a pitch near the geometric
mean of the two, off it by an amount that varies with the two frames' periods. The library plays it
at the geometric mean times `1 - 3.0e-4` (the factor, below) times `1 + (h(Pa) + h(Pb)) / 2`, a wave
in the two periods fitted to DVSI's decodes, with the agreement given below (**The wave**).
TIA-102.BABA-1 defines no mid-frame model, so neither is a departure from its text. The factor was
measured and adopted first; the paragraphs on it below, to the PESQ figures, are that measurement as
made then.

**The drift.** With [Phase regeneration](#phase-regeneration) in place, rate 33's aligned voiced
runs also grew with their length, less than P25's: 0.020, 0.030, 0.050, 0.065, 0.064 and 0.13 of
the signal over frames 0 to 3, 4 to 7 and so on. Per pitch code rate 33's fundamental is right
(0.0074% rms over 81 codes; 0.0001 samples a frame at an unchanged code, 1466 frame pairs). The
drift comes with pitch changes. Per frame pair, by the steps of the pitch code into the two frames
the measurement spans (development vectors, 16,033 frame pairs; the counts and medians here are over
all of them, and the fits over the 15,972 with a fitted shift under 0.5 samples), DVSI's pulse train
moves against ours by a median of: no step 0.000 samples (1466 pairs); one step of -1 or +1 with no
step before it -0.015 and -0.018 (665 and 974), of -2 or +2 -0.016 and -0.008, of -3 or +3 -0.031
and +0.005 (39 and 87); a step before and none now -0.011 (2286); a step before and a step of -1 or
+1 now -0.031 and -0.014 (1468 and 2615), of -2 or +2 -0.027 and -0.011, of -3 or +3 -0.056 and
-0.023 (398 and 1040). A least-squares fit of one shift for every mid-frame pitch taken between two
different codes gives -0.010 samples over all 15,972 pairs under the cut, of which 2,390 have a
step of more than three codes (larger steps were not tabulated), and -0.011 over the 13,582 whose
two steps are both of at most three codes; the six medians of the rows with one step and no step
before average -0.014 (-0.015 weighted by their pairs). The fit holds for periods from 20 to 95
samples alike: what a mid-frame pitch about 0.025% low gives, since it holds for a quarter of the
160 samples on either side. The same fit gives +0.002 samples for D-STAR (15,496 pairs, the fit over
the 15,423 under the cut; its medians by step run from -0.022, no step before and a step of -2 now,
173 pairs, to +0.015, with -0.037 on the 68 pairs stepping down by three) and 0.000 for P25 with
[P25 pitch](#p25-pitch) (14,833 pairs, the fit over the 14,774 under the cut; medians within ±0.005
for steps of one and two codes, -0.007 and -0.015 on the 61 and 113 pairs stepping by three with no
step before). The figures are the same over the pairs whose two steps are both of at most three
codes: +0.002 (12,172 pairs) and 0.000 (8,605). Neither shows an offset for every pair of unequal
pitches.

**The factor.** As a factor on the mid-frame pitch, the residual energy of the development
vectors against DVSI's decode is 0.0683 of the signal at 0, 0.0573 at -2e-4, 0.0564 at -3.3e-4,
0.0601 at -4.5e-4, 0.0706 at -6e-4 and 0.0941 at -8e-4; the parabola through the three best has its
minimum at -2.9e-4. The library takes `1 - 3.0e-4` (`MBE_AMBE2450_MID_PITCH_FACTOR`), between two
voiced frames of different pitch only, and applies it to every step, the larger steps that were
not tabulated included; every other case of the mid-frame pitch, and every other codec, is
unchanged. With the same factor D-STAR's correlation with DVSI's decode falls from 0.9844 to
0.9765 and P25 is best without one. Measured on the build with it (16,042 frame pairs), the same fit
gives +0.001 samples over all 15,990 pairs under the cut and +0.001 over the 13,594 whose two steps
are both of at most three codes, and the medians by step lie between -0.012 and +0.013 for steps of
one and two codes and between -0.033 and +0.010 for steps of three: no step 0.000; one step of -1 or
+1 with no step before -0.008 and -0.012, of -2 or +2 -0.009 and -0.004, of -3 or +3 -0.025 and
+0.010; a step before and none now +0.006; a step before and a step of -1 or +1 now -0.008 and
+0.010, of -2 or +2 -0.004 and +0.013, of -3 or +3 -0.033 and +0.000 (the pair counts as before,
within four). The factor removes the mean offset; it does not explain how the offset varies with the
step's direction and the step before, which the factor alone did not model. (Each of these rows, of
frame pairs that span two mid-frame pitches, sums the offsets of two frames; the wave below is
fitted to single frames, on the periods of the two frames each mid-frame pitch lies between.) Half
a step of a 2^-10 grid in log2 f0 is -3.4e-4, but the frame-pair data cannot tell a floor to such a
grid from a constant (nor do the grids below fit the offsets per code).

The subframe synthesis learns which codec renders through a third rendering,
`MBE_VOICED_SYNTHESIS_DVSI_AMBE2450`, which every rate-33 synthesis call passes (decoded frames,
frames from silence, unheard frames under a mute, repeats, the mute's fade-out and the speech under
a tone's first frame), so that a frame advanced without being rendered leaves the pulse train where
rendering it would. D-STAR, P25 and ProVoice kept `MBE_VOICED_SYNTHESIS_DVSI` then; since the
mid-frame voicing change P25 renders with `MBE_VOICED_SYNTHESIS_DVSI_IMBE` (and ProVoice did until
the pulsed state, since which it renders with `MBE_VOICED_SYNTHESIS_DVSI_PROVOICE`), and D-STAR
keeps `MBE_VOICED_SYNTHESIS_DVSI`. The public `mbe_synthesizeSpeechf()` is unchanged.

For the factor, on the scoreboard (joint lags, against the commit before), rate 33 decoding
DVSI's bits:

| Partition | Waveform corr. vs DVSI, sample for sample | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| development | 0.9659 → 0.9720 | 1.295 → 1.286 | 0.9919 → 0.9919 | +0.014 → +0.014 | 6.270 → 6.271 |
| validation | 0.9722 → 0.9788 | 1.115 → 1.106 | 0.9883 → 0.9883 | -0.041 → -0.040 | 4.833 → 4.834 |

The waveform correlation rises on 19 of the 20 development vectors and 4 of the 5 validation ones
(`tambf22a` falls from 0.9741 to 0.9712 and `t04` from 0.9869 to 0.9856; against main they rise,
from 0.9455 and 0.9709; the other 23 rise by 0.0010 or more), and the distance from DVSI's decode
falls by 0.009 dB on both partitions. The gate sets `r33-mid-pitch-development`, `-validation`, `-tones`
and `-channel` were declared before scoring. `-tones` and `-channel` pass. `-development` fails on one check, that the waveform correlation fall by no more than 0.002 on any
vector (`tambf22a`, -0.0029), and `-validation` on one, that it rise by at least 0.003 with its 95%
interval above 0: the mean rises by 0.0066, but over five vectors, one of which falls, the interval
reaches -0.0001. Their other checks pass.

- **Tones and channel errors.** The decodes of 5 rate-33 tone vectors change, and no tone metric
  moves on any vector. Every rate-33 channel-error decode changes; the mute decisions do not, and on the loud
  frames neither decoder mutes the distance from DVSI's decode falls by 0.007 dB (hard decisions;
  0.010 soft).
- **Other modes.** The D-STAR and P25 decodes are unchanged on every speech, tone and
  channel-error vector, and every encoder's bits are unchanged.

PESQ against the input falls by 0.0026 ± 0.0018 (9 of 25 vectors better), ESTOI moves by +0.0000 ±
0.0000 (14 of 25), and PESQ against DVSI's decode by -0.0002 ± 0.0002 (10 of 25).

**The drift, frame by frame** (measured later, on the mid-frame voicing change's head `dac6d14`, on
the development and validation vectors). At a frame boundary both decoders play the frame's own
pitch, so the lag between DVSI's decode and ours can be read there from the phases of the harmonics,
and the change of that lag over one frame reflects that frame's mid-frame pitch: for an ordinary
glide, a mid-frame pitch off by a relative d moves the pulse train by about 80 d samples over the
frame (a step of about an octave inside the frame changes that, through the octave step of the pulse
train; the data hold a step from code 39 to 119). The selection: a frame counts when its decoded
model has a voiced band and its V/UV row is below 16 and none of 13 and 15 (rows 0 to 12 and 14). A
boundary between two such frames is used when DVSI's decode has a standard deviation of at least 150
over the 160 samples centred on it; the lag there is the weighted least-squares slope of the
cross-spectrum phase over the voiced harmonics of the earlier frame's model up to 2.2 rad/sample (a
Hann window, each harmonic's band half a harmonic spacing either side, a harmonic kept when the
squared coherence is at least 0.9 and the phase lies within ±2 rad, at least four of them). A
frame's increment is the lag at its end less the lag at its start, where both exist and the frames
on both sides of both boundaries count; a file's first two frames and its last two are left out:
18,890 frames on the 25 vectors, 4,436 with the same pitch code as the frame before and 14,454 with
another. Over the frames whose two boundaries each kept nine harmonics or more, the median absolute
increment is 0.018 samples at an unchanged code (4,005 frames: the level of everything that is not
the mid-frame pitch, the estimator's noise and the per-code pitch differences among it), and at a
step of the code by -5, -4, -3, -2, -1, +1, +2, +3, +4 and +5 it is 0.064, 0.066, 0.059, 0.052,
0.052, 0.049, 0.048, 0.048, 0.047 and 0.040 samples (85, 186, 469, 1,100, 2,269, 3,702, 2,406,
1,135, 555 and 273 frames), with means between -0.011 and +0.018: at steps of one to five codes
DVSI's pulse train moves against ours by about 0.05 samples either way. Larger steps are fewer:
0.036 samples at steps of 6 to 9 codes (456 frames), 0.050 at 10 to 19 (61) and 0.155 at 20 or more
(21; the six steps of exactly 20 codes 0.167). The same measurement on P25 and D-STAR (19,598 and
20,394 frames with a lag at both boundaries, no row filter; the medians over the 16,468 and 17,429
of them with nine harmonics or more at both) gives 0.018 and 0.023 at an unchanged code, 0.021 and
0.022 at a step of one code, 0.026 and 0.028 at two or three, 0.033 and 0.033 at four to eight and
0.038 and 0.039 at nine or more: no jump at a step like rate 33's in these grouped medians, which do
not show that P25 and D-STAR have no effect of the pitch (the slow rise with the step's size was not
examined).

**One term per pitch code.** With the increment of a frame that steps from code a to code b taken as
40 (x[a] + x[b]) samples, one unknown x per code (the relative amount by which that code's pitch is
lower in the interpolation than the pitch it plays), least squares over the 12,683 stepping frames
with nine harmonics or more and an absolute increment under 0.4 samples, with unknowns for the 101
codes with 30 uses or more (41 of the frames have an end outside them, and stay in the fit), takes
the median absolute increment from 0.0496 to 0.0225 samples (its standard deviation from 0.0720 to
0.0548); fitted on every other row of those frames and applied to the rows between, the codes chosen
on all of them, from 0.0497 to 0.0230. The 101 terms have a standard deviation of 6.5e-4 and run
from -11.2e-4 to +14.3e-4. They are smooth along the codes (the differences between neighbouring
codes have a standard deviation of 3.1e-4, where unrelated values would give 9.2e-4) and periodic:
maxima near codes 24-26, 47-49, 71, 95 and 116, minima near 34-37, 58-59, 80-84 and 101-106, a
period of about 23 codes, which is half an octave of pitch (a code is 0.0220427 of an octave). They
are not the rounding of a code's period, frequency or log2 period to any grid tried: periods of 1 to
1/32 sample and frequencies of 2^-12 to 2^-17 cycle per sample, each with floor and nearest, and
log2 periods of 2^-5 to 2^-12 with floor, nearest and ceiling, give correlations of 0.2 or less (a
period grid of 1/64 sample, tried in a separate run with floor and nearest, 0.16); and they are only
partly related to Annex A's own table, whose f0 differs from the library's law by -2.2e-4 to
+24.1e-4 (correlation 0.41).

**The wave.** Fitted directly on the increments, as -40 (h(Pa) + h(Pb)) for a frame stepping between
codes of periods Pa and Pb, with h(P) = C0 + C1 cos t + S1 sin t and t = 4 pi log2(P), P a code's
period in samples (2 pi / w0, from the code's law), on the 10,045 stepping frames of the development
vectors under the same cut: C0 = +4.244e-5, C1 = +9.374e-5, S1 = +8.388e-4 (amplitude 8.44e-4, phase
83.6 degrees). This three-term fit, the one the library plays, takes the median absolute increment
from 0.0493 to 0.0226 samples on development; applied to the 2,638 such frames of the validation
vectors, from 0.0506 to 0.0266. With a second harmonic (C2 cos 2t + S2 sin 2t) the development fit
comes out at C2 = -0.99e-5 and S2 = +4.80e-5, with C0, C1 and S1 at +4.220e-5, +9.330e-5 and
+8.400e-4, and gives 0.0225 on development and 0.0263 on validation (five terms fitted on the
validation frames themselves: 0.0267); the library leaves the second harmonic out. With the period
of the wave left free (scanned from 0.30 to 1.20 octave in steps of 0.005) the fit chooses 0.500
octave. h is positive where the fit puts DVSI's mid-frame pitch above the library's before this
change; with the wave the library plays the geometric mean times `1 - 3.0e-4` times `1 + (h(Pa) +
h(Pb)) / 2`. The wave taken at the mid-frame pitch itself instead of at the two codes fits about as
well: a median of 0.0228 against 0.0226 samples over all the development steps, 0.0216 against
0.0217 on the 8,827 steps of one to three codes, and 0.0354 against 0.0338 on the 1,218 steps of
four codes or more.

**Beside the P25 pitch.** P25's decoder plays each code at eq 46 on a 2^-19 grid and, between two
voiced models, a floored parabola for the square root of the product of the two fundamentals
([P25 pitch](#p25-pitch)); its mid-frame deviation, `q(y) / sqrt(y) - 1` (y the product's
mantissa), is a wave in the same variable with the same half-octave period as rate 33's h, of about
the same size (+7.7e-4 to -1.2e-3; rate 33's two constant terms sum to -2.58e-4, with a single
harmonic of amplitude 8.44e-4). The law the library played for P25 before that change was a wave in
the same variable too (a constant of -1.98e-4, a first harmonic of amplitude 5.00e-4), fitted to the
per-code drift over runs at one code, which is the mean of a frame's own pitch and its mid-frame
pitch and so carried half the mid-frame deviation. Rate 33 shows no deviation at an unchanged code
(above), where P25 does; the parabola in place of rate 33's wave (unfloored, between two different
pitches only) measured no better on its present pitch law in the study's scratch build, and gains
there only once that law's two constants are refitted: a later change set's, not this one's. In
rate 33 the pitch a code plays agrees with the library's law to 0.0074% rms on the codes measured
([Rate 33: pitch](#rate-33-pitch)), which bounds a periodic part there and does not exclude one, and
the wave fitted here is in the mid-frame pitch. In this light the factor `1 - 3.0e-4` of the earlier
change stands for the offset's mean, and the variation from frame to frame that it did not follow is
what the wave fits. What in DVSI's decoder makes the wave is not known.

**How the library plays it.** Between two voiced frames of different pitch the rate-33 rendering
(`MBE_VOICED_SYNTHESIS_DVSI_AMBE2450`) plays the mid-frame pitch as the geometric mean of the two
frames' stored pitches times `MBE_AMBE2450_MID_PITCH_FACTOR` times `1 + (h(Pa) + h(Pb)) / 2`, with
C0, C1 and S1 as above (`MBE_AMBE2450_MID_PITCH_WAVE_C0`, `_C1` and `_S1` in
`src/internal/mbe_subframe.h`) and P each frame's period from its stored pitch, all formed in double
from the stored floats with the double-precision `log2`, `cos` and `sin` and rounded to float once
(`ambe2450_mid_pitch()` in `src/core/mbe_subframe.c`), as the P25 law was formed, so that the pitch
does not move with a libm's single-precision functions. h lies within C0 ± sqrt(C1^2 + S1^2), so the
two factors together lie between 0.998898 and 1.000587: for any two positive normal floats the
pitch is positive and finite in double. A product past `FLT_MAX` gives no pitch (the two largest
floats, in either order, reach about 1.0005 `FLT_MAX`), which the synthesis rejects as it rejects
every pitch that is not positive and normal, and the mid-frame model then has no voiced band. No
pair falls under `FLT_MIN`: both pitches would have to lie within about 0.23% of it, where the two
factors together are about 1.0005. Every other case of the mid-frame pitch (equal pitches, one frame without a voiced band,
the pitch-carrying mark, a held previous model) and every other rendering are unchanged. The
results are in [November 2026: rate 33's mid-frame pitch
wave](#november-2026-rate-33s-mid-frame-pitch-wave).

### D-STAR and P25: mid-frame voicing by columns

In D-STAR's and P25's renderings the library voices the mid-frame model ([Synthesis in
subframes](#synthesis-in-subframes)) by one rule, for every band of every frame but those it does
not cover (a tone model, which keeps its own model at mid-frame; a model held by
`MBE_FRAME_MARK_PREVIOUS_HELD`; a frame with no voiced band not marked as lacking a pitch, which
keeps its own voicing; invalid harmonic counts or pitches; below, **Where it does not apply**): a band in
the same state in both frames keeps it, and one whose state differs is voiced at mid-frame according
to fixed 500 Hz columns, read in the two frames together (**The rule**, below). The rule was fitted
to per-harmonic labels of DVSI's decodes of the development and validation vectors, each label
saying which of two test decodes is nearer DVSI's decode in a harmonic's band (below), and the
DVMEGA DVstick 30's single-band step curves of [#157](https://github.com/arancormonk/mbelib-neo/issues/157)
come closer to the unit's with it. The labels are not observations of DVSI's decoder's band
states, and what DVSI's decoders compute is not known. TIA-102.BABA and the D-STAR description hold
a frame's voicing for the whole frame and define no mid-frame model, so the rule is a departure from
no text ([Departures from the specification](#departures-from-the-specification), "Where voicing
starts and ends inside a frame").

**What the library did.** Before this change (at `c3c6d42`, whose decoder code is `d08ea39`'s) the
mid-frame model of a D-STAR or P25 frame took its voicing by three rules:

- between two voiced frames, this frame's voicing for every band;
- at a voiced frame after a frame with no voiced band, every voiced band of the frame in D-STAR, and
  in P25 the voiced harmonics in, or within 250 Hz below, a covered 500 Hz column ([P25: the bands
  that start with the frame](#p25-the-bands-that-start-with-the-frame));
- at a frame with no voiced band marked as lacking a pitch (`MBE_FRAME_MARK_PITCH_ABSENT`) after a
  voiced frame, all of the frame before's voicing when its harmonic at or next above 505 Hz is
  voiced, else none ([Voiced onsets and offsets](#voiced-onsets-and-offsets)); the pulse train ran
  on through the frame only in the first case.

Inside voiced runs the first of these and the rule give different states for 27,615 of the 57,411
labelled D-STAR bands (0.481) and 31,263 of the 67,645 labelled P25 bands (0.462): about half.

**How the bands were labelled.** The method of [P25: the bands that start with the
frame](#p25-the-bands-that-start-with-the-frame): two scratch decodes, outside the repository, that
differ in the mid-frame state of the bands in question; each runs through the whole stream, so a
difference in the pulse train's state from earlier frames can enter the window. DVSI's decode is
compared with each in the band of every harmonic concerned (a Hann window of 160 samples, a
1024-point spectrum, the band from (h - 0.5) f0 to (h + 0.5) f0, complex spectra). A harmonic is
labelled when the two decodes' squared spectral differences from DVSI's decode differ by at least 0.2
of DVSI's energy in the band (inclusive), and takes the state of the nearer decode. Frames whose DVSI
decode has an rms under 30 over the frame are left out, and so are a file's first two frames and its
last, and bands where DVSI's decode has no energy. Voicing is the model's as rendered.

| Situation | Frames | The two decodes | Window | Harmonics labelled |
|---|---|---|---|---|
| onset | voiced, the frame before with no voiced band | every onset from mid-frame / every voiced band from the frame's start | 160 samples from 40 before the frame | the frame's voiced harmonics, at its pitch |
| lacking a pitch | no voiced band, pitch code outside D-STAR b0 64..119 (and under 120) / P25 b0 128..183, the frame before voiced | never held / always held | the frame's 160 samples | the frame before's voiced harmonics, at its pitch |
| in-run | voiced, the frame before voiced | this frame's state (`c3c6d42`) / the frame before's state for every band | 160 samples from 40 before the frame | this frame's harmonics whose voicing differs from the frame before's at the nearest harmonic (`floor(h f0 / f0' + 0.5)`, kept within 1..L') |

Over the development and validation vectors together: D-STAR 34,555 harmonics labelled on 2,214
onset frames, 12,438 on 1,039 frames lacking a pitch and 57,411 on 8,963 in-run frames; P25 30,393 on
2,244, 11,905 on 1,230 and 67,645 on 12,196. Inside voiced runs, the labels of the bands that change:
of D-STAR's 30,308 labelled harmonics that turn unvoiced, 0.710 are labelled held through the first
10 ms, and of the 27,103 that turn voiced, 0.233 are labelled late; of P25's 34,548 that turn
unvoiced, 0.624 held, and of the 33,097 that turn voiced, 0.275 late.

**D-STAR by voicing word and column.** In-run labels by the pair of b1 words (bits for the 1 kHz bands
from 0 Hz up) and the 500 Hz column of the harmonic: "early / late" for bands turning voiced,
"held / let go" for bands turning unvoiced.

| b1 before → now | Column | Count |
|---|---|---|
| 1100 → 1110 | 4 | 1,163 early / 40 late |
| 1100 → 1110 | 5 (top of the run) | 34 / 1,168 |
| 1110 → 1111 | 6 | 2,315 / 119 |
| 1110 → 1111 | 7 (top) | 183 / 683 |
| 1100 → 1111 | 4, 5, 6 | 1,413 / 4, 1,484 / 7, 1,362 / 79 |
| 1100 → 1111 | 7 (top) | 115 / 364 |
| 1000 → 1100 | 2 / 3 (top) | 731 / 18, 31 / 725 |
| 1101 → 1111 | 4 / 5 (top, voiced above in both frames) | 709 / 37, 39 / 352 |
| 1100 → 1000 | 2 / 3 (top) | 1,693 held / 68 let go, 70 / 1,901 |
| 1110 → 1100 | 4 / 5 (top) | 1,119 / 71, 58 / 1,258 |
| 1111 → 1110 | 6 / 7 (top) | 1,912 / 170, 200 / 667 |
| 1111 → 1100 | 4, 5, 6 / 7 (top) | 754 / 16, 844 / 16, 788 / 50 / 80 / 263 |
| 1111 → 1101 | 4 / 5 (top, voiced above in both) | 785 / 37, 53 / 398 |

Where the column above changes the other way, the labels say voiced at mid-frame: 1000 → 0100
column 1, 159 voiced / 3 not; 1101 → 1110 column 5, 119 / 0; 1110 → 1101 column 5, 52 / 2; 1011 →
1101 column 3, 46 / 0.

**The rule.** For each band (harmonic l) of the mid-frame model, of pitch w0:

1. Its state in the frame before and in this frame, each read at the frame's harmonic nearest the
   band: the position l w0 / w0' among the frame's harmonics (w0' the frame's pitch), formed in
   double from the stored pitches, a half added and its whole part taken (a position exactly halfway
   between two harmonics takes the upper one), kept within 1..L'. A frame with no voiced band counts
   as unvoiced in every band.
2. The same state in both frames: that state.
3. Different: voiced at mid-frame if and only if `M(k)` or `M(k')`, `k` the 500 Hz column of the
   band's frequency and `k'` the column of the frequency 250 Hz higher, where

       M(k) = (Sp(k) and Sc(k)) or (Sp(k) != Sc(k) and Sp(k+1) != Sc(k+1))

   with `Sp` and `Sc` the two frames' column samples (below; the samples, not the bands' own states),
   the sample past the last column unvoiced in both. A changing band is so voiced where the sample of
   its column, or of the column 250 Hz higher, is voiced in both frames, or where that sample changes
   and the sample above it changes too. The last column (above 3505 Hz in P25, from 3500 Hz in
   D-STAR) has no sample above it that can change, but its own sample can be voiced in both frames:
   in P25's `ucarf15a`, frame 2385, harmonic 22 turns from unvoiced to voiced above 3505 Hz, both
   frames' column-7 samples are voiced, and `M(7)` voices it at mid-frame.

The column samples:

- **P25 and ProVoice**, as PR #198 read them: columns [0, 505), [505, 1005), ... [3005, 3505) Hz
  and one above 3505 Hz; `S(0)` is harmonic 1's voicing and `S(k)`, k = 1..7, the voicing at the
  harmonic at or next above 505 + 500 (k - 1) Hz, unvoiced past L. **Exception**: in a frame of
  more than 36 harmonics whose last band (harmonics 34 to L) is voiced while harmonic 33 is not,
  the last band reads as unvoiced in the samples (the frame's own model still voices it).
- **D-STAR**, its own columns, 500 k Hz for k = 0..7: both columns of 1 kHz band j (b1's bit for
  it) take the band's voicing, read from the model at its lowest harmonic (harmonic 1 for band 0),
  unvoiced past L; a band's column and the column 250 Hz higher are found as the decoder finds them
  to voice a harmonic (`mbe_ambe2400_harmonic_voiced()`: `l 16 f0`, rounded to float, its whole part
  and that of half a column more, at most 7).

What it reduces to:

- at a voiced frame after a frame with no voiced band, `M(k) = Sc(k) and Sc(k+1)`: P25's rule of
  PR #198, unchanged but for the exception, and in D-STAR a voiced band starts at mid-frame where
  neither its column nor the column 250 Hz higher has both its own and the next sample voiced;
- at a frame lacking a pitch after a voiced frame, `M(k) = Sp(k) and Sp(k+1)`: only the bands of the
  covered columns (and those within 250 Hz below one) hold through the first 10 ms; the pulse train
  runs on through the frame when one does (the mid-frame model then has a voiced band: no change to
  the train's code is needed), and starts afresh when none does;
- between two voiced frames, `M` as stated, at the column of the band's frequency and at the column
  250 Hz higher.

**Agreement with the labels**, the library's three rules at `c3c6d42` → the rule (the share of the
labelled harmonics whose state a reading gives; in parentheses the harmonics labelled). The readings
are scored at the labelled harmonics' frequencies (at onsets and inside runs the frame's own, at
frames lacking a pitch the frame before's), where the library evaluates the rule at the mid-frame
model's; the two differ across a pitch step and at an onset after a code that carries
its pitch (D-STAR: 577 labelled onsets, 10,207 labelled harmonics, after such a code with a change of
pitch).

| | Onsets | Lacking a pitch | In-run |
|---|---|---|---|
| D-STAR, development | 0.8386 → 0.9749 (26,334) | 0.7120 → 0.9914 (9,784) | 0.5155 → 0.9482 (43,167) |
| D-STAR, validation | 0.8466 → 0.9572 (8,221) | 0.7280 → 0.9864 (2,654) | 0.5143 → 0.9309 (14,244) |
| D-STAR, both | 0.8405 → 0.9707 | 0.7154 → 0.9904 | 0.5152 → 0.9439 |
| P25, development | 0.9511 → 0.9574 (23,161) | 0.7741 → 0.9600 (9,156) | 0.5461 → 0.8888 (51,424) |
| P25, validation | 0.9434 → 0.9491 (7,232) | 0.7745 → 0.9280 (2,749) | 0.5497 → 0.8779 (16,221) |
| P25, both | 0.9493 → 0.9555 | 0.7742 → 0.9526 | 0.5470 → 0.8862 |

The rule was found on these labels, the development and validation vectors together: neither
partition is out of sample for it.

**Alternatives scored** on the same labels (both partitions together, agreement as above), and why
each loses:

- A band voiced at mid-frame if it is covered in the frame that voices it (the onset rule applied
  both ways): in-run D-STAR 0.8964 (on its own columns), P25 0.7801 (0.7869 with the last-band
  exception at any L). It calls 1101 → 1111 column 5 early, where the labels say late.
- D-STAR read on P25's grid (505 + 500 k Hz, at the harmonic at or next above each edge): onsets
  0.9657, lacking a pitch 0.9846, in-run 0.9308, against 0.9707, 0.9904 and 0.9439 on its own
  columns; on the DVstick it takes the onset steps from 0.034 to 0.088 and the offset steps from
  0.061 to 0.138 ([Synthesis in subframes](#synthesis-in-subframes)).
- No grid (a band voiced at mid-frame if the two frames still differ X Hz above it): in-run D-STAR
  0.8565 (X 250 Hz) and 0.7925 (500 Hz), P25 0.7708 and 0.7037.
- P25's first edge moved (onsets / lacking a pitch, PR #198's reading at both ends of a column, no
  exception): 480 Hz 0.9439 / 0.9305, 500 Hz 0.9504 / 0.9394, 505 Hz 0.9493 / 0.9362, 510 Hz
  0.9466 / 0.9339, 530 Hz 0.9388 / 0.9213. 505 Hz is kept: PR #198's bracket of 502.8 to 508.1 Hz
  came from whole-run labels, and the per-harmonic difference between 500 and 505 Hz is 0.001 to
  0.003.
- P25 read at the top end of a column only, instead of both ends: 0.9318 / 0.9011 against 0.9493 /
  0.9362.

**P25's exception.** Of the labelled harmonics of an isolated voiced last band in a frame of more
than 36 harmonics that lie in, or within 250 Hz below, a column that is covered when the last band
is read as voiced, the labels say voiced at mid-frame for 26 and not for 208 at onsets, and 5 and
198 at frames lacking a pitch (with the shipped exception none of those 437 lies so: that is the
exception); of an isolated
last band in a frame of 36 harmonics or fewer, 9 and 0 at onsets, 6 and 3 at frames lacking a pitch.
With the exception agreement rises at onsets from 0.9493 to 0.9555, at frames lacking a pitch from
0.9362 to 0.9526 and in-run from 0.8806 to 0.8862. Applied whatever L (the form first tried) it gives
0.8771 in-run, and P25's `alltone` falls from 0.8253 to 0.8202 (whole-file correlation with DVSI's
decode): its single tones at b0 2 and 3 (L 9, voicing word `001`) then start their pulse train at
mid-frame. Reading the samples above harmonic 36 from band 11 instead (onset and lacking-a-pitch
labels together): against no exception it turns 161 predictions voiced (labels voiced 48, unvoiced
113) and 433 unvoiced (31, 402); against the shipped exception it changes 174 (48, 126).

**Other frames with no voiced band.** After a voiced frame, a frame with no voiced band at a code
DVSI's encoders keep for such frames is not held by the library. Labels comparing the unswitched
scratch build with one that holds on any such frame (`offset8b.py`: only frames where the two
decodes differ) say held in 12 of
953 labelled harmonics at D-STAR b0 64..71, in 0.02 to 0.11 of them at 72..103 (by eights of the
code) and 0.33 to 0.39 at 104..119; P25 0.01 to 0.13 at b0 128..167, 0.33 to 0.38 at 168..183; rate
33's row 16, 97 of 8,487. The codes with a third or more labelled "held" are ones whose level DVSI's
decoder shapes ([Known differences](#known-differences)); the rule does not hold them.

**Which onsets showed it first.** Voiced runs of three frames or more after a gap of one frame, the
onset's index f with 3 <= f < n - 3 (n the file's frames) and a voiced frame before the gap at index
1 or more, selected by DVSI's rms (at least 100 on the last frame before the gap and 200 on the run's second
frame), the run before the gap in line with DVSI's on that one frame (lag within 0.6 sample,
coherence 0.8), and coherence of at least 0.8 on the run's second and third frames with their lags
within 0.08 cycle of each other: of D-STAR's 119 such onsets two were out of line before this
change, `t03` frame 895 and `ucarm15a` frame 315, both after a frame lacking a pitch (b0 34) that
follows a frame voiced only above the first band (b1 0100 at b0 123; b1 0111 at b0 115). DVSI's pulse
train stood 0.214 and 0.192 cycle off the library's fresh start there, and a train running on through
the gap frame at the last voiced pitch predicts 0.218 and 0.197; by the rule that frame's covered
column holds. Of the 55 D-STAR onsets after a gap frame at a code DVSI's encoder keeps for frames
with no voiced band, all 55 agree with a fresh start (within 0.08 cycle) and 7 agree with a train
running on as well: 48 tell the two apart. In P25, of the onsets after a frame lacking a pitch,
eight tell the two apart, all in favour of a fresh start (others agree with both, `dam` frame 787
among them): the frame before each of the eight is voiced in isolated bands only (words
`10000000000`, `100010000`, `100000000100`, `10000000`, `100000000000` three times, `10010001`), so no
column is covered.

**Where it does not apply.** The mid-frame model takes this frame's voicing, as before, for a frame
with no voiced band that is not marked as lacking a pitch. No mark of the frame before is read by
the rule: a frame after a P25 tone model takes the rule as after any model. The first form of this
change took this frame's voicing after a tone model; nothing measured supported that, and the
previous model's mark is caller state: one forged float there changed 159 of 160 samples of a public
P25 decode. Measured on the scratch build (the rule after a
tone model against this frame's voicing there): no P25 speech decode of the 25 development and
validation vectors depends on it; three P25 tone decodes do (`dtmf60` 19 frames, `dtmf180` 8,
`dtmfvbad` 12), and over those frames the residual against DVSI's decode, over the signal, is
0.0328 → 0.0321, 6.0324 → 4.3615 and 0.0384 → 0.0379 with the rule, the whole-file correlations equal
to four decimals. Everything `mbe_subframe_mid_model()` decides before it fills the bands stays as
it was: a tone model's own model at mid-frame (`MBE_FRAME_MARK_TONE` on this frame, which the data
path clears on a frame that plays no tone), the previous model held by a rate-33 silence frame
(`MBE_FRAME_MARK_PREVIOUS_HELD`), invalid harmonic counts or pitches. In D-STAR's rendering
`MBE_FRAME_MARK_LATE_ONSET` after a frame with no voiced band (no D-STAR decoder writes it) still
leaves no band voiced at mid-frame. Rate 33's rendering does not use the rule: its V/UV rows carry
the mid-frame voicing (in its in-run labels this frame's state at mid-frame holds in 0.958 of the
bands that turn voiced and 0.989 of those that turn unvoiced; the column rule would score 0.56
there), and it keeps the hold the three renderings shared before for a frame marked as lacking a
pitch, which no rate-33 decoder writes; the TIA-102.BABA rendering of `mbe_synthesizeSpeechf()` has
no mid-frame model.

**How the library carries it.** `mbe_subframe_mid_model()` (`src/core/mbe_subframe.c`) reads each
frame's column samples from its model as rendered, after adaptive smoothing, and needs no new state
and no new mark. Each frame's nearest harmonic is found in double (step 1 above): a float quotient
rounds a position just under a half to the half itself, and at one pair of P25 codes, b0 134 before
161 (word `000000000010`, then `111111111101`), harmonic 36 of the mid-frame model lies at 33.4999986
of the first frame's harmonics, which float takes for 34 (unvoiced) instead of 33 (voiced), so that
the first form of this change silenced a band voiced in both frames; at the pitches the P25
decoder plays, that is the only pair of different codes, in either order, at which the two readings
differ (no decode of the files compared changes with it: the 106 D-STAR and 103 P25 bit files at the
top of each mode's directory of DVSI's vectors, `irstia` and `mark` left out). The rate-33 path, and
the magnitude and phase interpolation, keep their arithmetic: the band's frequency and its position
among a frame's harmonics in float, the interpolation's weights, logs and means in double.
D-STAR's positions among the columns are formed with the decoder's arithmetic, from the stored pitch
over 2 pi in float: at every D-STAR pitch code they give the decoder's columns, each at least 6.58e-5
of a column from a whole number and 2.00e-5 from a half (1.81e-5 where the decoder's stored pitch
rounds an ulp otherwise, with fast-math and in 32-bit x87 builds; `test_dstar_column_edges`): these
margins are for the ordinary pitch codes 0 to 123, at which no build moves a harmonic to another
column, and the model's voicing and the columns' reading cannot round apart. Code 124, the silence
model (L 14, a pitch of 250 Hz), puts its harmonics exactly on column positions 0.5 and 1.0, a
margin of zero; it decodes with no voiced band, and a harmonic exactly on an edge lies in the upper
column (`test_dstar_exact_column_edge`). P25's edges are those of PR #198 (`test_p25_onset_edges`).
Where the mid-frame pitch is a frame's own (an onset after a code that does not carry its pitch, a
frame lacking a pitch, two frames of one pitch code) the mid-frame harmonics lie as far from the
edges; at an onset after a pitch-carrying code the mid-frame pitch is the two frames' geometric
mean, which there, as between any two frames of different pitch, can put a harmonic anywhere. A band that plays nothing in either frame (a caller's
voicing that is neither 0 nor 1) plays nothing at mid-frame where the rule applies, so it still plays
nothing wherever it did before. P25's rendering still marks every voiced model it renders but a tone
model `MBE_FRAME_MARK_ONSET_COLUMNS`, which replaces any other mark there, and clears that mark from a
model with no voiced band; what ended is the mark's part in choosing bands. ProVoice, which shares the
data path and has no DVSI reference, takes P25's rule.

**Marks in a caller's state.** As before this change set, `MBE_FRAME_MARK_PREVIOUS_HELD` and
`MBE_FRAME_MARK_TONE` on the current model give the mid-frame model before any rendering is told
apart. `MBE_FRAME_MARK_PITCH_ABSENT` on a current model with no voiced band invokes the columns' hold
in D-STAR's and P25's renderings and the kept 505 Hz hold in rate 33's. A voiced model P25's rendering
renders is marked again. The pulse train runs on through a frame lacking a pitch only when a band is
held (after P25 b0 97, word `10000000000`, which covers no column, it starts afresh). No mark of the
previous model chooses the columns' bands; its `PITCH_CARRIES`, `PITCH_ABSENT` and `PREVIOUS_HELD`
still act on the mid-frame pitch and the pulse train, as before.

**P25's own models by columns** (the P25 decoder change set). DVSI's P25 decoder voices a frame's own
model, at the frame's end (and as the next frame's first model), by the same 500 Hz columns as the
mid-frame model: column 0 from harmonic 1, column k from the harmonic at or next above 505 + 500 (k
- 1) Hz (unvoiced past L; the isolated last band read as unvoiced: P25's samples below), and a
harmonic voiced if the column of its frequency, or of its frequency plus 250 Hz, is. Where that
differs from the bands of three (the first one or two harmonics above the voiced region, mostly) the
library played noise where DVSI's decoder plays a harmonic. Measured on the library before the
change (with the pitch rules and the phase kernel's reading in), on unvoiced harmonics of steady
mixed voiced frames in line with DVSI's decode (from a run's third frame, the frame's voicing equal
to the frames' before and after), by the coherence of DVSI's decode with ours in the harmonic's band
(the two decoders run the same noise generator, so an unvoiced band of DVSI's is coherent with ours,
a voiced one is not): of the unvoiced harmonics directly above the voiced region 49% are under 0.5
(375 harmonics, median 0.516), one further up 27%, two further up 6%, four to six 1%. By what the
column reading of the same decoded bands gives: ours unvoiced and voiced by columns, 511 harmonics,
median 0.352, 66.9% under 0.5; ours unvoiced and unvoiced by columns, 3,674, median 0.996, 2.7%;
ours voiced and voiced by columns, 5,740, median 1.000, 0.3%; ours voiced and unvoiced by columns,
12. With the change (the library's own voicing, the column reading, as the classifier) the
unvoiced harmonics directly above the voiced region are 6% under 0.5 (300 harmonics, median 0.994),
one further up 10%, two further up 5%, four to six 2%; the 511 that turned voiced are coherent with
DVSI's decode at a median 0.995 (3.5% under 0.5); and over all 101,444 unvoiced harmonics of the
mixed frames (any neighbours) 3.9% are under 0.5 (from 8.1% of 117,309 by the bands). Labels (450
frames, 9,491 harmonics labelled DVSI-voiced or DVSI-unvoiced by a coherence under 0.5 or of 0.9 or
more, 6,082 voiced; the study's `uvband10c.py` on the library before the change): the standard's
bands agree on 9,026 (0.9510; 440 of ours unvoiced where DVSI's is voiced, 25 the reverse), the
column reading with the library's constants on 9,341 (0.9842; 103 and 47); over first edges from
450 to 1,000 Hz, reaches from 0 to 500 Hz and sampling by ceiling, rounding or floor nothing agrees
better than 9,345 (510 Hz, 250 Hz, ceiling), and edges of 495 to 510 Hz agree within 4 harmonics of
it. What is left (about 150) is spread over the whole width of a column. No constant is fitted: the
reading is the mid-frame rule's.

A band that turns takes its magnitude across the unvoiced factor `0.2046 / sqrt(w0)` (divided by it
to the voiced domain, multiplied by it to the unvoiced one; `magnitude_across()`, which the pulsed
state's bands use too). Not for a model with no voiced band (its bands and its pulsed state as
before), a tone model, or a band with no state. A model whose voiced bands all lie off the samples
(no DVSI encoding of the development and validation vectors has one; adaptive smoothing under
channel errors makes 76 among the eight channel-error vectors' 10,200 frames) plays with no voiced
band; kept with its decoded voicing instead, the channel-error vectors' log-spectral distance from
DVSI's decode rises (hard decisions 4.7357 → 4.7377 dB, soft 5.6455 → 5.6477, the scoreboard's
measure on the library with the first form of the mid-frame clause), and the speech vectors do not
change. The decoded model, its voicing, its magnitudes and the prediction history keep the decoded
values; the public decoder's `cur_mp` too.

**What reads which voicing.** The models' bands, and so the voiced train and its pulses, the unvoiced
synthesis's bands and the pulse train's carrying into the next frame (`train_runs_on()`: the last
frame's own model as rendered) follow the column reading. The decoded (smoothed) voicing of the
parameter sets still decides the mark of the model rendered (`mark_rendered_model()`), the mid-frame
pitch's cases (an onset, an offset, two voiced models), whether the mid-frame rule applies
(`column_rule()`), the reserved codes of the pulsed state (a model with no voiced band, which the
reading leaves alone), the phase kernel's input ([Phase regeneration](#phase-regeneration): the
step on the column-read voicing measured worse) and the low band's tracker: a tracker that moves on a
model the columns read as fully voiced (345 of the 32,370 rendered models with a voiced band of the
25 vectors and the channel-error vectors differ in that) measured no better (the scoreboard's
log-spectral distance from DVSI's decode 1.1019 → 1.1024 dB development, 0.7466 → 0.7473
validation, channel hard decisions 4.7357 → 4.7367; the waveform correlation equal to the fifth
decimal), so the low band's rule stands as it was.

**The mid-frame clause.** With the own models read by columns, a harmonic can be voiced at both ends
of a steady frame and unvoiced at its middle, if the clause ("a band in the same state in both frames
keeps it") compares the decoded states. Three forms, measured on the library with the four rules
(the scoreboard's own measures, each build aligned freely; the lags came out the same in every
build):

| Clause | Waveform corr. vs DVSI, development / validation | LSD vs DVSI (dB), dev / val | LSD on the channel vectors (dB), hard / soft |
|---|---|---|---|
| on the decoded states (as it was) | 0.994636 / 0.996855 | 1.1019 / 0.7466 | 4.7357 / 5.6455 |
| none (every band by the columns' M) | 0.994673 / 0.996786 | 1.0745 / 0.7309 | 4.7369 / 5.6454 |
| on the column-read states (the library) | 0.994727 / 0.996757 | 1.0608 / 0.7251 | 4.7309 / 5.6406 |

The waveform correlation splits (the column-read clause the best on development, the decoded one on
validation, by about 1e-4 either way); the log-spectral distance from DVSI's decode, the
scoreboard's second measure, is lowest with the clause on the column-read states on both
partitions and on the channel-error vectors, by 0.041 / 0.022 dB against the clause as it was. The
library takes it: a band in the same state in both frames' own models, as the columns read them,
keeps it, and a frame whose two own models read alike plays that reading at mid-frame. The three
forms differ, for example, between b0 20 and b0 41 with the word 11011(0) in both: the clause as it
was leaves harmonic 8 of the mid-frame model unvoiced, voiced in both own models; no clause leaves
harmonic 13 unvoiced (`test_p25_own_columns`). At an onset, where the frame before has no voiced
band, the reading makes a harmonic the word leaves unvoiced in a covered column, or within 250 Hz
below one, sound with the frame (at b0 86 with the word 1000001100, harmonic 18 joins 19 to 23;
`test_p25_onset_columns`); ProVoice keeps its decoded voicing and the clause as it was.

**What is not known.** What in DVSI's decoders the rule stands for: why, in most labelled bands,
the top column of a run of changing columns takes the unvoiced state even under a column voiced in
both frames (1101 → 1111 column 5: 39 labelled early, 352 late; 1111 → 1101: 53 held, 398 let
go). Why P25's in-run agreement (0.886) is lower than D-STAR's (0.944) was not examined. A D-STAR
voice frame after a tone, a mute or a resume takes the rule against the last model rendered, as it
took this frame's voicing before; DVSI's vectors hold no tone inside speech.

### D-STAR and P25: the pulsed state

DVSI's P25 and D-STAR encoders give a voice frame with no voiced band either one default pitch code
(P25 b0 25, D-STAR b0 34) or one of 56 others (P25 b0 128..183, D-STAR b0 64..119; [Voiced onsets and
offsets](#voiced-onsets-and-offsets)). The library reads those 56 as side information for the voicing
around them (the voicing before ends where the frame starts, and some codes carry their pitch into
a voiced frame after them), and played the frame itself as plain noise. DVSI's decoders play it
otherwise: each of the 56 codes names, for the frame's own subframe and for its mid-frame subframe,
which of eight 500 Hz channels are **pulsed**, a third voicing state in which a band plays no noise
and, instead, the response to a single pulse. TIA-102.BABA plays every unvoiced band as noise
whatever the pitch code ([Departures from the specification](#departures-from-the-specification)).
In P25 the library reads a code's pattern only while a validity counter, as US 7,970,606 describes
its decoder, allows it (clause 8); D-STAR reads it whenever the frame's model qualifies.

How often, on DVSI's bit streams of the 20 development and 5 validation vectors (the study's count,
frames with no voiced band as decoded): 8,087 of P25's 28,384 frames with no voiced band carry a
reserved code (3,333 whose pattern pulses the own subframe only, 3,914 the mid-frame subframe only,
840 both), and the other 20,297 all carry b0 25; 8,222 of D-STAR's 28,555 (3,326 / 4,055 / 841), the
other 20,333 all b0 34. The library's own decodes of the same bit streams give the same counts
(frames its public decoder decodes as voice with no voiced band, from `mbe_quality_eval --params`),
each code on 27 frames or more in P25 and 28 or more in D-STAR.

**The patterns.** Index b0 - 128 (P25) or b0 - 64 (D-STAR); eight characters for the frame's own
subframe, then eight for its mid-frame subframe, the lowest channel first, channel k from 500 k to
500 (k + 1) Hz (the last open above 3500 Hz), `P` pulsed:

```
 0 UUUUUUUP UUUUUUUU   14 UUUUUUUU UUPPPUUU   28 PPUUUUUU UUUUUUUU   42 UUUUUUUU PPPUUUUU
 1 UUUUUUPU UUUUUUUU   15 UUUUUUUU UUUPPPUU   29 PPPPUUUU UUUUUUUU   43 UUUUUUUU UPPPPPPP
 2 UUUUUUPP UUUUUUUU   16 UUUUUUUU UUUUUUUP   30 UUUPPPPP UUUUUUUU   44 UUUUUUUU PPUUUUUU
 3 UUUUUPUU UUUUUUUU   17 UUUUUUUU UUUUUPPP   31 UUPPUUUU UUUUUUUU   45 UUUUUUUU UUUPPPPP
 4 UUUUUPPP UUUUUUUU   18 UUUUUUUU UUUUUUPU   32 PPPPPPPP UUUUUUUU   46 UUUUUUUU UUPPPPPP
 5 UUUUPPPP UUUUUUUU   19 UUUUUUUU UUUUUPUU   33 UUUUPUUU UUUUUUUU   47 UUUUUUUU PPPPPPUU
 6 UUUPUUUU UUUUUUUU   20 UUUUUUUU UUUUUUPP   34 UUPPPPPP UUUUUUUU   48 PPPPPPPP PPPPPPPP
 7 UUUPPUUU UUUUUUUU   21 UUUUUUUU UUUUPPPP   35 PPPPPUUU UUUUUUUU   49 PPPPPPPP UUPPPPPP
 8 UUUUUUUU PUUUUUUU   22 UUUUUUUU UUUUPUUU   36 UUPUUUUU UUUUUUUU   50 PPPUUUUU PPPUUUUU
 9 UUUUUUUU UUUPPUUU   23 UUUUUUUU UUUUPPUU   37 PPPPPPUU UUUUUUUU   51 UUUPPPPP UUUPPPPP
10 UUUUUUUU UUUPUUUU   24 PUUUUUUU UUUUUUUU   38 UPPPPPPP UUUUUUUU   52 UUPPPPPP UUPPPPPP
11 UUUUUUUU UUPPPPPU   25 UUPPPUUU UUUUUUUU   39 PPPPPPPU UUUUUUUU   53 PPPPPPPP PPPUUUUU
12 UUUUUUUU UUPUUUUU   26 PPPUUUUU UUUUUUUU   40 UUUUUUUU PPPPPPPP   54 PPUUUUUU PPUUUUUU
13 UUUUUUUU UUPPUUUU   27 UUUPPPUU UUUUUUUU   41 UUUUUUUU PPPPPUUU   55 UUPUUUUU UUPUUUUU
```

Provenance: the 56 patterns are Table 2 of US 7,970,606 ("Interoperable vocoder"), which was read
([legal/patents.md](../legal/patents.md)); it lists them for the P25 pitch values 128 to 183. That
D-STAR's codes 64..119 carry the same 56 patterns at index b0 - 64, and that DVSI's decoders play
them so, are measurements (below). The codes whose own subframe has a pulsed channel (0..7, 24..39,
48..55) are exactly D-STAR's `MBE_FRAME_MARK_PITCH_CARRIES` codes (b0 64..71, 88..103, 112..119).
A separate threshold measurement of each code's channels on DVSI's decodes (the study's, outside
the repository) selects frames with no voiced band whose two frames before and one after have none
either (four consecutive rendered frames; frames 2 to n - 3 of a file of n frames, so that the
penultimate frame is left out, which, admitted, adds no frame), reconstructs the plain noise rendering of the 240
samples from the frame's start (a numpy port of the library's), leaves out frames where that
reconstruction's standard deviation is under 30, and takes, per code, subframe and 500 Hz channel,
the median over the code's frames of DVSI's departure from that rendering in the channel over the
rendering's energy there, a channel called pulsed at 0.5 or more, for codes with five such frames
or more: all 56 codes in both codecs, 2,340 P25 and 2,376 D-STAR frames, 7 or more per code. It
gives the table's state for 851 of the 896 (code, subframe, channel) triples in P25 and 860 of 896
in D-STAR. The differences do not all lie near the threshold. P25's 45: 12 are the channel just
below a pulsed one, measured 0.51 to 0.97 (the 250 Hz reach of clause 2 makes its upper half
pulsed); 28 are pulsed channels of the table measured under the threshold, 0.26 to 0.50, 26 of them
in the mid-frame subframe; 5 are unpulsed channels measured 0.51 to 0.71 in one subframe's window
of a frame whose other subframe is pulsed in that channel (b0 144's own channel 7, 0.71; b0 167's
mid-frame channels 0, 1 and 3; b0 177's mid-frame channel 0: the two 160-sample windows overlap by
half). D-STAR's 36: 13 of the first kind (0.56 to 0.96), 20 of the second (0.32 to 0.50, all in the
mid-frame subframe), and 3 others at 0.56 to 0.65: b0 80's own channel 7 and b0 96's mid-frame
channel 1, of the third kind, and b0 78's mid-frame channel 7 (0.58, 10 frames), which is not
explained. So the threshold measurement supports the table only coarsely; the result per code (with
the table every code's frames come nearer DVSI's decode, but P25 b0 131's: [November 2026: the
pulsed state](#november-2026-the-pulsed-state)) is the stronger check. (The "27 frames or more" and
"28 or more" above count every decoded frame with no voiced band, another selection.)

**The rule** (D-STAR's and P25's renderings; not ProVoice's, rate 33's or the TIA-102.BABA one),
each clause with what the study measured it against. The figures in this list are the study's, on
its scratch build (outside the repository: a worktree at `47d8cf3` with environment switches): its
variant `pzx3` (the rule, with clause 5 applied to every pair of frames with no voiced band) against
the same build with the head's switches (`k1`, whose decodes of the 25 vectors equal those of the
library's head `74126a3` file for file in all three modes), each alternative one switch away from
`pzx3`; `pzx3`'s decodes of the 25 vectors differ from the library's in frame 0 of every P25 file
and, in D-STAR, in six files by one unit in one to three samples (`pzm7`, with clause 5 as shipped,
differs from the library in P25's frame 0 only; [November 2026: the pulsed
state](#november-2026-the-pulsed-state)). They are on the 20 development and 5 validation vectors
(in sample: the rule was found on them): for each of the 22 vectors whose decode by DVSI has an rms
of 300 or more (`dam20`, `dam40` and `dam80` left out: at their level the 16-bit quantisation
dominates such shares) the residual energy against
DVSI's decode over that vector's whole signal energy, then the mean of those ratios over the 22
(each vector weighted alike), P25 / D-STAR; the adopted rule gives 0.02099 / 0.01475, the head
0.02302 / 0.01689 (with the energies pooled over the 22 instead, 0.01898 / 0.01357 against 0.02143 /
0.01579). The alternatives named in the clauses were measured so, on scratch library builds over all
frames of those vectors, except the pulse's position, the reach and the gain: when the gate sets
were declared those rested on the study's earlier numpy model of P25 frames deep inside runs with
no voiced band. The position and the reach were measured on scratch library builds after the
declaration (clauses 2 and 3); the gain was not measured again.

1. **Which frames.** A model the synthesis renders (after adaptive smoothing, like every rule of the
   mid-frame voicing) that has no voiced band and whose pitch is the one a reserved code decodes to
   has that code; so does the frame before's model, the last frame synthesized. Nothing but the
   model selects a pattern: the codes come from the models' voicing and pitch, and no frame mark is
   read. The pitch is compared by its bit pattern with the decoder's own pitch function for the
   code (P25 `imbe_code_pitch()`; D-STAR `ambe2400_code_pitch()`, one body, not inlined, for the
   decode and this test, so that the two cannot round apart in any build). A frame whose smoothing
   voices a band has no code; a caller's model at a reserved pitch with no voiced band has the
   code (it holds nothing that could be forged). In P25 a code is read only while the validity
   counter allows it (clause 8), and the frame before's code only where that frame read it; a code
   that is not read is no code for clauses 2 to 5. In every state a stream reaches, a frame whose
   code is not read and whose frame before's code was not read renders as at the branch point: a
   frame stores the read flag only for a code it read with a counter of 0 or more, and then stores a
   counter of 1 or more, so no stream stores the flag with a counter below 0. A caller's state that
   holds the two together (a counter of -1 with the flag set, say) still lets the frame before's
   code pulse its own model and choose clause 5's pitch: there the current frame's pitch stands at
   mid-frame where the branch point takes the geometric mean (`test_pulsed_read_flag`).
2. **Which bands.** In a model with a code, an unvoiced harmonic l, at f = l w0 4000 / pi Hz, is
   pulsed if channel min(7, floor(f / 500)) or channel min(7, floor((f + 250) / 500)) is pulsed in
   the pattern that applies to the model: the code's own-subframe pattern for the frame's own model
   (as this frame's end model and as the next frame's first), its mid-frame pattern for its
   mid-frame model. A pulsed band plays no noise in that model. The study's earlier numpy model of
   frames deep inside runs with no voiced band pointed to 250 Hz too. Measured since the
   declaration on the scratch build (its variants `pzr0` to `pzr500`, the adopted rule with the reach
   alone changed, against `pzm7`), the totals at a reach of 0 / 125 / 250 /
   375 / 500 Hz: P25 0.02117 / 0.02104 / 0.02099 / 0.02107 / 0.02114, D-STAR 0.01494 / 0.01487 /
   0.01475 / 0.01481 / 0.01489.
3. **The pulse.** The three models of a frame stand at its sample 0 (the last frame's own), 80 (the
   mid-frame one) and 160 (this frame's own). A model with a pulsed band emits one pulse 24 samples
   after its point (the last frame's own model at sample 24 of the frame, the mid-frame model at
   104) unless the model before it emitted one; voiced models and models without pulsed bands never
   count as having emitted. The pulse plays through the emitting model's pulsed bands under that
   model's falling crossfade and through the next model's pulsed bands under its rising one, over
   the 80 samples of that interval and no further. The response is the voiced pulse response
   ([Synthesis in subframes](#synthesis-in-subframes): the window rising over 16 samples, whole
   until one period of that model after the pulse, falling over 16; phases regenerated from
   `log2Ml`; the low-band gain and phase, since the low-band change DVSI's fixed filter and the cut,
   [DVSI's low band](#dvsis-low-band); the voiced render gain), each band's magnitude taken as a
   voiced one (the unvoiced factor `0.2046 / sqrt(w0)` of that model taken out). Measured against:
   a pulsed model before a voiced one emitting nothing, 0.02106 / 0.01484 (with the excitations kept
   apart, clause 4); a voiced model counting as having emitted, 0.02120 / 0.01496; the pulse ringing
   on past the next model's point, 0.02109 / 0.01491. From the numpy model: the pulse at 24 samples
   (correlation with DVSI's 0.93 there, about 0 at 23 or 25), at the voiced gain (a free fit gives
   0.99 of it; not measured again). Measured since the declaration on the scratch build (its
   variants `pzd22` to `pzd26`, the adopted rule with the pulse's position alone changed, against
   `pzm7`), the share of the frames with a pulsed
   model in reach at a pulse 22 / 23 / 24 / 25 / 26 samples after the point: P25 0.00108 / 0.00080 /
   0.00033 / 0.00082 / 0.00110, D-STAR 0.00116 / 0.00092 / 0.00033 / 0.00093 / 0.00119 (totals
   0.02184 / 0.02160 / 0.02099 / 0.02162 / 0.02185 and 0.01569 / 0.01547 / 0.01475 / 0.01547 /
   0.01571).
4. **The two excitations do not cross.** The voiced pulse train never plays through pulsed bands,
   nor a pulsed model's pulse through voiced bands. The train itself is as it was (its phase, its
   pulses, when it runs on or starts afresh, the octave step): a model with pulsed bands and no
   voiced band counts as one with no voiced band for every rule of the train and of the mid-frame
   voicing. What holds is that the voiced train is the branch point's for the same state. On the
   decoders' own streams that puts the train of a voiced frame after a frame whose own subframe is
   pulsed where it was before the change: from its first sample when its mid-frame model has a
   voiced band, from mid-frame when it has none (P25 b0 160 with no voiced band, then b0 97 with the
   voicing word `10000000000`, which covers no onset column: the train's pulses at samples 80 and
   148.227). With a caller's `MBE_FRAME_MARK_PITCH_ABSENT` on the frame with the code the train runs
   on as that mark makes it (the review: voiced b0 97, unvoiced b0 160 so marked, voiced b0 86, the
   first new pulse at sample 19.4444). The pulse of the frame
   before's own model, at sample 24, sounds in that model's pulsed bands only. Measured against: both crossing, 0.02151 / 0.01545; the train kept out of
   pulsed bands only, 0.02117 / 0.01507; the pulse kept out of voiced bands only, 0.02134 / 0.01514;
   with both crossing the class of voiced frames right after a frame whose own subframe is pulsed
   goes from a residual share of 0.00008 to 0.00053 (P25). The voiced train of an onset starting
   from the pulsed model's pulse instead of where it starts (measured with the excitations
   still crossing) takes the development correlation from 0.987 to 0.770 (P25) and 0.992 to 0.789
   (D-STAR): the runs after such onsets fall out of line. The frame before's voicing held to
   mid-frame on a frame with a code (as on a default-code frame lacking a pitch) takes the
   correlation from 0.9875 to 0.9756 (P25) and 0.9926 to 0.9794 (D-STAR): the voicing before ends where the frame
   starts.
5. **Mid-frame pitch.** Between two frames with no voiced band of which at least one has a code that
   is read (clause 1), in
   P25 this frame's pitch; in D-STAR the frame before's, unless this frame's own-subframe pattern has
   a pulsed channel, and then this frame's, or the geometric mean if the frame before's has one too.
   Between two frames with no voiced band and no code the geometric mean stays (DVSI's default-code
   frames all share one pitch, so nothing was measured there); every other case is as it was.
   Measured against: the geometric mean everywhere, 0.02118 / 0.01489; P25's rule in D-STAR,
   0.01495; D-STAR's rule in P25, 0.02129; after a voiced frame, the mid-frame model of a frame with
   a code at that frame's own pitch, 0.02099 / 0.01477 (no gain, not adopted). The study's tag
   `pzx3` applied clause 5 to every pair of frames with no voiced band; as written here it is
   `pzm7`, whose P25 decodes equal `pzx3`'s on all 25 files and whose D-STAR decodes equal them on
   19, the other six differing in one to three samples by one unit (the geometric mean of two equal
   pitches against the pitch itself).
6. **State.** Whether the mid-frame model emitted is needed by the next frame (its first model is
   this frame's own); it travels in a free element of the subframe state (`noiseOverlap[34]`,
   `MBE_PULSED_SLOT`). In P25's rendering two more follow it: whether the frame's code was read
   (`noiseOverlap[35]`, `MBE_PULSED_READ_SLOT`) and the validity counter it left
   (`noiseOverlap[36]`, `MBE_PULSED_COUNT_SLOT`, a whole number -10..3 as a float); D-STAR's
   rendering and the other subframe renderings leave 0 in both and read neither. All three are
   written with the pulse train by rendering and by `mbe_advance_subframes()` alike and read by
   their bit patterns, only in a state with the subframe rendering's mark: exactly `1.0f` means
   emitted, or read; anything else, or a state without that mark, not; a counter that is not one of
   the fourteen whole numbers, or a state without the mark, reads as 0. Any value a frame does not
   store reads as a fresh stream's (not emitted, not read, a counter of 0). The TIA-102.BABA
   rendering keeps its noise overlap in those elements; the mark says whose they are. A frame from
   silence plays nothing of the last frame's own model, so that model emits nothing; a rejected call
   (`frame_accepted()`) leaves the state untouched.
7. **Where it does not apply.** Tone models (which have voiced bands), the model a frame marked
   `MBE_FRAME_MARK_TONE` or `MBE_FRAME_MARK_PREVIOUS_HELD` puts at mid-frame (that mid-frame model
   takes no pattern), a mid-frame model with a voiced band (a caller's `MBE_FRAME_MARK_PITCH_ABSENT`
   on a frame with a code can hold voiced bands there), rate 33's rendering, the TIA rendering and
   any other value of the rendering passed to the subframe synthesizer (only D-STAR's and P25's have
   the codes): no pulsed band. ProVoice, which shares P25's data path and has no DVSI material, renders P25's
   way without the pulsed state (`MBE_VOICED_SYNTHESIS_DVSI_PROVOICE`), so its frames decode
   exactly as before.
8. **P25's validity counter.** US 7,970,606 describes its decoder (the description of its FIG. 8,
   steps 815 to 850): on a frame whose voicing quantizer value is zero, if the pitch value is one of
   Table 2's (b0 25, unvoiced in both subframes, and 128 to 183) the table is read only while a state
   variable, ValidCount, is 0 or more, and the variable is then raised by one, up to 3; at a pitch
   value not in the table it is lowered by one, down to -10; where a test fails the voicing is the
   received one, as in the standard decoder. The document gives no initial value. P25's rendering
   keeps such a counter in the stream's subframe state (clause 6), -10 to 3, as the frame before left
   it: a model with a voiced band (a tone model included) leaves it and has no code that is read;
   one with no voiced band at the pitch b0 25 decodes to or at one of the 56 codes' (the decoder's
   own pitch function, by bit pattern) raises it by one, to 3 at most, and a code of the 56 is read
   if and only if the counter the frame before left is 0 or more (b0 25 has nothing to read); one
   with no voiced band at the pitch of any other P25 code (b0 0..207, the same pitch function)
   lowers it by one, to -10 at least, and reads nothing; one at a pitch no code decodes to (the
   repeat limit's reset model and the initial state, at eq 46's nominal pitch of b0 134, or a model a
   caller built) leaves it and reads nothing. A repeat within the limit replays a decoded model and so
   counts by that model's code again. The
   counter is taken from the document, not measured: on DVSI's own encodings it never falls below 0
   (simulated on the rendered models of the library's decodes of all 103 P25 bit files, `irstia`
   and `mark` left out, the codes matched by the decoder's pitch bit for bit: 10,918 frames with one
   of the 56 codes, none arriving with the counter below 0, the lowest value reached 0; the
   channel-error files at 5% and 10% errors hold 45 to 59 frames with no voiced band outside the
   table each, all at codes' pitches, never enough in a row), so every decode of those files is the
   same with it as without. Its starting value (0, for a stream with no marked state), its reading
   on the rendered model (after adaptive smoothing; a muted or repeated frame counts by the model it
   renders or advances) and that a model at no code's pitch leaves it are this library's choices:
   the description counts received frames, and a lost frame decodes nothing; what DVSI's P25
   decoder does with lost frames was not measured. Through a run of lost frames (b0 220, no usable
   pitch code) the first four replay the last decoded model, the 4th muted and advanced, and so count
   by that model's code; the 5th installs the repeat limit's reset model, at no code's pitch, with a
   repeat count of 0 and renders it with the mute flag clear, and from there the run goes on so,
   muting at the 9th and the 14th, the reset model leaving the counter. So the first code frame
   after the run is read only if the counter stood at 0 or more before it: after three frames at b0
   176 (code 48) it stands at 3 through fifteen losses and the next b0 176 frame is read; after
   twelve frames at b0 127 it stands at -10 through them, and the next b0 176 frame is not read and
   raises it to -9. Frames the error rate mutes are decoded and advanced as rendered and count by the
   model advanced, after adaptive smoothing (which can voice a harmonic and so leave the counter):
   three muted frames at b0 127 with level amplitudes, advanced unvoiced, lower it from 1 to -2, and
   the next two b0 176 frames are not read (`test_p25_counter_through_a_mute`). No decode of DVSI's
   vectors reaches these cases. What an AMBE-3000 does with a stream from an
   encoder that does not use the codes was not measured. D-STAR keeps the rule with no counter, for
   two reasons: the document describes the P25 decoder only, and DVSI's D-STAR decodes contradict a
   counter of the plain form. In D-STAR's three DTX vectors a counter lowered by every frame with no
   voiced band at a pitch outside the table (b0 34 or 64..119), simulated the same way, stands below
   0 at 253 frames with one of the 56 codes (`clean_dtx` 111, `dam_dtx` 53, `noisy_dtx` 89; the
   frames that lower it are the repeats of those vectors' b0 127 frames, at b0 120's pitch, and
   their b0 120 frames), and on those 253 DVSI's decode is nearer the pulsed rendering than the noise
   rendering in 247, nearer the noise in 5 and level in 1 (residual energy of the frame against
   DVSI's decode at lag 0; the library's builds of the branch before the counter, `d3fb70f`, and of
   `74126a3`). With only the pitches of the ordinary codes 0..119 outside the table lowering it, 2
   frames of all 106 D-STAR files would go unread, one each in `dam_e10_hd` and `dam_e10_sd`;
   nothing supports that form, and it is not adopted. Whether the chip's D-STAR decoder has a counter
   of another form is not known.

**On hardware.** Issue #157's DVMEGA DVstick 30 package (a real AMBE-3000 decoding constructed
D-STAR streams; full-band rms in 16-sample blocks of the difference between a stepped stream and its
control, relative curves; mean absolute error against the unit's curve over the step frame's ten
blocks), the study's scratch build with the head's switches (`hwk1b`) against the same build with
`pzx3`'s (`hwpz3`), decoding the package's streams: the three steps from frames with
`b1` 0 at `b0` 100 (index 36: the own subframe pulsed in channel 2) to voiced frames 0.0604 → 0.0241,
0.0483 → 0.0152 and 0.0886 → 0.0204 (mean 0.06577 → 0.01992); in the second block (samples 16 to 31,
which hold the rule's pulse at sample 24) the unit's curve reads 0.73 at one position, the head's
0.52 and the rule's 0.76.
Of the three steps into such frames one comes nearer the unit's curve (0.0692 → 0.0512) and two go
further from it (0.1652 → 0.1891, 0.1599 → 0.1653; mean 0.13143 → 0.13520). Every other row of the package is unchanged to the digits printed (gain
steps 0.006, band-voicing steps 0.189, the pitch step 0.091, double steps 0.256, the steps at `b0`
20 and at the stream's own pitch). What the curves establish: they are full-band rms in 16-sample
blocks, for one D-STAR code (index 36), and they corroborate the fit at the steps out of such frames.
They cannot identify a channel, place the pulse at sample 24 within its block, or say anything of
P25: the table's mapping onto D-STAR's codes, the pulse's position, the alternation, the separation
of the two excitations, the reach and the mid-frame pitches rest on DVSI's decodes of its vectors
(above), not on the hardware.

**How the library carries it.** `mbe_subframe.c` holds the patterns and gives the three models of
a frame their pulsed bands (`BAND_PULSED`, beside voiced, unvoiced and none) after the models are
built (`pulsed_frame()`), asking the codec whose rendering it serves for a model's code through
`mbe_pulsed.h` (`mbe_imbe4400_pulsed_code()`, `mbe_ambe2400_pulsed_code()`, each with its decoder's
own pitch function); the voiced and the pulsed bands are rendered apart (`prepare_render()`), and
the unvoiced part leaves pulsed bands out. `mbe_subframe_mid_model()` carries clause 5's pitch and no
pulsed band; it consults the read state for that pitch (`pulsed_mid_pitch()`, through `read_code()`
and `last_read_code()`). `pulsed_frame()` reads the same codes for the three models' bands and
computes the counter's step for the next frame (`next_count()`; a pitch's P25 code, if any, from
`mbe_imbe4400_pitch_code()`). Rendering and advancing a frame store the emission, the read flag and
the counter alike.

**What is not explained.** P25's first frame with a code after a voiced frame stays further from
DVSI's decode than D-STAR's when its mid-frame subframe alone is pulsed (0.177 against 0.034, 606 and
738 frames; with the frames whose two subframes are pulsed added, pooled over both kinds, 0.168743
over 775 and 0.034568 over 912). The library's builds (`74126a3`
against the branch before the counter, `d3fb70f`, whose decodes the counter does not change), all
25 vectors, frames with a code and no voiced band as rendered, "first" those right after a frame
with a voiced band, the residual energy against DVSI's decode over those frames' own signal with
the energies pooled: P25 first own / mid-frame / both pulsed 0.080 → 0.046 (331 frames), 0.256 →
0.177 (606), 0.755 → 0.138 (169), later 0.202 → 0.018 (3,002), 0.176 → 0.021 (3,308), 1.328 → 0.053
(671); D-STAR first 0.093 → 0.037 (330), 0.123 → 0.034 (738), 0.797 → 0.035 (174), later 0.165 →
0.013 (2,996), 0.162 → 0.019 (3,317), 1.288 → 0.046 (667). Of the DVstick's three steps into frames
at `b0` 100 two go further from the unit's curves; rate 33's V/UV rows 18 to 31, rate 33's form of
the same state, are played since the rate-33 pulsed rows change by a rule of their own, measured on
DVSI's rate-33 decodes (each pulsed model plays its own pulse, the pitch code its time, no
alternation; [Rate 33: the pulsed rows](#rate-33-the-pulsed-rows)); ProVoice
has no DVSI material; the mid-frame pitch between two default-code frames of different pitch was
not measured and is left as it was. What the library does with another encoder's stream: in P25 a
frame with no voiced band at one of these codes plays its pattern only while the validity counter
allows (clause 8), in D-STAR whenever its rendered model qualifies. What was measured is DVSI's own
encodings, on which P25's counter never falls below 0, and one D-STAR code on a DVstick; what an
AMBE-3000 does with a stream from an encoder that does not use the codes, and whether its D-STAR
decoder has a counter, are not known (this library's encoders send none of the codes: [Encoder pitch
code of unvoiced frames](#encoder-pitch-code-of-unvoiced-frames)). What DVSI's decoders compute is
known only through their output.

### Rate 33: the pulsed rows

TIA-102.BABA-1's Annex B decodes V/UV rows 16 to 31 (`b1`) as "no band voiced". Before this change
the library played rows 18 to 31 on the band grid ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)),
every band as noise. On DVSI's rate-33 decodes of its own encodings these rows sound like rate
33's form of the pulsed state of [D-STAR and P25: the pulsed state](#d-star-and-p25-the-pulsed-state):
in the 500 Hz channels each row names, of the frame's own subframe and of its mid-frame subframe, a
pulse in place of noise, its time moving with the pitch code `b0` (0..119 on a voice frame; two
times, on a coarser grid, for the two rows that pulse both subframes). The library plays them by a
rule fitted to those decodes (below), in which every pulsed model plays its own pulse through its
own bands only, with no alternation between neighbouring models and no pulse shared between two
models, where the rule differs from P25's and D-STAR's. It is a fit, not DVSI's decoder: what that
decoder computes is not known, the fit leaves a residual (**What is not known or not matched**,
below), and the repeat rule rests on six frames. All of it was measured on DVSI's decodes of its own
encodings; no document was read for the pulse times, the two-position code or the repeat rule.

**The masks**, measured on DVSI's rate-33 decodes (channel k = 500 k to 500 (k + 1) Hz, the last
open above 3500 Hz, lowest first, `P` pulsed; the frame's own subframe, then its mid-frame
subframe):

```
row  own       mid-frame      row  own       mid-frame
18   PUUUUUUU  UUUUUUUU       25   UUPPPPPP  UUPPPPPP
19   PPUUUUUU  UUUUUUUU       26   UUUUUUUU  PUUUUUUU
20   PPPPPPPP  UUUUUUUU       27   UUUUUUUU  UUPUUUUU
21   UUPUUUUU  UUUUUUUU       28   UUUUUUUU  PPPPPPPP
22   UUUUUUUP  UUUUUUUU       29   UUUUUUUU  PPUUUUUU
23   UUUUUPPP  UUUUUUUU       30   UUUUUUUU  UUUUUUUP
24   PPPPPPPP  PPPPPPPP       31   UUUUUUUU  UUUUUUPP
```

Provenance: the fourteen pairs are measurements of DVSI's rate-33 decodes (the study of the rows'
columns in [Known differences](#known-differences) as it stood, re-measured for this change).
Checked afterwards, as an observation and not as their source: each pair equals exactly one entry
of the 56 patterns the pulsed state of D-STAR and P25 takes from Table 2 of US 7,970,606
([legal/patents.md](../legal/patents.md)), and the fourteen entries are those at indices 0, 4, 8, ...
52 (row 18 index 24, 19 28, 20 32, 21 36, 22 0, 23 4, 24 48, 25 52, 26 8, 27 12, 28 40, 29 44, 30
16, 31 20; no pair equals a second entry). US 8,359,197 (the half-rate vocoder), read for an earlier
study, lists indices 16 to 31 of its 5-bit voicing codebook as all unvoiced and says that duplicate
indices leave room for decoders that "may decode these indices into related but different voicing
states" such as "pulsed" ones (paragraph [122] as Google Patents numbers it); it does not say what
each index means.

**The rule** (the rate-33 rendering, `MBE_VOICED_SYNTHESIS_DVSI_AMBE2450`, alone). A frame's
samples are counted from its first (0); the last frame's own model stands at 0 (its point), the
mid-frame model at 80, this frame's own model at 160; a model's crossfade window runs from its
point - 80 to its point + 80, rising, then falling, linearly.

1. **Which frames.** A voice frame decoded with V/UV row 18..31 whose model, as the synthesis
   renders it (after adaptive smoothing), has no voiced band is a pulsed-row frame with its row
   and its pitch code. If smoothing voiced a harmonic it is not, and renders as before.
2. **Bands.** In a model with no voiced band an unvoiced harmonic at f Hz is pulsed when the
   channel of f or the channel of f + 250 Hz is `P` in the mask that applies: the pulsed state's
   rule (`min(7, floor(f / 500))`, f in double from the float pitch; the unvoiced factor
   `0.2046 / sqrt(w0)` of that model taken out of the magnitude; no noise in a pulsed band). The
   mid-frame model takes the row's mid-frame mask, the frame's own model the row's own mask, and
   the last frame's own model the own mask of the row the last frame stored (clause 6).
3. **Times.** A model with a pulsed band emits exactly one pulse: for rows 18..23 and 26..31 at
   `point + 1.016 b0 - 80.89` samples; for rows 24 and 25 the code is two positions on a grid of
   11, the mid-frame model's pulse at `80 - 75.6 + 11.875 (b0 mod 11)` and the own model's at
   `160 - 75.6 + 11.875 (b0 div 11)`. Each time is formed in double from the integers and the
   constants and rounded to float once (`MBE_AMBE2450_PULSE_SLOPE`, `_OFFSET`, `_GRID`,
   `_GRID_STEP`, `_GRID_ORIGIN` in `src/internal/mbe_subframe.h`): every exact time is n / 1000
   samples, which lies far from any rounding edge of a float, so every build rounds it alike.
4. **Playing.** A pulse is played by its own model only: the voiced pulse response (the response
   window with that model's period, the voiced render gain, the regenerated phases the model
   carries) through that model's pulsed bands, under that model's crossfade weight, in the two
   intervals where the model has weight: rising in the one that ends at its point, falling in the
   one that starts there. So a pulse sounds from its time (or from the start of its model's window,
   if it falls before it, its response already running) until its response ends (one period and 16
   samples after the pulse) or its model's window ends, whichever is first. A mid-frame pulse can
   fall just before the frame (b0 0: -0.89 samples). An own-model pulse lies at 79.11 to 200.01
   samples (one subframe pulsed) or 84.4 to 203.15 (rows 24, 25); the part of it from sample 160 on
   is played by the next frame (clause 6).
5. **The voiced train** is as before: a model with pulsed bands counts as one with no voiced band
   for every rule of the train and of the mid-frame voicing; the train never plays a pulsed band
   and a row's pulse never a voiced one (a model has pulsed bands only if it has no voiced band).
   The mid-frame pitch is unchanged (the geometric mean between two frames with no voiced band).
6. **Across the frame boundary.** The next frame needs the row and the code of a frame whose own
   model had a pulsed band: they give the last frame's own model its pulsed bands again (clause 2)
   and that model's pulse, at its time less 160. The frame stores them in the stream's subframe
   state, rendering and `mbe_advance_subframes()` alike; a frame whose own model has no pulsed band
   stores none. Starting from silence (`MBE_SYNTHESIS_FROM_SILENCE`) the last frame's own model plays
   nothing, so it has no pulsed band and no pulse.
7. **A lost frame's repeat** (and the repeat under a tone's first frame, which loads the same way)
   plays the state of the last subframe again: the repeated model's own mask in both halves
   (mid-frame model and own model), each of the two models emitting its pulse at the own pulse's
   position relative to its point (rows 18..23: `point + 1.016 b0 - 80.89` for both; rows 24 and
   25: `point - 75.6 + 11.875 (b0 div 11)` for both). A repeat of rows 26..31, whose own mask is
   empty, has no pulsed band. The repeat reads the repeated model's word only when that model is
   the band-grid model of the word's code (clause 1's test, but for the voicing, which the rendering
   checks); on any other model it leaves the word as the branch point did. This rests on six frames
   of two files (below), repeats of rows 19, 26 and 27: for the other rows of 18..25 it is the
   rule's extension, unobserved, and rows 24 and 25's held position (`b0 div 11` at both points)
   above all, since there the repeat moves the mid-frame pulse to a position the frame's own code
   gives its own model.
8. **Other frames.** A silence frame (`MBE_FRAME_MARK_PREVIOUS_HELD`) is no pulsed-row frame: its
   mid-frame model, the copy of the model before, takes no pulsed band, while the last frame's own
   model keeps its own (clause 6). Tone frames, mutes and their fades change only through clauses 6
   and 7: the mute's fade, a repeat toward a silent copy that carries no word, plays the last own
   model's pulse out (clause 6), and the repeat a tone's first frame fades out under its first 80
   samples is loaded as a lost frame's repeat, held word and all (clause 7).
9. **Nothing else changes**: D-STAR, P25, ProVoice, the public `mbe_synthesizeSpeechf()`, every
   encoder, and rate 33 on any stream that holds no frame of rows 18..31. For rate 33 that holds of
   the states the decoder writes, and of a caller's word on a model that fails the grid test (its
   pitch, by bit pattern, and its harmonic count against the band-grid model of the word's code): in
   `PHIl[0]`, or in the subframe state with a last model that fails it, the samples and the
   returned states are the branch point's. On the grid model of the word's code a word is read
   whatever the voicing: a model there with a voiced band plays no pulse, and its samples and slot
   37 are the branch point's, but a lost frame's repeat gives its word the held form (rows 18..25)
   or clears it (26..31), as it does for a decoded frame whose adaptive smoothing voiced a band, so
   the returned `PHIl[0]` is not the branch point's (`test_r33_repeat_on_a_voiced_grid_model`). A
   caller's state holding a word on the grid model of its code with no voiced band, or that word
   in the subframe state with such a last model, is to the library the state of a stream that had
   such a frame, and plays as one; two codes of one harmonic count share a grid model (codes 40 and
   41, say), so the test cannot tell them apart.

**The constants are fits.** In the study's alternatives (below) slopes 1.012 to 1.024 with the
pulse of b0 87 between 7.25 and 8.0 samples after its model's point (the rule puts it at 7.50; the
eleven variants' offsets 80.29 to 81.34) measure alike (shares of class (a), below, 0.00070 to
0.00074): neither the slope nor the offset is known more finely than those ranges. For the grid, origins of 75.1 and 76.1 measured a
little worse than 75.6; no finer step of the origin or of the step (11.875 a position) was tried.
The two-position code of rows 24 and 25 was found by fitting a pulse time per pitch code on DVSI's
decodes: b0 66, 67 and 69, for instance, give mid-frame pulses at 4.75, 16.6 and 40.0 samples with
the own pulse at 155.4 each (the rule: 4.4, 16.3 and 40.0, and 155.65); b0 110, 111, 112, 114 and
115 give 2.25, 17.0, 28.75, 52.25 and 64.0 with the own pulse at 203.0 (the rule: 4.4, 16.3, 28.2,
51.9 and 63.8, and 203.15).

**How it was measured.** All on DVSI's rate-33 material, the held-out partition (`irstia`,
`mark`) never read. The study's scratch build is a worktree at `47d8cf3` with environment switches
(outside the repository): its variant `k1` is the baseline (its rate-33 decodes of all 88 vectors
below equal those of the tools built from the branch below, `2b6b6c2`, file for file), its variant
`rs1` the rule. The library's figures below are its builds': the tools built from `2b6b6c2` (the
baseline, whose decodes equal `k1`'s) and from this change (`a1398df`; Release, SIMD, IEEE
arithmetic, as the scoreboard builds them), decoding DVSI's bits. Selections and measures:

- "22 vectors": the development and validation vectors whose DVSI decode has an rms of 300 or more
  (all but `dam20`, `dam40` and `dam80`).
- A frame "with no voiced band": every `Vl` of the model as synthesized (the `--params` dump's
  `synth` record of the baseline's decode) is 0. A "voice frame": `flags & 0x1F0 == 0` in that dump.
  A frame's row is its `b1`.
- "Residual share" of a class of frames: per vector, the energy of (our decode - DVSI's decode) at
  lag 0 over the class's frames divided by the energy of DVSI's decode over the whole vector; then
  the mean of the 22 ratios. The classes (a) to (d), over frames 1 to n - 1 of each vector: (a)
  voice frames with no voiced band of rows 18..31; (b) other frames with no voiced band right after
  a voice frame of rows 18..25 (a pulsed own subframe); (c) frames with a voiced band right after
  one; (d) the other frames with no voiced band; (e) the rest, frame 0 of each vector included.
- "Residual over the frames' own signal": the residual energy over a set of frames divided by
  DVSI's energy over the same frames, the energies pooled over the 22 vectors (or, for a named
  file, over that file's frames).

Results, baseline → this change, the library's builds (each figure equals the study's `k1` → `rs1`
to the digits printed, except where said):

- Whole-file correlation with DVSI's decode at lag 0, mean over the 20 development / 5 validation
  vectors: 0.992438 / 0.992119 → 0.993226 / 0.993591. All 25 vectors rise: `dam10` +0.000007,
  `t02` +0.000056, `damcar` +0.000058, `dam80` +0.000083, `t01` +0.000084, `t03` +0.000093, `dam`
  +0.000140, `dam_10ov` +0.000147, `dam20` +0.000158, `dam40` +0.000176, `tia11` +0.000408, `noisy`
  +0.000687, `lvl` +0.000802, `fambf22c` +0.000821, `tambf32e` +0.000867, `ucarf15a` +0.001163,
  `tambf22b` +0.001172, `tambf22a` +0.001366, `tambf32b` +0.001435, `fambm22a` +0.001615,
  `ucarm20a` +0.001726, `clean` +0.001726, `ucarm15a` +0.001875, `p01mirs` +0.002653, `t04`
  +0.003794.
- Residual shares: (a) 0.00249 → 0.00071 (5,667 frames); (b) 0.00013 → 0.00004; (c) 0.00031 →
  0.00012; (d) 0.00019 → 0.00019; (e) 0.01294 → 0.01294; total 0.01606 → 0.01399.
- Residual over the frames' own signal (class (a)), by row: 18 0.28 → 0.13; 19 0.40 → 0.09; 20
  0.70 → 0.19; 21 0.21 → 0.13; 22 0.16 → 0.12; 23 0.36 → 0.18; 24 1.24 → 0.29; 25 0.38 → 0.12; 26
  0.36 → 0.22; 27 0.14 → 0.08; 28 0.73 → 0.29; 29 0.74 → 0.41; 30 0.22 → 0.07; 31 0.29 → 0.12.
- The same by the frame's neighbours (class (a), frames 1 to n - 2): both neighbours voice frames of
  row 16 with no voiced band: 1,417 frames, 0.19 → 0.08; at least one neighbour that is not a voice
  frame with no voiced band (a voiced frame, or a silence, tone or lost frame): 1,560 frames, 0.55 →
  0.22; the others, selected as both neighbours voice frames with no voiced band, at least one of
  another row than 16: 2,687 frames, 0.59 → 0.24 (30 of them have neighbours of rows 16 and 17 only,
  no pulsed-row frame among them). Of the three selections the frame between two frames of row 16
  is matched best.
- PESQ-NB (pesq 0.0.4; the mean change ± its standard error with the vectors clustered by recording,
  [Method](#method); the 25 vectors come from 12 recordings, so a 95% interval takes t with 11
  degrees of freedom, 2.201 standard errors) on the 25 vectors: against DVSI's decode 4.244 → 4.331
  (+0.0867 ± 0.0075, interval +0.0702 to +0.1032, 25 of 25 up); against the input 2.959 → 2.950
  (-0.0094 ± 0.0071, interval -0.0250 to +0.0062, 9 of 25 up); DVSI's own decode against the input
  scores 2.980. The change moves the decode toward DVSI's and costs a little spectral distance from
  the input (+0.008 and +0.016 dB, the scoreboard's, [November 2026: rate 33's pulsed
  rows](#november-2026-rate-33s-pulsed-rows)); PESQ against the input establishes neither a fall
  nor its absence (the interval holds both 0 and a fall of 0.025), and no equivalence was tested.
- Other vectors (hard-decision decodes unless said; "changed frames" are the frames whose samples
  differ between the two builds; the residual over the frames' own signal on them). Tone vectors:
  `alert` 22 of 483 frames, 0.6014 → 0.5939; `alltone` 320 of 950, 0.0635 → 0.0576; `cp0` 2 of 50,
  0.0016 → 0.0015; `cp1` 2 of 50, 0.0012 → 0.0012; `cpvbad` 56 of 418, 1.0065 → 0.2865; `dtmf` 416
  of 5,293, 0.1159 → 0.1048; `dtmf60` 582 of 1,920, 0.0189 → 0.0121; `dtmf180` 583 of 5,760, 0.0019
  → 0.0012; `dtmfvbad` 471 of 2,226, 0.1791 → 0.0849; the other 43 tone vectors have no changed
  frame. Channel errors: `dam_e1_hd` 175 of 1,275, 0.4142 → 0.2132; `dam_e2_hd` 175, 0.3804 →
  0.1891; `dam_e5_hd` 174, 0.4070 → 0.2215; `dam_e10_hd` 122, 0.6054 → 0.4797; soft decisions
  `dam_e1_sd` 175, 0.4142 → 0.2132; `dam_e2_sd` 175, 0.3804 → 0.1891; `dam_e5_sd` 176, 0.4033 →
  0.2286; `dam_e10_sd` 213, 1.1661 → 1.1599. DTX: `clean_dtx` 374 of 3,700, 0.6022 → 0.2260;
  `dam_dtx` 111 of 1,275, 0.4137 → 0.2101; `noisy_dtx` 154 of 3,588, 0.4124 → 0.1556. Of these the
  study's scratch build gave other figures for five tone vectors, `alltone` (325 frames, 0.0917 →
  0.0861), `cp0` and `cp1` (3 frames each), `dtmf` (461, 0.3419 → 0.3357) and `dtmfvbad` (475,
  0.2055 → 0.1142), and for `dam_e10_hd` (252, 1.7500 → 1.7340): it changed, besides, the voiced
  frames that start from silence after a tone or a mute, where the library changes nothing, and the
  frames those reach, and its `cp0` and `cp1` frame 0 (below).

Alternatives, each one switch away from the rule on the study's scratch build (correlation
development / validation; class (a)'s residual share; the rule itself 0.993226 / 0.993591, 0.00071):

| Alternative | Correlation | Share (a) |
|---|---|---|
| a pulsed model emits only if the model before it did not (P25's and D-STAR's alternation) | 0.993172 / 0.993528 | 0.00079 |
| a pulse played through the pulsed bands of both end models of each interval it reaches | 0.993196 / 0.993529 | 0.00079 |
| ... through its own model and the model after it | 0.993203 / 0.993565 | 0.00076 |
| ... through its own model and the model before it | 0.993218 / 0.993557 | 0.00073 |
| the response window's whole part twice the model's period | 0.993169 / 0.993487 | 0.00084 |
| ... half the period | 0.993197 / 0.993597 | 0.00074 |
| ... 40 samples | 0.993225 / 0.993617 | 0.00070 |
| ... 80 samples | 0.993164 / 0.993485 | 0.00084 |
| the pulse's amplitude 0.9 of the voiced response's | 0.993214 / 0.993591 | 0.00072 |
| ... 1.2 | 0.993234 / 0.993558 | 0.00071 |
| the pulse at full weight up to its model's point (falling after) | 0.993212 / 0.993525 | 0.00076 |
| no regenerated phase in a pulsed band | 0.993002 / 0.993186 | 0.00118 |
| half of it | 0.993177 / 0.993549 | 0.00079 |
| its negative | 0.992765 / 0.992603 | 0.00176 |
| slope 1.000 (the time at b0 87 kept) | 0.993179 / 0.993545 | 0.00080 |
| slope 1.040 (likewise) | 0.993124 / 0.993523 | 0.00089 |
| every single-row pulse one sample earlier | 0.993037 / 0.993432 | 0.00099 |
| ... one sample later | 0.993108 / 0.993456 | 0.00090 |
| rows 24 and 25 on the single law (no grid of 11) | 0.992952 / 0.993088 | 0.00129 |
| the grid's origin 75.1 | 0.993189 / 0.993547 | 0.00077 |
| ... 76.1 | 0.993181 / 0.993582 | 0.00079 |
| the mid-frame model's pulsed bands from this frame's spectrum alone | 0.993156 / 0.993428 | 0.00090 |
| ... from the frame before's alone | 0.993199 / 0.993523 | 0.00079 |
| the reach below a pulsed channel 0 / 125 / 375 / 500 Hz instead of 250 | 0.993206 / 0.993534; 0.993221 / 0.993549; 0.993222 / 0.993567; 0.993205 / 0.993505 | 0.00075; 0.00073; 0.00072; 0.00077 |

The 40-sample window and the amplitude 1.2 tie with the rule within what these figures resolve; the
rule keeps the model's period and the voiced response's amplitude because they are what the voiced
response already is, not because they measured better. These alternatives were measured on the
scratch build only, on the 25 development and validation vectors (in sample: the rule was found on
them); the scratch build's `rs1` equals the library on those vectors to the digits printed (below).

**Lost frames.** In the eight channel-error vectors (hard and soft decisions) six frames are
repeats, not muted, of a frame of rows 18..31 as the baseline decodes them (none in the
soft-decision vectors): three of row 19 (`dam_e10_hd` 350, 1025 and its second repeat 1026) and
three of rows 26 and 27 (`dam_e5_hd` 840, 1239; `dam_e10_hd` 647). Residual over the frame's own
signal (the study's scratch build; the baseline's, the word dropped and the rule's columns are also
the library's builds', to the digits printed: the baseline, a scratch copy of the library whose
repeat drops the word, and the library):

| Frame | baseline | repeat plays the model as plain noise (the word dropped) | literal repeat (both masks and pulses as coded) | own mask in both halves, own pulse only | the rule (clause 7) |
|---|---|---|---|---|---|
| `dam_e10_hd` 350 (row 19, b0 76) | 0.196 | 0.171 | 0.149 | 0.123 | 0.092 |
| `dam_e10_hd` 1025 (row 19, b0 82) | 0.490 | 0.328 | 0.260 | 0.233 | 0.116 |
| `dam_e10_hd` 1026 (its second repeat) | 0.741 | 0.741 | 0.447 | 0.302 | 0.214 |
| `dam_e5_hd` 840 (row 27, b0 110) | 0.188 | 0.188 | 0.287 | 0.188 | 0.188 |
| `dam_e5_hd` 1239 (row 26, b0 96) | 0.068 | 0.068 | 0.099 | 0.068 | 0.068 |
| `dam_e10_hd` 647 (row 26, b0 73) | 0.145 | 0.145 | 0.156 | 0.145 | 0.145 |

The rule is the nearest in the three row-19 frames and ties with the plain repeat in the three
others, where a literal repeat is worse. Six frames of two files at bit error rates of 5 and 10 per
cent: this is what was observed, not an established property of DVSI's decoder. The frame after
`dam_e10_hd` 1026 (frame 1027, a muted repeat: its first 80 samples fade the last frame out) goes
from 1.022 to 0.225 (the library's builds; the scratch build's literal repeat goes as far): what it
shows is clause 6, the last own model's pulsed bands and pulse in the fade. On the tone vectors the
same loader (the repeat under a tone's first frame) changes, against that scratch copy whose repeat
drops the word (the library's builds; the study's scratch build gave the same), 35 frames of
`alltone` (0.0040 → 0.0040), 78 of `dtmf` (0.0008 → 0.0007), 43 of `dtmfvbad` (0.0186 → 0.0140), 2
of `alert` and one each of `cp0`, `cp1` and `dtmf60` (residual over the changed frames' own signal);
in the channel-error vectors only 5 frames of `dam_e10_hd` (0.6029 → 0.2248: the three row-19
repeats, frame 1027 and frame 351, the one after the repeat at 350), none with soft decisions.

**The library against the study's scratch build**, 16-bit samples, file for file, the 88 rate-33
vectors the scratch build decoded (25 development and validation, 52 tone, 4 hard- and 4
soft-decision channel-error and 3 DTX vectors; the soft-decision ones decoded as soft): 60 are
identical. 43 of them the change leaves as the branch point decodes them (every KNOX and dual-tone
vector, and the tone vectors `cp2`, `cp31`, `dtmf4`, `dtmf4r`, `dtmf8`, `dtmf8r`, `dtmf15n`,
`dtmf15p`, `dtmf25`, `dtmf50` and `dtmf100`); the other 17 the change moves alike in both (`alert`,
`dam`, `dam10`, `dam20`, `dam40`, `dam80`, `dam_e1_hd`, `dam_e1_sd`, `dam_e2_hd`, `dam_e2_sd`,
`damcar`, `dtmf180`, `fambf22c`, `noisy_dtx`, `t01`, `t03`, `tambf32e`). 24,776 samples differ in
the remaining 28, for three reasons, each checked:

- The scratch build computed the pulse times in float; the library forms them in double and rounds
  once. That moves samples by one unit in frames of rows 18..31, repeats of them and frames right
  after one with a pulsed own subframe (counted below).
- From silence (the first frame heard after a mute or a tone) the scratch build let the voiced train
  of its first subframe glide from the kept model's pitch as if that model were voiced (it kept the
  model's voicing from before `MBE_SYNTHESIS_FROM_SILENCE` silenced it); the library's train there is
  the branch point's, as clause 5 asks. A difference of this kind starts in a voiced frame right
  after a mute or a tone and reaches the frames after it that the train or the state carries it
  to. Of `dam_e10_hd`'s 133 frames that differ by more than one unit, 22 are voiced frames right
  after a mute, where it starts; it reaches 83 voiced frames after a voiced frame, 13 muted frames,
  7 repeats and 8 frames with no voiced band. `alltone`'s six such frames are two voiced frames
  right after a tone (229, 235), where it starts, and frames 230 and 236 to 238 after them, 238 a
  tone frame after a voiced frame (302 units: the speech its first 80 samples fade out); `dtmf`'s
  52 are 20 voiced frames right after a tone and 12 voiced frames, 12 tone frames and 8 frames with
  no voiced band after those; `dtmfvbad`'s four are two voiced frames right after a tone and the
  frame with no voiced band after each. A scratch copy of the library with that one
  behaviour switched in decodes all 88 files as the scratch build does but for 117 samples of 26
  files that differ by one unit, in 104 frames, each a frame of rows 18..31, a repeat of one or the
  frame right after one that pulses its own subframe (the first reason), and the 118 samples of
  frame 0 of the next.
- The scratch build kept its pulsed-row state in process-wide variables, which `mbe_quality_eval`'s
  second (int16) pass inherited from the end of its first: in the four files whose last speech frame
  pulses its own subframe (`fambm22a` row 25, `noisy` row 18, and `cp0` and `cp1`, row 20 before
  their closing tones), frame 0 of the written decode carries that frame's pulse (at most 47 units;
  in `cp0` and `cp1` it is the scratch build's third changed frame, the only one that differs from
  the library by more than a unit). The library's frame 0 there equals the baseline's.

No difference is left unexplained; none is a difference in the rule as written.

**How the library carries it.** The rate-33 decoder (`ambe2450_decode_vuv()`) writes, on a voice
frame of rows 18..31, a word in `PHIl[0]` in place of a mark: `1024 + 2048 held + 128 (row - 18) +
b0` as a float (all exact; `mbe_set_pulsed_row_word()`, `mbe_pulsed_row_word()` in
`src/internal/mbe_adaptive.h`): 2,640 values are words, the 1,680 the decoder writes (1024 to 2807)
and the 960 held words of rows 18..25 the repeat writes (3072 to 4087); a held form of rows 26..31
is none, as no writer gives one. `mbe_frame_mark()` reads a word as no mark, so that every reader of
marks treats such a frame as before (rows 18..31 carried no mark). The reader works on the bit
pattern (a range check on the bits, the conversion, the float rebuilt and compared by its bits, the
fields' ranges): any value no frame carries, a neighbouring float included, is no word. The public
`mbe_decodeAmbe2450Parms()` returns the word; the library's encoder never sends these rows. One
function in `src/core/mbe_subframe.c` (`current_row()`) decides whether the current model is a
pulsed-row frame: the rate-33 rendering, a word, no voiced band, a positive normal pitch, and the
band-grid model of the word's code (`L` equal to the code's harmonic count and the pitch equal by
bit pattern to the grid pitch the decoder computes, through `mbe_ambe2450_pulsed_row_model()` in
`src/internal/mbe_pulsed.h`, implemented in `ambe3600x2450.c` on the non-inlined function the
decode itself stores the grid model with). A word that fails plays as no word. The word is read in
one other place, the repeat loader (clause 7), which applies the same grid test (pitch and harmonic
count, not voicing) before it writes the held form, so that a caller's word on a model off the grid
is copied as the branch point copied it, and one on the grid model is given the held form or
cleared whatever the model's voicing. The grid pitch
depends on L alone, so two codes of one L have the same model and the test cannot tell them apart: a
word of another code of the model's L is taken (with that code's times). The last frame's own model
takes its row from the new slot of the subframe state (`noiseOverlap[37]`, `MBE_PULSED_ROW_SLOT`:
the word, not held, of a frame whose own model had a pulsed band, 0 otherwise; every other subframe
rendering, D-STAR's, P25's and ProVoice's, writes 0 there and reads it not, while the public
TIA-102.BABA rendering keeps its noise overlap in this element as in the others, reading and writing
it), read by its bit pattern, only in a state with the subframe mark, only a word as a
frame stores it (not held, rows 18..25), and only when the last frame's model passes the same test
with it; no mark of the previous model is read for any of it. A lost frame's repeat
(`ambe2450_load_repeat()`) writes the held form of a word of rows 18..25 and clears one of rows
26..31; the mute's fade toward a silent copy clears it with the mark, the unheard frames under a
mute and the kept model through a tone keep it where they keep a mark. The pulsed bands reuse the
pulsed state's `BAND_PULSED`, `pulse_bands()` (given a mask) and `prepare_render()`;
`pulsed_rows()` gives the three models their bands and times beside `pulsed_frame()` (which stays
P25's and D-STAR's and finds no code in the rate-33 rendering), `render_row_pulses()` plays each
pulse through its own model in the two intervals where the model has weight (with `render_pulse()`,
reached through a const volatile function pointer, below), and `mbe_advance_subframes()` stores the
slot and renders nothing.

**The pointer.** With two more direct calls of `render_pulse()` GCC stopped inlining it into the voiced
train (it cloned it instead), and in the 32-bit Release build (x87, `-std=gnu99`, so excess precision
kept as the registers allow) that moved the voiced samples of every path: on a scratch harness's
random streams (3,000 frames each of D-STAR's, P25's and rate 33's data paths and of ProVoice's,
rate 33's holding no frame of rows 18..31) 2,791, 2,846, 2,735 and 2,939 frames, by 1 to 74,898
ulps of the sample's own value (the large counts at samples near zero), at most 3 ulps of the
frame's peak (9.8e-4 at the most, on samples up to 4,447), and so, rarely, a 16-bit sample by one
unit: through the int16 entry points on the same streams five samples changed against the branch
point, D-STAR frame 1626 sample 131 (-567 to -566), P25 frames 1352 and 1471 sample 123, rate 33
frame 436 sample 23 and ProVoice frame 768 sample 156, and the head restores all five (the review
found the same five); the other nine configurations did not differ, in float or int16. Reached through the pointer the voiced
train compiles as before: the harness gives the branch point's output byte for byte in all ten
configurations ([November 2026: rate 33's pulsed
rows](#november-2026-rate-33s-pulsed-rows)). The C arithmetic is the same either way, but the rows'
own samples are not: over all 1,680 pairs of row and code, each between two frames of row 16, through
the three entry points (parameter bits, hard and soft frames, alike), the pointer moves them against
direct calls in the 32-bit Release build (654 pairs, 1,624 samples per entry point, all in the
pair's frame) and in dev-release, GCC with fast-math, LTO and `-march=native` (1,544 pairs, 9,206
samples, 3,046 of them in the frame after), by at most 3 ulps of the peak of the sample's own
160-sample frame (2.0371e-7 of it in the 32-bit Release build, row 21 code 41 sample 141 of the
pair's frame, and 2.2316e-7 in dev-release, row 19 code 59 sample 155; 1.31e-7 and 1.80e-7 of the
peak of the three frames together, the figure given here before), no 16-bit sample; the other eight
configurations, the 32-bit Debug build and both Clang ones among them, are
identical. (The review that found it counted 644 pairs and 1,714 samples on its own streams.) Only
the ten configurations of one x86-64 machine were checked; MSVC and macOS were not, and the
construction binds no future compiler.

**What is not known or not matched.** Rows 29, 24, 28 and 26 stay at 0.41, 0.29, 0.29 and 0.22 of
their own signal and row 20 at 0.19 (the others 0.07 to 0.18), where D-STAR's pulsed frames reach
0.03 to 0.05 ([November 2026: the pulsed state](#november-2026-the-pulsed-state), its later frames,
another selection); what DVSI's decoder plays beyond this rule on the frames of rows 18..31 whose
neighbours are not both voice frames of row 16 with no voiced band (0.24 of their own signal for
those whose two neighbours are voice frames with no voiced band, at least one of another row than
16, 30 of whose 2,687 have only rows 16 and 17 for neighbours; 0.22 for those with a neighbour that
is not a voice frame with no voiced band; against 0.08 for a frame between two frames of row 16);
why the change costs against the
input; what the decoder does on a lost frame beyond the six frames above (repeats of rows 19,
26 and 27 only: the held position of the other rows, of rows 24 and 25 above all, is the rule's
extension); what it plays on a frame of these rows that starts from silence, which the rule plays
as any frame of its row but that the last frame's own model plays nothing (clause 6), not studied
as a case of its own, no alternative tried (below); and whether an
AMBE-3000 treats other encoders' rows 18..31 alike (DVSI's own encodings only were compared; this
library's encoder never sends these rows).

**From silence.** By the library's own rule (the frame before left its start-afresh word: a muted
frame, a played tone frame or a tone repeat, read from the decoder's state before each frame), 822
voice frames (as defined above: no flag of `0x1F0`) of rows 18..31 with no voiced band, frames 1 to
n - 1 of each vector, start from silence in the scoreboard's 85 rate-33 vectors, 817 right after a tone frame in the tone
vectors and 5 right after a muted frame in the channel-error vectors, and 29 in the three DTX
vectors, each right after a muted frame. Residual over those frames' own signal, the baseline →
this change: after a tone frame 0.0433 → 0.0397, after a muted frame 0.4243 → 0.3145 (channel
errors) and 0.9703 → 0.4265 (DTX). A count by the flags of the frame before alone, not the
library's rule, adds frames that do not start from silence: 47 after an unmuted repeat (8 in the
channel-error vectors, 0.2944 → 0.2821, and 39 in DTX, 1.1852 → 0.2641) and 30 after a silence
frame in DTX (0.2823 → 0.1628).

**The pulses' fit**, measured on the library's decodes (the 22 vectors; "isolated one-pulse frames":
voice frames with no voiced band of rows 18..23 or 26..31 whose two neighbours are voice frames of
row 16 with no voiced band, where DVSI's signal in the row's pulsed channels over the model's
window, the 160 samples around its point, has an rms of 20 or more; 371 frames): the residual in
those channels over DVSI's energy there, by the frame's code, mid-frame pulses (rows 26..31) 0.098
(b0 0..19, 12 frames), 0.262 (20..39, 15), 0.167 (40..59, 27), 0.161 (60..79, 44), 0.087 (80..99, 46)
and 0.046 (100..119, 43); own pulses (rows 18..23) 0.761 (14), 0.513 (12), 0.520 (21), 0.178 (37),
0.180 (60) and 0.116 (40). At own codes 0..19, whose pulse time lies 80.89 to 61.586 samples before
its model's point (the response's nominal peak, 16 samples after the pulse, 54.7 samples before it
at code 10, the range's middle), playing no pulse there leaves less, 0.674. Each range's least-squares gain of the rule's
pulse against DVSI's estimated pulse (DVSI's decode less the library's without its pulse; a numpy
rendering of the rule, which matches the library's decodes to 0.0001 of the pulse's energy, from
two scratch decodes that differ only in the pulse's gain) is 1.42, 1.37, 1.89, 1.19, 1.22 and 1.19
for the mid-frame pulses and, for the own pulses, -0.43 at codes 0..19 (at a correlation of
-0.113: a poor fit), then 1.46, 1.04, 1.22, 1.03 and 1.04: DVSI's pulses are louder than the
rule's, the mid-frame ones most, and what DVSI plays for an own pulse at the lowest codes is not
this pulse.

**Onsets.** On the voiced frame after a voice frame with no voiced band (the 22 vectors, the onset
frame's DVSI rms 100 or more), the residual over the onset frame's own signal, by the row of the
frame before, is 0.0186 after row 16 (955 onsets), 0.0747 after rows 18..23 (258; 0.3430 before
the change), 0.2027 after rows 24 and 25 (111; 0.3334) and 0.0689 after rows 26..31 (263, unchanged):
the residual energy of the 632 onsets after a pulsed row is 0.00024 of all of DVSI's signal energy
of the 22 vectors, against 0.000055 for the 955 after row 16; as shares of all the residual energy
of the 22 vectors (every frame, lag 0; itself 0.0168 of the signal energy), 1.45% and 0.33%. Onset
frames of row 15 make most of them: 772, 213, 63 and 192 at 0.018, 0.075, 0.257 and 0.093. Over the
next three frames of the runs whose first four frames are voiced (the onset frame excluded), the lag
of ours against DVSI's stays near 0 in every class (medians +0.007, -0.006, -0.022, -0.008 sample):
those frames are in line, which does not exclude a timing part in the onset frame's own residual or
in a single run (`ucarm15a`'s onset at frame 1177, after a frame of row 18, fits a lag of -0.838
sample over its next three frames). After rows 24 and 25
one onset holds most of it: `lvl` frame 460 (row 15 after row 25 at b0 92) 58% of the class's
residual energy, and `dam_10ov` frame 460, the same speech, 15% more.

**The noise beside a pulsed own subframe.** In a pulsed own mask's channels, over the samples where the
next frame's mid-frame noise dominates, DVSI's noise is lower than the library's after a frame of
rows 18..23 when the frame is of row 15 (least-squares gain 0.834 over 143 frames; residual 0.0856
at gain 1, 0.0478 at that gain and 0.0470 at the square root of the pulsed frame's unvoiced factor)
and about the same when it is of row 16 (0.986, 157 frames) or when a frame of row 16 comes before
the pulsed frame (1.054, 297 frames). Neither the gain nor a reason for it was looked for further.

### Tones as harmonics of a fundamental

TIA-102.BABA-1 clause 7.3 ("Tone Regeneration") says only that the decoder synthesizes a tone
"using the decoded tone amplitude and frequency parameters", and Table 9 lists each tone index's
nominal frequencies. Annex J ("Tone Frame Parameters", normative) gives each index a fundamental f0
and harmonic numbers l1 and l2, and clause 8 (parametric rate conversion) turns a half-rate tone
frame into full-rate voice parameters with them (eqs 65-68: `w0 = 2 pi f0 / 8000`,
`L = floor(3812.5 / f0)`, harmonics l1 and l2 voiced at `16384 * 10^(0.03555 (AD - 127))`, every
other magnitude 0). US 7,970,606's Table 1 is the same table. The library played Table 9's nominal
frequencies, in rate 33 and in D-STAR.

**What DVSI's decoders play.** Measured on runs of one tone index of eight frames or more in DVSI's
decode of its tone vectors, each spectral line's frequency where a Hann-windowed transform of the
run's interior peaks (the same method on the library's decode before this change, which played
Table 9's nominal pair, gives 697.000 and 1209.000 Hz for index 129, DTMF "1"; the library now plays
693.909 + 1214.340 Hz):

- **Dual tones** (AMBE+2 indices 128..163; D-STAR's 128..147 map onto them). Both lines of each of
  the 36 tones are harmonics l1 and l2 of one fundamental. Measured on DVSI's rate-33 decodes, in
  units of 2^-19 cycle per sample (0.01526 Hz), the 36 fundamentals lie within 0.003 of an integer.
  34 equal `round(f0 * 2^19 / 8000)` for Annex J's printed f0, and two do not: index 145 plays 9883
  units (150.80 Hz; the text prints 150.89, which would be 9889) and index 147 plays 5667 (86.47 Hz;
  the text prints 86.5, 5669). D-STAR has 20 of the 36 (16 DTMF and four call-progress tones; it has
  no KNOX indices), and its decodes give the same 20 values. Examples (the played frequencies are the
  committed steps times 8000 / 2^19):

  | Index | Table 9 (Hz) | Played (Hz) |
  |---|---|---|
  | 128, DTMF "0" | 941 + 1336 | 942.078 + 1334.610 |
  | 129, DTMF "1" | 697 + 1209 | 693.909 + 1214.340 |
  | 137 | 852 + 1477 | 848.022 + 1484.039 |
  | 160 | 350 + 440 | 351.135 + 438.919 |
  | 161 | 440 + 480 | 424.988 + 495.819 |
  | 162 | 480 + 620 | 487.976 + 609.970 |
  | 163 | 350 + 490 | 350.037 + 490.051 |

  Over all components the played frequency lies -0.47% to +0.48% from the nominal one for DTMF
  (median magnitude 0.13%; the largest, index 137's 1477 Hz at 1484.04), -0.77% to +1.55% for KNOX
  (median 0.37%; index 153's 1279 Hz at 1298.83) and -3.41% to +3.30% for call progress (median
  0.97%; index 161's 440 Hz at 424.99).
- **Single tones** (5..122): harmonic `l = index / 16 + 1` (integer division) of the fundamental
  `floor(2048 * index / l)` units, which keeps the fundamental under 500 Hz. The played frequency
  then lies 0 to 0.092 Hz below the nominal 31.25 times the index, and on it for 69 of the 118. This
  is a rule of the rate-33 and D-STAR decoders; the P25 tone table's single tones use Annex J's
  harmonic numbers ([P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)). The
  evidence: (a) the three single tones with long runs (indices 30, 46 and 57: nine runs of 17 frames
  per index in `alert` and `dtmf`, in both codecs' decodes, and one short run in `alltone`) measure
  937.4998, 1437.470 and 1781.2501 Hz on DVSI's decoders (nominal 937.5, 1437.5 and 1781.25), so
  index 46 two units low; the rule plays the grid values 937.500, 1437.469 (3 times 31402 units) and
  1781.250 Hz. Index 46 is not an Annex J harmonic here: the rule gives it harmonic 3, where Annex
  J lists harmonic 4, which would play it at the nominal frequency exactly; (b) `alert` plays 30, 46 and 57 back to back, and the new tone's
  phase follows the pulse train it continues (below) only for harmonic numbers in the ratio
  2 : 3 : 4, within 0.0001 cycle, which fixes the ratio of the three harmonic numbers, not the
  numbers themselves; (c) over the 118 single tones, each tone's frequency estimated from its runs in
  `alert`, `dtmf` and `alltone` (most have only `alltone`'s one run of 3 to 5 frames) pooled with
  weights of the run length cubed, the rule puts 61% within one unit of 2^-20 of the estimate (rms
  1.4 units), against 36% (rms 4.1) for the nominal frequency, the best of about 330 rules of
  harmonic number and grid tried; those short runs carry a bias of a unit or two that is the same in
  both codecs' decodes. So the rule is well supported for indices 30, 46 and 57, and plausible, not
  proven, for the others.

**It is the voiced synthesizer.** DVSI's tones carry the marks of its voiced synthesis, US
8,200,497's impulse-response synthesis, which the library's [subframe
synthesis](#synthesis-in-subframes) models: weak lines at `k f0 + m 100 Hz`; images of harmonic l at
`(256 - l) f0` and `(256 + l) f0`, folded into the band, about `(l / 256)^2` down (index 129:
3195.80 and 2375.52 Hz at -62.1 and -62.9 dB for l 7, 3716.25 and 2895.94 Hz at -72.0 and -72.3 dB
for l 4, where `(7 / 256)^2` is -62.5 dB and `(4 / 256)^2` -72.2 dB), as the patent's K = 256 table
read with linear interpolation leaves them; and a ratio of the two components' amplitudes equal to
`sinc^2(pi l2 / 256) / sinc^2(pi l1 / 256)` within 0.0004 for all 36 dual tones (0.9807 for index
140's harmonics 12 and 23).

**How a tone starts.** After a frame with no voiced band DVSI's decoders start a tone as a fresh
pulse train at the first tone frame's first sample (after a voiced frame they do not: they continue
that frame's train, below, and the library starts the tone fresh there too, a known difference): a
component of frequency F plays
`cos(2 pi F (t - 16) + p)`, p the phase of the fixed output high-pass at F: since the low-band
change DVSI's third-order filter at 75 Hz (`mbe_low_band_dvsi_phase()`; [DVSI's low
band](#dvsis-low-band)), before it the second-order one at 100 Hz (`mbe_low_band_phase()`). The
phase at onset is the same at every onset of a tone within 0.0002 cycle (20 to 23 onsets per DTMF
tone), and the measured phase minus this formula was +0.0016 cycle on average, standard deviation
0.0022, over the 64 DTMF and KNOX components with the second-order filter (as measured for this
change); measured again with the low-band change over 68 components and 1,540 onsets, +0.00178 with
the second-order filter and +0.00035 with the third-order one ([DVSI's low band](#dvsis-low-band)). The first 80 samples ramp in, as the library already played them
([D-STAR tone transitions](#d-star-tone-transitions), [Rate 33 tone
transitions](#rate-33-tone-transitions)). A model of all this outside the library (steady tones
with a fresh start, no rendering of the glide from tone to tone) follows DVSI's rate-33 decode at a
correlation of 0.9995 to 1.0000 on the 32 one-tone vectors (`dtone_1` to `dtone_16`, `knox_1` to
`knox_16`).

**What the library does.** Rate-33 and D-STAR tone frames play these frequencies exactly:
`mbe_tone_played_steps()` gives an AMBE+2 tone ID's components in units of 2^-19 cycle per sample
(the 36 fundamentals with their l1 and l2, and the single-tone rule) and
`mbe_tone_dstar_played_steps()` a D-STAR index's, through D-STAR's mapping to the AMBE+2 tone of the
same nominal pair (144..147 to 160..163), both in `src/internal/mbe_tone.h`. Each oscillator's
32-bit phase moves by the component's units times 2^13 a sample, so the frequency is exact. A tone
whose first frame follows a frame that was not a tone frame starts at the phase above
(`mbe_startToneAmbe2450()`, `mbe_startToneDstar()`), whatever that frame was: after a frame with no
voiced band as DVSI's decoders start it, after a voiced frame where they continue the voiced frame's
train instead (below). A held tone counts only with a tone the codec plays (the test its decoder
applies to a tone frame's index); a word in a caller's state that names any other is no held tone,
and the next tone frame starts as above. A tone after a tone, a held tone, lost frames
that repeat a tone and a tone's tail go on from the phases the state carries, as before. Levels,
ramps, crossfades, held tones and lost frames are as they were. Table 9's nominal lookups
(`mbe_tone_lookup_freqs()`, `mbe_tone_lookup_dstar_freqs()`) are unchanged and still serve the
encoders, tone detection and the validity of an index. The public `mbe_synthesizeTonef()` and
`mbe_synthesizeTonefdstar()` play the new frequencies and continue the phases they are given: they
have no frame history, so they never start a tone afresh. No public signature or structure changes.
On the tone vectors' steady frames the median correlation with DVSI's waveform rises from 0.0014 to
0.9997 in rate 33 and from 0.0005 to 0.9995 in D-STAR ([November 2026: tones and the octave
step](#november-2026-tones-and-the-octave-step)).

**Not done: one pulse train for speech and tones.** DVSI's decoders play a tone through the voiced
synthesizer, and the library still plays it from two oscillators beside it. So:

- a tone after a voiced voice frame continues that frame's pulse train on DVSI's decoder, the pulse
  rate crossing linearly to the tone's fundamental over the tone's first 80 samples (`cp2` and
  `cp31`: predicted pulse phase 0.4798 and 0.8068 cycle, measured 0.4805 and 0.807; `cp0` and
  `cp1`, whose tones follow an unvoiced frame, start fresh), and starts fresh here;
- a tone after another tone continues the train on DVSI's decoder: of 134 changes from tone to tone
  in nine vectors (`dtmf4`, `dtmf4r`, `dtmf8`, `dtmf8r`, `dtmf15n`, `dtmf15p`, `dtmf25`, `dtmf50`,
  `dtmf100`), the same in both codecs, 66 follow the plain linear crossing, 53 the octave up and 10
  the octave down of [The pulse train across an octave](#the-pulse-train-across-an-octave), each
  within 0.0012 cycle of the prediction, and 5 cannot be told (index 138, whose harmonics 6 and 14
  leave the pulse phase open by half a cycle); here the oscillators go on from their phases;
- a voiced frame after a tone continues the tone's train on DVSI's decoders. Here, in rate 33 it
  starts from silence with a new train (`dtmf`: 40 voiced frames right after a tone, the second half
  of the frame at a median correlation of -0.03 with DVSI's on `fac0111`; not measured again, since
  nothing in this change touches it); in D-STAR, when the frame before the tone was heard (muted
  neither by the error level nor by the repeat limit, `ambe2400_last_frame_audible()`), that frame
  is played once more under the tone's first frame and leaves its pulse train in the last
  synthesized frame, the tone frames after it leave that frame as it is, and a voiced frame after
  the tone continues that train, as before the tone model change; when the frame before the tone was
  not heard, no
  repeat plays under the tone, and the voiced frame after it goes on from what the last synthesized
  frame held;
- the two components of a dual tone differ by the sinc^2 droop on DVSI's decoder (up to 0.17 dB,
  index 140), and play at one level here.

12 of the 50 rate-33 tone vectors with three steady frames or more are cases of the first two.

### The pulse train across an octave

Through voiced speech the pulse train of [Synthesis in subframes](#synthesis-in-subframes) moves its
rate linearly from one model's pitch to the next's across a subframe. Where the step is about an
octave DVSI's decoders do not glide. With R the new pitch period over the old: at 1.6 < R < 2.4 the
old train's phase (in [0, 1)) and its rate are halved first, and at 0.4 < R < 0.6 its phase is
doubled (modulo 1) and its rate doubled; the linear crossing then follows as at any step.
TIA-102.BABA and TIA-102.BABA-1 define no pulse train, so this is not a departure from their text.

**On tones.** Of the 134 changes from tone to tone of [Tones as harmonics of a
fundamental](#tones-as-harmonics-of-a-fundamental), where DVSI's decoder carries the train from one
tone to the next, 129 are classified, each within 0.0012 cycle of its prediction: the 53 at R 0.528
to 0.584 follow the doubled train, the 10 at R 1.712 the halved one, and the 66 at R 0.73 to 1.353
and 2.475 to 2.98 the plain crossing. The other 5 cannot be told (index 138, whose harmonics 6 and
14 leave the pulse phase open by half a cycle).

**On speech**, `fac0111`'s decode of the 25 development and validation vectors against DVSI's:

- **Onsets.** Of the D-STAR voiced runs that start after a frame with no voiced band (three voiced
  frames or more, loud; out of line where the two frames after the first correlate with DVSI's under
  0.8), 25 of 1,945 are out of line, 20 of them after a pitch-carrying code (11 of 421 after codes
  that pulse the frame's own subframe, 9 of 69 after codes that pulse both subframes; 3 of 1,044
  after a plain unvoiced frame). After a pitch-carrying code the first subframe runs at the
  geometric mean of the carried and the new pitch ([Voiced onsets and
  offsets](#voiced-onsets-and-offsets)), so the second subframe steps by the square root of the new
  period over the carried one. By that ratio, of 490 such onsets: 0.50-0.55, 1 of 4 out of line;
  0.55-0.60, 18 of 21; 0.60-0.65, 0 of 39; 0.65-1.6, 1 of 426.
- **Inside voiced runs** (frames in line over the two frames before; the mid-frame model sits at
  the geometric mean, so each subframe steps by the square root of the frame-to-frame ratio). The
  frames whose frame-to-frame pitch ratio is 0.30 to 0.42 (per-subframe period ratios 1.62, 1.62,
  1.62 and 1.72) are one in D-STAR, three in rate 33 and none in P25, and all four lose the
  alignment (correlation of the two frames after: 0.65; 0.40, 0.38 and -0.23). The frames nearest
  the two windows on the outside, at per-subframe period ratios of 1.534 and 0.611, keep it. Losses
  of alignment at other steps exist and this rule does not explain them: on `fac0111`, at
  frame-to-frame pitch ratios of 0.425 to 2.68, 31 D-STAR, 71 P25 and 96 rate-33 frames that were in
  line lose the alignment by this section's criterion (counted for the review of the tone model
  change).

**The edges.** The inner ones are bracketed, not measured: a per-subframe period ratio of 1.534
glides and 1.618 halves; 0.611 glides and 0.584 doubles. Nothing brackets the outer ones: 2.475
glides, but no step was observed between the largest halved one (1.72) and 2.475; the smallest step
observed to double is 0.528 (tones), and the four D-STAR onsets at 0.50 to 0.55 do not tell doubling
from gliding (the two end at nearly the same phase there). The library takes 20% on either side of
2 and of 1/2 (`octave_step()` in `src/core/mbe_subframe.c`): the outer edges, 0.4 and 2.4, are its
choice, by symmetry.

**The library** takes the step wherever both models of an interval are voiced, in every codec
path's subframe synthesis, in rendering and in `mbe_advance_subframes()` alike. The train keeps its
pulses' launch times, not their models: a pulse launched before a step is rendered, as before this
change, through both end models of each later interval it reaches; the step changes only the phase
and the rate the interval starts from.
The public `mbe_synthesizeSpeechf()` has no pulse train and is unchanged.

**Effect.** D-STAR's onsets out of line fall from 25 to 4 of 1,945 (after codes that pulse the own
subframe 11 to 0, after codes that pulse both 9 to 0, after a plain unvoiced frame 3 to 2, after
codes that pulse the mid-frame subframe 2 to 2). The four frames inside voiced runs come back in
line (D-STAR 0.647 to 0.939; rate 33, median, 0.383 to 0.873). P25 has no pitch-carrying codes; 60
of its speech frames change by this rule alone, and their residual against DVSI's, as a share of
all the vectors' signal energy, falls from 0.00176 to 0.00127.

### Amplitude limit

TIA-102.BABA eqs 115-116 scale a frame's amplitudes down whenever they sum past a threshold, 20480
with error-free bits. DVSI's P25 and rate-33 decoders do not: on identical bits, the frames that eq
116 scales (1.7% of rate 33's loud voice frames on the development vectors, 1% of P25's) play
within 0.01-0.13 dB (medians) of their unscaled level on DVSI's decoder, where eq 116 left them 2.95
and 1.05 dB quieter. Those paths, ProVoice's with P25's, now skip the scaling; eqs 111-114 (the
smoothed energy and the voicing threshold under errors) are unchanged. As a numeric guard they
still scale a frame whose amplitudes sum past 10^9, about 94 dB over 20480 and far beyond
full-scale output. D-STAR keeps it: there its
frames sit 0.22 dB under DVSI's with the scaling and 0.52 dB over it without, and the loudest
D-STAR frames already play louder than a DVstick 30 (see [Known differences](#known-differences)).
The public `mbe_synthesizeSpeechf()` keeps eq 116.

On the scoreboard, against the commit before it: the level offset from DVSI's decoder falls from
a mean of 0.332 to 0.065 dB on rate-33 development (`lvl` -3.55 to +0.22 dB, `dam_10ov` -2.02 to
-0.02) and from 0.195 to 0.131 on P25; log-spectral distance against DVSI falls 0.082 and 0.016 dB
(development), 0.009 and 0.012 (validation). D-STAR's output, every encoder's bits, the channel and
tone results are unchanged. The declared gate that rate 33's level move by at most 0.2 dB fails, on
that move toward DVSI (`lvl`).

### Encoder voicing over background noise

A recording encoded by this library's D-STAR encoder and played on DVSI hardware buzzed at the
lowest pitch in its pauses ([#149](https://github.com/arancormonk/mbelib-neo/issues/149)), where
DVSI's own encoding plays noise.

**Cause.** The pauses hold 60 Hz mains hum (about 75 LSB, -53 dBFS, with weak 120 and 180 Hz
harmonics) over about 13 LSB of background. In them all three encoders sent the lowest pitch with
every band voiced: 62 voiced pause frames in D-STAR, all 56 harmonics from 0 to 0.4 s. Two things
combine:

- The hum makes the error function E(P) small at the longest periods.
- At such a pitch a harmonic band spans about two bins of the 256-point spectrum, and fitting one
  window-shaped line there (eq 25) leaves noise little residual. White noise's expected fit error,
  computed from `w_R`'s spectral covariance and checked by simulation, is about 0.15 at 62 Hz, 0.3
  near 90 Hz and 0.7 at 250 Hz, against an eq (37) threshold near 0.4. So the noise bands passed
  as voiced.

It is not specific to hum. On DVSI's development and validation vectors, DVSI's encoder sent 81%
of the frames the D-STAR encoder voiced below 70 Hz with every band unvoiced.

**What DVSI's encoders do.** DVSI documents noise suppression in its encoders, on by default, and
its vectors were made with it. Its effect shows in DVSI's bits decoded by this library's decoder:

- In the pauses of `dam10` and `damcar` (`dam` plus white-like noise at -34 dBFS and car noise at
  -46 dBFS), the known noise comes out 10-12 dB down from 1 kHz and 23-29 dB down below 500 Hz.
- In `t01`/`t02`, whose pauses carry 60.1 Hz hum, the 60 Hz line is 50 dB down.
- On the final partition's `mark` the pauses come down only about 6 dB.
- Steady pulse trains (-20 to -39 dBFS) keep their level.

How much of the low-frequency loss comes from suppression and how much from DVSI sending most fully
unvoiced frames at one pitch (250.3 Hz D-STAR, 248.1 Hz P25, which leaves nothing below 125 Hz)
cannot be separated from the bits.

**This change.** The analysis estimates the background noise and uses it in the voicing decisions
only; magnitudes are unchanged (`src/core/mbe_speech_analysis.c`, "Background noise").

- **Noise estimate.** Per spectrum bin: the lowest smoothed power of the last 48 or so noise-like
  frames, doubled and at most their mean (after minimum statistics, R. Martin, IEEE Trans. Speech
  Audio Process. 9(5), 2001). A frame is noise-like when all of these hold:
  - it does not repeat in the voice band: E at the tracked pitch, on the input high-passed at
    100 Hz by a 6th-order Butterworth, is above the periodicity gate;
  - its power is within 6 dB of the previous frame's;
  - it is quiet (under -30 dBFS, where DVSI's encoder treats noise as noise from the start of a
    stream) or 15 dB under the loudest sound tracked so far (`xi_max` of eq 41); background lies 18 dB
    and more under the speech in DVSI's noisy vectors.

  The estimate is forgotten 50 frames (one second) after the last noise-like frame.
- **Voicing.** A 500 Hz column is treated as unvoiced before the codec quantizes its voicing when
  its energy is under four times the noise estimate and white noise's expected fit error there is
  under 0.3. That is mostly at pitches below about 90 Hz, where noise passes the voicing test;
  the expected error is weighted by the column's energy, so a column dominated by one low harmonic
  can fall under it at somewhat higher pitches (up to about 115 Hz).
  - D-STAR applies this to its own columns, P25 and rate 33 to the column of each band's middle
    harmonic.
  - D-STAR's 1 kHz bit can still voice such a column when the other column outweighs it.
- **What was tried and dropped.** Bringing the magnitudes down as DVSI does, with a floor fitted on
  `dam10` and `damcar`, made the held-out `mark` 0.4-0.5 dB further from DVSI's encoding: there
  DVSI attenuates far less, and low-level speech came out about 3 dB too quiet.

**Gates.** The `noise-voicing` gates were declared before scoring this design, and the final
partition had already been scored once for that suppressing variant, so it is not held out here.

**Results** (joint lags, against main; both encodings decoded by the same decoder):

- "Voicing agreement" is the scoreboard's per-1 kHz-band agreement with DVSI's encoding.
- "False voicing" is the share of DVSI's fully unvoiced frames that we voice.
- "Recall" is the share of DVSI's voiced frames that we voice.
- "Band recall" is the share of DVSI's voiced 1 kHz bands that we voice.

| Codec | Partition | Voicing agreement | False voicing | Recall | Band recall | LSD vs DVSI's encoding (dB) |
|---|---|---|---|---|---|---|
| D-STAR | development | 0.865 → 0.881 | 11.0% → 6.1% | 0.826 → 0.825 | 0.794 → 0.792 | 5.733 → 5.731 |
| D-STAR | validation | 0.821 → 0.853 | 23.0% → 11.4% | 0.824 → 0.823 | 0.796 → 0.795 | 4.612 → 4.613 |
| D-STAR | final | 0.837 → 0.840 | 11.9% → 7.5% | 0.783 → 0.781 | 0.767 → 0.756 | 5.794 → 5.802 |
| P25 | development | 0.864 → 0.879 | 11.7% → 7.0% | 0.820 → 0.820 | 0.795 → 0.790 | 5.669 → 5.665 |
| P25 | validation | 0.828 → 0.855 | 24.0% → 13.1% | 0.825 → 0.823 | 0.803 → 0.801 | 4.456 → 4.456 |
| P25 | final | 0.842 → 0.844 | 12.6% → 8.1% | 0.787 → 0.787 | 0.770 → 0.756 | 5.486 → 5.494 |
| Rate 33 | development | 0.844 → 0.860 | 14.4% → 8.5% | 0.791 → 0.790 | 0.795 → 0.790 | 5.716 → 5.711 |
| Rate 33 | validation | 0.832 → 0.860 | 23.0% → 9.4% | 0.794 → 0.789 | 0.799 → 0.795 | 4.508 → 4.508 |
| Rate 33 | final | 0.832 → 0.835 | 11.5% → 6.8% | 0.749 → 0.744 | 0.764 → 0.750 | 5.682 → 5.691 |

Levels against DVSI's encoding are unchanged (within 0.04 dB), and so is the pitch. On `t01` and
`t02` D-STAR's false voicing falls from 0.31 and 0.42 to under 0.1. On #149's recording the voiced
pause frames fall from 62 to 1 (D-STAR), 3 (P25) and 0 (rate 33).

Of the 121 checks, 16 fail:

- **Band recall on single vectors.** 2-3 points of DVSI's voiced 1 kHz bands are lost on `mark`
  (all modes) and on `ucarm15a` and `ucarf15a` (P25, rate 33), against a 2-point limit; on average
  0.4-1.4 points.
- **Low-voice recall over hum.** DVSI-voiced frames below 80 Hz over the hum of `t02` and `t03`
  (rate 33 and D-STAR): 15 of 34 and 4 of 12 in rate 33, 1 of 9 in D-STAR. DVSI voices some
  frames there that are low voice over hum.
- **Rate-33 voiced recall on single vectors.** 1.5 points on `t04` and 1.0 on `mark`, against a
  1-point limit.
- **Clipping.** More clipped samples on inputs that are overloaded or reach the rails: P25
  `dam_10ov` +3, `lvl` +1 and `t04` +3, and rate-33 `t04` +2.

Variants that changed the periodicity gate instead were rejected for real low voices:

- Cutting the top period to 122, with E taken from the high-passed input, never voiced a steady
  62.5 Hz voice and voiced 776 of DVSI's 816 sub-80 Hz voiced frames that main voiced.
- A gate passing on either a high-passed or a 50/60 Hz-notched input kept low voices, but the
  notches take 0.3 s to settle, so every stream still began with a buzz over hum.

Still different:

- Pauses keep their level: DVSI's encoders bring background noise down, this library's do not
  (a median 9.2 dB between the two encodings in the quietest quarter of a vector's frames, see
  [Known differences](#known-differences)). A model fitted to DVSI's suppression of added noise did
  not carry over to its real recordings. A second attempt in November 2026, a per-band gain over
  a tracked noise floor, did no better: its pauses stayed 3.1-7.9 dB from DVSI's encodings on the
  development vectors and 4.0-11.5 dB on the validation vectors, and it took `dam40` and `dam80`
  down 13 dB, where DVSI's encoding lies no lower than an unsuppressed one.
- A quiet voice whose harmonics are weak beside its fundamental, low and in noise, is periodic in
  the voice band only faintly and can lose voicing like hum. A synthetic 66 Hz voice at about
  -50 dBFS with harmonics falling as 1/h², at about 10 dB SNR, keeps voicing in 33-43% of frames,
  against all of them before. Loud (-20 dBFS), such a voice keeps all of them, also after a vowel
  13 dB louder, since only quiet sound or sound 15 dB under the loudest so far is taken for noise.
  Right after a sound 17 dB louder (a tone near full scale), it loses voicing for about a second,
  until that sound's level has decayed in the tracker or the estimate expires.
- After noise that stops under an unbroken vowel, the estimate takes 50 frames (one second) to
  go. Until then a column of that vowel whose noise-expected fit error is under 0.3 (mostly
  below about 90 Hz, in places up to about 115 Hz) can lose voicing.

### Encoder pitch code of unvoiced frames

A voice frame with no voiced band still carries a pitch code, and DVSI's D-STAR and P25 decoders,
and this library's, read it as side information ([Voiced onsets and
offsets](#voiced-onsets-and-offsets)): D-STAR b0 64..119 and P25 b0 128..183 end the voicing of the
frame before where the frame starts (and D-STAR's 64..71, 88..103 and 112..119 carry their pitch
into a voiced frame after them), any other code lets that voicing last through the frame's first
10 ms if it reaches 505 Hz. The encoders sent the tracked pitch on such frames, which on noise is
arbitrary, so what a DVSI decoder read there was arbitrary too, and so was the frame's harmonic grid
(9 to 56 bands for the same noise). On 17 clean development and validation vectors, the share of
frames with no voiced band at DVSI's default code, inside DVSI's range and elsewhere:

| Codec | Default code | Range | DVSI's encoder | This library's encoder (2.6.11) |
|---|---|---|---|---|
| D-STAR | b0 34 | 64..119 | 65.7% / 34.3% / 0.0% | 0.4% / 46.7% / 52.9% |
| P25 | b0 25 | 128..183 | 66.0% / 34.0% / 0.0% | 0.3% / 21.0% / 78.7% |

Rate 33 has no such side channel (V/UV row 16 ends the voicing where the frame starts), but its
grid was arbitrary in the same way.

**The rule.** The voicing is decided as before, at the tracked pitch. A frame whose result has no
voiced band is then sent at one fixed pitch code, and its amplitudes are measured again from the
same frame's spectrum on that code's harmonic grid (`mbe_analysis_measure_at()`), then quantized
against the decoder's history as before:

| Codec | No voiced band | Code | Bands | Amplitudes |
|---|---|---|---|---|
| D-STAR | b1 = 0 | b0 123 | 56 | voicing-independent estimate |
| Rate 33 | V/UV row 16 | b0 119 | 56 | eq 44 (unvoiced) |
| P25 | every bit of b1 0 | b0 127 | 37 | eq 44 (unvoiced) |

TIA-102.BABA-1 clause 4.1 codes every voice frame's own fundamental, the Annex A entry closest to
the estimate, and TIA-102.BABA 6.1 (eq 45) every frame's own pitch; neither treats a frame with no
voiced band differently. DVSI's encoders do: in D-STAR and P25 they never send the tracked pitch on
such a frame (table above), and in rate 33 they send b0 119 on 20.7% of their row-16 frames, more
than any other code. The encoders now depart from the text as DVSI's do, with codes of their own
([Departures from the specification](#departures-from-the-specification)). D-STAR b0 123 and P25 b0
127 lie outside DVSI's unvoiced ranges, so a DVSI decoder reads no carried pitch or energy pattern
into the frame, and every such frame says "the voicing before lasts through the first 10 ms" (if it
reaches 505 Hz), which is also what the codes outside the range sent before (53% and 79% of these
frames) said. That reading of D-STAR b0 120 to 123 rests on this library's decoder and on the
DVstick 30 measurement of [#157](https://github.com/arancormonk/mbelib-neo/issues/157) for other
codes outside the range, and that of P25 b0 127 on DVSI's decodes of its default code b0 25 ([Voiced
onsets and offsets](#voiced-onsets-and-offsets)); no P25 code has been measured on hardware. A
ThumbDV plays D-STAR b0 120 to 123 at b1 0 as low-pitch voice ([D-STAR b0 124](#d-star-b0-124),
"Checked on the hardware"), but in streams that start with them: whether b0 123 lets the voicing of
the frame before last through its first 10 ms has not been measured on hardware. Rate 33's b0 119 is
DVSI's most frequent code there, and the published AMBE-3000R measurements of the
[mbevocoder](https://github.com/srcejon/mbevocoder) project (`docs/encoder-quality.md`, the item
"Unvoiced frames on a fine grid") report that the AMBE-3000R codes fully unvoiced rate-33 frames on
at least 50 harmonics.

**The codes.** Chosen by a coarse sweep before this change, every unvoiced frame at one code: PESQ
against the input speech ([Method](#method)), mean paired change from the encoder before on the 25
development and validation vectors (and on the second set of 70 recordings where given; decoded by
this library's 2.6.11):

- D-STAR: b0 34 (DVSI's default, 14 bands) -0.025; 63 (22 bands) +0.037; 111 (47, inside DVSI's
  range, no carried pitch) +0.011; 119 (53, inside the range, pitch carried) -0.238; 122 (55)
  +0.073 / +0.090; 123 (56) +0.074.
- P25: b0 25 (DVSI's default, 14 bands) -0.079 / -0.085; 60 (23) +0.004 / +0.001; 100 (32) +0.038
  / +0.044; 127 (37) +0.058 / +0.050; 150 (inside the range) -0.025; 184 (51) +0.031 / +0.034; 207
  (56) +0.032 / +0.019.
- Rate 33: at least 50 bands +0.036 / +0.043; at least 56 bands +0.035.

DVSI's own default codes, on 14 bands, score below the tracked pitch's arbitrary grid, and a carried
pitch (D-STAR 119) far below it.

**Results.** PESQ and ESTOI ([Method](#method)) for this change against the commit before it, both
encodings decoded by this library's 2.6.11 decoder (this change touches no decoder):

| Codec | PESQ, 25 vectors | ESTOI, 25 vectors | PESQ, 70 recordings | ESTOI, 70 recordings |
|---|---|---|---|---|
| D-STAR | +0.074 ± 0.008 (24/25) | +0.0036 ± 0.0005 | +0.093 ± 0.011 (62/70) | +0.0064 ± 0.0012 |
| Rate 33 | +0.046 ± 0.006 (23/25) | +0.0030 ± 0.0007 | +0.055 ± 0.007 (59/70) | +0.0067 ± 0.0012 |
| P25 | +0.058 ± 0.010 (23/25) | +0.0039 ± 0.0020 | +0.050 ± 0.009 (53/70) | +0.0037 ± 0.0011 |

On the scoreboard (joint lags, against the commit before it; both runs' encodings decoded by this
change's decoder, which it does not touch). LSD is log-spectral distance; partition means:

| Codec | Partition | LSD vs DVSI's encoding (dB) | LSD vs input (dB) | Level gap to DVSI's encoding (dB) | Envelope corr. vs DVSI's encoding |
|---|---|---|---|---|---|
| D-STAR | development | 5.605 → 5.531 | 5.496 → 5.428 | 0.435 → 0.484 | 0.925 → 0.919 |
| D-STAR | validation | 4.377 → 4.309 | 4.689 → 4.623 | 0.398 → 0.455 | 0.938 → 0.936 |
| P25 | development | 5.494 → 5.443 | 5.313 → 5.314 | 0.229 → 0.241 | 0.939 → 0.935 |
| P25 | validation | 4.227 → 4.180 | 4.535 → 4.493 | 0.151 → 0.160 | 0.949 → 0.947 |
| Rate 33 | development | 5.670 → 5.598 | 5.307 → 5.265 | 0.212 → 0.190 | 0.923 → 0.921 |
| Rate 33 | validation | 4.347 → 4.289 | 4.715 → 4.666 | 0.132 → 0.130 | 0.938 → 0.935 |

Every encoding of every speech vector comes closer to DVSI's encoding in LSD (by 0.006 to 0.225
dB). Against the input, P25's development mean is unchanged (-0.06 to +0.30 dB per vector, `dam10`
the largest rise). The voicing decisions are unchanged on every vector, and so is the pitch and
voicing of every frame with a voiced band (none of 24,051 D-STAR, 24,506 rate-33 and 24,354 P25
voiced frames differs). Every frame with no voiced band is now at the fixed code (31,270 D-STAR,
30,815 rate-33 and 30,967 P25 frames on the 25 vectors), where the encoders before sent 50.2% of
D-STAR's and 24.9% of P25's inside DVSI's ranges. The scoreboard's encoder pitch error rises on
some vectors (D-STAR development median 6.5 to 10.1 cents, rate 33 9.9 to 17.2) only because it
measures frames whose pitch is steady over two frames either side, which now admits voiced frames
at the lowest pitch code (65 Hz, the fixed code's pitch) next to unvoiced ones (327 D-STAR and 302
rate-33 frames more).

The gate sets `unvoiced-code-development`, `-validation`, `-tones` and `-channel` were declared
before scoring and pass: no decode of DVSI's bits changes in any mode, the speech vectors' by their
samples' hashes and the tone and channel-error vectors' by theirs since a second scoring ([Scoring
history](#scoring-history)), on no tone vector does the encoder's tone agreement fall or its rate of
tone frames where DVSI sends voice rise (both unchanged on every vector; that rate is not zero, at
most 2.6%, `dtmfvbad`, in D-STAR and rate 33 alike), and the encoders' LSD from DVSI's encoding and
from the input does not rise by more than 0.01 dB nor the level gap by more than 0.1 dB (largest
mean rise: D-STAR's level gap, +0.058 dB on validation; on single vectors up to +0.17 dB, `noisy`).
The final partition was not scored for this change alone ([November 2026: encoders and rate-33 rows
18 to 31](#november-2026-encoders-and-rate-33-rows-18-to-31)).

Not done: choosing per frame between "ends at the frame's start" and "lasts to mid-frame" (the
codes inside DVSI's ranges), sending DVSI's default codes, and sending rate-33 rows 17 to 31.

### Joint choice of the AMBE spectral vectors

After the gain, the D-STAR and rate-33 encoders code the frame's prediction residual `T_l` (l = 1..L)
with six vectors: b3 and b4 (the PRBA vector, 512 and 128 rows) and b5 to b8 (one higher-order
vector per block). TIA-102.BABA-1 4.3.1 and 4.3.2 take each as the row of its table "closest (i.e.
has the minimum mean-squared error)" to its own target in the transform domain, and the encoders did
exactly that, each index on its own, every coefficient weighted alike; D-STAR, which has no public
specification, used the same search. Two things are lost that way: an error in a transform
coefficient lands on the harmonics unevenly (the blocks differ in length, and the decoder removes the
mean of the reconstructed residual), and an error on a spectral peak counts no more than the same
error in a valley, although the peaks are what is heard.

**The rule.** The nearest-row search stays as it was; the six indices are then refined together
against the error the decoder will produce (`src/ambe/ambe_enc_search.c`):

- E = sum over l of w_l (That_l - mean(That) - (T_l - mean(T)))², where That is the decoder's
  residual for the chosen rows (the decoder uses it less its mean; the gain carries the mean), and
  w_l = min(10, max(0.1, 2^((a_l - mean(a)) / 2))) with a_l the target log2 amplitudes: the square
  root of the harmonic's target amplitude relative to the frame's geometric mean, bounded.
- Two passes visit b3, b4, b5, b6, b7 and b8 in turn. At each visit, the other five held, the index
  becomes the row with the smallest E among those the decoder can receive (D-STAR's b8: every second
  row); the current row stays unless another is strictly better, and of equally good rows the lower
  index wins. A block of at most two harmonics has no higher-order vector and keeps index 0. E never
  rises.
- The decoder's residual is linear in 23 numbers (G2..G8 of the PRBA vector and C(i,3..6) of each
  block), so E is a quadratic in them, built from each block's weighted moments of its cosines, and
  one visit scores a row by a quadratic in at most four numbers. A row is scored by its difference
  from the current one, so a row equal to it there scores exactly 0 whatever the build's float
  evaluation order.

The refinement leaves the gain index b2, the pitch, the voicing and tone frames as they were, and
the IMBE encoder is unchanged; packing, `cur_mp` and the prediction history use the refined indices.
The recalibrated magnitude scales (Level, below) do move the gain: against the commit before, b2
differs on 42.0% of the D-STAR frames and 15.0% of the rate-33 frames of the 25 development and
validation vectors (9 to 53% and 4 to 20% per vector), and with the refinement alone on none. The
bits are ordinary D-STAR and rate-33 bits: the refinement only chooses among rows every decoder
receives, so they remain valid for any decoder.

**DVSI.** DVSI's re-quantizations of a decoded model pick the PRBA vectors this library's nearest-row
quantizer picks on 84 to 91% of frames ([Rate 33: field assignment of
u3](#rate-33-field-assignment-of-u3)), so DVSI's encoder appears to follow the specification here.
This change departs from both, for the gain measured below ([Departures from the
specification](#departures-from-the-specification)).

**The weighting.** A sweep of the exponent in w_l, measured on a prototype before [Encoder pitch
code of unvoiced frames](#encoder-pitch-code-of-unvoiced-frames) (PESQ against the input speech,
[Method](#method), mean paired change on the 25 development and validation vectors, D-STAR / rate
33): exponent 1/4 +0.020 / +0.026, 1/2 +0.031 / +0.041, 1 +0.013 / +0.027. With every weight 1 (the
decoder's unweighted error) the joint search does worse than the nearest rows: PESQ -0.010 / -0.014
then, and on this change's base -0.010 ± 0.007 / -0.009 ± 0.004 (1.5 and 2.0 standard errors below
no change), while ESTOI falls clearly, by 0.0011 ± 0.0001 / 0.0021 ± 0.0003, lower on 20 and 22 of
the 25 vectors. The weighting, not the joint search, makes the gain. Only the exponent 1/2 is used.

**Level.** Fitting the peaks better raises the decoded level against DVSI's encoding of the same
input, by 0.12 to 0.16 dB (D-STAR) and 0.13 to 0.17 dB (rate 33) on average over a partition. `AMBE2400_ENC_MAG_SCALE` and
`AMBE2450_ENC_MAG_SCALE` were therefore re-derived by their procedure, the median over the
development vectors of the level that puts our decoded encoding at DVSI's: D-STAR 1.02 to 1.075,
rate 33 1.08 to 1.06. D-STAR's encoding sat about 0.5 dB below DVSI's before this change (a median
-0.45 dB on the development vectors before the previous change, -0.49 dB after it); the scale now
takes that out as well. Level of our decoded encoding against DVSI's, median (mean) over each
partition's vectors, in dB:

| Codec | Partition | Before | Joint search | Recalibrated |
|---|---|---|---|---|
| D-STAR | development | -0.49 (-0.46) | -0.38 (-0.34) | -0.01 (+0.06) |
| D-STAR | validation | -0.46 (-0.46) | -0.27 (-0.29) | +0.16 (+0.14) |
| Rate 33 | development | -0.01 (+0.04) | +0.20 (+0.21) | -0.02 (+0.04) |
| Rate 33 | validation | -0.09 (-0.03) | +0.13 (+0.10) | -0.01 (-0.07) |

Neighbouring D-STAR scales put the development median at -0.04 dB (1.07) and +0.04 dB (1.08); 1.075
puts it nearest 0.

**Results.** PESQ and ESTOI ([Method](#method)) against the commit before this change, both
encodings decoded by this library's 2.6.11 decoder (this change touches no decoder):

| Codec | PESQ, 25 vectors | ESTOI, 25 vectors | PESQ, 70 recordings | ESTOI, 70 recordings |
|---|---|---|---|---|
| D-STAR | +0.040 ± 0.010 (23/25) | +0.0099 ± 0.0015 (24/25) | +0.038 ± 0.008 (51/70) | +0.0114 ± 0.0010 (67/70) |
| Rate 33 | +0.040 ± 0.008 (24/25) | +0.0142 ± 0.0020 (25/25) | +0.040 ± 0.009 (51/70) | +0.0149 ± 0.0009 (69/70) |

Before recalibration the search alone gave D-STAR +0.039 / +0.041 and rate 33 +0.042 / +0.045 PESQ
(25 vectors / 70 recordings), and the same ESTOI.

On the scoreboard (joint lags, against the commit before it; both runs' encodings decoded by the
same build), partition means:

| Codec | Partition | LSD vs DVSI's encoding (dB) | LSD vs input (dB) | Level gap to DVSI's encoding (dB) | Envelope corr. vs input | Envelope corr. vs DVSI's encoding |
|---|---|---|---|---|---|---|
| D-STAR | development | 5.529 → 5.503 | 5.423 → 5.434 | 0.484 → 0.184 | 0.920 → 0.923 | 0.919 → 0.921 |
| D-STAR | validation | 4.309 → 4.295 | 4.613 → 4.619 | 0.455 → 0.141 | 0.930 → 0.933 | 0.936 → 0.938 |
| Rate 33 | development | 5.577 → 5.527 | 5.262 → 5.274 | 0.190 → 0.188 | 0.930 → 0.935 | 0.921 → 0.925 |
| Rate 33 | validation | 4.266 → 4.226 | 4.661 → 4.660 | 0.130 → 0.123 | 0.930 → 0.933 | 0.934 → 0.937 |

The decoded encoding comes closer to DVSI's encoding (log-spectral distance down on every partition
and on 47 of the 50 encodings) and its envelope follows the input more closely (on 48 of the 50; the
other two lose at most 0.0022). Its log-spectral distance from the input rises slightly on the
development vectors, by 0.011 dB (D-STAR) and 0.012 dB (rate 33):
the search alone does that (0.009 and 0.014; the recalibration adds 0.001 and takes 0.001 off). The
weighting trades error on weak harmonics for error on strong ones, and log-spectral distance counts
every spectrum bin's dB alike. Split by the bin's level in the input frame, the RMS error of the
bins within 10 dB of the frame's strongest falls by 0.09 dB (D-STAR) and 0.10 dB (rate 33), while
that of the bins 20 dB and more below it, 62% of the bins measured, rises by 0.04 to 0.06 dB; with
every weight 1 neither moves. The rise is largest on `dam_10ov` and `lvl` (D-STAR, +0.13 and +0.12
dB); on the validation vectors the change is +0.006 and -0.001 dB.

The gate sets `joint-search-development`, `-validation`, `-tones` and `-channel` were declared
before scoring. `-validation`, `-tones` and `-channel` pass: no decode of DVSI's bits changes in any
mode (by the samples' hashes, those of the tone and channel-error vectors since a second scoring,
[Scoring history](#scoring-history)), the IMBE encoder's bits are unchanged on every vector, the
encoders' tone results are unchanged on every tone vector, and the D-STAR and rate-33 encodings'
log-spectral distance from DVSI's encoding and from the input, envelope correlation with the input
and level gap to DVSI's encoding all stay within their limits. `-development` fails on one check:
the log-spectral distance from the input rises by 0.0105 (D-STAR) and 0.0123 dB (rate 33) against a
limit of 0.01, for the reason above; its other checks pass. The final partition was not scored for
this change alone ([November 2026: encoders and rate-33 rows 18 to
31](#november-2026-encoders-and-rate-33-rows-18-to-31)).

**Cost.** The refinement adds 3.6 to 3.9 µs a frame in `dev-release` on the machine measured
(`bench_encode`, medians: D-STAR 19.8 to 23.6 µs, rate 33 20.9 to 24.8 µs).

## Departures from the specification

| Item | Specification text | DVSI's AMBE-3000 | This library |
|---|---|---|---|
| Rate-33 u3 bits 8-5 | Table 8: b3(0), b4(2), b4(1), b4(0) | b4(2), b4(1), b4(0), b3(0) | DVSI |
| Rate-33 pitch | Annex A table | smooth log2 law | DVSI |
| P25 pitch | eq 46: `w0 = 4 pi / (b0 + 39.5)`; the synthesis of eq 127-141 has no mid-frame model | eq 46's fundamental cut to a grid of 2^-19 cycle per sample, `floor(2^21 / (2 b0 + 79)) / 2^19` (0 to 2.26e-4 under eq 46), for speech and tones alike; between two voiced models a mid-frame pitch from a floored parabola for the square root of the product of the two fundamentals | DVSI ([P25 pitch](#p25-pitch)) |
| P25 tones sent as voice frames | TIA-102.BABA has no tone frames: every voice frame plays its decoded harmonics at eq 46's pitch | a voice frame whose pitch code and voicing word are those its encoder sends for a DTMF, KNOX or single tone, and whose tone harmonics dominate, plays a tone model: the tone's harmonics at one level, the others at least 27 dB down, eq 46's fundamental cut to a 2^-19 grid; call-progress tones, which DVSI's decoder does not play as equal-amplitude tone models, are left as decoded | DVSI, by a test with two measured limits ([P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)) |
| Rate-33 harmonic count, b0 17 | Annex A: 11 | 12 | DVSI |
| Rate-33 repeat | eq 57-58: `eps_0 >= 4`, or `eps_0 >= 2` and `eps_T >= 6` | `eps_0 >= 4` | DVSI |
| P25 output under channel errors | eq 112-116 | attenuates progressively | specification |
| Unvoiced synthesis (P25, and rate 33, whose synthesis TIA-102.BABA-1 takes from TIA-102.BABA) | eq 117-126: the 11.2 noise generator, one 256-sample window per frame | a 16-bit generator, `x <- 173 x + 13849`, and 10 ms subframes under a square-root triangular window | DVSI ([Unvoiced synthesis](#unvoiced-synthesis)) |
| Voiced synthesis (P25, rate 33) | eq 127-141: a frame in one piece; harmonics 1-7 at steady pitch interpolate linearly, the others crossfade with `wS`; the phase advances by `l 80 (w0(-1) + w0(0))` per frame, through unvoiced frames too | two 10 ms subframes around a mid-frame model, linear crossfades, a train of pitch pulses that starts afresh after an unvoiced subframe, harmonic phases regenerated from the decoded log2 amplitudes | DVSI ([Synthesis in subframes](#synthesis-in-subframes), [Phase regeneration](#phase-regeneration)) |
| Frames with no voiced band at a reserved pitch code (P25 b0 128..183; D-STAR b0 64..119, which no TIA text covers) | TIA-102.BABA eq 117-126: every unvoiced band plays noise, whatever the pitch code | each code names, for the frame's own subframe and its mid-frame subframe, pulsed 500 Hz channels, whose bands play the response to one pulse in place of noise | DVSI ([D-STAR and P25: the pulsed state](#d-star-and-p25-the-pulsed-state)) |
| Where voicing starts and ends inside a frame (D-STAR, P25, rate 33) | the frame's V/UV decisions hold for the whole frame | the pitch code of a frame with no voiced band (D-STAR, P25) moves it to mid-frame; in D-STAR and P25 a band whose state differs between two frames is voiced at mid-frame by 500 Hz columns read in both frames (at an onset only the voiced harmonics in or within 250 Hz below the columns this frame covers start with the frame, through a frame lacking a pitch only those the frame before covers hold); V/UV rows 13, 15 and 17 (rate 33) voice the first half of the frame otherwise than their end | DVSI ([Voiced onsets and offsets](#voiced-onsets-and-offsets), [D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns), [P25: the bands that start with the frame](#p25-the-bands-that-start-with-the-frame), [Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)) |
| Rate-33 tone frames | clause 7 gives the tone a frame carries, nothing on the frames around it | the tone crossfades with the speech over 80 samples at both ends | DVSI ([Rate 33 tone transitions](#rate-33-tone-transitions)) |
| Rate-33 tone frequencies | 7.3: the tone synthesized from its decoded frequency parameters; Table 9: the nominal frequencies | harmonics l1 and l2 of Annex J's fundamental on a grid of 2^-19 cycle per sample (indices 145 and 147 at fundamentals of its own), a single tone as harmonic `ID / 16 + 1` of a fundamental under 500 Hz; DTMF within 0.48% of nominal, call progress within 3.41% | DVSI ([Tones as harmonics of a fundamental](#tones-as-harmonics-of-a-fundamental)) |
| Rate-33 silence frames | 4.1 eq 2: L = 14; 4.3-4.4: amplitudes decoded as a voice frame's | 15 bands, amplitudes a further 13.8 dB down, the previous frame's model held to mid-frame | DVSI ([Rate 33 silence frames](#rate-33-silence-frames)) |
| Rate-33 mute | 5.7: bypass the synthesis, output noise uniform over [-5, 5] | noise of RMS about 0.6 (0, +-1 and a few +-2 at 16 bits); the first muted frame fades the last frame out over 80 samples; the model goes on underneath | DVSI's level, as 0 and +-1 ([Rate 33 mutes](#rate-33-mutes)) |
| Rate-33 tone ID 255 | 7.3 and Table 9: a tone of zero amplitude, whatever AD | a lost frame, as an unusable tone index | DVSI ([Rate 33 lost frames and tones](#rate-33-lost-frames-and-tones)) |
| Rate-33 frame repeat after a tone | 7.3 with 5.6: the previous frame's decoded model; silent on a tone before it | the tone repeats at its level until the repeat limit's mute fades it out | DVSI ([Rate 33 lost frames and tones](#rate-33-lost-frames-and-tones)) |
| Amplitude limit (P25, rate 33) | eqs 115-116 scale a frame whose amplitudes sum past 20480 (error-free) | no scaling | DVSI ([Amplitude limit](#amplitude-limit)) |
| Rate-33 V/UV vectors | Annex B rows and eq 5's columns, one voicing per frame; rows 17-31 all unvoiced | each voiced column reaches 250 Hz into an unvoiced one below; each row's own voicing, at the frame's end, is Annex B's; rows 13, 15 and 17 voice the first half of the frame otherwise (row 13 columns 0 and 1, row 15 none, row 17 column 0) | DVSI ([Rate 33: V/UV vectors](#rate-33-vuv-vectors), [Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17)) |
| Rate-33 V/UV rows 18-31 | Annex B: every band unvoiced, as row 16; 4.1 (Annex A): the fundamental of the frame's pitch code | on DVSI's encodings, a band grid whose L bands span 3.58 kHz, `w0 = 2.81 / L`, the pitch code setting only L; the row names pulsed 500 Hz channels of one subframe (rows 24 and 25: both), whose bands play one pulse each at a time the pitch code gives | the band grid and the pulsed rows as measured on DVSI's decodes ([Rate 33: V/UV vectors](#rate-33-vuv-vectors), [Rate 33: the pulsed rows](#rate-33-the-pulsed-rows)) |
| Encoder pitch code of a frame with no voiced band | TIA-102.BABA-1 4.1 and TIA-102.BABA 6.1 (eq 45): the code of the frame's own estimated fundamental, for every voice frame | D-STAR and P25: never the tracked pitch (a default code or one of 56 codes kept for such frames); rate 33: b0 119 on 20.7% of its row-16 frames | one fixed code per codec: D-STAR b0 123, rate 33 b0 119, P25 b0 127, amplitudes measured on its grid ([Encoder pitch code of unvoiced frames](#encoder-pitch-code-of-unvoiced-frames)) |
| Encoder spectral vectors (D-STAR, rate 33) | TIA-102.BABA-1 4.3.1 and 4.3.2: each of b3..b8 the row closest to its own target (minimum mean-squared error) | appears to follow the text (its re-quantizations pick this library's nearest PRBA rows on 84-91% of frames) | the nearest rows, then the six refined together to lower the amplitude-weighted error of the residual the decoder reconstructs ([Joint choice of the AMBE spectral vectors](#joint-choice-of-the-ambe-spectral-vectors)) |
| Encoder voicing over background noise | eq (37) decides voicing from the fit error and E(P) alone | suppresses background noise by default (documented); voices far less of it at low pitches | at pitches below about 90 Hz a 500 Hz column holding little more than the learned background noise is unvoiced ([Encoder voicing over background noise](#encoder-voicing-over-background-noise)) |

The P25 attenuation is not adopted: no rule tried (eq 112-116 as printed, and instantaneous,
clamped or leaky variants of it) explains it, and under channel errors our output stays closer to
the input speech than DVSI's (log-spectral distance from the input at 10% bit errors 10.4 against
12.1 dB).

## Known differences

- Soft-decision D-STAR streams are not muted on the error level, where DVSI mutes some of them
  (see above).
- Rate-33 V/UV rows 1 and 3: as measured for 2.6.7, DVSI's decoder plays some of their columns
  between voiced and unvoiced; no reading tried resolves them, and the library plays Annex B's
  rows ([Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17), where rows 13, 15 and 17, read so
  before, turned out to voice each half of the frame apart).
- Rate 33 inside voiced runs, with rows 13, 15 and 17 read as DVSI's decoder plays them and the
  pulse train running on through a row-17 frame. With the mid-frame pitch wave, a fit whose cause in DVSI's decoder is not known
  ([Rate 33: mid-frame pitch](#rate-33-mid-frame-pitch)), the residual inside the runs that stay in
  line no longer grows with the run's length on the development vectors; on the validation vectors
  it is lower in every group of frames into the run, and the group from frame 20 on is still the
  largest ([November 2026: rate 33's mid-frame pitch
  wave](#november-2026-rate-33s-mid-frame-pitch-wave)).
- What is left of rate 33's waveform difference with the pulse train running on through a row-17
  frame, in the accounting of [November 2026: mid-frame voicing](#november-2026-mid-frame-voicing):
  0.0163 of the signal on the development vectors, 0.0138 of it in the voiced runs that stay in
  line, 0.0000 in the runs under a correlation of 0.8 and 0.0025 in the frames with no voiced band,
  and 0.0163 on the validation vectors (0.0124, 0.0000 and 0.0039), whose largest group of frames
  into a run in line is that from frame 20 on (as measured with the mid-frame pitch wave); the parts
  are in [November 2026: the pulse train through a row-17
  frame](#november-2026-the-pulse-train-through-a-row-17-frame) and [November 2026: rate 33's
  mid-frame pitch wave](#november-2026-rate-33s-mid-frame-pitch-wave).
- Rate-33 V/UV rows 18 to 31. Since the rate-33 pulsed rows change they play rate 33's form of the
  pulsed state on their band grid ([Rate 33: the pulsed rows](#rate-33-the-pulsed-rows); before it,
  noise on the grid, [Rate 33: V/UV vectors](#rate-33-vuv-vectors); the study of the rows before the
  change is kept below). DVSI's encoder
  sends rows 17 to 31 on 12.2% of the voice frames of the speech vectors (row 16 on 37.3%). What is
  left: rows 29, 24, 28 and 26 stay at 0.41, 0.29, 0.29 and 0.22 of their own signal and row 20 at
  0.19 (the others 0.07 to 0.18; residual over the frames' own signal on the 22 vectors, the
  library's builds), where D-STAR's pulsed frames reach 0.03 to 0.05 (another selection); what
  DVSI's decoder plays beyond the rule on the frames of these rows whose neighbours are not both
  voice frames of row 16 with no voiced band (0.24 of their own signal where both neighbours are
  voice frames with no voiced band, at least one of another row than 16, a group that holds 30
  frames whose neighbours are rows 16 and 17 only; 0.22 where a neighbour is not such a voice
  frame; against 0.08 for a frame between two frames of row 16); why the change costs a little against the input (+0.008 and +0.016 dB of log-spectral
  distance from it on development and validation; PESQ against the input -0.0094 ± 0.0071, a 95%
  interval of -0.0250 to +0.0062 over the 12 recordings, which establishes neither a cost nor its
  absence); the pulses' level (DVSI's louder than the rule's, least-squares gains 1.19 to 1.89 for
  the mid-frame pulses, 1.03 to 1.46 for the own ones) and DVSI's own pulse at the lowest codes
  (rows 18..23 at b0 0..19 leave 0.761 of DVSI's in-channel energy, more than no pulse, 0.674, the
  least-squares gain there -0.43 at a correlation of -0.113);
  the onsets after a pulsed row (0.0747, 0.2027 and 0.0689 of the onset frame's own signal after
  rows 18..23, 24 and 25, and 26..31, against 0.0186 after row 16; the runs' next three frames
  stay in line, which does not exclude timing in the onset frame itself); what the
  decoder does on a lost frame beyond the six frames of the channel-error vectors the repeat rule
  rests on (rows 19, 26 and 27: the held position of the other rows, of rows 24 and 25 above all,
  is the rule's extension); what it plays on a frame of these rows that starts from silence (822
  such frames in the scoreboard's 85 rate-33 vectors, 817 after tone frames, and 29 in the DTX
  vectors; not studied as a case of its own, no alternative tried); and whether an AMBE-3000 treats
  other encoders' rows 18..31 alike (DVSI's own encodings only were compared; this library's encoder
  never sends these rows). The figures are in [Rate 33: the pulsed rows](#rate-33-the-pulsed-rows).
  The study of these rows before the change, on the band grid's decodes, kept as it was measured
  (what follows describes DVSI's output against noise on the grid, and its "not determined" is what
  the pulsed rows now model): putting DVSI's samples in place of ours on the frames of rows 17 to 31
  would have raised PESQ against the input by 0.019 ± 0.003 on the 25 development and validation
  vectors, nearly all of the 0.021 ± 0.004 by which DVSI's decode of the rate-33 bits scored above
  ours before the grid; the grid took about 0.003 of it. On plain voice frames (before the grid) our
  frames correlated with DVSI's 0.98 on row 16, 0.93 on row 17, 0.89 on rows 18, 19 and 21 to 23,
  0.87 on rows 26, 27, 30 and 31, 0.61 on rows 28 and 29, 0.75 on row 20 and 0.50 on rows 24 and 25.
  - **Where they act.** Each row replaces a fixed set of 500 Hz columns in one 10 ms subframe (rows
    24 and 25 in both, below) and leaves the other bands as row 16 plays them (coherence with ours
    0.93-0.97); in the affected ones of rows 18 to 23 and 26 to 31 (24 and 25 below) DVSI's signal
    has a coherence of 0.4-0.5 with ours and lies 1.5 to 8 dB below it within the subframe. Rows 18
    to 23 act on the frame's own (second) subframe, rows 26 to 31 on the mid-frame (first) one, and
    the two halves mirror each other; each partly affected band is the 250 Hz reach into the column
    below, as in the voicing rule:

    | Columns | Own subframe | Mid-frame subframe |
    |---|---|---|
    | 0 | 18 | 26 |
    | 0 and 1 | 19 | 29 |
    | 2 | 21 | 27 |
    | 7 | 22 | 30 |
    | 5 to 7 (own), 6 and 7 (mid) | 23 | 31 |
    | all, impulsive | 20 | 28 |

    Rows 24 and 25 act on both subframes, 24 on every band and 25 above about 750 Hz, as attenuated
    noise (1.5 to 13 dB down) that does not correlate with ours (coherence 0.05-0.19). Their level
    over the frame follows no law of time: rows 18 to 31 other than 20, 24, 25 and 28 play within
    about 1 dB of our level in every 20-sample block.
  - **What plays there is not determined.** It is not periodic at the coded pitch or at the grid's
    (normalized autocorrelation at the period 0.08-0.44, as for our noise in the same bands,
    0.09-0.48), and it is not DVSI's noise generator restarted from 0 at the subframe (correlation
    0.09 or less at offsets -3 to 3 samples).
  - **Rules tried and not adopted** (each ± here is over the vectors, [Method](#method)). Scaling
    the unvoiced bands of those columns in that subframe (rows 20, 24, 25 and 28 aside) brings our
    waveform closer to DVSI's (correlation +0.02 and log-spectral distance up to 0.24 dB lower on
    the frames concerned, 0.39 with the grid; on the scoreboard -0.005 to -0.007 dB alone and -0.013
    to -0.015 dB with the grid) but not closer to the input: by 0.5, PESQ -0.0013 ± 0.0014 and ESTOI
    -0.0003 ± 0.0001 (3 of 25 vectors better), log-spectral distance from the input +0.004 dB; by
    0.3, PESQ -0.0023 ± 0.0016; with the grid, the grid's PESQ gain falls from 0.0027 to 0.0007.
    Voicing those columns in that subframe instead is worse on every waveform measure: correlation
    with DVSI's -0.03 to -0.16 on these rows and -0.11 on rows 0 to 15 (the pulse trains it starts
    move later pulses), the scoreboard's log-spectral distance from DVSI's decode +0.07 to +0.09 dB
    and its waveform correlation -0.08 to -0.12; ESTOI falls by 0.0003 ± 0.0001 (8 of 25 better),
    while PESQ rises by 0.005 ± 0.004, within its error.
  - **The impulsive rows.** Row 20 plays a quiet second half (1.8 to 18 dB below ours by
    20-sample block) and then a burst about 108 samples plus one coded pitch period into the frame
    (43% of frames within 4 samples of it, against 7% by chance). Row 28 plays a quiet first half
    (4 to 13 dB down), with a peak related to the pitch on a minority of frames. Rows 24 and 25 are
    peaky (kurtosis about 8.5 against 3.8 for ours), and on a given pitch code row 25's peak often
    falls at the same place (60% of frames within 3 samples of their code's median); no single law
    was found.
- The pitch code of a P25 or D-STAR frame with no voiced band. Since the pulsed state ([D-STAR and
  P25: the pulsed state](#d-star-and-p25-the-pulsed-state)) the 56 codes kept for such frames (P25 b0
  128..183, D-STAR b0 64..119) play the patterns of pulsed 500 Hz channels they name, and the
  pitch code is no longer a side channel the library leaves unread. Measured before that state:
  frames at the encoders' default codes (P25 b0 25, D-STAR b0 34) were in close waveform agreement
  with DVSI's output, across level steps too (median correlation 0.999), and those at the 56 codes correlated
  0.87-0.98, where the level stepped up DVSI playing one half of such a frame, or the start of the
  next, 1 to 9 dB quieter than this library, which half and by how much going with the block and the
  code's low bits; a D-STAR frame of this kind right before a tone played its second half a median
  7.5 dB louder here than on DVSI's decoder. What the pulsed state leaves: P25's first frame with a
  code after a voiced frame, when its mid-frame subframe alone is pulsed, stays further from DVSI's
  decode than D-STAR's (on the library's builds, residual over those frames' own signal 0.177 over
  606 frames, D-STAR's 0.034 over 738; [D-STAR and P25: the pulsed
  state](#d-star-and-p25-the-pulsed-state)); two of the DVstick
  30's three steps into frames at `b0` 100 go further from the unit's curves; ProVoice, with no DVSI
  material, plays such frames as noise; and the mid-frame pitch between two default-code frames of
  different pitch is unmeasured and left as it was. A frame with no voiced band at one of the 56
  codes from another encoder plays its pattern here in D-STAR whenever its model qualifies and in
  P25 while the validity counter allows; what an AMBE-3000 does with such a stream is not known.
  This library's encoders send neither those codes nor the default ones on such a frame, but one
  fixed code outside the ranges ([Encoder pitch code of unvoiced frames](#encoder-pitch-code-of-unvoiced-frames)),
  so they do not use this channel: every frame with no voiced band tells a DVSI decoder, as far as
  the codes measured show, that the voicing before lasts through its first 10 ms (in the bands whose
  500 Hz column, or the column 250 Hz higher, that voicing covers, by the rule of [D-STAR and P25:
  mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns)), where DVSI's encoders
  end it at the frame's start on about a third of them, and plays as noise, where DVSI's encoders
  send a pattern on 28.5% (P25) and 28.8% (D-STAR) of them (8,087 of P25's 28,384 frames with no
  voiced band on the development and validation vectors, 8,222 of D-STAR's 28,555).
- P25 tones sent as voice frames ([P25 tones sent as voice
  frames](#p25-tones-sent-as-voice-frames)). The first frame or two of a tone can fail the test and
  play as speech: a first frame passes the dual-tone limit with as little as a third of its energy
  in the tone's harmonics, and a tone's first frames still differ from DVSI's more than its steady
  ones (a median residual of -13.0 dB on 1,254 first tone frames, where rate 33's reach -31.6 dB).
  Call-progress tones play as their decoded harmonics: DVSI's decoder does not play them as
  equal-amplitude tone models (0 of 50 selected frames), and the library leaves them as decoded. Of the single tones
  only `alert`'s 30, 46 and 57 have long runs to measure, and the single-tone limit rests on the
  isolated speech frames of the 25 development and validation vectors that pass a single-tone row
  (five below it, about 65 above). Between DVSI's phases and the regenerated ones mean offsets of
  0.00 to 0.07 cycle remain (median magnitude 0.015), largest for harmonics 4 and 7 at short
  periods.
- Rate-33 and D-STAR tones beside the synthesizer ([Tones as harmonics of a
  fundamental](#tones-as-harmonics-of-a-fundamental)). DVSI's decoders play a tone through the
  voiced synthesizer and the library from two oscillators beside it. A tone after a voiced frame
  takes its phase from that frame's pulse train on DVSI's decoder and starts fresh here (as DVSI's
  start one after a frame with no voiced band); a tone after another tone keeps the train there
  (plain or [across an octave](#the-pulse-train-across-an-octave)) and continues the oscillators
  here; a voiced frame after a tone continues the tone's train there, while here rate 33 starts it
  from silence with a new train (`dtmf`: the second half of the 40 voiced frames right after a tone
  correlates with DVSI's at a median -0.03 on `fac0111`; not measured again) and D-STAR, when the
  frame before the tone was heard (muted neither by the error level nor by the repeat limit), goes
  on from the speech pulse train it had before the tone (the repeat of that frame under the tone's
  first frame leaves it, and the tone frames keep it); and the two components of a dual tone differ by the
  sinc^2 droop there (up to 0.17 dB, index 140) and play at one level here. 12 of the 50 rate-33 tone vectors
  with three steady frames or more hold the first two cases, and they are the 12 whose steady
  frames' median correlation with DVSI's stays under 0.99. The single tones other than 30, 46 and 57
  rest on a rule that is plausible, not proven.
- The `sinc^2(pi l / 256)` droop on voiced speech. On steady, fully voiced speech frames the median
  gain of DVSI's decoder over the library's, by harmonic number, follows it in rate 33 up to l 29
  (-0.33 dB at l 25-29, against -0.31) and partly in D-STAR, not in P25, and a fit takes nothing off
  the 0.6 to 0.75 dB rms scatter. It is modelled for P25's tone frames only.
- DVSI's low band ([DVSI's low band](#dvsis-low-band)). The fundamental's level under the corner and
  the low harmonics' phase, which differed from DVSI's by up to 12 dB and 0.5 rad under 140 Hz before
  the change, are modelled; what is left: DVSI's harmonic 1 over ours at an rms of 0.600 dB over the
  81,342 readings (0.390 where the model predicts no deficit, the measure's own level; 0.65 to 1.95 dB
  where it predicts one, rising with the deficit), the most on the car-noise and filtered versions of
  one recording (`ucarm20a`, `fambm22a`, `ucarm15a`, `ucarf15a`, 0.9 to 1.2 dB in rate 33 and P25)
  and P25's `noisy` (1.20 dB, a mean of -0.50); D-STAR's pulse trains' harmonic-1 level, off by -0.9 to
  +1.0 dB on three of them with or without the cut, which the cut does not explain; a phase of DVSI's
  voiced harmonics over ours in P25 under 300 Hz, +0.02 to +0.05 rad at fundamentals over 75 Hz,
  that is not the filter's (its correlation alone would take a corner of 80 Hz); in the unvoiced bands near a high corner, +1.0 dB
  (P25) and +3.5 dB (D-STAR) more in DVSI's decode at 164 Hz than the library's whole bins let
  through; DVSI's DTX vectors, which come no nearer (0.17761 → 0.17784 D-STAR, 0.17401 → 0.17464
  rate 33, on the frames the change moves). Not known: the tracker's law beyond the fit (one rate in
  frequency, two in log frequency and a share of voiced bands fit alike), whether it steps once a
  frame or a subframe, the ceiling's exact value, what DVSI's decoders do with it on repeated, muted,
  silence and tone frames (the library's choices), whether they cut the pulsed state's bands, and
  ProVoice.
- P25 onsets the column rule does not fit ([P25: the bands that start with the
  frame](#p25-the-bands-that-start-with-the-frame)): the five onsets of the development and
  validation vectors that its whole-frame form does not explain (three that DVSI's decoder starts
  early with no column covered, two it starts late with a column covered) and the 14 onsets that
  neither start fits. Of the voiced runs of two frames or more under a correlation of 0.8 left on
  the development vectors with the rule (20 by this count; the accounting of [November 2026:
  mid-frame voicing](#november-2026-mid-frame-voicing) counts 36 its own way), three are runs of
  `noisy` whose onsets are among those 14 (the runs from frames 2795, 2620 and 1944; the lag of our
  decode against DVSI's that lines up the two frames after the onset is -4, +19 and +2 samples), and
  one is `dam40`'s 27-frame run from frame 758. Frame 758 is a voiced onset that the study does not
  select (the two frames after it are quiet, a standard deviation of 10.8 in DVSI's decode); the
  run is in line (frame correlation above 0.8) for its first 11 frames, falls to 0.26 and comes
  back to 0.80 (0.65 over the whole run): the loss is inside the run.
- P25 voiced runs out of line with DVSI's decode after a speech frame DVSI's decoder plays as a tone
  model: with the four rules of the P25 decoder change set, about half of P25's remaining whole-file
  residual on the study's scratch build lay in a dozen voiced runs out of line (`noisy` from frames
  2796 and 1946, `t04` from 684, ...); in the ones examined, DVSI's decoder plays the one speech frame
  just before the run as a tone model (two harmonics, the others and the noise silent: `noisy` 2795
  harmonics 7 and 11, `t04` 683 harmonics 5 and 7, `noisy` 1944 harmonics 4 and 5), at a pitch code
  one or two away from a tone's own, which the library's tone test (exact rows of code and voicing
  word) does not pass, and the frame after it takes its own pitch at mid-frame. The tone test's
  tolerance is a later study's.
- P25 frames between tone runs: the floored parabola between a tone model and the voice frame after
  it moves the pulse train's phase into the next tone, and eight DTMF vectors lose whole-file
  correlation with DVSI's decode (0.976-0.979 → 0.958-0.966); a geometric mean there fits six and
  not the other two ([P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)).
- P25's own model by columns ([D-STAR and P25: mid-frame voicing by
  columns](#d-star-and-p25-mid-frame-voicing-by-columns)): whether the column reading is DVSI's rule
  or the nearest the library's constants come to it (it leaves about a third of the bands'
  disagreements with the labels, about 150 of 9,491, spread over the columns), and what DVSI's
  decoder does with a model whose voiced bands all lie off the samples (channel errors' smoothing
  makes a few; the library plays it with no voiced band). Why P25's phases come from the half-rate's
  magnitude domain ([Phase regeneration](#phase-regeneration)): the simplest reading is one phase
  routine for all the decoder's modes; the step's size is the half-rate standard's constant,
  measured to about 0.1 in log2.
- A single band whose voicing changes between two frames, on a DVMEGA DVstick 30
  ([#157](https://github.com/arancormonk/mbelib-neo/issues/157)). With the rule of
  [D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns) the mean error of the 18 single-band steps' curves falls from 0.346 to 0.189
  ([Synthesis in subframes](#synthesis-in-subframes)). Still far from the unit's: the steps where
  band 2 turns voiced (1101 → 1111, a mean of 0.384 over three positions from 0.391, and at position
  34 worse, 0.426 → 0.478) and 1011 → 1111 (0.222 in both builds to three decimals, though the decodes
  differ inside the step frame). The curves are full-band rms in 16-sample blocks and do not show
  which band starts when.
- The mid-frame voicing of D-STAR and P25 by columns ([D-STAR and P25: mid-frame voicing by
  columns](#d-star-and-p25-mid-frame-voicing-by-columns)) is a fit to labels of DVSI's decodes
  (which of two test decodes is nearer DVSI's in each harmonic's band), with the DVstick 30's curves
  coming closer with it; what in DVSI's decoders it stands for is not known, nor why, in most
  labelled bands, the top column of a run of changing columns takes the unvoiced state even under a
  column voiced in both frames (1101 → 1111 column 5: 39 labelled early, 352 late; 1111 → 1101: 53
  held, 398 let go). Why P25's in-run agreement with the labels (0.886) is lower than D-STAR's (0.944) was
  not examined. ProVoice has no DVSI reference and takes P25's rule: its
  rendering (`MBE_VOICED_SYNTHESIS_DVSI_PROVOICE`) shares P25's column rule. A D-STAR voice frame after a tone, a mute or a resume takes
  the rule against the last model rendered, as it took this frame's voicing before; DVSI's vectors
  hold no tone inside speech.
- Encoder: pauses keep their level, where DVSI's encoders suppress background noise
  ([Encoder voicing over background noise](#encoder-voicing-over-background-noise)). In the
  quietest quarter of each speech vector's frames, our rate-33 encoding decoded, as measured for
  2.6.8, a median 9.2 dB above DVSI's encoding of the same input (0.2 to 13.6 dB per vector, `dam40`
  and `dam80` aside, whose input lies 40 and 80 dB down), and the loudest quarter sat within 2.0 dB
  on every vector but `dam10` (+4.2). The encoder changes of 2.6.12 move the level of the quietest
  quarter by a few tenths of a dB.
- DVSI hard-clips at full scale; this library clips at 0.95 of full scale. Clipping at full
  scale moved no metric, so it is not adopted.
- A frame at the largest D-STAR gain and b3 (b2 63, b3 511; `9E B0 6F EE F8 7B 9F FE 74`) held
  for 50 frames ([#158](https://github.com/arancormonk/mbelib-neo/issues/158)): a DVMEGA DVstick 30
  plays -18.4 dBFS on the first frame and -9.0 to -11.4 dBFS after it, this library -13.0 and -6.5
  to -8.4 dBFS. It used to play silence there from the 5th frame, when the spectral enhancement's
  weights overflowed float.
- After a D-STAR b0 124 frame (the null frame among them), a DVMEGA DVstick 30 plays the lost
  frames of a run from the 4th or 5th on at a comfort-noise level that follows the gain field b2
  of the b0 124 frames (up to -2 dBFS after b2 63, b3 511), and plays all lost frames at level
  directly after one ([#156](https://github.com/arancormonk/mbelib-neo/issues/156)); this library
  mutes from the 5th, at about -95 dBFS. Of those frames radios send only the null frame (b2 4).
  After it, on the measured passage, the unit's 4th lost frame follows the voice before the run, as
  ours does (within 6 dB), while the run starts up to 16 voice frames after the null frame; with 17
  to 20 voice frames between, it plays the 4th at -73.9 to -75.9 dBFS where ours plays -87.7 to
  -90.4, and a run of 30 (20 voice frames after) from the 5th on at about -74 dBFS (frame RMS about
  6) where ours plays about -95 (RMS about 0.6). It is not modelled: that would take the comfort
  noise's spectrum (not measured), a level law over b2 to b5 fitted to one sweep, a rule for how
  several b0 124 frames combine (not determined), and more hidden state in `mbe_parms`, whose layout
  is fixed.
- Pitch per code: P25 plays eq 46 on DVSI's 2^-19 grid, which leaves 0.00338% rms (weighted) of the
  per-code drift on the 49 rows measured, and between two voiced models DVSI's floored parabola
  ([P25 pitch](#p25-pitch)); frames at a pitch step keep a little more of the lag's change over a
  frame than frames at one code (0.0101 against 0.0086 samples at a step of one or two codes, 0.0185
  at seven or more), not explained. D-STAR within 0.013% rms (at most 0.04% on a code with 10 or more frame
  pairs) and rate 33 within 0.0074% rms (at most 0.02%), measured as in [P25 pitch](#p25-pitch). A
  law with log2 f0 rounded to a 2^-12 grid explains part of rate 33's difference (0.0074% to 0.0048%
  rms; its slope then comes out at 1 / 45.368) but gains 0.0007 of waveform correlation only, and a
  2^-11 grid for D-STAR measured worse; neither is adopted.
- What is left of the waveform difference. On the development vectors our decode of DVSI's bits
  still differs from DVSI's decode by a residual energy of 0.037 (D-STAR), 0.069 (P25) and 0.065
  (rate 33) of DVSI's, sample for sample at lag 0, pooled over the 20 vectors with each vector
  scaled to the same RMS level of DVSI's decode; these figures and their parts below are pooled so.
  The means over the vectors of each vector's own ratio are 0.0317, 0.0664 and 0.0561 ([November
  2026: phase regeneration and pitch](#november-2026-phase-regeneration-and-pitch)), and with the
  tone model change 0.0191, 0.0633 and 0.0527 ([November 2026: tones and the octave
  step](#november-2026-tones-and-the-octave-step)), and with the mid-frame voicing change 0.0191
  (D-STAR's decode does not change), 0.0309 and 0.0388 ([P25: the bands that start with the
  frame](#p25-the-bands-that-start-with-the-frame), [Rate 33: rows 13, 15 and
  17](#rate-33-rows-13-15-and-17)), and for rate 33 with its mid-frame pitch wave 0.0214 ([November
  2026: rate 33's mid-frame pitch wave](#november-2026-rate-33s-mid-frame-pitch-wave)), and with the
  mid-frame columns change 0.0172 (D-STAR) and 0.0273 (P25), rate 33's unchanged ([November 2026:
  mid-frame voicing by columns](#november-2026-mid-frame-voicing-by-columns), which also gives the
  accounting by kind of frame before and with that change). The pooled
  figures and their parts below were measured for
  PR #196's decoder changes, before the tone model change; an accounting by kind of frame measured
  before and with the mid-frame voicing change is in [November 2026: mid-frame
  voicing](#november-2026-mid-frame-voicing). Of it:
  - Voiced runs that line up with DVSI's (run correlation 0.8 or more) hold 0.018, 0.021 and 0.033.
    P25's no longer grow with the run's length (0.019 to 0.028 over frames 0 to 3, 4 to 7 and so
    on); D-STAR's are 0.016 to 0.018 over a run's first twenty frames and 0.037 past them; rate 33's
    still grew then (0.018, 0.024, 0.036, 0.053, 0.067 and 0.067); with the mid-frame pitch wave
    they no longer grow on the development vectors and are lower in every group on the validation
    vectors, the group from frame 20 on still the largest there (above).
  - Voiced runs that do not line up hold most of the rest: D-STAR 45 of 1911 runs, holding 0.017 of
    the signal energy as residual, 44 of them off from their first frames (19 after a frame whose
    pitch code is one of the 56 kept for frames with no voiced band and names a pattern for the
    frame's own subframe, 12 after one naming both subframes, 4 after one naming the mid-frame
    subframe and 9 after a plain unvoiced frame; counted for the tone model change on the
    development and
    validation vectors, three voiced frames or more, the octave step brings every out-of-line onset
    after a code that pulses the own subframe or both in line, [The pulse train across an
    octave](#the-pulse-train-across-an-octave)); P25 85 of 1959 runs holding 0.045, 83 off from
    their first frames, which no single lag lines up (of the 77 whose first three frames correlate
    under 0.8, 22 reach 0.85 at any lag within 130 samples); rate 33 43 of 1825 runs holding 0.030,
    31 of them off from their first frames (0.012) and 12 starting aligned and losing it inside the
    run (0.018; median length 21 frames).
  - Frames with no voiced band hold 0.003 in each codec (measured before the pulsed state; in that
    change's accounting, another selection, the frames with no voiced band that a pulsed model reaches
    fall from 0.00203 to 0.00033 of the signal in P25 and from 0.00204 to 0.00033 in D-STAR, the mean
    over 22 vectors of each vector's share, on the study's scratch build, `pzx3` against `k1`:
    [November 2026: the pulsed state](#november-2026-the-pulsed-state)):
    the pitch codes and V/UV rows that tell
    DVSI's decoder more than is modelled (above).
  - Inside voiced runs rate 33 loses the alignment with DVSI's waveform on 0.6% of its frames
    (D-STAR 0.2%, P25 0.3%; frames whose two frames before correlate 0.9 or more and whose two
    frames after fall under 0.8), more often across large pitch steps (over 20%: 6 of 65 in rate
    33). Of the frames at the largest steps, a frame-to-frame pitch ratio of 0.30 to 0.42 (three in
    rate 33, one in D-STAR), all four come back in line with the octave step ([The pulse train
    across an octave](#the-pulse-train-across-an-octave)). The rate-33 runs under 0.8 looked at
    lost the alignment at sequences of V/UV rows (row 14 then row 15, row 17 then row 13 or 15),
    which the mid-frame voicing change reads as DVSI's decoder plays them, the pulse train running
    on through a row-17 frame since the row-17 train change ([Rate 33: rows 13, 15 and
    17](#rate-33-rows-13-15-and-17)); the 0.6% was measured before those changes and before the
    mid-frame pitch wave. What moves the pulse train at the others is not known, and nor is what
    in DVSI's decoder sets rate 33's mid-frame pitch, which the library follows by a fit ([Rate 33:
    mid-frame pitch](#rate-33-mid-frame-pitch)).
- DVSI's rate-33 encoder sends V/UV indices 17-31 on some voice frames, where clause 4.2
  searches 0-16 and keeps the lowest index on ties, and this library's encoder never sends them.
  Rows 18-31 decode all unvoiced here, as row 16 does, but on their band grid ([Rate 33: V/UV
  vectors](#rate-33-vuv-vectors)), and play rate 33's pulsed state ([Rate 33: the pulsed
  rows](#rate-33-the-pulsed-rows)); row 17 decodes all unvoiced, as Annex B lists it, with its lowest
  column voiced in the frame's first half ([Rate 33: rows 13, 15 and
  17](#rate-33-rows-13-15-and-17)).
- Rate-33 silence frames play on 15 bands on the strength of the published AMBE-3000R measurements;
  DVSI's vectors do not tell 14 bands from 15. A run of consecutive silence frames is not verified,
  since DVSI's vectors hold none ([Rate 33 silence frames](#rate-33-silence-frames)).
- Rate-33 mute noise: DVSI's muted samples include +-2 (2.0% each on `dam_e10_hd`); ours are 0 and
  +-1 only ([Rate 33 mutes](#rate-33-mutes)).
- After a rate-33 tone, the pulse train of a voiced frame does not line up with DVSI's (median
  correlation of the next frame with DVSI's waveform 0.26, against 0.9 and more in running speech):
  DVSI's decoder continues the tone's pulse train into it ([Tones as harmonics of a
  fundamental](#tones-as-harmonics-of-a-fundamental)), and the library plays its tones beside the
  synthesizer. What DVSI's decoder plays after a tone inside speech (speech, a tone, then different
  speech) is not known: its tone vectors hold loud frames before a tone only between tones, where
  the tone signal is coded as voice frames ([Rate 33 frames after a mute or a
  tone](#rate-33-frames-after-a-mute-or-a-tone)).
- A rate-33 silence frame right after a mute or a tone rises from silence to its own model. No DVSI
  vector shows that case, and the rule is a design choice ([Rate 33 frames after a mute or a
  tone](#rate-33-frames-after-a-mute-or-a-tone)).
- Rate-33 tone frames with an unusable index (ID 255 included) and lost frames next to a tone
  follow the published AMBE-3000R measurements; no DVSI vector shows how DVSI's decoder plays them
  ([Rate 33 lost frames and tones](#rate-33-lost-frames-and-tones)).
- P25's 16-bit output holds more zero samples on DVSI's decoder than truncation gives (`dam80`:
  44.4% against 24.4%; rounding to nearest gives 12.3%). This library keeps truncating, which comes
  closer; DVSI's rule is not determined ([16-bit output](#16-bit-output)).

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

### November 2026: synthesis in subframes

Main (`8875729`) against the three synthesis changes of November 2026 together
([Unvoiced synthesis](#unvoiced-synthesis), [Synthesis in subframes](#synthesis-in-subframes) with
[Voiced onsets and offsets](#voiced-onsets-and-offsets), [Rate 33 tone
transitions](#rate-33-tone-transitions)), both scored at their joint lags. No partition was held
out ([Method](#method)). Decoding DVSI's bits:

| Codec | Partition | LSD vs DVSI (dB) | Waveform corr. vs DVSI, sample for sample | Envelope corr. vs DVSI | LSD vs input (dB) |
|---|---|---|---|---|---|
| D-STAR | development | 3.20 → 1.55 | 0.673 → 0.946 | 0.969 → 0.988 | 6.53 → 6.41 |
| D-STAR | validation | 2.72 → 1.37 | 0.655 → 0.930 | 0.963 → 0.985 | 4.95 → 4.85 |
| D-STAR | final | 3.87 → 1.81 | 0.499 → 0.932 | 0.946 → 0.978 | 5.70 → 5.56 |
| P25 | development | 3.27 → 1.76 | 0.831 → 0.918 | 0.973 → 0.988 | 6.44 → 6.32 |
| P25 | validation | 2.71 → 1.40 | 0.851 → 0.942 | 0.969 → 0.987 | 4.83 → 4.73 |
| P25 | final | 3.56 → 1.73 | 0.818 → 0.934 | 0.959 → 0.984 | 5.54 → 5.38 |
| Rate 33 | development | 2.95 → 1.47 | 0.896 → 0.943 | 0.979 → 0.991 | 6.40 → 6.29 |
| Rate 33 | validation | 2.51 → 1.26 | 0.923 → 0.952 | 0.973 → 0.988 | 4.94 → 4.84 |
| Rate 33 | final | 3.59 → 1.73 | 0.867 → 0.949 | 0.961 → 0.977 | 5.69 → 5.54 |

Every one of the 81 decodes comes closer to DVSI's in log-spectral distance, by 0.41 dB (P25
`dam80`) to 3.86 dB (D-STAR `dam10`). The level against DVSI's stays within 0.11 dB in every codec
and partition (means of the vectors' absolute offsets; up to 0.18 dB on main).

- **Channel errors.** The mute decisions are unchanged at every error rate, hard and soft. On
  the loud frames neither decoder mutes, the log-spectral distance from DVSI's decode falls by
  1.06 dB for D-STAR (hard decisions; 0.89 soft), 0.60 dB for P25 (0.49) and 0.92 dB for rate 33
  (0.81), means over the four error rates.
- **Tones.** On the 104 D-STAR and rate-33 tone vectors the tone detector agrees with DVSI's
  output more often on ten and less often on one (rate 33's `dtmfvbad`, 99.38% to 98.96% of its
  windows, see [Rate 33 tone transitions](#rate-33-tone-transitions)); no window plays a wrong
  tone, and the level error of the steady tones changes by at most 0.09 dB on a vector. P25,
  which carries tones as voice, is in [Known differences](#known-differences).
- **Encoders.** Every encoder's bits are unchanged, and so are the encoder comparisons: both
  runs' encodings were decoded by the same decoder.
- **Cost.** A complete frame decode takes 17-24% longer, about 1 µs a frame
  ([benchmarks](testing.md#subframe-synthesis)): the pitch pulses are rendered by scalar code.

The gate sets (`unvoiced-subframes-*`, `voiced-subframes-*`, `ambe2450-tone-fade-*` in
`tools/quality/dvsi_gates.json`) were declared per change before that change was scored, each
against the commit before it; their failing checks are named in the three sections.

### November 2026: rate-33 frame handling

Main (`75544c9`, 2.6.8) against the five changes of 2.6.9: [Rate 33 silence
frames](#rate-33-silence-frames), [Rate 33 mutes](#rate-33-mutes), [Rate 33 frames after a mute or
a tone](#rate-33-frames-after-a-mute-or-a-tone), [Rate 33 lost frames and
tones](#rate-33-lost-frames-and-tones) and [16-bit output](#16-bit-output). The first four change
mostly frames the scoreboard leaves out (tone transitions, and muted frames with their neighbours):
those frames are judged by the unit tests and by the measurements in their sections, and the
scoreboard's gate sets check that the rest does not move.

Two comparisons, each with both builds aligning freely first and then scored at their joint lags,
both runs' encodings decoded by the candidate's decoder. No partition was held out: the
measurements that found these rules read all the vectors, so every figure is retrospective.

**Frame handling** (main against the four frame-handling changes; `ambe2450-frames-*`). All five
gate sets pass.

- **Speech.** No decode changes in any mode on any partition (DVSI's speech vectors hold no silence
  frame, lost frame or tone), and no encoder's bits change.
- **Tones.** D-STAR and P25 are unchanged in every tone metric. On rate 33 the tone detector agrees
  with DVSI's output on 99.45% of `dtmfvbad`'s windows, from 98.96% ([Rate 33 tone
  transitions](#rate-33-tone-transitions)), and as before on the other 51 vectors; no window plays a
  wrong tone, and the level, frequency and balance errors of the steady tones are unchanged.
- **Channel errors** (`dam` at 1.1, 2.2, 4.9 and 10.6% bit errors, hard and soft decisions). The
  mute decisions are unchanged in every mode at every error rate (rate 33's agreement with DVSI's
  on loud frames: 100%, 100%, 100% and 92.5% hard, 96.9% soft at 10.6%), and D-STAR's and P25's
  output is unchanged. On the loud frames neither decoder mutes, rate 33's log-spectral distance
  from DVSI's decode falls by 0.13 dB at 10.6% with hard decisions (0.03 dB as the mean of the four
  rates; coverage 0.75) and by 0.015 dB at 10.6% with soft ones, and does not change at the lower
  rates. That metric leaves out muted frames and their neighbours, which are the frames these
  changes are about; the sections above measure them.

**Rounding** (the commit before it against the rounding; `pcm-rounding-*`).
`pcm-rounding-validation`, `-final` and `-channel` pass; `pcm-rounding-development` and `-tones`
fail:

- **Development, the level check** (the level against DVSI's decode may move by 0.05 dB on any
  vector). On `dam80` it moves by +0.18 dB for D-STAR (-0.02 to +0.16 dB) and +0.15 dB for rate 33
  (-0.14 to +0.01 dB), and on `dam40` by +0.09 dB (D-STAR, -0.04 to +0.05) and +0.08 dB (rate 33,
  -0.05 to +0.03). Those two inputs lie 80 and 40 dB down, where truncation's dead zone removed part
  of the signal. With rounding, `dam80` sits where the same decoder puts `dam` against DVSI's level
  (D-STAR +0.16 against +0.10 dB, rate 33 +0.01 against +0.00 dB); with truncation it sat 0.12 and
  0.15 dB below that. Every other vector moves by less than 0.01 dB.
- **Tones**, two checks declared with no tolerance ("not higher on any vector"), by amounts of
  rounding noise: the relative frequency error rises by at most 3.9e-7 (D-STAR) and 6.8e-7 (rate
  33) on a vector, and the balance of a dual tone's components by at most 0.0008 and 0.0003 dB. The
  detector's agreement and the wrong-tone fraction are unchanged on every vector, and the level
  error of the steady tones changes by at most 0.0012 dB.

The other checks pass. The log-spectral distance from DVSI's decode falls by 1.04 dB on D-STAR
`dam80` (3.42 to 2.38), 0.70 dB on D-STAR `dam40` (2.49 to 1.79), 0.40 dB on rate-33 `dam80` (3.03
to 2.63) and 0.19 dB on rate-33 `dam40` (1.90 to 1.71); on the other vectors it moves by less than
0.01 dB (it rises on 21 D-STAR and 19 rate-33 vectors, by at most 0.0017 dB). The sample-for-sample
waveform correlation with DVSI's decode is unchanged to the fourth decimal on every vector but
`dam80` (+0.0007 D-STAR, +0.0001 rate 33). P25's decodes are byte-identical, every encoder's bits
are unchanged, and under channel errors the mute decisions are unchanged and the distance on the
loud frames neither decoder mutes moves by 0.0002 dB at most. Partition means:

| Codec | Partition | LSD vs DVSI (dB) | Waveform corr. vs DVSI, sample for sample | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| D-STAR | development | 1.548 → 1.461 | 0.946 → 0.946 | -0.005 → +0.010 | 6.413 → 6.431 |
| D-STAR | validation | 1.371 → 1.372 | 0.930 → 0.930 | -0.020 → -0.018 | 4.853 → 4.854 |
| D-STAR | final | 1.805 → 1.806 | 0.932 → 0.932 | +0.040 → +0.042 | 5.560 → 5.561 |
| Rate 33 | development | 1.470 → 1.440 | 0.943 → 0.943 | -0.013 → +0.000 | 6.285 → 6.264 |
| Rate 33 | validation | 1.257 → 1.257 | 0.952 → 0.952 | -0.055 → -0.053 | 4.837 → 4.838 |
| Rate 33 | final | 1.728 → 1.728 | 0.949 → 0.949 | -0.089 → -0.088 | 5.539 → 5.539 |

P25 is unchanged in every column. The development rows move because `dam40` and `dam80` are
development vectors. D-STAR's distance from the input rises on `dam80`, 6.37 to 6.64 dB, toward
DVSI's own decode's 6.71 dB; rate 33's falls, 6.23 to 5.98 dB (DVSI's 6.03 dB).

**Byte identity.** Decodes of DVSI's vectors (27 speech, 52 tone and 8 channel vectors, and for
D-STAR and rate 33 the three DTX vectors), PCM and per-frame flags, hard and soft:

- Main against the four frame-handling changes (before the rounding): every D-STAR stream (90 of
  90) and every P25 stream (87 of 87) is identical in PCM and flags. Rate 33's flags are identical
  on all 90 streams and its PCM on 78; the other twelve are the seven tone vectors with a voice
  frame right after a tone ([Rate 33 frames after a mute or a
  tone](#rate-33-frames-after-a-mute-or-a-tone)), `dam_e10_hd`, `dam_e10_sd` and the three DTX
  vectors. No rate-33 speech vector and no channel vector up to 4.9% bit errors changes.
- The rounding against the commit before it: P25 is identical on all 87 streams, PCM and flags;
  every D-STAR and rate-33 stream differs in PCM, and none in flags.

**DTX and channel-error vectors**, our decode against DVSI's decode of the same bits, main and now
(details in each section):

| Frames | Measure | Main | Now |
|---|---|---|---|
| Silence frames of the DTX vectors, after a voiced frame (53) | frame correlation, median | 0.631 | 0.979 |
| Silence frames of the DTX vectors, after a frame with no voiced band (84) | frame correlation, median | 0.852 | 0.994 |
| `dam_e10_hd`, muted frames inside muted runs (185) | RMS on the 16-bit scale, median (DVSI 0.59) | 2.46 | 0.58 |
| `dam_e10_hd`, mute onsets (13) | waveform correlation over the first 60 samples, median | -0.08 | 0.94 |
| `dam_e10_hd`, frames leaving a mute (14) | correlation, first half / second half / next frame, medians | 0.00 / 0.04 / 0.02 | 0.95 / 0.95 / 0.95 |
| `dam_e10_hd`, frames leaving a mute (14) | first half's level against DVSI's, median | -54.2 dB | +0.95 dB |
| Tone vectors, voice frames right after a tone (1206) | mean absolute level error, whole frame | 0.58 dB | 0.37 dB |

The share of samples identical to DVSI's, with the rounding, rises on every codec it applies to
(mean over the 27 speech vectors 0.140 to 0.177 for D-STAR and 0.174 to 0.190 for rate 33,
[16-bit output](#16-bit-output)). The golden PCM hashes that moved are the int16 hashes of the
D-STAR and rate-33 sequences, which now hash what their `short` entry points return; the float
hashes are unchanged.

### November 2026: encoders and rate-33 rows 18 to 31

Main (`68b110e`, 2.6.11) against the three changes of 2.6.12:

- [Encoder pitch code of unvoiced frames](#encoder-pitch-code-of-unvoiced-frames): every encoder
  sends a voice frame with no voiced band at one fixed pitch code, its amplitudes measured on that
  code's harmonic grid.
- [Joint choice of the AMBE spectral vectors](#joint-choice-of-the-ambe-spectral-vectors): the
  D-STAR and rate-33 encoders refine their six spectral vectors together, weighted toward the strong
  harmonics, with their magnitude scales recalibrated.
- [Rate 33: V/UV vectors](#rate-33-vuv-vectors): the rate-33 decoder plays V/UV rows 18 to 31 on the
  published band grid, `w0 = 2.81 / L`.

Both builds aligned freely first and were then scored at their joint lags, both runs' encodings
decoded by the candidate's decoder (our encoders never send rows 18 to 31, so the grid changes no
decode of their bits). The development and validation partitions were measured while the changes
were made; the final partition was held out and scored once, after all three were fixed
([Method](#method)).

**Encoders.** Our decoded encoding of the input against DVSI's encoding of it and against the input,
partition means (LSD is log-spectral distance):

| Codec | Partition | LSD vs DVSI's encoding (dB) | LSD vs input (dB) | Level gap to DVSI's encoding (dB) | Envelope corr. vs input | Envelope corr. vs DVSI's encoding |
|---|---|---|---|---|---|---|
| D-STAR | development | 5.606 → 5.504 | 5.492 → 5.434 | 0.435 → 0.184 | 0.921 → 0.923 | 0.925 → 0.921 |
| D-STAR | validation | 4.377 → 4.295 | 4.680 → 4.619 | 0.398 → 0.141 | 0.932 → 0.933 | 0.938 → 0.938 |
| D-STAR | final | 5.543 → 5.471 | 5.462 → 5.398 | 0.649 → 0.100 | 0.898 → 0.899 | 0.911 → 0.910 |
| P25 | development | 5.494 → 5.443 | 5.313 → 5.314 | 0.229 → 0.241 | 0.932 → 0.931 | 0.939 → 0.935 |
| P25 | validation | 4.227 → 4.180 | 4.535 → 4.493 | 0.151 → 0.160 | 0.935 → 0.934 | 0.949 → 0.947 |
| P25 | final | 5.237 → 5.211 | 5.298 → 5.291 | 0.399 → 0.363 | 0.904 → 0.905 | 0.918 → 0.915 |
| Rate 33 | development | 5.669 → 5.545 | 5.305 → 5.276 | 0.213 → 0.188 | 0.931 → 0.935 | 0.924 → 0.925 |
| Rate 33 | validation | 4.324 → 4.224 | 4.710 → 4.660 | 0.132 → 0.123 | 0.931 → 0.933 | 0.936 → 0.936 |
| Rate 33 | final | 5.466 → 5.388 | 5.468 → 5.452 | 0.168 → 0.143 | 0.900 → 0.905 | 0.909 → 0.912 |

Every encoding of every speech vector comes closer to DVSI's encoding in log-spectral distance, and
every codec's comes closer to the input on every partition but P25's development one (+0.0015 dB,
`dam10` +0.30 dB, from the pitch-code change). D-STAR's encoding, a median 0.45 dB (development) and
0.43 dB (validation) below DVSI's on main, now sits at -0.01, +0.16 and -0.10 dB on the three
partitions. The envelope correlation with DVSI's encoding falls by up to 0.004 for D-STAR and P25,
from the pitch-code change. The share of 1 kHz bands voiced as DVSI's encoder voices them is
unchanged on every vector, and no encoding holds a tone or silence frame.

**Rate 33 decoding DVSI's bits.** The D-STAR and P25 decodes are unchanged on every speech vector.
Rate 33, partition means and the mean paired change in log-spectral distance from DVSI's decode with
its 95% interval:

| Partition | LSD vs DVSI (dB) | Waveform corr. vs DVSI, sample for sample | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) | Change in LSD vs DVSI (dB) |
|---|---|---|---|---|---|---|
| development | 1.440 → 1.432 | 0.943 → 0.943 | 0.992 → 0.992 | +0.000 → +0.000 | 6.264 → 6.263 | -0.0081 [-0.0130, -0.0053] |
| validation | 1.257 → 1.249 | 0.952 → 0.952 | 0.988 → 0.988 | -0.052 → -0.053 | 4.838 → 4.836 | -0.0079 [-0.0101, -0.0043] |
| final | 1.728 → 1.710 | 0.949 → 0.949 | 0.977 → 0.977 | -0.088 → -0.089 | 5.539 → 5.536 | -0.0183 [-0.0211, -0.0156] |

The distance from DVSI's decode falls on all 27 vectors; on development and validation these are
the band-grid change's own figures, since the other two changes touch no decoder.

**PESQ and ESTOI** ([Method](#method)) of the three encoders, main against the three changes, both
encodings decoded by this library's 2.6.11 decoder (mean paired change ± standard error, vectors
improved; main's encoding scores PESQ 2.98, 3.03 and 3.13 on the 25 vectors and 2.94, 2.94 and
3.12 on the 70 recordings, D-STAR, rate 33 and P25):

| Codec | PESQ, 25 vectors | ESTOI, 25 vectors | PESQ, 70 recordings | ESTOI, 70 recordings |
|---|---|---|---|---|
| D-STAR | +0.114 ± 0.009 (25/25) | +0.0135 ± 0.0018 (25/25) | +0.131 ± 0.010 (66/70) | +0.0178 ± 0.0011 (68/70) |
| Rate 33 | +0.086 ± 0.011 (24/25) | +0.0172 ± 0.0026 (25/25) | +0.095 ± 0.009 (64/70) | +0.0216 ± 0.0010 (70/70) |
| P25 | +0.058 ± 0.010 (23/25) | +0.0039 ± 0.0020 (17/25) | +0.050 ± 0.009 (53/70) | +0.0037 ± 0.0011 (40/70) |

Each is the sum of the encoder changes' own figures. Decoded by this branch's own decoder, the
rate-33 figures are the same (+0.086 ± 0.011 and +0.0172 ± 0.0026 on the 25 vectors, +0.095 ± 0.009
and +0.0216 ± 0.0010 on the 70 recordings). On our encodings this branch's decoder gives exactly the
samples of main's decoder built the same way, in every codec (our encoders never send rate-33 rows
18 to 31); the 2.6.11 build used for the scores, built with fast-math and link-time optimization,
differs from such a plain Release build by one 16-bit step on 0.07 to 0.08% of the samples. The
rate-33 decode of DVSI's bits gains 0.0035 ± 0.0016 PESQ against the input ([Rate 33: V/UV
vectors](#rate-33-vuv-vectors)).

**The final partition.** The gate set `encoder-quality-final` was declared before the final scoring
and built from the three changes' own gate sets: their metrics and rules, their limits added where
two changes bound the same metric (the D-STAR and rate-33 encodings' distance from DVSI's encoding
and from the input may rise by 0.02 dB and their level gap by 0.2 dB, P25's by 0.01 and 0.1 dB, the
envelope correlation with the input fall by 0.002), "unchanged" kept only for the D-STAR and P25
decodes of DVSI's bits, and the band-grid change's checks for the rate-33 decode. It passes on every
check: the D-STAR and P25 decodes are unchanged; the encoders' distance from DVSI's encoding falls by
0.071 (D-STAR), 0.077 (rate 33) and 0.026 dB (P25), from the input by 0.064, 0.016 and 0.007 dB, and
their level gap by 0.549, 0.024 and 0.036 dB, while the D-STAR and rate-33 envelope correlations with
the input rise by 0.001 and 0.004; the rate-33 decode comes 0.018 dB closer to DVSI's (95% interval
-0.021 to -0.016), its waveform and envelope correlations with DVSI's rise by 0.00006 and 0.00002,
its level moves by -0.0007 dB and its distance from the input falls by 0.003 dB.

**The gate that failed.** `joint-search-development` failed on one check, the log-spectral distance
of the D-STAR and rate-33 encodings from the input on the development vectors (+0.0105 and +0.0123
dB against a limit of 0.01), because the search's weighting trades error on weak harmonics for error
on strong ones and that distance counts every spectrum bin's dB alike ([Joint choice of the AMBE
spectral vectors](#joint-choice-of-the-ambe-spectral-vectors)). Over the three changes that distance
falls for both codecs on every partition (development: -0.058 and -0.030 dB).

### November 2026: phase regeneration and pitch

Main (`7fb5ca6`, 2.6.12) against the three decoder changes of 2.6.13:

- [Phase regeneration](#phase-regeneration): the codec paths regenerate voiced phases from each
  frame's decoded log2 amplitudes, as DVSI's decoders do.
- [P25 pitch](#p25-pitch): P25 plays each pitch code at a fit to the pitch DVSI's decoder plays on
  speech, eq 46's moved by at most 0.071%. DVSI's decoder plays a tone sent as voice frames at
  another pitch, and there the change lowers the steady frames' median correlation with DVSI's
  waveform from 0.880 to 0.656 (0.847 on main). The tone model change plays such frames at the pitch
  DVSI's decoder
  plays them at, eq 46's fundamental cut to a 2^-19 grid, and that cost is gone ([November 2026:
  tones and the octave step](#november-2026-tones-and-the-octave-step)).
- [Rate 33: mid-frame pitch](#rate-33-mid-frame-pitch): between two voiced rate-33 frames of
  different pitch the mid-frame pitch is the geometric mean times `1 - 3.0e-4`.

Both builds aligned freely first and were then scored at their joint lags, both runs' encodings
decoded by the candidate's decoder. The development and validation partitions were measured while
the changes were made; the final partition was held out and scored once, after the three changes
and the review's fixes were complete ([Method](#method); its result is below).

**Decoding DVSI's bits**, partition means:

| Codec | Partition | Waveform corr. vs DVSI, sample for sample | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|---|
| D-STAR | development | 0.9461 → 0.9844 | 1.461 → 1.357 | 0.9877 → 0.9881 | +0.010 → +0.025 | 6.428 → 6.426 |
| D-STAR | validation | 0.9299 → 0.9612 | 1.372 → 1.248 | 0.9846 → 0.9849 | -0.018 → -0.006 | 4.854 → 4.847 |
| P25 | development | 0.9182 → 0.9672 | 1.762 → 1.649 | 0.9882 → 0.9883 | +0.032 → +0.046 | 6.320 → 6.320 |
| P25 | validation | 0.9422 → 0.9855 | 1.402 → 1.270 | 0.9871 → 0.9873 | -0.025 → -0.014 | 4.734 → 4.728 |
| Rate 33 | development | 0.9428 → 0.9720 | 1.432 → 1.286 | 0.9916 → 0.9919 | +0.000 → +0.014 | 6.263 → 6.272 |
| Rate 33 | validation | 0.9524 → 0.9788 | 1.249 → 1.106 | 0.9880 → 0.9883 | -0.053 → -0.040 | 4.836 → 4.834 |

Every one of the 75 decodes correlates more closely with DVSI's waveform and comes closer to DVSI's
decode in log-spectral distance, and so do the six decodes of the final partition (81 in all;
below). On development and validation the partition means of the distance from the input move by
at most 0.009 dB (rate-33 development, +0.0084), and those of the level against DVSI's rise by 0.011
to 0.015 dB. On the same vectors the residual energy of our decode against DVSI's, sample for
sample at lag 0, falls from 0.1085 to 0.0317 (D-STAR), 0.1642 to 0.0664 (P25) and 0.1142 to 0.0561
(rate 33) of DVSI's on development, and from 0.1403 to 0.0779, 0.1155 to 0.0289 and 0.0949 to
0.0424 on validation. These are means over the vectors of each vector's own ratio; pooled over the
development vectors, each scaled to the same RMS level of DVSI's decode, the residual now is
0.037, 0.069 and 0.065 ([Known differences](#known-differences)). Where a figure here differs
from the same build's in the sections' tables (P25's development distance and envelope correlation
on main, 1.759 dB and 0.9886 there; rate 33's development distance from the input now, 6.271 dB
there), the two comparisons' joint lags differ on some vectors.

- **Tones and channel errors** (scored per change, in the sections). No change alters a mute
  decision at any error rate, hard or soft, or makes a D-STAR or rate-33 tone window play a wrong tone, and at every step
  the distance from DVSI's decode on the loud frames neither decoder mutes falls in every mode the
  step changes. The one tone check that fails, D-STAR `knox_12`'s level error in the phase change
  (+0.070 dB), is on a KNOX vector, which D-STAR codes as voice. On P25's tone vectors, whose tone
  rows no gate set of the P25 pitch checks, the P25 pitch lowers the steady tone frames' median
  correlation with DVSI's waveform from 0.880 to 0.656, since DVSI's decoder plays a tone sent as
  voice frames at another pitch ([P25 pitch](#p25-pitch); since the tone model change the library
  plays it so:
  [P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)).
- **Encoders.** Every encoder's bits are unchanged: the three encoders give the same parameter rows
  in all four builds (main and each change) for the 95 DVSI inputs other than `irstia` and `mark`.
- **Cost.** A complete hard-decision frame decode takes 3.2-3.4% longer, measured at `3e0c6ce`
  against main ([benchmarks](testing.md#phase-regeneration-and-pitch)).

**Gate sets**, each declared before its change was scored, and evaluated against the commit before
it (main, the phase change, the P25 pitch, the mid-frame pitch, in turn) in one run over the
development, validation, tone and channel partitions ([Scoring history](#scoring-history)):

| Gate set | Result | Failing check |
|---|---|---|
| `phase-source-development` | fails | P25 waveform correlation on every vector (`t01` -0.0087, allowed -0.002) |
| `phase-source-validation` | passes | |
| `phase-source-tones` | fails | D-STAR steady-tone level error on every vector (`knox_12` +0.070 dB, allowed 0.05) |
| `phase-source-channel` | passes | |
| `p25-pitch-development` | passes | |
| `p25-pitch-validation` | passes | |
| `p25-pitch-tones` | passes | |
| `p25-pitch-channel` | passes | |
| `r33-mid-pitch-development` | fails | rate-33 waveform correlation on every vector (`tambf22a` -0.0029, allowed -0.002) |
| `r33-mid-pitch-validation` | fails | rate-33 waveform correlation rising by 0.003 with its 95% interval above 0 (+0.0066, interval -0.0001 to +0.0119) |
| `r33-mid-pitch-tones` | passes | |
| `r33-mid-pitch-channel` | passes | |

Every other check of the failing sets passes, and the failures are kept as declared. The gate set of
the final partition, `decoder-phase-pitch-final`, scored for the three changes together against
main, passes (below). Over the three changes the two vectors that fall in the waveform checks end
closer to DVSI's waveform than on main (`t01` 0.9476 to 0.9956, `tambf22a` 0.9455 to 0.9712);
`knox_12`'s steady-tone level error stays further from DVSI's (0.117 to 0.187 dB with the phase
change, which the later changes, leaving D-STAR's decodes as they are, do not move).

**PESQ and ESTOI** ([Method](#method)) of the decode of DVSI's bits, main against the three
changes (mean paired change ± standard error, vectors improved, of 25; DVSI's own decode of the same
bits scores PESQ 2.933, 3.026 and 2.980 against the input for D-STAR, P25 and rate 33):

| Codec | PESQ vs input | ESTOI vs input | PESQ vs DVSI's decode |
|---|---|---|---|
| D-STAR | 2.9476 → 2.9468: -0.0008 ± 0.0020 (11/25) | +0.0004 ± 0.0001 (22/25) | +0.0034 ± 0.0004 (19/25) |
| P25 | 3.0399 → 3.0333: -0.0066 ± 0.0046 (13/25) | +0.0001 ± 0.0002 (17/25) | +0.0042 ± 0.0013 (21/25) |
| Rate 33 | 2.9615 → 2.9530: -0.0085 ± 0.0032 (8/25) | -0.0001 ± 0.0005 (18/25) | +0.0048 ± 0.0008 (23/25) |

The three changes bring the waveform to DVSI's and leave the input measures nearly where they were.
D-STAR's ESTOI against the input rises on 22 of 25 vectors, by +0.0004 ± 0.0001 (four standard
errors); P25's and rate 33's ESTOI do not move by a standard error. PESQ against the input falls by
less than a standard error for D-STAR, by about one and a half for P25 (-0.0066 ± 0.0046) and by
about two and a half for rate 33 (-0.0085 ± 0.0032). The project's rule is the specification first
and DVSI's decoder where the specification is silent. The phase regeneration belongs to the voiced
synthesis already recorded as a departure, the rate-33 mid-frame pitch lies where TIA-102.BABA-1 is
silent (it defines no mid-frame model), and the P25 pitch is a recorded departure from TIA-102.BABA
eq 46, by up to 0.071%, taken for the measured waveform agreement with DVSI's decoder on speech
([Departures from the specification](#departures-from-the-specification)). The changes are judged by
the waveform against DVSI's decoder, not by a gain against the input.

**Byte identity**, by the scoreboard's hashes of our decode of DVSI's bits (25 speech vectors, 52
tone vectors and 8 channel-error vectors per mode, hard and soft decisions among the last):

- The phase change changes every speech and channel-error decode in every mode, and 25 D-STAR, all
  52 P25 and 10 rate-33 tone decodes.
- The P25 pitch leaves every D-STAR and rate-33 decode unchanged, speech, tone and channel-error
  vectors alike, and changes every P25 one.
- The rate-33 mid-frame pitch leaves every D-STAR and P25 decode unchanged, speech, tone and
  channel-error vectors alike. It changes every rate-33 speech and channel-error
  decode and 5 rate-33 tone decodes.

**The final partition.** The gate set `decoder-phase-pitch-final` was declared before the final
scoring and built from the development gate sets of the three changes: their metrics and rules,
their limits added where two changes bound one metric of a mode (P25: the phase change and the P25
pitch; rate 33: the phase change and the mid-frame pitch; D-STAR: the phase change alone). So the
waveform correlation with DVSI's decode must rise by at least 0.015 (D-STAR), 0.025 (P25) and 0.018
(rate 33), with its 95% interval above 0, and fall by no more than 0.002, 0.004 and 0.004 on any
vector; the distance from DVSI's decode may rise by 0.01, 0.02 and 0.02 dB (0.05, 0.1 and 0.1 on any
vector), the envelope correlation fall by 0.002, 0.004 and 0.004, the distance from the input rise
by 0.02, 0.03 and 0.03 dB and the level move by 0.1, 0.15 and 0.15 dB, and every encoder's bits stay
unchanged. The final partition was scored once, for the three changes together, against it: after
the review's fixes, main (`7fb5ca6`) against this branch at `6491fe2` (the commits after it change
tests, comments and documentation only), on `irstia` and `mark`, as above (both builds aligned
freely first and then scored at their joint lags, both runs' encodings decoded by the candidate's
decoder). It passes all 24 checks (16 declared, the per-mode ones counted per mode). Partition means
(two vectors):

| Codec | Waveform corr. vs DVSI, sample for sample | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| D-STAR | 0.9316 → 0.9656 | 1.922 → 1.781 | 0.9753 → 0.9766 | +0.042 → +0.058 | 5.561 → 5.544 |
| P25 | 0.9342 → 0.9798 | 1.731 → 1.570 | 0.9835 → 0.9840 | -0.052 → -0.034 | 5.379 → 5.355 |
| Rate 33 | 0.9489 → 0.9766 | 1.710 → 1.559 | 0.9771 → 0.9781 | -0.089 → -0.072 | 5.536 → 5.521 |

The waveform correlation rises on all six decodes: D-STAR `irstia` 0.8980 → 0.9426 and `mark`
0.9652 → 0.9886, P25 0.9134 → 0.9727 and 0.9550 → 0.9869, rate 33 0.9347 → 0.9697 and 0.9632 →
0.9834. The distance from DVSI's decode falls on all six, by 0.10 to 0.22 dB, and the distance from
the input does not rise on five (rate 33 `irstia`: 5.562 → 5.563 dB). Against the gate's limits: the
correlation rises by 0.0340 (D-STAR, 0.015 required), 0.0456 (P25, 0.025) and 0.0276 (rate 33,
0.018), on every vector; the distance from DVSI's decode falls by 0.141, 0.162 and 0.150 dB; the
envelope correlation rises by 0.0013, 0.0005 and 0.0010; the distance from the input falls by 0.017,
0.025 and 0.015 dB; the level moves by +0.016, +0.018 and +0.017 dB; and every encoder's bits are
unchanged on both vectors.

### November 2026: tones and the octave step

The head of PR #196 (`fac0111`) against the three decoder changes of the tone model change set
(`c736edf` to `c2ba220`):

- [The pulse train across an octave](#the-pulse-train-across-an-octave): the pulse train is taken an
  octave down or up before a subframe step whose period ratio lies within 20% of 2 or of 1/2.
- [Tones as harmonics of a fundamental](#tones-as-harmonics-of-a-fundamental): rate-33 and D-STAR
  tone frames play the harmonics DVSI's decoders play, and a tone after a frame that was not a tone
  frame starts as a fresh pulse train.
- [P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames): a P25 voice frame that carries
  a DTMF, KNOX or single tone as DVSI's encoder sends it plays the tone model DVSI's decoder plays.

Both builds aligned freely first and were then scored at their joint lags, both runs' encodings
decoded by the candidate's decoder. The study that found these rules read DVSI's tone vectors and
the development and validation vectors, never `irstia` or `mark`; the change set was measured and
gated on the development, validation, tone and channel partitions, and the final partition was held
out and scored once, after the review ([Method](#method); below).

**Tones**, against DVSI's decode of the same bits, sample for sample at lag 0, measured with a local
script on DVSI's 52 tone vectors of each codec. The frames it counts:

- A **tone frame**: in rate 33 and D-STAR a frame the library's decoder flags as a tone
  (`MBE_PROCESS_FLAG_TONE`), keyed by its tone ID (D-STAR's index mapped to the AMBE+2 ID of the
  same pair); in P25 a frame that passes the script's own implementation of the library's tone test
  (the table's pairs of b0 and b1 without the call-progress tones, the tone harmonics within L, the
  rms of the other decoded magnitudes under 0.375 of the tone's for a dual tone and under 0.0235 for
  a single tone), keyed by its pitch code and tone harmonics. That implementation is not a copy of
  the library's classifier (it lacks the positive-normal test of the tone magnitudes and adds 1e-12
  to a denominator); its frame selection agrees with the library's test on every one of the 26,200
  frames of the P25 tone vectors.
- A **steady frame**: a tone frame between two tone frames of the same key, whose DVSI output has a
  standard deviation above 100 (16-bit scale). A vector counts if it holds three steady frames or
  more (rate 33 50, D-STAR 34, P25 47), and the other frames below are counted on those vectors
  only.
- A **first tone frame**: a tone frame after a frame that is not a tone frame (a change from one
  tone to another is not one), whose DVSI output is above the same level.
- The **frame after a tone**: a frame that is not a tone frame after a tone frame, over its first 80
  samples, whose DVSI output over those 80 samples is above the same level.

A frame's correlation is over its 160 samples, not centred; the residual is the energy of the
difference over that of DVSI's output, in dB, over the frame (the frame after a tone: its first 80
samples).

| Codec | Steady frames | Correlation, median | Correlation, mean | First tone frames: median residual |
|---|---|---|---|---|
| Rate 33 | 12,555 | +0.0014 → +0.9997 | -0.015 → +0.647 | +3.0 → -31.6 dB (1,021 frames) |
| D-STAR | 10,570 | +0.0005 → +0.9995 | -0.017 → +0.584 | +3.0 → -31.3 dB (989) |
| P25 | 12,721 | +0.3876 → +0.9862 | +0.2871 → +0.9628 | -2.5 → -13.0 dB (1,254) |

- In rate 33, 38 of the 50 vectors reach a median steady-frame correlation of 0.99 or more. The 12
  that do not are the cases not modelled ([Known differences](#known-differences)): tones after a
  voiced frame (`cp2`, `cp31`), tones changed back to back (`dtmf4`, `dtmf4r`, `dtmf8`, `dtmf8r`,
  `dtmf15n`, `dtmf15p`, `dtmf25`, `dtmf50`, `dtmf100`) and `alert`'s single tones back to back.
- In P25 the mean over the 47 vectors of each vector's median rises from +0.1138 to +0.9802. These
  P25 figures are on this measurement's selection, which takes the frames P25's tone test passes;
  PR #196 quoted 0.880 to 0.656 for its P25 pitch on another selection and calculation (the 17,334
  steady loud frames of PR #196's measurement, [P25 pitch](#p25-pitch)). The like-for-like
  comparison is the row "tone model" under **Before this change** in [P25 tones sent as voice
  frames](#p25-tones-sent-as-voice-frames): on those frames the median correlation is 0.986 (from
  0.656 with the P25 pitch change, `470374c`), the mean 0.940 (from 0.408), the steady frames' level
  error 0.279 dB (from
  0.560), and the mean over the vectors of the whole-file correlation 0.918 (from 0.127).
- Over the first 80 samples after a tone the median residual falls from +2.9 to -30.5 dB in rate 33
  (976 frames), from +3.0 to -23.4 dB in D-STAR (960) and from -6.2 to -10.3 dB in P25 (708).
- On the scoreboard (means over the vectors that pair), the steady tones' relative frequency error
  against DVSI's output falls by 0.001579 (D-STAR) and 0.003214 (rate 33), to about zero from means
  of 0.0016 and 0.0033, and rises on no vector; their level error falls by 0.027 and 0.012 dB, and
  the balance error of a dual tone's components rises by at most 0.002 dB on a vector and falls by
  up to 0.051 dB (`cp0`: 0.0508 dB in D-STAR, 0.0458 in rate 33). No vector plays a wrong tone more
  often, and the encoders' tone frames agree with DVSI's exactly as before. In P25, on the 47 or 48
  vectors that pair, the balance error falls from 2.0999 to 0.0159 dB, the level error from 0.2587
  to 0.0959 dB and the wrong-tone fraction from 0.0169 to 0; the detector's agreement with DVSI's
  output (not gated), over all 52 vectors, rises from 0.9477 to 0.9994.

**Speech**, decoding DVSI's bits: the waveform correlation with DVSI's decode, sample for sample,
partition means (vectors up and down at four decimals), and the residual energy against DVSI's
decode, means over the vectors of each vector's own ratio:

| Codec | Partition | Waveform corr. vs DVSI | Vectors up / down | Residual energy vs DVSI |
|---|---|---|---|---|
| D-STAR | development | 0.9844 → 0.9907 | 10 / 0 | 0.0317 → 0.0191 |
| D-STAR | validation | 0.9612 → 0.9820 | 3 / 0 | 0.0779 → 0.0363 |
| P25 | development | 0.9672 → 0.9687 | 4 / 1 (`ucarm15a` -0.0003) | 0.0664 → 0.0633 |
| P25 | validation | 0.9855 → 0.9855 | 1 / 1 (`tambf22b` -0.0001) | 0.0289 → 0.0289 |
| Rate 33 | development | 0.9720 → 0.9737 | 4 / 0 | 0.0561 → 0.0527 |
| Rate 33 | validation | 0.9788 → 0.9788 | 2 / 0 | 0.0424 → 0.0423 |

The distance from DVSI's decode falls on both partitions (development by 0.0041 dB for D-STAR,
0.0008 for P25 and 0.0015 for rate 33, rising on no vector; validation by 0.0104, 0.0008 and 0.0005
dB, rising by at most 0.0006 dB on a vector). The envelope correlation with DVSI's decode moves by
less than 0.0001, the distance from the input and the level against DVSI's by less than 0.0004 and
0.0005 dB (partition means), and no encoder's bits change.

**PESQ** ([Method](#method)) of the decode of DVSI's bits on the 25 development and validation
vectors (mean paired change ± standard error clustered by recording): against the input +0.0008 ±
0.0007 (D-STAR), -0.0000 ± 0.0002 (P25) and +0.0003 ± 0.0002 (rate 33); against DVSI's decode
+0.0104 ± 0.0045 (D-STAR, 16 of 25 vectors better), +0.0003 ± 0.0008 (P25) and +0.0009 ± 0.0012
(rate 33).

**Channel errors.** No mute decision changes at any error rate, hard or soft. On the loud frames
neither decoder mutes the distance from DVSI's decode falls by 0.030 dB for D-STAR (hard decisions;
0.008 soft), 0.010 dB for P25 (it rises by 0.0004 with soft decisions) and 0.023 dB for rate 33
(0.007), means over the four error rates.

**Cost.** A complete hard-decision frame decode, medians of seven invocations
([benchmarks](testing.md#tones-and-the-octave-step)): with the index every codec's median is
within 0.6% of `fac0111`'s (-0.55% D-STAR, +0.3% AMBE+2, +0.1% IMBE) and the ranges of the
invocations overlap. Before the P25 tone test looked a frame's pitch code up in an index, IMBE's
median was 1.4% higher, with ranges that also overlapped (6.696-6.801 against 6.776-6.914 us), so
this measurement does not isolate the cost of the scan.

**Gate sets**, declared in `c736edf` before any build of the change set was scored, and evaluated
in one run over the development, validation, tone and channel partitions, `fac0111` against the
change set:

| Gate set | Result | Failing checks |
|---|---|---|
| `tone-model-development` | passes | |
| `tone-model-validation` | passes | |
| `tone-model-tones` | fails 5 of its 20 checks | P25: the balance error falling by 1.0 dB, the level error not rising, the frequency error not rising by 0.0005 on any vector and the wrong-tone fraction not rising, all four for `knox_9`, which the comparison cannot pair; rate 33: the detector's agreement not falling on any vector (`alltone` -0.003175, allowed 0) |
| `tone-model-channel` | passes | |
| `tone-model-final` | passes (scored once, after the review, at `c2ba220`; below) | |

The four P25 checks fail because the scoreboard cannot pair `knox_9`: `fac0111`'s P25 decode of it
has no steady tone window the scoreboard's tone detector recognises, so its tone metrics are null
there and a paired check refuses the vector. On the vectors that pair, each of the four moves the
way its check asks (the figures under **Tones** above; the frequency error rises by at most 8.7e-05
on a vector, `cp31`), and `knox_9`'s detector agreement goes from 0 to 1. The rate-33 failure is
one of the 315 detector windows of `alltone` (whose tones are 156 runs, 120 of them single tones);
the mean agreement over the vectors rises by 0.003886. Both failures are kept as declared, and the
gate sets were not changed. The description text of `tone-model-tones` in
`tools/quality/dvsi_gates.json` names the two baseline frequency-error means in the wrong order
(they are 0.001634 for D-STAR and 0.003307 for rate 33) and quotes the P25 balance error of the
free alignment (2.05 dB; 2.0999 dB at the joint lags); the checks and their limits are as declared,
and the file was not changed.

The first build of the P25 tone path (`d245efd`, one limit of 0.375 for single and dual tones) was
scored in the same way first. D-STAR and rate 33 came out as above, `tone-model-tones` failed the
same five checks and `tone-model-channel` passed, but the P25 checks of `tone-model-development`
failed five times (the waveform correlation's mean, -0.004377, and every vector, `noisy` -0.0585;
the distance from DVSI's decode, mean +0.02598 dB and `noisy` +0.1005; the envelope correlation,
-0.001189) and those of `tone-model-validation` four times (the waveform correlation's mean,
-0.007197, and every vector, `t03` -0.0275; the distance, +0.01615 and +0.0318 dB). The single-tone
limit of [P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames) followed from that, as
`b58d91d`, before any of this documentation, and was scored against the same gate sets.

**The final partition.** The gate set `tone-model-final` was declared with the others in `c736edf`,
before any build was scored: the checks of `tone-model-development` on `irstia` and `mark` (in
every mode the waveform correlation with DVSI's decode falls by no more than 0.0002 on the mean and
0.001 on any vector, the distance from DVSI's decode rises by no more than 0.01 dB on the mean and
0.03 dB on any vector, the envelope correlation falls by no more than 0.001, the distance from the
input rises by no more than 0.01 dB, the level against DVSI's moves by at most 0.05 dB, and every
encoder's bits are unchanged). The final partition was scored once for this change set, on
2026-11-26, after the review: `fac0111` against the branch at `c2ba220`, the reviewed code (the
commits after it change documents only), on `irstia` and `mark`, both builds aligned freely first
and then scored at their joint lags, both runs' encodings decoded by the candidate's decoder. It
passes all 24 checks of `tone-model-final` (none fails). The library's decode of DVSI's bits,
`fac0111` → `c2ba220`, partition means over the two vectors and the waveform correlation of each:

| Codec | Waveform corr. vs DVSI, sample for sample, mean | `irstia` | `mark` | LSD vs DVSI (dB), mean | LSD vs input (dB), mean |
|---|---|---|---|---|---|
| D-STAR | 0.9656 → 0.9830 | 0.9426 → 0.9774 | 0.9886 → 0.9886 | 1.7814 → 1.7794 | 5.5438 → 5.5432 |
| P25 | 0.9798 → 0.9798 | 0.9727 → 0.9727 | 0.9869 → 0.9869 | 1.5700 → 1.5696 | 5.3702 → 5.3703 |
| Rate 33 | 0.9766 → 0.9766 | 0.9697 → 0.9697 | 0.9834 → 0.9834 | 1.5592 → 1.5592 | 5.5208 → 5.5208 |

The envelope correlation with DVSI's decode stays at 0.9766 (D-STAR), 0.9840 (P25) and 0.9781
(rate 33), means. The gate's paired deltas: the waveform correlation +0.01738 for D-STAR on the
mean (`irstia` +0.0348, `mark` 0), +8.4e-07 for P25 (`irstia` +1.7e-06) and 0 for rate 33; the
distance from DVSI's decode -0.0020 dB (D-STAR), -0.0003 dB (P25) and 0; the distance from the
input -0.0006 dB, +0.0001 dB and 0; the level against DVSI's +0.0010 dB, -0.000006 dB and 0; and
every encoder's bits unchanged (no vector with changed rows in any mode). What this shows: on the
two held-out vectors no metric of the gate set moves the wrong way past its limit, and D-STAR's
decode follows DVSI's more closely, a gain that is all `irstia`'s (+0.0348 in waveform correlation;
`mark` does not move). What it does not show: it is two vectors, and P25 and rate 33 do not move
there (rate 33's two decodes are unchanged, P25's correlation by at most 1.7e-06), so it says
nothing for or against the change set in those two codecs; every `mark` figure is the same at four
decimals.

### November 2026: mid-frame voicing

The tone model change's code against the two decoder changes of the mid-frame voicing change set
(branch `fix/mid-frame-voicing`, from `6acba79`):

- [P25: the bands that start with the frame](#p25-the-bands-that-start-with-the-frame): after a
  frame with no voiced band, a P25 (and ProVoice) frame sounds at mid-frame only the voiced
  harmonics in, or within 250 Hz below, a 500 Hz column its voicing covers.
- [Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17): these V/UV rows voice the first half
  of the frame otherwise than its end, and row 17's own voicing is Annex B's, all unvoiced.

The change set was written on `cc272ef`, the tone model change before its review's fixes, and the
gate sets name `cc272ef`, the branch point when they were declared, as the baseline. It was then
moved onto the tone branch's final code (`c2ba220`; the branch's head, `adecb29`, has the same code)
and scored there, `c2ba220` against the change set on it. The two baselines decode alike: all 255
decode hashes and 75 encoder-row hashes of `c2ba220` equal `cc272ef`'s, and every partition mean
and vector hash of that run equals the run made before the move. After its own review the change
set took two fixes: each frame mark is read only in the rendering of the codec whose decoder writes
it, P25 and ProVoice rendering with one of their own ([Rate 33: rows 13, 15 and
17](#rate-33-rows-13-15-and-17), [P25: the bands that start with the
frame](#p25-the-bands-that-start-with-the-frame)), and the P25 rendering marks the model it renders,
after adaptive smoothing. The run reported here was made with them, `adecb29` against the change
set with the fixes (`201bafd`'s code). Both builds aligned freely first and were then scored at
their joint lags, both runs' encodings decoded by the candidate's decoder. The study
that found these rules read the development and validation vectors and DVSI's tone vectors, never
`irstia` or `mark`; the change set was measured and gated on the development, validation, tone and
channel partitions, and the final partition was held out and scored once, after the review
([Method](#method); below).

**Speech**, decoding DVSI's bits: the waveform correlation with DVSI's decode, sample for sample
(`vs_dvsi.waveform_corr_unshifted`), partition means with the mean paired change and its 95%
interval, the change in log-spectral distance from DVSI's decode (`vs_dvsi.lsd_db`), and the vectors
whose correlation rises and falls:

| Codec | Partition | Waveform corr. vs DVSI | LSD vs DVSI (dB) | Vectors up / down |
|---|---|---|---|---|
| P25 | development | 0.9687 → 0.9849 (+0.0161, 95% interval 0.0088 to 0.0192) | -0.035 | 20 / 0 (largest rise `lvl` +0.088) |
| P25 | validation | 0.9855 → 0.9896 (+0.0041, interval 0.00004 to 0.0091) | -0.023 | 5 / 0 |
| Rate 33 | development | 0.9737 → 0.9806 (+0.0069, interval 0.0041 to 0.0112) | -0.023 | 19 / 1 (`ucarm20a` -0.0006; largest rise `fambm22a` +0.040) |
| Rate 33 | validation | 0.9788 → 0.9834 (+0.0046, interval 0.00006 to 0.0097) | -0.026 | 5 / 0 |

The partition means of the envelope correlation with DVSI's decode and of the log-spectral distance
from the input improve in both codecs and both partitions. Per vector, the envelope correlation
falls on six rate-33 vectors (development `dam` -6.2e-05, `dam_10ov` -4.1e-05, `lvl` -2.1e-05,
`dam80` -1.2e-05 and `t01` -1e-06; validation `t03` -1.0e-05), and the distance from the input
rises on P25's `lvl` (+0.0014 dB) and `t04` (+0.0007 dB) and on rate 33's `t04` (+0.0029 dB),
`dam_10ov` (+0.0009 dB), `dam80` (+0.0008 dB) and `dam40` (+0.0002 dB). The level against DVSI's
moves by under 0.0006 dB on any vector. The correlation of the script in [P25: the bands that start
with the frame](#p25-the-bands-that-start-with-the-frame) and the scoreboard's differ (0.9848 and
0.9849 for P25 on development) because the scoreboard correlates only the frames where DVSI's
decode lies within 40 dB of its loudest frame. D-STAR's decodes are identical, and every encoder's
bits are unchanged in all three modes.

**Where the residual sits**, on the development vectors: the energy of the difference between our
decode and DVSI's over the signal's, by kind of frame, at `cc272ef` and with the change (the
scratch builds, whose decodes of these vectors the library's equal). The accounting is its own, not
the mean over the vectors of each vector's ratio, which [P25: the bands that start with the
frame](#p25-the-bands-that-start-with-the-frame) and [Rate 33: rows 13, 15 and
17](#rate-33-rows-13-15-and-17) give: every vector is scaled to unit power and the residual and the
signal are summed over the vectors; a voiced run is a stretch of frames with a voiced band in the
decode before the change, in both rows of a codec (so a rate-33 row-17 frame counts as voiced in
both, and the run counts are not the change's own runs), and a run is under 0.8 by its whole-run
correlation; a file's first frame and a run under way at its start are in the total and in no
part, so the parts need not add up to it.

| Codec | Build | Residual | Voiced runs in line | Voiced runs under a correlation of 0.8 | Frames with no voiced band |
|---|---|---|---|---|---|
| D-STAR | `cc272ef` and with the change | 0.0204 | 0.0172 | 0.0008 | 0.0024 |
| P25 | `cc272ef` | 0.0647 | 0.0211 | 0.0407 (82 of 1,959 runs) | 0.0029 |
| P25 | with the change | 0.0332 | 0.0211 | 0.0094 (36 runs) | 0.0027 |
| Rate 33 | `cc272ef` | 0.0603 | 0.0331 | 0.0246 | 0.0025 |
| Rate 33 | with the change | 0.0449 | 0.0332 | 0.0091 (25 of 1,825 runs) | 0.0025 |

On the validation vectors P25's residual with the change is 0.0205 (0.0155 + 0.0024 + 0.0027; 12 of
542 runs under 0.8). In both P25 and rate 33 the voiced runs in line hold about what they held
before, and the runs out of line give up most of theirs. Inside rate 33's runs in line the residual
grows with the run's length, from 0.0175 of the signal in its frames 0 to 3 to 0.0675 from frame 20
on (the same grouping).

**Tones.** The decodes of 9 P25 tone vectors and 1 rate-33 tone vector change. The five tone
measures the gates check (the frequency, level and balance errors, the share of windows playing a
wrong tone and the detector's agreement) do not move on any vector; the scoreboard's
`tone.alignment_corr`, which no gate checks, changes on nine P25 vectors and one rate-33 vector, by
at most 0.000113 (P25 `cp0`, -0.000113). Frame for frame against DVSI's decode, as measured in
[November
2026: tones and the octave step](#november-2026-tones-and-the-octave-step), the steady tone frames'
median correlation stays at 0.9997 (rate 33) and 0.9862 (P25) and the first tone frames' median
residual at -31.6 and -13.0 dB (a P25 tone frame keeps its own model at mid-frame: [P25: the bands
that start with the frame](#p25-the-bands-that-start-with-the-frame)).

**Channel errors.** No mute decision changes at any error rate, hard or soft, frame for frame: the
mute flag of every frame of the 24 channel-error decodes (8 vectors in each of the three codecs,
30,600 frames, 1,678 of them muted) is the same in both builds. Every rate-33 channel-error decode
changes. The scoreboard's distance from DVSI's decode on the loud frames is taken over the frames
that neither decoder nor DVSI mutes, away from a muted frame, and measured in both runs; there it
falls on every vector, means over the error rates, by 0.0137 dB with hard decisions (by 0.0095 to
0.0177 dB per vector) and 0.0140 dB with soft ones (0.0026 to 0.0178 dB). Those frames are 0.750
(hard) and 0.991 (soft) of the loud frames DVSI plays on `dam_e10`, the lowest coverage, and all of
them on the other six vectors. No D-STAR or P25 channel-error decode changes.

**PESQ** ([Method](#method)) of the decode of DVSI's bits on the 25 development and validation
vectors (mean paired change ± standard error clustered by recording):

| Codec | Against the input | Against DVSI's decode |
|---|---|---|
| P25 | 3.033 → 3.035, +0.0017 ± 0.0010 (16 of 25 vectors better; DVSI's own decode 3.026) | 4.274 → 4.278, +0.0033 ± 0.0007 (23 of 25) |
| Rate 33 | 2.953 → 2.956, +0.0023 ± 0.0013 (19 of 25; DVSI's own decode 2.980) | 4.238 → 4.244, +0.0052 ± 0.0009 (24 of 25) |

**Cost.** No measurable cost: every codec's median complete frame decode moves by -1.3% to +0.3%,
`cc272ef` against the change set on it, measured on a busy machine (load 1.5 to 2.5;
[benchmarks](testing.md#mid-frame-voicing)).

**Gate sets**, declared in the change set's first commit, `6acba79`, before any build was scored,
and evaluated in one run over the development, validation, tone and channel partitions, `adecb29`
against the change set with the review's fixes:

| Gate set | Result |
|---|---|
| `mid-voicing-development` | passes every check |
| `mid-voicing-validation` | passes every check |
| `mid-voicing-tones` | passes every check |
| `mid-voicing-channel` | passes every check |
| `mid-voicing-final` | passes (scored once, after the review, at `504740c`; below) |

Every check passed with the same values before the fixes, `c2ba220` against the change set, but
one: P25's soft-decision distance on the loud frames of the channel set (above).

Notes on the declarations: what each check establishes, and what was verified beside it. A check
that the mean rises by at least a value is the scoreboard's `increase` rule, which also requires
the 95% interval above 0, stricter than a rise of the mean alone; the sets' descriptions say so.
D-STAR's unchanged decode is checked as no vector with changed PCM, over the speech vectors and,
in the tone and channel sets, over their own vectors too. The
encoder checks of the development, validation and channel sets count the speech vectors whose
encoder rows change, and those of the tones set compare the D-STAR and rate-33 encoders' tone
classification (agreement with DVSI's tone frames and the rate of tones where DVSI sends voice), not
their bits; beside the gates every encoder's bits were compared directly: all 285 encodings (95
input files, `irstia` and `mark` not among them, three encoders) are identical between the two
builds, and so are all 106 D-STAR decodes of DVSI's bits (4 of them soft-decision). In the tones
set, `cp1`, `cpvbad` and `dtmfvbad` have no steady tone window in DVSI's decode, so neither run
measures the frequency, level, balance and wrong-tone figures there, and `alert` has no balance,
being a single tone (it has steady windows: 128 in P25); the set's description records the coverage,
from the baseline's own scoreboard, before scoring. The channel set's "no mute decision changes" is
expressed as zero change, on every vector, of the mute agreement, the mute rate and the false and
missed mute rates on the loud frames, hard and soft: statistics per vector, which two changes of
opposite sign could leave as they were; the frame-level comparison above, beside it, finds no mute
decision changed. The sets' descriptions name `cc272ef` as the baseline (above).

Two statements of the declarations are wrong, and the file was not changed (the checks and their
limits are as declared). The description of `mid-voicing-tones` puts `alert` among the vectors
whose DVSI decode has no steady window holding a tone; P25's `alert` has steady single-tone
windows and lacks only a dual tone's balance (as above). And every set's description says DVSI's
decoders sound "a P25 (and ProVoice) voice frame" so; there is no ProVoice reference ([P25
pitch](#p25-pitch)),
and ProVoice follows P25's rule because it shares P25's data path.

**The final partition.** The gate set `mid-voicing-final` was declared with the others in `6acba79`,
before any build was scored: on `irstia` and `mark`, D-STAR's decode of DVSI's bits unchanged, every
encoder's bits unchanged, and for P25 and rate 33 the waveform correlation with DVSI's decode
falling by no more than 0.0005 on the mean and 0.002 on any vector, the distance from DVSI's decode
rising by no more than 0.01 dB on the mean and 0.03 dB on any vector, the envelope correlation
falling by no more than 0.001, the distance from the input rising by no more than 0.01 dB, and the
level against DVSI's moving by at most 0.05 dB. The final partition was scored once for this change
set, on 2026-11-26, after its review: the tone model change's code (the branch's build of
`c2ba220`, whose decoder code `adecb29` shares) against `504740c`, the reviewed code and documents,
on `irstia` and `mark`, both builds aligned freely first and then scored at their joint lags, both
runs' encodings decoded by the candidate's decoder. It passes all 18 checks of `mid-voicing-final`
(none fails). The library's decode of DVSI's bits, partition means over the two vectors and the
waveform correlation of each:

| Codec | Waveform corr. vs DVSI, sample for sample, mean | `irstia` | `mark` | LSD vs DVSI (dB), mean | LSD vs input (dB), mean |
|---|---|---|---|---|---|
| P25 | 0.9798 → 0.9837 | 0.9727 → 0.9795 | 0.9869 → 0.9880 | 1.5696 → 1.5248 | 5.3703 → 5.3690 |
| Rate 33 | 0.9766 → 0.9783 | 0.9697 → 0.9715 | 0.9834 → 0.9850 | 1.5592 → 1.5328 | 5.5208 → 5.5211 |
| D-STAR | 0.9830 → 0.9830 | 0.9774 → 0.9774 | 0.9886 → 0.9886 | 1.7794 → 1.7794 | 5.5432 → 5.5432 |

The gate's paired deltas: the waveform correlation +0.003939 for P25 on the mean (`irstia`
+0.006790, `mark` +0.001088) and +0.001702 for rate 33 (+0.001795, +0.001608); the distance from
DVSI's decode -0.04488 dB for P25 (-0.0443 and -0.0455 dB) and -0.02643 dB for rate 33 (-0.0228 and
-0.0301 dB); the envelope correlation +0.00031 for P25 (0.9840 → 0.9843) and +0.00022 for rate 33
(0.9781 → 0.9783); the distance from the input -0.0013 dB for P25 (`irstia` 5.3272 → 5.3236, `mark`
5.4134 → 5.4145) and +0.0003 dB for rate 33 (`irstia` 5.5626 → 5.5614, `mark` 5.4790 → 5.4807); the
level against DVSI's -0.0009 dB (P25) and -0.0011 dB (rate 33); D-STAR's decode unchanged (no
vector with changed PCM), and every encoder's bits unchanged (no vector with changed rows in any
mode). What this shows: on the two held-out vectors the waveform correlation with DVSI's decode
rises on both vectors in both codecs, more on `irstia` (P25 +0.0068) than on `mark`, and the
distance from DVSI's decode falls on both, as on the development and validation partitions. What
it does not show: it is two vectors, a check that nothing moves the wrong way past its limit
rather than a measure of the change's size; and the distance from the input rises on `mark` in
both codecs, by 0.001 to 0.002 dB, and falls on `irstia`, within the gate's limits but not an
improvement.

### November 2026: rate 33's mid-frame pitch wave

The mid-frame voicing change's head (`dac6d14`, PR #198) against rate 33's mid-frame pitch wave
([Rate 33: mid-frame pitch](#rate-33-mid-frame-pitch)): between two voiced rate-33 frames of
different pitch the mid-frame pitch is the geometric mean times `1 - 3.0e-4` times `1 + (h(Pa) +
h(Pb)) / 2`. The study that found the wave read the development and validation vectors (the wave's
constants were fitted on the development vectors alone) and, for the checks of the change, DVSI's
rate-33 tone vectors, never `irstia` or `mark`. Its scratch build of the change was measured on the
development and validation vectors (the whole-file figures, the accounting and PESQ below) before
the gate sets were declared, and their limits were chosen knowing those results, so the development
and validation partitions are not out-of-sample here; the final partition was held out and scored
once, after the review (below). The candidate is the change set's code at `d08ea39` (the
commits after it change tests, comments and documents only). Both builds aligned freely first and
were then scored at their joint lags, both runs' encodings decoded by the candidate's decoder.

**Speech**, rate 33 decoding DVSI's bits, on the scoreboard: the waveform correlation with DVSI's
decode, sample for sample (`vs_dvsi.waveform_corr_unshifted`, over the frames where DVSI's decode
lies within 40 dB of its loudest), with the mean paired change and its 95% interval and the
vectors whose correlation rises and falls; the log-spectral distance from DVSI's decode
(`vs_dvsi.lsd_db`), the envelope correlation with it (`vs_dvsi.env_corr`), the level against it
(`vs_dvsi.level_offset_db`) and the log-spectral distance from the input (`ours.lsd_db`), partition
means:

| Partition | Waveform corr. vs DVSI | Vectors up / down | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|---|
| development | 0.9806 → 0.9893 (+0.0087, 95% interval 0.0052 to 0.0135) | 20 / 0 (smallest rise `dam10` +0.0030, largest `fambm22a` +0.0267) | 1.2618 → 1.2299 | 0.99205 → 0.99207 | +0.0139 → +0.0139 | 6.2696 → 6.2692 |
| validation | 0.9834 → 0.9921 (+0.0087, interval 0.0031 to 0.0145) | 5 / 0 (smallest `tambf32e` +0.0027, largest `tambf22b` +0.0167) | 1.0797 → 1.0642 | 0.98839 → 0.98840 | -0.0403 → -0.0403 | 4.8304 → 4.8303 |

Per vector, the distance from DVSI's decode falls on all 25 (by 0.0063 to 0.1361 dB on development,
`dam40` and `lvl`, and 0.0108 to 0.0179 dB on validation, `t03` and `t04`); the envelope
correlation falls on one, `t04` (-7.7e-06); the distance from the input rises on 13, by at most
0.0017 dB (`dam80`; `t02` +0.0011, `tambf32e` +0.0011, `fambf22c` +0.0010, `tia11` +0.0010, the
others +0.0004 or less), and falls on 12, by at most 0.0046 dB (`dam_10ov`); the level against
DVSI's moves by at most 0.0013 dB on any vector (`dam_10ov`). D-STAR's and P25's partition means
and decodes are unchanged.

Over whole files, the measure of the study: the correlation of our decode with DVSI's, sample for
sample over the whole file, and the residual energy of the difference over the signal's, means over
the vectors, the vectors whose correlation moves by more than 1e-6:

| Rate 33 | Correlation | Residual | Up / down |
|---|---|---|---|
| development, `dac6d14` | 0.9806 | 0.0388 | |
| development, the wave | 0.9893 | 0.0214 | 20 / 0 (smallest `dam10` +0.0030, largest `fambm22a` +0.0267) |
| validation, `dac6d14` | 0.9834 | 0.0331 | |
| validation, the wave | 0.9921 | 0.0157 | 5 / 0 (smallest `tambf32e` +0.0027, largest `tambf22b` +0.0167) |

These were measured on the scratch build of the study and again on the library's builds of `dac6d14`
and `d08ea39`, which give the same figures to four decimals. The scratch build forms the geometric
mean in float before the factor, the library in double: of the 25 rate-33 decodes none is identical
between the two, and they differ in 2 (`dam80`) to 3,767 (`lvl`) samples of a file, by one 16-bit
step but on `lvl` (up to 4), `dam_10ov` and `fambm22a` (up to 2). On the scratch build, the
five-term development fit (with the second harmonic) gives the same means at four decimals (0.9893
and 0.0214 on development, 0.9921 and 0.0157 on validation); C1 and S1 times 1.19 give 0.9891 /
0.0219 and 0.9924 / 0.0153, times 0.83 0.9890 / 0.0220 and 0.9914 / 0.0171, and the wave's opposite
sign 0.9562 and 0.9664.

**Where the residual sits**, in the accounting of [November 2026: mid-frame
voicing](#november-2026-mid-frame-voicing): every vector scaled to unit power, the residual and the
signal summed over the vectors, the voiced runs as that section defines them, a run counting as in
line or as under a correlation of 0.8 by its own correlation in each decode (so the two rows need
not classify the same runs), and the frames of the runs in line grouped by their position in the
run. Measured on the scratch builds, and again on the library's builds of `dac6d14` and `d08ea39`,
which give the same figures (the validation row with the mid-frame voicing change comes from the
library's build alone):

| Rate 33 | Residual | Voiced runs in line | Voiced runs under a correlation of 0.8 | Frames with no voiced band |
|---|---|---|---|---|
| development, with the mid-frame voicing change | 0.0449 | 0.0332 | 0.0091 | 0.0025 |
| development, with the wave | 0.0248 | 0.0138 | 0.0085 | 0.0025 |
| validation, with the mid-frame voicing change | 0.0358 | 0.0182 | 0.0137 | 0.0039 |
| validation, with the wave | 0.0163 | 0.0124 | 0.0000 | 0.0039 |

By frames into a run in line (frames 0 to 3, 4 to 7, 8 to 11, 12 to 15, 16 to 19, 20 on), on the
development vectors the residual was 0.0175, 0.0243, 0.0358, 0.0517, 0.0614 and 0.0675 of the signal
and is 0.0155, 0.0128, 0.0118, 0.0148, 0.0157 and 0.0185: it no longer grows with the run's length
there. On the validation vectors it was 0.0164, 0.0125, 0.0172, 0.0269, 0.0297 and 0.0375 and is
0.0157, 0.0075, 0.0076, 0.0122, 0.0125 and 0.0333: lower in every group, the group from frame 20 on
still the largest.

**PESQ** ([Method](#method)) of the decode of DVSI's bits on the 25 development and validation
vectors, the scratch builds (mean paired change ± standard error clustered by recording): against
the input 2.956 → 2.959, +0.0035 ± 0.0018 (13 of 25 vectors better), within about two standard
errors; against DVSI's decode 4.244 → 4.244, +0.0006 ± 0.0004 (14 of 25).

**Tones.** The decodes of 5 rate-33 tone vectors change (`alltone`, `cpvbad`, `dtmf`, `dtmf180`,
`dtmfvbad`). Of the five tone measures the gates check, only the detector's agreement with DVSI's
tone output moves, on `cpvbad`, by one of its 359 windows (0.9944 → 0.9916); the scoreboard's
`tone.alignment_corr`, which no gate checks, moves on the five by under 1e-6. Frame for frame
against DVSI's decode, as measured in [November 2026: tones and the octave
step](#november-2026-tones-and-the-octave-step), the steady tone frames' median correlation stays
at 0.9997 and the first tone frames' median residual at -31.6 dB (the scratch builds). No D-STAR or
P25 tone decode changes.

**Channel errors.** No mute decision changes at any error rate, hard or soft, frame for frame: the
mute flag of every frame of the 24 channel-error decodes (8 vectors in each of the three codecs,
30,600 frames, 1,678 of them muted) is the same in both builds. Every rate-33 channel-error decode
changes, and on the loud frames neither decoder mutes the log-spectral distance from DVSI's decode
falls, means over the error rates, by 0.0137 dB with hard decisions (on every vector, by 0.0095 to
0.0177 dB; on 0.750 of the loud frames DVSI plays) and 0.0140 dB with soft ones (0.0026 to 0.0178
dB; 0.991). No D-STAR or P25 channel-error decode changes.

**Other modes and the encoders.** Beside the scoreboard's hashes (25 speech, 52 tone and 8
channel-error vectors per mode), every D-STAR and P25 bit file of DVSI's vectors at the top of each
mode's directory, `irstia` and `mark` left out (106 D-STAR and 103 P25 files, the DTX, synthetic and
channel-error ones among them, the soft-decision ones decoded as soft), decodes identically in four
builds of `dac6d14` and `d08ea39` each: the scoreboard's (Release, SIMD, IEEE arithmetic), GCC
Release with fast-math and LTO for the build host, Clang Release with fast-math and LTO for
x86-64-v3, and 32-bit Release with x87 arithmetic. In the same four builds all 285 encodings (95
input files, `irstia` and `mark` not among them, three encoders) are identical. Of the 106 rate-33
bit files, 52 decode otherwise.

**Cost.** No measurable cost: the median complete hard-decision frame decode moves by +0.9%
(AMBE+2), -0.4% (D-STAR) and +0.4% (IMBE, neither of which the change touches), `dac6d14` against
`d08ea39`, with the ranges of the invocations overlapping in every codec, measured on a busy
machine (load 1.7 to 2.0; [benchmarks](testing.md#rate-33s-mid-frame-pitch-wave)).

**Gate sets**, declared in the branch's first commit, `8942924`, before any build of the branch was
scored (the study's scratch build had been scored on the development and validation vectors before
that, above), and evaluated in one combined run over the development, validation, tone and channel
partitions, `dac6d14` against `d08ea39`:

| Gate set | Result | What it establishes |
|---|---|---|
| `r33-pitch-wave-development` | passes every check | on the vectors the constants were fitted on: rate 33's correlation rises (mean +0.0087, interval above 0, no vector falling by more than 0.0005), the spectral, envelope, input and level measures stay within their limits, D-STAR's and P25's decodes of the speech vectors and every encoder's rows on them are unchanged |
| `r33-pitch-wave-validation` | passes every check | the same on vectors held out from the fit of the constants (read in the study) |
| `r33-pitch-wave-tones` | passes every check | rate 33's tone measures within their limits, D-STAR's and P25's speech and tone decodes unchanged, the D-STAR and rate-33 encoders' per-vector tone statistics unchanged |
| `r33-pitch-wave-channel` | passes every check | the per-vector mute statistics unchanged in every mode, rate 33's mean distance from DVSI's decode on the loud frames rising by no more than 0.01 dB, D-STAR's and P25's decodes and every encoder's rows unchanged |
| `r33-pitch-wave-final` | passes all 12 checks (scored once, after the review; below) | on `irstia` and `mark`, held out: nothing moves the wrong way past its limit (below) |

Notes on the declarations. A check that the mean rises by at least a value is the scoreboard's
`increase` rule, which also requires the 95% interval above 0. D-STAR's and P25's unchanged decodes
are checked as no vector with changed PCM, over the speech vectors and, in the tone and channel
sets, over their own vectors too. These identity counters, and the encoder-row counter, are counted
over the speech vectors present in the run: they were evaluated on the combined run, and a tone or
channel partition scored alone leaves some of them without a value. The tones group holds no
encoder-row hashes, so the tones set compares the D-STAR and rate-33 encoders' per-vector tone
statistics (agreement with DVSI's tone frames, tones where DVSI sends voice, and the level
difference), which do not show that the same frames or bits were sent; the bit identity is the
separate comparison of the 285 encodings above. The channel set's mute checks are statistics per
vector, and the frame-level comparison above, beside them, finds no mute decision changed. The
per-vector checks of the distance from DVSI's decode (`vs_dvsi.lsd_db`, `all_no_increase`) do not
require every vector, so they would pass a vector that neither run measures; in the recorded run all
20 development and all 5 validation vectors were measured, and in the final partition's run both
vectors were measured in both runs for every gated rate-33 measure (below). The tones set's description records, from the baseline's own
scoreboard of the tone partition before the set was declared, which tone measures each vector can be
measured on in rate 33: the detector agreement and the extra-tone rate on all 52; the frequency,
level and wrong-tone figures on all but `cp1`, `cpvbad` and `dtmfvbad`; the balance on all but those
and `alert`, a single tone; the rate-33 encoder agreement and level error on all but `cpvbad`.

**The final partition.** The gate set `r33-pitch-wave-final` was declared with the others in
`8942924`, before any build of the branch was scored: on `irstia` and `mark`, D-STAR's and P25's
decodes of DVSI's bits unchanged, every encoder's bits unchanged, and for rate 33 the waveform
correlation with DVSI's decode falling by no more than 0.0005 on the mean and 0.002 on any vector,
the distance from DVSI's decode rising by no more than 0.01 dB on the mean and 0.03 dB on any
vector, the envelope correlation falling by no more than 0.001, the distance from the input rising
by no more than 0.01 dB, and the level against DVSI's moving by at most 0.05 dB. The final
partition was scored once for this change set, on 2026-11-26, after its review: `dac6d14` against
the head's code (the build of `d08ea39`, whose code the head shares apart from comments), on
`irstia` and `mark`, both builds aligned freely first and then scored at their joint lags, both
runs' encodings decoded by the candidate's decoder. It passes all 12 checks of
`r33-pitch-wave-final` (none fails). The library's decode of DVSI's bits, partition means over the
two vectors and the waveform correlation of each:

| Codec | Waveform corr. vs DVSI, sample for sample, mean | `irstia` | `mark` | LSD vs DVSI (dB), mean | LSD vs input (dB), mean |
|---|---|---|---|---|---|
| Rate 33 | 0.9783 → 0.9860 | 0.9715 → 0.9828 | 0.9850 → 0.9891 | 1.5328 → 1.5074 | 5.5211 → 5.5216 |
| P25 | 0.9837 → 0.9837 | 0.9795 → 0.9795 | 0.9880 → 0.9880 | 1.5248 → 1.5248 | 5.3690 → 5.3690 |
| D-STAR | 0.9830 → 0.9830 | 0.9774 → 0.9774 | 0.9886 → 0.9886 | 1.7794 → 1.7794 | 5.5432 → 5.5432 |

The gate's paired deltas for rate 33: the waveform correlation +0.007689 on the mean (`irstia`
+0.011249, `mark` +0.004128); the distance from DVSI's decode -0.02537 dB (`irstia` 1.7053 → 1.6848,
`mark` 1.3603 → 1.3301 dB); the envelope correlation +0.00002 (0.9783 → 0.9783); the distance from
the input +0.0005 dB (`irstia` 5.5614 → 5.5624, `mark` 5.4807 → 5.4808 dB); the level against DVSI's
-0.0004 dB (-0.0727 → -0.0731); D-STAR's and P25's decodes unchanged (no vector with changed PCM in
either), and every encoder's bits unchanged (no vector with changed rows in any mode). Both vectors
were measured in both runs for every gated rate-33 measure, so the per-vector check of the distance
from DVSI's decode, which does not itself require every vector, compared both.

What this shows: it is the out-of-sample check of the change, since the constants were fitted on the
development vectors and the validation vectors were read in the study; on the two held-out vectors
the waveform correlation with DVSI's decode rises on both, more on `irstia` than on `mark`, by
+0.0077 on the mean, a little under the development and validation partitions' +0.0087 each, and the
distance from DVSI's decode falls on both. What it does not show: it is two vectors, a check that
nothing moves the wrong way past its limit rather than a measure of the change's size; and the
distance from the input rises on both vectors, by 0.0001 dB (`mark`) to 0.0010 dB (`irstia`), within
the gate's limit but not an improvement.

### November 2026: mid-frame voicing by columns

Rate 33's mid-frame pitch wave's head (`c3c6d42`, PR #199, whose decoder code is `d08ea39`'s)
against the mid-frame columns change ([D-STAR and P25: mid-frame voicing by
columns](#d-star-and-p25-mid-frame-voicing-by-columns)): the mid-frame model of a D-STAR or P25
(and ProVoice) frame voices each band by one rule of 500 Hz columns, read in both frames, in place
of this frame's voicing between voiced frames, the onset rules and the hold through a frame lacking
a pitch. The study that found the rule read the development and validation vectors (its labels),
DVSI's tone, channel-error and DTX vectors and the DVstick 30 curves of
[#157](https://github.com/arancormonk/mbelib-neo/issues/157), never `irstia` or `mark`. Its scratch
build was measured on the development and validation vectors (the whole-file figures, the
accounting and PESQ below) before the gate sets were declared, and the rule was found on those
vectors' labels, so the development and validation partitions are not out of sample here; the final
partition was held out and scored once, after the review (below). The library's D-STAR decodes equal the scratch build's with the rule,
and its P25 decodes those of the scratch build with the rule applied after a tone model too, on all
25 development and validation vectors (16-bit samples, every file identical), and on every D-STAR
and P25 bit file of DVSI's vectors at the top of each mode's directory, `irstia` and `mark` left out
(106 D-STAR and 103 P25 files, the tone, channel-error, DTX and synthetic ones among them, the
soft-decision ones decoded as soft); so the scratch build's figures below are the library's. The
review of the first form of the change (`6f8ad6d`) found two defects, fixed in `e86afd0`: each
frame's nearest harmonic is now read in double, and the rule applies after a P25 tone model
([D-STAR and P25: mid-frame voicing by
columns](#d-star-and-p25-mid-frame-voicing-by-columns)). The first fix changes no decode of DVSI's
vectors; the second changes three P25 tone decodes (`dtmf60`, `dtmf180`, `dtmfvbad`) and no speech
decode, on the files compared (the 106 D-STAR and 103 P25 bit files above). The scoreboard figures
below are those of the second run, `c3c6d42` against `e86afd0`; the
first run, against `e9cf77c`, gave the same figures but for those three tone vectors. Both builds
aligned freely first and were then scored at their joint lags, both runs' encodings decoded by the
candidate's decoder.

**Speech**, decoding DVSI's bits, on the scoreboard: the waveform correlation with DVSI's decode,
sample for sample (`vs_dvsi.waveform_corr_unshifted`, over the frames where DVSI's decode lies within
40 dB of its loudest), with the mean paired change and its 95% interval and the vectors whose
correlation rises and falls; the log-spectral distance from DVSI's decode (`vs_dvsi.lsd_db`), the
envelope correlation with it (`vs_dvsi.env_corr`), the level against it (`vs_dvsi.level_offset_db`)
and the log-spectral distance from the input (`ours.lsd_db`), partition means:

| Codec | Partition | Waveform corr. vs DVSI | Vectors up / down | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|---|---|
| D-STAR | development | 0.9907 → 0.9916 (+0.0010, 95% interval 0.0005 to 0.0016) | 20 / 0 (smallest rise `dam40` +0.0001, largest `ucarm15a` +0.0051) | 1.3533 → 1.1607 | 0.98816 → 0.98861 | +0.0247 → +0.0255 | 6.4257 → 6.4229 |
| D-STAR | validation | 0.9820 → 0.9935 (+0.0115, interval 0.0008 to 0.0400) | 5 / 0 (smallest `tambf32b` +0.0007, largest `t03` +0.0530) | 1.2378 → 1.0564 | 0.98501 → 0.98550 | -0.0053 → -0.0043 | 4.8467 → 4.8433 |
| P25 | development | 0.9849 → 0.9867 (+0.0018, interval 0.0008 to 0.0045) | 20 / 0 (smallest `t02` +0.00004, largest `noisy` +0.0148) | 1.6111 → 1.4750 | 0.98937 → 0.98965 | +0.0453 → +0.0465 | 6.3174 → 6.3169 |
| P25 | validation | 0.9896 → 0.9908 (+0.0012, interval 0.0005 to 0.0018) | 5 / 0 (smallest `t03` +0.0001, largest `t04` +0.0021) | 1.2462 → 1.1252 | 0.98751 → 0.98768 | -0.0143 → -0.0125 | 4.7259 → 4.7209 |

Per vector, the distance from DVSI's decode falls on all 25 in both codecs (D-STAR by 0.073 to
0.331 dB on development, `dam40` and `lvl`, and 0.126 to 0.238 dB on validation, `t03` and `t04`;
P25 by 0.030 to 0.247 dB, `dam80` and `lvl`, and 0.051 to 0.154 dB, `t03` and `tambf32e`), and the
envelope correlation rises on all 25. The distance from the input rises on 4 D-STAR development
vectors, by at most 0.0025 dB (`t02`; `fambm22a` +0.0017, `ucarm15a` +0.0013, `tambf22a` +0.0002),
and on `t03` of validation (+0.0033 dB), and falls on the other 20; in P25 it rises on 10
development vectors, by at most 0.0039 dB (`tia11`; `dam_10ov` +0.0026, `ucarm15a` +0.0025,
`fambm22a` +0.0023, `noisy` +0.0022, `dam10` +0.0018, the others +0.0009 or less), and on `t03`
(+0.0032 dB), and falls on the other 14. The level against DVSI's moves by at most 0.0038 dB on a
development vector (P25 `clean`) and 0.0051 dB on a validation vector (P25 `t04`). Rate 33's
partition means and decodes are unchanged, and so is every encoder's bits.

Over whole files, the measure of the study (the correlation of our decode with DVSI's, sample for
sample over the whole file at lag 0, and the residual energy of the difference over the signal's,
means over the vectors; the vectors whose correlation moves by more than 1e-6), measured on the
scratch build, whose decodes are the library's:

| Codec | Partition | Correlation | Residual | Up / down |
|---|---|---|---|---|
| D-STAR | development | 0.9907 → 0.9916 | 0.0191 → 0.0172 | 20 / 0 (smallest `dam40` +0.0001, largest `ucarm15a` +0.0051) |
| D-STAR | validation | 0.9819 → 0.9934 | 0.0363 → 0.0132 | 5 / 0 (smallest `tambf32b` +0.0007, largest `t03` +0.0530) |
| P25 | development | 0.9848 → 0.9866 | 0.0309 → 0.0273 | 20 / 0 (smallest `t02`, under +0.00005, largest `noisy` +0.0148) |
| P25 | validation | 0.9896 → 0.9907 | 0.0208 → 0.0185 | 5 / 0 (smallest `t03` +0.0001, largest `t04` +0.0021) |

Scratch variants each switched on alone, same measure, development / validation, the change of
the partition mean against the unswitched scratch build (to six decimals, then four), D-STAR
read on P25's grid (which end to end gives the same means as its own columns to four decimals):
frames lacking a pitch held by the covered columns of the frame before (and within 250 Hz below
one), D-STAR +0.000192 / +0.010679 (+0.0002 / +0.0107), P25 +0.000123 / +0.000011 (+0.0001 /
+0.0000); D-STAR onsets read in columns, +0.000141 / +0.000047 (+0.0001 / +0.0000). Inside voiced
runs the variants were not the shipped rule's terms: a changing band voiced at mid-frame when the
sample above its column, at either reach, changes, without the both-voiced term and without a
condition on the column's own sample: in both directions D-STAR +0.000586 / +0.000741 (+0.0006 /
+0.0007), P25 +0.000948 / +0.001095 (+0.0009 / +0.0011); for bands turning unvoiced only (those
turning voiced starting with the frame) D-STAR +0.000489 / +0.000653 (+0.0005 / +0.0007), P25
+0.000886 / +0.001086 (+0.0009 / +0.0011); for bands turning voiced only (those turning unvoiced
let go) D-STAR +0.000097 / +0.000089 (+0.0001 / +0.0001), P25 +0.000062 / +0.000009 (+0.0001 /
+0.0000). Variants that lose, on development:
steady voiced bands eroded too (a band voiced at mid-frame when covered in either frame), D-STAR
0.9902 (2 vectors up, 18 down), P25 0.9842 (4 up, 16 down); bands turning voiced sounding only if
covered in this frame, D-STAR 0.9906 (8 up, 11 down), P25 0.9847 (3 up, 17 down).

**Where the residual sits**, in the accounting of [November 2026: mid-frame
voicing](#november-2026-mid-frame-voicing) (every vector scaled to unit power, the residual and the
signal summed over the vectors; voiced runs with a whole-run correlation under 0.8 out of line),
the head and with the change, on the scratch builds:

| Codec | Partition | Residual | Voiced runs in line | Voiced runs under a correlation of 0.8 | Frames with no voiced band |
|---|---|---|---|---|---|
| D-STAR | development | 0.0204 → 0.0183 | 0.0172 → 0.0155 | 0.0008 → 0.0004 (29 → 19 of 1,911 runs) | 0.0024 → 0.0024 |
| D-STAR | validation | 0.0301 → 0.0138 | 0.0129 → 0.0115 | 0.0148 → 0.0000 (11 → 8 of 527) | 0.0024 → 0.0023 |
| P25 | development | 0.0332 → 0.0286 | 0.0211 → 0.0188 | 0.0094 → 0.0071 (36 → 27 of 1,959) | 0.0027 → 0.0027 |
| P25 | validation | 0.0205 → 0.0182 | 0.0155 → 0.0134 | 0.0024 → 0.0022 (12 → 11 of 542) | 0.0027 → 0.0027 |

D-STAR validation's 0.0148 is the sum over 11 runs out of line, of which the `t03` run from frame
895 (55 frames, correlation 0.08 with DVSI's) is 0.0147; that run's residual is 0.107 of `t03`'s
signal, and the run holds 0.057 of it. With the change that run is in line (its onset is one of the two
of [D-STAR and P25: mid-frame voicing by columns](#d-star-and-p25-mid-frame-voicing-by-columns),
"Which onsets showed it first").

**PESQ** ([Method](#method)) of the decode of DVSI's bits on the 25 development and validation
vectors, the scratch builds (mean paired change ± standard error clustered by recording):

| Codec | Against the input | Against DVSI's decode |
|---|---|---|
| D-STAR | 2.948 → 2.952, +0.0041 ± 0.0014 (19 of 25 vectors better; DVSI's own decode 2.933) | 4.283 → 4.323, +0.0397 ± 0.0089 (24 of 25) |
| P25 | 3.035 → 3.040, +0.0046 ± 0.0021 (16 of 25; DVSI's own decode 3.026) | 4.278 → 4.306, +0.0284 ± 0.0024 (25 of 25) |

**Tones.** The decodes of 25 D-STAR and 8 P25 tone vectors change, and no rate-33 one. The five tone
measures the gates check (the frequency, level and balance errors, the share of windows playing a
wrong tone and the detector's agreement) do not move on any vector in any mode; the scoreboard's
`tone.alignment_corr`, which no gate checks, moves by at most +0.0012 (D-STAR `cpvbad`) and in P25
by -0.00008 to +0.00016. Applying the rule after a P25 tone model (the second fix) changes three P25
tone decodes, `dtmf60` (19 frames), `dtmf180` (8) and `dtmfvbad` (12); over those frames the residual
against DVSI's decode, over the signal, is 0.0328, 6.0324 and 0.0384 with this frame's voicing there
and 0.0321, 4.3615 and 0.0379 with the rule (the scratch build), their whole-file correlations equal
to four decimals; their `tone.alignment_corr` against the head moves from +0.0000061, 0 and -0.0000650
in the first run to -0.0000188, +0.0000013 and -0.0000754. Whole-file correlation with DVSI's decode
(the scratch build): `cpvbad` P25 0.9465 → 0.9615 and D-STAR 0.9494 → 0.9750, D-STAR `dtmfvbad`
0.9602 → 0.9610, no other vector moving in the fourth decimal; P25 `alltone` stays at 0.8253 (12 of
its frames differ).

**Channel errors.** No mute decision changes at any error rate, hard or soft, frame for frame: the
mute flag of every frame of the 24 channel-error decodes (8 vectors in each of the three codecs,
30,600 frames, 1,678 of them muted) is the same in both builds. Every D-STAR and P25 channel-error
decode changes, and no rate-33 one. The scoreboard's log-spectral distance from DVSI's decode
(`lsd_common_db`: over the frames loud in DVSI's decode of the clean bits, not muted by DVSI,
measured in both runs, with no muted frame at or next to them; the set asks for a coverage of at
least 0.1 of the loud frames DVSI plays) falls, means over the error rates: D-STAR by 0.2051 dB with
hard decisions (0.140 to 0.272 dB per vector, a coverage of 1.0 on every vector) and 0.2294 dB with
soft ones (0.117 to 0.319 dB; a coverage of 0.937 at the least, on `dam_e10_sd`), P25 by 0.0501 dB
(0.032 to 0.084 dB; at the least 0.974, on `dam_e10_hd`) and 0.0401 dB (0.013 to 0.085 dB; at the
least 0.971, on `dam_e10_sd`). The waveform correlation, which no channel check uses, over whole
files, head → rule (the library's decodes): with hard decisions P25 `dam_e1_hd` 0.818979 → 0.804567,
`dam_e2_hd` 0.468670 → 0.453890, `dam_e5_hd` 0.034481 → 0.033988 and `dam_e10_hd` 0.010412 → 0.009247
(it falls at all four rates), D-STAR 0.916266 → 0.936127, 0.898719 → 0.917346, 0.715943 → 0.737629
and 0.330184 → 0.329617 (it rises at three and falls on `dam_e10_hd`); with soft decisions P25
`dam_e1_sd` 0.665945 → 0.666056, `dam_e2_sd` 0.458545 → 0.458726, `dam_e5_sd` 0.103938 → 0.104091 and
`dam_e10_sd` 0.052836 → 0.054864 (it rises at all four rates), D-STAR 0.973799 → 0.995100, 0.972339 →
0.994828, 0.706225 → 0.717810 and 0.085516 → 0.082729 (it rises at three and falls on
`dam_e10_sd`). With hard decisions, over all the frames whose decode changes, the head's
residual against DVSI's decode, over the signal, is 0.1942, 1.1237, 2.2521 and 3.1105 at the four
rates for P25 (D-STAR 0.3170, 0.3516, 0.5891 and 1.2930). A magnitude-spectrum distance (the rms
over 0.1 to 3.8 kHz of the dB difference of 20 ms spectra smoothed over 125 Hz), over the loud frames
that change (DVSI's rms 100 or more), falls on average at every rate: P25 `dam_e1_hd` 3.040 → 2.603
dB, `dam_e2_hd` 4.697 → 4.447, `dam_e5_hd` 7.660 → 7.523, `dam_e10_hd` 10.296 → 10.231 (1,324 frames,
nearer in 803, farther in 521), D-STAR 1.960 → 1.368, 1.951 → 1.396, 2.150 → 1.637 and 3.241 →
2.924 (1,043 frames, 930 nearer, 113 farther); on the clean validation vectors the same measure gives
P25 2.380 → 1.858 dB (2,150 frames nearer / 428 farther) and D-STAR 2.257 → 1.594 (2,373 / 279).
That distance does not see the phase: it does not show why P25's waveform correlation falls with
hard decisions (and rises with soft ones), nor which band's voicing is right, and its frame
selection is not that of the scoreboard's `lsd_common_db`. The cause of P25's fall with hard
decisions was not established. D-STAR's DTX vectors, which the
scoreboard does not score as speech: waveform correlation `clean_dtx` 0.8581 → 0.8482, `noisy_dtx`
0.9535 → 0.9328, `dam_dtx` 0.9580 → 0.9584, the same magnitude-spectrum distance 2.213 → 1.590,
1.996 → 1.432 and 1.926 → 1.375 dB (1,442 frames nearer, 167 farther).

**On hardware**, the DVstick 30's steps of
[#157](https://github.com/arancormonk/mbelib-neo/issues/157): the mean error of the 18 single-band
steps' curves falls from 0.346 to 0.189 (two of the 18 get worse, one at position 34 of 1101 → 1111
by 0.426 → 0.478), the onset steps' mean moves in the fifth decimal (0.034 in both to three
decimals), and the offset, gain, pitch and double steps' means are equal to nine decimals (the
table to five decimals, and the six means of three positions each, are in [Synthesis in
subframes](#synthesis-in-subframes)); the curves are full-band rms in 16-sample blocks and show the
steps as a whole, not which band starts when. The library's decodes of these streams are the same after the review's fixes.

**Other modes and the encoders.** Every rate-33 bit file of DVSI's vectors at the top of its mode's
directory (106 files, `irstia` and `mark` left out) decodes identically in the two builds; of the
D-STAR and P25 files, 78 of 106 and 56 of 103 decode otherwise (55 before the second fix). All 285
encodings (95 input files,
`irstia` and `mark` not among them, three encoders) are identical between the two builds. These
were compared in the scoreboard's builds (Release, SIMD, IEEE arithmetic) only.

**Cost.** The median complete hard-decision frame decode moves by +2.0% (D-STAR) and +0.9% (IMBE),
and by -1.6% for AMBE+2, which the change does not reach (a control), `c3c6d42` against `e9cf77c`,
with the ranges of the invocations overlapping in every codec, measured on a machine busy with
other work (one-minute load 1.7 to 2.1; [benchmarks](testing.md#mid-frame-voicing-by-columns)): the
cost, if any, is within what this measurement resolves, the D-STAR change the nearest to it. The
review's fixes (both frames' nearest harmonics in double for every band the columns evaluate,
before the two states are compared, and one test fewer per frame) were
not timed.

**Gate sets**, declared in the branch's first commit, `3593e8c`, before any build of the branch was
scored (the study's scratch build had been measured on the development and validation vectors
before that, above), and evaluated in one combined run over the development, validation, tone and
channel partitions, `c3c6d42` against `e9cf77c` and again, after the review's fixes, `c3c6d42`
against `e86afd0`, with the same outcome; the observed values are the second run's:

| Gate set | Result | What it establishes |
|---|---|---|
| `mid-columns-development` | passes all 18 checks | on the vectors the rule was found on, within the limits declared: D-STAR's and P25's correlation means rising by at least 0.0005 and 0.001 with intervals above 0 (observed +0.0010 and +0.0018) and no vector falling by more than 0.0005 (observed: none falls), the distance from DVSI's decode within its limits on the mean and on every vector, and the envelope correlation, the distance from the input and the level within theirs on the mean only (one vector may pass a mean limit and the set still pass), each check measured on all 20 vectors (`every_vector`), rate 33's speech decodes and every encoder's speech rows unchanged |
| `mid-columns-validation` | passes all 18 checks | the same on the validation vectors, the mean limits again on the mean only (limits 0.003 and 0.0005; observed D-STAR +0.0115, P25 +0.0012, none falling), also read in the study |
| `mid-columns-tones` | passes all 18 checks | the D-STAR and P25 tone measures within their limits on the vectors both runs measure, rate 33's speech and tone decodes unchanged, the D-STAR and rate-33 encoders' per-vector tone statistics (agreement, extra tones, level) unchanged; not that the same tone frames or bits were sent |
| `mid-columns-channel` | passes all 33 checks | the per-vector mute statistics unchanged in every mode (not by itself that no frame's mute decision changed), D-STAR's and P25's mean distance from DVSI's decode on the loud frames not rising past 0.01 dB (observed: falling), rate 33's decodes and every encoder's speech rows unchanged |
| `mid-columns-final` | passes all 18 checks (scored once, on 2026-11-26, after the review; below) | on `irstia` and `mark`, held out: D-STAR's and P25's correlation means rising (+0.0026 and +0.0032; no vector falling), the distance from DVSI's decode falling on the mean and on both vectors, the envelope correlation, the distance from the input and the level within their limits on the mean only, rate 33's decodes and every encoder's bits unchanged, each gated measure measured on both vectors in both runs |

Notes on the declarations. A check that the mean rises by at least a value is the scoreboard's
`increase` rule, which also requires the 95% interval above 0. Every paired check of the
development and validation sets requires every vector (`every_vector`), and the baseline measured
each of them on all 20 and all 5 vectors; the tones set requires every vector for the detector
agreement and the encoders' extra-tone rates, which cover all 52, and its description records the
coverage of the others from the baseline's own scoreboard (D-STAR's frequency, level and
wrong-tone figures on 44 of the 52 vectors, its balance on 43, its encoder agreement and volume
error on 35; P25's on 49 and 48; the rate-33 encoder's agreement and level error on 51, `cpvbad`
absent). A check on a metric that applies to part of a partition cannot
require every vector, so a vector on which the metric went missing in both runs would pass
unnoticed; in this run each such check covered the vectors listed (and the channel checks the four
hard- and the four soft-decision vectors each). `every_vector` makes a check need a measurement on
every vector; it does not make a limit on the mean a limit on each vector. The tones set's description says that neither run
can measure a figure where the baseline lacks it: that is not so, the candidate can measure it
alone, and the comparator then fails the check as unpaired (no check here was unpaired). The identity
counters (`speech.pcm_changed_vectors`, the encoder-row counter) count the speech vectors present in
the run, so they were evaluated on the combined run. Bit identity of the encodings is the separate
comparison of the 285 encodings above, and frame-for-frame mute identity the separate comparison of
the 30,600 channel-error frames: the gates establish neither. The channel set gates no waveform
correlation; its description gives as the reason the residual "0.18 to 0.26" in P25's changed frames
and runs "out of phase ... whichever voicing is right": the figure was two subgroups of `dam_e1_hd`
alone, its frames voiced after a voiced frame whose rendered voicing equals the decoded one, 136
with corrected errors (0.176) and 68 without (0.264); the explanation was not established, and
P25's correlation falls with hard decisions only, rising with soft ones (**Channel errors**,
above). The tones, channel and
final sets' descriptions do not give the order of events: the study's scratch build was measured on
the development and validation vectors before the sets were declared, and the baseline's own
scoreboard was read for coverage. The final set's description, "scored once, after the review", states the
plan, which was kept: it was scored once, after the review and its two rounds of fixes (below). All five sets were declared for the rule as first written, which took
this frame's voicing after a P25 tone model, and their descriptions say so ("nor in P25's
rendering after a tone model"); the reviewed change applies the rule after a tone model as after any, three P25 tone
decodes differ by it (`dtmf60` 19 frames, `dtmf180` 8, `dtmfvbad` 12), and the scores reported
here are the final code's. (A tone model itself still keeps its own model at mid-frame; the frame
after it follows the rule like any frame.)

**The final partition.** The gate set `mid-columns-final` was declared with the others in `3593e8c`,
before any build of the branch was scored: on `irstia` and `mark`, rate 33's decodes of DVSI's bits
unchanged, every encoder's bits unchanged, and for D-STAR and P25 the waveform correlation with
DVSI's decode falling by no more than 0.0005 on the mean and 0.002 on any vector, the distance from
DVSI's decode rising by no more than 0.01 dB on the mean and 0.03 dB on any vector, and, on the
mean only, the envelope correlation falling by no more than 0.001, the distance from the input
rising by no more than 0.01 dB and the level against DVSI's moving by at most 0.05 dB. The final
partition was scored once for this change set, on 2026-11-26, after the review and its two rounds of
fixes: `c3c6d42` against the head's code (the build of `e86afd0`, whose code the head shares apart
from comments, as candidate and decode tools), on `irstia` and `mark`, both builds aligned freely
first and then scored at their joint lags, both runs' encodings decoded by the candidate's decoder.
It passes all 18 checks of `mid-columns-final` (none fails). The library's decode of DVSI's bits,
partition means over the two vectors and the waveform correlation of each:

| Codec | Waveform corr. vs DVSI, sample for sample, mean | `irstia` | `mark` | LSD vs DVSI (dB), mean | Envelope corr., mean | Level vs DVSI (dB), mean | LSD vs input (dB), mean |
|---|---|---|---|---|---|---|---|
| D-STAR | 0.9830 → 0.9856 | 0.9774 → 0.9806 | 0.9886 → 0.9906 | 1.7794 → 1.5116 | 0.9766 → 0.9774 | +0.0592 → +0.0608 | 5.5432 → 5.5380 |
| P25 | 0.9837 → 0.9869 | 0.9795 → 0.9840 | 0.9880 → 0.9898 | 1.5248 → 1.3346 | 0.9843 → 0.9846 | -0.0349 → -0.0328 | 5.3690 → 5.3658 |
| Rate 33 | 0.9860 → 0.9860 | 0.9828 → 0.9828 | 0.9891 → 0.9891 | 1.5074 → 1.5074 | 0.9783 → 0.9783 | -0.0731 → -0.0731 | 5.5216 → 5.5216 |

The gate's paired deltas: the waveform correlation, D-STAR +0.002598 on the mean (`irstia`
+0.003181, `mark` +0.002015), P25 +0.003178 (`irstia` +0.004545, `mark` +0.001810); the distance
from DVSI's decode, D-STAR -0.2678 dB (`irstia` 1.9671 → 1.7470, `mark` 1.5917 → 1.2762 dB), P25
-0.1901 dB (`irstia` 1.4894 → 1.2989, `mark` 1.5601 → 1.3704 dB); the envelope correlation +0.00083
and +0.00028; the distance from the input, D-STAR -0.0052 dB (`irstia` 5.5420 → 5.5401, `mark`
5.5444 → 5.5359 dB), P25 -0.0032 dB (`irstia` 5.3236 → 5.3243, `mark` 5.4145 → 5.4074 dB); the level
against DVSI's +0.0017 dB and +0.0021 dB; rate 33's decodes unchanged (no vector with changed PCM),
and every encoder's bits unchanged (no vector with changed rows in any mode). Both vectors were
measured in both runs for every gated measure, as the checks require (`every_vector`), so none
failed for a missing measurement.

This is the change's out-of-sample check: the rule was found on the labels of the development and
validation vectors together, and its scratch build was measured on both, so neither of those
partitions is out of sample for it. On the final partition the correlation gains, D-STAR +0.0026
and P25 +0.0032, are larger than development's (+0.0010 and +0.0018) and P25's on validation
(+0.0012), and smaller than D-STAR's on validation (+0.0115, most of it one run of `t03`). They rest
on two vectors, so the partition means say little beyond the sign: the correlation rises and the
distance from DVSI's decode falls on both vectors of both codecs. The distance from the input rises
on one vector of one codec (P25 `irstia`, by 0.0007 dB) and falls on the other three. The final
partition shows nothing about the channel-error or tone behaviour, nor which band's voicing is
right in a given frame.

### November 2026: the pulse train through a row-17 frame

The mid-frame columns change's head (`e446950`, whose rate-33 decodes are those of `c3c6d42`, PR
#199) against the row-17 train change ([Rate 33: rows 13, 15 and 17](#rate-33-rows-13-15-and-17),
**The pulse train through a row-17 frame**): in the rate-33 rendering the train a stream carries
into a frame is also taken from a frame marked `MBE_FRAME_MARK_MID_COLUMN_0`. The study that found
it read the development and validation vectors and DVSI's rate-33 tone vectors, never `irstia` or
`mark`. Its scratch build was measured on the development and validation vectors (the cases of the
section above, the whole-file figures and the accounting below) and on the rate-33 tone vectors
before the gate sets were declared, so the development and validation partitions are not out of
sample here; the final partition was held out and scored once, after the review (below). The
library's rate-33 decodes of the 25 development and validation vectors equal the scratch build's
with the change, file for file, so the scratch build's rate-33 figures below are the library's. The
scratch build's switch is not limited to the rate-33 rendering, as the library's condition is; no
D-STAR or P25 decoder writes the mark, so its D-STAR and P25 decodes cannot show that difference,
and D-STAR's and P25's identity below comes from the comparison of the library's own two builds. The
candidate is the branch's head at `2a3a64f`, whose library code is `8cbef45`'s (the commits after
it change tests and documents only). Both builds aligned freely first and were then scored at their
joint lags, both runs' encodings decoded by the candidate's decoder.

**Speech**, rate 33 decoding DVSI's bits, on the scoreboard: the waveform correlation with DVSI's
decode, sample for sample (`vs_dvsi.waveform_corr_unshifted`, over the frames where DVSI's decode
lies within 40 dB of its loudest), with the mean paired change and its 95% interval and the vectors
whose correlation rises and falls by more than 1e-6; the log-spectral distance from DVSI's decode
(`vs_dvsi.lsd_db`), the envelope correlation with it (`vs_dvsi.env_corr`), the level against it
(`vs_dvsi.level_offset_db`) and the log-spectral distance from the input (`ours.lsd_db`), partition
means:

| Partition | Waveform corr. vs DVSI | Vectors up / down | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|---|
| development | 0.989337 → 0.992449 (+0.003112, 95% interval 0.000123 to 0.011608) | 4 / 0 (`t02` +0.029258, `clean` +0.016964, `fambf22c` +0.014837, `ucarm20a` +0.001187) | 1.2299 → 1.2291 | 0.99207 → 0.99209 | +0.0139 → +0.0139 | 6.2692 → 6.2692 |
| validation | 0.992142 → 0.992142 (+3.7718e-09, 95% interval 0 to 8.0824e-09) | 0 / 0 (`tambf22b` +1.886e-08, under 1e-6) | 1.0642 → 1.0641 | 0.98840 → 0.98840 | -0.0403 → -0.0403 | 4.8303 → 4.8303 |

Per vector, the correlation of `ucarm15a` moves by +1.7e-08 and of validation's `tambf22b` by +1.9e-08,
and no other vector's moves. The distance from DVSI's decode falls on the five development vectors
that change (`fambf22c` by 0.0077 dB, `t02` 0.0033, `clean` 0.0027, `ucarm20a` 0.0022, `ucarm15a`
0.00002) and on `tambf22b` (0.0001 dB); the envelope correlation rises on those six; the distance from
the input rises on `fambf22c` (+0.00025 dB), `t02` (+0.00015) and `ucarm20a` (+0.00013) and falls on
`clean` (-0.00011), `ucarm15a` (-0.00003) and `tambf22b` (-0.000005); the level against DVSI's
moves by at most 0.00006 dB (`clean`). D-STAR's and P25's partition means and decodes are unchanged.
Validation's one changed decode is `tambf22b`'s: a row-17 frame at 1641 that follows a frame with no
voiced band (so outside the twelve cases' selection) and is followed by a voiced frame; frames 1642
and 1643 change, its correlation moves from 0.991451890 to 0.991451909 (+1.886e-08) and its distance
from DVSI's decode from 1.08087003 to 1.08075404 dB (-0.000116). Validation so shows no gain at the
precision reported (four decimals, and under the 1e-6 per vector counted above).

The decodes that change, and their frames that change (16-bit samples, the library's builds): `t02`
31 frames, from frame 1418, `fambf22c` 28 from 1213, `clean` 25 from 1417, `ucarm20a` 6 from 610,
`ucarm15a` 2 (610 and 611) and, on validation, `tambf22b` 2 (1642 and 1643); the other 19 rate-33
decodes are identical.

Over whole files, the measure of the study (the correlation of our decode with DVSI's, sample for
sample over the whole file at lag 0, and the residual energy of the difference over the signal's,
means over the vectors; the scratch builds): development 0.989326 → 0.992438 (+0.003112; 0.9893 →
0.9924), residual 0.0214 → 0.0152, four vectors rising by more than 1e-6 and none falling (`t02`
0.9660 → 0.9953, +0.029257; `clean` 0.9723 → 0.9893, +0.016962; `fambf22c` 0.9774 → 0.9922,
+0.014836; `ucarm20a` 0.9891 → 0.9903, +0.001187; `ucarm15a` moving by less than 1e-7); validation
0.992119 → 0.992119, `tambf22b` moving by less than 1e-7. A scratch variant that never drops a
carried train gives the same 25 rate-33 decodes ([Rate 33: rows 13, 15 and
17](#rate-33-rows-13-15-and-17)).

**Where the residual sits**, in the accounting of [November 2026: mid-frame
voicing](#november-2026-mid-frame-voicing) (every vector scaled to unit power, the residual and the
signal summed over the vectors; voiced runs with a whole-run correlation under 0.8 out of line), the
baseline and with the change, on the scratch builds:

| Rate 33 | Residual | Voiced runs in line | Voiced runs under a correlation of 0.8 | Frames with no voiced band |
|---|---|---|---|---|
| development | 0.0248 → 0.0163 | 0.0138 → 0.0138 | 0.0085 → 0.0000 | 0.0025 → 0.0025 |
| validation | 0.0163 → 0.0163 | 0.0124 → 0.0124 | 0.0000 → 0.0000 | 0.0039 → 0.0039 |

The runs are grouped by the voicing of the decode the accounting has always used, of the 2.6.12
era, which still voiced row 17's own model, so that a run spans its row-17 frame. On development, 23
of the 1,825 runs lie under 0.8 with the baseline, together 0.00846548 of the signal, and four of
them hold 0.00845526 of it: `t02` from frame 1597 (49 frames, correlation -0.112), 0.00386282;
`clean` from 1391 (50 frames, +0.398), 0.00285607; `fambf22c` from 1197 (29, +0.412), 0.00106158;
and `fambf22c` from 1468 (21, +0.774), 0.00067479; the other 19 hold 0.00001 together. With the
change 19 of the 1,825 runs lie under 0.8, together 0.00001023 (the table's 0.0000); the four runs
correlate +0.993, +0.988, +0.987 and +0.997 and hold 0.00002560, 0.00005467, 0.00002330 and
0.00000965, counted among the runs in line. At six decimals the development cells are 0.024787 →
0.016305 (residual), 0.013846 → 0.013822 (runs in line: the four runs come in with 0.000113, and
other runs that span a row-17 frame and stayed above 0.8 as wholes lose 0.000137), 0.008465 →
0.000010 (under 0.8) and 0.002475 → 0.002473 (no voiced band); the validation cells 0.016323,
0.012390, 0.0000014 (4 of 495 runs) and 0.003932, unchanged. Each frame's correlation with DVSI's
decode after the row-17 frame, baseline → change, frames f + 1 to f + 4 of the row-17 frame f:
`clean` 1416 +0.45071 → +0.97698, +0.14528 → +0.98479, -0.04680 → +0.99088, +0.11331 → +0.99547;
`fambf22c` 1212 -0.20528 → +0.99686, -0.06268 → +0.99559, -0.15887 → +0.98261, -0.50298 →
+0.98837; `t02` 1618 +0.41759 → +0.95308, -0.84452 → +0.99460, -0.82241 → +0.99516, -0.68819 →
+0.99351; `fambf22c` 1474 +0.99269 → +0.99636, +0.97436 → +0.99476, +0.75344 → +0.99291, +0.78892
→ +0.99551 (the baseline falls over three frames there); `ucarm20a` 609 +0.50423 → +0.91756,
+0.85541 → +0.99329, +0.47183 → +0.99463, +0.03832 → +0.97346. (The run accounting and these frame
correlations were measured by the study's scripts on its scratch builds, whose rate-33 decodes are
the library's.)

**Tones.** No rate-33 tone decode changes (52 of 52 by the scoreboard's hashes), and none of
D-STAR's or P25's; the five tone measures the gates check and the encoders' tone statistics do not
move on any vector.

**Channel errors.** No mute decision changes at any error rate, hard or soft, frame for frame: the
mute flag of every frame of the 24 channel-error decodes (8 vectors in each of the three codecs,
30,600 frames, 1,678 of them muted) is the same in both builds. Two rate-33 channel-error decodes
change, `dam_e10_hd` in 3 frames (387, 413 and 861) and `dam_e10_sd` in 10 ([Rate 33: rows 13, 15
and 17](#rate-33-rows-13-15-and-17), "Where else the change acts"). The scoreboard's log-spectral
distance from DVSI's decode, `lsd_common_db`, is taken over the frames loud in DVSI's decode of the
clean bits that DVSI does not mute, with no muted frame at or next to them (neither a loud frame
DVSI mutes nor a frame either run muted, mapped at that run's lag) and with a spectrum measured in
both runs (`channel_lsd()`); its coverage is their share of the loud frames DVSI leaves unmuted
(rate 33 `dam_e10_hd`: 273 of 364, 0.750; with soft decisions 0.991 at the least). It does not move
on any of the eight vectors (a paired change of 0 on each). No D-STAR or P25 channel-error decode
changes.

**Other modes and the encoders.** Every D-STAR and P25 bit file of DVSI's vectors at the top of each
mode's directory, `irstia` and `mark` left out (106 D-STAR and 103 P25 files, the tone,
channel-error, DTX and synthetic ones among them, the soft-decision ones decoded as soft), decodes
identically in the two builds, and so do 97 of the 106 rate-33 files: `clean`, `clean_dtx` (from
frame 1417), `dam_e10_hd`, `dam_e10_sd`, `fambf22c`, `t02`, `tambf22b`, `ucarm15a` and `ucarm20a`
decode otherwise. All 285 encodings (95 input files, `irstia` and `mark` not among them, three
encoders) are identical between the two builds. These were compared in the scoreboard's builds
(Release, SIMD, IEEE arithmetic) only.

**Cost.** Not measured: the change adds one comparison of the previous frame's mark per frame, in
the rate-33 rendering, and only after a frame with no voiced band.

**Gate sets**, declared in the branch's first commit, `6a32fb8`, before any build of the branch was
scored (the study's scratch build had been measured on the development, validation and rate-33 tone
vectors before that, above, and the baseline's own scoreboard of the four partitions, `e446950`'s
build decoding its own encodings, was read to record each measure's coverage), and evaluated in one
combined run over the development, validation, tone and channel partitions, `e446950` against the
head's code:

| Gate set | Result | What it establishes |
|---|---|---|
| `r33-row17-train-development` | passes all 12 checks | on the vectors the change was found on, within the limits declared: rate 33's correlation mean rising by at least 0.002 with its interval above 0 (observed +0.003112, interval 0.000123 to 0.011608) and no vector falling by more than 0.0005 (observed: none falls), the distance from DVSI's decode within its limits on the mean and on every vector, and the envelope correlation, the distance from the input and the level within theirs on the mean only, each check measured on all 20 vectors (`every_vector`); D-STAR's and P25's speech decodes and every encoder's speech rows unchanged |
| `r33-row17-train-validation` | passes all 12 checks | the same limits for harm on the validation vectors, also read in the study (the correlation's mean falling by no more than 0.0005, observed +3.7718e-09); one decode changes, in two frames, and the partition shows no gain at the precision reported |
| `r33-row17-train-tones` | passes all 15 checks | rate 33's tone measures within their limits on the vectors both runs measure, D-STAR's and P25's speech and tone decodes unchanged, the D-STAR and rate-33 encoders' per-vector tone statistics (agreement, extra tones, level) unchanged; not that the same tone frames or bits were sent |
| `r33-row17-train-channel` | passes all 33 checks | the per-vector mute statistics unchanged in every mode (not by itself that no frame's mute decision changed), rate 33's mean distance from DVSI's decode on the loud frames not rising past 0.01 dB (observed: no change), D-STAR's and P25's decodes and every encoder's speech rows unchanged |
| `r33-row17-train-final` | passes all 12 checks (scored once, on 2026-11-26, after the review and its fixes; below) | on `irstia` and `mark`, held out: no harm (no decode of either vector changes, in any mode, and every encoder's rows are unchanged); it cannot show the gain, since nothing in it changes |

Notes on the declarations. A check that the mean rises by at least a value is the scoreboard's
`increase` rule, which also requires the 95% interval above 0. Every paired check of the
development, validation and final sets requires every vector (`every_vector`), and the baseline
measured each of the development and validation ones on all 20 and all 5 vectors; `every_vector`
makes a check need a measurement on every vector, it does not make a limit on the mean a limit on
each vector. The tones set requires every vector for rate 33's detector agreement and the D-STAR and
rate-33 encoders' extra-tone rates, which cover all 52, and its description records the coverage of
the others from the baseline's own scoreboard (rate 33's frequency, level and wrong-tone figures on
49 of the 52 vectors, its balance on 48; D-STAR's encoder agreement and volume error on 35; the
rate-33 encoder's agreement and level error on 51, `cpvbad` absent); a check on a metric that
applies to part of a partition cannot require every vector, so a vector on which the metric went
missing in both runs would pass it unnoticed, and in this run each such check covered the vectors
listed (the channel checks the four hard- and the four soft-decision vectors each). The identity
counters (`speech.pcm_changed_vectors`, the encoder-row counter) count the speech vectors present in
the run, so they were evaluated on the combined run. Bit identity of the encodings is the separate
comparison of the 285 encodings above, and frame-for-frame mute identity the separate comparison of
the 30,600 channel-error frames: the gates establish neither. No channel-error vector is gated on
rate 33's waveform. The sets' descriptions give the order of events (the scratch build measured
before the declaration) and the coverage read from the baseline's scoreboard. Where the
descriptions, which stay as declared, say more than the checks or the measurements:

- All five open with a row-17 frame "which voices column 0 at mid-frame and no band at its end" and
  a frame after it that "started a new train at its first sample": true of a row-17 frame as
  decoded without channel errors. The two other cases of [Rate 33: rows 13, 15 and
  17](#rate-33-rows-13-15-and-17), "Where else the change acts", are not said: a row-17 frame whose
  end the adaptive smoothing voices, where nothing changes, and another row-17 frame after it, where
  the change acts too.
- The channel set's description says the distance is measured "on the loud frames neither run nor
  DVSI mutes"; what is checked is narrower (`lsd_common_db`, above: no muted frame at or next to
  them, each run's mutes mapped at its lag, a spectrum measured in both runs). The figure is formed
  at comparison time from both runs' frames: the baseline's own scoreboard, read before the
  declaration, holds `lsd_db` and `lsd_coverage`, not it. And the set does not require the four
  vectors per kind that the run measured: a figure missing in both runs on a vector would pass.
- The tones set's description says the encoder's level error is taken "where both send a tone"; it
  is taken over the frames where both send the same tone (`encoder_tone_agreement()`; rate-33
  `dtmfvbad`: 115 aligned pairs of frames where both send a tone, 114 of them the same tone). And
  its "say nothing of the voice frames in these vectors" is too strong: the detector's agreement
  counts the steady windows where both detect no tone, and the extra-tone rate counts the tone
  frames sent where DVSI sends voice.
- The encoder-row counter (`encoder_row_changes()`) skips a vector whose baseline holds no row hash
  (the PCM counters count a hash missing on either side as changed), so "no speech vector's
  encoder rows change" rests on every speech vector carrying a row hash in both runs. In this run
  all do: all 75 speech entries (25 per mode) carry `encoder.rows_sha256` in the baseline's and the
  candidate's reports, and all 255 entries a PCM hash in both. The final scoring is to be checked
  the same way, and the count recorded with its result.
- The validation set's description says the partition "can show no gain here": it shows none at
  the precision reported (above), with one decode changing in two frames.

**The final partition.** The gate set `r33-row17-train-final` was declared with the others in
`6a32fb8`, before any build of the branch was scored: on `irstia` and `mark`, D-STAR's and P25's
decodes of DVSI's bits unchanged, every encoder's bits unchanged, and for rate 33 the waveform
correlation with DVSI's decode falling by no more than 0.0005 on the mean and 0.002 on any vector,
the distance from DVSI's decode rising by no more than 0.01 dB on the mean and 0.03 dB on any
vector, and, on the mean only, the envelope correlation falling by no more than 0.001, the distance
from the input rising by no more than 0.01 dB and the level against DVSI's moving by at most 0.05
dB. Neither the change nor its scratch build was measured on these vectors. The final partition was
scored once for this change set, on 2026-11-26, after the review and its fixes: `e446950` (the
build of `git archive e446950`) against the head's code (the build of `3840d2b`, whose library code
is `8cbef45`'s; its decodes of the 25 development and validation vectors are those of the build
scored above, in all three modes, 25 of 25 each), on `irstia` and `mark`, both builds aligned freely
first and then scored at their joint lags, both runs' encodings decoded by the candidate's decoder.
It passes all 12 checks of `r33-row17-train-final` (none fails). No decode of the two vectors
changes: the PCM hash of each of the six entries (`irstia` and `mark` in each mode) is the same in
both runs, rate 33's included, and so is every encoder-row hash; the two vectors hold no frame the
change touches. As the notes above ask of the encoder-row counter, all six entries carry an
encoder-row hash and a PCM hash in both reports. The library's decode of DVSI's bits, the same in
both runs:

| Codec | Waveform corr. vs DVSI, sample for sample, mean | `irstia` | `mark` | LSD vs DVSI (dB), mean | Envelope corr., mean | Level vs DVSI (dB), mean | LSD vs input (dB), mean |
|---|---|---|---|---|---|---|---|
| Rate 33 | 0.985956 (0.9860) | 0.982792 | 0.989120 | 1.5074 (`irstia` 1.6848, `mark` 1.3301) | 0.9783 | -0.0731 | 5.5216 |
| P25 | 0.986919 (0.9869) | 0.984030 | 0.989808 | 1.3346 | 0.9846 | -0.0328 | 5.3658 |
| D-STAR | 0.985597 (0.9856) | 0.980559 | 0.990634 | 1.5116 | 0.9774 | +0.0608 | 5.5380 |

Every paired delta of the gate is 0. What this establishes: the held-out partition shows no harm,
and it cannot show the gain, since nothing in it changes. The change's evidence is the development
vectors', in sample (the twelve cases, the phase measurements, the one row-17 pair without channel
errors), with validation moving by +3.8e-09.

### November 2026: the pulsed state

The row-17 train change's head (`74126a3`, whose library code is `8cbef45`'s) against the pulsed
state ([D-STAR and P25: the pulsed state](#d-star-and-p25-the-pulsed-state)): a P25 or D-STAR frame
with no voiced band at one of the 56 reserved pitch codes plays the pattern of pulsed 500 Hz
channels its code names. The study that found the rule read the development and validation vectors,
DVSI's tone, hard-decision channel-error and D-STAR DTX vectors and the DVstick 30 package of
[#157](https://github.com/arancormonk/mbelib-neo/issues/157), never `irstia` or `mark`. Its scratch
build (outside the repository) was measured on those before the gate sets were declared, so the
development and validation partitions are not out of sample here; the final partition was held out
and scored once, after the two reviews and their fixes ([the final partition](#the-final-partition),
below). The candidate is the branch's head at `d3fb70f`, whose
library code is `834d001`'s (the commit after it changes tests only). Both builds aligned freely
first and were then scored at their joint lags, both runs' encodings decoded by the candidate's
decoder (Release, SIMD, IEEE arithmetic).

**The library against the study's scratch build**, 16-bit samples, file for file. D-STAR: the
library's decodes of the 25 development and validation vectors equal the scratch build's with the
rule as written (its tag `pzm7`) on all 25. P25: they equal them on every frame of all 25 but the
first, where the library leaves 79 samples (at most 11 units) as the head plays them: the initial
state (`mbe_initMbeParms()`, all unvoiced) has eq 46's nominal pitch of b0 134, which is not the
pitch the P25 decoder stores for b0 134, so by clause 1 it has no code, where the scratch build's
test took the nearest code to any pitch (rounding, with no comparison of the bits) and gave it code
6 (its own subframe pulsed in channel 3), a pulse at sample 24 of every P25 decode's first frame.
Rate 33: the library's decodes equal the head's (`k1`) on all 25. On the other files the scratch build
decoded (115 tone, hard-decision channel-error and D-STAR DTX files, under its tag `pzx3`, which
applies clause 5 to every pair of frames with no voiced band) the library's decodes differ from
the scratch build's in 48 of the 59 D-STAR files and all 56 P25 files, in frame 0 and in some
elsewhere; a scratch copy of the library with the study's two tests switched
in (the nearest code, and clause 5 on every pair) decodes 162 of the 165 files that tag holds (the
50 speech files among them) exactly as the scratch build, and the other three, P25 `dtmf60`,
`dtmf180` and `dtmfvbad`, in all but 19, 8 and 12 frames: the frames after a P25 tone model where
the columns rule applies since the review of the mid-frame columns change (`e86afd0`), which the
scratch build, made at `47d8cf3` with switches, leaves out for a tone model. With the nearest-code
switch alone 54 of the 84 D-STAR files that tag holds still differ, and with clause 5 on every
pair alone all 84 match: D-STAR's differences from `pzx3` are all clause 5's extent. In P25, with
the nearest code alone 7 of its 81 files still differ (the three tone files above and the four
hard-decision channel-error files, which need both switches), and with clause 5 on every pair
alone all 81 differ (in frame 0 at least): P25's differences are mostly the nearest code's (in
frame 0, the initial state's; where else it acts in the channel-error files was not traced). No difference is left
unexplained; none is a difference in the rule as written.

**Speech**, decoding DVSI's bits, on the scoreboard: the waveform correlation with DVSI's decode,
the two waveforms correlated sample against sample at lag 0 (`vs_dvsi.waveform_corr_unshifted`, over the frames where DVSI's decode lies within
40 dB of its loudest), with the mean paired change and its 95% interval and the vectors whose
correlation rises and falls by more than 1e-6; the log-spectral distance from DVSI's decode
(`vs_dvsi.lsd_db`), the envelope correlation with it (`vs_dvsi.env_corr`), the level against it
(`vs_dvsi.level_offset_db`) and the log-spectral distance from the input (`ours.lsd_db`), partition
means, the library's builds:

| Codec | Partition | Waveform corr. vs DVSI | Vectors up / down | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|---|---|
| D-STAR | development | 0.991645 → 0.992583 (+0.000938, 95% interval 0.000529 to 0.001492) | 20 / 0 (smallest `dam10` +0.000026, largest `p01mirs` +0.004036) | 1.1607 → 0.9915 | 0.98861 → 0.99727 | +0.0255 → +0.0233 | 6.4197 → 6.4243 |
| D-STAR | validation | 0.993452 → 0.994654 (+0.001201, interval 0.000621 to 0.002134) | 5 / 0 (smallest `t03` +0.000373, largest `t04` +0.002722) | 1.0564 → 0.8656 | 0.98550 → 0.99410 | -0.0043 → -0.0099 | 4.8394 → 4.8628 |
| P25 | development | 0.986679 → 0.987520 (+0.000841, interval 0.000467 to 0.001283) | 20 / 0 (smallest `damcar` +0.000036, largest `p01mirs` +0.002896) | 1.4750 → 1.3145 | 0.98965 → 0.99478 | +0.0465 → +0.0431 | 6.3169 → 6.3290 |
| P25 | validation | 0.990761 → 0.991984 (+0.001223, interval 0.000644 to 0.002001) | 5 / 0 (smallest `t03` +0.000389, largest `t04` +0.002538) | 1.1252 → 0.9375 | 0.98768 → 0.99282 | -0.0125 → -0.0181 | 4.7209 → 4.7467 |

Per vector, the distance from DVSI's decode falls on all 25 in both codecs (D-STAR by 0.027 to 0.292
dB on development and 0.083 to 0.351 dB on validation; P25 by 0.031 to 0.323 dB and 0.082 to 0.304
dB), and the envelope correlation rises on all 25. **The distance from the input rises**: on 15 of
the 20 D-STAR development vectors (by at most 0.0254 dB, `tambf22a`; it falls by up to 0.0543 dB,
`dam80`) and on all 5 validation vectors (0.0114 to 0.0442 dB, `t04` the most), and on 17 of the 20
P25 development vectors (at most 0.0342 dB, `p01mirs`; down by up to 0.0532 dB, `dam80`) and all 5
validation vectors (0.0145 to 0.0500 dB, `t04`). The level against DVSI's moves by at most 0.0081 dB
on a D-STAR development vector and 0.0173 dB on a validation vector (P25 0.0150 and 0.0155 dB). Rate
33's partition means and decodes are unchanged, and so is every encoder's bits.

Over whole files, the measure of the study (the correlation of our decode with DVSI's, sample against
sample over the whole file at lag 0, and the residual energy of the difference over the signal's,
means over the vectors), on the library's decodes: P25 development 0.986625 → 0.987487 (+0.000862),
residual 0.02733 → 0.02558, validation 0.990749 → 0.991983 (+0.001234), residual 0.01848 → 0.01602;
D-STAR development 0.991634 → 0.992583 (+0.000949), 0.01717 → 0.01527, validation 0.993435 →
0.994653 (+0.001218), 0.01315 → 0.01072; every one of the 25 vectors rises in both codecs (P25
+0.000036 `dam10` to +0.002910 `p01mirs` on development; D-STAR +0.000026 `dam10` to +0.004056
`p01mirs`). These are the study's figures (`pzx3` against `k1`) to six decimals: the library's
decodes differ from `pzx3`'s in frame 0 of every P25 file and in six D-STAR files by one unit in one
to three samples (from `pzm7`'s in P25's frame 0 alone), and none of that moves them there.

**By class of frame** (the study's accounting, its scratch build's `pzx3` against `k1`; for each of
the 22 vectors whose decode by DVSI has an rms of 300 or more the class's residual energy over that vector's whole signal
energy, then the mean of those ratios over the 22, each vector weighted alike; classes by the head's
decode, frame 0 of every file counted with all other frames): P25 frames with no voiced band and no pulsed model in reach 0.00024 → 0.00020; with one
(the frame's own code, or the frame before's own-subframe pattern) 0.00203 → 0.00033; voiced frames
right after a frame whose own subframe is pulsed 0.00038 → 0.00008; all other frames 0.02037 →
0.02037; total 0.02302 → 0.02099. D-STAR 0.00019 → 0.00012; 0.00204 → 0.00033; 0.00050 → 0.00014;
0.01416 → 0.01416; total 0.01689 → 0.01475 (the energies pooled over the 22 instead, P25 0.02143 →
0.01898 and D-STAR 0.01579 → 0.01357). **By code** (the study's, `pzx3` against `k1`: frames with no voiced
band between two such frames, 23 vectors without `dam40` and `dam80`, residual over the frames' own
signal per code): the median over the 56 codes 0.277 → 0.019 (P25) and 0.290 → 0.015 (D-STAR);
every code improves but P25 b0 131 (0.101 → 0.106, 129 frames); codes left above 0.10: P25 181
(0.132, 26 frames) and 131, D-STAR 118 (0.128, 108 frames).

**Against the input and by ear-model** (the cost of this change). PESQ-NB (pesq 0.0.4) on the 25
vectors, the library's decodes (`74126a3` against `d3fb70f`; the review's fixes change no decode of
them; mean paired change ± its standard error, the vectors clustered by recording), against DVSI's
decode: P25 4.306 → 4.367 (+0.0605 ± 0.0056, 24 of 25 rise), D-STAR 4.323 →
4.407 (+0.0838 ± 0.0093, 25 of 25); against the **input**: P25 3.040 → 3.017 (-0.0230 ± 0.0039, 2 of
25 rise), D-STAR 2.952 → 2.932 (-0.0198 ± 0.0027, 4 of 25), where DVSI's own decode scores 3.026
(P25) and 2.933 (D-STAR). The library's score against the input ends within 0.01 of DVSI's own
decode's, below the score its noise rendering had; no equivalence was tested. The library's build
changes 12,261 P25 and 12,511 D-STAR frames of the 55,296 per codec against the head (the scratch
`pzx3`, which applies clause 5 to every pair of frames with no voiced band and takes the nearest
code by rounding: 12,286 and 12,517). In `mbe_quality_eval`'s
measures as the study took them, its scratch build `pzx3` against the head's decodes `k1` (each
file aligned freely; development / validation changes of the means; the scoreboard's on the
library's builds are in the table above): the
log-spectral distance from DVSI's decode P25 -0.1606 / -0.1878 dB (all 25 fall), D-STAR -0.1693 /
-0.1908 (all 25 fall); the envelope correlation with it +0.0051 / +0.0051 (P25), +0.0087 / +0.0086
(D-STAR); the level against it -0.0034 / -0.0056 dB (P25), -0.0023 / -0.0056 (D-STAR); the
log-spectral distance from the input P25 5.9977 → 6.0092 dB (+0.0080 / +0.0258; 22 of 25 rise),
D-STAR 6.1068 → 6.1091 (-0.0007 / +0.0140; 17 rise); the envelope correlation with the input
-0.0023 / -0.0119 (P25), +0.0000 / -0.0114 (D-STAR). The project's rule is to match DVSI where
DVSI's output and the input disagree; these are the figures where they do.

**Tones.** By the scoreboard's hashes, 15 of the 52 D-STAR and 14 of the 52 P25 tone decodes change,
and no rate-33 one; they change in 2 to 796 frames per file in D-STAR (`dtmfvbad` 796, `dtmf60`
730, `dtmf180` 729, `dtmf` 662, `alltone` 418, `cpvbad` 150, `alert` 35, eight files 2) and 1 to 772
in P25 (`dtmfvbad` 772, `dtmf180` 451, `dtmf60` 449, `dtmf` 414, `alltone` 299, `cpvbad` 98, `alert`
35, five files 2, `cp0` and `cp1` 1; 16-bit samples, the library's builds). The five tone measures the gates check (the frequency, level and balance
errors, the share of windows playing a wrong tone and the detector's agreement) do not move on any
vector in any mode (a paired change of 0 on each), nor do the encoders' tone statistics; the
scoreboard's `tone.alignment_corr`, which no gate checks, moves on 29 vectors, most up (D-STAR
`alltone` 0.954084 → 0.998272, P25 `alltone` 0.968384 → 0.987903; 20 up, 9 down by more than
1e-6), the largest fall -0.000993 (P25 `dtmf`). Of DVSI's 32 KNOX and dual-tone vectors (`knox_*`, `dtone_*`) none changes in either
codec: the study's single-frame changes there (frame 0 of `knox_13` to `knox_16`, among others)
are its tests' (above), not the rule's.

**Channel errors.** No mute decision changes at any error rate, hard or soft, frame for frame: the
mute flag of every frame of the 24 channel-error decodes (8 vectors in each of the three codecs,
30,600 frames, 1,678 of them muted) is the same in both builds. Every D-STAR and P25 channel-error
decode changes, and no rate-33 one. The scoreboard's log-spectral distance from DVSI's decode,
`lsd_common_db` (`channel_lsd()`: over the frames loud in DVSI's decode of the clean bits that DVSI
does not mute, with a spectrum measured in both runs, and with no frame at or next to them that
either run mutes (each run's mutes mapped at its lag) or that is loud and muted by DVSI; a quiet
frame DVSI mutes next to them does not exclude them, and the review found seven selected D-STAR and
P25 frames with one, zero-based frame 625 of D-STAR `dam_e5_hd` among them; coverage their share of
the loud frames DVSI leaves unmuted), falls on all 16 D-STAR and
P25 vectors: D-STAR by 0.0658, 0.0549, 0.0461 and 0.0223 dB with hard decisions (`dam_e1_hd` to
`dam_e10_hd`) and 0.0658, 0.0583, 0.0411 and 0.0095 dB with soft ones, P25 by 0.0479, 0.0262, 0.0160
and 0.0011 dB and 0.0434, 0.0231, 0.0092 and 0.0033 dB; the coverage is 0.937 at the least
(D-STAR `dam_e10_sd`). On the frames whose 16-bit samples differ between the library's builds
(`74126a3` against the final code, whose decodes equal `d3fb70f`'s file for file), the residual
energy against DVSI's decode over DVSI's energy there, at lag 0, `dam_e1` to `dam_e10`: D-STAR hard
decisions 0.0876 → 0.0212, 0.0888 → 0.0241, 0.1519 → 0.0515 and 2.1739 → 2.1727 (249, 237, 217 and
20 of 1,275 frames differ), soft decisions 0.0876 → 0.0212, 0.0870 → 0.0226, 0.1274 → 0.0537 and
2.4670 → 2.3639 (249, 240, 219 and 160 frames); P25 hard decisions 0.1916 → 0.1211, 0.6758 →
0.6325, 0.8990 → 0.8640 and 0.9921 → 0.9522 (262, 254, 221 and 95 frames), soft decisions 0.1645 →
0.0910, 0.5472 → 0.5260, 0.9513 → 0.9455 and 1.0911 → 1.0714 (268, 268, 249 and 134 frames).
(The figures the gate set's description gives, hard decisions only, are the study's scratch build
`pzx3` against `k1` on the frames that differ between those two.)

**Other modes and the encoders.** Of the bit files at the top of each mode's directory of DVSI's
vectors (`irstia` and `mark` left out: 106 D-STAR, 103 P25 and 106 rate-33 files, the tone,
channel-error, DTX and synthetic ones among them, the soft-decision ones decoded as soft), every
rate-33 file decodes identically in the two builds; 62 of the 106 D-STAR files and 58 of the 103
P25 files decode otherwise (the list, with the frames that change, is with the scoring run's
records, outside the repository). D-STAR's DTX vectors, which the scoreboard does not score as
speech: `clean_dtx` 595 of 3,700 frames change, `dam_dtx` 163 of 1,275, `noisy_dtx` 247 of 3,588;
the residual energy against DVSI's decode over DVSI's energy in the frames that change, the
library's builds as above: 0.4088 → 0.1417, 0.0981 → 0.0236 and 0.3683 → 0.1688. All 285 encodings (95 input files,
`irstia` and `mark` not among them, three encoders) are identical between the two builds. These
were compared in the scoreboard's builds (Release, SIMD, IEEE arithmetic). In each of the ten local
configurations (below), a differential harness outside the repository, built against the
baseline's library and the head's, gives the same output, byte for byte, for the pitch every D-STAR
and P25 pitch code decodes to and for rate 33's and ProVoice's data paths, `mbe_synthesizeSpeechf()`
and the three encoders over constructed streams whose caller states hold reserved pitches with no
voiced band, every frame mark, values in the state's slots and states without the subframe
rendering's mark (hashes of the output and of all three returned states).

**On hardware**, the DVstick 30's steps of [#157](https://github.com/arancormonk/mbelib-neo/issues/157)
(the study's scratch build with the head's switches, `hwk1b`, against the same build with `pzx3`'s,
`hwpz3`, decoding the package's streams): the three steps from frames with `b1` 0 at `b0` 100 to voiced
frames come nearer the unit's curves (mean error 0.06577 → 0.01992); of the three steps into such
frames one comes nearer (0.0692 → 0.0512) and two go further (0.1652 → 0.1891, 0.1599 → 0.1653; mean
0.13143 → 0.13520); every other row of the package is unchanged to the digits printed. The curves
are full-band levels in 16-sample blocks of one D-STAR code: they corroborate the fit at the steps
out of such frames and establish no clause of the rule ([D-STAR and P25: the pulsed
state](#d-star-and-p25-the-pulsed-state)).

**Cost.** The median complete hard-decision frame decode moves by +5.6% in D-STAR (6.454 → 6.816
us; the seven invocations' medians 6.445 to 6.482 and 6.791 to 6.880 us) and +5.6% in IMBE (6.864
→ 7.246 us; 6.850 to 6.959 and 7.211 to 7.307), and by +0.9% for AMBE+2 (7.032 → 7.092 us; 7.013 to
7.056 and 7.088 to 7.241), `74126a3` against `d3fb70f`, `bench_replay hard <codec> 60 5` on the `dam`
stream, seven alternating invocations pinned to one CPU, a portable Release build with SIMD,
fast-math and LTO, on a lightly loaded machine (one-minute load 1.5 falling to 1.1; [benchmarks](testing.md#the-pulsed-state)).
AMBE+2's output is identical in the two builds (its checksums are the same), but it is not an
untouched control: its timed path runs the changed code too (`prepare_frame()`, `pulsed_frame()`,
three more `prepare_render()` calls and the state store), and the cause of its +0.9% (ranges that do
not overlap) was not examined. Timed twice more after the review, the three builds alternating
(`74126a3`, `d3fb70f` and the final code, with the validity counter; the same method, one-minute load
about 1.0): D-STAR 6.442 → 6.695 → 6.682 us and 6.454 → 6.697 → 6.692 us (+3.7% and +3.7% for the
final code against `74126a3`), IMBE 6.851 → 7.121 → 7.227 us and 6.849 → 7.146 → 7.180 us (+5.5% and
+4.8%), AMBE+2 7.022 → 7.032 → 7.020 us and 7.014 → 7.014 → 7.016 us (0.0% and 0.0%). Against
`74126a3` the ranges do not overlap in D-STAR or IMBE and do in AMBE+2, whose +0.9% of the first run
does not recur. Against `d3fb70f` the final code's ranges overlap in D-STAR and AMBE+2 in both runs;
in IMBE the first run puts it 1.5% above with ranges apart (7.100 to 7.146 against 7.175 to 7.252)
and the second 0.5% above with ranges overlapping (7.109 to 7.190 against 7.162 to 7.318): these
give no bound on the counter's own cost; the change against `74126a3`
moves between runs (IMBE +5.6%, +5.5%, +4.8%). Under
callgrind (the first run's builds, the P25 `dam` stream once) the head executes 2.5% more
instructions: the rendering of the emitted pulses (`render_pulse()` for them, and the setup of the
pulsed bands' oscillators, `atan2` and `sincos`) and the code tests (`pulsed_code()`,
`imbe_code_pitch()`, `pulse_bands()`) account for most of the increase, and `prepare_render()`,
inlined in the baseline's build, is a function of its own in the head's (its instructions moved, not
added). Why the time rises more than the instructions was not examined.

**After the review.** The review's fixes add P25's validity counter (clause 8, `6efedb7`; moving only
for a model at a P25 code's pitch since `43850cc`) and
restrict the reserved codes to D-STAR's and P25's renderings; on DVSI's material they change
nothing. Every decode of the 315 bit files at the top of the three modes' directories (`irstia` and
`mark` left out: 106 D-STAR, 103 P25 and 106 rate-33 files, the soft-decision ones decoded as soft)
is identical, sample for sample and in every frame's mute flag, to the decode of `d3fb70f` (the
library code of `33e7e58`); all 285 encodings are identical to `74126a3`'s; the scoreboard, run
again in one combined run, `74126a3` against the final code, passes the four sets with the same 87
checks and the same figures, every measurement of the candidate's scoreboard equal to the run
against `d3fb70f`; the differential harness, which now also forges the read flag and the counter,
gives `74126a3`'s output in the ten configurations; and the golden hashes are unchanged and pass
under glibc 2.39 (scalar, SSE2 and AVX2). The counter's effect is shown by the tests alone (a
stream from an encoder that does not use the codes, the repeat limit's mute: clause 8).

**Gate sets**, declared in the branch's first commit, `525d017`, before any build of the branch was
scored (the study's scratch build had been measured on the development and validation vectors and
on the other material above before that, and the baseline's own scoreboard of the four partitions,
a build of `git archive 74126a3` decoding its own encodings, was read to record each measure's
coverage), and evaluated in one combined run over the development, validation, tone and channel
partitions, `74126a3` against `d3fb70f`:

| Gate set | Result | What it establishes |
|---|---|---|
| `pulsed-state-development` | passes all 18 checks | on the vectors the rule was found on, within the limits declared: D-STAR's and P25's correlation means rising by at least 0.0005 with intervals above 0 (observed +0.000938 and +0.000841) and no vector falling by more than 0.0005 (observed: none falls), the distance from DVSI's decode falling by at least 0.10 dB with its interval below 0 (observed -0.1692 and -0.1605 dB) and rising by no more than 0.02 dB on any vector (observed: falling on every vector), and the envelope correlation, the distance from the input (allowed to rise by 0.02 dB; observed +0.0046 and +0.0121 dB) and the level within theirs on the mean only, each check measured on all 20 vectors (`every_vector`); rate 33's speech decodes and every encoder's speech rows unchanged |
| `pulsed-state-validation` | passes all 18 checks | the same on the validation vectors, the distance from the input allowed to rise by 0.04 dB on the mean (observed +0.0234 and +0.0258 dB), also read in the study |
| `pulsed-state-tones` | passes all 18 checks | rate 33's speech and tone decodes unchanged, the D-STAR and rate-33 encoders' per-vector tone statistics (agreement, extra tones, level over the frames where both send the same tone) unchanged, and D-STAR's and P25's five tone measures within their limits on the vectors both runs measure (observed: unchanged); not that the same tone frames or bits were sent |
| `pulsed-state-channel` | passes all 33 checks | the per-vector mute statistics unchanged in every mode (not by itself that no frame's mute decision changed), D-STAR's and P25's mean distance from DVSI's decode on the loud frames (`lsd_common_db`) not rising past 0.01 dB with a coverage of at least 0.1 (observed: falling on every vector), rate 33's decodes and every encoder's speech rows unchanged |
| `pulsed-state-final` | passes all 18 checks | scored once, on 2026-11-27, after the two reviews and their fixes, on `irstia` and `mark` ([the final partition](#the-final-partition)): D-STAR's and P25's correlation means rising by at least 0.0003 with intervals above 0 and neither vector falling by more than 0.0005 (observed +0.003707 and +0.003896, both vectors rising in both codecs), the distance from DVSI's decode falling by at least 0.05 dB with its interval below 0 and rising by no more than 0.02 dB on either vector (observed -0.3686 and -0.3406 dB, falling on both), the envelope correlation, the distance from the input (allowed to rise by 0.05 dB; observed +0.0179 and +0.0266 dB) and the level within theirs on the mean; rate 33's decodes and every encoder's rows unchanged |

Notes on the declarations. A check that the mean rises (falls) by at least a value is the
scoreboard's `increase` (`decrease`) rule, which also requires the 95% interval above (below) 0.
Every paired check of the development, validation and final sets requires every vector
(`every_vector`), and the baseline measured each of the development and validation ones on all 20
and all 5 vectors; `every_vector` makes a check need a measurement on every vector, it does not make
a limit on the mean a limit on each vector. The tones set requires every vector for the detector
agreement and the D-STAR and rate-33 encoders' extra-tone rates, which cover all 52, and its
description records the coverage of the others from the baseline's own scoreboard; a check on a
metric that applies to part of a partition cannot require every vector, so a vector on which the
metric went missing in both runs would pass it unnoticed, and in this run each such check covered
the vectors listed there (the channel checks the four hard- and the four soft-decision vectors
each). The identity counters (`speech.pcm_changed_vectors`, `tone.pcm_changed_vectors`,
`channel.pcm_changed_vectors` and the encoder-row counter) count the vectors present in the run, so
they were evaluated on the combined run; the encoder-row counter skips a vector whose baseline holds
no row hash, and in this run all 75 speech entries (25 per mode) carry an encoder-row hash and a PCM
hash in both reports, and all 156 tone and 24 channel entries a PCM hash. Bit identity of the
encodings is the separate comparison of the 285 encodings above, and frame-for-frame mute identity
the separate comparison of the 30,600 channel-error frames: the gates establish neither. The
distance from the input rising is allowed by the limits, and observed: it is this change's cost
(above). Where the descriptions, which stay as declared, say more than the checks or the
measurements:

- All five state the rule without its exceptions for the mid-frame model: the model a frame marked
  `MBE_FRAME_MARK_TONE` or `MBE_FRAME_MARK_PREVIOUS_HELD` puts at mid-frame, and a mid-frame model
  with a voiced band (a caller's `MBE_FRAME_MARK_PITCH_ABSENT` on a frame with a code can hold voiced
  bands there), take no pattern ([D-STAR and P25: the pulsed state](#d-star-and-p25-the-pulsed-state),
  clause 7). No decode of DVSI's bits reaches these but at a P25 tone model, which has voiced bands
  and so no code.
- The tones set's description says the scratch build changed some D-STAR and P25 tone decodes "in
  frames at a tone's edge". The library changes 15 and 14 tone decodes in 1 to 796 frames each
  (above); the scratch build's own caches show the same order of change against the head (`pzx3`
  against the head's build: P25 `dtmfvbad` in 785 frames, D-STAR `dtmfvbad` in 797), and its
  single-frame changes lie in frame 0 (its tests', above). The set gates the tone measures, not
  where the decodes change.
- The channel set's description gives the scratch build's residual "at all four hard-decision rates
  in both codecs"; those are the study's figures, `pzx3` against `k1` on the frames that differ
  between those two; the library's builds' figures, hard and soft, are above.
- All five say that each clause "was chosen against its alternatives on DVSI's decodes of the
  development and validation vectors by a scratch build of the library". For the pulse's position,
  the reach and the gain that is wrong: when the sets were declared those rested on the study's
  earlier numpy model of P25 frames deep inside runs with no voiced band. The position and the reach
  were measured on scratch library builds after the declaration ([D-STAR and P25: the pulsed
  state](#d-star-and-p25-the-pulsed-state), clauses 2 and 3); the gain was not measured again.
- All five describe the rule as it stood before the review: they do not have P25's validity
  counter (clause 8), added after it, nor do they say that the slots of the subframe state carry a
  read flag and the counter. With the counter no decode of a non-held-out file changes (below), so
  the scores above stand for the final code; `pulsed-state-final` was scored with it (below).
- The channel set's description says the frames `lsd_common_db` takes have "no muted frame at or
  next to them"; its parenthesis names what the check excludes (a loud frame DVSI mutes, or a frame
  either run muted), and a quiet frame DVSI mutes next to them is admitted (the review found seven
  selected D-STAR and P25 frames with one, zero-based frame 625 of D-STAR `dam_e5_hd` among them).
- The tones set's description says the tone measures "do not see a tone's first and last frames".
  That is not guaranteed: the selection is three adjacent acoustic classifications of overlapping
  320-sample windows, and the review passed a synthetic step from 3500 to 3531.25 Hz through it as
  steady in the window spanning the step and in the windows beside it. The checks constrain the
  detector's results on the windows it calls steady, no more.

#### The final partition

`pulsed-state-final` was scored once, on 2026-11-27, after the two reviews and their fixes. The
baseline is the build of `3840d2b`, whose library code is `74126a3`'s (its decodes of the 25
development and validation vectors are those of the stack's head, in all three modes, 25 of 25
each); the candidate is the build of `5cfdf37`, whose library code is `43850cc`'s (its decodes of
the 25 vectors are those of the build scored above, in all three modes, 25 of 25 each). Both builds
aligned freely first and were then scored at their joint lags on `irstia` and `mark`, both runs'
encodings decoded by the candidate's decoder. Neither the rule nor any build of it had been measured
on these two vectors. The set passes all 18 checks (none fails); it required a gain, and the gain is
there on both vectors in both codecs. The library's decode of DVSI's bits, baseline → candidate:

| Codec | Waveform corr. vs DVSI, sample against sample at lag 0, mean | `irstia` | `mark` | LSD vs DVSI (dB), mean (`irstia`, `mark`) | Envelope corr. vs DVSI, mean | Level vs DVSI (dB), change of the mean | LSD vs input (dB), mean |
|---|---|---|---|---|---|---|---|
| D-STAR | 0.985589 → 0.989296 (+0.003707) | 0.980544 → 0.984752 | 0.990634 → 0.993839 | 1.3658 → 0.9971 (-0.3686; 1.4553 → 1.0756, 1.2762 → 0.9187) | 0.97917 → 0.99589 (+0.01673) | -0.0104 | 5.5172 → 5.5351 (+0.0179) |
| P25 | 0.986919 → 0.990815 (+0.003896) | 0.984030 → 0.988222 | 0.989808 → 0.993407 | 1.3346 → 0.9941 (-0.3406; 1.2989 → 0.9609, 1.3704 → 1.0273) | 0.98461 → 0.99593 (+0.01132) | -0.0118 | 5.3576 → 5.3842 (+0.0266) |
| Rate 33 | 0.985956, unchanged | 0.982792 | 0.989120 | 1.5074, unchanged | 0.9783, unchanged | 0 | 5.5216, unchanged |

Rate 33's two decodes are the same in both runs (equal PCM hashes), and no encoder row changes in
any mode; all six entries (`irstia` and `mark` in each mode) carry an encoder-row hash and a PCM
hash in both reports, as the notes above ask of the identity counters. The baseline's D-STAR figures
differ from those of the row-17 train change's held-out run ([November 2026: the pulse train through
a row-17 frame](#november-2026-the-pulse-train-through-a-row-17-frame): 0.980559 for `irstia`, a mean
of 0.985597, and a distance from DVSI's decode of 1.5116 dB against 1.3658 here) although the decode
is the same (equal PCM hashes): each run scores both builds at its own joint lags, and the joint lag
of D-STAR's `irstia` against DVSI's decode is -8 samples there and 0 here (that vector's distance
from DVSI's decode 1.747 dB there, 1.4553 dB here). The two runs' baselines are not one number.

What this establishes: on two vectors neither the rule nor any build of it was measured on before,
the gain declared in advance is there, larger than on development and validation (+0.0037 and
+0.0039 of correlation against +0.0008 to +0.0012), and so is the cost against the input (+0.018 and
+0.027 dB of log-spectral distance from it; PESQ was not measured on the held-out vectors). Two
vectors give no interval beyond what the gate's own rule computes.

### November 2026: rate 33's pulsed rows

The pulsed state's head (`2b6b6c2`) against rate 33's pulsed rows ([Rate 33: the pulsed
rows](#rate-33-the-pulsed-rows)): a rate-33 voice frame of V/UV rows 18 to 31 plays the row's
pulsed 500 Hz channels, each pulsed model its own pulse at the time the pitch code gives. The study
that found the rule read the development and validation vectors, DVSI's rate-33 tone, channel-error
(hard and soft decisions) and DTX vectors, never `irstia` or `mark`. Its scratch build (outside the
repository) was measured on those, with the scoreboard's own measures for rate 33 too, before the
gate sets were declared, so none of those partitions is out of sample here; the final partition was
held out and scored once, at the end (below). The baseline is a build of `git archive 2b6b6c2`, the
candidate the branch's library code at `a1398df`; the branch's later library commit, `1bc657f`,
changes how the compiler builds the voiced train in the 32-bit Release build (below, **Other modes
and the encoders**) and with it moves rate 33's own float samples on frames of rows 18..31 in two
builds, the 32-bit Release build and dev-release (GCC with fast-math, LTO and `-march=native`), by
up to 3 ulps of a frame's peak ([Rate 33: the pulsed rows](#rate-33-the-pulsed-rows), **The
pointer**); in the scoreboard's configuration its rate-33 decodes of the 88 vectors of [Rate 33:
the pulsed rows](#rate-33-the-pulsed-rows), the scoreboard's 85 rate-33 vectors among them, equal
`a1398df`'s file for file. The review's fix, `326c89c` (a lost frame's repeat reads the word only
on its code's band-grid model, and a held form of rows 26..31 is no word), decodes all 315 bit
files below as `1bc657f` does: the 16-bit WAV files are equal byte for byte, and so are the
`--params` JSONL files, which hold every frame's flags, error count and bits and its decoded,
predicted and synthesized model as `mbe_quality_eval` serializes it (`w0`, `L`, `Vl`, `Ml`,
`log2Ml`, `gamma`, and the synthesized `PHIl` and `PSIl` from harmonic 1), but not `PHIl[0]`, slot
37 or every other field of `mbe_parms`: evidence about the samples and the serialized model, not
about caller state, for which the differential harness below is the evidence (the first review's
defect kept the samples and changed exactly a field the JSONL omits). The scores below stand for the
final code, `1502dc0`, which the held-out run scored (later commits change comments, tests and
documents only); both
aligned freely first and were then scored at their joint lags, both runs' encodings decoded by the
candidate's decoder (Release, SIMD, IEEE arithmetic).

**The library against the study's scratch build** is in [Rate 33: the pulsed
rows](#rate-33-the-pulsed-rows): 60 of the 88 vectors identical, the others equal but for one-unit
differences from the pulse times' rounding, the scratch build's voiced train from silence (and the
frames it reaches) and its frame 0 carried over between `mbe_quality_eval`'s two passes; the figures
of the development and validation vectors are the study's to the digits printed, those of five tone
vectors (`alltone`, `cp0`, `cp1`, `dtmf`, `dtmfvbad`) and of `dam_e10_hd` are not.

**Speech**, rate 33 decoding DVSI's bits, on the scoreboard: the waveform correlation with DVSI's
decode sample against sample at lag 0 (`vs_dvsi.waveform_corr_unshifted`, over the frames where
DVSI's decode lies within 40 dB of its loudest), with the mean paired change, its 95% interval and
the smallest change; the log-spectral distance from DVSI's decode (`vs_dvsi.lsd_db`), the envelope
correlation with it (`vs_dvsi.env_corr`), the level against it (`vs_dvsi.level_offset_db`) and the
log-spectral distance from the input (`ours.lsd_db`), partition means, the library's builds:

| Partition | Waveform corr. vs DVSI | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| development | 0.992449 → 0.993228 (+0.000778, 95% interval 0.000300 to 0.001396; all 20 rise, the least `dam10` +0.0000073, the most `p01mirs` +0.002618) | 1.2291 → 1.1747 (-0.0545; all 20 fall, by 0.0036 at the least, `dam10`, and 0.1149 at the most, `fambm22a`) | 0.99209 → 0.99703 (+0.00494) | +0.0139 → +0.0106 | 6.2709 → 6.2789 (+0.0080; all 20 rise, the most `ucarm20a` +0.0253) |
| validation | 0.992142 → 0.993595 (+0.001452, interval 0.000473 to 0.002828; all 5 rise, the least `t03` +0.000090, the most `t04` +0.003740) | 1.0641 → 0.9878 (-0.0764; all 5 fall, by 0.0309 at the least, `t03`, and 0.1252 at the most, `t04`) | 0.98840 → 0.99172 (+0.00331) | -0.0403 → -0.0468 | 4.8352 → 4.8515 (+0.0163; all 5 rise, the most `t04` +0.0357) |

These are the figures the gate sets' descriptions give for the study's scratch build, to the digits
given there. **The distance from the input rises** on every vector (by +0.000408 at the least on
development, `lvl`, and +0.0076 on validation, `tambf32b`): it is this change's cost. PESQ against
the input moves by -0.0094 ± 0.0071 (9 of 25 vectors up; a 95% interval over the 12 recordings, t
with 11 degrees of freedom, -0.0250 to +0.0062, so neither a cost nor its absence is established)
and against DVSI's decode by +0.0867 ± 0.0075 (+0.0702 to +0.1032, 25 of 25 up; [Rate 33: the pulsed
rows](#rate-33-the-pulsed-rows)). D-STAR's
and P25's partition means and decodes are unchanged, and so is every encoder's bits.

**Tones.** By the scoreboard's hashes 9 of the 52 rate-33 tone decodes change (`alert`, `alltone`,
`cp0`, `cp1`, `cpvbad`, `dtmf`, `dtmf60`, `dtmf180`, `dtmfvbad`) and no D-STAR or P25 one; the
frames that change and their residuals are in [Rate 33: the pulsed rows](#rate-33-the-pulsed-rows).
Of the five tone measures the gates check, the frequency, level and balance errors and the share of
windows playing a wrong tone do not move on any vector; the detector's agreement moves on one,
`alltone` (+0.003175; +0.00006 in the mean, no vector down); the encoders' tone statistics do not
move. `tone.alignment_corr`, which no gate checks, moves on the nine vectors (from -0.003136, `cpvbad`,
to +0.005765, `cp1`). The frequency, level, balance and wrong-tone checks see only the windows where
both decodes hold a tone the detector recognises. The detector's agreement is taken over every
window whose DVSI classification equals its two neighbours' (`compare_tones()`), with or without a
tone: such a window agrees when both decodes hold no tone or the same tone, and disagrees when one
holds a tone and the other none or another (a window where DVSI's decode holds a steady tone and
ours none counts as checked and disagreeing: the saved `dtmf` comparison has one, at hop 3278);
three consecutive windows of DVSI's decode without a tone, voice frames between and around the tones
among them, count too, so the agreement does see some of the frames where the change acts.

**Channel errors.** No mute decision changes, frame for frame: the flags of every frame of the 24
channel-error decodes (8 vectors in each of the three codecs) are the same in both builds, as are
those of every frame of the 315 bit files below. Every rate-33 channel-error decode changes, and no
D-STAR or P25 one. The scoreboard's log-spectral distance from DVSI's decode on the loud frames that
neither run nor DVSI mutes (`lsd_common_db`; [November 2026: the pulsed
state](#november-2026-the-pulsed-state) gives its selection) falls on all eight: `dam_e1_hd` to
`dam_e10_hd` by 0.0125, 0.0144, 0.0146 and 0.0229 dB, `dam_e1_sd` to `dam_e10_sd` by 0.0126,
0.0145, 0.0148 and 0.0073 dB, the means -0.0161 (hard, coverage 0.750 at the least, `dam_e10_hd`)
and -0.0123 dB (soft, coverage 0.991). The gate sets' description gives -0.0122 for the hard
decisions, the study's scratch build's: that build also changed the voiced frames after the mutes of
`dam_e10_hd`, which the library leaves as they were.

**Other modes and the encoders.** Of the bit files at the top of each mode's directory of DVSI's
vectors (`irstia` and `mark` left out: 106 D-STAR, 103 P25 and 106 rate-33 files, the tone,
channel-error, DTX and synthetic ones among them, the soft-decision ones decoded as soft), every
D-STAR and P25 file decodes identically in the two builds, sample for sample and in every frame's
flags (113,596 and 105,033 frames); 56 of the rate-33 files decode otherwise (every speech,
channel-error and DTX vector, the nine tone vectors above, and the synthetic `grbge`, `ptrain1`,
`ptrain3` to `ptrain5`, `ptrain7` to `ptrain9`, `sine0_1k`, `sine0_4k` and `xfer`), the flags of
all 106 the same. All 285 encodings (95 input files, `irstia` and `mark` not among them, three
encoders) are identical between the two builds; the 95 rate-33 encodings hold no voice frame of
rows 18 to 31 (94,928 frames) and decode identically with both builds' decoders.
In each of the ten local configurations (dev-debug; dev-release, GCC Release with SIMD, fast-math,
LTO and `-march=native`, the preset's `MBELIB_ENABLE_NATIVE=ON`; asan-ubsan-debug; notones-debug; a
Debug build with SIMD; a Debug `-march=x86-64-v3` build; Clang Release with SIMD, fast-math and LTO,
once with `MBELIB_ENABLE_NATIVE=ON`, so `-march=native`, and once with `-march=x86-64-v3` and
`MBELIB_ENABLE_NATIVE=OFF`; and 32-bit Release and Debug), a differential harness outside the repository,
built against the baseline's library and the head's, gives the same output byte for byte: the pitch
and harmonic count every D-STAR and P25 code decodes to and every rate-33 code at rows 16 and 18,
and hashes of the output and of all three returned states of D-STAR's, P25's, rate 33's and
ProVoice's paths (3,000 frames each) and `mbe_synthesizeSpeechf()`, and of the three encoders' bits
and `cur_mp` (1,500 frames), over constructed streams whose caller states hold reserved pitches,
every frame mark, pulsed-row words (held or not, with and without an own mask, and values beside
them), values in slots 34 to 37 and states without the subframe mark; rate 33's own stream holds no
voice frame of rows 18 to 31. As first run the harness gave rate 33's caller states no word (a word
there on a band-grid model would be such a frame's), and its report did not say so: that hid the
repeat loader's reading of a word on any model, which the review found and the fix round
corrected (a word is read only on its band-grid model, clause 7). Run again with the exclusion
removed, rate 33's states holding the same words as the other paths' (on these coded pitches no
word is a band-grid model's), against `2b6b6c2` with the head's library, every configuration gives
the same output byte for byte. At `a1398df` the 32-bit Release build differed there in D-STAR's,
P25's, rate 33's and ProVoice's paths: the two new direct calls of `render_pulse()` made GCC clone it
out of the voiced train, and x87's excess precision then moved most voiced frames' samples (by up to
3 ulps of the frame's peak, and five 16-bit samples by one unit through the int16 entry points; the
states were equal); `1bc657f` reaches it through a
const volatile pointer, and the ten configurations agree since (the pointer's own effect on the
rows' samples: [Rate 33: the pulsed rows](#rate-33-the-pulsed-rows), **The pointer**). The golden PCM hashes that moved (rate 33's, scalar, SSE2 and AVX2; the
dispatched AVX2 equal to a `-march=x86-64-v3` build's) pass under glibc 2.39 (an `ubuntu:24.04`
container) as under the local 2.44.

**Cost.** The median complete hard-decision frame decode of AMBE+2 moves by +1.5% (7.308 → 7.420
us; the seven invocations' medians 7.162 to 7.592 and 7.369 to 7.533 us) and, run again, +3.6%
(7.153 → 7.411 us; 7.126 to 7.544 and 7.358 to 7.792), the ranges overlapping both times; D-STAR
and IMBE, whose timed paths run the changed code too (`pulsed_rows()` and the slot's store) and whose
outputs are the same (equal checksums), move by +0.9% and -1.2%, then -0.2% and -0.4%. `2b6b6c2`
against `1bc657f`, `bench_replay hard <codec> 60 5` on DVSI's `dam` stream of each codec (its bits
reframed), seven alternating invocations pinned to one CPU, a portable Release build with SIMD,
fast-math and LTO, on a machine with a browser open (one-minute load 1.2 to 1.5): the change's cost
is within that spread ([benchmarks](testing.md#rate-33s-pulsed-rows)).

**Gate sets**, declared in the branch's first commit, `5af6391`, before any build of the branch was
scored (the study's scratch build had been measured on the development, validation, tone, channel and
rate-33 DTX vectors before that, and the baseline's own scoreboard of the four partitions, a build of
`git archive 2b6b6c2` decoding its own encodings, was read to record each measure's coverage), and
evaluated in one combined run over the development, validation, tone and channel partitions,
`2b6b6c2` against `a1398df`:

| Gate set | Result | What it establishes |
|---|---|---|
| `r33-pulsed-rows-development` | passes all 12 checks | on the vectors the rule was found on, within the limits declared: rate 33's correlation mean rising by at least 0.0005 with its interval above 0 (observed +0.000778) and no vector falling by more than 0.0001 (observed: none falls), the distance from DVSI's decode falling by at least 0.03 dB with its interval below 0 (observed -0.0545 dB) and rising by no more than 0.005 dB on any vector (observed: falling on every vector), and the envelope correlation, the distance from the input (allowed to rise by 0.02 dB on the mean; observed +0.0080 dB) and the level within theirs on the mean only, each check measured on all 20 vectors (`every_vector`); D-STAR's and P25's speech decodes and every encoder's speech rows unchanged |
| `r33-pulsed-rows-validation` | passes all 12 checks | the same on the validation vectors, the distance from the input allowed to rise by 0.04 dB on the mean (observed +0.0163 dB), also read in the study |
| `r33-pulsed-rows-tones` | passes all 15 checks | D-STAR's and P25's speech and tone decodes unchanged, the D-STAR and rate-33 encoders' per-vector tone statistics (agreement, extra tones, level where both send the same tone) unchanged, and rate 33's five tone measures within their limits on the vectors both runs measure (observed: unchanged but for the detector's agreement on `alltone`, +0.003175); not that the frames between the tones, where the change acts, are right |
| `r33-pulsed-rows-channel` | passes all 33 checks | the per-vector mute statistics unchanged in every mode (not by itself that no frame's mute decision changed), rate 33's mean distance from DVSI's decode on the loud frames (`lsd_common_db`) not rising at all, with a coverage of at least 0.1 (observed: falling on every vector, -0.0161 dB hard and -0.0123 dB soft), D-STAR's and P25's decodes and every encoder's speech rows unchanged |
| `r33-pulsed-rows-final` | passes all 12 checks, scored once (below) | on `irstia` and `mark`, held out, a prediction: rate 33's correlation mean rising by at least 0.0003 with its interval above 0, which with two recordings asks for the gain on both vectors (observed +0.003599, both up), no vector falling by more than 0.0001, the distance from DVSI's decode falling by at least 0.02 dB with its interval below 0 (observed -0.151 dB) and rising on neither vector by more than 0.005 dB, the envelope correlation (falling by no more than 0.001), the distance from the input (rising by no more than 0.05 dB) and the level (moving by no more than 0.05 dB) within theirs on the mean; D-STAR's and P25's decodes and every encoder row unchanged |

Notes on the declarations. The identity counters were evaluated on the combined run: all 75 speech
entries (25 per mode) carry an encoder-row hash and a PCM hash in both reports, and all 156 tone and
24 channel entries a PCM hash. Bit identity of the encodings is the separate comparison of the 285
encodings above, and frame-for-frame mute identity the separate comparison of the frames' flags: the
gates establish neither. Where the descriptions, which stay as declared, say more than the checks
or the measurements:

- All five give the study's scratch build's figures for the scoreboard's measures; on the library's
  builds they are the same to the digits given but one, `channel.hd.lsd_common_db`: -0.0161 dB
  where the descriptions give -0.0122. The scratch build changed, besides the rule, the voiced
  frames that start from silence after a mute or a tone (its first interval's train), which the
  library leaves as the baseline plays them ([Rate 33: the pulsed rows](#rate-33-the-pulsed-rows)):
  in `dam_e10_hd`'s frames after its mutes that moved the hard-decision distance.
- The tones set's description says the tone measures "cover the scoreboard's steady tone windows,
  not the voice frames between and around the tones, which is where this change acts in these
  vectors". As the pulsed state's notes on its declarations say, the selection is three adjacent
  acoustic classifications of overlapping 320-sample windows, which does not guarantee that a
  window holds no voice frame; and the detector's agreement counts three consecutive windows with
  no tone as steady too, so it covers some voice frames by design. The checks constrain the
  detector's results on the windows it calls steady, no more.
- All five say "on the tone partition the five tone metrics unchanged (tone.detector_agreement
  +0.00006 in the mean, no vector down)": four were unchanged; the detector's agreement moved, on
  `alltone`.
- All five say "The distance from the input speech rising is the change's measured cost, which these
  checks allow within their limits and do not hide". The development and validation sets check it
  (`ours.lsd_db` on the mean); the tones and channel sets have no check on the distance from the
  input, so there the sentence bounds nothing.
- "No build was measured on the held-out partition" is true of this change's builds (the study's
  scratch build and the branch's); earlier changes' builds were scored on it, once each, as their
  sections record.

**The final partition**, scored once, on 2026-11-27, after the two reviews and their fixes: the
build of `git archive 2b6b6c2` against the build of `1502dc0` (`build/tools-final`, the
scoreboard's Release configuration), both aligning freely first and then scored at their joint lags,
both runs' encodings decoded by the candidate's decoder, gate set `r33-pulsed-rows-final`: passed,
12 of 12 checks. The gate lines:

```
PASS dstar/final speech.pcm_changed_vectors max: value=0
PASS p25/final speech.pcm_changed_vectors max: value=0
PASS dstar/final encoder.rows_changed_vectors max: value=0
PASS p25/final encoder.rows_changed_vectors max: value=0
PASS r33/final encoder.rows_changed_vectors max: value=0
PASS r33/final vs_dvsi.waveform_corr_unshifted increase: mean delta=+0.003599 ... min=+0.00319 (and all_no_decrease)
PASS r33/final vs_dvsi.lsd_db decrease: mean delta=-0.151 ... max=-0.1462 (and all_no_increase)
PASS r33/final vs_dvsi.env_corr no_decrease: mean delta=+0.01371
PASS r33/final ours.lsd_db no_increase: mean delta=+0.003067 ci95=(-0.00661, 0.01274) worst=('irstia', 0.012739)
PASS r33/final vs_dvsi.level_offset_db abs_delta_max: mean delta=-0.01399 worst=('irstia', -0.015191)
```

Rate 33, per vector:

| Vector | Waveform corr. vs DVSI | LSD vs DVSI (dB) | Envelope corr. vs DVSI | LSD vs input (dB) | Level vs DVSI (dB) |
|---|---|---|---|---|---|
| `irstia` | 0.982794 → 0.985984 (+0.003190) | 1.409731 → 1.263516 (-0.146215) | 0.975280 → 0.994251 | 5.562402 → 5.575142 (+0.012739) | -0.031362 → -0.046553 |
| `mark` | 0.989120 → 0.993127 (+0.004008) | 1.330127 → 1.174357 (-0.155770) | 0.987956 → 0.996401 | 5.511256 → 5.504649 (-0.006606) | -0.114772 → -0.127559 |

The partition's mean correlation goes 0.985957 → 0.989556. D-STAR and P25: every measurement of
both vectors is equal in the two runs; all six entries carry both PCM hashes and both encoder-row
hashes in both reports, the four D-STAR and P25 PCM hashes and all six row hashes equal, the two
rate-33 PCM hashes different. The prediction, a gain on both held-out vectors, is met. The mean gain
in correlation there, +0.003599, is larger than on development (+0.000778) and validation
(+0.001452); nothing was measured that explains the difference. The distance from the input rises
on `irstia` (+0.012739 dB) and falls on `mark` (-0.006606 dB) at the lag both runs are scored at;
against the baseline's own alignment of `mark` it would read as a rise of 0.023873 dB (below), so the
sign of the change on `mark` is not established. The gate's limit, 0.05 dB on the mean, holds either
way: +0.0031 dB at the joint lags, +0.0183 dB with the baseline's own figure for `mark`. It was one
scoring, after the two reviews and their fixes, of the code of `1502dc0`.

The baseline's rate-33 figures here are not those the row-17 run's table above gives for the same
decode (rate 33's two PCM hashes, 0x4F2FBCFE and 0x47687202, are those of rate 33's held-out run of
the mid-frame pitch wave too), and the lags account for the difference. `joint_lags()`
(`tools/quality/dvsi_scoreboard.py`) scores both runs at one lag per comparison: their common lag,
or, where the two alignment searches disagree, the lag of the run that aligned with the higher
correlation (the baseline's on a tie). On `irstia`, our decode against DVSI's (`lags.vs_dvsi`), the
baseline aligning freely takes -8 samples (its free report: distance from DVSI's decode 1.684767 dB,
envelope correlation 0.968730, correlation 0.982792, level -0.031336 dB) and the candidate 0, and the
joint lag is 0, where the baseline scores 1.409731, 0.975280, 0.982794 and -0.031362. In the row-17
run the two decodes of the vector were one file, both aligned at -8, and its table gives the -8
figures (1.6848; its envelope mean, 0.9783, is the mean of 0.968730 and `mark`'s 0.987956, and its
correlation mean, 0.985956, that of 0.982792 and `mark`'s 0.989120). On `mark`, our decode against
the input (`lags.ours`), the baseline's free alignment takes 360 samples, the candidate's 352, and the
joint lag is 352. The baseline's distance from the input is 5.480776 dB at 360 (the row-17 table's
partition mean, 5.5216, is the mean of that and `irstia`'s 5.562402) and 5.511256 dB at 352, so the
fall of 0.006606 dB on `mark` is at the joint lag; against the baseline's own alignment the
candidate's 5.504649 is 0.023873 dB higher. Every other figure of the table above is equal between
the free and joint reports of each run: the candidate's two reports are equal in every field, D-STAR's
and P25's entries in the baseline's too, and the baseline's rate-33 entries differ only in those two
comparisons and in the fields taken with them (pitch statistics, band and crest deltas).

### November 2026: DVSI's low band

Rate 33's pulsed rows' head (`9a101f8`) against DVSI's low band ([DVSI's low
band](#dvsis-low-band)): a third-order Butterworth high-pass at 75 Hz in place of the second-order
one at 100 Hz, and the adaptive low cut, in D-STAR's, P25's and rate 33's renderings. The study that
found the rule read the development and validation vectors and DVSI's synthetic pulse-train vectors,
never `irstia` or `mark`; its scratch build (outside the repository, a worktree at `47d8cf3` with
environment switches) was measured on the development, validation, tone, channel and pulse-train
vectors, with the scoreboard's own measures too, before the gate sets were declared in the branch's
first commit, `4677a41`, so none of those partitions is out of sample here. The baseline is a build
of `git archive 9a101f8`, the candidate one of the branch's library commit, `0e3463e` (the later
commits change tests and documents only); both aligned freely first and were then scored at their
joint lags, both runs' encodings decoded by the candidate's decoder (Release, SIMD, IEEE
arithmetic).

**The library against the study's scratch build.** The 25 development and validation vectors,
decoded by the scoreboard's tools: D-STAR and P25 against the scratch decodes of the one-rate rule on
the stack's head (the study's `pll1`), rate 33 against those on its scratch pulsed rows (`rll1`).
D-STAR and P25: every P25 file differs in its frame 0 (79 samples, up to 11 units: the scratch build's
base differs from the library there too, as the pulsed state's results record); apart from that frame,
16 of the 50 decodes are equal sample for sample, and the other 34 differ in 135 samples of 132 frames,
each by one unit. A
copy of the library that forms the cut as the scratch build did (`f / fc` in float, its square
multiplied into the model's magnitude) equals the scratch decodes on all 50 files but P25's frame 0:
the 135 samples are the cut's arithmetic, the library's square formed in double and rounded once.
Rate 33: 3 of 25 equal; 209 samples differ, 133 with that copy: frame 0 of `fambm22a` (34 samples) and
`noisy` (31), the scratch build's row state carried between `mbe_quality_eval`'s two passes, and 68
one-unit samples in 65 frames, 43 of them frames of V/UV rows 18 to 31 and 22 the frame right after a
row whose own subframe pulses: the scratch build's pulse times in float, as the pulsed rows' results
record. Nothing else differs. The library's tones, channel-error, DTX and pulse-train decodes were not
compared with the scratch build's (it kept the tracker in elements its other paths did not carry as
the library does); the figures above and below are the library's.

**Speech**, decoding DVSI's bits, on the scoreboard: the waveform correlation with DVSI's decode at
lag 0 (`vs_dvsi.waveform_corr_unshifted`), the log-spectral distance from DVSI's decode
(`vs_dvsi.lsd_db`), the envelope correlation with it, the level against it and the log-spectral
distance from the input (`ours.lsd_db`), partition means, the mean paired change and its extremes:

| Codec, partition | Waveform corr. vs DVSI | LSD vs DVSI (dB) | Envelope corr. | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| D-STAR, development | 0.992583 → 0.993593 (+0.001010; 18 up, 2 down: `t02` -0.000095, `noisy` -0.000081; the most `lvl` +0.004003) | 0.9915 → 0.9414 (-0.0501; 19 fall, `dam10` rises +0.0092) | 0.99727 → 0.99763 (+0.00036) | +0.0233 → +0.0380 (+0.0147) | 6.4234 → 6.4456 (+0.0222; the most `ucarf15a` +0.0781, 5 fall) |
| D-STAR, validation | 0.994654 → 0.996149 (+0.001496; 5 up, the least `t03` +0.000041) | 0.8656 → 0.8084 (-0.0572; all fall, the least `t04` 0.0263) | 0.99410 → 0.99462 | -0.0099 → +0.0473 (+0.0572; `t04` +0.1565) | 4.8572 → 4.8792 (+0.0220; `t03` +0.0446) |
| P25, development | 0.987520 → 0.989085 (+0.001565; 20 up, the least `ucarm20a` +0.000481) | 1.3145 → 1.2510 (-0.0635; 19 fall, `dam10` +0.0019) | 0.99478 → 0.99549 (+0.00071) | +0.0431 → +0.0522 (+0.0091) | 6.3249 → 6.3475 (+0.0227; `ucarf15a` +0.0707, 6 fall) |
| P25, validation | 0.991984 → 0.993818 (+0.001833; 5 up, the least `t03` +0.000972) | 0.9375 → 0.8783 (-0.0592; all fall, the least `t04` 0.0167) | 0.99282 → 0.99340 | -0.0181 → +0.0298 (+0.0479; `t04` +0.1200) | 4.7467 → 4.7763 (+0.0295; all rise, `t03` +0.0462) |
| Rate 33, development | 0.993228 → 0.993793 (+0.000565; 16 up, 4 down: `t02` -0.000330, `noisy` -0.000163, `clean` -0.000124, `fambm22a` -0.000026) | 1.1747 → 1.1190 (-0.0557; 19 fall, `dam10` +0.0105) | 0.99703 → 0.99733 (+0.00030) | +0.0106 → +0.0289 (+0.0183) | 6.2809 → 6.3063 (+0.0254; `ucarf15a` +0.0880, 6 fall) |
| Rate 33, validation | 0.993595 → 0.994696 (+0.001102; 4 up, `t03` -0.000085) | 0.9878 → 0.9345 (-0.0533; all fall, the least `t04` 0.0368) | 0.99172 → 0.99226 | -0.0468 → +0.0088 (+0.0557; `t04` +0.1391) | 4.8515 → 4.8769 (+0.0254; all rise, `t03` +0.0478) |

These are the gate sets' measures; the study's scratch build gave the same to the digits the brief of
this change quoted, the vectors that fall the same high voices (`t02`, `noisy`, `clean`, `t03`), where
the cut does nothing and only the fixed filter changes. The whole-file correlation at lag 0 (the
study's measure, mean over each partition's vectors): P25 0.987487 / 0.991983 → 0.989053 / 0.993816,
D-STAR 0.992583 / 0.994653 → 0.993593 / 0.996148, rate 33 0.993226 / 0.993591 → 0.993791 / 0.994693,
each the scratch build's to the sixth decimal. **The distance from the input rises**: by +0.022 to
+0.030 dB in the mean in all six, the most on `ucarf15a` (+0.071 to +0.088), falling on `p01mirs`
(-0.024 to -0.030), `clean`, `fambm22a`, `noisy` and some of the car-noise versions: DVSI's decoders
remove fundamentals the input speech has, and a decode that follows them moves away from the input.
It is this change's cost. PESQ (the 25 vectors; a 95% interval over the 12 recordings, t with 11
degrees of freedom, from the standard error clustered by recording): against the input D-STAR +0.0023 ±
0.0017 (-0.0014 to +0.0060, 14 of 25 up), P25 -0.0024 ± 0.0017 (-0.0061 to +0.0013, 12 up), rate 33
+0.0001 ± 0.0015 (-0.0032 to +0.0034, 12 up): neither a cost nor a gain established; against DVSI's
decode +0.0029 ± 0.0018 (-0.0011 to +0.0069, 20 up), +0.0055 ± 0.0018 (+0.0015 to +0.0095, 20 up),
+0.0048 ± 0.0018 (+0.0008 to +0.0088, 23 up). The level against DVSI's moves by +0.009 to +0.057 dB
in the mean, the most on `t04` (+0.120 to +0.156), the lowest voice of the set (median f0 86 Hz): its
energy under 150 Hz (3.4 to 4.9% of it) rises by 1.6 dB and that from 150 to 300 Hz (23 to 27%) by
0.21 dB, the third-order filter at 75 Hz passing more there than the second-order one at 100 Hz
while its corner follows its own low pitch; that adds 0.12 to 0.16 dB.

**Tones.** By the scoreboard's hashes every tone decode of the three codecs changes (the AMBE tones
start at the new filter's phase; P25's tones are voice frames and its tone models play through the
new filter). No tone measure of D-STAR or rate 33 moves in the third decimal (the worst vector's
level error 0.2126 → 0.2137 and 0.2829 → 0.2847 dB, balance 0.1687 → 0.1688 and 0.1696 → 0.1697);
P25's level error 0.5897 → 0.5890 and balance error 0.3392 → 0.3601 (the worst vector, `cp31`:
+0.021); the share of wrong tones and the detector's agreement do not move on any vector; the
encoders' tone statistics do not move. In DVSI's tone vectors no P25 tone model has a component under
the corner (the lowest, ID 5, lies at 156.25 Hz); the voice frames between and around the tones are cut
as speech is.

**Channel errors.** No mute decision changes, frame for frame: the flags of every frame of the 24
channel-error decodes, and of every frame of the 315 bit files below, are the same in both builds.
Every decode changes. The scoreboard's log-spectral distance from DVSI's decode on the loud frames that
neither run nor DVSI mutes (`lsd_common_db`) falls in all six: hard decisions D-STAR -0.0583, P25 -0.0377,
rate 33 -0.0668 dB (coverage 1.000, 0.974, 0.750 at the least), soft decisions -0.0697, -0.0424, -0.0571
(0.937, 0.971, 0.991); `channel.hd.lsd_db` / `sd` 1.666 → 1.608 and 2.036 → 1.966 (D-STAR), 4.824 →
4.787 and 5.723 → 5.681 (P25), 1.565 → 1.498 and 1.563 → 1.506 (rate 33). Whether the tracker should
follow the decoded voicing or the rendered one (after the adaptive smoothing channel errors bring) can
hardly be told on these vectors: a copy of the library whose tracker reads the decoded voicing differs
from the library on 320 D-STAR, 297 P25 and 188 rate-33 frames that we do not mute, nearly all of them in
`dam_e5_hd` and the 10% vectors where neither decode follows DVSI's (a residual larger than DVSI's
energy); on those that do (residual under half of DVSI's energy: 3, 9 and 144 frames) the rendered
voicing is nearer DVSI's on 3 of 3, 7 of 9 and 81 of 144 (pooled residual 0.30992 against 0.31122,
0.25320 against 0.25325, 0.06007 against 0.06037). The library keeps the rendered model, as every rule
of its synthesis does.

**DTX.** D-STAR's and rate 33's DTX vectors (`clean_dtx`, `dam_dtx`, `noisy_dtx`), on the frames the
change moves (5,930 D-STAR, 5,766 rate 33), residual energy against DVSI's decode at lag 0, pooled, over
DVSI's: D-STAR 0.17761 → 0.17784 (the library's decode nearer on 3,530 frames, the branch point's on
2,351), rate 33 0.17401 → 0.17464 (2,867 against 2,846): no gain there.

**Other modes and the encoders.** Of the bit files at the top of each mode's directory (`irstia` and
`mark` left out: 106 D-STAR, 103 P25 and 106 rate-33 files, the soft-decision ones decoded as soft), every
decode changes and every frame's flags are the same (113,596, 105,033 and 113,596 frames). All 285
encodings (95 input files, three encoders) are identical between the two builds. In each of the ten local
configurations (dev-debug; dev-release, GCC Release with SIMD, fast-math, LTO and `-march=native`;
asan-ubsan-debug; notones-debug; a Debug build with SIMD; a Debug `-march=x86-64-v3` build; Clang Release
with SIMD, fast-math and LTO, once with `MBELIB_ENABLE_NATIVE=ON` and once with `-march=x86-64-v3`; and
32-bit Release and Debug), a differential harness outside the repository, built against the baseline's
library and the head's, gives the same output byte for byte for: the pitch and harmonic count every
D-STAR and P25 code decodes to and every rate-33 code at rows 16 and 18; hashes of the output and of all
three returned states of ProVoice's path (its float and 16-bit frame entry points, 3,000 frames), of
`mbe_synthesizeSpeechf()` and `mbe_synthesizeSpeech()` (3,000 frames), and of the three encoders' bits
and `cur_mp` (1,500 frames), over constructed streams whose caller states hold reserved pitches, every
frame mark, pulsed-row words, values in slots 34 to 37, every float class and both neighbours of each
bound in the low cut's slots 38 and 39 (and in-range trackers and corners), and states without the
subframe mark; and ProVoice, on the states a D-STAR, P25 and rate-33 stream of the same build left, plays
the same samples and returns the same states, but for the two slots, whatever those slots hold. D-STAR's,
P25's and rate 33's own streams change, as they are to. A first form of the change, before its commit,
failed the harness in dev-release (ProVoice's samples moved): the unvoiced shaping, written as one loop
for both filters, was compiled otherwise there; ProVoice's loop is kept as it was (`shape_subframe()`
apart from `shape_subframe_dvsi()`), and the ten configurations agree since, at `0e3463e`.
The golden PCM hashes that moved (D-STAR's and rate 33's, scalar, SSE2 and AVX2; the dispatched AVX2
equal to a `-march=x86-64-v3` build's) pass under glibc 2.39 (an `ubuntu:24.04` container) as under the
local 2.44.

**Cost.** The median complete hard-decision frame decode, `bench_replay hard <codec> 60 5` on DVSI's
`dam` stream of each codec, seven alternating invocations pinned to one CPU, a portable Release build
with SIMD, fast-math and LTO, `9a101f8` against `0e3463e`: AMBE+2 7.345 → 7.650 us (+4.2%; the
invocations' medians 7.327 to 7.682 and 7.582 to 7.970), D-STAR 6.809 → 7.100 (+4.3%; 6.791 to 7.078
and 7.037 to 7.192), IMBE 7.279 → 7.528 (+3.4%; 7.231 to 7.574 and 7.492 to 7.560), the ranges barely
overlapping; the new filter's gain and phase are three double-precision functions per voiced harmonic
of each model, against one before.

**Gate sets**, declared in the branch's first commit, `4677a41`, before any build of the branch was
scored (the study's scratch build had been measured on the development, validation, tone, channel and
pulse-train vectors before that, and the baseline's own scoreboard of the four partitions, a build of
`git archive 9a101f8` scoring itself, was read for each measure's coverage), and evaluated in one
combined run over the development, validation, tone and channel partitions, `9a101f8` against
`0e3463e`; all 75 speech entries carry an encoder-row hash in both reports, every one equal:

| Gate set | Result | What it establishes |
|---|---|---|
| `low-band-development` | passes all 24 checks | on the vectors the rule was found on: in each codec a mean gain in correlation with DVSI's decode with its interval above 0 (+0.001010, +0.001565, +0.000565), no vector falling by more than the stated little (the least -0.000095, +0.000481, -0.000330), the distance from DVSI's decode falling with its interval below 0 (-0.0501, -0.0635, -0.0557 dB) and rising by no more than 0.03 dB on any vector (the most +0.0092, +0.0019, +0.0105, all `dam10`), the envelope correlation, the distance from the input (allowed to rise by 0.05 dB in the mean; +0.022 to +0.025) and the level (0.08 dB in the mean, 0.25 on any vector) within theirs; every encoder's speech rows unchanged |
| `low-band-validation` | passes all 24 checks | the same on the validation vectors (the level allowed 0.1 dB in the mean: +0.048 to +0.057, `t04` +0.156 at the most), also read in the study |
| `low-band-tones` | passes all 21 checks | the D-STAR and rate-33 encoders' tone statistics unchanged and the decoders' tone measures within their limits on the windows the detector calls steady; not that the frames between and around the tones are right |
| `low-band-channel` | passes all 33 checks | the per-vector mute statistics unchanged in every mode (not by themselves that no frame's mute decision changed: the flags say that), the distance from DVSI's decode on the loud frames falling in every mode and decision type, every encoder's speech rows unchanged |
| `low-band-final` | not scored | the prediction for `irstia` and `mark` |

### November 2026: four rules of DVSI's P25 decoder

DVSI's low band's head (`f725a66`) against four rules of DVSI's P25 decoder, committed one at a
time: each frame's own pitch eq 46 on a 2^-19 grid (`6fb4ddc`), the mid-frame pitch between two
voiced models a floored parabola for the square root (`bfe59d4`; [P25 pitch](#p25-pitch)), the phase
kernel's reading of unvoiced bands with the unvoiced factor taken out (`072dde2`; [Phase
regeneration](#phase-regeneration)), and the frame's own models voiced by the 500 Hz columns, with
the mid-frame clause on the column-read states (`aff61e5`; [D-STAR and P25: mid-frame voicing by
columns](#d-star-and-p25-mid-frame-voicing-by-columns)). P25's rendering alone; ProVoice keeps its
nominal pitch, the geometric mean, the plain phase source and its decoded voicing, a decision (it
has no DVSI material). The study that found the rules read the development and validation vectors
and DVSI's synthetic pulse-train vectors, never `irstia` or `mark`; its scratch build (outside the
repository, a worktree at `47d8cf3` with environment switches) was measured on the development,
validation, tone and channel vectors, with the scoreboard's own measures too, before the gate sets
were declared in the branch's first commit, `f1be05d`, so none of those partitions is out of sample
here. The baseline is a build of `git archive f725a66`, the candidate one of `999647e`, an earlier
form of `aff61e5` whose `mid_pitch()` nested its branches otherwise (rewritten for the complexity
limit; the two builds decode all 98 P25 speech, channel-error, tone and pulse-train vectors byte for
byte alike, and the differential harness gives the same for every protected stream; the later
commits change documents only); both aligned freely first and were then scored at their joint lags,
both runs' encodings decoded by the candidate's decoder (Release, SIMD, IEEE arithmetic). The
measurements of the sections above that are not the scoreboard's were made on a copy of the library
with environment switches (outside the repository) that decodes every P25 speech, channel-error and
tone vector byte for byte as the library with all four rules on and as the baseline with all four
off.

**The library against the study's scratch build.** The scratch build's low band is its own two-rate
scratch form, not the library's, so the decodes cannot be compared sample for sample (the library's
baseline and the scratch build's, `plu1`, differ in 1.70 million of 8.85 million samples, a residual
of 1.5e-5 of the scratch decode's energy). Compared as the rules' effects (each build's decode with
the rules less its decode without them, over the 25 development and validation vectors): rules 1 to
3 against the scratch build's (`pk2` less `plu1`), the two effects differ by 0.036% of the scratch
effect's energy, pooled; with the own models by columns under the clause as it was (`pk3`), by
0.027%; rule 4's own effect against the scratch build's (`pk3` less `pk2`) by 0.001%. The largest part
of the difference lies in `dam40` frames 768 and 769, a frame played as a P25 tone model (a single
tone at harmonic 1) and the frame after it: the scratch build put the phase step on the tone model's
unvoiced bands, the library reads a tone model as decoded; the rest is spread at the level of the
two baselines' own difference, more in `dam_10ov`'s overloaded frames, where a small difference
moves the output limit. The correlation with DVSI's decode, development / validation means: the
library's base 0.989053 / 0.993816 (scratch 0.989053 / 0.993810), with rules 1 to 3 0.993692 /
0.996434 (0.993689 / 0.996426), with rule 4 and the clause as it was 0.994604 / 0.996853 (0.994600 /
0.996845). The library's tones, channel-error and pulse-train decodes were not compared with the
scratch build's.

**Speech**, decoding DVSI's bits, on the scoreboard (P25; the gate sets' measures, partition means,
the mean paired change and its extremes):

| Partition | Waveform corr. vs DVSI | LSD vs DVSI (dB) | Envelope corr. | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| development | 0.989085 → 0.994727 (+0.005642; 20 up, the least `t01` +0.002228, the most `dam10` +0.013283) | 1.2510 → 1.0608 (-0.1902; all fall, the least `dam80` -0.0221, the most `lvl` -0.3258) | 0.99549 → 0.99565 (+0.00017; `dam40` -0.000007) | +0.0522 → +0.0549 (+0.0027) | 6.3446 → 6.3501 (+0.0055; 17 rise, the most `dam10` +0.0222) |
| validation | 0.993818 → 0.996757 (+0.002939; 5 up, the least `t04` +0.001747) | 0.8783 → 0.7251 (-0.1532; all fall, the least `t03` -0.0870) | 0.99340 → 0.99347 | +0.0298 → +0.0324 | 4.7763 → 4.7745 (-0.0018; 4 fall, `t03` +0.0074) |

The scratch build's run gave +0.005548 / +0.003036 of correlation and -0.1491 / -0.1317 dB of
distance from DVSI's decode: the library's clause on the column-read states lowers the distance by a
further 0.04 / 0.02 dB (above). **The distance from the input rises on development**, by 0.0055 dB in
the mean (17 of 20 vectors, the most `dam10`, +0.022): the own models' column reading voices
harmonics that the analysis of the input left unvoiced, and a decode that follows DVSI's there moves
away from the input. It is this change's cost. The whole-file correlation at lag 0 (the study's
measure), the rules in turn: 0.989053 / 0.993816 → the grid 0.965205 / 0.972297 → the parabola
0.991623 / 0.995199 → the phase step 0.993692 / 0.996434 → the own models 0.994695 / 0.996755; each
rule alone on the baseline: the grid -0.023848 / -0.021519, the parabola -0.033616 / -0.025837, the
phase step +0.002747 / +0.001540, the own models +0.000890 / +0.000267 (validation 4 up, `tambf32b`
-0.000397). Of the residual over DVSI's energy pooled over the 22 vectors with a DVSI rms of 300 or
more (the study's `shares10.py`), 0.01468 → 0.00518: in voiced frames in line with DVSI's 0.01118 →
0.00201 (0.0115 → 0.0021 of their own signal), in frames with no voiced band 0.00069 → 0.00041,
and the 59 to 71 voiced frames out of line (0.0028) do not move. In steady voiced frames in line
(the study's `resid10.py`), residual over DVSI's energy 0.01043 = amplitude 0.00062 + phase 0.00796
+ incoherent 0.00185 → 0.00168 = 0.00046 + 0.00079 + 0.00043 (after the grid and the parabola
0.00677, phase 0.00433; with the phase step 0.00378, phase 0.00153); the voiced harmonics' phase
against DVSI's, energy-weighted rms, 0.0593 rad after the pitch rules → 0.0270, their level 0.167 →
0.164 dB rms. The onsets' lead (the study's `p25on10b.py`, 1,854 onsets followed by two voiced
frames: DVSI's decode ahead of ours by up to 0.35 samples in the windows 160 and 240 samples into
the onset frame at short periods) is gone with the pitch rules and the phase step: within about
0.05 samples in every period range with the first half sounding. About half of P25's remaining
whole-file residual lies in a dozen voiced runs out of line with DVSI's decode (`noisy` from frames
2796 and 1946, `t04` from 684, ...), in the ones examined after a speech frame DVSI's decoder plays
as a tone model at a pitch code one or two away from a tone's own, which the library's tone test
does not pass: a later study (Known differences); `noisy` and `t04` gain least on validation and
development for it.

PESQ (the 25 vectors; a 95% interval over the 12 recordings, t with 11 degrees of freedom, from the
standard error clustered by recording): against the input +0.0029 ± 0.0026 (-0.0027 to +0.0086, 16
of 25 up), neither a cost nor a gain established; against DVSI's decode +0.0468 ± 0.0049 (+0.0360 to
+0.0575, 25 of 25 up).

**Tones.** All 52 P25 tone decodes change (the D-STAR and rate-33 ones do not). The detector's
agreement moves on one vector, `cp31` (a call-progress vector P25 plays as decoded speech: one
steady window of 37 lost, 0.973 → 0.946); the level error 0.0967 → 0.0953 dB in the mean over the 49
vectors that measure it (`cp31` -0.044, `cp2` -0.024), the balance error 0.0158 → 0.0161 dB over 48
(`cp2` +0.041, `cp31` -0.028), the frequency error 0.00023 → 0.00024 over 49 (`cp31` +0.00018), the
share of wrong tones 0 on every vector. The DTMF vectors whose tones are separated by short voiced
runs lose whole-file correlation ([P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)).
The encoders' tone statistics do not move.

**Channel errors.** No mute decision changes, frame for frame: the flags of every frame of the 103
P25 bit files (`irstia` and `mark` left out; the soft-decision files decoded as soft; 105,033
frames) are the same in both builds; every decode changes. The scoreboard's log-spectral distance
from DVSI's decode on the loud frames that neither run nor DVSI mutes falls: hard decisions -0.0559
dB in the mean (every vector, the least `dam_e10_hd` -0.0164, coverage 0.974), soft decisions
-0.0400 (coverage 0.971); `channel.hd.lsd_db` 4.7868 → 4.7309, `channel.sd.lsd_db` 5.6805 →
5.6406. Whole-file correlation with DVSI's decode on the hard-decision vectors (the study's
`parts8.py`): 0.8056 → 0.8098, 0.4554 → 0.4657, 0.0340 → 0.0502, 0.0080 → 0.0121.

**Pulse-train vectors** ([P25 pitch](#p25-pitch)): the lag increments fall (0.0185 → 0.0163 samples),
the whole-file correlation rises on 8 of 13 and falls on 5.

**Other modes and the encoders.** All 285 encodings (95 input files, three encoders) are identical
between the two builds; the IMBE encoder's `cur_mp` differs in its `w0` only, which the public
decoder sets. In each of the ten local configurations (dev-debug; dev-release, GCC Release with
SIMD, fast-math, LTO and `-march=native`; asan-ubsan-debug; notones-debug; a Debug build with SIMD;
a Debug `-march=x86-64-v3` build; Clang Release with SIMD, fast-math and LTO, once with
`MBELIB_ENABLE_NATIVE=ON` and once with `-march=x86-64-v3`; and 32-bit Release and Debug), after
each of the four rule commits, a differential harness outside the repository, built against the
baseline's library and the commit's, gives the same output byte for byte for: the pitch and harmonic
count every D-STAR code decodes to and every rate-33 code at rows 16 and 18; hashes of the output
and of all three returned states of the D-STAR and rate-33 data paths (3,000 random frames each) and
of ProVoice's path (its float and 16-bit frame entry points, 3,000 frames; and 500 frames each from
the states the D-STAR and rate-33 streams left), of `mbe_synthesizeSpeechf()` and
`mbe_synthesizeSpeech()` (3,000 frames), of the two AMBE encoders' bits and `cur_mp` and of the IMBE
encoder's bits (1,500 frames), over constructed streams whose caller states hold reserved pitches,
P25 pitches of every code in three forms (the earlier law's, the grid's and eq 46's nominal, from
the harness's own formulas, so that both builds see the same inputs) with voiced bands (a low voiced
region, alternate bands, every band, a lone band) and unvoiced ones and shaped `log2Ml`, every frame
mark, pulsed-row words, values in slots 34 to 37, every float class in the low cut's slots, and
states without the subframe mark; and ProVoice plays the same samples and returns the same states
whatever the low cut's slots hold. P25's stream, ProVoice from P25's states and the IMBE encoder's
`cur_mp` change, as they are to. A first form of rules 3 and 4, before their final commits, moved
D-STAR's and rate 33's streams in the 32-bit Release build: P25's mid-frame filler and column state,
written as branches of the shared ones, gave a column helper a second caller, which GCC then no
longer inlined, and the shared code spilled its x87 intermediates otherwise; P25's rendering has its
own copies of the three (`p25_fill_mid_bands()`, `p25_column_band_state()`,
`p25_columns_voice()`), and the ten configurations agree since. No golden PCM sequence is P25's: no
golden hash moves.

**Cost.** The median complete hard-decision frame decode, `bench_replay hard <codec> 60 5` on DVSI's
`dam` stream of each codec, seven alternating invocations pinned to one CPU, a portable Release build
with SIMD, fast-math and LTO, taken under load (one-minute load 9.4 to 12.6, another agent's builds
and tests running): IMBE `f725a66` 8.558 us, the grid 8.358 (-2.3%), the parabola 9.191 (+7.4%), the
phase step 8.213 (-4.0%), the own models 8.698 (+1.6%; the invocations' medians 8.117 to 8.777 and
8.457 to 9.293); the controls, whose decodes do not change, moved by as much: D-STAR -4.9% to +6.0%,
AMBE+2 -0.7% to +1.5%. Under that load the runs give no bound on the change's cost tighter than about
6%.

**Gate sets**, declared in the branch's first commit, `f1be05d`, before any build of the branch was
scored (the study's scratch build had been measured on the development, validation, tone and channel
vectors before that, and the baseline's own scoreboard of the four partitions, a build of `git
archive f725a66` scoring itself, was read for each measure's coverage; two checks of the brief the
sets were written from could not be evaluated as written and were corrected before the declaration:
the IMBE encoder's distance from the input is `encoder.lsd_db`, and the channel vectors, one
recording, give no bootstrap interval, so their hard-decision distance is checked on the mean and on
each vector), and evaluated in one combined run over the development, validation, tone and channel
partitions, `f725a66` against `999647e` (above); all 75 speech entries carry an encoder-row hash in both
reports, every one equal:

| Gate set | Result | What it establishes |
|---|---|---|
| `p25-decoder-development` | passes all 15 checks | on the vectors the rules were found on: a mean gain in P25's correlation with DVSI's decode with its interval above 0 (+0.005642), every vector up (the least +0.002228), the distance from DVSI's decode falling with its interval below 0 (-0.1902 dB) and on every vector (the most -0.0221), the envelope correlation, the distance from the input (+0.0055 dB in the mean against 0.012 allowed, +0.022 on `dam10` against 0.04) and the level within theirs; D-STAR's and rate 33's decodes and every encoder's rows unchanged |
| `p25-decoder-validation` | passes all 15 checks | the same on the validation vectors (+0.002939, the least +0.001747; -0.1532 dB), also read in the study |
| `p25-decoder-tones` | passes all 13 checks | D-STAR's and rate 33's decodes and tone statistics unchanged, P25's tone measures within their limits on the windows the detector calls steady (`cp31`'s one window inside the 0.03 allowed); not that the frames between and around the tones are right |
| `p25-decoder-channel` | passes all 33 checks | the per-vector mute statistics unchanged in every mode (the flags say no frame's decision changed), D-STAR's and rate 33's channel decodes unchanged, P25's distance from DVSI's decode on the loud frames falling on the mean and on every hard-decision vector and not rising with soft decisions |
| `p25-decoder-final` | not scored | the prediction for `irstia` and `mark` |

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

The November 2026 synthesis changes were each scored once against the commit before them, with
two exceptions. The voiced synthesis was scored at a voiced gain of +0.45 dB first, failed its level
check, and was scored again after the gain was set on development frame medians ([Synthesis in
subframes](#synthesis-in-subframes)). And the tone partition was scored again for each change
after a late correction to D-STAR's unvoiced overlap across a tone, which changes seven D-STAR
tone vectors and no other stream: every check came out as before. The D-STAR checks of
`ambe2450-tone-fade-tones` were tightened after its first scoring, from the mean level and the
detector agreement to no change on any vector in the five metrics its description names, and
pass. The final partition was not held out for these changes.

The three changes of 2.6.12 were each scored once against the commit before them, on the
development, validation, tone and channel partitions; `joint-search-development` failed and was kept
as declared. The final partition was then scored once for the three together, against a gate set
declared before it ([November 2026: encoders and rate-33 rows 18 to
31](#november-2026-encoders-and-rate-33-rows-18-to-31)), and passed. After that, each change was
scored a second time on the four other partitions, with the decode counts below.

`r33-row-grid-tones` and `r33-row-grid-channel` state that the D-STAR and P25 decodes of DVSI's bits
are unchanged, but checked it through `speech.pcm_changed_vectors`, a count over the speech vectors
that the scoreboard copied into every group, so the check did not cover the tone and channel-error
vectors those sets judge. The scoreboard now hashes our decode of every tone and channel-error
vector and counts the changed ones (`tone.pcm_changed_vectors`, `channel.pcm_changed_vectors`;
[DVSI scoreboard](testing.md#dvsi-scoreboard)). After their first scoring the two sets were given
those counts, at most 0 for D-STAR and P25, and the band-grid change was scored again on the
development, validation, tone and channel partitions, with the same builds and protocol: the
development and validation results are the same to the last digit, every check passes as before,
and none of the D-STAR and P25 decodes of the 52 tone and 8 channel-error vectors changes (the
grid changes the rate-33 decodes of 9 tone vectors and of all 8 channel-error vectors).

`unvoiced-code-tones`, `unvoiced-code-channel`, `joint-search-tones` and `joint-search-channel`
state that every decode of DVSI's bits is unchanged, but checked it through tone, mute and distance
metrics, which a changed decode can leave as they were: on the band-grid change's runs both tone
sets pass, with nine rate-33 tone decodes changed. After their first scoring the four sets were
given the same counts, at most 0 in all three modes, and both changes were scored again on the
development, validation, tone and channel partitions, with the same builds and protocol: the
development and validation results are the same to the last digit (`joint-search-development`
failing as before: the distance from the input rises by 0.0105 and 0.0123 dB, against a limit of
0.01), every other check passes as before, and for either change none of the decodes of the 52 tone
and 8 channel-error vectors changes in any mode.

The three decoder changes of 2.6.13 were each scored against the commit before them, on the
development, validation, tone and channel partitions, against gate sets declared before any of them
was scored. The tone and channel partitions were scored twice. First each was scored alone, in a run
of its own; there the `speech.pcm_changed_vectors` checks of `p25-pitch-tones`, `-channel`,
`r33-mid-pitch-tones` and `-channel` find no speech vector to count and fail for want of data. Then,
as the protocol has it ([DVSI scoreboard](testing.md#dvsi-scoreboard)), each step was scored in a
combined run, the four partitions in one run and its four gate sets together, with the same builds:
those checks pass, and every other check and figure is the same as when the tone and channel
partitions were scored alone. The outcomes reported above are the combined run's. The final
partition was scored once, for the three changes together, after they and the review's fixes were
complete, and passed ([November 2026: phase regeneration and
pitch](#november-2026-phase-regeneration-and-pitch)).

The three changes of the tone model change set were scored together, `fac0111` against the change
set, on the development, validation, tone and channel partitions in one combined run with the four
gate sets
declared in `c736edf`, before any build of the change set was scored. The first build of the P25
tone path, `d245efd`, failed `tone-model-development` and `-validation` on P25; its single-tone
limit was corrected (`b58d91d`) and the change set scored again with the same gate sets, builds of
the baseline and protocol. `tone-model-tones` failed the same five checks both times and is kept as
declared ([November 2026: tones and the octave step](#november-2026-tones-and-the-octave-step)).
The final partition was scored once for the tone model change, on 2026-11-26, after the review:
`fac0111` against
`c2ba220`, the reviewed code, against `tone-model-final`, declared in `c736edf`. It passes all 24
checks ([November 2026: tones and the octave step](#november-2026-tones-and-the-octave-step)).

The two decoder changes of the mid-frame voicing change set were scored together, on the
development, validation, tone and channel partitions in one combined run with the four gate sets
declared in `6acba79`, before any build of the change set was scored: first `cc272ef` against the
change set on it, then, after the change set was moved onto the tone branch's final code, `c2ba220`
against it, with the same partition means and vector hashes, and a third time after the review's
fixes, `adecb29` (the same code as `c2ba220`) against the fixed change set, where only P25's
soft-decision distance on the loud frames of the channel set moved (-0.0106 to -0.0109 dB); all
four pass each time ([November 2026: mid-frame voicing](#november-2026-mid-frame-voicing)). The
final partition was scored once for the mid-frame voicing change, on 2026-11-26, after the review:
the tone model change's code against `504740c`, the reviewed code and documents, against
`mid-voicing-final`, declared in `6acba79`. It passes all 18 checks.

Rate 33's mid-frame pitch wave was scored once, `dac6d14` against `d08ea39`, on the development,
validation, tone and channel partitions in one combined run with the four gate sets declared in
`8942924`, the branch's first commit, before any build of the branch was scored; all four pass
([November 2026: rate 33's mid-frame pitch wave](#november-2026-rate-33s-mid-frame-pitch-wave)).
Before that declaration the study's scratch build had been measured on the development and
validation vectors, and the baseline alone was scored on the tone partition, to record each tone
measure's coverage in the tones set's description. The identity counters of the four sets were
evaluated on the combined run: a tone or channel partition scored alone leaves some of them without
a value. The final partition was scored once for it, on 2026-11-26, after the review: `dac6d14`
against the head's code, against `r33-pitch-wave-final`, declared in `8942924`. It passes all 12
checks ([November 2026: rate 33's mid-frame pitch wave](#november-2026-rate-33s-mid-frame-pitch-wave)).

The mid-frame columns change was scored first `c3c6d42` against `e9cf77c`, on the development,
validation, tone and channel partitions in one combined run with the four gate sets declared in
`3593e8c`, the branch's first commit, before any build of the branch was scored; all four pass
([November 2026: mid-frame voicing by columns](#november-2026-mid-frame-voicing-by-columns)). After the
review's two fixes (`e86afd0`, which change three P25 tone decodes) it was scored a second time,
`c3c6d42` against `e86afd0`, with the same sets, builds of the baseline and protocol: all four pass
again, every figure but those of the three tone vectors as before. Before
that declaration the study's scratch build had been measured on the development and validation
vectors, and the baseline's own scoreboard of the four partitions (from the rate-33 pitch wave's
run) was read to record each measure's coverage in the sets' descriptions. The final partition was
scored once for it, on 2026-11-26, after the review and its two rounds of fixes: `c3c6d42` against
the head's code (the build of `e86afd0`), against `mid-columns-final`, declared in `3593e8c`. It
passes all 18 checks ([November 2026: mid-frame voicing by
columns](#november-2026-mid-frame-voicing-by-columns)).

The row-17 train change was scored once, `e446950` against the head's code (the build of `2a3a64f`,
whose library code is `8cbef45`'s), on the development, validation, tone and channel partitions in
one combined run with the four gate sets declared in `6a32fb8`, the branch's first commit, before
any build of the branch was scored; all four pass, 72 of 72 checks ([November 2026: the pulse train
through a row-17 frame](#november-2026-the-pulse-train-through-a-row-17-frame)). Before that
declaration the study's scratch build had been measured on the development, validation and rate-33
tone vectors, and the baseline alone was scored on the four partitions (its own encodings decoded
by its own decoder) to record each measure's coverage in the sets' descriptions. The final partition
was scored once for it, on 2026-11-26, after the review and its fixes: `e446950` against the head's
code (the build of `3840d2b`), against `r33-row17-train-final`, declared in `6a32fb8`. It passes all
12 checks; no decode of `irstia` or `mark` changes in any mode, so it shows no harm and cannot show
the gain ([November 2026: the pulse train through a row-17
frame](#november-2026-the-pulse-train-through-a-row-17-frame)).

The pulsed state was scored once, `74126a3` against the head's code (the build of `d3fb70f`, whose
library code is `834d001`'s), on the development, validation, tone and channel partitions in one
combined run with the four gate sets declared in `525d017`, the branch's first commit, before any
build of the branch was scored; all four pass, 87 of 87 checks ([November 2026: the pulsed
state](#november-2026-the-pulsed-state)). Before that declaration the study's scratch build had been
measured on the development, validation, tone, hard-decision channel-error and D-STAR DTX vectors,
and the baseline alone was scored on the four partitions (a build of `git archive 74126a3`, its own
encodings decoded by its own decoder) to record each measure's coverage in the sets' descriptions.
After the review the same four sets were scored again in one
combined run, `74126a3` against the final code, the validity counter in: 87 of 87 checks, with the
same figures ([November 2026: the pulsed state](#november-2026-the-pulsed-state)). The final
partition was scored once for it, on 2026-11-27, after the two reviews and their fixes: the build of
`3840d2b` (library code `74126a3`'s) against the build of `5cfdf37` (library code `43850cc`'s),
against `pulsed-state-final`, declared in `525d017`. It passes all 18 checks; the gain it required is
there on both vectors in both codecs ([the final partition](#the-final-partition)).

Rate 33's pulsed rows were scored once, `2b6b6c2` against the branch's library code at `a1398df`, on
the development, validation, tone and channel partitions in one combined run with the four gate sets
declared in `5af6391`, the branch's first commit, before any build of the branch was scored; all four
pass, 72 of 72 checks ([November 2026: rate 33's pulsed rows](#november-2026-rate-33s-pulsed-rows)).
Before that declaration the study's scratch build had been measured on the development, validation,
tone, channel (hard and soft decisions) and rate-33 DTX vectors, with the scoreboard's own measures
for rate 33 too, and the baseline alone was scored on the four partitions (a build of `git archive
2b6b6c2`, its own encodings decoded by its own decoder) to record each measure's coverage in the
sets' descriptions. The final partition was scored once, on 2026-11-27, after the two reviews and
their fixes, the build of `2b6b6c2` against the build of `1502dc0`: `r33-pulsed-rows-final` passed,
12 of 12 checks, rate 33's correlation rising on both held-out vectors (mean 0.985957 → 0.989556),
D-STAR's and P25's decodes and every encoder row unchanged ([November 2026: rate 33's pulsed
rows](#november-2026-rate-33s-pulsed-rows)).

DVSI's low band was scored once, `9a101f8` against the branch's library code at `0e3463e`, on the
development, validation, tone and channel partitions in one combined run with the four gate sets
declared in `4677a41`, the branch's first commit, before any build of the branch was scored; all four
pass, 102 of 102 checks ([November 2026: DVSI's low band](#november-2026-dvsis-low-band)). Before that
declaration the study's scratch build had been measured on the development, validation, tone, channel
and pulse-train vectors, with the scoreboard's own measures too, and the baseline alone was scored on
the four partitions (a build of `git archive 9a101f8`, its own encodings decoded by its own decoder) to
record each measure's coverage in the sets' descriptions. The final partition is not yet scored.
