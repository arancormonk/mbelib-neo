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
  pitch](#november-2026-phase-regeneration-and-pitch)). The three decoder changes of 2.6.14 ([The
  pulse train across an octave](#the-pulse-train-across-an-octave), [Tones as harmonics of a
  fundamental](#tones-as-harmonics-of-a-fundamental) and [P25 tones sent as voice
  frames](#p25-tones-sent-as-voice-frames)) hold it out as well: the study that found them read
  DVSI's tone vectors and the development and validation vectors, never `irstia` or `mark`; the
  change set was measured and gated on the development, validation, tone and channel partitions,
  and the final partition was scored once for it, after the review, against the gate set
  `tone-model-final`, declared before any build was scored, which it passed ([November 2026: tones
  and the octave step](#november-2026-tones-and-the-octave-step)).
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
5) is voiced or when the column 250 Hz above it is: each voiced column reaches half a column
(250 Hz) into an unvoiced one below it, as D-STAR's 1 kHz bands reach 250 Hz into theirs ([D-STAR voicing bands](#d-star-voicing-bands)).
It also reads row 17 with its lowest column voiced, where Annex B lists rows 17 to 31 all
unvoiced like row 16. The decoder, and the encoder's eqs 4 and 8, now read the voicing that way;
the encoder still searches rows 0 to 16 (clause 4.2), so row 17 only changes how DVSI's bits play.
On DVSI's speech vectors our decode of frames coded 17 comes 0.5 dB of log-spectral distance closer
to DVSI's (27 frames), and of voiced frames generally 0.03 dB.

Rows 1, 3, 13, 15 and 17 play some columns neither voiced nor unvoiced on DVSI's decoder (see
[Known differences](#known-differences)); voiced is the nearer binary value for each (log-spectral
distance against DVSI's decode on frames coded 13 or 15: 6.87 dB voiced, 7.56 dB unvoiced), so the
table keeps Annex B's voiced columns there.

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
on 17 of the 25 vectors; this ± is over the vectors, [Method](#method)). Row 17 keeps its lowest
column voiced.

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
31](#november-2026-encoders-and-rate-33-rows-18-to-31)). What these rows do beyond the grid is not
modelled ([Known differences](#known-differences)).

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
DVSI's noise generator runs as for any frame. The rate-33 path uses it for the first valid voice or
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

### Output high-pass phase

DVSI's decoders, in every mode, roll off the band below 200 Hz like a second-order high-pass
near 100 Hz, and this library already rendered the magnitude of a second-order Butterworth
high-pass at 100 Hz for it. It now renders that filter's phase too, on every voiced harmonic of
both frames being interpolated, at render time, so the model, its regenerated phase and the
stored state keep the decoded values. On the validation vectors the phase DVSI renders differs
from ours by a median 8.9, 6.6 and 7.7 degrees (D-STAR, P25, rate 33), 4-5 below 300 Hz, where it
differed by 12.1, 10.0 and 11.5 (12-14 below 300 Hz): per harmonic of steady voiced frames, once
each frame's time shift and common phase are removed. The filter's corner and shape come from the
magnitude already measured; nothing in it was fitted to the phase. A rate-33 or D-STAR tone that
starts after a frame that was not a tone frame starts at this filter's phase at each component's
frequency, as DVSI's decoders start a tone after a frame with no voiced band; after a voiced frame
they continue that frame's pulse train into the tone, which the library does not model ([Tones as
harmonics of a fundamental](#tones-as-harmonics-of-a-fundamental)).

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
  unvoiced render gain; the output high-pass acts on every bin, in magnitude and in phase.

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
  between two voiced frames of different pitch, a little below it: [Rate 33: mid-frame
  pitch](#rate-33-mid-frame-pitch)), its harmonic count that pitch's, each band's voicing this
  frame's at the band's frequency, its
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
  model's. A subframe with no voiced band at either end starts the train afresh.

DVSI describes a synthesis of this shape in US 8,200,497: square-root triangular windows over
10 ms subframes, and an impulse sequence convolved with a pitch-dependent impulse response. The
rules and constants here were measured on DVSI's decoders, from its vectors and, for D-STAR, on
the hardware below.

The codec paths render this, with a voiced gain of +0.55 dB: fully voiced frames then play
within 0.02 dB of DVSI's (frame medians on the development vectors, all three modes; 0.09-0.12 dB
low at the +0.45 dB of the rendering before). The pulse phase and the pulses still sounding
travel in `noiseOverlap`, which these paths no longer use for noise, and `PSIl` and `PHIl` hold
the harmonic phases at the frame's end as before. The public `mbe_synthesizeSpeechf()` keeps
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

The single-band steps come out worse: on the unit a band above 1 kHz that turns voiced starts
with the frame and one that turns unvoiced lasts to mid-frame, while the lowest band does the
opposite on both; the mid-frame model here takes this frame's voicing for every band.

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
  D-STAR b0 64..119). Rate 33 has its V/UV rows 16 to 31 for this.
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
  reaching 507.9 Hz do. Rate 33's row 16 ends the voicing where the frame starts.
- **D-STAR b0 64..71, 88..103 and 112..119 carry their pitch** into a voiced frame after them:
  its first subframe runs at the geometric mean of the two pitches (253 onsets of 253 on DVSI's
  vectors), and at the voiced frame's own pitch after any other unvoiced frame (925 of 925).
- **Late onsets.** A rate-33 frame coded with V/UV row 15 after an unvoiced one starts its
  voiced bands in its second subframe, and so does a P25 frame whose voicing does not reach
  505 Hz (all 40 such onsets on DVSI's vectors).

With these, voiced runs after two or more unvoiced frames start where DVSI's do in 95% of cases for
P25 (791 of 831 on ten vectors), 96% for D-STAR (812 of 846) and 96% for rate 33 (775 of 810).
Measured again for the octave step of 2.6.14 on all 25 development and validation vectors (D-STAR
voiced runs of three voiced frames or more after a frame with no voiced band, loud; out of line
where the two frames after the first correlate with DVSI's under 0.8), 25 of 1,945 D-STAR runs were
out of line before that step and 4 of 1,945 are with it: most of the misses after a pitch-carrying
code were a step of about an octave into the second subframe ([The pulse train across an
octave](#the-pulse-train-across-an-octave)). Not resolved: P25 onsets that DVSI's decoder starts in
the frame's second subframe and the library starts with the frame. Measured for 2.6.14, 29 of the
1,940 P25 onsets of the development and validation vectors are such onsets; in 28 of them the
harmonic at or next above 505 Hz is unvoiced, but so it is in 42 onsets that DVSI's decoder starts
with the frame, and no rule is adopted (the measurement is in [Known
differences](#known-differences)). See also [Known differences](#known-differences) for what the
unvoiced pitch codes and rate 33's rows 17 to 31 further tell DVSI's decoder.

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

On speech, DVSI's P25 decoder does not play pitch code b0 at TIA-102.BABA eq 46's fundamental
`w0 = 4 pi / (b0 + 39.5)` exactly: it plays that value times `1 + g`, with g a smooth function of
the nominal pitch period `P = (b0 + 39.5) / 2` that repeats every half octave and lies between
-0.071% and +0.034%. A tone its encoder sent as voice frames it plays at eq 46's fundamental cut
to a grid of 2^-19 cycle per sample ([P25 tones sent as voice
frames](#p25-tones-sent-as-voice-frames)).

**The drift.** With [Phase regeneration](#phase-regeneration) in place, the residual of P25's
voiced runs that line up with DVSI's grew with the run's length (development vectors, residual
energy over DVSI's): 0.023 of the signal over a run's first four frames, 0.043 over the next four,
then 0.085, 0.098, 0.111, and 0.21 past the twentieth frame. D-STAR's was 0.016 to 0.018 over a
run's first twenty frames and 0.037 past them.

**Per pitch code.** On two consecutive frames at one pitch code, measured only where four
consecutive frames have that code (each frame of the pair needs the two frames before it at its
code), the phase of DVSI's output against ours on each voiced harmonic advances by `l w0` times the
time by which DVSI's pulse train gains on ours; divided by the frame length that is DVSI's
fundamental against ours. It is the same on every harmonic (b0 40: -0.086 to -0.096 samples a frame
on harmonics 1 to 18; b0 80, on the pulse-train vectors: -0.100 to -0.113 on harmonics 1 to 27; b0
29: +0.049 to +0.060 on harmonics 1 to 15), so it is a pitch difference. The measurement keeps one
row per pitch code and corpus: 58 rows from DVSI's synthetic pulse-train vectors (`ptrain0` to
`ptrain12`) and 24 from the speech vectors of the development partition, 5,453 frame pairs. Over the
39 rows with 10 pairs or more DVSI's fundamental runs from -0.069% (b0 42) to +0.036% (b0 57) of eq
46's (a row with 5 pairs, b0 193, reads -0.106%), as a function of P that repeats every half octave:
zero at P = 32, 45.3 and 64 (2^(k/2)), about +0.036% a fifth of the way through each half octave,
zero again at 0.39 of it and about -0.067% at 0.72. Examples (b0: DVSI's fundamental against eq
46's, frame pairs): 21: -0.063% (72), 29: +0.036% (35), 30: +0.031% (74), 40: -0.058% (80), 51:
-0.004% (80), 61: +0.011% (70), 71: -0.045% (74), 80: -0.064% (2931), 81: -0.068% (227), 165:
-0.010% (609). The best single sinusoid in log2 f0 has exactly 2.00 cycles an octave.

**The rule.** `w0 = 4 pi / (b0 + 39.5) (1 + g)`, with t = 4 pi log2 P and
`g = a0 + a1 cos t + b1 sin t + a2 cos 2t + b2 sin 2t` (`imbe_dvsi_pitch_offset()`: a0 = -1.980564e-4,
a1 = 1.263798e-4, b1 = 4.835086e-4, a2 = 2.352386e-5, b2 = 7.623621e-5), fitted by least squares
over the rows with five frame pairs or more: 69 rows (50 from the pulse-train vectors, 19 from
speech) holding 5,401 frame pairs on 61 different pitch codes (eight codes are in both corpora: 30,
40, 46, 54, 57, 63, 67 and 79), each row's weight capped at 40 pairs, so a code in both corpora can
weigh more (b0 40 has 57). It leaves 0.0064% rms with those weights, at most 0.014% on a row with 10
or more pairs; a third harmonic changes nothing end to end. On the validation vectors' own codes (9
codes, 115 pairs, noisy) the rms difference from DVSI's, weighted by pairs, falls from 0.046% to
0.018%. The nominal pitch of eq 46 is computed in float as the decoder always has (rounded to float
before the correction also where the arithmetic keeps excess precision, as on 32-bit x87), and L, K
and so the bit layout follow it as before. The pitch stored for synthesis is that float times
`1 + g`, g computed in double with the double-precision libm functions (from -7.06e-4 at b0 127 to
+3.38e-4 at b0 96), the product formed in double and rounded to float once. That is not the law
evaluated wholly in double, from which it differs by an ulp at 79 of the 208 codes. For all 208
codes and any harmonic as the top voiced one, neither the late-onset mark the decoder sets ([Voiced
onsets and offsets](#voiced-onsets-and-offsets)) nor the harmonic whose voicing the synthesis reads
to tell whether a frame's voicing reaches 505 Hz (`ceil(505 Hz / f0)`) moves under the correction:
the nearest harmonic to 505 Hz is b0 119's 5th, 0.053% below it, and the correction moves it further
below. Other discrete steps of synthesis follow the pitch played and do move at some codes (a band
edge's transform bin, the mid-frame model's harmonic count, the harmonic where the high-band shelf
starts). ProVoice, which reaches the same data path, has no DVSI reference and keeps eq 46; the
public `mbe_decodeImbe4400Parms()`, which has no frame context, decodes P25's. The P25 encoder
quantizes from its own pitch code, so its bits do not change.

With it, the aligned runs no longer grow: 0.023, 0.019, 0.023, 0.028, 0.024 and 0.024 over frames 0
to 3, 4 to 7 and so on, 20 on. Where the curve comes from is not established; a low-order
approximation of a logarithm or an exponential over half-octave segments in the decoder's
fixed-point arithmetic would give one of this shape.

**Tones sent as voice frames.** DVSI's P25 decoder does not play a tone its encoder sent as voice
frames at the pitch it plays speech at. On DVSI's 52 tone vectors, over pairs of consecutive voiced
frames at one pitch code (measured where four consecutive frames have that code, as above), the
phase of DVSI's output against the phase-change build's, which plays eq 46's pitch, measured per
harmonic (the strong voiced harmonics, within 10 dB of the frame's strongest, coherence 0.98 or
more), moves by a median over the frame pairs of -44 to +37 mrad a frame for each of the 61
combinations of pitch code and harmonic that have 60 frame pairs or more (medians, not a bound on
single frame pairs or on harmonics not measured), not in proportion to the harmonic number and with
either sign on one code: b0 36, harmonics 4 and 7, +7 and +37 mrad a frame; b0 92, harmonics 7 and
11, +16 and -11; b0 144, harmonics 8 and 17, -15 and +16; b0 189, harmonics 10, 11, 19 and 21, -21,
-39, -19 and -44; b0 127, harmonics 8 and 17, -32 and -14, where the law gives -68 and -145. For
2.6.13 this was read as each measured harmonic playing at eq 46's frequency to within 0.35 Hz (44
mrad a frame) in the median over its frame pairs, with a small offset of its own that follows no
rule found. That reading is corrected: those offsets came from the short runs of the DTMF vectors
and the sidebands of their frame-to-frame amplitude changes. The long runs of the KNOX vectors,
measured since, are clean: DVSI's decoder plays such a tone at eq 46's fundamental cut to a grid of
2^-19 cycle per sample, both harmonics of a tone at the same relative offset from eq 46 ([P25 tones
sent as voice frames](#p25-tones-sent-as-voice-frames)). Fitted
as one time shift per frame pair to the phases of the voiced harmonics whose coherence is 0.9 or
more (at least two of them, DVSI's frame standard deviation 100 or more: 11,647 frame pairs, of
which the 11,535 with a fitted shift under 0.5 samples are used; 28 codes with 100 pairs or more),
these offsets come out at a median of -0.0006 (b0 127) to -0.041 samples a frame (b0 144) on a code,
whatever the law gives for the code (-0.113 to +0.054 samples a frame on those codes), and the slope
of that fit on the law's value is 0.12 over those 11,535 pairs. The same fit on the development
speech vectors (298 pairs; three frames at one code are rare in speech) gives a slope of 0.92, and
1.03 on the 105 pairs with every band voiced at codes where the law moves the train by more than
0.04 samples a frame (rms difference from the law 0.037 samples a frame); the per-code figures of
"Per pitch code" above came from speech and DVSI's pulse-train vectors (5,453 pairs). At b0 80
DVSI's pulse-train vectors drift by -0.10 samples a frame, its tone vectors by -0.010 in the
one-time-shift fit.

For 2.6.13 what tells DVSI's decoder the two apart was not found. Tried: the tone frames' voicing
patterns vary (the bands that
hold the tone's harmonics are voiced: `01100`, `1110000`, `0111110010000`, ...; 4 of the 28 codes
have the lowest band voiced, as speech has, and do not drift either); the share of a frame's decoded
energy outside its two strongest harmonics does not separate them (tone frames: median 1.1%, 90th
percentile 26%; speech frames: 1st percentile 0.17%, 10th percentile 2.4%; 16,543 and 25,081 voiced
frames above an RMS of 100); no floor or rounding of a harmonic's frequency `l / P`, of the
fundamental's `1 / P` times l, or of `w0` times l, to a grid of 2^-10 to 2^-25 (cycles a sample, or
radians) explains the per-harmonic offsets (the best leaves 77% of their rms); and for the
one-time-shift fit, no floor, rounding or ceiling of w0, f0, the period or log2 f0 to a grid of 2^-8
to 2^-23 does (the best leaves 63%). Those searches ran on the contaminated offsets above; on the
long runs the floor of the fundamental to 2^-19 cycle per sample does explain the pitch. What tells
DVSI's decoder a tone from speech was found later: the pitch code and the voicing word its encoder
sends for the tone, with the tone's harmonics dominating the others ([P25 tones sent as voice
frames](#p25-tones-sent-as-voice-frames)); the share of energy alone does not separate them.

2.6.13 applied the pitch of this section to every P25 frame. On speech that is a fit to the pitch
DVSI's decoder plays (above). On a tone sent as voice frames it moved harmonic l by `l w0 160 g` a
frame against eq 46's: by up to 145 mrad a frame over the 61 measured combinations of pitch code and
harmonic (b0 127's 17th harmonic, -145), and by -307 mrad a frame on the ninth harmonic of b0 2,
which `alltone` holds for a few frames. So the median correlation of the steady tone frames with
DVSI's waveform fell from 0.880 to 0.656 (2.6.13's selection: the 17,334 loud frames of the 52
vectors whose neighbours are loud too), while their level error moved from 0.555 to 0.560 dB and
the scoreboard's tone detector's agreement with DVSI's output moved by at most 0.0083 on a vector
(`knox_12`) (below). That cost belonged to the 2.6.13 change and is gone: since 2.6.14 a frame
recognised as a tone plays the grid pitch, not the speech pitch. Measured like for like, on
2.6.13's own selection of 17,334 frames and its own calculation, the median correlation is 0.986
(from 0.656; the row 2.6.14 of the table under **Before this change** in [P25 tones sent as voice
frames](#p25-tones-sent-as-voice-frames)); on the selection of 2.6.14's measurement, the frames
its tone test passes between two such frames ([November 2026: tones and the octave
step](#november-2026-tones-and-the-octave-step)), it is 0.986 as well, from 0.388.

On the scoreboard (joint lags, against the commit before), P25 decoding DVSI's bits:

| Partition | Waveform corr. vs DVSI, sample for sample | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| development | 0.9386 → 0.9672 | 1.668 → 1.647 | 0.9887 → 0.9888 | +0.046 → +0.046 | 6.318 → 6.320 |
| validation | 0.9619 → 0.9855 | 1.291 → 1.270 | 0.9872 → 0.9873 | -0.014 → -0.014 | 4.724 → 4.728 |

The waveform correlation rises and the distance from DVSI's decode falls on every vector. The
distance from the input rises by 0.0024 dB (development) and 0.0033 dB (validation). The gate sets
`p25-pitch-development`, `-validation`, `-tones` and `-channel` were declared before scoring;
all four pass. `-tones` checks that the D-STAR and rate-33 decodes are unchanged: P25 has no tone
frames, and no gate set of this step checks P25's tone rows (`r33-mid-pitch-tones`, of the next
step, checks that the mid-frame pitch leaves the P25 tone decodes unchanged).

- **Tones and channel errors.** P25 has no tone frames, and no gate set of this step checks its tone
  rows; on its 52 tone vectors, whose decodes all change, the detector's agreement and wrong-tone
  fraction move by less than 0.01 and the steady tones' level error by at most 0.051 dB on a vector
  (0.555 to 0.560 dB over the steady frames of all 52). DVSI's decoder plays a tone sent as voice
  frames at another pitch than speech's (eq 46's fundamental cut to a 2^-19 grid, as found since),
  so the steady tone frames' median correlation with DVSI's waveform falls from 0.880 to 0.656
  (above; 2.6.14 plays recognised tone frames at that grid pitch). Every P25 channel-error decode
  changes; the mute decisions do not, and on the loud
  frames neither decoder mutes the distance from DVSI's decode falls by 0.007 dB (hard decisions;
  0.006 soft).
- **Other modes.** The D-STAR and rate-33 decodes are unchanged on every speech, tone and
  channel-error vector, and every encoder's bits are unchanged.

PESQ against the input falls by 0.0043 ± 0.0022 (8 of 25 vectors better), ESTOI by 0.0001 ± 0.0000
(6 of 25), and PESQ against DVSI's decode rises by 0.0002 ± 0.0003 (9 of 25). The pitch of a P25 code
now departs from eq 46 by up to 0.071% to follow DVSI's decoder ([Departures from the
specification](#departures-from-the-specification)).

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

**The first implementation** (`d245efd`) used one limit, 0.375, for single and dual tones. P25
speech then moved away from DVSI's decode: its waveform correlation with DVSI's fell by 0.0044 on
the development vectors (16 of 20 lower, `noisy` by 0.058) and by 0.0072 on the validation ones,
all of it from the 142 speech frames that passed single-tone rows and the frames after them in the
same voiced stretch (1,104 changed frames, their residual 0.0017 to 0.0140 of the signal energy),
and its development and validation gate sets failed ([November 2026: tones and the octave
step](#november-2026-tones-and-the-octave-step)). The single-tone limit above, set on the speech
frames DVSI's decoder plays as speech or as tones, replaced it (`b58d91d`).

**Before this change.** 2.6.13 played every P25 frame as speech. Over the 52 tone vectors decoded as
P25, frame for frame at lag 0 (a frame is loud when DVSI's frame RMS is above 30; 17,334 steady
frames, loud with every neighbour loud: 17,286 between two loud frames and 48 a file's first or last
frame, with one neighbour only; 791 first loud frames of a tone; 743 first quiet frames after one),
against DVSI's decode (level: mean absolute difference, dB):

| Build | Steady frames: level | Steady frames: correlation with DVSI's waveform, median | First frame of a tone: level | Last frame: level | Frame after: level |
|---|---|---|---|---|---|
| main (`7fb5ca6`, 2.6.12) | 0.555 | 0.847 | 4.27 | 0.46 | 4.82 |
| phase change (`468e582`) | 0.555 | 0.880 | 4.27 | 0.44 | 4.81 |
| P25 pitch (`470374c`), 2.6.13 | 0.560 | 0.656 | 4.27 | 0.45 | 4.82 |
| 2.6.14 | 0.279 | 0.986 | 0.37 | 0.21 | 3.60 |

(The mean of the steady frames' correlation: 0.582, 0.590, 0.408 and 0.940; the mean over the
vectors of the whole-file correlation: 0.576, 0.550, 0.127 and 0.918.) The row 2.6.14 is this change
set, measured on the same frames with the same calculation: the like-for-like comparison with the
rows above. The measurement is the one made at 2.6.8, which
main reproduces. On 2.6.7, before the subframe synthesis, it gave 1.3 dB for the first frame of a
tone (its faster onset is 2 dB too loud on speech), 0.62 dB and a median correlation of 0.76 for the
steady frames, 1.24 dB for a tone's last frame and 6.7 dB for the frame after it. DVSI's decoder
brings such a tone in at its full level over the first 80 samples of its first frame (0.13, 0.40,
0.61, 0.90 of the tone's level in 20-sample blocks, 65 onsets of `dtmf`), and plays the unvoiced
frame before it, whose amplitudes already hold part of the tone, much quieter than 2.6.13 did; there
the first frame followed the mid-frame model like any voiced onset (0.01, 0.03, 0.04, 0.04, then
0.18, 0.47, 0.71, 0.95). The measurement of [November 2026: tones and the octave
step](#november-2026-tones-and-the-octave-step) selects other frames (those the tone test passes,
between two such frames, with DVSI's frame above a level); on that selection the median correlation
rises from 0.388 to 0.986, which is not the comparison with the rows above.

### Rate 33: mid-frame pitch

When both frames have a voiced band and their pitches differ, DVSI's rate-33 decoder plays the
mid-frame model ([Synthesis in subframes](#synthesis-in-subframes)) at the geometric mean of the
two pitches times a factor a little below 1. TIA-102.BABA-1 defines no mid-frame model, so this is
not a departure from its text.

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
step's direction and the step before, which is not modelled. How DVSI's decoder comes to it is not
known. Half a step of a 2^-10 grid in log2 f0 is -3.4e-4, but the frame-pair data cannot tell a
floor to such a grid from a constant.

The subframe synthesis learns which codec renders through a third rendering,
`MBE_VOICED_SYNTHESIS_DVSI_AMBE2450`, which every rate-33 synthesis call passes (decoded frames,
frames from silence, unheard frames under a mute, repeats, the mute's fade-out and the speech under
a tone's first frame), so that a frame advanced without being rendered leaves the pulse train where
rendering it would. D-STAR, P25 and ProVoice keep `MBE_VOICED_SYNTHESIS_DVSI`, and the public
`mbe_synthesizeSpeechf()` is unchanged.

On the scoreboard (joint lags, against the commit before), rate 33 decoding DVSI's bits:

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
`cos(2 pi F (t - 16) + p)`, p the phase of the output high-pass at F (`mbe_low_band_phase()`,
second order at 100 Hz; [Output high-pass phase](#output-high-pass-phase)). The phase at onset is
the same at every onset of a tone within 0.0002 cycle (20 to 23 onsets per DTMF tone), and the
measured phase minus this formula is +0.0016 cycle on average, standard deviation 0.0022, over the
64 DTMF and KNOX components. The first 80 samples ramp in, as the library already played them
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
  the tone continues that train, as before 2.6.14; when the frame before the tone was not heard, no
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
  line lose the alignment by this section's criterion (counted for the review of 2.6.14).

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
| P25 pitch | eq 46: `w0 = 4 pi / (b0 + 39.5)` | on speech and its pulse-train vectors, that pitch times `1 + g`, g between -0.071% and +0.034%, repeating every half octave of the period; a tone its encoder sent as voice frames at eq 46's fundamental cut to a grid of 2^-19 cycle per sample | a fit to DVSI's pitch on speech, for every P25 frame but a recognised tone frame, which plays the grid pitch ([P25 pitch](#p25-pitch), [P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)) |
| P25 tones sent as voice frames | TIA-102.BABA has no tone frames: every voice frame plays its decoded harmonics at eq 46's pitch | a voice frame whose pitch code and voicing word are those its encoder sends for a DTMF, KNOX or single tone, and whose tone harmonics dominate, plays a tone model: the tone's harmonics at one level, the others at least 27 dB down, eq 46's fundamental cut to a 2^-19 grid; call-progress tones, which DVSI's decoder does not play as equal-amplitude tone models, are left as decoded | DVSI, by a test with two measured limits ([P25 tones sent as voice frames](#p25-tones-sent-as-voice-frames)) |
| Rate-33 harmonic count, b0 17 | Annex A: 11 | 12 | DVSI |
| Rate-33 repeat | eq 57-58: `eps_0 >= 4`, or `eps_0 >= 2` and `eps_T >= 6` | `eps_0 >= 4` | DVSI |
| P25 output under channel errors | eq 112-116 | attenuates progressively | specification |
| Unvoiced synthesis (P25, and rate 33, whose synthesis TIA-102.BABA-1 takes from TIA-102.BABA) | eq 117-126: the 11.2 noise generator, one 256-sample window per frame | a 16-bit generator, `x <- 173 x + 13849`, and 10 ms subframes under a square-root triangular window | DVSI ([Unvoiced synthesis](#unvoiced-synthesis)) |
| Voiced synthesis (P25, rate 33) | eq 127-141: a frame in one piece; harmonics 1-7 at steady pitch interpolate linearly, the others crossfade with `wS`; the phase advances by `l 80 (w0(-1) + w0(0))` per frame, through unvoiced frames too | two 10 ms subframes around a mid-frame model, linear crossfades, a train of pitch pulses that starts afresh after an unvoiced subframe, harmonic phases regenerated from the decoded log2 amplitudes | DVSI ([Synthesis in subframes](#synthesis-in-subframes), [Phase regeneration](#phase-regeneration)) |
| Where voicing starts and ends inside a frame (P25, rate 33) | the frame's V/UV decisions hold for the whole frame | the pitch code of a frame with no voiced band (P25), V/UV row 15 (rate 33) and how far up the voicing reaches (P25) move it to mid-frame | DVSI ([Voiced onsets and offsets](#voiced-onsets-and-offsets)) |
| Rate-33 tone frames | clause 7 gives the tone a frame carries, nothing on the frames around it | the tone crossfades with the speech over 80 samples at both ends | DVSI ([Rate 33 tone transitions](#rate-33-tone-transitions)) |
| Rate-33 tone frequencies | 7.3: the tone synthesized from its decoded frequency parameters; Table 9: the nominal frequencies | harmonics l1 and l2 of Annex J's fundamental on a grid of 2^-19 cycle per sample (indices 145 and 147 at fundamentals of its own), a single tone as harmonic `ID / 16 + 1` of a fundamental under 500 Hz; DTMF within 0.48% of nominal, call progress within 3.41% | DVSI ([Tones as harmonics of a fundamental](#tones-as-harmonics-of-a-fundamental)) |
| Rate-33 silence frames | 4.1 eq 2: L = 14; 4.3-4.4: amplitudes decoded as a voice frame's | 15 bands, amplitudes a further 13.8 dB down, the previous frame's model held to mid-frame | DVSI ([Rate 33 silence frames](#rate-33-silence-frames)) |
| Rate-33 mute | 5.7: bypass the synthesis, output noise uniform over [-5, 5] | noise of RMS about 0.6 (0, +-1 and a few +-2 at 16 bits); the first muted frame fades the last frame out over 80 samples; the model goes on underneath | DVSI's level, as 0 and +-1 ([Rate 33 mutes](#rate-33-mutes)) |
| Rate-33 tone ID 255 | 7.3 and Table 9: a tone of zero amplitude, whatever AD | a lost frame, as an unusable tone index | DVSI ([Rate 33 lost frames and tones](#rate-33-lost-frames-and-tones)) |
| Rate-33 frame repeat after a tone | 7.3 with 5.6: the previous frame's decoded model; silent on a tone before it | the tone repeats at its level until the repeat limit's mute fades it out | DVSI ([Rate 33 lost frames and tones](#rate-33-lost-frames-and-tones)) |
| Amplitude limit (P25, rate 33) | eqs 115-116 scale a frame whose amplitudes sum past 20480 (error-free) | no scaling | DVSI ([Amplitude limit](#amplitude-limit)) |
| Rate-33 V/UV vectors | Annex B rows and eq 5's columns; rows 17-31 all unvoiced | each voiced column reaches 250 Hz into an unvoiced one below; row 17 voices its lowest column | DVSI ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)) |
| Rate-33 V/UV rows 18-31 | Annex B: every band unvoiced, as row 16; 4.1 (Annex A): the fundamental of the frame's pitch code | a band grid whose L bands span 3.58 kHz, `w0 = 2.81 / L`, the pitch code setting only L; in one subframe (rows 24 and 25: in both) some 500 Hz columns play a signal not determined ([Known differences](#known-differences)) | the band grid ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)) |
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
- Rate-33 V/UV rows 1, 3, 13, 15 and 17: DVSI's decoder plays some of their columns between voiced
  and unvoiced, a state the binary voicing of `mbe_parms` cannot hold; this library plays them
  voiced, the nearer value.
- Rate-33 V/UV rows 18 to 31 play here as unvoiced noise on their band grid
  ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)); DVSI's decoder plays more than that. DVSI's
  encoder sends rows 17 to 31 on 12.2% of the voice frames of the speech vectors (row 16 on 37.3%).
  Putting DVSI's samples in place of ours on the frames of rows 17 to 31 would raise PESQ against
  the input by 0.019 ± 0.003 on the 25 development and validation vectors, nearly all of the
  0.021 ± 0.004 by which DVSI's decode of the rate-33 bits scored above ours before the grid; the
  grid takes about 0.003 of it. On plain voice frames (before the grid) our frames correlate with DVSI's
  0.98 on row 16, 0.93 on row 17, 0.89 on rows 18, 19 and 21 to 23, 0.87 on rows 26, 27, 30 and
  31, 0.61 on rows 28 and 29, 0.75 on row 20 and 0.50 on rows 24 and 25.
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
- The pitch code of a P25 or D-STAR frame with no voiced band tells DVSI's decoder more than the
  marks of [Voiced onsets and offsets](#voiced-onsets-and-offsets). Frames at the encoders' default
  codes (P25 b0 25, D-STAR b0 34) match DVSI's output sample for sample, across level steps too
  (median correlation 0.999). At the 56 codes kept for unvoiced frames (P25 b0 128..183, D-STAR b0
  64..119, seven blocks of eight in both) they correlate 0.87-0.98: where the level steps up,
  DVSI plays one half of such a frame, or the start of the next, 1 to 9 dB quieter than this
  library, and which half and by how much goes with the block and with the code's low bits. A
  D-STAR frame of this kind right before a tone plays its second half a median 7.5 dB louder here
  than on DVSI's decoder. It reads as if the code said how the frame's energy lies in its two
  halves; the law is not determined. This library's encoders send neither those codes nor the
  default ones on such a frame, but one fixed code outside the ranges ([Encoder pitch code of
  unvoiced frames](#encoder-pitch-code-of-unvoiced-frames)), so they do not use this channel: every
  frame with no voiced band tells a DVSI decoder, as far as the codes measured show, that the
  voicing before lasts through its first 10 ms if it reaches 505 Hz, where DVSI's encoders end it at
  the frame's start on about a third of them.
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
- P25 voiced onsets that DVSI's decoder starts a subframe late and the library starts with the
  frame ([Voiced onsets and offsets](#voiced-onsets-and-offsets)). Measured on 2.6.14's decode with
  two test builds outside the repository, one that starts every onset late (no voiced band in the
  mid-frame model after a frame with no voiced band) and one that starts none late, over the 1,940
  onsets of P25 voiced runs on the 25 development and validation vectors (a voiced frame after a
  frame with no voiced band, followed by two voiced frames, the two frames after the onset frame
  together, 320 samples, with a standard deviation of at least 300 in DVSI's decode; a file's first
  three and last four frames left out), each compared with DVSI's decode over the 400 samples from
  the onset frame's mid-point: DVSI's decode follows the early start in 1,771 (correlation above 0.8
  and more than 0.1 above the other build's), the late start in 72, both in 83 (both above 0.8,
  neither more than 0.1 above the other) and neither in 14. Of the 72 late ones the library starts 43 late
  (the frames voiced only below 505 Hz) and 29 early; those 29 are the open item. 1,898 onsets are in
  line now (above 0.8), and 1,926 would be with the better start for each. In 28 of the 29 the
  harmonic at or next above 505 Hz is unvoiced, as in all 43 the library starts late, but so it is
  in 42 onsets that DVSI's decoder starts early; a rule is under study, and none is adopted here.
- A single band whose voicing changes between two frames, on a DVMEGA DVstick 30
  ([#157](https://github.com/arancormonk/mbelib-neo/issues/157)): a band above 1 kHz that turns
  voiced starts with the frame and one that turns unvoiced lasts to mid-frame, and the lowest band
  does the opposite on both. The mid-frame model here takes the current frame's voicing for every
  band ([Synthesis in subframes](#synthesis-in-subframes)).
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
- Pitch per code: P25 now plays within 0.0064% rms of DVSI's decoder over the 69 rows (61 pitch
  codes) it was fitted on, weighted as fitted, and 0.018% over the validation vectors' 9 codes ([P25
  pitch](#p25-pitch)), D-STAR within 0.013% rms (at most 0.04% on a code with 10 or more frame
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
  changes of 2.6.14 0.0191, 0.0633 and 0.0527 ([November 2026: tones and the octave
  step](#november-2026-tones-and-the-octave-step)); the pooled figures and their parts below were
  measured for 2.6.13. Of it:
  - Voiced runs that line up with DVSI's (run correlation 0.8 or more) hold 0.018, 0.021 and 0.033.
    P25's no longer grow with the run's length (0.019 to 0.028 over frames 0 to 3, 4 to 7 and so
    on); D-STAR's are 0.016 to 0.018 over a run's first twenty frames and 0.037 past them; rate 33's
    still grow (0.018, 0.024, 0.036, 0.053, 0.067 and 0.067).
  - Voiced runs that do not line up hold most of the rest: D-STAR 45 of 1911 runs, holding 0.017 of
    the signal energy as residual, 44 of them off from their first frames (19 after a frame whose
    pitch code is one of the 56 kept for frames with no voiced band and names a pattern for the
    frame's own subframe, 12 after one naming both subframes, 4 after one naming the mid-frame
    subframe and 9 after a plain unvoiced frame; counted for 2.6.14 on the development and
    validation vectors, three voiced frames or more, the octave step brings every out-of-line onset
    after a code that pulses the own subframe or both in line, [The pulse train across an
    octave](#the-pulse-train-across-an-octave)); P25 85 of 1959 runs holding 0.045, 83 off from
    their first frames, which no single lag lines up (of the 77 whose first three frames correlate
    under 0.8, 22 reach 0.85 at any lag within 130 samples); rate 33 43 of 1825 runs holding 0.030,
    31 of them off from their first frames (0.012) and 12 starting aligned and losing it inside the
    run (0.018; median length 21 frames).
  - Frames with no voiced band hold 0.003 in each codec: the pitch codes and V/UV rows that tell
    DVSI's decoder more than is modelled (above).
  - Inside voiced runs rate 33 loses the alignment with DVSI's waveform on 0.6% of its frames
    (D-STAR 0.2%, P25 0.3%; frames whose two frames before correlate 0.9 or more and whose two
    frames after fall under 0.8), more often across large pitch steps (over 20%: 6 of 65 in rate
    33). Of the frames at the largest steps, a frame-to-frame pitch ratio of 0.30 to 0.42 (three in
    rate 33, one in D-STAR), all four come back in line with the octave step of 2.6.14 ([The pulse
    train across an octave](#the-pulse-train-across-an-octave)). What moves the pulse train at the
    others, as what sets rate 33's mid-frame pitch ([Rate 33: mid-frame
    pitch](#rate-33-mid-frame-pitch)), is not known.
- DVSI's rate-33 encoder sends V/UV indices 17-31 on some voice frames, where clause 4.2
  searches 0-16 and keeps the lowest index on ties, and this library's encoder never sends them.
  Rows 18-31 decode all unvoiced here, as row 16 does, but on their band grid; row 17 voices its
  lowest column ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)).
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
  waveform from 0.880 to 0.656 (0.847 on main). 2.6.14 plays such frames at the pitch DVSI's decoder
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
  voice frames at another pitch ([P25 pitch](#p25-pitch); since 2.6.14 the library plays it so:
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

The head of PR #196 (`fac0111`, 2.6.13) against the changes of 2.6.14:

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
  steady loud frames of 2.6.13's measurement, [P25 pitch](#p25-pitch)). The like-for-like comparison
  is the row 2.6.14 under **Before this change** in [P25 tones sent as voice
  frames](#p25-tones-sent-as-voice-frames): on those frames the median correlation is 0.986 (from
  0.656 on 2.6.13), the mean 0.940 (from 0.408), the steady frames' level error 0.279 dB (from
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

The three changes of 2.6.14 were scored together, `fac0111` against the change set, on the
development, validation, tone and channel partitions in one combined run with the four gate sets
declared in `c736edf`, before any build of the change set was scored. The first build of the P25
tone path, `d245efd`, failed `tone-model-development` and `-validation` on P25; its single-tone
limit was corrected (`b58d91d`) and the change set scored again with the same gate sets, builds of
the baseline and protocol. `tone-model-tones` failed the same five checks both times and is kept as
declared ([November 2026: tones and the octave step](#november-2026-tones-and-the-octave-step)).
The final partition was scored once for 2.6.14, on 2026-11-26, after the review: `fac0111` against
`c2ba220`, the reviewed code, against `tone-model-final`, declared in `c736edf`. It passes all 24
checks ([November 2026: tones and the octave step](#november-2026-tones-and-the-octave-step)).
