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
  retrospective.
- **Hardware.** Where a contributor measured a real AMBE-3000 (a DVMEGA DVstick 30 or a ThumbDV)
  on constructed frames, those measurements decide between rules DVSI's vectors leave open; each
  section names the issue that carries the data.

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
[-1.5, 1.5] on the int16 scale, which the int16 conversion truncates to the same mix
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
from silence. The gain is a fit to frame levels; how the chip does this is not known. The mean
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

### Rate 33 tone transitions

DVSI's rate-33 decoder joins a tone to the frames around it over 80 samples, as its D-STAR decoder
does ([D-STAR tone transitions](#d-star-tone-transitions)); TIA-102.BABA-1 says nothing on it. The
tone fades in linearly over the first 80 samples of its first frame while the last frame's speech,
played once more, fades out; and it goes on from its phases through the first 80 samples of the
frame after its last, fading out while that frame fades in. A tone that changes while it sounds
changes at once. The library cut the tone in and out at the frame edges.

Over DVSI's 52 tone vectors, in 20-sample blocks relative to the tone's level, DVSI's medians are
0.15, 0.37, 0.63, 0.87 going in (1102 starts) and 0.90, 0.65, 0.39, 0.16 going out (1099 ends); ours
are now 0.14, 0.38, 0.63, 0.88 and 0.89, 0.65, 0.40, 0.15. Frame for frame, the level of a tone's
first frame comes from 1.83 to 0.14 dB of DVSI's and that of the frame after a tone from 22.1 to
0.58 dB; where a tone starts on loud speech the first 40 samples come from 9.7 to 0.9 dB (123
frames). The tone a frame leaves sounding is kept in the unused element 0 of `Vl` of the last
synthesized frame; a frame the error rate mutes plays its noise alone. Declared gate that fails
(`ambe2450-tone-fade-tones`): on `dtmfvbad` the tone detector agrees with DVSI's output on 98.96%
of its windows, from 99.93%, while the frames of that vector with a 40-sample block more than 6 dB
from DVSI's fall from 459 to 256; the other 51 vectors keep their agreement, and no speech vector
changes.

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
clears the overlap with its last synthesized frame, and the voice after a tone fades in
([Rate 33 tone transitions](#rate-33-tone-transitions)). The two renderings keep different things in `previousUw` and `noiseOverlap`: a
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
  frame's and this one's stands at mid-frame: its pitch the geometric mean of the two, its
  harmonic count that pitch's, each band's voicing this frame's at the band's frequency, its
  magnitudes the two frames' at that frequency, met in the log domain. Unvoiced magnitudes enter
  that mean without their factor `0.2046 / sqrt(w0)` and take the mid-frame pitch's on the way
  out. With one frame voiced and the other not, the mid-frame pitch is the voiced frame's.
- **Linear crossfades.** Each subframe plays the model at its start with a falling weight and
  the model at its end with a rising one.
- **Pitch pulses.** The voiced part is a train of pulses. Each pulse plays a model's voiced
  harmonics, their phases the regenerated ones, peaking 16 samples after the pulse, through a
  window that rises over 16 samples (a raised cosine), stays whole until one pitch period after
  the pulse and falls over the 16 samples after that. A response so lasts a period and 16
  samples and overlaps the next by 16, where the two windows add up to one: the responses of a
  steady train add up to the harmonics. The pulse rate moves linearly from one model's
  pitch to the next's across a subframe; into or out of voicing it keeps the voiced model's.
  A subframe with no voiced band at either end starts the train afresh.

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

With these, voiced runs after two or more unvoiced frames start where DVSI's do in 95% of cases
for P25 (791 of 831 on ten vectors), 96% for D-STAR (812 of 846) and 96% for rate 33 (775 of
810). Not resolved: of the 64 P25 onsets voiced in the lowest band and in bands above it, but not
at 505 Hz, 43 start with the frame, 17 in its second subframe and 4 elsewhere, and no rule
separates them (the width of the widest upper voiced region does for 53 of 57). See also
[Known differences](#known-differences) for what the unvoiced pitch codes and rate 33's rows 17
to 31 further tell DVSI's decoder.

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

## Departures from the specification

| Item | Specification text | DVSI's AMBE-3000 | This library |
|---|---|---|---|
| Rate-33 u3 bits 8-5 | Table 8: b3(0), b4(2), b4(1), b4(0) | b4(2), b4(1), b4(0), b3(0) | DVSI |
| Rate-33 pitch | Annex A table | smooth log2 law | DVSI |
| Rate-33 harmonic count, b0 17 | Annex A: 11 | 12 | DVSI |
| Rate-33 repeat | eq 57-58: `eps_0 >= 4`, or `eps_0 >= 2` and `eps_T >= 6` | `eps_0 >= 4` | DVSI |
| P25 output under channel errors | eq 112-116 | attenuates progressively | specification |
| Unvoiced synthesis (P25, and rate 33, whose synthesis TIA-102.BABA-1 takes from TIA-102.BABA) | eq 117-126: the 11.2 noise generator, one 256-sample window per frame | a 16-bit generator, `x <- 173 x + 13849`, and 10 ms subframes under a square-root triangular window | DVSI ([Unvoiced synthesis](#unvoiced-synthesis)) |
| Voiced synthesis (P25, rate 33) | eq 127-141: a frame in one piece; harmonics 1-7 at steady pitch interpolate linearly, the others crossfade with `wS`; the phase advances by `l 80 (w0(-1) + w0(0))` per frame, through unvoiced frames too | two 10 ms subframes around a mid-frame model, linear crossfades, a train of pitch pulses that starts afresh after an unvoiced subframe | DVSI ([Synthesis in subframes](#synthesis-in-subframes)) |
| Where voicing starts and ends inside a frame (P25, rate 33) | the frame's V/UV decisions hold for the whole frame | the pitch code of a frame with no voiced band (P25), V/UV row 15 (rate 33) and how far up the voicing reaches (P25) move it to mid-frame | DVSI ([Voiced onsets and offsets](#voiced-onsets-and-offsets)) |
| Rate-33 tone frames | clause 7 gives the tone a frame carries, nothing on the frames around it | the tone crossfades with the speech over 80 samples at both ends | DVSI ([Rate 33 tone transitions](#rate-33-tone-transitions)) |
| Amplitude limit (P25, rate 33) | eqs 115-116 scale a frame whose amplitudes sum past 20480 (error-free) | no scaling | DVSI ([Amplitude limit](#amplitude-limit)) |
| Rate-33 V/UV vectors | Annex B rows and eq 5's columns; rows 17-31 all unvoiced | each voiced column reaches 250 Hz into an unvoiced one below; row 17 voices its lowest column | DVSI ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)) |
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
- Rate-33 V/UV rows 18 to 31 play as plain noise here, as Annex B gives them and as row 16 does.
  DVSI's decoder plays them otherwise, and DVSI's encoder sends rows 17 to 31 on 12.2% of the voice
  frames of the speech vectors (row 16 on 37.3%). Against a frame of row 16, rows 18 to 23 change
  the frame's second half and the start of the next frame, rows 26 to 31 its first half, rows 24
  and 25 both; and part of what DVSI plays there is impulsive (kurtosis 13.7 on row 20 and about 7
  on rows 24 and 25, against 2.7 for noise) without being periodic at the frame's pitch. Our
  frames correlate with DVSI's 0.99 on row 16, 0.91-0.96 on rows 18, 19, 21 to 23, 26, 27, 30 and
  31, 0.82-0.83 on rows 28 and 29, 0.74 on row 20 and 0.34-0.47 on rows 24 and 25. What each row
  tells the decoder is not determined.
- The pitch code of a P25 or D-STAR frame with no voiced band tells DVSI's decoder more than the
  marks of [Voiced onsets and offsets](#voiced-onsets-and-offsets). Frames at the encoders' default
  codes (P25 b0 25, D-STAR b0 34) match DVSI's output sample for sample, across level steps too
  (median correlation 0.999). At the 56 codes kept for unvoiced frames (P25 b0 128..183, D-STAR b0
  64..119, seven blocks of eight in both) they correlate 0.87-0.98: where the level steps up,
  DVSI plays one half of such a frame, or the start of the next, 1 to 9 dB quieter than this
  library, and which half and by how much goes with the block and with the code's low bits. A
  D-STAR frame of this kind right before a tone plays its second half a median 7.5 dB louder here
  than on DVSI's decoder. It reads as if the code said how the frame's energy lies in its two
  halves; the law is not determined.
- P25 tones. P25 has no tone frames, and DVSI's encoder sends a tone it detects as voice frames
  at a pitch whose harmonics carry it, with one frame of that pitch and no amplitude after it.
  DVSI's decoder brings such a tone in at its full level over the first 80 samples of its first
  frame (0.13, 0.40, 0.61, 0.90 of the tone's level in 20-sample blocks, 65 onsets of `dtmf`), as
  its rate-33 decoder does a tone frame, and plays the unvoiced frame before it, whose amplitudes
  already hold part of the tone, much quieter than this library. Here the first frame follows the mid-frame model like
  any voiced onset (0.01, 0.03, 0.04, 0.04, then 0.18, 0.47, 0.71, 0.95), which on speech is what
  DVSI's decoder does at every size of step (2106 P25 onsets, also those whose two strongest
  harmonics hold 95% of the energy). What tells its decoder the two apart was not found. Over
  the 52 tone vectors decoded as P25, frame for frame, the first frame of a tone lies 4.3 dB from
  DVSI's level (791 frames; 1.3 dB on main, whose faster onset is 2 dB too loud on speech), the
  steady frames 0.56 dB (0.62 on main; their median correlation with DVSI's waveform 0.85, from
  0.76), a tone's last frame 0.46 dB (1.24) and the frame after it 4.8 dB (6.7).
- P25 voiced onsets whose voicing covers the lowest band and bands above it but not 505 Hz: 17 of
  64 start in the frame's second subframe on DVSI's decoder and no rule found separates them
  ([Voiced onsets and offsets](#voiced-onsets-and-offsets)); all start with the frame here.
- A single band whose voicing changes between two frames, on a DVMEGA DVstick 30
  ([#157](https://github.com/arancormonk/mbelib-neo/issues/157)): a band above 1 kHz that turns
  voiced starts with the frame and one that turns unvoiced lasts to mid-frame, and the lowest band
  does the opposite on both. The mid-frame model here takes the current frame's voicing for every
  band ([Synthesis in subframes](#synthesis-in-subframes)).
- Encoder: pauses keep their level, where DVSI's encoders suppress background noise
  ([Encoder voicing over background noise](#encoder-voicing-over-background-noise)). In the
  quietest quarter of each speech vector's frames, our rate-33 encoding decodes a median 9.2 dB
  above DVSI's encoding of the same input (0.2 to 13.6 dB per vector, `dam40` and `dam80` aside,
  whose input lies 40 and 80 dB down); the loudest quarter sits within 2.0 dB on every vector but
  `dam10` (+4.2).
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
- Small pitch differences remain: P25 codes up to 0.07% from DVSI's (rounding-like), D-STAR codes
  within 0.03%.
- DVSI's rate-33 encoder sends V/UV indices 17-31 on some voice frames, where clause 4.2
  searches 0-16 and keeps the lowest index on ties, and this library's encoder never sends them.
  Rows 18-31 decode all unvoiced here, as row 16 does; row 17 voices its lowest column
  ([Rate 33: V/UV vectors](#rate-33-vuv-vectors)).

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
