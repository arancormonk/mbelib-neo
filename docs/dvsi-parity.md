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
frame still plays 12 dB low, as the voice after any repeat-limit mute fades in. The reporter's patch
was the starting point. Only frame levels were measured, so the fade's shape and the continuing phases
are the library's choice; other dual indices and volumes, runs of 2 or more than 4 and single tones
other than 32 were not measured. No DVSI vector changes.

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

### Voiced amplitude

Within a frame DVSI's decoders move a voiced harmonic from one frame's amplitude to the next's
differently from TIA-102.BABA eq 133-135. On a DVMEGA DVstick 30
([#157](https://github.com/arancormonk/mbelib-neo/issues/157)), D-STAR gain steps of 7.0 to 41.6 dB,
up and down, at three positions in a stream, reach half their settled change in dB within one
16-sample block of sample 80 of the step frame at every size (24 of 30 at 80 exactly). The
specification's rendering, eq 135's linear amplitude for harmonics 1-7 and eq 133's crossfade with
the window `wS` above them, reaches it at 16 to 144 samples depending on the step. On every codec
path, a harmonic voiced in both frames at steady pitch (eq 134's pitch condition) is now
interpolated whatever its number, with the fourth root of its amplitude moving linearly over the
frame:

```
a_l(n) = ((1 - n/N) M_l(-1)^(1/4) + (n/N) M_l(0)^(1/4))^4
```

Over the 30 gain steps, the mean error of the level curve (the level change of each 16-sample
block of the step frame as a fraction of the settled change) against the unit's falls from 0.103 to
0.021; interpolating every harmonic with linear amplitude gives 0.110, a geometric mean 0.048 and
the fourth root, the best of the power laws tried, 0.021. The curves of the unit's b0 and b1 steps
stay as close as before. A fade from or to a silent model (after a mute, a tone or the repeat
limit) keeps the specification's rendering: there the law would fade in as (n/N)^4, and the
repeat-limit fade was matched to DVSI's with the specification's.

The law carries less energy between two levels than linear amplitude: with it, DVSI's voiced
frames played 0.11-0.16 dB louder than ours again (per-frame medians, development vectors). The
voiced render gain on the codec paths is therefore +0.45 dB, chosen from +0.40, +0.45 and +0.50 dB
on development, which puts those medians within 0.03 dB of where +0.3 dB put them with the
specification's rendering. The public `mbe_synthesizeSpeechf()` keeps the specification's
rendering and +0.3 dB. Measured on the unit's 179 streams from #157, which played no part in
choosing the gain, frame levels come closer to its output in every group of streams (mean
difference 0.62 to 0.53 dB).

On DVSI's bits, each of the 25 development and validation vectors comes closer to DVSI's decode in
each of the three modes, and on average closer to the input speech; 8 of the 75 move away from the
input by at most 0.03 dB (`dam10` in every mode). The final partition, scored once, comes closer
to both on all six of its vectors. Scored at joint lags (`voiced-amplitude` gate sets):

| Codec | Partition | LSD vs DVSI (dB) | Envelope corr. vs DVSI | Level vs DVSI (dB) | LSD vs input (dB) |
|---|---|---|---|---|---|
| D-STAR | development | 3.705 → 3.347 | 0.9523 → 0.9535 | +0.049 → +0.070 | 6.578 → 6.521 |
| D-STAR | validation | 3.191 → 2.899 | 0.9483 → 0.9496 | +0.025 → +0.055 | 4.990 → 4.921 |
| D-STAR | final | 4.449 → 4.162 | 0.9218 → 0.9222 | +0.168 → +0.176 | 5.759 → 5.698 |
| P25 | development | 3.798 → 3.455 | 0.9555 → 0.9561 | +0.027 → +0.045 | 6.514 → 6.454 |
| P25 | validation | 3.180 → 2.904 | 0.9557 → 0.9565 | -0.034 → -0.013 | 4.900 → 4.820 |
| P25 | final | 4.146 → 3.886 | 0.9286 → 0.9289 | +0.123 → +0.100 | 5.640 → 5.548 |
| Rate 33 | development | 3.622 → 3.260 | 0.9527 → 0.9555 | -0.244 → -0.247 | 6.513 → 6.454 |
| Rate 33 | validation | 3.055 → 2.747 | 0.9521 → 0.9542 | -0.074 → -0.074 | 5.021 → 4.940 |
| Rate 33 | final | 4.204 → 3.895 | 0.9321 → 0.9331 | +0.074 → +0.034 | 5.793 → 5.704 |

Most of the gain in log-spectral distance comes from interpolating every harmonic (linear
amplitude alone: -0.31 dB on development in every mode); the fourth-root law adds the rest and
the gains in envelope correlation and in distance from the input. Under channel errors, the
frames that neither decoder mutes come 0.39, 0.24 and 0.48 dB closer to DVSI's (D-STAR, P25,
rate 33, hard decisions) with mute decisions unchanged, and the encoders' bits are unchanged. On
the tone vectors, rate 33's `alltone` agrees with DVSI's tone detector on two more windows and
`dtmf180` on one fewer: at -49 dBFS, the voice frame right after a run of tone frames (frame
1504, rendered with the crossfade as before, as its pitch jumps) now reads as a dual tone where
DVSI's reads none, from the 0.15 dB of voiced gain, so `dvsi-parity-tones` and
`voiced-amplitude-tones` fail on that window.

## Departures from the specification

| Item | Specification text | DVSI's AMBE-3000 | This library |
|---|---|---|---|
| Rate-33 u3 bits 8-5 | Table 8: b3(0), b4(2), b4(1), b4(0) | b4(2), b4(1), b4(0), b3(0) | DVSI |
| Rate-33 pitch | Annex A table | smooth log2 law | DVSI |
| Rate-33 harmonic count, b0 17 | Annex A: 11 | 12 | DVSI |
| Rate-33 repeat | eq 57-58: `eps_0 >= 4`, or `eps_0 >= 2` and `eps_T >= 6` | `eps_0 >= 4` | DVSI |
| P25 output under channel errors | eq 112-116 | attenuates progressively | specification |
| Voiced amplitude within a frame (P25, and rate 33, whose synthesis TIA-102.BABA-1 takes from TIA-102.BABA) | eq 133-135: harmonics 1-7 at steady pitch interpolate linearly, the others crossfade with `wS` | every harmonic at steady pitch, its fourth root linear | DVSI ([Voiced amplitude](#voiced-amplitude)) |

The P25 attenuation is not adopted: no rule tried (eq 112-116 as printed, and instantaneous,
clamped or leaky variants of it) explains it, and under channel errors our output stays closer to
the input speech than DVSI's (log-spectral distance from the input at 10% bit errors 10.4 against
12.1 dB).

## Known differences

- Soft-decision D-STAR streams are not muted on the error level, where DVSI mutes some of them
  (see above).
- After a fully unvoiced gap DVSI repositions the voiced pulse train in a way no simple rule
  reproduced; within voiced runs, including pitch steps, our pulse positions follow DVSI's.
- DVSI hard-clips at full scale; this library clips at 0.95 of full scale. Clipping at full
  scale moved no metric, so it is not adopted.
- A frame at the largest D-STAR gain and b3 (b2 63, b3 511; `9E B0 6F EE F8 7B 9F FE 74`) held
  for 50 frames ([#158](https://github.com/arancormonk/mbelib-neo/issues/158)): a DVMEGA DVstick 30
  plays -18.4 dBFS on the first frame and -9.0 to -11.4 dBFS after it, this library -13.0 and -6.5
  to -8.4 dBFS. It used to play silence there from the 5th frame, when the spectral enhancement's
  weights overflowed float.
- After a D-STAR b0 124 frame (the null frame among them), a DVMEGA DVstick 30 plays the lost
  frames of a run from the 4th or 5th on at a comfort-noise level that follows the gain field b2
  of the b0 124 frames (-75 dBFS after the null frame, up to -2 dBFS after b2 63), and plays all
  lost frames at level directly after one ([#156](https://github.com/arancormonk/mbelib-neo/issues/156));
  this library mutes from the 5th, at about -95 dBFS.
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
