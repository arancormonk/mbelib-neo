// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Speech synthesis in 10 ms subframes, as DVSI's decoders render it.
 *
 * The codec paths (MBE_VOICED_SYNTHESIS_DVSI for D-STAR,
 * MBE_VOICED_SYNTHESIS_DVSI_IMBE for P25, MBE_VOICED_SYNTHESIS_DVSI_PROVOICE
 * for ProVoice,
 * MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 for rate 33) synthesize a frame as two subframes of
 * MBE_SUBFRAME samples (docs/dvsi-parity.md, "Synthesis in subframes"):
 *  - the frame's own model stands at its end; a mid-frame model, between the
 *    previous frame's and this one's, stands at sample 80;
 *  - each interval between two models crossfades linearly from one to the next;
 *  - the voiced part is a train of pitch pulses, each playing a model's
 *    impulse response: its voiced harmonics, peaking MBE_PULSE_PEAK samples
 *    after the pulse, through a window that rises over MBE_PULSE_PEAK samples
 *    (a raised cosine), stays whole until one pitch period after the pulse
 *    and falls over MBE_PULSE_PEAK more, so a response lasts a period and
 *    MBE_PULSE_PEAK samples and the responses of a steady train add up to
 *    the harmonics; the pulse rate moves linearly from one model's pitch to
 *    the next's across an interval, after taking the train an octave down or
 *    up first where the new period is within 20% of twice or of half the old;
 *  - the unvoiced part is mbe_synthesizeUnvoicedSubframes().
 */

#ifndef MBEINT_MBE_SUBFRAME_H
#define MBEINT_MBE_SUBFRAME_H

#include "mbe_adaptive.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

/* Samples from a pitch pulse to the peak of the response it plays, and the
 * length of each raised-cosine end of the response's window. */
#define MBE_PULSE_PEAK                 16

/* An unvoiced band's amplitude is its magnitude times this over sqrt(w0); the
 * mid-frame model's magnitudes are interpolated with that factor taken out,
 * so that a band voiced in one frame and unvoiced in the other moves between
 * like quantities. */
#define MBE_UNVOICED_MAGNITUDE_FACTOR  0.2046f

/* Pulses a stream keeps from one frame to the next (those whose responses can
 * still sound), and where: element 0 of noiseOverlap holds their count, the
 * elements after it their times in samples before the frame's end, the one
 * after those the pulse phase reached at the frame's end, in cycles, and the
 * one after that (MBE_PULSED_SLOT) whether the frame's mid-frame model emitted
 * the pulse of the pulsed state: exactly MBE_PULSED_EMITTED if it did, 0 if
 * not (D-STAR's and P25's renderings write it; the other subframe renderings
 * leave 0). In P25's rendering alone the two after that hold the pulsed
 * state's validity counter: MBE_PULSED_READ_SLOT whether the frame's reserved
 * code was read (exactly MBE_PULSED_READ if it was, 0 if not), and
 * MBE_PULSED_COUNT_SLOT the counter as the frame left it, a whole number
 * MBE_PULSED_COUNT_MIN..MBE_PULSED_COUNT_MAX as a float; every other subframe
 * rendering leaves 0 in both, and reads neither. Each is read by its bit
 * pattern, only in a state with the subframe rendering's mark
 * (MBE_SUBFRAME_STATE_MARK); any value a frame does not store reads as a fresh
 * stream's (not emitted, not read, a count of 0). The TIA-102.BABA rendering
 * keeps its noise overlap in all these elements instead. */
#define MBE_PULSE_HISTORY              32
#define MBE_PULSED_SLOT                (MBE_PULSE_HISTORY + 2)
#define MBE_PULSED_EMITTED             1.0f
#define MBE_PULSED_READ_SLOT           (MBE_PULSE_HISTORY + 3)
#define MBE_PULSED_READ                1.0f
#define MBE_PULSED_COUNT_SLOT          (MBE_PULSE_HISTORY + 4)
#define MBE_PULSED_COUNT_MIN           (-10)
#define MBE_PULSED_COUNT_MAX           3
/* In the rate-33 rendering alone, the element after those holds rate 33's
 * pulsed-row state (docs/dvsi-parity.md, "Rate 33: the pulsed rows"): the
 * word (mbe_pulsed_row_value(), mbe_adaptive.h, not held: a row 18..25 and a
 * pitch code) of a frame whose own model had a pulsed band, 0 otherwise, so
 * that the next frame plays that model's pulsed bands and pulse again. Every
 * other subframe rendering (D-STAR's, P25's, ProVoice's) writes 0 there and
 * reads it not; the TIA-102.BABA rendering keeps its noise overlap in this
 * element as in all the others, reading and writing it. It is read by its bit
 * pattern, only in a state with the subframe rendering's mark, and only when
 * the last frame's model is the band-grid model of the word's code; any value
 * no frame stores (a held word, a row 26..31, anything that is no word) reads
 * as none. A caller's word of rows 18..25 there with such a last model is, to
 * the library, the state of a stream that had such a frame (two codes of one
 * harmonic count share a grid model). */
#define MBE_PULSED_ROW_SLOT            (MBE_PULSE_HISTORY + 5)
/* In D-STAR's, P25's and rate 33's renderings alone, the two elements after
 * that hold DVSI's adaptive low cut (mbe_low_band.h; docs/dvsi-parity.md,
 * "DVSI's low band"): MBE_LOW_CUT_TRACK_SLOT the tracker as the frame left it,
 * in Hz, and MBE_LOW_CUT_CORNER_SLOT the corner the frame applied to its own
 * model (min of the tracker it started from and MBE_LOW_CUT_CEILING_HZ), so
 * that the next frame renders that model with the corner it had. ProVoice's
 * rendering writes 0 in both and reads neither. Each is read by its bit
 * pattern, only in a state with the subframe rendering's mark, and only within
 * the range a frame stores (the tracker MBE_LOW_CUT_TRACK_MIN_HZ to
 * MBE_LOW_CUT_TRACK_MAX_HZ, the corner MBE_LOW_CUT_TRACK_MIN_HZ to
 * MBE_LOW_CUT_CEILING_HZ); any other value (0, a negative or subnormal value,
 * an infinity, not a number) reads as a fresh stream's: the tracker at
 * MBE_LOW_CUT_START_HZ, the last corner the one this frame applies. The
 * TIA-102.BABA rendering keeps its noise overlap in these elements too. */
#define MBE_LOW_CUT_TRACK_SLOT         (MBE_PULSE_HISTORY + 6)
#define MBE_LOW_CUT_CORNER_SLOT        (MBE_PULSE_HISTORY + 7)

/* Samples from a model's point (0, 80 or 160 in its frame) to the one pulse a
 * model with pulsed bands emits (docs/dvsi-parity.md, "D-STAR and P25: the
 * pulsed state"). */
#define MBE_PULSED_DELAY               24

/* Rate 33's pulse times (docs/dvsi-parity.md, "Rate 33: the pulsed rows"): a
 * model of a frame of V/UV rows 18 to 23 or 26 to 31 with a pulsed band emits
 * its pulse at point + MBE_AMBE2450_PULSE_SLOPE b0 - MBE_AMBE2450_PULSE_OFFSET
 * samples, point being 80 for the mid-frame model and 160 for the frame's own;
 * rows 24 and 25, which pulse both subframes, code two positions on a grid of
 * MBE_AMBE2450_PULSE_GRID in b0, the mid-frame model's pulse at 80 -
 * MBE_AMBE2450_PULSE_GRID_ORIGIN + MBE_AMBE2450_PULSE_GRID_STEP (b0 mod
 * MBE_AMBE2450_PULSE_GRID) and the own model's at 160 - ORIGIN + STEP (b0 div
 * MBE_AMBE2450_PULSE_GRID). Each time is formed in double from the integers and
 * these constants and rounded to float once. All are fits to DVSI's rate-33
 * decodes of the development and validation vectors, not values from a
 * document: slopes from 1.012 to 1.024 with the pulse of b0 87 between 7.25 and
 * 8.0 samples after its model's point (here 7.50), offsets from 80.29 to 81.34,
 * measured alike, so neither is known more finely than those ranges; for the grid,
 * origins of 75.1 and 76.1 measured a little worse than 75.6, and no finer
 * variation of the origin or of the step (11.875) was tried. */
#define MBE_AMBE2450_PULSE_SLOPE       1.016
#define MBE_AMBE2450_PULSE_OFFSET      80.89
#define MBE_AMBE2450_PULSE_GRID        11
#define MBE_AMBE2450_PULSE_GRID_STEP   11.875
#define MBE_AMBE2450_PULSE_GRID_ORIGIN 75.6

/** One subframe's model. */
struct mbe_subframe_model {
    float w0;                 /**< Fundamental, radians per sample; 0 for none. */
    int L;                    /**< Harmonic count. */
    int any_voiced;           /**< Whether a band is voiced. */
    unsigned char voiced[57]; /**< [1..L]: the band's state: unvoiced (0), voiced (1), none (2) or pulsed (3). */
    float magnitude[57];      /**< [1..L]: the band's amplitude, as mbe_parms.Ml holds it (a pulsed band's as a
                                   voiced one's: the unvoiced factor taken out). */
    float phase[57];          /**< [1..L]: regenerated phase of a voiced or pulsed band (from log2Ml). */
    int any_pulsed;           /**< Whether a band is pulsed (only in a model with no voiced band). */
};

/**
 * @brief The model a frame's parameters give at its end.
 *
 * Its phases are regenerated from the frame's decoded log2 amplitudes
 * (log2Ml), not from the amplitudes it renders, as DVSI's decoders regenerate
 * them (docs/dvsi-parity.md, "Phase regeneration"); a log2Ml that is not a
 * finite number counts as the kernel's floor.
 * @param mp    Frame parameters.
 * @param model Receives the model; it has no band when w0 is not a pitch
 *              (positive, finite and not subnormal: a subnormal value is zero
 *              to the arithmetic where the caller's floating-point environment
 *              flushes denormals) or the harmonic count is not 1..56.
 */
void mbe_subframe_own_model(const mbe_parms* mp, struct mbe_subframe_model* model);

/**
 * @brief The model a frame's parameters give at its end in @p rendering.
 *
 * mbe_subframe_own_model() but in P25's rendering
 * (MBE_VOICED_SYNTHESIS_DVSI_IMBE), whose phase kernel reads each band the
 * frame holds as unvoiced with the unvoiced factor taken out: its log2Ml minus
 * log2(MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt(w0)), the domain of
 * TIA-102.BABA-1's log magnitudes in which the AMBE codecs' decoded unvoiced
 * bands already are (docs/dvsi-parity.md, "Phase regeneration"); a tone model
 * (MBE_FRAME_MARK_TONE) and a band with no state are read as decoded.
 * @param mp        Frame parameters.
 * @param rendering The rendering.
 * @param model     Receives the model, as mbe_subframe_own_model() says.
 */
void mbe_subframe_rendered_own_model(const mbe_parms* mp, enum mbe_voiced_synthesis rendering,
                                     struct mbe_subframe_model* model);

/* The rate-33 rendering plays the mid-frame pitch between two voiced frames of
 * different pitch at the geometric mean of their pitches times this, times the
 * half-octave wave below (docs/dvsi-parity.md, "Rate 33: mid-frame pitch").
 * A measured offset, not a known mechanism: a least-squares fit of one shift
 * of DVSI's rate-33 pulse train against ours for each mid-frame pitch taken
 * between two different pitch codes gives -0.010 samples over the frame pairs
 * measured and -0.011 over those whose two code steps are of at most three
 * codes (development vectors); as a factor on the mid-frame pitch the
 * development vectors fit 1 - 2.9e-4 (1 - 3.3e-4 scored a residual energy against DVSI's
 * of 0.0564, against 0.0683 at 1). Implemented as 1 - 3.0e-4 for every step.
 * Measured later frame by frame, the offset varies with the two frames' pitch
 * codes; the wave fitted on top of this factor has a constant term of
 * +4.2e-5 (the two constant terms sum to -2.58e-4). */
#define MBE_AMBE2450_MID_PITCH_FACTOR  (1.0f - 3.0e-4f)

/* The half-octave wave of rate 33's mid-frame pitch (docs/dvsi-parity.md,
 * "Rate 33: mid-frame pitch"): between two voiced frames of different pitch,
 * with periods Pa and Pb in samples (2 pi / w0), the rate-33 rendering plays
 * the geometric mean of the two pitches times MBE_AMBE2450_MID_PITCH_FACTOR
 * times 1 + (h(Pa) + h(Pb)) / 2, where, with t = 4 pi log2(P),
 *   h(P) = C0 + C1 cos(t) + S1 sin(t),
 * a wave in the period that repeats every half octave (amplitude 8.44e-4,
 * phase 83.6 degrees). A fit to DVSI's rate-33 decodes of the development
 * vectors: least squares on the change over a frame of the lag between
 * DVSI's decode and ours, read at the frame's two boundaries (where both play
 * the frame's own pitch), on 10,045 frames between two different pitch codes
 * (the selection is in docs/dvsi-parity.md); it takes the median absolute
 * change from 0.0493 to 0.0226 samples there, and from 0.0506 to 0.0266 on
 * the 2,638 such frames of the validation vectors.
 * The law P25's decoder played before its pitch was found to lie on a 2^-19
 * grid was a wave in the same variable with the same period, fitted to the
 * pulse train's drift over runs at one pitch code, which cannot tell a
 * frame's own pitch from its mid-frame pitch (docs/dvsi-parity.md, "P25
 * pitch"). What in DVSI's decoder makes this one is not known. */
#define MBE_AMBE2450_MID_PITCH_WAVE_C0 4.244e-5
#define MBE_AMBE2450_MID_PITCH_WAVE_C1 9.374e-5
#define MBE_AMBE2450_MID_PITCH_WAVE_S1 8.388e-4

/**
 * @brief The mid-frame model between two frames.
 *
 * Its pitch is the geometric mean of the two frames' when both have a voiced
 * band, the voiced frame's when one has (the geometric mean again when the
 * unvoiced frame before a voiced one is marked MBE_FRAME_MARK_PITCH_CARRIES),
 * and the geometric mean when neither has, but in D-STAR's and P25's
 * renderings where one of two frames with no voiced band has a reserved pitch
 * code of the pulsed state (mbe_pulsed.h) that the rendering reads (in P25's
 * this frame's code only while the validity counter that @p prev_mp's state
 * holds allows, and the frame before's only where that frame read it;
 * mbe_synthesize_subframes()): there P25's is this frame's
 * pitch, and D-STAR's the frame before's unless this frame's own-subframe
 * pattern pulses a channel (then this frame's, or the geometric mean if the
 * frame before's pattern does too). The model itself carries no pulsed band:
 * mbe_synthesize_subframes() gives the three models of a frame their pulsed
 * bands. In the rate-33 rendering
 * (MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) the geometric mean of two voiced
 * frames whose pitches differ is multiplied by MBE_AMBE2450_MID_PITCH_FACTOR
 * and by the half-octave wave (MBE_AMBE2450_MID_PITCH_WAVE_C0), in double,
 * rounded to float once; no other case changes.
 * In D-STAR's and P25's renderings (MBE_VOICED_SYNTHESIS_DVSI,
 * MBE_VOICED_SYNTHESIS_DVSI_IMBE) a band's voicing comes from 500 Hz columns
 * (docs/dvsi-parity.md, "D-STAR and P25: mid-frame voicing by columns"): a
 * band in the same state in both frames, each read at its harmonic nearest the
 * band's frequency (l w0 / w0' in double, ties to the upper harmonic), keeps
 * it; one that differs is voiced if and only if M holds at the column of its
 * frequency or at the column of 250 Hz more, M(k) = (Sp(k) and Sc(k)) or
 * (Sp(k) != Sc(k) and Sp(k + 1) != Sc(k + 1)) over the two frames' column
 * samples Sp and Sc, the sample past the last column unvoiced (the
 * samples: D-STAR's 500 k Hz, both columns of each 1 kHz band taking the band's
 * lowest harmonic; P25's from 505 Hz, read at harmonic 1 and at the harmonic
 * at or next above each edge, an isolated voiced last band of a frame of more
 * than 36 harmonics read as unvoiced). A frame with no voiced band counts as
 * unvoiced; a band with no state (a voicing that is neither 0 nor 1) in either
 * frame has none. The columns do not apply, and the current frame's voicing at
 * the band's frequency stands, when the current frame has no voiced band and is
 * not marked MBE_FRAME_MARK_PITCH_ABSENT; no mark of the previous frame is
 * read for them. In the rate-33 rendering a band takes the
 * current frame's voicing at its frequency, or the previous frame's when the
 * current frame has no voiced band, is marked MBE_FRAME_MARK_PITCH_ABSENT
 * (which no rate-33 decoder writes) and the previous frame is voiced at its
 * harmonic at or next above MBE_VOICING_REACH_HZ. With @p late_onset no band is
 * voiced (the caller passes it for a frame marked MBE_FRAME_MARK_LATE_ONSET: in
 * the rate-33 rendering whatever the frame before, in the others after a frame
 * with no voiced band; mbe_synthesize_subframes()).
 * In the rate-33 rendering a current frame
 * marked MBE_FRAME_MARK_MID_COLUMN_0 or MBE_FRAME_MARK_MID_COLUMNS_0_1
 * (V/UV rows 17 and 13) instead voices the bands the 500 Hz column 0, or
 * columns 0 and 1, voice at the mid-frame model's own pitch
 * (mbe_ambe2450_harmonic_voiced()), whatever the two frames' voicing, and
 * leaves the others unvoiced. In any other rendering these marks change
 * nothing.
 * Its magnitudes are the two frames' at its frequency, interpolated in the log
 * domain. Its phases are regenerated from the mean of the two frames' log2Ml
 * at the band's frequency, read as the magnitudes are, each log2Ml used as
 * the decoder stored it: where the magnitudes' interpolation takes the
 * unvoiced factor out of an unvoiced band's Ml, nothing is taken out of or put
 * into its log2Ml, which carries no such factor (an AMBE decoder builds an
 * ordinary voice frame's unvoiced band as Ml = MBE_UNVOICED_MAGNITUDE_FACTOR
 * / sqrt(w0) times 2^log2Ml; Ml then also carries what log2Ml does not, such
 * as rate 33's further silence scale and spectral enhancement); in P25's
 * rendering, whose decoded log2Ml of an unvoiced band is the IMBE magnitude
 * itself, each frame's unvoiced bands are read with the unvoiced factor taken
 * out (mbe_subframe_rendered_own_model()). A current
 * frame marked MBE_FRAME_MARK_PREVIOUS_HELD
 * (a rate-33 silence frame) instead gets the previous frame's own model
 * (mbe_subframe_rendered_own_model()), whatever the two frames' voicing: DVSI's
 * decoder plays it through the silence frame's first 10 ms and crossfades to
 * the silence model over the second. A current frame marked
 * MBE_FRAME_MARK_TONE (a P25 tone sent as voice frames, played as its tone
 * model) gets its own model (mbe_subframe_rendered_own_model()): DVSI's decoder brings
 * the tone's amplitude from the last frame's to its own over the first 10 ms
 * and holds it over the second. Harmonic counts outside 1..56 give a model
 * with no band.
 *
 * @param prev_mp    Previous frame parameters.
 * @param cur_mp     Current frame parameters.
 * @param late_onset Nonzero when the frame's voiced bands start late (no band
 *                   voiced at mid-frame).
 * @param rendering  MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 for rate 33,
 *                   MBE_VOICED_SYNTHESIS_DVSI_IMBE for P25 and
 *                   MBE_VOICED_SYNTHESIS_DVSI_PROVOICE for ProVoice (P25's
 *                   without the pulsed state); any other value renders as
 *                   MBE_VOICED_SYNTHESIS_DVSI (D-STAR) but for the pulsed
 *                   state, which only MBE_VOICED_SYNTHESIS_DVSI and
 *                   MBE_VOICED_SYNTHESIS_DVSI_IMBE have (its mid-frame pitch;
 *                   mbe_synthesize_subframes()).
 * @param model      Receives the model.
 */
void mbe_subframe_mid_model(const mbe_parms* prev_mp, const mbe_parms* cur_mp, int late_onset,
                            enum mbe_voiced_synthesis rendering, struct mbe_subframe_model* model);

/**
 * @brief Synthesize one frame in subframes.
 *
 * Adds the voiced and unvoiced parts to @p aout_buf (the caller clears it,
 * scales and clips the sum). Reads the pulse train (noiseOverlap: the pulses
 * kept and the pulse phase) and the unvoiced overlap (previousUw) of
 * @p prev_mp, when that state is marked as this rendering's, and writes
 * cur_mp's, with the harmonic phases at the frame's end (PSIl, PHIl) and the
 * noise state. The train runs on from a previous frame that has a voiced
 * band, or is marked MBE_FRAME_MARK_PITCH_ABSENT or
 * MBE_FRAME_MARK_PREVIOUS_HELD, or, in the rate-33 rendering only,
 * MBE_FRAME_MARK_MID_COLUMN_0 (V/UV row 17); after any other it starts
 * afresh.
 *
 * A frame whose harmonic count, or the previous frame's, is not 1..56 renders
 * nothing and leaves the state as it is.
 *
 * A frame that starts from silence (@p from_silence, as DVSI's rate-33
 * decoder renders the first frame heard after a mute or a tone) does not play
 * the previous frame's own model, which counts as having no voiced band: its
 * first subframe rises from silence to the mid-frame model, which is the
 * ordinary one between the previous frame's model and this one's (a late
 * onset keeps its meaning for it), save for a frame marked
 * MBE_FRAME_MARK_PREVIOUS_HELD: the previous model was not heard and is not
 * held, so such a frame (a rate-33 silence frame) rises from silence to its
 * own model, which also stands at mid-frame. The pulse train the previous state
 * carries is dropped, whatever its mark, and a new one starts at the frame's
 * start with phase 0; no unvoiced overlap comes from the previous state. The
 * noise generator runs as for any frame.
 *
 * In D-STAR's and P25's renderings (MBE_VOICED_SYNTHESIS_DVSI,
 * MBE_VOICED_SYNTHESIS_DVSI_IMBE) alone, not in ProVoice's
 * (MBE_VOICED_SYNTHESIS_DVSI_PROVOICE), rate 33's, the TIA-102.BABA
 * rendering's (MBE_VOICED_SYNTHESIS_SPEC) or any other value's, a model with
 * no voiced band at the pitch of one of the 56 reserved pitch codes
 * (mbe_pulsed.h) has that code, the frame before's model as this frame's,
 * whatever their marks (docs/dvsi-parity.md, "D-STAR and P25: the pulsed
 * state"). D-STAR's rendering reads every code a model has. P25's reads this
 * frame's code only while the validity counter the frame before left stands
 * at 0 or more (US 7,970,606, FIG. 8), and the frame before's only where
 * that frame read it (MBE_PULSED_READ_SLOT); the counter
 * (MBE_PULSED_COUNT_SLOT, 0 for a fresh stream) stays through a frame whose
 * model has a voiced band, rises by one, to MBE_PULSED_COUNT_MAX at most,
 * through one with no voiced band at the pitch of the default code b0 25 or
 * of one of the 56, falls by one, to MBE_PULSED_COUNT_MIN at least, through
 * one with no voiced band at another P25 code's pitch, and stays through one
 * at a pitch no P25 code decodes to (the repeat limit's reset model). A code that is not read
 * is no code: everything below (the patterns, the emission, the mid-frame
 * pitch) applies to codes that are read. Its unvoiced
 * harmonics in the 500 Hz channels the code's pattern pulses (the channel of
 * the harmonic's frequency or of 250 Hz more), its own-subframe pattern for a
 * frame's own model and its mid-frame pattern for the mid-frame model (none
 * for the model a frame marked MBE_FRAME_MARK_PREVIOUS_HELD or
 * MBE_FRAME_MARK_TONE puts at mid-frame), are pulsed: they play no noise, and
 * a model with a pulsed band emits one pulse MBE_PULSED_DELAY samples after
 * its point unless the model before it emitted one, played through its
 * pulsed bands under its falling crossfade and the next model's under its
 * rising one, within that interval, as the voiced pulse response with the
 * unvoiced factor taken out of the magnitude. The voiced train plays no
 * pulsed band, and is as it was: a model with pulsed bands counts as one with
 * no voiced band for it. Whether the mid-frame model emitted is stored with
 * the train (MBE_PULSED_SLOT), and in P25's rendering whether the code was
 * read and the counter, each read exactly by its bits.
 *
 * In the rate-33 rendering (MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) alone, rate
 * 33's form of the pulsed state (docs/dvsi-parity.md, "Rate 33: the pulsed
 * rows"): a model with no voiced band (after the core's adaptive smoothing)
 * carrying the pulsed-row word of a V/UV row 18..31 and pitch code b0
 * (mbe_adaptive.h) in PHIl[0], on the band-grid model of that code
 * (mbe_ambe2450_pulsed_row_model(), mbe_pulsed.h), is a pulsed-row frame: the
 * row's mid-frame mask applies to the mid-frame model (its own mask for a
 * repeat's held word), its own mask to the frame's own model, and the own mask
 * of the row the frame before stored (MBE_PULSED_ROW_SLOT, read only when
 * @p prev_mp is the grid model of that word's code) to the frame before's own
 * model. Pulsed bands are chosen as above (the channel of the harmonic's
 * frequency or of 250 Hz more, in a model with no voiced band) and play no
 * noise. Each model with a pulsed band emits exactly one pulse, at the time
 * the row and code give (MBE_AMBE2450_PULSE_SLOPE and the constants after it),
 * played by that model alone, as the voiced pulse response through its pulsed
 * bands under its crossfade weight, rising in the interval that ends at its
 * point and falling in the one that starts there; no alternation, no pulse
 * shared between two models. The voiced train, the mid-frame voicing and the
 * mid-frame pitch are as without the word. From silence the frame before's own
 * model has no pulsed band. A frame whose own model has a pulsed band stores
 * its word (not held) for the next; any other stores 0, as every other
 * rendering does.
 *
 * In P25's rendering (MBE_VOICED_SYNTHESIS_DVSI_IMBE, and ProVoice's) the frame is first
 * marked by the voicing it is rendered with (after the adaptive smoothing of
 * the synthesis core): MBE_FRAME_MARK_ONSET_COLUMNS when it has a voiced band
 * and is not a tone model (MBE_FRAME_MARK_TONE), that mark cleared when it
 * has none; mbe_advance_subframes() marks it alike. The mark chooses no band
 * (the columns apply to every voiced model of that rendering, and to a frame
 * after a tone model as after any); it replaces any other mark a voiced model
 * carries.
 *
 * @param aout_buf     160 samples.
 * @param cur_mp       Current frame: enhanced, smoothed amplitudes without render gains.
 * @param prev_mp      Previous frame as synthesized.
 * @param plan         Transform plan, or NULL (no unvoiced part).
 * @param voiced_gain  Render gain of the voiced harmonics.
 * @param from_silence Nonzero for a frame that starts from silence.
 * @param rendering    The codec's subframe rendering (mbe_subframe_mid_model()).
 */
void mbe_synthesize_subframes(float* aout_buf, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                              float voiced_gain, int from_silence, enum mbe_voiced_synthesis rendering);

/**
 * @brief Move a stream's synthesis state one frame on without rendering it.
 *
 * For a frame muted inside the synthesis core: the pulse train and the noise
 * generator run on (the pulses still sounding are kept for the next frame),
 * and no unvoiced overlap is left for it. @p from_silence starts the train as
 * mbe_synthesize_subframes() does, so that the pulse train, the noise
 * generator and the pulsed state (the stored emission, MBE_PULSED_SLOT, in
 * P25's rendering the code read and the validity counter, and in rate 33's the
 * pulsed-row word, MBE_PULSED_ROW_SLOT) end where
 * rendering the frame leaves them; so does @p rendering,
 * which must be the one the stream renders with (and marks the frame as
 * rendering it does). The state is not the one a
 * rendering leaves in full: a rendered frame with unvoiced bands leaves its own
 * subframe's second half in previousUw, an advanced one leaves none.
 */
void mbe_advance_subframes(mbe_parms* cur_mp, const mbe_parms* prev_mp, int from_silence,
                           enum mbe_voiced_synthesis rendering);

#endif /* MBEINT_MBE_SUBFRAME_H */
