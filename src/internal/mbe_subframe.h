// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Speech synthesis in 10 ms subframes, as DVSI's decoders render it.
 *
 * The codec paths (MBE_VOICED_SYNTHESIS_DVSI) synthesize a frame as two
 * subframes of MBE_SUBFRAME samples (docs/dvsi-parity.md, "Synthesis in
 * subframes"):
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
 *    the next's across an interval;
 *  - the unvoiced part is mbe_synthesizeUnvoicedSubframes().
 */

#ifndef MBEINT_MBE_SUBFRAME_H
#define MBEINT_MBE_SUBFRAME_H

#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

/* Samples from a pitch pulse to the peak of the response it plays, and the
 * length of each raised-cosine end of the response's window. */
#define MBE_PULSE_PEAK                16

/* An unvoiced band's amplitude is its magnitude times this over sqrt(w0); the
 * mid-frame model's magnitudes are interpolated with that factor taken out,
 * so that a band voiced in one frame and unvoiced in the other moves between
 * like quantities. */
#define MBE_UNVOICED_MAGNITUDE_FACTOR 0.2046f

/* Pulses a stream keeps from one frame to the next (those whose responses can
 * still sound), and where: element 0 of noiseOverlap holds their count, the
 * elements after it their times in samples before the frame's end, and the
 * one after those the pulse phase reached at the frame's end, in cycles. The
 * TIA-102.BABA rendering keeps its noise overlap there instead. */
#define MBE_PULSE_HISTORY             32

/** One subframe's model. */
struct mbe_subframe_model {
    float w0;                 /**< Fundamental, radians per sample; 0 for none. */
    int L;                    /**< Harmonic count. */
    int any_voiced;           /**< Whether a band is voiced. */
    unsigned char voiced[57]; /**< [1..L]: the band is voiced (and below Nyquist). */
    float magnitude[57];      /**< [1..L]: the band's amplitude, as mbe_parms.Ml holds it. */
    float phase[57];          /**< [1..L]: regenerated phase of a voiced band. */
};

/**
 * @brief The model a frame's parameters give at its end.
 * @param mp    Frame parameters.
 * @param model Receives the model; it has no band when w0 is not a pitch
 *              (positive, finite and not subnormal: a subnormal value is zero
 *              to the arithmetic where the caller's floating-point environment
 *              flushes denormals) or the harmonic count is not 1..56.
 */
void mbe_subframe_own_model(const mbe_parms* mp, struct mbe_subframe_model* model);

/**
 * @brief The mid-frame model between two frames.
 *
 * Its pitch is the geometric mean of the two frames' when both have a voiced
 * band, the voiced frame's when one has (the geometric mean again when the
 * unvoiced frame before a voiced one is marked MBE_FRAME_MARK_PITCH_CARRIES),
 * and the geometric mean when neither has. A band takes the current frame's
 * voicing at its frequency, or the previous frame's when the current frame has
 * no voiced band, is marked MBE_FRAME_MARK_PITCH_ABSENT and the previous
 * frame's voicing reaches MBE_VOICING_REACH_HZ; with @p late_onset no band is
 * voiced (the voiced bands of a frame so marked start in its own subframe).
 * Its magnitudes are the two frames' at its frequency, interpolated in the log
 * domain. A current frame marked MBE_FRAME_MARK_PREVIOUS_HELD (a rate-33
 * silence frame) instead gets the previous frame's own model
 * (mbe_subframe_own_model()), whatever the two frames' voicing: DVSI's
 * decoder plays it through the silence frame's first 10 ms and crossfades to
 * the silence model over the second. Harmonic counts outside 1..56 give a
 * model with no band.
 *
 * @param prev_mp    Previous frame parameters.
 * @param cur_mp     Current frame parameters.
 * @param late_onset Nonzero when the frame's voiced bands start late.
 * @param model      Receives the model.
 */
void mbe_subframe_mid_model(const mbe_parms* prev_mp, const mbe_parms* cur_mp, int late_onset,
                            struct mbe_subframe_model* model);

/**
 * @brief Synthesize one frame in subframes.
 *
 * Adds the voiced and unvoiced parts to @p aout_buf (the caller clears it,
 * scales and clips the sum). Reads the pulse train (noiseOverlap: the pulses
 * kept and the pulse phase) and the unvoiced overlap (previousUw) of
 * @p prev_mp, when that state is marked as this rendering's, and writes
 * cur_mp's, with the harmonic phases at the frame's end (PSIl, PHIl) and the
 * noise state.
 *
 * A frame whose harmonic count, or the previous frame's, is not 1..56 renders
 * nothing and leaves the state as it is.
 *
 * @param aout_buf     160 samples.
 * @param cur_mp       Current frame: enhanced, smoothed amplitudes without render gains.
 * @param prev_mp      Previous frame as synthesized.
 * @param plan         Transform plan, or NULL (no unvoiced part).
 * @param voiced_gain  Render gain of the voiced harmonics.
 */
void mbe_synthesize_subframes(float* aout_buf, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                              float voiced_gain);

/**
 * @brief Move a stream's synthesis state one frame on without rendering it.
 *
 * For a frame muted inside the synthesis core: the pulse train and the noise
 * generator run on (the pulses still sounding are kept for the next frame),
 * and no unvoiced overlap is left for it.
 */
void mbe_advance_subframes(mbe_parms* cur_mp, const mbe_parms* prev_mp);

#endif /* MBEINT_MBE_SUBFRAME_H */
