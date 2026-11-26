// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Speech synthesis in 10 ms subframes, as DVSI's decoders render it.
 *
 * The codec paths (MBE_VOICED_SYNTHESIS_DVSI and, for rate 33,
 * MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) synthesize a frame as two subframes of
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
 *    the next's across an interval;
 *  - the unvoiced part is mbe_synthesizeUnvoicedSubframes().
 */

#ifndef MBEINT_MBE_SUBFRAME_H
#define MBEINT_MBE_SUBFRAME_H

#include "mbe_adaptive.h"
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
    float phase[57];          /**< [1..L]: regenerated phase of a voiced band (from log2Ml). */
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

/* The rate-33 rendering plays the mid-frame pitch between two voiced frames of
 * different pitch at the geometric mean of their pitches times this
 * (docs/dvsi-parity.md, "Rate 33: mid-frame pitch"). A measured offset, not a
 * known mechanism: DVSI's rate-33 pulse train falls behind ours by 0.013
 * samples a frame for each mid-frame pitch taken between two different pitch
 * codes, and by nothing between equal ones, measured on steps of one to three
 * codes in either direction; as a factor on the mid-frame pitch the
 * development vectors fit 1 - 2.9e-4 (1 - 3.3e-4 scored a residual energy
 * against DVSI's of 0.0564, against 0.0683 at 1). Implemented as 1 - 3.0e-4
 * for every step; larger steps were not measured. */
#define MBE_AMBE2450_MID_PITCH_FACTOR (1.0f - 3.0e-4f)

/**
 * @brief The mid-frame model between two frames.
 *
 * Its pitch is the geometric mean of the two frames' when both have a voiced
 * band, the voiced frame's when one has (the geometric mean again when the
 * unvoiced frame before a voiced one is marked MBE_FRAME_MARK_PITCH_CARRIES),
 * and the geometric mean when neither has. In the rate-33 rendering
 * (MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) the geometric mean of two voiced
 * frames whose pitches differ is multiplied by MBE_AMBE2450_MID_PITCH_FACTOR;
 * no other case changes. A band takes the current frame's
 * voicing at its frequency, or the previous frame's when the current frame has
 * no voiced band, is marked MBE_FRAME_MARK_PITCH_ABSENT and the previous
 * frame's voicing reaches MBE_VOICING_REACH_HZ; with @p late_onset no band is
 * voiced (the voiced bands of a frame so marked start in its own subframe).
 * Its magnitudes are the two frames' at its frequency, interpolated in the log
 * domain. Its phases are regenerated from the mean of the two frames' log2Ml
 * at the band's frequency, read as the magnitudes are, each log2Ml used as
 * the decoder stored it: where the magnitudes' interpolation takes the
 * unvoiced factor out of an unvoiced band's Ml, nothing is taken out of or put
 * into its log2Ml, which carries no such factor (an AMBE decoder builds an
 * ordinary voice frame's unvoiced band as Ml = MBE_UNVOICED_MAGNITUDE_FACTOR
 * / sqrt(w0) times 2^log2Ml; Ml then also carries what log2Ml does not, such
 * as rate 33's further silence scale and spectral enhancement). A current
 * frame marked MBE_FRAME_MARK_PREVIOUS_HELD
 * (a rate-33 silence frame) instead gets the previous frame's own model
 * (mbe_subframe_own_model()), whatever the two frames' voicing: DVSI's
 * decoder plays it through the silence frame's first 10 ms and crossfades to
 * the silence model over the second. Harmonic counts outside 1..56 give a
 * model with no band.
 *
 * @param prev_mp    Previous frame parameters.
 * @param cur_mp     Current frame parameters.
 * @param late_onset Nonzero when the frame's voiced bands start late.
 * @param rendering  MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 for rate 33; any other
 *                   value renders as MBE_VOICED_SYNTHESIS_DVSI.
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
 * noise state.
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
 * mbe_synthesize_subframes() does, so that the pulse train and the noise
 * generator end where rendering the frame leaves them; so does @p rendering,
 * which must be the one the stream renders with. The state is not the one a
 * rendering leaves in full: a rendered frame with unvoiced bands leaves its own
 * subframe's second half in previousUw, an advanced one leaves none.
 */
void mbe_advance_subframes(mbe_parms* cur_mp, const mbe_parms* prev_mp, int from_silence,
                           enum mbe_voiced_synthesis rendering);

#endif /* MBEINT_MBE_SUBFRAME_H */
