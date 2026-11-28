// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 *
 * Internal header for adaptive smoothing algorithms.
 * Implements TIA-102.BABA equations 111-116 for error-based audio smoothing.
 */

#ifndef MBEINT_MBE_ADAPTIVE_H
#define MBEINT_MBE_ADAPTIVE_H

#include <stdint.h>

#include "mbelib-neo/mbelib.h"

/* Algorithm constants from TIA-102.BABA sections 8 and 9 */

/** Default local energy value (Algorithm #111). */
#define MBE_DEFAULT_LOCAL_ENERGY        75000.0f

/** Minimum local energy threshold (Algorithm #111). */
#define MBE_MIN_LOCAL_ENERGY            10000.0f

/** Energy smoothing coefficient alpha (Algorithm #111). */
#define MBE_ENERGY_SMOOTH_ALPHA         0.95f

/** Energy smoothing coefficient beta (Algorithm #111). */
#define MBE_ENERGY_SMOOTH_BETA          0.05f

/** Default amplitude threshold (Algorithm #115). */
#define MBE_DEFAULT_AMPLITUDE_THRESHOLD 20480

/** Error rate threshold for smoothing entry (Algorithm #112). */
#define MBE_ERROR_THRESHOLD_ENTRY       0.0125f

/** Low error rate threshold (Algorithm #112). */
#define MBE_ERROR_THRESHOLD_LOW         0.005f

/** Adaptive gain constant (Algorithm #112). */
#define MBE_ADAPTIVE_GAIN               45.255f

/** Adaptive exponent constant (Algorithm #112). */
#define MBE_ADAPTIVE_EXPONENT           277.26f

/** Alternative adaptive multiplier (Algorithm #112). */
#define MBE_ADAPTIVE_ALT                1.414f

/** Amplitude penalty per error (Algorithm #115). */
#define MBE_AMPLITUDE_PENALTY_PER_ERROR 300

/** Amplitude base constant (Algorithm #115). */
#define MBE_AMPLITUDE_BASE              6000

/**
 * @brief Apply adaptive smoothing with caller-supplied pre-enhancement RM0.
 *
 * TIA-102.BABA equation 111 uses local energy derived from pre-enhanced spectral
 * amplitudes. Internal decode/synthesize paths pass that value explicitly so
 * the handoff cannot retain caller stack addresses.
 *
 * @param cur_mp Current frame parameters (modified in-place).
 * @param prev_mp Previous frame parameters (for local energy).
 * @param rm0 Sum of squared pre-enhancement amplitudes.
 */
void mbe_applyAdaptiveSmoothingWithRm0(mbe_parms* cur_mp, const mbe_parms* prev_mp, float rm0);

/**
 * Whether adaptive smoothing scales a frame's amplitudes down to the threshold
 * of TIA-102.BABA eqs 115-116; chosen by the codec path.
 */
enum mbe_amplitude_limit {
    /** Eq 116 as specified: D-STAR and the public synthesis API. */
    MBE_AMPLITUDE_LIMIT_SPEC,
    /** No eq 116 scaling: the IMBE (P25 and ProVoice) and rate-33 paths. On identical bits
     * DVSI's decoders play the frames it would scale (1.7% of rate 33's voice
     * frames, 1% of P25's) 0.0-0.1 dB from their unscaled level, where eq 116
     * left them 2.9 and 1.0 dB quieter (docs/dvsi-parity.md, "Amplitude
     * limit"). */
    MBE_AMPLITUDE_LIMIT_NONE,
};

/* The amplitude sum MBE_AMPLITUDE_LIMIT_NONE still scales down to: a numeric
 * guard only, about 94 dB over eq 115's 20480 and far past full-scale output,
 * so an extreme (finite) prediction history cannot overflow the synthesis. */
#define MBE_AMPLITUDE_SAFETY_LIMIT 1.0e9

/** @brief mbe_applyAdaptiveSmoothingWithRm0() with the eq 116 scaling chosen by @p limit. */
void mbe_applyAdaptiveSmoothingWithRm0Limit(mbe_parms* cur_mp, const mbe_parms* prev_mp, float rm0,
                                            enum mbe_amplitude_limit limit);

/**
 * @brief Apply spectral amplitude enhancement and return pre-enhancement RM0.
 * @param cur_mp In/out parameter set to enhance.
 * @param high_band_gain Height of the high-band shelf applied afterwards (1 for none).
 * @return Sum of squared amplitudes before enhancement.
 */
float mbe_spectralAmpEnhanceWithRm0(mbe_parms* cur_mp, float high_band_gain);

/*
 * High-band shelf heights per codec (linear amplitude gain above 3 kHz; see
 * mbe_apply_high_band_gain), each matched to DVSI's AMBE-3000 decoder for that
 * mode on its test vectors. ProVoice shares the IMBE synthesis path.
 */
#define MBE_HIGH_BAND_GAIN_IMBE     1.25892541f /* +2.0 dB */
#define MBE_HIGH_BAND_GAIN_AMBE2450 1.33352143f /* +2.5 dB */
#define MBE_HIGH_BAND_GAIN_AMBE2400 1.33352143f /* +2.5 dB */

/**
 * @brief Low-band roll-off synthesis renders at a harmonic frequency.
 * @param hz Harmonic frequency in Hz.
 * @return Amplitude gain of a second-order high-pass at 100 Hz.
 */
float mbe_low_band_gain(float hz);

/**
 * @brief Phase (radians) of the same high-pass at a harmonic frequency.
 * @param hz Harmonic frequency in Hz.
 * @return arg H: pi toward DC, pi/2 at 100 Hz, near 0 far above; 0 for hz <= 0.
 */
float mbe_low_band_phase(float hz);

/* Restore the library default voice model while retaining stream continuity. */
void mbe_resetImbeModel(mbe_parms* mp);

/** Noise a frame muted inside the synthesis core is replaced with; chosen by the codec path. */
enum mbe_mute_noise {
    MBE_MUTE_NOISE_COMFORT, /**< Library comfort noise: D-STAR, ProVoice and the public synthesis API. */
    MBE_MUTE_NOISE_SPEC,    /**< MBE_SPEC_MUTE_NOISE_AMPLITUDE: the P25 IMBE path (the AMBE 3600x2450
                               path passes it too, but mutes before synthesis and never reaches the
                               core's mute: it plays MBE_AMBE_MUTE_NOISE_AMPLITUDE itself). */
};

/** How a frame's model is rendered; chosen by the codec path. */
enum mbe_voiced_synthesis {
    /** TIA-102.BABA: at steady pitch the amplitudes of harmonics 1..7 move
     * linearly over the frame and the others crossfade with the window Ws;
     * the unvoiced part is one transform per frame (eq 117-126). The public
     * mbe_synthesizeSpeechf() renders this way. */
    MBE_VOICED_SYNTHESIS_SPEC,
    /** The codec paths, as DVSI's AMBE-3000 decoders: two 10 ms subframes
     * per frame, the voiced part a train of pitch pulses and the unvoiced
     * part DVSI's noise (mbe_subframe.h; docs/dvsi-parity.md, "Synthesis in
     * subframes"), the mid-frame model voiced by D-STAR's 500 Hz columns
     * (mbe_subframe_mid_model(); docs/dvsi-parity.md, "D-STAR and P25:
     * mid-frame voicing by columns"), and D-STAR's pulsed state on a frame
     * with no voiced band at b0 64..119 (mbe_pulsed.h; docs/dvsi-parity.md,
     * "D-STAR and P25: the pulsed state"). D-STAR, and with the three below
     * the base of every codec path's rendering. */
    MBE_VOICED_SYNTHESIS_DVSI,
    /** MBE_VOICED_SYNTHESIS_DVSI as DVSI's rate-33 (AMBE 3600x2450) decoder
     * renders it, whose mid-frame pitch between two voiced frames of
     * different pitch is their geometric mean times a factor and a
     * half-octave wave in their periods, within -0.11% and +0.06% of it
     * (mbe_subframe_mid_model(); docs/dvsi-parity.md, "Rate 33: mid-frame
     * pitch"), and which honours the marks of the rate-33 rows 13, 15 and 17
     * (MBE_FRAME_MARK_MID_COLUMNS_0_1, MBE_FRAME_MARK_LATE_ONSET whatever the
     * frame before, MBE_FRAME_MARK_MID_COLUMN_0, which on the frame before
     * also carries the pulse train). Every rate-33 synthesis call. */
    MBE_VOICED_SYNTHESIS_DVSI_AMBE2450,
    /** MBE_VOICED_SYNTHESIS_DVSI as DVSI's P25 (IMBE 7200x4400) decoder
     * renders it: the mid-frame model voiced by P25's 500 Hz columns, from
     * 505 Hz (mbe_subframe_mid_model()), after a tone model as after any
     * (a tone model itself, MBE_FRAME_MARK_TONE, keeps its own model at
     * mid-frame, and a frame with no voiced band not marked
     * MBE_FRAME_MARK_PITCH_ABSENT takes its own voicing); the rendering marks every voiced model it
     * renders that is not a tone model MBE_FRAME_MARK_ONSET_COLUMNS, after
     * adaptive smoothing (mbe_synthesize_subframes(),
     * mbe_advance_subframes()), and P25's pulsed state on a frame with no
     * voiced band at b0 128..183, read while the stream's validity counter
     * allows (mbe_pulsed.h; docs/dvsi-parity.md, "D-STAR and P25: the pulsed
     * state"). The P25 data path. */
    MBE_VOICED_SYNTHESIS_DVSI_IMBE,
    /** The ProVoice (IMBE 7100x4400) frames of the shared IMBE 4400 data
     * path, which has no DVSI reference for them (MBE_PROCESS_FLAG_PROVOICE):
     * P25's rendering as it stood before the pulsed state and before the four
     * rules of DVSI's P25 decoder (docs/dvsi-parity.md), so that ProVoice
     * plays as it did. It shares P25's mid-frame voicing by the columns (P25's
     * samples, of the decoded voicing), the marking of the model it renders
     * and the marks it honours; it has no pulsed state and no DVSI low band
     * (the second-order filter at 100 Hz, mbe_low_band.h), its mid-frame pitch
     * between two voiced models is the geometric mean, its phase kernel reads
     * every band as decoded and its own models keep their decoded voicing. */
    MBE_VOICED_SYNTHESIS_DVSI_PROVOICE,
};

/** How a frame's rendering starts; chosen per call, it has a meaning for
 * the subframe renderings (MBE_VOICED_SYNTHESIS_DVSI,
 * MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_VOICED_SYNTHESIS_DVSI_IMBE and
 * MBE_VOICED_SYNTHESIS_DVSI_PROVOICE) only. */
enum mbe_synthesis_start {
    /** From the previous frame as it was synthesized: its model, the pulse
     * train it carries and its unvoiced overlap. */
    MBE_SYNTHESIS_CONTINUES,
    /** From silence, as DVSI's rate-33 decoder renders the first frame heard
     * after a mute or a tone (docs/dvsi-parity.md, "Rate 33 frames after a
     * mute or a tone"): the previous frame's own model is not played (it
     * counts as having no voiced band), the pulse train starts anew at the
     * frame's start, and no unvoiced overlap comes from before. The
     * mid-frame model is the ordinary one between the previous model and this
     * frame's (mbe_subframe_mid_model()), save that a frame marked
     * MBE_FRAME_MARK_PREVIOUS_HELD (a rate-33 silence frame) does not hold a
     * previous model that was not heard: it rises from silence to its own. */
    MBE_SYNTHESIS_FROM_SILENCE,
};

/** Whether the synthesis core limits the frame it renders; chosen per call. */
enum mbe_output_limit {
    /** Limited to the float output limit (mbe_clipFloatBuffer()), as every
     * frame the core renders otherwise is. */
    MBE_OUTPUT_LIMITED,
    /** Left unlimited: the caller mixes another signal into the frame (a
     * tone's tail, a crossfade with a tone, mute noise) and limits the sum
     * once, after the mix. */
    MBE_OUTPUT_UNLIMITED,
};

/* The onset delay mbe_set_voicing_onset_delay() can give a model (samples):
 * its voiced bands then start in the frame's second subframe
 * (MBE_FRAME_MARK_LATE_ONSET: in the rate-33 rendering whatever the frame
 * before, in the others after a frame with no voiced band). */
#define MBE_VOICING_ONSET_LATE 80

/**
 * @brief Mark a decoded model whose voiced harmonics start late:
 * @p late nonzero gives it MBE_VOICING_ONSET_LATE samples, zero none. The
 * mark lives in the otherwise unused element 0 of PHIl.
 */
void mbe_set_voicing_onset_delay(mbe_parms* mp, int late);

/** @brief The onset delay mbe_set_voicing_onset_delay() gave @p mp (0 if none). */
int mbe_voicing_onset_delay(const mbe_parms* mp);

/* The first edge of P25's 500 Hz columns, in which P25's rendering reads the
 * mid-frame voicing (mbe_subframe_mid_model()): on DVSI's P25 decodes a first
 * band of three voiced harmonics reaching 502.8 Hz at the pitch played (pitch
 * code 56; 502.6 Hz at eq 46's nominal pitch) covers no column and starts at
 * mid-frame after a frame with no voiced band, one reaching 508.1 Hz (code 55;
 * 507.9 Hz nominal) covers the first and starts with the frame. The rate-33
 * rendering, whose decoder writes no MBE_FRAME_MARK_PITCH_ABSENT, keeps the
 * hold every rendering applied before the columns: the voicing of the frame
 * before a frame so marked held through that frame's first half when a voiced
 * harmonic lies at or next above this frequency. */
#define MBE_VOICING_REACH_HZ 505.0f

/**
 * What a codec's decoder tells the synthesizer about a frame beyond its model
 * (docs/dvsi-parity.md, "Voiced onsets and offsets"). A model carries one mark
 * at most, in the otherwise unused element 0 of PHIl, stored as the mark's
 * value in float. The marks of the mid-frame voicing (MID_COLUMN_0,
 * MID_COLUMNS_0_1, and LATE_ONSET's reading whatever the frame before) act
 * only in the rendering of the codec whose decoder writes them: in D-STAR's
 * rendering 5, 6 and 7 play as none and LATE_ONSET acts as it did before them,
 * only after a frame with no voiced band; in rate 33's 5 plays as none. Other
 * marks in a caller's state are read as marks in every rendering:
 * PREVIOUS_HELD and TONE on the current model give the mid-frame model before
 * any rendering is told apart (the previous model, or the frame's own);
 * PITCH_ABSENT on a current model with no voiced band gives the columns' hold
 * in D-STAR's and P25's renderings and the 505 Hz hold in rate 33's;
 * PITCH_CARRIES, PITCH_ABSENT and PREVIOUS_HELD on the previous model act on
 * the mid-frame pitch and the pulse train; in rate 33's rendering
 * MID_COLUMN_0 on the previous model carries the pulse train too (in the
 * others it plays there as no mark does). These marks of the previous model
 * are read as the caller's state holds them, unchecked; each decides only a
 * pitch or whether the train that state stored is loaded, which the loader
 * bounds. No mark of the previous model
 * chooses the columns' bands. P25's rendering gives the model it renders a
 * mark of its own (MBE_FRAME_MARK_ONSET_COLUMNS on any model with a voiced
 * band that is not a tone model), so there a caller's mark on a voiced model,
 * LATE_ONSET and PREVIOUS_HELD included, is overwritten and has no effect; on a
 * model with no voiced band PITCH_ABSENT, PREVIOUS_HELD and PITCH_CARRIES still
 * act, and LATE_ONSET, 6 and 7 change nothing.
 */
enum mbe_frame_mark {
    MBE_FRAME_MARK_NONE = 0,
    /* No voiced band, and the pitch enters the mid-frame pitch of a voiced
     * frame after it. */
    MBE_FRAME_MARK_PITCH_CARRIES = 1,
    /* No voiced band, at a pitch code other than those DVSI's encoders keep
     * for such frames: the voicing of the frame before it lasts to mid-frame,
     * in D-STAR's and P25's renderings in the bands whose 500 Hz column, or
     * the column 250 Hz higher, that frame covers (mbe_subframe_mid_model()),
     * the pulse train running on through the frame
     * when a band is held and starting afresh when none is (a P25 frame voiced
     * in isolated bands only, b0 97 word 10000000000, covers no column); in the
     * rate-33 rendering all of it when it is voiced at its harmonic at or next
     * above MBE_VOICING_REACH_HZ. */
    MBE_FRAME_MARK_PITCH_ABSENT = 2,
    /* No voiced band, and the previous frame's model, as it stands at that
     * frame's end, stands at this frame's mid-frame too: DVSI's rate-33
     * decoder plays it through the first 10 ms of a silence frame and
     * crossfades to the silence model over the second, the pulse train
     * running on through the frame. */
    MBE_FRAME_MARK_PREVIOUS_HELD = 3,
    /* A P25 tone sent as voice frames, played as DVSI's decoder plays it: the
     * model is the tone model (imbe7200x4400.c), and the frame's own model
     * stands at mid-frame too, so the tone's amplitude moves from the last
     * frame's to its own over the first 10 ms and holds over the second. Only
     * imbe4400_play_tone() sets it, on the frame it plays; no decoder does,
     * and the P25 and ProVoice data path clears it on any other frame. */
    MBE_FRAME_MARK_TONE = 4,
    /* A P25 or ProVoice voice frame with a voiced band, as decoded
     * (imbe7200x4400.c) and as rendered: P25's rendering marks the model it
     * renders so whenever, after adaptive smoothing, it has a voiced band and
     * is not a tone model, whatever its decoded mark (smoothing can voice
     * harmonics of a frame decoded with none), and clears the mark from one
     * with no voiced band. Since the mid-frame columns change it chooses no
     * band (mark_rendered_model() still reads it, to clear it from a model with
     * no voiced band): P25's rendering reads the mid-frame voicing in
     * 500 Hz columns, whatever the current model's mark, for every model but a
     * tone model, a model held by MBE_FRAME_MARK_PREVIOUS_HELD on one with no
     * voiced band, and one with no voiced band not marked
     * MBE_FRAME_MARK_PITCH_ABSENT (columns ([0, 505), [505, 1005), ... [3005, 3505) Hz and one
     * above; docs/dvsi-parity.md, "D-STAR and P25: mid-frame voicing by
     * columns"), and the marking only replaces any other
     * mark a voiced model carries. (From PR #198 to that change the mark
     * gave the onset rule of docs/dvsi-parity.md, "P25: the bands that start
     * with the frame", which the columns' rule contains.) */
    MBE_FRAME_MARK_ONSET_COLUMNS = 5,
    /* A rate-33 voice frame of V/UV row 17: the mid-frame model voices the
     * 500 Hz column 0 (through mbe_ambe2450_harmonic_voiced(), at its own
     * harmonics) and no other, whatever the two frames' voicing; the frame's
     * own voicing is Annex B's row 17, all unvoiced. On the previous model
     * the pulse train runs on from the frame, as from one marked
     * MBE_FRAME_MARK_PITCH_ABSENT or MBE_FRAME_MARK_PREVIOUS_HELD: the train
     * the row-17 frame stored, which ran through its voiced first half and
     * out of its mid-frame model. Read only in the rate-33 rendering
     * (MBE_VOICED_SYNTHESIS_DVSI_AMBE2450). */
    MBE_FRAME_MARK_MID_COLUMN_0 = 6,
    /* A rate-33 voice frame of V/UV row 13: the mid-frame model voices the
     * 500 Hz columns 0 and 1 (as MBE_FRAME_MARK_MID_COLUMN_0 does column 0)
     * and no other; the frame's own voicing is row 13's, columns 0 to 2.
     * Read only in the rate-33 rendering. */
    MBE_FRAME_MARK_MID_COLUMNS_0_1 = 7,
    /* In the rate-33 rendering, the mid-frame model voices no band, whatever
     * the frame before: the voiced bands start MBE_VOICING_ONSET_LATE samples
     * late (V/UV row 15; after a voiced frame the pulse train runs on through
     * the unvoiced mid-frame, at the previous frame's pitch into it and this
     * frame's out of it). The rate-33 rows 13, 15 and 17 give each half of
     * the frame its own voicing (docs/dvsi-parity.md, "Rate 33: rows 13, 15
     * and 17"); a lost frame that repeats such a frame plays the repeated
     * model's own voicing through (the rate-33 path clears these marks on a
     * repeat). In D-STAR's rendering the mark keeps its earlier meaning: no
     * band voiced at mid-frame only after a frame with no voiced band (no
     * decoder of D-STAR writes it). In P25's rendering it is overwritten on a
     * voiced model (MBE_FRAME_MARK_ONSET_COLUMNS) and changes nothing on one
     * with no voiced band. */
    MBE_FRAME_MARK_LATE_ONSET = MBE_VOICING_ONSET_LATE,
};

/** @brief Give @p mp a mark, replacing any it had. */
void mbe_set_frame_mark(mbe_parms* mp, enum mbe_frame_mark mark);

/** @brief The mark of @p mp; MBE_FRAME_MARK_NONE for any other stored value. */
enum mbe_frame_mark mbe_frame_mark(const mbe_parms* mp);

/*
 * The word a rate-33 voice frame of V/UV rows 18 to 31 carries in PHIl[0] in
 * place of a mark (docs/dvsi-parity.md, "Rate 33: the pulsed rows"): its row,
 * its pitch code b0 (0..119) and whether a lost frame's repeat carries it
 * (held), as the float MBE_PULSED_ROW_WORD_BASE + MBE_PULSED_ROW_WORD_HELD held +
 * MBE_PULSED_ROW_WORD_ROW (row - MBE_PULSED_ROW_FIRST) + b0, every one a
 * whole number a float holds exactly: 2,640 values, the 1,680 words of rows 18
 * to 31 (1024 to 2807), which the decoder writes, and the 960 held words of
 * rows 18 to 25 (3072 to 4087), which the repeat loader writes; a held form of
 * rows 26 to 31 is written by nothing and is no word. It is no enum
 * mbe_frame_mark value: mbe_frame_mark() reads it as MBE_FRAME_MARK_NONE, so
 * that every reader of marks treats such a frame as it did before the word
 * (rows 18 to 31 carried no mark). The rate-33 decoder writes it
 * (ambe2450_decode_vuv()). It is read in two places, each against the model it
 * sits in (the band-grid model of the word's code, mbe_pulsed.h): the repeat
 * loader (ambe2450_load_repeat()), which turns it into the held form or clears
 * it, and the rate-33 rendering (mbe_subframe.c), which also requires no voiced
 * band. The same value without the held part is what that rendering stores in
 * the stream's subframe state for the next frame (MBE_PULSED_ROW_SLOT,
 * mbe_subframe.h).
 */
#define MBE_PULSED_ROW_FIRST     18
#define MBE_PULSED_ROW_LAST      31
/* Rows 18 to 25 pulse their own subframe (26 to 31 only the mid-frame one). */
#define MBE_PULSED_ROW_OWN_LAST  25
#define MBE_PULSED_ROW_CODES     120
#define MBE_PULSED_ROW_WORD_BASE 1024
#define MBE_PULSED_ROW_WORD_HELD 2048
#define MBE_PULSED_ROW_WORD_ROW  128

/** A rate-33 pulsed-row word's fields: row 18..31, b0 0..119, held 0 or 1. */
struct mbe_pulsed_row {
    int row;
    int b0;
    int held;
};

/** @brief The word of @p row (18..31), @p b0 (0..119) and @p held as the float stored. */
float mbe_pulsed_row_value(const struct mbe_pulsed_row* word);

/**
 * @brief Read a pulsed-row word from a stored float, by its bit pattern.
 *
 * Only a value a writer gives (a whole number whose row is 18..31 and whose
 * code is 0..119, held only for rows 18..25: 2,640 values from 1024 to 4087)
 * is a word; anything else (not a number, an infinity, a negative value, a
 * fraction, a neighbouring float, a code of 120 to 127, a held form of rows 26
 * to 31) is none.
 * @return 1 with @p word set for a word, 0 for none.
 */
int mbe_pulsed_row_read(const float* stored, struct mbe_pulsed_row* word);

/** @brief Give @p mp a pulsed-row word in PHIl[0], replacing any mark it had. */
void mbe_set_pulsed_row_word(mbe_parms* mp, const struct mbe_pulsed_row* word);

/** @brief The pulsed-row word of @p mp (mbe_pulsed_row_read() on PHIl[0]). */
int mbe_pulsed_row_word(const mbe_parms* mp, struct mbe_pulsed_row* word);

/**
 * @brief Regenerate harmonic phases from a spectral envelope (US 5,701,390).
 * @param Ml  Amplitudes, [1..L] read.
 * @param L   Harmonic count, 1..56.
 * @param phi Receives the phases of harmonics 1..56 (those past L from the
 *            envelope's extrapolation); element 0 is not written.
 */
void mbe_regenerate_phases(const float Ml[57], int L, float phi[57]);

/** log2 of the smallest amplitude the phase kernel tells apart (1e-6). */
#define MBE_PHASE_LOG2_FLOOR (-19.931568f)

/**
 * @brief Regenerate harmonic phases from log2 amplitudes (US 5,701,390).
 *
 * The kernel of mbe_regenerate_phases() applied to log2 amplitudes given
 * directly, as the codec paths apply it to a frame's decoded log2Ml
 * (docs/dvsi-parity.md, "Phase regeneration"); a value below
 * MBE_PHASE_LOG2_FLOOR counts as that floor.
 * @param log2_amplitude Finite log2 amplitudes, [1..L] read.
 * @param L              Harmonic count, 1..56.
 * @param phi            Receives the phases of harmonics 1..56 (those past L
 *                       from the envelope's extrapolation); element 0 is not
 *                       written.
 */
void mbe_regenerate_phases_log2(const float log2_amplitude[57], int L, float phi[57]);

/**
 * @brief Synthesize speech using a captured pre-enhancement RM0 value.
 * @param aout_buf Output buffer of 160 float samples.
 * @param cur_mp Current frame parameters.
 * @param prev_mp Previous enhanced frame parameters.
 * @param rm0 Sum of squared pre-enhancement amplitudes.
 * @param mute_noise Noise output if the frame is muted.
 * @param voiced How the frame is rendered.
 * @param limit Whether adaptive smoothing applies eq 116's amplitude scaling.
 */
void mbe_synthesizeSpeechWithPreEnhRm0f(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, float rm0,
                                        enum mbe_mute_noise mute_noise, enum mbe_voiced_synthesis voiced,
                                        enum mbe_amplitude_limit limit);

/**
 * @brief mbe_synthesizeSpeechWithPreEnhRm0f() with the way its rendering starts
 *        and whether its output is limited.
 * @param start  How the frame's rendering starts (the subframe renderings only;
 *               a frame the core mutes moves its state on the same way).
 * @param output Whether the rendered frame is limited (a frame the core mutes
 *               is noise either way). The state the call leaves is the same.
 */
void mbe_synthesizeSpeechWithPreEnhRm0Startf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, float rm0,
                                             enum mbe_mute_noise mute_noise, enum mbe_voiced_synthesis voiced,
                                             enum mbe_amplitude_limit limit, enum mbe_synthesis_start start,
                                             enum mbe_output_limit output);

/**
 * @brief Synthesize a frame-repeat model as-is, without adaptive smoothing.
 *
 * TIA-102.BABA-1 5.6 synthesizes the repeated model parameters unchanged;
 * smoothing them again could change their voicing and amplitudes.
 *
 * @param aout_buf Output buffer of 160 float samples.
 * @param cur_mp Repeated (already enhanced) frame parameters.
 * @param prev_mp Previous enhanced frame parameters.
 * @param voiced How the frame is rendered.
 */
void mbe_synthesizeRepeatedSpeechf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp,
                                   enum mbe_voiced_synthesis voiced);

/**
 * @brief mbe_synthesizeRepeatedSpeechf() with the way its rendering starts and
 *        whether its output is limited.
 * @param start  How the frame's rendering starts (the subframe renderings only).
 * @param output Whether the rendered frame is limited.
 */
void mbe_synthesizeRepeatedSpeechStartf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp,
                                        enum mbe_voiced_synthesis voiced, enum mbe_synthesis_start start,
                                        enum mbe_output_limit output);

/**
 * @brief Render 160 samples of an AMBE 3600x2450 tone, continuing the tone phases in @p mp.
 *
 * mbe_synthesizeTonef() without the frame's bits, for a tone known by its
 * fields. Both play the frequencies DVSI's decoders play for the ID
 * (mbe_tone_played_steps(), mbe_tone.h), not Table 9's nominal ones. An ID
 * with no tone renders silence. Limited, it equals mbe_synthesizeTonef() bit
 * for bit; unlimited, it is the same samples before
 * their limit, for a caller that fades or mixes the tone and limits the
 * result. The oscillator phases it leaves are the same either way.
 *
 * @param aout_buf Output buffer of 160 float samples.
 * @param mp Parameter set carrying the tone phases (updated).
 * @param id Tone ID (TIA-102.BABA-1 Table 9).
 * @param ad Tone amplitude field, 0..127.
 * @param output Whether the samples are limited.
 */
void mbe_renderToneAmbe2450f(float* aout_buf, mbe_parms* mp, int id, int ad, enum mbe_output_limit output);

/**
 * @brief Set the tone phases in @p mp to start an AMBE+2 tone as DVSI's decoders
 * start a tone after a frame with no voiced band: a fresh pulse train from the
 * frame's start (mbelib.c, tone_start_phase()). The codec paths call it for a
 * tone after any frame that was not a tone frame; after a voiced frame DVSI's
 * decoders continue that frame's pulse train instead (not modelled). The next
 * mbe_renderToneAmbe2450f() of that ID then plays the tone from there. An ID
 * with no tone leaves them as they are.
 */
void mbe_startToneAmbe2450(mbe_parms* mp, int id);

/**
 * @brief Render 160 samples of a D-STAR tone, continuing the tone phases in @p mp.
 *
 * mbe_synthesizeTonefdstar() without the bit validation, for a tone known by
 * its index and volume. An index with no D-STAR tone renders silence.
 *
 * @param aout_buf Output buffer of 160 float samples.
 * @param mp Parameter set carrying the tone phases (updated).
 * @param index D-STAR tone index.
 * @param volume D-STAR tone volume, 0..255.
 */
void mbe_renderToneDstarf(float* aout_buf, mbe_parms* mp, int index, int volume);

/**
 * @brief mbe_renderToneDstarf() with the tone's level moving from @p from_volume
 * to @p volume over the first @p ramp samples, before the output is clipped.
 *
 * @param ramp Samples of the level ramp, 1..160 (anything else: none).
 */
void mbe_renderToneDstarFromf(float* aout_buf, mbe_parms* mp, int index, int from_volume, int volume, int ramp);

/** @brief mbe_startToneAmbe2450() for a D-STAR tone index. */
void mbe_startToneDstar(mbe_parms* mp, int index);

/**
 * Mute-noise amplitude of TIA-102.BABA 7.8 (IMBE) and TIA-102.BABA-1 5.7
 * (AMBE 3600x2450): uniform in [-5, 5] on the synthesized-speech scale s(n).
 * P25's IMBE path plays it; the AMBE 3600x2450 path plays DVSI's level
 * (MBE_AMBE_MUTE_NOISE_AMPLITUDE, ambe_common.h).
 */
#define MBE_SPEC_MUTE_NOISE_AMPLITUDE 5.0f

/**
 * Gain from the synthesized-speech scale s(n) to this library's float samples,
 * which mbe_floattoshort() multiplies by 7. s(n) is on the 16-bit scale of the
 * TIA-102 reference, so 1/7 makes the short output match DVSI's
 * AMBE-3000 decoders on identical bits (within 0.5 dB). Without it, decoded
 * speech ran 17 dB hotter than the original and the soft clip limited 5-7% of
 * the samples of normally recorded speech. Speech and the 7.8 mute noise, both
 * on s(n), take this gain; comfort noise and tones already have absolute
 * output levels.
 */
#define MBE_SPEECH_OUTPUT_GAIN        (1.0f / 7.0f)

/*
 * Render gains of voiced and unvoiced harmonics (+0.3 dB and -0.4 dB), shared
 * by every codec and applied with the low-band roll-off, so the model and
 * adaptive smoothing see the decoded amplitudes. On identical bits DVSI's
 * AMBE-3000 decoders play voiced frames 0.32-0.38 dB louder and unvoiced frames
 * 0.32-0.40 dB quieter than this synthesis otherwise does, in all three modes.
 */
#define MBE_RENDER_VOICED_GAIN        1.03514217f /* +0.3 dB */
/* The voiced gain on the codec paths (the subframe renderings), whose pitch
 * pulses in 10 ms subframes (mbe_subframe.h) carry a little less energy
 * between two levels than the rendering before them: at +0.55 dB DVSI's
 * decoders play fully voiced frames within 0.02 dB of ours on identical bits
 * (frame medians on development vectors, all three modes); at +0.45 dB they
 * played them 0.09-0.12 dB louder. */
#define MBE_RENDER_VOICED_GAIN_DVSI   1.06536876f /* +0.55 dB */
#define MBE_RENDER_UNVOICED_GAIN      0.95499259f /* -0.4 dB */

/**
 * @brief Fill 160 float samples with uniform noise in [-amplitude, +amplitude].
 *
 * Draws from the same thread-local Java Random-compatible generator as
 * mbe_synthesizeComfortNoisef(), which delegates here with the library comfort gain.
 * The float domain is the synthesized-speech scale s(n) of TIA-102.BABA
 * (int16 output is 7x this), so amplitude 5 is the spec's [-5, 5] mute noise.
 *
 * @param aout_buf Output buffer of 160 float samples.
 * @param amplitude Peak amplitude in the library float domain.
 */
void mbe_synthesizeUniformNoisef(float* aout_buf, float amplitude);

/**
 * @brief Limit float output samples to the library's headroom (95% of full
 *        scale), as the synthesis core limits a frame; NaN passes through.
 *
 * For codec paths that add two signals the core rendered separately (a tone's
 * tail and speech, a fading frame and mute noise): the sum is limited.
 *
 * @param samples Samples, limited in place.
 * @param count   Number of samples.
 */
void mbe_clipFloatBuffer(float* samples, int count);

/**
 * @brief Seed the comfort-noise RNG used by mbe_synthesizeComfortNoisef().
 *
 * Uses Java Random-compatible 48-bit state initialization.
 *
 * @param seed Thread-local seed value (0 maps to a fixed non-zero default).
 */
void mbe_seedComfortNoiseRng(uint32_t seed);

#endif /* MBEINT_MBE_ADAPTIVE_H */
