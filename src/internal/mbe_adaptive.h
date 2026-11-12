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
    MBE_MUTE_NOISE_SPEC,    /**< MBE_SPEC_MUTE_NOISE_AMPLITUDE: the P25 IMBE and AMBE 3600x2450 paths. */
};

/** How a frame's model is rendered; chosen by the codec path. */
enum mbe_voiced_synthesis {
    /** TIA-102.BABA: at steady pitch the amplitudes of harmonics 1..7 move
     * linearly over the frame and the others crossfade with the window Ws;
     * the unvoiced part is one transform per frame (eq 117-126). The public
     * mbe_synthesizeSpeechf() renders this way. */
    MBE_VOICED_SYNTHESIS_SPEC,
    /** Every codec path, as DVSI's AMBE-3000 decoders: two 10 ms subframes
     * per frame, the voiced part a train of pitch pulses and the unvoiced
     * part DVSI's noise (mbe_subframe.h; docs/dvsi-parity.md, "Synthesis in
     * subframes"). */
    MBE_VOICED_SYNTHESIS_DVSI,
};

/* The onset delay mbe_set_voicing_onset_delay() can give a model (samples):
 * its voiced bands then start in the frame's second subframe. */
#define MBE_VOICING_ONSET_LATE 80

/**
 * @brief Mark a decoded model whose voiced harmonics start late at an onset:
 * @p late nonzero gives it MBE_VOICING_ONSET_LATE samples, zero none. The
 * mark lives in the otherwise unused element 0 of PHIl.
 */
void mbe_set_voicing_onset_delay(mbe_parms* mp, int late);

/** @brief The onset delay mbe_set_voicing_onset_delay() gave @p mp (0 if none). */
int mbe_voicing_onset_delay(const mbe_parms* mp);

/* Where a frame's voicing must reach for it to hold at a subframe the frame
 * does not own: DVSI's decoders keep the voicing of the frame before a frame
 * marked MBE_FRAME_MARK_PITCH_ABSENT through that frame's first half, and
 * start a P25 onset at once, only when a voiced harmonic lies at or above this
 * frequency. On DVSI's P25 decodes three voiced harmonics reaching 502.6 Hz
 * (pitch code 56) do not count and three reaching 507.9 Hz (code 55) do. */
#define MBE_VOICING_REACH_HZ 505.0f

/**
 * What a codec's decoder tells the synthesizer about a frame beyond its model
 * (docs/dvsi-parity.md, "Voiced onsets and offsets"). A model carries one mark
 * at most, in the otherwise unused element 0 of PHIl, stored as the mark's
 * value in float.
 */
enum mbe_frame_mark {
    MBE_FRAME_MARK_NONE = 0,
    /* No voiced band, and the pitch enters the mid-frame pitch of a voiced
     * frame after it. */
    MBE_FRAME_MARK_PITCH_CARRIES = 1,
    /* No voiced band, at a pitch code other than those DVSI's encoders keep
     * for such frames: the voicing of the frame before it lasts to mid-frame. */
    MBE_FRAME_MARK_PITCH_ABSENT = 2,
    /* No voiced band, and the previous frame's model, as it stands at that
     * frame's end, stands at this frame's mid-frame too: DVSI's rate-33
     * decoder plays it through the first 10 ms of a silence frame and
     * crossfades to the silence model over the second, the pulse train
     * running on through the frame. */
    MBE_FRAME_MARK_PREVIOUS_HELD = 3,
    /* The voiced bands start MBE_VOICING_ONSET_LATE samples late at an onset. */
    MBE_FRAME_MARK_LATE_ONSET = MBE_VOICING_ONSET_LATE,
};

/** @brief Give @p mp a mark, replacing any it had. */
void mbe_set_frame_mark(mbe_parms* mp, enum mbe_frame_mark mark);

/** @brief The mark of @p mp; MBE_FRAME_MARK_NONE for any other stored value. */
enum mbe_frame_mark mbe_frame_mark(const mbe_parms* mp);

/**
 * @brief Regenerate harmonic phases from a spectral envelope (US 5,701,390).
 * @param Ml  Amplitudes, [1..L] read.
 * @param L   Harmonic count, 1..56.
 * @param phi Receives the phases of harmonics 1..56 (those past L from the
 *            envelope's extrapolation); element 0 is not written.
 */
void mbe_regenerate_phases(const float Ml[57], int L, float phi[57]);

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
 * @brief Render 160 samples of an AMBE 3600x2450 tone, continuing the tone phases in @p mp.
 *
 * mbe_synthesizeTonef() without the frame's bits, for a tone known by its
 * fields. An ID with no tone renders silence.
 *
 * @param aout_buf Output buffer of 160 float samples.
 * @param mp Parameter set carrying the tone phases (updated).
 * @param id Tone ID (TIA-102.BABA-1 Table 9).
 * @param ad Tone amplitude field, 0..127.
 */
void mbe_renderToneAmbe2450f(float* aout_buf, mbe_parms* mp, int id, int ad);

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

/**
 * Mute-noise amplitude of TIA-102.BABA 7.8 (IMBE) and TIA-102.BABA-1 5.7
 * (AMBE 3600x2450): uniform in [-5, 5] on the synthesized-speech scale s(n).
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
/* The voiced gain on the codec paths (MBE_VOICED_SYNTHESIS_DVSI), whose pitch
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
 * @brief Seed the comfort-noise RNG used by mbe_synthesizeComfortNoisef().
 *
 * Uses Java Random-compatible 48-bit state initialization.
 *
 * @param seed Thread-local seed value (0 maps to a fixed non-zero default).
 */
void mbe_seedComfortNoiseRng(uint32_t seed);

#endif /* MBEINT_MBE_ADAPTIVE_H */
