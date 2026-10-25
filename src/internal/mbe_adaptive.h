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

/**
 * @brief Synthesize speech using a captured pre-enhancement RM0 value.
 * @param aout_buf Output buffer of 160 float samples.
 * @param cur_mp Current frame parameters.
 * @param prev_mp Previous enhanced frame parameters.
 * @param rm0 Sum of squared pre-enhancement amplitudes.
 * @param mute_noise Noise output if the frame is muted.
 */
void mbe_synthesizeSpeechWithPreEnhRm0f(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, float rm0,
                                        enum mbe_mute_noise mute_noise);

/**
 * @brief Synthesize a frame-repeat model as-is, without adaptive smoothing.
 *
 * TIA-102.BABA-1 5.6 synthesizes the repeated model parameters unchanged;
 * smoothing them again could change their voicing and amplitudes.
 *
 * @param aout_buf Output buffer of 160 float samples.
 * @param cur_mp Repeated (already enhanced) frame parameters.
 * @param prev_mp Previous enhanced frame parameters.
 */
void mbe_synthesizeRepeatedSpeechf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp);

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
