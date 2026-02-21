// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 *
 * Adaptive smoothing implementation.
 * Implements JMBE Algorithms #111-116 for error-based audio quality improvement.
 */

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "mbe_adaptive.h"
#include "mbe_compiler.h"
#include "mbe_math.h"
#include "mbelib-neo/mbelib.h"

/* Thread-local storage for comfort noise RNG to avoid cross-thread interference.
 * JMBE uses per-synthesizer java.util.Random instances. */
#define MBE_JAVA_RNG_MULT      0x5DEECE66DULL
#define MBE_JAVA_RNG_ADD       0xBULL
#define MBE_JAVA_RNG_MASK      ((1ULL << 48) - 1ULL)
#define MBE_JAVA_RNG_INIT_SEED 0x12345678ULL
static MBE_THREAD_LOCAL uint64_t mbe_comfort_noise_seed48 = 0;
static MBE_THREAD_LOCAL int mbe_comfort_noise_seeded = 0;
/* Thread-local pre-enhancement RM0 handoff from mbe_spectralAmpEnhance(). */
static MBE_THREAD_LOCAL float mbe_pre_enh_rm0 = 0.0f;
static MBE_THREAD_LOCAL int mbe_pre_enh_rm0_valid = 0;
static MBE_THREAD_LOCAL const mbe_parms* mbe_pre_enh_owner = NULL;

void
mbe_seedComfortNoiseRng(uint32_t seed) {
    if (seed == 0u) {
        seed = 0x6d25357bu;
    }
    mbe_comfort_noise_seed48 = (((uint64_t)seed) ^ MBE_JAVA_RNG_MULT) & MBE_JAVA_RNG_MASK;
    mbe_comfort_noise_seeded = 1;
}

void
mbe_setPreEnhRm0(const mbe_parms* owner, float rm0) {
    if (owner && rm0 >= 0.0f) {
        mbe_pre_enh_rm0 = rm0;
        mbe_pre_enh_rm0_valid = 1;
        mbe_pre_enh_owner = owner;
    } else {
        mbe_pre_enh_rm0 = 0.0f;
        mbe_pre_enh_rm0_valid = 0;
        mbe_pre_enh_owner = NULL;
    }
}

/**
 * @brief Java Random-compatible next(bits) generator for comfort noise.
 *
 * Replicates java.util.Random's 48-bit LCG:
 *   seed = (seed * 0x5DEECE66D + 0xB) & ((1<<48)-1)
 *   return seed >>> (48 - bits)
 *
 * @param bits Number of high-order bits requested (1..32).
 * @return Pseudorandom value with the requested number of bits.
 */
static inline uint32_t
mbe_java_random_next_bits(int bits) {
    if (!mbe_comfort_noise_seeded) {
        mbe_comfort_noise_seed48 = (MBE_JAVA_RNG_INIT_SEED ^ MBE_JAVA_RNG_MULT) & MBE_JAVA_RNG_MASK;
        mbe_comfort_noise_seeded = 1;
    }

    mbe_comfort_noise_seed48 = (mbe_comfort_noise_seed48 * MBE_JAVA_RNG_MULT + MBE_JAVA_RNG_ADD) & MBE_JAVA_RNG_MASK;
    return (uint32_t)(mbe_comfort_noise_seed48 >> (48 - bits));
}

/**
 * @brief Check if adaptive smoothing is required based on error rates.
 *
 * Smoothing is required when error rate exceeds 1.25% or total errors exceed 4.
 *
 * @param mp Parameter set to check.
 * @return Non-zero if smoothing should be applied.
 */
int
mbe_requiresAdaptiveSmoothing(const mbe_parms* mp) {
    if (MBE_UNLIKELY(!mp)) {
        return 0;
    }
    return (mp->errorRate > MBE_ERROR_THRESHOLD_ENTRY) || (mp->errorCountTotal > 4);
}

/**
 * @brief Check if frame should be muted due to excessive errors.
 *
 * Uses the codec-specific muting threshold stored in mp->mutingThreshold.
 * IMBE uses 8.75% (0.0875), AMBE uses 9.6% (0.096).
 *
 * @param mp Parameter set to check.
 * @return Non-zero if frame should be muted.
 */
int
mbe_requiresMuting(const mbe_parms* mp) {
    if (MBE_UNLIKELY(!mp)) {
        return 0;
    }
    return mp->errorRate > mp->mutingThreshold;
}

/**
 * @brief Check if max repeat threshold has been exceeded.
 *
 * @param mp Parameter set to check.
 * @return Non-zero if repeatCount >= MBE_MAX_FRAME_REPEATS.
 */
int
mbe_isMaxFrameRepeat(const mbe_parms* mp) {
    if (MBE_UNLIKELY(!mp)) {
        return 0;
    }
    return mp->repeatCount >= MBE_MAX_FRAME_REPEATS;
}

/**
 * @brief Generate comfort noise for muted frames (float version).
 *
 * Generates low-level uniform white noise to fill gaps during frame muting.
 *
 * @param aout_buf Output buffer of 160 float samples.
 */
void
mbe_synthesizeComfortNoisef(float* aout_buf) {
    if (MBE_UNLIKELY(!aout_buf)) {
        return;
    }

    /* JMBE muted-noise model: uniform white noise in [-1, +1] with gain 0.003.
     * Translate to this library's float-domain scale (short path multiplies by 7). */
    const float gain = (0.003f * 32767.0f) / 7.0f;

    for (int i = 0; i < 160; i++) {
        /* JMBE parity: use Java Random-like 24-bit float generation. */
        float u = ((float)mbe_java_random_next_bits(24) / 16777216.0f) * 2.0f - 1.0f;
        aout_buf[i] = u * gain;
    }
}

/**
 * @brief Generate comfort noise for muted frames (16-bit version).
 *
 * @param aout_buf Output buffer of 160 16-bit samples.
 */
void
mbe_synthesizeComfortNoise(short* aout_buf) {
    if (MBE_UNLIKELY(!aout_buf)) {
        return;
    }

    float float_buf[160];
    mbe_synthesizeComfortNoisef(float_buf);

    /* Reuse float->short scaling so noise amplitude matches JMBE's 0.003 model. */
    mbe_floattoshort(float_buf, aout_buf);
}

/**
 * @brief Apply adaptive smoothing to parameters based on error rates.
 *
 * Implements JMBE Algorithms #111-116:
 * - Algorithm #111: Local energy tracking with IIR smoothing
 * - Algorithm #112: Adaptive threshold calculation
 * - Algorithm #113: Apply threshold to voicing decisions
 * - Algorithm #114: Calculate amplitude measure
 * - Algorithm #115: Calculate amplitude threshold
 * - Algorithm #116: Scale enhanced spectral amplitudes
 *
 * @param cur_mp Current frame parameters (modified in-place).
 * @param prev_mp Previous frame parameters (for local energy).
 */
void
mbe_applyAdaptiveSmoothing(mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    if (MBE_UNLIKELY(!cur_mp || !prev_mp)) {
        return;
    }

    float* M = cur_mp->Ml;
    int* V = cur_mp->Vl;
    int L = cur_mp->L;
    float errorRate = cur_mp->errorRate;
    int errorTotal = cur_mp->errorCountTotal;
    int errorCount4 = cur_mp->errorCount4;

    /* Algorithm #111: Calculate local energy with IIR smoothing */
    float RM0 = 0.0f;
    if (mbe_pre_enh_rm0_valid && mbe_pre_enh_owner == cur_mp) {
        /* Use pre-enhancement RM0 provided by mbe_spectralAmpEnhance(). */
        RM0 = mbe_pre_enh_rm0;
        mbe_pre_enh_rm0_valid = 0;
        mbe_pre_enh_owner = NULL;
    } else {
        /* Fallback for direct callers that bypass spectral enhancement. */
        for (int l = 1; l <= L; l++) {
            RM0 += M[l] * M[l];
        }
    }

    float prevEnergy = prev_mp->localEnergy;
    if (prevEnergy < MBE_MIN_LOCAL_ENERGY) {
        prevEnergy = MBE_DEFAULT_LOCAL_ENERGY;
    }

    cur_mp->localEnergy = MBE_ENERGY_SMOOTH_ALPHA * prevEnergy + MBE_ENERGY_SMOOTH_BETA * RM0;
    if (cur_mp->localEnergy < MBE_MIN_LOCAL_ENERGY) {
        cur_mp->localEnergy = MBE_MIN_LOCAL_ENERGY;
    }

    /* Algorithm #112: Calculate adaptive threshold VM */
    float VM;
    if (errorRate <= MBE_ERROR_THRESHOLD_LOW && errorTotal <= 4) {
        VM = FLT_MAX; /* No smoothing at very low error rates */
    } else {
        /* x^(3/8) = (x^(1/8))^3, where x^(1/8) = sqrtf(sqrtf(sqrtf(x)))
         * Faster than powf() and maintains adequate precision for adaptive smoothing */
        float x8 = sqrtf(sqrtf(sqrtf(cur_mp->localEnergy)));
        float energy = x8 * x8 * x8;
        if (errorRate <= MBE_ERROR_THRESHOLD_ENTRY && errorCount4 == 0) {
            /* Formula 1: exponential decay based on error rate */
            VM = (MBE_ADAPTIVE_GAIN * energy) / expf(MBE_ADAPTIVE_EXPONENT * errorRate);
        } else {
            /* Formula 2: simple scaling for higher error conditions */
            VM = MBE_ADAPTIVE_ALT * energy;
        }
    }

    /* Algorithm #113: Apply threshold to voicing decisions */
    for (int l = 1; l <= L; l++) {
        if (M[l] > VM) {
            V[l] = 1; /* Force voiced when amplitude exceeds threshold */
        }
    }

    /* Algorithm #114: Calculate amplitude measure */
    float Am = 0.0f;
    for (int l = 1; l <= L; l++) {
        Am += M[l];
    }

    /* Algorithm #115: Calculate amplitude threshold */
    int Tm;
    int prevThreshold = prev_mp->amplitudeThreshold;
    if (prevThreshold <= 0) {
        prevThreshold = MBE_DEFAULT_AMPLITUDE_THRESHOLD;
    }

    if (errorRate <= MBE_ERROR_THRESHOLD_LOW && errorTotal <= 6) {
        Tm = MBE_DEFAULT_AMPLITUDE_THRESHOLD;
    } else {
        Tm = MBE_AMPLITUDE_BASE - (MBE_AMPLITUDE_PENALTY_PER_ERROR * errorTotal) + prevThreshold;
    }
    cur_mp->amplitudeThreshold = Tm;

    /* Algorithm #116: Scale enhanced spectral amplitudes if exceeded */
    if (Am > (float)Tm && Am > 0.0f) {
        float scale = (float)Tm / Am;
        for (int l = 1; l <= L; l++) {
            M[l] *= scale;
        }
    }
}
