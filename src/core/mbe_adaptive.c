// SPDX-License-Identifier: GPL-2.0-or-later
/* Independently implemented from TIA-102.BABA equations 111--116 and the policies in legal/implementation-provenance.md. */
#include "mbe_adaptive.h"
#include "mbe_compiler.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#define MBE_COMFORT_MULTIPLIER UINT64_C(0x5deece66d)
#define MBE_COMFORT_MASK       UINT64_C(0xffffffffffff)
#define MBE_ZERO_SEED          UINT32_C(0x6d25357b)

static MBE_THREAD_LOCAL uint64_t comfort_state = UINT64_C(0x12345678) ^ MBE_COMFORT_MULTIPLIER;

static int
valid_harmonics(const mbe_parms* mp) {
    return mp != NULL && mp->L >= 1 && mp->L <= 56;
}

int
mbe_requiresAdaptiveSmoothing(const mbe_parms* mp) {
    return mp != NULL && (mp->errorRate > MBE_ERROR_THRESHOLD_ENTRY || mp->errorCountTotal > 4);
}

int
mbe_requiresMuting(const mbe_parms* mp) {
    return mp != NULL && mp->errorRate > mp->mutingThreshold;
}

int
mbe_isMaxFrameRepeat(const mbe_parms* mp) {
    return mp != NULL && mp->repeatCount >= MBE_MAX_FRAME_REPEATS;
}

static float
voicing_threshold(const mbe_parms* mp) {
    if (mp->errorRate <= MBE_ERROR_THRESHOLD_LOW && mp->errorCountTotal <= 4) {
        return FLT_MAX;
    }
    /* x^(3/8) = (sqrt(sqrt(sqrt(x))))^3, equation 112. */
    const float root = sqrtf(sqrtf(sqrtf(mp->localEnergy)));
    const float energy_power = (root * root) * root;
    if (mp->errorRate <= MBE_ERROR_THRESHOLD_ENTRY && mp->errorCount4 == 0) {
        return (MBE_ADAPTIVE_GAIN * energy_power) / expf(MBE_ADAPTIVE_EXPONENT * mp->errorRate);
    }
    return MBE_ADAPTIVE_ALT * energy_power;
}

static int
amplitude_threshold(const mbe_parms* cur, int prior) {
    if (cur->errorRate <= MBE_ERROR_THRESHOLD_LOW && cur->errorCountTotal <= 6) {
        return MBE_DEFAULT_AMPLITUDE_THRESHOLD;
    }
    if (prior <= 0) {
        prior = MBE_DEFAULT_AMPLITUDE_THRESHOLD;
    }
    const int64_t updated =
        (int64_t)prior + MBE_AMPLITUDE_BASE - (int64_t)MBE_AMPLITUDE_PENALTY_PER_ERROR * cur->errorCountTotal;
    if (updated > INT_MAX) {
        return INT_MAX;
    }
    if (updated < INT_MIN) {
        return INT_MIN;
    }
    return (int)updated;
}

void
mbe_applyAdaptiveSmoothingWithRm0(mbe_parms* cur_mp, const mbe_parms* prev_mp, float rm0) {
    if (!valid_harmonics(cur_mp) || !valid_harmonics(prev_mp)) {
        return;
    }
    float prior_energy = prev_mp->localEnergy;
    const int prior_threshold = prev_mp->amplitudeThreshold;
    if (prior_energy < MBE_MIN_LOCAL_ENERGY) {
        prior_energy = MBE_DEFAULT_LOCAL_ENERGY;
    }
    float energy = MBE_ENERGY_SMOOTH_ALPHA * prior_energy + MBE_ENERGY_SMOOTH_BETA * rm0;
    if (energy < MBE_MIN_LOCAL_ENERGY) {
        energy = MBE_MIN_LOCAL_ENERGY;
    }
    cur_mp->localEnergy = energy;
    const float threshold = voicing_threshold(cur_mp);
    /* In double: a float sum of loud amplitudes can overflow, and the scale
     * below would then silence every harmonic. */
    double amplitude_sum = 0.0;
    for (int l = 1; l <= cur_mp->L; ++l) {
        if (cur_mp->Ml[l] > threshold) {
            cur_mp->Vl[l] = 1;
        }
        amplitude_sum += (double)cur_mp->Ml[l];
    }
    cur_mp->amplitudeThreshold = amplitude_threshold(cur_mp, prior_threshold);
    float scale = 1.0f;
    if (amplitude_sum > 0.0 && (double)cur_mp->amplitudeThreshold <= amplitude_sum) {
        scale = (float)((double)cur_mp->amplitudeThreshold / amplitude_sum);
    }
    if (scale < 0.0f) {
        scale = 0.0f;
    }
    for (int l = 1; l <= cur_mp->L; ++l) {
        cur_mp->Ml[l] *= scale;
    }
}

void
mbe_applyAdaptiveSmoothing(mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    if (!valid_harmonics(cur_mp) || !valid_harmonics(prev_mp)) {
        return;
    }
    float rm0 = 0.0f;
    for (int l = 1; l <= cur_mp->L; ++l) {
        rm0 += cur_mp->Ml[l] * cur_mp->Ml[l];
    }
    mbe_applyAdaptiveSmoothingWithRm0(cur_mp, prev_mp, rm0);
}

void
mbe_seedComfortNoiseRng(uint32_t seed) {
    if (seed == 0) {
        seed = MBE_ZERO_SEED;
    }
    comfort_state = ((uint64_t)seed ^ MBE_COMFORT_MULTIPLIER) & MBE_COMFORT_MASK;
}

void
mbe_synthesizeUniformNoisef(float* aout_buf, float amplitude) {
    if (aout_buf == NULL) {
        return;
    }
    for (int i = 0; i < 160; ++i) {
        comfort_state = (comfort_state * MBE_COMFORT_MULTIPLIER + 11u) & MBE_COMFORT_MASK;
        const float unit = (float)(comfort_state >> 24) * (1.0f / 16777216.0f);
        aout_buf[i] = ((unit * 2.0f) - 1.0f) * amplitude;
    }
}

void
mbe_synthesizeComfortNoisef(float* aout_buf) {
    mbe_synthesizeUniformNoisef(aout_buf, (0.003f * 32767.0f) / 7.0f);
}

void
mbe_synthesizeComfortNoise(short* aout_buf) {
    if (aout_buf == NULL) {
        return;
    }
    float samples[160];
    mbe_synthesizeComfortNoisef(samples);
    mbe_floattoshort(samples, aout_buf);
}
