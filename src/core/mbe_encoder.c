// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Analysis front end shared by the speech encoders.
 */

#include "mbe_encoder.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "mbe_validation.h"

#define MBE_ENCODER_PCM_SCALE 32768.0f /* analysis runs on the 16-bit scale */

int
mbe_encoder_frontend_open(struct mbe_encoder_frontend* fe) {
    fe->fft = mbe_fft_plan_alloc();
    fe->acf = mbe_acf_plan_alloc();
    if (fe->fft == NULL || fe->acf == NULL) {
        mbe_encoder_frontend_close(fe);
        return -1;
    }
    mbe_analysis_init_tables(&fe->tables);
    mbe_encoder_frontend_reset(fe);
    return 0;
}

void
mbe_encoder_frontend_close(struct mbe_encoder_frontend* fe) {
    mbe_fft_plan_free(fe->fft);
    mbe_acf_plan_free(fe->acf);
    fe->fft = NULL;
    fe->acf = NULL;
}

void
mbe_encoder_frontend_reset(struct mbe_encoder_frontend* fe) {
    mbe_analysis_reset(&fe->analysis);
}

/*
 * Every test below reads a float's bit pattern from memory and compares its
 * magnitude against a bound, which rejects NaN and infinity, whose patterns
 * lie above every bound used here. Both halves matter under fast-math. Clang
 * marks a by-value float parameter nofpclass(nan inf), so a test on one may be
 * folded away (an exponent-mask test, bits & 0x7F800000 == 0x7F800000, is
 * folded to "finite"), while a load carries no such assumption. And the
 * optimizer reads an exponent-mask test as a NaN/Inf class test, which a value
 * computed by fast-math arithmetic is assumed to pass; a range test against an
 * ordinary bound is left alone.
 */
#define MBE_ENCODER_LIMIT_BITS     0x49800000u /* 2^20 */
#define MBE_ENCODER_LOG2_HIGH_BITS 0x42F00000u /* 120 */
#define MBE_ENCODER_LOG2_LOW_BITS  0x4A800000u /* 2^22 */
#define MBE_ENCODER_AMPLITUDE_BITS 0x7F000000u /* 2^127 */

static uint32_t
mbe_encoder_bits(const float* x) {
    uint32_t bits;
    memcpy(&bits, x, sizeof(bits));
    return bits;
}

static int
mbe_encoder_magnitude_at_most(const float* x, uint32_t limit_bits) {
    return (mbe_encoder_bits(x) & 0x7FFFFFFFu) <= limit_bits;
}

/* Finite and within +-2^20. */
static int
mbe_encoder_value_valid(const float* x) {
    return mbe_encoder_magnitude_at_most(x, MBE_ENCODER_LIMIT_BITS);
}

/* At most 120, so 2^log2Ml times any codec's unvoiced scale (below 1.2) is far
 * from overflow, and at least -2^22 (a history within +-2^20 gives no less
 * than about -2.1 * 2^20). */
static int
mbe_encoder_log2_valid(const float* x) {
    const uint32_t bits = mbe_encoder_bits(x);
    return (bits & 0x7FFFFFFFu) <= (((bits >> 31) != 0u) ? MBE_ENCODER_LOG2_LOW_BITS : MBE_ENCODER_LOG2_HIGH_BITS);
}

int
mbe_encoder_samples_valid(const float samples[MBE_ENCODER_SAMPLES]) {
    for (int i = 0; i < MBE_ENCODER_SAMPLES; i++) {
        if (!mbe_encoder_value_valid(&samples[i])) {
            return 0;
        }
    }
    return 1;
}

int
mbe_encoder_history_valid(const mbe_parms* prev_mp) {
    for (int l = 0; l <= MBE_MAX_HARMONIC_BANDS; l++) {
        if (!mbe_encoder_value_valid(&prev_mp->log2Ml[l])) {
            return 0;
        }
    }
    return mbe_encoder_value_valid(&prev_mp->gamma);
}

/* log2Ml is finite by construction from a validated history, so it is judged
 * first; Ml[l], which overflows exactly when log2Ml[l] is too large, is read
 * only after its log2Ml[l] passed. That order matters even though the bits
 * come from memory: the optimizer may forward a value a fast-math decoder just
 * stored, assumptions included. gamma is left out: it stays within +-2^20
 * from a validated history, and the IMBE decoder does not write it. */
int
mbe_encoder_model_valid(const mbe_parms* mp) {
    const int L = mbe_clamp_harmonic_count(mp->L);
    int ok = mbe_encoder_value_valid(&mp->w0);
    for (int l = 1; ok && l <= L; l++) {
        ok = mbe_encoder_log2_valid(&mp->log2Ml[l])
             && mbe_encoder_magnitude_at_most(&mp->Ml[l], MBE_ENCODER_AMPLITUDE_BITS);
    }
    return ok;
}

void
mbe_encoder_short_to_float(const short in[MBE_ENCODER_SAMPLES], float out[MBE_ENCODER_SAMPLES]) {
    for (int i = 0; i < MBE_ENCODER_SAMPLES; i++) {
        out[i] = (float)in[i] / 32768.0f;
    }
}

int
mbe_encoder_frontend_analyze(struct mbe_encoder_frontend* fe, const float samples[MBE_ENCODER_SAMPLES],
                             struct mbe_analysis_result* result) {
    float scaled[MBE_ENCODER_SAMPLES];
    for (int i = 0; i < MBE_ENCODER_SAMPLES; i++) {
        scaled[i] = samples[i] * MBE_ENCODER_PCM_SCALE;
    }
    mbe_analysis_push(&fe->analysis, scaled);
    return mbe_analysis_frame(&fe->tables, &fe->analysis, fe->fft, fe->acf, result);
}

/* Energy of the envelope mean + s * (a - mean_a), in the log2 amplitude domain. */
static float
mbe_encoder_envelope_energy(const float* a, int L, float mean_a, float mean, float s) {
    float energy = 0.0f;
    for (int l = 1; l <= L; l++) {
        energy += exp2f(2.0f * (mean + (s * (a[l] - mean_a))));
    }
    return energy;
}

/* The same energy as log2, summed relative to the envelope's largest term (each
 * term at most 1, the sum between 1 and L), so it never overflows. */
static float
mbe_encoder_envelope_log2_energy(const float* a, int L, float mean_a, float mean, float s) {
    float top = mean + (s * (a[1] - mean_a));
    for (int l = 2; l <= L; l++) {
        top = fmaxf(top, mean + (s * (a[l] - mean_a)));
    }
    float sum = 0.0f;
    for (int l = 1; l <= L; l++) {
        sum += exp2f(2.0f * ((mean + (s * (a[l] - mean_a))) - top));
    }
    /* The largest term is 2^0, so the sum is at least 1. */
    return (2.0f * top) + log2f(fmaxf(sum, 1.0f));
}

/* The largest log2 amplitude any envelope of the fit reaches: max a for the
 * target, floor_mean + s * (a - mean_a) with s in [0, 1] for the floor. */
static float
mbe_encoder_envelope_peak(const float* a, int L, float mean_a, float floor_mean) {
    float a_max = a[1];
    for (int l = 2; l <= L; l++) {
        a_max = fmaxf(a_max, a[l]);
    }
    return fmaxf(a_max, floor_mean + fmaxf(a_max - mean_a, 0.0f));
}

/* 56 harmonics at 2^(2 * 60) still sum inside the float range. */
#define MBE_ENCODER_ENERGY_PEAK_MAX 60.0f

static float
mbe_encoder_fit_energy(const float* a, int L, float mean_a, float mean, float s, int log2_domain) {
    return log2_domain ? mbe_encoder_envelope_log2_energy(a, L, mean_a, mean, s)
                       : mbe_encoder_envelope_energy(a, L, mean_a, mean, s);
}

float
mbe_encoder_fit_floor(float* a, int L, float mean_a, float floor_mean) {
    /* An input level or prediction history far beyond speech could overflow
     * the energies, which fast-math code must never do; such a frame compares
     * them as log2. Every other frame takes the linear sums. */
    const int log2_domain = mbe_encoder_envelope_peak(a, L, mean_a, floor_mean) > MBE_ENCODER_ENERGY_PEAK_MAX;
    const float target = mbe_encoder_fit_energy(a, L, mean_a, mean_a, 1.0f, log2_domain);
    float lo = 0.0f;
    float hi = 1.0f;
    if (mbe_encoder_fit_energy(a, L, mean_a, floor_mean, 0.0f, log2_domain) >= target) {
        hi = 0.0f;
    }
    for (int i = 0; i < 30 && hi > 0.0f; i++) {
        float mid = 0.5f * (lo + hi);
        if (mbe_encoder_fit_energy(a, L, mean_a, floor_mean, mid, log2_domain) > target) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    for (int l = 1; l <= L; l++) {
        a[l] = floor_mean + (hi * (a[l] - mean_a));
    }
    return floor_mean;
}
