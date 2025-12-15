// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 *
 * FFT-based unvoiced synthesis implementation.
 * Implements JMBE Algorithms #117-126 for high-quality unvoiced audio.
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "mbe_compiler.h"
#include "mbe_unvoiced_fft.h"

/* LCG integer constants matching the float #defines in the header */
#define MBE_LCG_A_INT 171u
#define MBE_LCG_B_INT 11213u
#define MBE_LCG_M_INT 53125u

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief 211-element synthesis window (indices -105 to +105).
 *
 * Trapezoidal window with linear ramps at edges and flat region in center.
 * Matches JMBE specification for WOLA synthesis.
 */
static const float Ws_synthesis[211] = {
    /* Indices -105 to -56: linear ramp up from 0.0 to ~0.98 */
    0.000f, 0.020f, 0.040f, 0.060f, 0.080f, 0.100f, 0.120f, 0.140f, 0.160f, 0.180f, 0.200f, 0.220f, 0.240f, 0.260f,
    0.280f, 0.300f, 0.320f, 0.340f, 0.360f, 0.380f, 0.400f, 0.420f, 0.440f, 0.460f, 0.480f, 0.500f, 0.520f, 0.540f,
    0.560f, 0.580f, 0.600f, 0.620f, 0.640f, 0.660f, 0.680f, 0.700f, 0.720f, 0.740f, 0.760f, 0.780f, 0.800f, 0.820f,
    0.840f, 0.860f, 0.880f, 0.900f, 0.920f, 0.940f, 0.960f, 0.980f,
    /* Indices -55 to +55: flat region at 1.0 (111 values) */
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f, 1.000f,
    /* Indices +56 to +105: linear ramp down from ~0.98 to 0.0 */
    0.980f, 0.960f, 0.940f, 0.920f, 0.900f, 0.880f, 0.860f, 0.840f, 0.820f, 0.800f, 0.780f, 0.760f, 0.740f, 0.720f,
    0.700f, 0.680f, 0.660f, 0.640f, 0.620f, 0.600f, 0.580f, 0.560f, 0.540f, 0.520f, 0.500f, 0.480f, 0.460f, 0.440f,
    0.420f, 0.400f, 0.380f, 0.360f, 0.340f, 0.320f, 0.300f, 0.280f, 0.260f, 0.240f, 0.220f, 0.200f, 0.180f, 0.160f,
    0.140f, 0.120f, 0.100f, 0.080f, 0.060f, 0.040f, 0.020f, 0.000f};

/**
 * @brief Fast inline window lookup for hot paths.
 *
 * Directly indexes the window table without function call overhead.
 * Caller must ensure n is in [-105, +105] for valid results.
 */
static inline float
mbe_synthesisWindow_fast(int n) {
    return Ws_synthesis[n + 105];
}

/**
 * @brief FFT plan structure wrapping KISS FFT configs and scratch buffers.
 *
 * The scratch buffers are allocated once per thread and reused across frames
 * to reduce stack traffic and improve cache locality.
 */
struct mbe_fft_plan {
    kiss_fftr_cfg fwd; /**< Forward FFT config */
    kiss_fftr_cfg inv; /**< Inverse FFT config */

    /* Scratch buffers for unvoiced synthesis (reused across frames) */
    float Uw[MBE_FFT_SIZE];                    /**< Windowed noise buffer */
    kiss_fft_cpx Uw_fft[MBE_FFT_SIZE / 2 + 1]; /**< FFT output bins */
    float dftBinScalor[MBE_FFT_SIZE / 2 + 1];  /**< Per-bin scaling factors */
    float Uw_out[MBE_FFT_SIZE];                /**< IFFT output buffer */
    int a_min[57];                             /**< Band lower bin edges */
    int b_max[57];                             /**< Band upper bin edges */
};

mbe_fft_plan*
mbe_fft_plan_alloc(void) {
    mbe_fft_plan* plan = (mbe_fft_plan*)malloc(sizeof(mbe_fft_plan));
    if (!plan) {
        return NULL;
    }

    plan->fwd = kiss_fftr_alloc(MBE_FFT_SIZE, 0, NULL, NULL);
    plan->inv = kiss_fftr_alloc(MBE_FFT_SIZE, 1, NULL, NULL);

    if (!plan->fwd || !plan->inv) {
        if (plan->fwd) {
            kiss_fft_free(plan->fwd);
        }
        if (plan->inv) {
            kiss_fft_free(plan->inv);
        }
        free(plan);
        return NULL;
    }

    return plan;
}

void
mbe_fft_plan_free(mbe_fft_plan* plan) {
    if (plan) {
        if (plan->fwd) {
            kiss_fft_free(plan->fwd);
        }
        if (plan->inv) {
            kiss_fft_free(plan->inv);
        }
        free(plan);
    }
}

float
mbe_synthesisWindow(int n) {
    if (n < -105 || n > 105) {
        return 0.0f;
    }
    return Ws_synthesis[n + 105];
}

void
mbe_generate_noise_lcg(float* restrict buffer, int count, float* restrict seed) {
    if (MBE_UNLIKELY(!buffer || !seed)) {
        return;
    }

    /* Use integer arithmetic for the LCG to avoid expensive fmodf() calls.
     * The state is always in [0, 53124] which is exactly representable in float.
     * This optimization replaces per-sample fmodf() with integer modulo. */
    uint32_t state = (uint32_t)(*seed) % 53125u;
    for (int i = 0; i < count; i++) {
        /* Write current state to buffer BEFORE updating (preserves JMBE sequence) */
        buffer[i] = (float)state;
        state = (171u * state + 11213u) % 53125u;
    }
    *seed = (float)state;
}

void
mbe_generate_noise_with_overlap(float* restrict buffer, float* restrict seed, float* restrict overlap) {
    if (MBE_UNLIKELY(!buffer || !seed || !overlap)) {
        return;
    }

    /* Copy overlap from previous frame (first 96 samples) */
    memcpy(buffer, overlap, MBE_NOISE_OVERLAP * sizeof(float));

    /* Generate 160 new samples (256 - 96 = 160) */
    mbe_generate_noise_lcg(buffer + MBE_NOISE_OVERLAP, MBE_FFT_SIZE - MBE_NOISE_OVERLAP, seed);

    /* Save new overlap for next frame (last 96 samples) */
    memcpy(overlap, buffer + (MBE_FFT_SIZE - MBE_NOISE_OVERLAP), MBE_NOISE_OVERLAP * sizeof(float));
}

void
mbe_wola_combine(float* restrict output, const float* restrict prevUw, const float* restrict currUw, int N) {
    if (MBE_UNLIKELY(!output || !prevUw || !currUw)) {
        return;
    }

    for (int n = 0; n < N; n++) {
        /* Window positions relative to frame center
         * Indices may be outside [-105, +105], use bounds-checked version */
        float w_prev = mbe_synthesisWindow(n);     /* w(n) for previous frame */
        float w_curr = mbe_synthesisWindow(n - N); /* w(n-160) for current frame */

        /* Extract samples from buffers at correct positions
         * prevUw is centered at sample 128, currUw is also centered at 128
         * For output sample n:
         *   - prev contribution: prevUw[n + 128] if in range
         *   - curr contribution: currUw[n - 32] if in range (shifted by N-128=32)
         */
        float prev_sample = 0.0f;
        float curr_sample = 0.0f;

        int prev_idx = n + 128;
        if (MBE_LIKELY(prev_idx >= 0 && prev_idx < MBE_FFT_SIZE)) {
            prev_sample = prevUw[prev_idx];
        }

        int curr_idx = n + 128 - N; /* n - 32 for N=160 */
        if (MBE_LIKELY(curr_idx >= 0 && curr_idx < MBE_FFT_SIZE)) {
            curr_sample = currUw[curr_idx];
        }

        /* Normalized weighted overlap-add */
        float w_prev_sq = w_prev * w_prev;
        float w_curr_sq = w_curr * w_curr;
        float denom = w_prev_sq + w_curr_sq;

        if (MBE_LIKELY(denom > 1e-10f)) {
            output[n] += ((w_prev * prev_sample) + (w_curr * curr_sample)) / denom;
        }
    }
}

void
mbe_synthesizeUnvoicedFFTWithNoise(float* restrict output, mbe_parms* restrict cur_mp, mbe_parms* restrict prev_mp,
                                   mbe_fft_plan* restrict plan, const float* restrict noise_buffer) {
    if (MBE_UNLIKELY(!output || !cur_mp || !prev_mp || !plan || !noise_buffer)) {
        return;
    }

    /* Use plan's scratch buffers instead of stack-allocated arrays */
    float* Uw = plan->Uw;
    kiss_fft_cpx* Uw_fft = plan->Uw_fft;
    float* dftBinScalor = plan->dftBinScalor;
    float* Uw_out = plan->Uw_out;
    int* a_min = plan->a_min;
    int* b_max = plan->b_max;

    int L = cur_mp->L;
    float w0 = cur_mp->w0;

    /* Initialize scalors to zero (voiced bands stay zeroed) */
    memset(dftBinScalor, 0, (MBE_FFT_SIZE / 2 + 1) * sizeof(float));

    /* Copy pre-generated noise buffer and apply synthesis window (Algorithm #118 prep)
     * Window is centered at sample 128, indices go from -128 to +127
     * Use fast inline lookup - all indices are in valid range [-105, 105] after offset
     */
    for (int i = 0; i < MBE_FFT_SIZE; i++) {
        int win_idx = i - 128;
        /* Indices outside [-105, +105] get zero from the window table edges */
        float w = (win_idx >= -105 && win_idx <= 105) ? mbe_synthesisWindow_fast(win_idx) : 0.0f;
        Uw[i] = noise_buffer[i] * w;
    }

    /* Algorithm #118: 256-point real FFT */
    kiss_fftr(plan->fwd, Uw, Uw_fft);

    /* Algorithms #122-123: Calculate frequency band edges for each harmonic */
    float multiplier = MBE_256_OVER_2PI * w0;

    for (int l = 1; l <= L; l++) {
        a_min[l] = (int)ceilf((l - 0.5f) * multiplier);
        b_max[l] = (int)ceilf((l + 0.5f) * multiplier);
        /* Clamp to valid bin range */
        if (a_min[l] < 0) {
            a_min[l] = 0;
        }
        if (b_max[l] > MBE_FFT_SIZE / 2) {
            b_max[l] = MBE_FFT_SIZE / 2;
        }
    }

    /* Algorithm #120: Calculate band-level scaling for unvoiced bands */
    for (int l = 1; l <= L; l++) {
        if (cur_mp->Vl[l] == 0) { /* Unvoiced band */
            float numerator = 0.0f;
            int bin_count = 0;

            for (int bin = a_min[l]; bin < b_max[l]; bin++) {
                float re = Uw_fft[bin].r;
                float im = Uw_fft[bin].i;
                numerator += (re * re) + (im * im);
                bin_count++;
            }

            if (bin_count > 0 && numerator > 1e-10f) {
                float denominator = (float)bin_count;
                float scalor = MBE_UNVOICED_SCALE_COEFF * cur_mp->Ml[l] / sqrtf(numerator / denominator);

                /* Apply scaling factor to all bins in this band */
                for (int bin = a_min[l]; bin < b_max[l]; bin++) {
                    dftBinScalor[bin] = scalor;
                }
            }
        }
        /* Voiced bands: scalor remains 0, effectively zeroing those bins */
    }

    /* Algorithms #119, #120, #124: Apply scaling to FFT bins
     * Voiced bins get scaled by 0 (zeroed), unvoiced bins get proper scaling
     */
    for (int bin = 0; bin <= MBE_FFT_SIZE / 2; bin++) {
        Uw_fft[bin].r *= dftBinScalor[bin];
        Uw_fft[bin].i *= dftBinScalor[bin];
    }

    /* Algorithm #125: Inverse FFT */
    kiss_fftri(plan->inv, Uw_fft, Uw_out);

    /* Normalize IFFT output (KISS FFT doesn't normalize) */
    float scale = 1.0f / (float)MBE_FFT_SIZE;
    for (int i = 0; i < MBE_FFT_SIZE; i++) {
        Uw_out[i] *= scale;
    }

    /* Algorithm #126: WOLA combine with previous frame */
    mbe_wola_combine(output, prev_mp->previousUw, Uw_out, 160);

    /* Save current output for next frame's WOLA */
    memcpy(cur_mp->previousUw, Uw_out, MBE_FFT_SIZE * sizeof(float));
}

void
mbe_synthesizeUnvoicedFFT(float* restrict output, mbe_parms* restrict cur_mp, mbe_parms* restrict prev_mp,
                          mbe_fft_plan* restrict plan) {
    if (MBE_UNLIKELY(!output || !cur_mp || !prev_mp || !plan)) {
        return;
    }

    /* Algorithm #117: Generate 256 white noise samples with LCG and overlap */
    float noise_buffer[MBE_FFT_SIZE];
    mbe_generate_noise_with_overlap(noise_buffer, &cur_mp->noiseSeed, cur_mp->noiseOverlap);

    /* Use the shared implementation */
    mbe_synthesizeUnvoicedFFTWithNoise(output, cur_mp, prev_mp, plan, noise_buffer);
}
