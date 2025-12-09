// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 *
 * Copyright (C) 2010 mbelib Author
 * GPG Key ID: 0xEA5EFE2C (9E7A 5527 9CDC EBF7 BF1B  D772 4F98 E863 EA5E FE2C)
 *
 * Portions were originally under the ISC license; this mbelib-neo
 * distribution is provided under GPL-2.0-or-later. See LICENSE for details.
 */

/**
 * @file
 * @brief Core MBE parameter utilities and speech/tone synthesis.
 *
 * @defgroup mbe_internal Internal Helpers
 * @brief Private helpers for voiced/unvoiced mixing, RNG, and SIMD paths.
 *
 * These functions and utilities are used internally by the library’s synthesis
 * implementation and are not part of the public API. Names and semantics may
 * change between releases. Where appropriate, helpers note the impact of
 * `MBELIB_ENABLE_SIMD` and `MBELIB_STRICT_ORDER` on determinism and ordering.
 *
 * @{ 
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#if defined(MBELIB_ENABLE_SIMD)
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
#include <emmintrin.h>
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#endif
#if defined(__i386__) || defined(_M_IX86)
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif
#endif
#endif

#include "mbe_math.h"
#include "mbelib-neo/mbelib.h"
#include "mbelib_const.h"

/* Thread-local PRNG state and helpers (xorshift32) */
#if defined(_MSC_VER)
#define MBE_THREAD_LOCAL __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
#define MBE_THREAD_LOCAL __thread
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define MBE_THREAD_LOCAL _Thread_local
#else
/* Fallback for non-C11 compilers: most support __thread as an extension */
#define MBE_THREAD_LOCAL __thread
#endif
static MBE_THREAD_LOCAL uint32_t mbe_rng_state = 0x12345678u;

void
mbe_setThreadRngSeed(uint32_t seed) {
    if (seed == 0u) {
        seed = 0x6d25357bu; /* avoid zero state */
    }
    mbe_rng_state = seed;
}

/**
 * @brief Thread-local xorshift32 PRNG step.
 * @return New 32-bit pseudo-random state value (never 0).
 */
static inline uint32_t
mbe_xorshift32(void) {
    uint32_t x = mbe_rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    mbe_rng_state = x ? x : 0x6d25357bu;
    return mbe_rng_state;
}

/**
 * @brief Generate a uniform float in [0, 1).
 * @return Pseudo-random float in [0, 1).
 */
/** @internal @ingroup mbe_internal */
static float
mbe_rand(void) {
    /* Map upper 24 bits to [0,1) */
    return ((mbe_xorshift32() >> 8) * (1.0f / 16777216.0f));
}

/**
 * @brief Generate a pseudo-random phase in [-pi, +pi].
 * @return Random float in [-pi, +pi].
 */
/** @internal @ingroup mbe_internal */
static float
mbe_rand_phase(void) {
    return (mbe_rand() * (((float)M_PI) * 2.0F)) - ((float)M_PI);
}

/**
 * @brief Accumulate current unvoiced cosine components and advance states.
 *
 * Sums the current cosine terms over a bank of unvoiced oscillators for one
 * sample, adds scaled noise per oscillator when requested, and advances every
 * oscillator state by one step using a sine/cosine recurrence.
 *
 * @param c           In/out cosine register array (size >= count).
 * @param s           In/out sine register array (size >= count).
 * @param cd          Per-oscillator cos(Δθ) increments.
 * @param sd          Per-oscillator sin(Δθ) increments.
 * @param count       Number of oscillators.
 * @param noise_scale Additive random noise scale, 0 to disable.
 * @return Sum of cosine terms for the current sample.
 */
/** @internal @ingroup mbe_internal */
#if defined(MBELIB_STRICT_ORDER)
static inline float
mbe_unvoiced_mix_accum_update(float* restrict c, float* restrict s, const float* restrict cd, const float* restrict sd,
                              int count, float noise_scale) {
    float sum = 0.0f;
    int i = 0;
#if defined(MBELIB_ENABLE_SIMD)
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
    for (; i + 3 < count; i += 4) {
        __m128 vc = _mm_loadu_ps(c + i);
        __m128 vs = _mm_loadu_ps(s + i);
        __m128 vcd = _mm_loadu_ps(cd + i);
        __m128 vsd = _mm_loadu_ps(sd + i);
        /* accumulate current cosines before update */
        float buf[4];
        _mm_storeu_ps(buf, vc);
        sum += buf[0] + buf[1] + buf[2] + buf[3];
        /* advance oscillators */
        __m128 cpn = _mm_sub_ps(_mm_mul_ps(vc, vcd), _mm_mul_ps(vs, vsd));
        __m128 spn = _mm_add_ps(_mm_mul_ps(vs, vcd), _mm_mul_ps(vc, vsd));
        _mm_storeu_ps(c + i, cpn);
        _mm_storeu_ps(s + i, spn);
        if (noise_scale != 0.0f) {
            /* preserve RNG sequence: four scalar draws */
            sum += noise_scale * mbe_rand();
            sum += noise_scale * mbe_rand();
            sum += noise_scale * mbe_rand();
            sum += noise_scale * mbe_rand();
        }
    }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    for (; i + 3 < count; i += 4) {
        float32x4_t vc = vld1q_f32(c + i);
        float32x4_t vs = vld1q_f32(s + i);
        float32x4_t vcd = vld1q_f32(cd + i);
        float32x4_t vsd = vld1q_f32(sd + i);
        float buf[4];
        vst1q_f32(buf, vc);
        sum += buf[0] + buf[1] + buf[2] + buf[3];
        float32x4_t cpn = vsubq_f32(vmulq_f32(vc, vcd), vmulq_f32(vs, vsd));
        float32x4_t spn = vaddq_f32(vmulq_f32(vs, vcd), vmulq_f32(vc, vsd));
        vst1q_f32(c + i, cpn);
        vst1q_f32(s + i, spn);
        if (noise_scale != 0.0f) {
            sum += noise_scale * mbe_rand();
            sum += noise_scale * mbe_rand();
            sum += noise_scale * mbe_rand();
            sum += noise_scale * mbe_rand();
        }
    }
#endif
#endif /* MBELIB_ENABLE_SIMD */
    /* scalar tail (or whole loop if SIMD not enabled) */
    for (; i < count; ++i) {
        sum += c[i];
        if (noise_scale != 0.0f) {
            sum += noise_scale * mbe_rand();
        }
        float cpn = (c[i] * cd[i]) - (s[i] * sd[i]);
        float spn = (s[i] * cd[i]) + (c[i] * sd[i]);
        c[i] = cpn;
        s[i] = spn;
    }
    return sum;
}
#endif /* MBELIB_STRICT_ORDER */

/**
 * @brief Compute four successive unvoiced sums and advance oscillator states.
 *
 * Produces four sums (one per successive sample) of the unvoiced cosine
 * oscillators and advances their states by four steps in total.
 *
 * When MBELIB_STRICT_ORDER is defined, preserves legacy sample-major noise
 * and accumulation ordering; otherwise, uses a faster oscillator-major path
 * (vectorized when available), with noise added after the math.
 *
 * @param c           In/out cosine registers (size >= count).
 * @param s           In/out sine registers   (size >= count).
 * @param cd          Per-oscillator cos(Δθ) increments.
 * @param sd          Per-oscillator sin(Δθ) increments.
 * @param count       Number of oscillators.
 * @param noise_scale Additive random noise scale, 0 to disable.
 * @param sums        Output four sums for the four samples (sums[0..3]).
 */
/** @internal @ingroup mbe_internal */
static inline void
mbe_unvoiced_mix_block4(float* restrict c, float* restrict s, const float* restrict cd, const float* restrict sd,
                        int count, float noise_scale, float sums[4]) {
#if defined(MBELIB_STRICT_ORDER)
    /* Preserve legacy sample-major RNG/accumulation order: 4 passes */
    sums[0] = 0.0f;
    sums[1] = 0.0f;
    sums[2] = 0.0f;
    sums[3] = 0.0f;
    for (int t = 0; t < 4; ++t) {
        sums[t] += mbe_unvoiced_mix_accum_update(c, s, cd, sd, count, noise_scale);
    }
#else
    /* Oscillator-major: faster by processing 4 oscillators per SIMD vector */
    float sum0 = 0.0f, sum1 = 0.0f, sum2 = 0.0f, sum3 = 0.0f;
    int i = 0;
#if defined(MBELIB_ENABLE_SIMD)
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
    __m128 vS0 = _mm_setzero_ps();
    __m128 vS1 = _mm_setzero_ps();
    __m128 vS2 = _mm_setzero_ps();
    __m128 vS3 = _mm_setzero_ps();
    for (; i + 3 < count; i += 4) {
        __m128 vc = _mm_loadu_ps(c + i);
        __m128 vs = _mm_loadu_ps(s + i);
        const __m128 vcd = _mm_loadu_ps(cd + i);
        const __m128 vsd = _mm_loadu_ps(sd + i);
        /* sample 0 */
        vS0 = _mm_add_ps(vS0, vc);
        /* step -> sample 1 */
        __m128 c1 = _mm_sub_ps(_mm_mul_ps(vc, vcd), _mm_mul_ps(vs, vsd));
        __m128 s1 = _mm_add_ps(_mm_mul_ps(vs, vcd), _mm_mul_ps(vc, vsd));
        vS1 = _mm_add_ps(vS1, c1);
        /* step -> sample 2 */
        __m128 c2 = _mm_sub_ps(_mm_mul_ps(c1, vcd), _mm_mul_ps(s1, vsd));
        __m128 s2 = _mm_add_ps(_mm_mul_ps(s1, vcd), _mm_mul_ps(c1, vsd));
        vS2 = _mm_add_ps(vS2, c2);
        /* step -> sample 3 */
        __m128 c3 = _mm_sub_ps(_mm_mul_ps(c2, vcd), _mm_mul_ps(s2, vsd));
        __m128 s3 = _mm_add_ps(_mm_mul_ps(s2, vcd), _mm_mul_ps(c2, vsd));
        vS3 = _mm_add_ps(vS3, c3);
        /* step -> advance state to sample 4 and store */
        __m128 c4 = _mm_sub_ps(_mm_mul_ps(c3, vcd), _mm_mul_ps(s3, vsd));
        __m128 s4 = _mm_add_ps(_mm_mul_ps(s3, vcd), _mm_mul_ps(c3, vsd));
        _mm_storeu_ps(c + i, c4);
        _mm_storeu_ps(s + i, s4);
    }
    /* horizontal add */
    float tmp[4];
    _mm_storeu_ps(tmp, vS0);
    sum0 += tmp[0] + tmp[1] + tmp[2] + tmp[3];
    _mm_storeu_ps(tmp, vS1);
    sum1 += tmp[0] + tmp[1] + tmp[2] + tmp[3];
    _mm_storeu_ps(tmp, vS2);
    sum2 += tmp[0] + tmp[1] + tmp[2] + tmp[3];
    _mm_storeu_ps(tmp, vS3);
    sum3 += tmp[0] + tmp[1] + tmp[2] + tmp[3];
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    float32x4_t vS0 = vdupq_n_f32(0.0f);
    float32x4_t vS1 = vdupq_n_f32(0.0f);
    float32x4_t vS2 = vdupq_n_f32(0.0f);
    float32x4_t vS3 = vdupq_n_f32(0.0f);
    for (; i + 3 < count; i += 4) {
        float32x4_t vc = vld1q_f32(c + i);
        float32x4_t vs = vld1q_f32(s + i);
        const float32x4_t vcd = vld1q_f32(cd + i);
        const float32x4_t vsd = vld1q_f32(sd + i);
        /* sample 0 */
        vS0 = vaddq_f32(vS0, vc);
        /* step -> sample 1 */
        float32x4_t c1 = vsubq_f32(vmulq_f32(vc, vcd), vmulq_f32(vs, vsd));
        float32x4_t s1 = vaddq_f32(vmulq_f32(vs, vcd), vmulq_f32(vc, vsd));
        vS1 = vaddq_f32(vS1, c1);
        /* step -> sample 2 */
        float32x4_t c2 = vsubq_f32(vmulq_f32(c1, vcd), vmulq_f32(s1, vsd));
        float32x4_t s2 = vaddq_f32(vmulq_f32(s1, vcd), vmulq_f32(c1, vsd));
        vS2 = vaddq_f32(vS2, c2);
        /* step -> sample 3 */
        float32x4_t c3 = vsubq_f32(vmulq_f32(c2, vcd), vmulq_f32(s2, vsd));
        float32x4_t s3 = vaddq_f32(vmulq_f32(s2, vcd), vmulq_f32(c2, vsd));
        vS3 = vaddq_f32(vS3, c3);
        /* step -> advance state to sample 4 and store */
        float32x4_t c4 = vsubq_f32(vmulq_f32(c3, vcd), vmulq_f32(s3, vsd));
        float32x4_t s4 = vaddq_f32(vmulq_f32(s3, vcd), vmulq_f32(c3, vsd));
        vst1q_f32(c + i, c4);
        vst1q_f32(s + i, s4);
    }
    float32x4_t vtmp;
    vtmp = vS0;
    float tbuf0[4];
    vst1q_f32(tbuf0, vtmp);
    sum0 += tbuf0[0] + tbuf0[1] + tbuf0[2] + tbuf0[3];
    vtmp = vS1;
    float tbuf1[4];
    vst1q_f32(tbuf1, vtmp);
    sum1 += tbuf1[0] + tbuf1[1] + tbuf1[2] + tbuf1[3];
    vtmp = vS2;
    float tbuf2[4];
    vst1q_f32(tbuf2, vtmp);
    sum2 += tbuf2[0] + tbuf2[1] + tbuf2[2] + tbuf2[3];
    vtmp = vS3;
    float tbuf3[4];
    vst1q_f32(tbuf3, vtmp);
    sum3 += tbuf3[0] + tbuf3[1] + tbuf3[2] + tbuf3[3];
#endif
#endif /* MBELIB_ENABLE_SIMD */
    /* scalar tail */
    for (; i < count; ++i) {
        float ci = c[i];
        float si = s[i];
        const float cdi = cd[i];
        const float sdi = sd[i];
        sum0 += ci;
        float c1 = (ci * cdi) - (si * sdi);
        float s1 = (si * cdi) + (ci * sdi);
        sum1 += c1;
        float c2 = (c1 * cdi) - (s1 * sdi);
        float s2 = (s1 * cdi) + (c1 * sdi);
        sum2 += c2;
        float c3 = (c2 * cdi) - (s2 * sdi);
        float s3 = (s2 * cdi) + (c2 * sdi);
        sum3 += c3;
        float c4 = (c3 * cdi) - (s3 * sdi);
        float s4 = (s3 * cdi) + (c3 * sdi);
        c[i] = c4;
        s[i] = s4;
    }
    /* add noise contributions after math (different deterministic sequence) */
    if (noise_scale != 0.0f) {
        for (int k = 0; k < count; ++k) {
            sum0 += noise_scale * mbe_rand();
        }
        for (int k = 0; k < count; ++k) {
            sum1 += noise_scale * mbe_rand();
        }
        for (int k = 0; k < count; ++k) {
            sum2 += noise_scale * mbe_rand();
        }
        for (int k = 0; k < count; ++k) {
            sum3 += noise_scale * mbe_rand();
        }
    }
    sums[0] = sum0;
    sums[1] = sum1;
    sums[2] = sum2;
    sums[3] = sum3;
#endif /* MBELIB_STRICT_ORDER */
}

/*
 * Vectorized helper: add voiced oscillator contribution for 4 consecutive
 * samples to the output buffer, and advance oscillator state by 4 steps.
 * This optimizes the per-sample multiply-add with window weights.
 */
/**
 * @brief Add four voiced samples to the output using oscillator recurrence.
 *
 * Generates four consecutive cosine samples from the current oscillator state,
 * multiplies by window `W[0..3]` and amplitude `amp`, and accumulates into
 * `Ss[0..3]`. Advances the oscillator state by four steps.
 *
 * @param Ss Output sample pointer (adds to Ss[0..3]).
 * @param W  Window values for these four samples.
 * @param amp Amplitude multiplier for this band/component.
 * @param c  In/out cosine register.
 * @param s  In/out sine register.
 * @param sd sin(Δθ) per step.
 * @param cd cos(Δθ) per step.
 */
/** @internal @ingroup mbe_internal */
static inline void
mbe_add_voiced_block4(float* restrict Ss, const float* restrict W, float amp, float* restrict c, float* restrict s,
                      float sd, float cd) {
    /* produce c[0..3] from current state */
    float cblk[4];
    float cc = *c, ss = *s;
    for (int k = 0; k < 4; ++k) {
        cblk[k] = cc;
        float cpn = (cc * cd) - (ss * sd);
        float spn = (ss * cd) + (cc * sd);
        cc = cpn;
        ss = spn;
    }
    *c = cc;
    *s = ss;
#if defined(MBELIB_ENABLE_SIMD)
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
    __m128 vC = _mm_loadu_ps(cblk);
    __m128 vW = _mm_loadu_ps(W);
    __m128 vA = _mm_set1_ps(amp);
    __m128 vS = _mm_loadu_ps(Ss);
    vS = _mm_add_ps(vS, _mm_mul_ps(_mm_mul_ps(vC, vW), vA));
    _mm_storeu_ps(Ss, vS);
    return;
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    float32x4_t vC = vld1q_f32(cblk);
    float32x4_t vW = vld1q_f32(W);
    float32x4_t vA = vdupq_n_f32(amp);
    float32x4_t vS = vld1q_f32(Ss);
    vS = vaddq_f32(vS, vmulq_f32(vmulq_f32(vC, vW), vA));
    vst1q_f32(Ss, vS);
    return;
#endif
#endif
    /* scalar fallback */
    Ss[0] += W[0] * amp * cblk[0];
    Ss[1] += W[1] * amp * cblk[1];
    Ss[2] += W[2] * amp * cblk[2];
    Ss[3] += W[3] * amp * cblk[3];
}

/**
 * @brief Write the library version string into the provided buffer.
 * @param str Output buffer receiving a NUL-terminated version string.
 */
void
mbe_printVersion(char* str) {
    if (!str) {
        return;
    }
    /* Ensure we never overrun caller buffer; tests pass a 32B buffer. */
    (void)snprintf(str, 32, "%s", MBELIB_VERSION);
}

/* @} end of mbe_internal */

/* Convenience accessor that avoids buffer management. */
const char*
mbe_versionString(void) {
    return MBELIB_VERSION;
}

/**
 * @brief Copy MBE parameter set from input to output.
 * @param in  Source parameter set.
 * @param out Destination parameter set.
 */
void
mbe_moveMbeParms(mbe_parms* in, mbe_parms* out) {

    int l;
    out->swn = in->swn;
    out->w0 = in->w0;
    out->L = in->L;
    out->K = in->K;
    out->Ml[0] = (float)0;
    out->gamma = in->gamma;
    out->repeat = in->repeat;
    for (l = 0; l <= 56; l++) {
        out->Ml[l] = in->Ml[l];
        out->Vl[l] = in->Vl[l];
        out->log2Ml[l] = in->log2Ml[l];
        out->PHIl[l] = in->PHIl[l];
        out->PSIl[l] = in->PSIl[l];
    }
}

/**
 * @brief Replace current parameters with the last known parameters.
 * @param out Destination parameter set to fill.
 * @param in  Source parameter set from previous frame.
 */
void
mbe_useLastMbeParms(mbe_parms* out, mbe_parms* in) {

    int l;
    out->swn = in->swn;
    out->w0 = in->w0;
    out->L = in->L;
    out->K = in->K;
    out->Ml[0] = (float)0;
    out->gamma = in->gamma;
    out->repeat = in->repeat;
    for (l = 0; l <= 56; l++) {
        out->Ml[l] = in->Ml[l];
        out->Vl[l] = in->Vl[l];
        out->log2Ml[l] = in->log2Ml[l];
        out->PHIl[l] = in->PHIl[l];
        out->PSIl[l] = in->PSIl[l];
    }
}

/**
 * @brief Initialize MBE parameter state for decoding and synthesis.
 * @param cur_mp Output: current parameter state.
 * @param prev_mp Output: previous parameter state (reset to defaults).
 * @param prev_mp_enhanced Output: enhanced previous parameter state.
 */
void
mbe_initMbeParms(mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {

    int l;
    prev_mp->swn = 0;
    prev_mp->w0 = 0.09378;
    prev_mp->L = 30;
    prev_mp->K = 10;
    prev_mp->gamma = (float)0;
    for (l = 0; l <= 56; l++) {
        prev_mp->Ml[l] = (float)0;
        prev_mp->Vl[l] = 0;
        prev_mp->log2Ml[l] = (float)0; // log2 of 1 == 0
        prev_mp->PHIl[l] = (float)0;
        prev_mp->PSIl[l] = (M_PI / (float)2);
    }
    prev_mp->repeat = 0;
    mbe_moveMbeParms(prev_mp, cur_mp);
    mbe_moveMbeParms(prev_mp, prev_mp_enhanced);
}

/**
 * @brief Apply spectral amplitude enhancement to the current parameters.
 * @param cur_mp In/out parameter set to enhance.
 */
void
mbe_spectralAmpEnhance(mbe_parms* cur_mp) {

    float Rm0, Rm1, R2m0, R2m1, Wl[57];
    int l;
    float sum, gamma, M;

    // Precompute cos(w0 * l) table via recurrence to avoid repeated cosf
    float cos_tab[57];
    {
        const float w0 = cur_mp->w0;
        float s_step, c_step;
        mbe_sincosf(w0, &s_step, &c_step);
        float c = 1.0f, s = 0.0f; // angle = 0
        for (l = 1; l <= cur_mp->L; ++l) {
            // advance by w0
            float cn = (c * c_step) - (s * s_step);
            float sn = (s * c_step) + (c * s_step);
            c = cn;
            s = sn;
            cos_tab[l] = c;
        }
    }

    Rm0 = 0.0f;
    Rm1 = 0.0f;
    for (l = 1; l <= cur_mp->L; l++) {
        const float Ml2 = cur_mp->Ml[l] * cur_mp->Ml[l];
        Rm0 += Ml2;
        Rm1 += Ml2 * cos_tab[l];
    }

    R2m0 = (Rm0 * Rm0);
    R2m1 = (Rm1 * Rm1);

    for (l = 1; l <= cur_mp->L; l++) {
        if (cur_mp->Ml[l] != 0.0f) {
            const float cos_w0l = cos_tab[l];
            Wl[l] = sqrtf(cur_mp->Ml[l])
                    * sqrtf(sqrtf(((float)0.96 * (float)M_PI * ((R2m0 + R2m1) - ((float)2 * Rm0 * Rm1 * cos_w0l)))
                                  / (cur_mp->w0 * Rm0 * (R2m0 - R2m1))));

            if ((8 * l) <= cur_mp->L) {
                // no-op
            } else if (Wl[l] > 1.2f) {
                cur_mp->Ml[l] = 1.2f * cur_mp->Ml[l];
            } else if (Wl[l] < 0.5f) {
                cur_mp->Ml[l] = 0.5f * cur_mp->Ml[l];
            } else {
                cur_mp->Ml[l] = Wl[l] * cur_mp->Ml[l];
            }
        }
    }

    // generate scaling factor
    sum = 0.0f;
    for (l = 1; l <= cur_mp->L; l++) {
        M = cur_mp->Ml[l];
        if (M < 0.0f) {
            M = -M;
        }
        sum += (M * M);
    }
    if (sum == 0.0f) {
        gamma = 1.0f;
    } else {
        gamma = sqrtf(Rm0 / sum);
    }

    // apply scaling factor
    for (l = 1; l <= cur_mp->L; l++) {
        cur_mp->Ml[l] = gamma * cur_mp->Ml[l];
    }
}

// Tone synthesis mapping adapted from OP25 (Boatbod)
/**
 * @brief Synthesize a tone frame into 160 float samples at 8 kHz.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49) providing tone indices.
 * @param cur_mp   Current parameter set (tone synthesis state).
 */
void
mbe_synthesizeTonef(float* aout_buf, char* ambe_d, mbe_parms* cur_mp) {
    int i, n;
    float* aout_buf_p;

    int u0, u1, u2, u3;
    u0 = u1 = u2 = u3 = 0;

    for (i = 0; i < 12; i++) {
        u0 = u0 << 1;
        u0 = u0 | (int)ambe_d[i];
    }

    for (i = 12; i < 24; i++) {
        u1 = u1 << 1;
        u1 = u1 | (int)ambe_d[i];
    }

    for (i = 24; i < 35; i++) {
        u2 = u2 << 1;
        u2 = u2 | (int)ambe_d[i];
    }

    for (i = 35; i < 49; i++) {
        u3 = u3 << 1;
        u3 = u3 | (int)ambe_d[i];
    }

    int AD, ID0, ID1, ID2, ID3, ID4;
    AD = ((u0 & 0x3f) << 1) + ((u3 >> 4) & 0x1);
    ID0 = 0;
    ID1 = ((u1 & 0xfff) >> 4);
    ID2 = ((u1 & 0xf) << 4) + ((u2 >> 7) & 0xf);
    ID3 = ((u2 & 0x7f) << 1) + ((u2 >> 13) & 0x1);
    ID4 = ((u3 & 0x1fe0) >> 5);

    // TODO: Cross-validate related ID fields (per OP25). For now, rely on error counts.

    float step1, step2, amplitude;
    float freq1 = 0, freq2 = 0;
    (void)ID0;
    (void)ID2;
    (void)ID3;
    (void)ID4; // parsed for potential validation

#ifdef DISABLE_AMBE_TONES // generate silence if tones disabled
    aout_buf_p = aout_buf;
    for (n = 0; n < 160; n++) {
        *aout_buf_p = (float)0;
        aout_buf_p++;
    }
    return;
#endif

    // Current implementation selects tones solely by ID1
    switch (ID1) {
        // single tones, set frequency
        case 5:
            freq1 = 156.25;
            freq2 = freq1;
            break;
        case 6:
            freq1 = 187.5;
            freq2 = freq1;
            break;
        // DTMF
        case 128:
            freq1 = 1336;
            freq2 = 941;
            break;
        case 129:
            freq1 = 1209;
            freq2 = 697;
            break;
        case 130:
            freq1 = 1336;
            freq2 = 697;
            break;
        case 131:
            freq1 = 1477;
            freq2 = 697;
            break;
        case 132:
            freq1 = 1209;
            freq2 = 770;
            break;
        case 133:
            freq1 = 1336;
            freq2 = 770;
            break;
        case 134:
            freq1 = 1477;
            freq2 = 770;
            break;
        case 135:
            freq1 = 1209;
            freq2 = 852;
            break;
        case 136:
            freq1 = 1336;
            freq2 = 852;
            break;
        case 137:
            freq1 = 1477;
            freq2 = 852;
            break;
        case 138:
            freq1 = 1633;
            freq2 = 697;
            break;
        case 139:
            freq1 = 1633;
            freq2 = 770;
            break;
        case 140:
            freq1 = 1633;
            freq2 = 852;
            break;
        case 141:
            freq1 = 1633;
            freq2 = 941;
            break;
        case 142:
            freq1 = 1209;
            freq2 = 941;
            break;
        case 143:
            freq1 = 1477;
            freq2 = 941;
            break;
        // KNOX
        case 144:
            freq1 = 1162;
            freq2 = 820;
            break;
        case 145:
            freq1 = 1052;
            freq2 = 606;
            break;
        case 146:
            freq1 = 1162;
            freq2 = 606;
            break;
        case 147:
            freq1 = 1279;
            freq2 = 606;
            break;
        case 148:
            freq1 = 1052;
            freq2 = 672;
            break;
        case 149:
            freq1 = 1162;
            freq2 = 672;
            break;
        case 150:
            freq1 = 1279;
            freq2 = 672;
            break;
        case 151:
            freq1 = 1052;
            freq2 = 743;
            break;
        case 152:
            freq1 = 1162;
            freq2 = 743;
            break;
        case 153:
            freq1 = 1279;
            freq2 = 743;
            break;
        case 154:
            freq1 = 1430;
            freq2 = 606;
            break;
        case 155:
            freq1 = 1430;
            freq2 = 672;
            break;
        case 156:
            freq1 = 1430;
            freq2 = 743;
            break;
        case 157:
            freq1 = 1430;
            freq2 = 820;
            break;
        case 158:
            freq1 = 1052;
            freq2 = 820;
            break;
        case 159:
            freq1 = 1279;
            freq2 = 820;
            break;
        // dual tones
        case 160:
            freq1 = 440;
            freq2 = 350;
            break;
        case 161:
            freq1 = 480;
            freq2 = 440;
            break;
        case 162:
            freq1 = 620;
            freq2 = 480;
            break;
        case 163:
            freq1 = 490;
            freq2 = 350;
            break;
        // zero amplitude
        case 255:
            freq1 = 0;
            freq2 = 0;
            break;
        // single tones, calculated frequency
        default:
            if ((ID1 >= 7) && (ID1 <= 122)) {
                freq1 = 31.25 * ID1;
                freq2 = freq1;
            }
    }

    // Zero amplitude or unimplemented tone IDs
    if ((freq1 == 0) && (freq2 == 0)) {
        aout_buf_p = aout_buf;
        for (n = 0; n < 160; n++) {
            *aout_buf_p = (float)0;
            aout_buf_p++;
        }
        return;
    }

    // Debug: uncomment to inspect tone ID and amplitude
    // fprintf(stderr, "TONE ID = %d AD = %d\n", ID1, AD);

    // Synthesize tones
    step1 = 2.0f * (float)M_PI * freq1 / 8000.0f;
    step2 = 2.0f * (float)M_PI * freq2 / 8000.0f;
    amplitude = AD * 75.0f; //
    aout_buf_p = aout_buf;
    for (n = 0; n < 160; n++) {
        *aout_buf_p = amplitude * (sinf((cur_mp->swn) * step1) / 2.0f + sinf((cur_mp->swn) * step2) / 2.0f);
        *aout_buf_p = *aout_buf_p / 6.0f;
        aout_buf_p++;
        cur_mp->swn++;
    }
}

// Simplified D-STAR single-frequency tone synthesis based on existing approximations
/**
 * @brief Synthesize a D-STAR style tone into 160 float samples.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49).
 * @param cur_mp   Current parameter set.
 * @param ID1      Tone index selector.
 */
void
mbe_synthesizeTonefdstar(float* aout_buf, char* ambe_d, mbe_parms* cur_mp, int ID1) {
    int n;
    float* aout_buf_p;

    int AD = 103; // nominal amplitude aligned with other tone cases
    float step1, step2, amplitude;
    float freq1 = 0, freq2 = 0;
    (void)ambe_d;

#ifdef DISABLE_AMBE_TONES // generate silence if tones disabled
    aout_buf_p = aout_buf;
    for (n = 0; n < 160; n++) {
        *aout_buf_p = (float)0;
        aout_buf_p++;
    }
    return;
#endif

    switch (ID1) {
        // single tones, set frequency
        case 5:
            freq1 = 156.25;
            freq2 = freq1;
            break;
        case 6:
            freq1 = 187.5;
            freq2 = freq1;
            break;
        // single tones, calculated frequency
        default:
            if ((ID1 >= 7) && (ID1 <= 122)) {
                freq1 = 31.25 * ID1;
                freq2 = freq1;
            }
    }

    // Debug: uncomment to inspect tone ID and amplitude
    // fprintf(stderr, "TONE ID = %d AD = %d\n", ID1, AD);

    // Synthesize tones
    step1 = 2.0f * (float)M_PI * freq1 / 8000.0f;
    step2 = 2.0f * (float)M_PI * freq2 / 8000.0f;
    amplitude = AD * 75.0f; //
    aout_buf_p = aout_buf;
    for (n = 0; n < 160; n++) {
        *aout_buf_p = amplitude * (sinf((cur_mp->swn) * step1) / 2.0f + sinf((cur_mp->swn) * step2) / 2.0f);
        *aout_buf_p = *aout_buf_p / 6.0f;
        aout_buf_p++;
        cur_mp->swn++;
    }
}

/**
 * @brief Write 160 float samples of silence.
 * @param aout_buf Output buffer of 160 float samples.
 */
void
mbe_synthesizeSilencef(float* aout_buf) {
    memset(aout_buf, 0, 160 * sizeof(*aout_buf));
}

/**
 * @brief Write 160 16-bit samples of silence.
 * @param aout_buf Output buffer of 160 16-bit samples.
 */
void
mbe_synthesizeSilence(short* aout_buf) {
    memset(aout_buf, 0, 160 * sizeof(*aout_buf));
}

/**
 * @brief Synthesize one speech frame into 160 float samples at 8 kHz.
 * @param aout_buf Output buffer of 160 float samples.
 * @param cur_mp   Current parameter set.
 * @param prev_mp  Previous parameter set.
 * @param uvquality Unvoiced synthesis quality (1..64).
 */
void
mbe_synthesizeSpeechf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, int uvquality) {

    int i, l, n, maxl;
    float *Ss, loguvquality;
    int numUv;
    float cw0, pw0, cw0l, pw0l;
    float uvsine, uvrand, uvthreshold, uvthresholdf;
    float uvstep, uvoffset;
    float qfactor;
    float rphase[64], rphase2[64];

    const int N = 160;

    uvthresholdf = (float)2700;
    uvthreshold = ((uvthresholdf * M_PI) / (float)4000);

    // voiced/unvoiced/gain settings
    uvsine = (float)1.3591409 * M_E;
    uvrand = (float)2.0;

    if ((uvquality < 1) || (uvquality > 64)) {
        fprintf(stderr, "\nmbelib: Error - uvquality must be within the range 1 - 64, setting to default value of 3\n");
        uvquality = 3;
    }

    // calculate loguvquality
    if (uvquality == 1) {
        loguvquality = (float)1 / M_E;
    } else {
        loguvquality = logf((float)uvquality) / (float)uvquality;
    }

    // calculate unvoiced step and offset values
    uvstep = (float)1.0 / (float)uvquality;
    qfactor = loguvquality;
    uvoffset = (uvstep * (float)(uvquality - 1)) / (float)2;

    // count number of unvoiced bands
    numUv = 0;
    for (l = 1; l <= cur_mp->L; l++) {
        if (cur_mp->Vl[l] == 0) {
            numUv++;
        }
    }

    cw0 = cur_mp->w0;
    pw0 = prev_mp->w0;

    // init aout_buf
    Ss = aout_buf;
    for (n = 0; n < N; n++) {
        *Ss = (float)0;
        Ss++;
    }

    // eq 128 and 129
    if (cur_mp->L > prev_mp->L) {
        maxl = cur_mp->L;
        for (l = prev_mp->L + 1; l <= maxl; l++) {
            prev_mp->Ml[l] = (float)0;
            prev_mp->Vl[l] = 1;
        }
    } else {
        maxl = prev_mp->L;
        for (l = cur_mp->L + 1; l <= maxl; l++) {
            cur_mp->Ml[l] = (float)0;
            cur_mp->Vl[l] = 1;
        }
    }

    // update phil from eq 139,140
    for (l = 1; l <= 56; l++) {
        cur_mp->PSIl[l] = prev_mp->PSIl[l] + ((pw0 + cw0) * ((float)(l * N) / (float)2));
        if (l <= (cur_mp->L / 4)) {
            cur_mp->PHIl[l] = cur_mp->PSIl[l];
        } else {
            cur_mp->PHIl[l] = cur_mp->PSIl[l] + ((numUv * mbe_rand_phase()) / cur_mp->L);
        }
    }

    for (l = 1; l <= maxl; l++) {
        cw0l = (cw0 * (float)l);
        pw0l = (pw0 * (float)l);
        if ((cur_mp->Vl[l] == 0) && (prev_mp->Vl[l] == 1)) {
            Ss = aout_buf;
            // init random phase
            for (i = 0; i < uvquality; i++) {
                rphase[i] = mbe_rand_phase();
            }
            /* Recurrence for voiced (prev) component */
            const float amp_prev = prev_mp->Ml[l];
            float sd_prev, cd_prev;
            mbe_sincosf(pw0l, &sd_prev, &cd_prev);
            float s_prev, c_prev;
            mbe_sincosf(prev_mp->PHIl[l], &s_prev, &c_prev);
            /* Precompute unvoiced oscillator steps for current (cw0) */
            float sd_u[64] = {0.0f}, cd_u[64] = {0.0f};
            float su[64] = {0.0f}, cu[64] = {0.0f};
            const float base = (float)l - uvoffset;
            for (i = 0; i < uvquality; i++) {
                float inc = cw0 * (base + ((float)i * uvstep));
                mbe_sincosf(inc, &sd_u[i], &cd_u[i]);
                mbe_sincosf(rphase[i], &su[i], &cu[i]);
            }
            const float noise_c = (cw0l > uvthreshold) ? ((cw0l - uvthreshold) * uvrand) : 0.0f;
            for (n = 0; n < N; n += 4) {
                /* add 4-sample voiced block (prev) */
                mbe_add_voiced_block4(Ss, Ws + n + N, amp_prev, &c_prev, &s_prev, sd_prev, cd_prev);
                /* unvoiced 4-sample block */
                float sumu[4];
                mbe_unvoiced_mix_block4(cu, su, cd_u, sd_u, uvquality, noise_c, sumu);
                Ss[0] += (uvsine * Ws[n + 0] * cur_mp->Ml[l] * qfactor) * sumu[0];
                Ss[1] += (uvsine * Ws[n + 1] * cur_mp->Ml[l] * qfactor) * sumu[1];
                Ss[2] += (uvsine * Ws[n + 2] * cur_mp->Ml[l] * qfactor) * sumu[2];
                Ss[3] += (uvsine * Ws[n + 3] * cur_mp->Ml[l] * qfactor) * sumu[3];
                Ss += 4;
            }
        } else if ((cur_mp->Vl[l] == 1) && (prev_mp->Vl[l] == 0)) {
            Ss = aout_buf;
            // init random phase
            for (i = 0; i < uvquality; i++) {
                rphase[i] = mbe_rand_phase();
            }
            /* Recurrence for voiced (current) component */
            const float amp_cur = cur_mp->Ml[l];
            float sd_cur, cd_cur;
            mbe_sincosf(cw0l, &sd_cur, &cd_cur);
            float s_cur, c_cur;
            mbe_sincosf(cur_mp->PHIl[l] - (cw0l * (float)N), &s_cur, &c_cur);
            /* Precompute unvoiced oscillator steps for previous (pw0) */
            float sd_u[64] = {0.0f}, cd_u[64] = {0.0f};
            float su[64] = {0.0f}, cu[64] = {0.0f};
            const float base = (float)l - uvoffset;
            for (i = 0; i < uvquality; i++) {
                float inc = pw0 * (base + ((float)i * uvstep));
                mbe_sincosf(inc, &sd_u[i], &cd_u[i]);
                mbe_sincosf(rphase[i], &su[i], &cu[i]);
            }
            const float noise_p = (pw0l > uvthreshold) ? ((pw0l - uvthreshold) * uvrand) : 0.0f;
            for (n = 0; n < N; n += 4) {
                /* add 4-sample voiced block (cur) */
                mbe_add_voiced_block4(Ss, Ws + n, amp_cur, &c_cur, &s_cur, sd_cur, cd_cur);
                /* unvoiced 4-sample block (prev) */
                float sumu[4];
                mbe_unvoiced_mix_block4(cu, su, cd_u, sd_u, uvquality, noise_p, sumu);
                Ss[0] += (uvsine * Ws[n + 0 + N] * prev_mp->Ml[l] * qfactor) * sumu[0];
                Ss[1] += (uvsine * Ws[n + 1 + N] * prev_mp->Ml[l] * qfactor) * sumu[1];
                Ss[2] += (uvsine * Ws[n + 2 + N] * prev_mp->Ml[l] * qfactor) * sumu[2];
                Ss[3] += (uvsine * Ws[n + 3 + N] * prev_mp->Ml[l] * qfactor) * sumu[3];
                Ss += 4;
            }
        }
        //      else if (((cur_mp->Vl[l] == 1) || (prev_mp->Vl[l] == 1)) && ((l >= 8) || (fabsf (cw0 - pw0) >= ((float) 0.1 * cw0))))
        else if ((cur_mp->Vl[l] == 1) || (prev_mp->Vl[l] == 1)) {
            Ss = aout_buf;
            const float amp_prev = prev_mp->Ml[l];
            const float amp_cur = cur_mp->Ml[l];
            float sd_prev, cd_prev;
            mbe_sincosf(pw0l, &sd_prev, &cd_prev);
            float sd_cur, cd_cur;
            mbe_sincosf(cw0l, &sd_cur, &cd_cur);
            float s_prev, c_prev;
            mbe_sincosf(prev_mp->PHIl[l], &s_prev, &c_prev);
            float s_cur, c_cur;
            mbe_sincosf(cur_mp->PHIl[l] - (cw0l * (float)N), &s_cur, &c_cur);
            for (n = 0; n < N; n += 4) {
                /* voiced prev */
                mbe_add_voiced_block4(Ss, Ws + n + N, amp_prev, &c_prev, &s_prev, sd_prev, cd_prev);
                /* voiced cur (additive) */
                mbe_add_voiced_block4(Ss, Ws + n, amp_cur, &c_cur, &s_cur, sd_cur, cd_cur);
                Ss += 4;
            }
        }
        /*
      // Alternate phase/amplitude interpolation path (disabled for performance)
      else if ((cur_mp->Vl[l] == 1) || (prev_mp->Vl[l] == 1))
        {
          Ss = aout_buf;
          // eq 137
          deltaphil = cur_mp->PHIl[l] - prev_mp->PHIl[l] - (((pw0 + cw0) * (float) (l * N)) / (float) 2);
          // eq 138
          deltawl = ((float) 1 / (float) N) * (deltaphil - ((float) 2 * M_PI * (int) ((deltaphil + M_PI) / (M_PI * (float) 2))));
          for (n = 0; n < N; n++)
            {
              // eq 136
              thetaln = prev_mp->PHIl[l] + ((pw0l + deltawl) * (float) n) + (((cw0 - pw0) * ((float) (l * n * n)) / (float) (2 * N)));
              // eq 135
              aln = prev_mp->Ml[l] + (((float) n / (float) N) * (cur_mp->Ml[l] - prev_mp->Ml[l]));
              // eq 134
              *Ss = *Ss + (aln * cosf (thetaln));
              Ss++;
            }
        }
*/
        else {
            Ss = aout_buf;
            // init random phase
            for (i = 0; i < uvquality; i++) {
                rphase[i] = mbe_rand_phase();
            }
            // init random phase
            for (i = 0; i < uvquality; i++) {
                rphase2[i] = mbe_rand_phase();
            }
            /* Precompute unvoiced oscillator steps for both prev (pw0) and cur (cw0) */
            float sd_p[64] = {0.0f}, cd_p[64] = {0.0f}, sp[64] = {0.0f}, cp[64] = {0.0f};
            float sd_c[64] = {0.0f}, cd_c[64] = {0.0f}, sc[64] = {0.0f}, cc[64] = {0.0f};
            const float base = (float)l - uvoffset;
            for (i = 0; i < uvquality; i++) {
                float incp = pw0 * (base + ((float)i * uvstep));
                mbe_sincosf(incp, &sd_p[i], &cd_p[i]);
                mbe_sincosf(rphase[i], &sp[i], &cp[i]);
                float incc = cw0 * (base + ((float)i * uvstep));
                mbe_sincosf(incc, &sd_c[i], &cd_c[i]);
                mbe_sincosf(rphase2[i], &sc[i], &cc[i]);
            }
            const float noise_p = (pw0l > uvthreshold) ? ((pw0l - uvthreshold) * uvrand) : 0.0f;
            const float noise_c = (cw0l > uvthreshold) ? ((cw0l - uvthreshold) * uvrand) : 0.0f;
            for (n = 0; n < N; n += 4) {
                float sump[4], sumc[4];
                mbe_unvoiced_mix_block4(cp, sp, cd_p, sd_p, uvquality, noise_p, sump);
                mbe_unvoiced_mix_block4(cc, sc, cd_c, sd_c, uvquality, noise_c, sumc);
                Ss[0] += (uvsine * Ws[n + 0 + N] * prev_mp->Ml[l] * qfactor) * sump[0]
                         + (uvsine * Ws[n + 0] * cur_mp->Ml[l] * qfactor) * sumc[0];
                Ss[1] += (uvsine * Ws[n + 1 + N] * prev_mp->Ml[l] * qfactor) * sump[1]
                         + (uvsine * Ws[n + 1] * cur_mp->Ml[l] * qfactor) * sumc[1];
                Ss[2] += (uvsine * Ws[n + 2 + N] * prev_mp->Ml[l] * qfactor) * sump[2]
                         + (uvsine * Ws[n + 2] * cur_mp->Ml[l] * qfactor) * sumc[2];
                Ss[3] += (uvsine * Ws[n + 3 + N] * prev_mp->Ml[l] * qfactor) * sump[3]
                         + (uvsine * Ws[n + 3] * cur_mp->Ml[l] * qfactor) * sumc[3];
                Ss += 4;
            }
        }
    }
}

/**
 * @brief Synthesize one speech frame into 160 16-bit samples at 8 kHz.
 * @param aout_buf Output buffer of 160 16-bit samples.
 * @param cur_mp   Current parameter set.
 * @param prev_mp  Previous parameter set.
 * @param uvquality Unvoiced synthesis quality (1..64).
 */
void
mbe_synthesizeSpeech(short* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, int uvquality) {
    float float_buf[160];

    mbe_synthesizeSpeechf(float_buf, cur_mp, prev_mp, uvquality);
    mbe_floattoshort(float_buf, aout_buf);
}

/**
 * @brief Convert 160 float samples to clipped/scaled 16-bit PCM.
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 16-bit samples.
 */
/*
 * Runtime-dispatched float->short conversion with SIMD specializations.
 * Keeps public API unchanged while selecting the best implementation at runtime.
 */
/**
 * @brief Portable scalar fallback for float→int16 conversion.
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 int16 samples.
 */
#if defined(MBELIB_ENABLE_SIMD)
static void
mbe_floattoshort_scalar(float* restrict float_buf, short* restrict aout_buf) {
    const float again = 7.0f;
    for (int i = 0; i < 160; i++) {
        float audio = again * float_buf[i];
        if (audio > 32760.0f) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = 32760.0f;
        } else if (audio < -32760.0f) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = -32760.0f;
        }
        aout_buf[i] = (short)(audio);
    }
}

/**
 * @brief SSE2 specialization for float→int16 conversion.
 */
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
static void
mbe_floattoshort_sse2(float* restrict float_buf, short* restrict aout_buf) {
    const __m128 vscale = _mm_set1_ps(7.0f);
    const __m128 vmaxv = _mm_set1_ps(32760.0f);
    const __m128 vminv = _mm_set1_ps(-32760.0f);
    for (int i = 0; i < 160; i += 8) {
        __m128 a = _mm_mul_ps(_mm_loadu_ps(float_buf + i), vscale);
        __m128 b = _mm_mul_ps(_mm_loadu_ps(float_buf + i + 4), vscale);
        a = _mm_min_ps(_mm_max_ps(a, vminv), vmaxv);
        b = _mm_min_ps(_mm_max_ps(b, vminv), vmaxv);
        __m128i ia = _mm_cvttps_epi32(a);
        __m128i ib = _mm_cvttps_epi32(b);
        __m128i packed = _mm_packs_epi32(ia, ib);
        _mm_storeu_si128((__m128i*)(aout_buf + i), packed);
    }
}
#endif

/**
 * @brief NEON specialization for float→int16 conversion.
 */
#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
static void
mbe_floattoshort_neon(float* restrict float_buf, short* restrict aout_buf) {
    const float32x4_t vscale = vdupq_n_f32(7.0f);
    const float32x4_t vmaxv = vdupq_n_f32(32760.0f);
    const float32x4_t vminv = vdupq_n_f32(-32760.0f);
    for (int i = 0; i < 160; i += 8) {
        float32x4_t a = vmulq_f32(vld1q_f32(float_buf + i), vscale);
        float32x4_t b = vmulq_f32(vld1q_f32(float_buf + i + 4), vscale);
        a = vminq_f32(vmaxq_f32(a, vminv), vmaxv);
        b = vminq_f32(vmaxq_f32(b, vminv), vmaxv);
        int32x4_t ia = vcvtq_s32_f32(a);
        int32x4_t ib = vcvtq_s32_f32(b);
        int16x4_t na = vqmovn_s32(ia);
        int16x4_t nb = vqmovn_s32(ib);
        int16x8_t packed = vcombine_s16(na, nb);
        vst1q_s16(aout_buf + i, packed);
    }
}
#endif

typedef void (*mbe_floattoshort_fn)(float*, short*);
static mbe_floattoshort_fn mbe_floattoshort_impl = NULL; /**< Runtime-selected impl pointer. */

/**
 * @brief Initialize runtime dispatch by probing CPU features.
 */
static void
mbe_init_runtime_dispatch(void) {
    if (mbe_floattoshort_impl) {
        return;
    }
    /* Default to scalar */
    mbe_floattoshort_impl = mbe_floattoshort_scalar;

#if defined(__x86_64__) || defined(_M_X64)
    /* SSE2 is guaranteed on x86_64 */
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64)
    mbe_floattoshort_impl = mbe_floattoshort_sse2;
#endif
#elif defined(__i386__) || defined(_M_IX86)
    /* Probe SSE2 on 32-bit x86 */
#if defined(_MSC_VER)
    int regs[4] = {0};
    __cpuid(regs, 1);
    int edx = regs[3];
    if (edx & (1 << 26)) {
/* SSE2 supported */
#if defined(__SSE2__)
        mbe_floattoshort_impl = mbe_floattoshort_sse2;
#endif
    }
#else
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        if (edx & (1u << 26)) {
#if defined(__SSE2__)
            mbe_floattoshort_impl = mbe_floattoshort_sse2;
#endif
        }
    }
#endif
#elif defined(__aarch64__) || defined(_M_ARM64)
    /* NEON is mandatory on AArch64 */
#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    mbe_floattoshort_impl = mbe_floattoshort_neon;
#endif
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    /* Assume NEON if building with NEON intrinsics */
    mbe_floattoshort_impl = mbe_floattoshort_neon;
#endif
}

void
mbe_floattoshort(float* restrict float_buf, short* restrict aout_buf) {
    if (!mbe_floattoshort_impl) {
        mbe_init_runtime_dispatch();
    }
    mbe_floattoshort_impl(float_buf, aout_buf);
}

#else /* MBELIB_ENABLE_SIMD not set: keep scalar implementation */
void
mbe_floattoshort(float* restrict float_buf, short* restrict aout_buf) {
    const float again = 7.0f;
    for (int i = 0; i < 160; i++) {
        float audio = again * float_buf[i];
        if (audio > 32760.0f) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = 32760.0f;
        } else if (audio < -32760.0f) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = -32760.0f;
        }
        aout_buf[i] = (short)(audio);
    }
}
#endif /* MBELIB_ENABLE_SIMD */
