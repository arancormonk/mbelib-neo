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

#include "mbe_adaptive.h"
#include "mbe_compiler.h"
#include "mbe_math.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"
#include "mbelib_const.h"

#define MBE_TWO_PI (2.0f * (float)M_PI)

/* Thread-local PRNG state and helpers (xorshift32) */
static MBE_THREAD_LOCAL uint32_t mbe_rng_state = 0x12345678u;

/* Thread-local FFT plan for unvoiced synthesis */
static MBE_THREAD_LOCAL mbe_fft_plan* mbe_fft_plan_instance = NULL;

/**
 * @brief Get or create the thread-local FFT plan.
 * @return FFT plan for unvoiced synthesis, or NULL on allocation failure.
 */
static mbe_fft_plan*
mbe_get_fft_plan(void) {
    if (MBE_LIKELY(mbe_fft_plan_instance != NULL)) {
        return mbe_fft_plan_instance;
    }
    mbe_fft_plan_instance = mbe_fft_plan_alloc();
    return mbe_fft_plan_instance;
}

void
mbe_setThreadRngSeed(uint32_t seed) {
    if (seed == 0u) {
        seed = 0x6d25357bu; /* avoid zero state */
    }
    mbe_rng_state = seed;
    mbe_seedComfortNoiseRng(seed);
    mbe_seedUnvoicedNoiseLcg(seed);
}

/**
 * @brief Thread-local xorshift32 PRNG step.
 * @return New 32-bit pseudo-random state value (never 0).
 * @note Retained as a lightweight utility; current synthesis paths use
 *       codec-specific generators seeded via mbe_setThreadRngSeed().
 */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((unused))
#endif
static inline uint32_t
mbe_xorshift32(void) {
    uint32_t x = mbe_rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    mbe_rng_state = x ? x : 0x6d25357bu;
    return mbe_rng_state;
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
    __m128 vA = _mm_set1_ps(2.0f * amp); /* JMBE multiplies voiced output by 2.0 */
    __m128 vS = _mm_loadu_ps(Ss);
    vS = _mm_add_ps(vS, _mm_mul_ps(_mm_mul_ps(vC, vW), vA));
    _mm_storeu_ps(Ss, vS);
    return;
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    float32x4_t vC = vld1q_f32(cblk);
    float32x4_t vW = vld1q_f32(W);
    float32x4_t vA = vdupq_n_f32(2.0f * amp); /* JMBE multiplies voiced output by 2.0 */
    float32x4_t vS = vld1q_f32(Ss);
    vS = vaddq_f32(vS, vmulq_f32(vmulq_f32(vC, vW), vA));
    vst1q_f32(Ss, vS);
    return;
#endif
#endif
    /* scalar fallback (JMBE multiplies voiced output by 2.0) */
    Ss[0] += 2.0f * W[0] * amp * cblk[0];
    Ss[1] += 2.0f * W[1] * amp * cblk[1];
    Ss[2] += 2.0f * W[2] * amp * cblk[2];
    Ss[3] += 2.0f * W[3] * amp * cblk[3];
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
 *
 * Uses struct assignment for efficiency - the compiler generates optimal
 * code (typically a single memcpy-like operation) rather than many individual
 * field copies. The mbe_parms struct contains no pointers, so shallow copy
 * is correct and safe.
 */
void
mbe_moveMbeParms(mbe_parms* in, mbe_parms* out) {
    *out = *in;
}

/**
 * @brief Replace current parameters with the last known parameters.
 * @param out Destination parameter set to fill.
 * @param in  Source parameter set from previous frame.
 *
 * Uses struct assignment for efficiency. See mbe_moveMbeParms() for details.
 */
void
mbe_useLastMbeParms(mbe_parms* out, mbe_parms* in) {
    *out = *in;
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
    prev_mp->un = 0;
    prev_mp->w0 = 0.09378;
    prev_mp->L = 30;
    prev_mp->K = 10;
    prev_mp->gamma = (float)0;
    for (l = 0; l <= 56; l++) {
        prev_mp->Ml[l] = (float)0;
        prev_mp->Vl[l] = 0;
        prev_mp->log2Ml[l] = (float)0; // log2 of 1 == 0
        prev_mp->PHIl[l] = (float)0;
        /* JMBE previous-phase arrays start at 0.0f. */
        prev_mp->PSIl[l] = 0.0f;
    }
    prev_mp->repeat = 0;

    /* Initialize adaptive smoothing state */
    prev_mp->localEnergy = MBE_DEFAULT_LOCAL_ENERGY;
    prev_mp->amplitudeThreshold = MBE_DEFAULT_AMPLITUDE_THRESHOLD;
    prev_mp->errorRate = 0.0f;
    prev_mp->errorCountTotal = 0;
    prev_mp->errorCount4 = 0;

    /* Initialize frame repeat state */
    prev_mp->repeatCount = 0;
    prev_mp->mutingThreshold = MBE_MUTING_THRESHOLD_IMBE;

    /* Initialize FFT-based unvoiced synthesis state.
     * Use seed<0 sentinel so first nextBuffer() returns JMBE's initial all-zero
     * 256-sample buffer, then primes with MBE_LCG_DEFAULT_SEED. */
    prev_mp->noiseSeed = -1.0f;
    memset(prev_mp->noiseOverlap, 0, sizeof(prev_mp->noiseOverlap));
    memset(prev_mp->previousUw, 0, sizeof(prev_mp->previousUw));

    mbe_moveMbeParms(prev_mp, cur_mp);
    mbe_moveMbeParms(prev_mp, prev_mp_enhanced);
}

/**
 * @brief Apply spectral amplitude enhancement to the current parameters.
 * @param cur_mp In/out parameter set to enhance.
 *
 * Uses SIMD optimizations for accumulation and scaling loops when available.
 */
void
mbe_spectralAmpEnhance(mbe_parms* cur_mp) {

    float Rm0, Rm1, R2m0, R2m1, Wl[57];
    int l;
    float sum, gamma;
    const int L = cur_mp->L;

    /* Precompute cos(w0 * l) table via recurrence to avoid repeated cosf */
    float cos_tab[57];
    {
        const float w0 = cur_mp->w0;
        float s_step, c_step;
        mbe_sincosf(w0, &s_step, &c_step);
        float c = 1.0f, s = 0.0f;
        for (l = 1; l <= L; ++l) {
            float cn = (c * c_step) - (s * s_step);
            float sn = (s * c_step) + (c * s_step);
            c = cn;
            s = sn;
            cos_tab[l] = c;
        }
    }

    /* Compute Rm0 = sum(Ml^2) and Rm1 = sum(Ml^2 * cos)
     * Uses SIMD accumulation when available */
    Rm0 = 0.0f;
    Rm1 = 0.0f;
#if defined(MBELIB_ENABLE_SIMD)
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
    if (L >= 4) {
        __m128 vRm0 = _mm_setzero_ps();
        __m128 vRm1 = _mm_setzero_ps();
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            __m128 vMl = _mm_loadu_ps(&cur_mp->Ml[l4]);
            __m128 vCos = _mm_loadu_ps(&cos_tab[l4]);
            __m128 vMl2 = _mm_mul_ps(vMl, vMl);
            vRm0 = _mm_add_ps(vRm0, vMl2);
            vRm1 = _mm_add_ps(vRm1, _mm_mul_ps(vMl2, vCos));
        }
        /* Horizontal sum */
        __m128 shuf = _mm_shuffle_ps(vRm0, vRm0, _MM_SHUFFLE(2, 3, 0, 1));
        __m128 sums = _mm_add_ps(vRm0, shuf);
        shuf = _mm_movehl_ps(shuf, sums);
        sums = _mm_add_ss(sums, shuf);
        Rm0 = _mm_cvtss_f32(sums);

        shuf = _mm_shuffle_ps(vRm1, vRm1, _MM_SHUFFLE(2, 3, 0, 1));
        sums = _mm_add_ps(vRm1, shuf);
        shuf = _mm_movehl_ps(shuf, sums);
        sums = _mm_add_ss(sums, shuf);
        Rm1 = _mm_cvtss_f32(sums);

        /* Handle remaining elements */
        for (l = l4; l <= L; l++) {
            const float Ml2 = cur_mp->Ml[l] * cur_mp->Ml[l];
            Rm0 += Ml2;
            Rm1 += Ml2 * cos_tab[l];
        }
    } else
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    if (L >= 4) {
        float32x4_t vRm0 = vdupq_n_f32(0.0f);
        float32x4_t vRm1 = vdupq_n_f32(0.0f);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            float32x4_t vMl = vld1q_f32(&cur_mp->Ml[l4]);
            float32x4_t vCos = vld1q_f32(&cos_tab[l4]);
            float32x4_t vMl2 = vmulq_f32(vMl, vMl);
            vRm0 = vaddq_f32(vRm0, vMl2);
            vRm1 = vmlaq_f32(vRm1, vMl2, vCos);
        }
        /* Horizontal sum */
        float32x2_t sum2 = vadd_f32(vget_low_f32(vRm0), vget_high_f32(vRm0));
        Rm0 = vget_lane_f32(vpadd_f32(sum2, sum2), 0);
        sum2 = vadd_f32(vget_low_f32(vRm1), vget_high_f32(vRm1));
        Rm1 = vget_lane_f32(vpadd_f32(sum2, sum2), 0);

        /* Handle remaining elements */
        for (l = l4; l <= L; l++) {
            const float Ml2 = cur_mp->Ml[l] * cur_mp->Ml[l];
            Rm0 += Ml2;
            Rm1 += Ml2 * cos_tab[l];
        }
    } else
#endif
#endif
    {
        /* Scalar fallback */
        for (l = 1; l <= L; l++) {
            const float Ml2 = cur_mp->Ml[l] * cur_mp->Ml[l];
            Rm0 += Ml2;
            Rm1 += Ml2 * cos_tab[l];
        }
    }

    /* JMBE Algorithm #111 uses pre-enhancement RM0 for local-energy tracking. */
    mbe_setPreEnhRm0(cur_mp, Rm0);

    R2m0 = (Rm0 * Rm0);
    R2m1 = (Rm1 * Rm1);

    /* Compute spectral weight Wl and apply clipping
     * This loop has complex conditionals that are hard to vectorize effectively */
    for (l = 1; l <= L; l++) {
        if (cur_mp->Ml[l] != 0.0f) {
            const float cos_w0l = cos_tab[l];
            Wl[l] = sqrtf(cur_mp->Ml[l])
                    * sqrtf(sqrtf(((float)0.96 * (float)M_PI * ((R2m0 + R2m1) - ((float)2 * Rm0 * Rm1 * cos_w0l)))
                                  / (cur_mp->w0 * Rm0 * (R2m0 - R2m1))));

            if ((8 * l) <= L) {
                /* no-op for low harmonics */
            } else if (Wl[l] > 1.2f) {
                cur_mp->Ml[l] = 1.2f * cur_mp->Ml[l];
            } else if (Wl[l] < 0.5f) {
                cur_mp->Ml[l] = 0.5f * cur_mp->Ml[l];
            } else {
                cur_mp->Ml[l] = Wl[l] * cur_mp->Ml[l];
            }
        }
    }

    /* Compute sum of squared magnitudes for scaling factor
     * Uses SIMD when available */
    sum = 0.0f;
#if defined(MBELIB_ENABLE_SIMD)
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
    if (L >= 4) {
        __m128 vsum = _mm_setzero_ps();
        __m128 sign_mask = _mm_set1_ps(-0.0f);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            __m128 vMl = _mm_loadu_ps(&cur_mp->Ml[l4]);
            vMl = _mm_andnot_ps(sign_mask, vMl); /* fabs */
            vsum = _mm_add_ps(vsum, _mm_mul_ps(vMl, vMl));
        }
        /* Horizontal sum */
        __m128 shuf = _mm_shuffle_ps(vsum, vsum, _MM_SHUFFLE(2, 3, 0, 1));
        __m128 sums = _mm_add_ps(vsum, shuf);
        shuf = _mm_movehl_ps(shuf, sums);
        sums = _mm_add_ss(sums, shuf);
        sum = _mm_cvtss_f32(sums);

        /* Handle remaining elements */
        for (l = l4; l <= L; l++) {
            float M = cur_mp->Ml[l];
            if (M < 0.0f) {
                M = -M;
            }
            sum += (M * M);
        }
    } else
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    if (L >= 4) {
        float32x4_t vsum = vdupq_n_f32(0.0f);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            float32x4_t vMl = vld1q_f32(&cur_mp->Ml[l4]);
            vMl = vabsq_f32(vMl); /* fabs */
            vsum = vmlaq_f32(vsum, vMl, vMl);
        }
        /* Horizontal sum */
        float32x2_t sum2 = vadd_f32(vget_low_f32(vsum), vget_high_f32(vsum));
        sum = vget_lane_f32(vpadd_f32(sum2, sum2), 0);

        /* Handle remaining elements */
        for (l = l4; l <= L; l++) {
            float M = cur_mp->Ml[l];
            if (M < 0.0f) {
                M = -M;
            }
            sum += (M * M);
        }
    } else
#endif
#endif
    {
        /* Scalar fallback */
        for (l = 1; l <= L; l++) {
            float M = cur_mp->Ml[l];
            if (M < 0.0f) {
                M = -M;
            }
            sum += (M * M);
        }
    }

    if (sum == 0.0f) {
        gamma = 1.0f;
    } else {
        gamma = sqrtf(Rm0 / sum);
    }

    /* Apply scaling factor to all bands
     * Uses SIMD when available */
#if defined(MBELIB_ENABLE_SIMD)
#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
    if (L >= 4) {
        __m128 vgamma = _mm_set1_ps(gamma);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            __m128 vMl = _mm_loadu_ps(&cur_mp->Ml[l4]);
            vMl = _mm_mul_ps(vMl, vgamma);
            _mm_storeu_ps(&cur_mp->Ml[l4], vMl);
        }
        /* Handle remaining elements */
        for (l = l4; l <= L; l++) {
            cur_mp->Ml[l] = gamma * cur_mp->Ml[l];
        }
    } else
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64)
    if (L >= 4) {
        float32x4_t vgamma = vdupq_n_f32(gamma);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            float32x4_t vMl = vld1q_f32(&cur_mp->Ml[l4]);
            vMl = vmulq_f32(vMl, vgamma);
            vst1q_f32(&cur_mp->Ml[l4], vMl);
        }
        /* Handle remaining elements */
        for (l = l4; l <= L; l++) {
            cur_mp->Ml[l] = gamma * cur_mp->Ml[l];
        }
    } else
#endif
#endif
    {
        /* Scalar fallback */
        for (l = 1; l <= L; l++) {
            cur_mp->Ml[l] = gamma * cur_mp->Ml[l];
        }
    }
}

/* JMBE float-domain soft clip translated to this library's float scale. */
#define MBE_AUDIO_SOFT_CLIP_FLOAT   ((32767.0f * 0.95f) / 7.0f)
#define MBE_TONE_PHASE_SCALE        4294967296.0
#define MBE_TONE_PHASE_RAD_PER_TICK ((2.0 * M_PI) / MBE_TONE_PHASE_SCALE)

static inline float
mbe_clipFloatSample(float sample) {
    if (sample > MBE_AUDIO_SOFT_CLIP_FLOAT) {
        return MBE_AUDIO_SOFT_CLIP_FLOAT;
    }
    if (sample < -MBE_AUDIO_SOFT_CLIP_FLOAT) {
        return -MBE_AUDIO_SOFT_CLIP_FLOAT;
    }
    return sample;
}

static void
mbe_clipFloatBuffer(float* samples, int count) {
    for (int i = 0; i < count; i++) {
        samples[i] = mbe_clipFloatSample(samples[i]);
    }
}

static inline uint32_t
mbe_tonePhaseStep(double freq_hz) {
    double step = (freq_hz / 8000.0) * MBE_TONE_PHASE_SCALE;
    if (step <= 0.0) {
        return 0u;
    }
    return (uint32_t)(step + 0.5);
}

static inline float
mbe_toneSampleFromPhase(uint32_t phase) {
    float angle = (float)(((double)phase * MBE_TONE_PHASE_RAD_PER_TICK) - (M_PI / 2.0));
    return sinf(angle);
}

static int
mbe_lookupToneFreqs(int tone_id, float* freq1, float* freq2) {
    *freq1 = 0.0f;
    *freq2 = 0.0f;

    switch (tone_id) {
        case 5:
            *freq1 = 156.25f;
            *freq2 = *freq1;
            break;
        case 6:
            *freq1 = 187.5f;
            *freq2 = *freq1;
            break;
        case 128:
            *freq1 = 1336.0f;
            *freq2 = 941.0f;
            break;
        case 129:
            *freq1 = 1209.0f;
            *freq2 = 697.0f;
            break;
        case 130:
            *freq1 = 1336.0f;
            *freq2 = 697.0f;
            break;
        case 131:
            *freq1 = 1477.0f;
            *freq2 = 697.0f;
            break;
        case 132:
            *freq1 = 1209.0f;
            *freq2 = 770.0f;
            break;
        case 133:
            *freq1 = 1336.0f;
            *freq2 = 770.0f;
            break;
        case 134:
            *freq1 = 1477.0f;
            *freq2 = 770.0f;
            break;
        case 135:
            *freq1 = 1209.0f;
            *freq2 = 852.0f;
            break;
        case 136:
            *freq1 = 1336.0f;
            *freq2 = 852.0f;
            break;
        case 137:
            *freq1 = 1477.0f;
            *freq2 = 852.0f;
            break;
        case 138:
            *freq1 = 1633.0f;
            *freq2 = 697.0f;
            break;
        case 139:
            *freq1 = 1633.0f;
            *freq2 = 770.0f;
            break;
        case 140:
            *freq1 = 1633.0f;
            *freq2 = 852.0f;
            break;
        case 141:
            *freq1 = 1633.0f;
            *freq2 = 941.0f;
            break;
        case 142:
            *freq1 = 1209.0f;
            *freq2 = 941.0f;
            break;
        case 143:
            *freq1 = 1477.0f;
            *freq2 = 941.0f;
            break;
        case 144:
            *freq1 = 1162.0f;
            *freq2 = 820.0f;
            break;
        case 145:
            *freq1 = 1052.0f;
            *freq2 = 606.0f;
            break;
        case 146:
            *freq1 = 1162.0f;
            *freq2 = 606.0f;
            break;
        case 147:
            *freq1 = 1279.0f;
            *freq2 = 606.0f;
            break;
        case 148:
            *freq1 = 1052.0f;
            *freq2 = 672.0f;
            break;
        case 149:
            *freq1 = 1162.0f;
            *freq2 = 672.0f;
            break;
        case 150:
            *freq1 = 1279.0f;
            *freq2 = 672.0f;
            break;
        case 151:
            *freq1 = 1052.0f;
            *freq2 = 743.0f;
            break;
        case 152:
            *freq1 = 1162.0f;
            *freq2 = 743.0f;
            break;
        case 153:
            *freq1 = 1279.0f;
            *freq2 = 743.0f;
            break;
        case 154:
            *freq1 = 1430.0f;
            *freq2 = 606.0f;
            break;
        case 155:
            *freq1 = 1430.0f;
            *freq2 = 672.0f;
            break;
        case 156:
            *freq1 = 1430.0f;
            *freq2 = 743.0f;
            break;
        case 157:
            *freq1 = 1430.0f;
            *freq2 = 820.0f;
            break;
        case 158:
            *freq1 = 1052.0f;
            *freq2 = 820.0f;
            break;
        case 159:
            *freq1 = 1279.0f;
            *freq2 = 820.0f;
            break;
        case 160:
            *freq1 = 440.0f;
            *freq2 = 350.0f;
            break;
        case 161:
            *freq1 = 480.0f;
            *freq2 = 440.0f;
            break;
        case 162:
            *freq1 = 620.0f;
            *freq2 = 480.0f;
            break;
        case 163:
            *freq1 = 490.0f;
            *freq2 = 350.0f;
            break;
        case 255: break;
        default:
            if ((tone_id >= 7) && (tone_id <= 122)) {
                *freq1 = 31.25f * (float)tone_id;
                *freq2 = *freq1;
            }
            break;
    }

    return ((*freq1 > 0.0f) || (*freq2 > 0.0f));
}

static void
mbe_renderTonef(float* aout_buf, mbe_parms* cur_mp, float freq1, float freq2, int amplitude_id) {
    if (!aout_buf || !cur_mp || freq1 <= 0.0f) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }

    const int dual_tone = (freq2 > 0.0f) && (fabsf(freq2 - freq1) > 1e-6f);
    const float gain = (((amplitude_id < 0) ? 0.0f : (float)amplitude_id) / 127.0f) * MBE_AUDIO_SOFT_CLIP_FLOAT;
    const uint32_t step1 = mbe_tonePhaseStep((double)freq1);
    const uint32_t step2 = dual_tone ? mbe_tonePhaseStep((double)freq2) : 0u;
    uint32_t phase1 = (uint32_t)cur_mp->swn;
    uint32_t phase2 = (uint32_t)cur_mp->un;

    for (int n = 0; n < 160; n++) {
        phase1 += step1;
        float s1 = mbe_toneSampleFromPhase(phase1);

        if (dual_tone) {
            phase2 += step2;
            float s2 = mbe_toneSampleFromPhase(phase2);
            aout_buf[n] = (0.5f * gain * s1) + (0.5f * gain * s2);
        } else {
            aout_buf[n] = gain * s1;
        }
    }

    cur_mp->swn = (int)phase1;
    cur_mp->un = (int)phase2;
}

/**
 * @brief Synthesize a tone frame into 160 float samples at 8 kHz.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49) providing tone indices.
 * @param cur_mp   Current parameter set (tone synthesis state).
 */
void
mbe_synthesizeTonef(float* aout_buf, char* ambe_d, mbe_parms* cur_mp) {
    int i;

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

    float freq1, freq2;
    (void)ID0;
    (void)ID2;
    (void)ID3; /* parsed for potential validation */
    (void)ID4;

#ifdef DISABLE_AMBE_TONES // generate silence if tones disabled
    mbe_synthesizeSilencef(aout_buf);
    return;
#endif

    if (!mbe_lookupToneFreqs(ID1, &freq1, &freq2)) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }

    mbe_renderTonef(aout_buf, cur_mp, freq1, freq2, AD);
}

/**
 * @brief Synthesize a D-STAR style tone into 160 float samples.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49).
 * @param cur_mp   Current parameter set.
 * @param ID1      Tone index selector.
 */
void
mbe_synthesizeTonefdstar(float* aout_buf, char* ambe_d, mbe_parms* cur_mp, int ID1) {
    int AD = 103; /* JMBE nominal D-STAR tone amplitude */
    float freq1 = 0, freq2 = 0;
    (void)ambe_d;

#ifdef DISABLE_AMBE_TONES // generate silence if tones disabled
    mbe_synthesizeSilencef(aout_buf);
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
                freq1 = 31.25f * (float)ID1;
                freq2 = freq1;
            }
    }

    if (freq1 <= 0.0f) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }

    mbe_renderTonef(aout_buf, cur_mp, freq1, freq2, AD);
}

/**
 * @brief Write 160 float samples of silence.
 * @param aout_buf Output buffer of 160 float samples.
 */
void
mbe_synthesizeSilencef(float* aout_buf) {
    if (aout_buf == NULL) {
        return;
    }
    memset(aout_buf, 0, 160 * sizeof(*aout_buf));
}

/**
 * @brief Write 160 16-bit samples of silence.
 * @param aout_buf Output buffer of 160 16-bit samples.
 */
void
mbe_synthesizeSilence(short* aout_buf) {
    if (aout_buf == NULL) {
        return;
    }
    memset(aout_buf, 0, 160 * sizeof(*aout_buf));
}

/**
 * @brief Synthesize one speech frame into 160 float samples at 8 kHz.
 *
 * Uses FFT-based unvoiced synthesis (JMBE Algorithms #117-126) for
 * high-quality unvoiced audio with proper WOLA frame blending.
 *
 * @param aout_buf Output buffer of 160 float samples.
 * @param cur_mp   Current parameter set.
 * @param prev_mp  Previous parameter set.
 * @param uvquality Unvoiced synthesis quality (ignored, kept for API compatibility).
 */
/* JMBE-compatible white noise scalar for phase calculation: 2*PI / 53125 */
#define MBE_WHITE_NOISE_SCALAR (2.0f * (float)M_PI / 53125.0f)

void
mbe_synthesizeSpeechf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, int uvquality) {

    int l, n, maxl;
    float* Ss;
    int numUv;
    float cw0, pw0, cw0l, pw0l;

    const int N = 160;

    /* Silence unused parameter warning - uvquality is kept for API compatibility */
    (void)uvquality;

    /* Apply adaptive smoothing (Algorithms #111-116) before muting checks.
     * JMBE computes/update smoothing state during parameter preparation even
     * for frames that are subsequently muted. */
    mbe_applyAdaptiveSmoothing(cur_mp, prev_mp);

    /* Frame muting:
     * - JMBE IMBE path mutes on max repeats OR error-rate threshold.
     * - JMBE AMBE path mutes on max repeats only (error-rate muting is not applied in AMBE synth path). */
    int mute_on_error_rate = (fabsf(cur_mp->mutingThreshold - MBE_MUTING_THRESHOLD_AMBE) > 1e-6f);
    if (mbe_isMaxFrameRepeat(cur_mp) || (mute_on_error_rate && mbe_requiresMuting(cur_mp))) {
        /* Muted frames output comfort noise while preserving model progression. */
        mbe_synthesizeComfortNoisef(aout_buf);
        return;
    }

    /* Algorithm #117: Generate 256 white noise samples FIRST (JMBE-compatible)
     * This buffer is used for both phase calculation and unvoiced synthesis */
    float noise_buffer[256];
    mbe_generate_noise_with_overlap(noise_buffer, &cur_mp->noiseSeed, cur_mp->noiseOverlap);

    /* Count number of unvoiced bands */
    numUv = 0;
    for (l = 0; l <= cur_mp->L; l++) {
        if (cur_mp->Vl[l] == 0) {
            numUv++;
        }
    }

    cw0 = cur_mp->w0;
    pw0 = prev_mp->w0;

    /* Initialize output buffer to zero */
    memset(aout_buf, 0, N * sizeof(float));

    /* eq 128 and 129: Handle different L values between frames */
    if (cur_mp->L > prev_mp->L) {
        maxl = cur_mp->L;
        for (l = prev_mp->L + 1; l <= maxl; l++) {
            prev_mp->Ml[l] = 0.0f;
            prev_mp->Vl[l] = 1;
        }
    } else {
        maxl = prev_mp->L;
        for (l = cur_mp->L + 1; l <= maxl; l++) {
            cur_mp->Ml[l] = 0.0f;
            cur_mp->Vl[l] = 1;
        }
    }

    /* Update phase from eq 139, 140
     * JMBE-compatible: use noise_buffer[l] for phase randomization instead of separate RNG */
    for (l = 1; l <= 56; l++) {
        /* JMBE parity: wrap previous voiced phase to [0, 2PI) each frame before advancing. */
        float prev_psil_wrapped = fmodf(prev_mp->PSIl[l], MBE_TWO_PI);
        if (prev_psil_wrapped < 0.0f) {
            prev_psil_wrapped += MBE_TWO_PI;
        }
        prev_mp->PSIl[l] = prev_psil_wrapped;

        cur_mp->PSIl[l] = prev_psil_wrapped + ((pw0 + cw0) * ((float)(l * N) / 2.0f));
        if (l <= (cur_mp->L / 4)) {
            cur_mp->PHIl[l] = cur_mp->PSIl[l];
        } else {
            /* JMBE Algorithm #140: Use noise sample for phase jitter
             * pl = (2*PI/53125) * u[l] - PI maps noise to [-PI, +PI] */
            float pl = (MBE_WHITE_NOISE_SCALAR * noise_buffer[l]) - (float)M_PI;
            cur_mp->PHIl[l] = cur_mp->PSIl[l] + (((float)numUv * pl) / (float)cur_mp->L);
        }
    }

    /* Synthesize voiced components
     * Use phase/amplitude interpolation (Algorithms #134-138) for low harmonics
     * when pitch is stable, otherwise use windowed oscillator approach
     */
    for (l = 1; l <= maxl; l++) {
        cw0l = cw0 * (float)l;
        pw0l = pw0 * (float)l;

        /* Check if this band has any voiced component */
        int cur_voiced = (cur_mp->Vl[l] == 1);
        int prev_voiced = (prev_mp->Vl[l] == 1);

        if (cur_voiced || prev_voiced) {
            /* Use interpolation for low harmonics (l < 8) with stable pitch
             * This produces smoother transitions per JMBE Algorithms #134-138 */
            int use_interpolation = (l < 8) && cur_voiced && prev_voiced && (fabsf(cw0 - pw0) < (0.1f * cw0));

            if (use_interpolation) {
                Ss = aout_buf;

                /* Algorithm #137: Phase deviation */
                float deltaphil = cur_mp->PHIl[l] - prev_mp->PHIl[l] - (((pw0 + cw0) * (float)(l * N)) / 2.0f);

                /* Algorithm #138: Phase deviation rate (wrap to [-pi, pi]) */
                float deltawl =
                    (1.0f / (float)N)
                    * (deltaphil - (2.0f * (float)M_PI * floorf((deltaphil + (float)M_PI) / (2.0f * (float)M_PI))));

                for (n = 0; n < N; n++) {
                    /* Algorithm #136: Phase function with quadratic frequency interpolation */
                    float thetaln = prev_mp->PHIl[l] + ((pw0l + deltawl) * (float)n)
                                    + (((cw0 - pw0) * (float)(l * n * n)) / (float)(2 * N));

                    /* Algorithm #135: Linear amplitude interpolation */
                    float aln = prev_mp->Ml[l] + (((float)n / (float)N) * (cur_mp->Ml[l] - prev_mp->Ml[l]));

                    /* Algorithm #134: Synthesize sample (JMBE multiplies by 2.0) */
                    *Ss += 2.0f * aln * cosf(thetaln);
                    Ss++;
                }
            } else {
                /* Windowed oscillator approach for higher harmonics or pitch changes */
                Ss = aout_buf;

                /* Synthesize previous voiced component (fading out) */
                if (prev_voiced) {
                    const float amp_prev = prev_mp->Ml[l];
                    float sd_prev, cd_prev;
                    mbe_sincosf(pw0l, &sd_prev, &cd_prev);
                    float s_prev, c_prev;
                    mbe_sincosf(prev_mp->PHIl[l], &s_prev, &c_prev);

                    for (n = 0; n < N; n += 4) {
                        mbe_add_voiced_block4(Ss, Ws + n + N, amp_prev, &c_prev, &s_prev, sd_prev, cd_prev);
                        Ss += 4;
                    }
                }

                /* Synthesize current voiced component (fading in) */
                if (cur_voiced) {
                    Ss = aout_buf;
                    const float amp_cur = cur_mp->Ml[l];
                    float sd_cur, cd_cur;
                    mbe_sincosf(cw0l, &sd_cur, &cd_cur);
                    float s_cur, c_cur;
                    mbe_sincosf(cur_mp->PHIl[l] - (cw0l * (float)N), &s_cur, &c_cur);

                    for (n = 0; n < N; n += 4) {
                        mbe_add_voiced_block4(Ss, Ws + n, amp_cur, &c_cur, &s_cur, sd_cur, cd_cur);
                        Ss += 4;
                    }
                }
            }
        }
    }

    /* Synthesize unvoiced components using FFT method (JMBE Algorithms #117-126)
     * Use the same noise buffer that was used for phase calculation */
    mbe_fft_plan* plan = mbe_get_fft_plan();
    if (plan) {
        mbe_synthesizeUnvoicedFFTWithNoise(aout_buf, cur_mp, prev_mp, plan, noise_buffer);
    }

    /* Match JMBE float-path soft clipping semantics for synthesized speech. */
    mbe_clipFloatBuffer(aout_buf, N);
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
    /* JMBE-compatible soft clipping at 95% of maximum amplitude */
    const float again = 7.0f;
    const float max_amplitude = 32767.0f * 0.95f; /* ~31128.65 */
    for (int i = 0; i < 160; i++) {
        float audio = again * float_buf[i];
        if (audio > max_amplitude) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = max_amplitude;
        } else if (audio < -max_amplitude) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = -max_amplitude;
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
    /* JMBE-compatible soft clipping at 95% of maximum amplitude */
    const __m128 vscale = _mm_set1_ps(7.0f);
    const __m128 vmaxv = _mm_set1_ps(32767.0f * 0.95f);
    const __m128 vminv = _mm_set1_ps(-32767.0f * 0.95f);
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
    /* JMBE-compatible soft clipping at 95% of maximum amplitude */
    const float32x4_t vscale = vdupq_n_f32(7.0f);
    const float32x4_t vmaxv = vdupq_n_f32(32767.0f * 0.95f);
    const float32x4_t vminv = vdupq_n_f32(-32767.0f * 0.95f);
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
    if (MBE_UNLIKELY(!mbe_floattoshort_impl)) {
        mbe_init_runtime_dispatch();
    }
    mbe_floattoshort_impl(float_buf, aout_buf);
}

#else /* MBELIB_ENABLE_SIMD not set: keep scalar implementation */
void
mbe_floattoshort(float* restrict float_buf, short* restrict aout_buf) {
    /* JMBE-compatible soft clipping at 95% of maximum amplitude
     * This provides headroom and prevents harsh clipping artifacts */
    const float again = 7.0f;
    const float max_amplitude = 32767.0f * 0.95f; /* ~31128.65 */
    for (int i = 0; i < 160; i++) {
        float audio = again * float_buf[i];
        if (audio > max_amplitude) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = max_amplitude;
        } else if (audio < -max_amplitude) {
#ifdef MBE_DEBUG
            fprintf(stderr, "audio clip: %f\n", audio);
#endif
            audio = -max_amplitude;
        }
        aout_buf[i] = (short)(audio);
    }
}
#endif /* MBELIB_ENABLE_SIMD */
