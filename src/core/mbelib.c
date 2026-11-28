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
 * `MBELIB_ENABLE_SIMD` on determinism and ordering.
 *
 * @{ 
 */

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "mbe_compiler.h"
#include "mbe_voiced.h"
#if defined(MBELIB_ENABLE_SIMD)
#include "mbe_simd.h"
#if defined(MBE_SIMD_TARGET_SSE2)
#include <emmintrin.h>
#endif
#if defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h>
#endif
#endif

#include "mbe_adaptive.h"
#include "mbe_low_band.h"
#include "mbe_math.h"
#include "mbe_result.h"
#include "mbe_subframe.h"
#include "mbe_tone.h"
#include "mbe_unvoiced_fft.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"
#include "mbelib-neo/version.h"
#include "mbelib_const.h"

#define MBE_TWO_PI (2.0f * (float)M_PI)

/* Thread-local PRNG state and helpers (xorshift32) */
static MBE_THREAD_LOCAL uint32_t mbe_rng_state = 0x12345678u;

/* Thread-local FFT plan for unvoiced synthesis */
static MBE_THREAD_LOCAL mbe_fft_plan* mbe_fft_plan_instance = NULL;

void
mbe_initProcessResult(mbe_process_result* result) {
    if (!result) {
        return;
    }
    memset(result, 0, sizeof(*result));
}

void
mbe_formatProcessResult(char* str, size_t size, const mbe_process_result* result) {
    static const struct {
        unsigned flag;
        char marker;
    } markers[] = {
        {MBE_PROCESS_FLAG_ERASURE, 'E'},
        {MBE_PROCESS_FLAG_TONE, 'T'},
        {MBE_PROCESS_FLAG_REPEAT, 'R'},
        {MBE_PROCESS_FLAG_MUTE, 'M'},
    };

    if (!str || size == 0u) {
        return;
    }

    size_t pos = 0u;
    int total_errors = result ? result->total_errors : 0;
    if (total_errors < 0) {
        total_errors = 0;
    }

    for (int i = 0; i < total_errors && pos + 1u < size; ++i) {
        str[pos++] = '=';
    }

    if (result) {
        const unsigned flags = result->flags;
        for (size_t i = 0u; i < sizeof(markers) / sizeof(markers[0]) && pos + 1u < size; ++i) {
            if ((flags & markers[i].flag) != 0u) {
                str[pos++] = markers[i].marker;
            }
        }
    }
    str[pos] = '\0';
}

static uint8_t
mbe_clamp_u8_int(int value) {
    if (value < 0) {
        return 0u;
    }
    if (value > 255) {
        return 255u;
    }
    return (uint8_t)value;
}

mbe_soft_bit
mbe_softBitFromHard(int bit, uint8_t reliability) {
    mbe_soft_bit soft;
    soft.bit = (uint8_t)(bit ? 1u : 0u);
    soft.reliability = reliability;
    return soft;
}

mbe_soft_bit
mbe_softBitFromLlr(int16_t llr) {
    int mag = (llr < 0) ? -(int)llr : (int)llr;
    mbe_soft_bit soft;
    soft.bit = (uint8_t)((llr > 0) ? 1u : 0u);
    soft.reliability = mbe_clamp_u8_int(mag);
    return soft;
}

int
mbe_softBitsFromHard(const char* bits, mbe_soft_bit* soft, size_t count, uint8_t reliability) {
    if (!soft) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    int ret = mbe_validate_bits(bits, count);
    if (ret < 0) {
        return ret;
    }
    for (size_t i = 0u; i < count; ++i) {
        soft[i] = mbe_softBitFromHard(bits[i], reliability);
    }
    return 0;
}

int
mbe_softBitsFromLlr(const int16_t* llr, mbe_soft_bit* soft, size_t count) {
    if (!llr || !soft) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (size_t i = 0u; i < count; ++i) {
        soft[i] = mbe_softBitFromLlr(llr[i]);
    }
    return 0;
}

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

#if !defined(MBE_VOICED_KERNELS)
/*
 * Targets without vector voiced kernels keep the historical voiced synthesis
 * unchanged: these helpers and the render loops below are the original scalar
 * code, so its rounding points and expression shapes, and thus its PCM, stay
 * the same under x87 excess precision and fast-math reassociation too.
 */
static inline void
mbe_fill_voiced_cos_block4(float* restrict cblk, float* restrict c, float* restrict s, float sd, float cd) {
    float cc = *c;
    float ss = *s;

    for (int k = 0; k < 4; ++k) {
        cblk[k] = cc;
        float cpn = (cc * cd) - (ss * sd);
        float spn = (ss * cd) + (cc * sd);
        cc = cpn;
        ss = spn;
    }

    *c = cc;
    *s = ss;
}

/**
 * @brief Add four voiced samples to the output using oscillator recurrence.
 *
 * Generates four consecutive cosine samples from the current oscillator state,
 * multiplies by window `W[0..3]` and gain `gain`, and accumulates into
 * `Ss[0..3]`. `gain` already includes the TIA-102.BABA eq. 136 `2.0f` voiced multiplier.
 *
 * @param Ss Output sample pointer (adds to Ss[0..3]).
 * @param W  Window values for these four samples.
 * @param gain Voiced gain for this band/component, including the `2.0f` factor.
 * @param c  In/out cosine register.
 * @param s  In/out sine register.
 * @param sd sin(Δθ) per step.
 * @param cd cos(Δθ) per step.
 */
/** @internal @ingroup mbe_internal */
static inline void
mbe_add_voiced_block4(float* restrict Ss, const float* restrict W, float gain, float* restrict c, float* restrict s,
                      float sd, float cd) {
    float cblk[4];
    mbe_fill_voiced_cos_block4(cblk, c, s, sd, cd);

    Ss[0] += gain * W[0] * cblk[0];
    Ss[1] += gain * W[1] * cblk[1];
    Ss[2] += gain * W[2] * cblk[2];
    Ss[3] += gain * W[3] * cblk[3];
}

/**
 * @brief Add previous/current voiced contributions in one 4-sample block.
 *
 * This preserves the historical add order (`prev`, then `cur`) while avoiding
 * a second load/store of the output buffer when both voiced components are
 * active for the same harmonic.
 */
typedef struct mbe_voiced_block_component {
    const float* W;
    float gain;
    float* c;
    float* s;
    float sd;
    float cd;
} mbe_voiced_block_component;

/** @internal @ingroup mbe_internal */
static inline void
mbe_add_voiced_dual_block4(float* restrict Ss, const mbe_voiced_block_component* restrict prev,
                           const mbe_voiced_block_component* restrict cur) {
    float prev_cblk[4];
    float cur_cblk[4];

    mbe_fill_voiced_cos_block4(prev_cblk, prev->c, prev->s, prev->sd, prev->cd);
    mbe_fill_voiced_cos_block4(cur_cblk, cur->c, cur->s, cur->sd, cur->cd);

    Ss[0] += prev->gain * prev->W[0] * prev_cblk[0];
    Ss[1] += prev->gain * prev->W[1] * prev_cblk[1];
    Ss[2] += prev->gain * prev->W[2] * prev_cblk[2];
    Ss[3] += prev->gain * prev->W[3] * prev_cblk[3];
    Ss[0] += cur->gain * cur->W[0] * cur_cblk[0];
    Ss[1] += cur->gain * cur->W[1] * cur_cblk[1];
    Ss[2] += cur->gain * cur->W[2] * cur_cblk[2];
    Ss[3] += cur->gain * cur->W[3] * cur_cblk[3];
}
#endif

/** @} */ /* end of mbe_internal */

const char*
mbe_versionString(void) {
    return MBELIB_VERSION;
}

/**
 * @brief Copy MBE parameter set from input to output.
 * @param source_mp      Source parameter set.
 * @param destination_mp Destination parameter set.
 *
 * Uses struct assignment for efficiency - the compiler generates optimal
 * code (typically a single memcpy-like operation) rather than many individual
 * field copies. The mbe_parms struct contains no pointers, so shallow copy
 * is correct and safe.
 */
void
mbe_moveMbeParms(const mbe_parms* source_mp, mbe_parms* destination_mp) {
    if (!source_mp || !destination_mp) {
        return;
    }
    *destination_mp = *source_mp;
}

/**
 * @brief Replace current parameters with the last known parameters.
 * @param cur_mp  Destination parameter set to fill.
 * @param prev_mp Source parameter set from previous frame.
 *
 * Uses struct assignment for efficiency. See mbe_moveMbeParms() for details.
 */
void
mbe_useLastMbeParms(mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    if (!cur_mp || !prev_mp) {
        return;
    }
    *cur_mp = *prev_mp;
}

/**
 * @brief Initialize MBE parameter state for decoding and synthesis.
 * @param cur_mp Output: current parameter state.
 * @param prev_mp Output: previous parameter state (reset to defaults).
 * @param prev_mp_enhanced Output: enhanced previous parameter state.
 */
void
mbe_resetImbeModel(mbe_parms* mp) {
    if (mp == NULL) {
        return;
    }
    mp->w0 = (float)(4.0 * M_PI / (134.0 + 39.5));
    mp->L = (int)floor(0.9254 * floor(M_PI / (double)mp->w0 + 0.25));
    mp->K = 12;
    mp->gamma = 0.0f;
    for (int l = 0; l <= 56; ++l) {
        mp->Ml[l] = 1.0f;
        mp->Vl[l] = 0;
        mp->log2Ml[l] = 0.0f;
    }
    mp->swn = 0;
    mp->tonePhase = 0;
    mp->repeatCount = 0;
    /* The mark belongs to the frame the model came from. */
    mbe_set_frame_mark(mp, MBE_FRAME_MARK_NONE);
    mp->localEnergy = MBE_DEFAULT_LOCAL_ENERGY;
    mp->amplitudeThreshold = MBE_DEFAULT_AMPLITUDE_THRESHOLD;
    mp->mutingThreshold = MBE_MUTING_THRESHOLD_IMBE;
}

void
mbe_initMbeParms(mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    if (cur_mp == NULL || prev_mp == NULL || prev_mp_enhanced == NULL) {
        return;
    }
    mbe_parms initial;
    memset(&initial, 0, sizeof(initial));
    mbe_resetImbeModel(&initial);
    initial.noiseSeed = -1.0f;
    memcpy(cur_mp, &initial, sizeof(initial));
    memcpy(prev_mp, &initial, sizeof(initial));
    memcpy(prev_mp_enhanced, &initial, sizeof(initial));
}

static void
mbe_precompute_harmonic_cosines(float cos_tab[57], int L, float w0) {
    float s_step, c_step;
    mbe_sincosf(w0, &s_step, &c_step);
    float c = 1.0f, s = 0.0f;
    for (int l = 1; l <= L; ++l) {
        float cn = (c * c_step) - (s * s_step);
        float sn = (s * c_step) + (c * s_step);
        c = cn;
        s = sn;
        cos_tab[l] = c;
    }
}

static void
mbe_accumulate_spectral_energy(const mbe_parms* cur_mp, const float cos_tab[57], float* out_rm0, float* out_rm1) {
    const int L = cur_mp->L;
    float Rm0 = 0.0f;
    float Rm1 = 0.0f;
    int l;

#if defined(MBELIB_ENABLE_SIMD)
#if defined(MBE_SIMD_TARGET_SSE2)
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

        for (l = l4; l <= L; l++) {
            const float Ml2 = cur_mp->Ml[l] * cur_mp->Ml[l];
            Rm0 += Ml2;
            Rm1 += Ml2 * cos_tab[l];
        }
    } else
#elif defined(MBE_SIMD_TARGET_NEON)
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
        float32x2_t sum2 = vadd_f32(vget_low_f32(vRm0), vget_high_f32(vRm0));
        Rm0 = vget_lane_f32(vpadd_f32(sum2, sum2), 0);
        sum2 = vadd_f32(vget_low_f32(vRm1), vget_high_f32(vRm1));
        Rm1 = vget_lane_f32(vpadd_f32(sum2, sum2), 0);

        for (l = l4; l <= L; l++) {
            const float Ml2 = cur_mp->Ml[l] * cur_mp->Ml[l];
            Rm0 += Ml2;
            Rm1 += Ml2 * cos_tab[l];
        }
    } else
#endif
#endif
    {
        for (l = 1; l <= L; l++) {
            const float Ml2 = cur_mp->Ml[l] * cur_mp->Ml[l];
            Rm0 += Ml2;
            Rm1 += Ml2 * cos_tab[l];
        }
    }

    *out_rm0 = Rm0;
    *out_rm1 = Rm1;
}

/*
 * W_l = sqrt(M_l) (0.96 pi (Rm0^2 + Rm1^2 - 2 Rm0 Rm1 cos(w0 l)) / (w0 Rm0 (Rm0^2 - Rm1^2)))^(1/4),
 * evaluated with r = Rm1 / Rm0 and M_l / sqrt(Rm0), which are at most 1: the
 * terms in Rm0^3 overflow float for loud models (a D-STAR frame held at b2 63
 * reaches Rm0 1e19 and played silence) and underflow for quiet ones. The
 * numerator is at least (1 - |r|)^2. A model whose r or a numerator comes
 * within float's error of 0 (nearly all its energy at cos(w0 l) = +-1, which
 * no decoded model reaches) is enhanced in double precision instead
 * (mbe_spectral_float_unresolved()).
 */
static void
mbe_apply_spectral_weights(mbe_parms* cur_mp, const float cos_tab[57], float Rm0, float Rm1) {
    float Wl[57];
    const int L = cur_mp->L;

    if (!(Rm0 > 0.0f) || !(cur_mp->w0 > 0.0f)) {
        return;
    }
    const float r = Rm1 / Rm0;
    const float shape = ((float)0.96 * (float)M_PI) / (cur_mp->w0 * ((float)1 - (r * r)));
    const float inv_sqrt_rm0 = (float)1 / sqrtf(Rm0);

    for (int l = 1; l <= L; l++) {
        if (cur_mp->Ml[l] != 0.0f) {
            const float cos_w0l = cos_tab[l];
            const float spread = fmaxf(((float)1 + (r * r)) - ((float)2 * r * cos_w0l), 0.0f);
            Wl[l] = sqrtf(cur_mp->Ml[l] * inv_sqrt_rm0) * sqrtf(sqrtf(shape * spread));

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
}

static float
mbe_sum_spectral_magnitudes_squared(const mbe_parms* cur_mp) {
    const int L = cur_mp->L;
    float sum = 0.0f;
    int l;

#if defined(MBELIB_ENABLE_SIMD)
#if defined(MBE_SIMD_TARGET_SSE2)
    if (L >= 4) {
        __m128 vsum = _mm_setzero_ps();
        __m128 sign_mask = _mm_set1_ps(-0.0f);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            __m128 vMl = _mm_loadu_ps(&cur_mp->Ml[l4]);
            vMl = _mm_andnot_ps(sign_mask, vMl);
            vsum = _mm_add_ps(vsum, _mm_mul_ps(vMl, vMl));
        }
        __m128 shuf = _mm_shuffle_ps(vsum, vsum, _MM_SHUFFLE(2, 3, 0, 1));
        __m128 sums = _mm_add_ps(vsum, shuf);
        shuf = _mm_movehl_ps(shuf, sums);
        sums = _mm_add_ss(sums, shuf);
        sum = _mm_cvtss_f32(sums);

        for (l = l4; l <= L; l++) {
            float M = cur_mp->Ml[l];
            if (M < 0.0f) {
                M = -M;
            }
            sum += M * M;
        }
    } else
#elif defined(MBE_SIMD_TARGET_NEON)
    if (L >= 4) {
        float32x4_t vsum = vdupq_n_f32(0.0f);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            float32x4_t vMl = vld1q_f32(&cur_mp->Ml[l4]);
            vMl = vabsq_f32(vMl);
            vsum = vmlaq_f32(vsum, vMl, vMl);
        }
        float32x2_t sum2 = vadd_f32(vget_low_f32(vsum), vget_high_f32(vsum));
        sum = vget_lane_f32(vpadd_f32(sum2, sum2), 0);

        for (l = l4; l <= L; l++) {
            float M = cur_mp->Ml[l];
            if (M < 0.0f) {
                M = -M;
            }
            sum += M * M;
        }
    } else
#endif
#endif
    {
        for (l = 1; l <= L; l++) {
            float M = cur_mp->Ml[l];
            if (M < 0.0f) {
                M = -M;
            }
            sum += M * M;
        }
    }

    return sum;
}

static void
mbe_scale_spectral_magnitudes(mbe_parms* cur_mp, float gamma) {
    const int L = cur_mp->L;
    int l;

#if defined(MBELIB_ENABLE_SIMD)
#if defined(MBE_SIMD_TARGET_SSE2)
    if (L >= 4) {
        __m128 vgamma = _mm_set1_ps(gamma);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            __m128 vMl = _mm_loadu_ps(&cur_mp->Ml[l4]);
            vMl = _mm_mul_ps(vMl, vgamma);
            _mm_storeu_ps(&cur_mp->Ml[l4], vMl);
        }
        for (l = l4; l <= L; l++) {
            cur_mp->Ml[l] = gamma * cur_mp->Ml[l];
        }
    } else
#elif defined(MBE_SIMD_TARGET_NEON)
    if (L >= 4) {
        float32x4_t vgamma = vdupq_n_f32(gamma);
        int l4;
        for (l4 = 1; l4 + 3 <= L; l4 += 4) {
            float32x4_t vMl = vld1q_f32(&cur_mp->Ml[l4]);
            vMl = vmulq_f32(vMl, vgamma);
            vst1q_f32(&cur_mp->Ml[l4], vMl);
        }
        for (l = l4; l <= L; l++) {
            cur_mp->Ml[l] = gamma * cur_mp->Ml[l];
        }
    } else
#endif
#endif
    {
        for (l = 1; l <= L; l++) {
            cur_mp->Ml[l] = gamma * cur_mp->Ml[l];
        }
    }
}

/*
 * High-band compensation, applied after the enhancement's energy
 * renormalization. On DVSI's AMBE-3000 test vectors, with identical bits, the
 * TIA-102.BABA synthesis renders 2.5-3.75 kHz below DVSI's decoder relative to
 * 0-1 kHz, for voiced and unvoiced bands alike and at every pitch. A shelf
 * (0 dB below 2.5 kHz, rising to its height at 3 kHz, flat to 3.6 kHz, back to
 * 0 dB at 3.8 kHz) closes that gap; its height per codec (mbe_adaptive.h) was
 * chosen on one set of DVSI vectors and checked on another. Repeats that
 * resynthesize already enhanced amplitudes do not pass through here, so it is
 * never applied twice.
 */
static float
mbe_high_band_gain(float hz, float gain) {
    float t;
    if (hz <= 2500.0f || hz >= 3800.0f) {
        return 1.0f;
    }
    if (hz < 3000.0f) {
        t = (hz - 2500.0f) / 500.0f;
    } else if (hz <= 3600.0f) {
        t = 1.0f;
    } else {
        t = (3800.0f - hz) / 200.0f;
    }
    return 1.0f + ((gain - 1.0f) * t);
}

static void
mbe_apply_high_band_gain(mbe_parms* cur_mp, float gain) {
    const float hz_per_harmonic = cur_mp->w0 * (8000.0f / (2.0f * (float)M_PI));
    for (int l = 1; l <= cur_mp->L; l++) {
        cur_mp->Ml[l] *= mbe_high_band_gain(hz_per_harmonic * (float)l, gain);
    }
}

/*
 * Low-band roll-off, rendered by synthesis. On the same test vectors DVSI's
 * decoders roll off everything below 200 Hz like a second-order high-pass at
 * 100 Hz (-8 dB at 62 Hz, -2 dB at 112 Hz), identically in every mode, where
 * the TIA-102.BABA synthesis passes it flat. DVSI's harmonics also carry that
 * filter's phase: with it, the phase DVSI renders on held-out speech differs
 * from ours by a median 6.6-8.9 degrees (4-5 below 300 Hz) where it differed by
 * 10.0-12.1 (12.1-14.3), per harmonic of steady voiced frames once each
 * frame's time shift and common phase are removed. It is a property of the
 * output, not of the spectral envelope, so it changes only the harmonics being
 * rendered: the model, its regenerated phase and the stored state keep the
 * envelope and phases as decoded.
 */
#define MBE_LOW_BAND_CORNER_HZ 100.0f

float
mbe_low_band_gain(float hz) {
    /* 1 / sqrt(1 + (corner / hz)^4), written so that hz near or at zero gives
     * zero instead of overflowing. */
    if (!(hz > 0.0f)) {
        return 0.0f;
    }
    const float hz2 = hz * hz;
    const float corner2 = MBE_LOW_BAND_CORNER_HZ * MBE_LOW_BAND_CORNER_HZ;
    return hz2 / sqrtf((hz2 * hz2) + (corner2 * corner2));
}

float
mbe_low_band_phase(float hz) {
    /* arg of H(s) = s^2 / (s^2 + sqrt2 s + 1) at s = j hz / corner, in double:
     * atan2f's last bit differs between C libraries (glibc 2.39 and 2.44 render
     * different PCM from it), while double atan2 rounded to float agrees. */
    if (!(hz > 0.0f)) {
        return 0.0f;
    }
    const double x = (double)hz / (double)MBE_LOW_BAND_CORNER_HZ;
    return (float)(M_PI - atan2(M_SQRT2 * x, 1.0 - (x * x)));
}

/*
 * DVSI's fixed high-pass (mbe_low_band.h): a third-order Butterworth at
 * MBE_LOW_BAND_DVSI_CORNER_HZ, H(s) = s^3 / ((s + 1)(s^2 + s + 1)) at
 * s = j hz / corner, which the D-STAR, P25 and rate-33 renderings render in
 * place of the second-order filter above. In double with the double-precision
 * functions and rounded to float once, for the reason mbe_low_band_phase()
 * gives. Total over every float, by the class of its bit pattern
 * (dvsi_filter_class(): under fast-math a float test of a NaN or an infinity
 * folds away): zero, a subnormal (zero under flush-to-zero and
 * denormals-are-zero, which audio hosts set, as the synthesis treats a
 * subnormal pitch), a negative value and a NaN give 0 for both (the
 * second-order pair gives 0 for both at a NaN too, and is left as it is);
 * +infinity gives the limits, a gain of 1 and a phase of 0; every positive
 * normal value the formulas below. Near zero the gain is x^3 over a
 * denominator of about 1 (never a division by zero): at FLT_MIN x^3 is about
 * 3.85e-120 in double, which converts to 0 in float; the phase tends to
 * 3 pi / 2. Bin 0 of the transform's table (mbe_unvoiced_fft.c) calls both at
 * zero and takes 0 from them, as from the functions before.
 */
enum dvsi_filter_class { DVSI_FILTER_NONE, DVSI_FILTER_FINITE, DVSI_FILTER_INFINITE };

static enum dvsi_filter_class
dvsi_filter_class(float hz) {
    uint32_t bits;
    memcpy(&bits, &hz, sizeof(bits));
    if (bits == UINT32_C(0x7F800000)) {
        return DVSI_FILTER_INFINITE;
    }
    /* Under FLT_MIN (+0 and the subnormals), above +infinity (the positive
     * NaNs) or with the sign bit set (every negative value, -0 and the
     * negative NaNs among them). */
    if (bits < UINT32_C(0x00800000) || bits > UINT32_C(0x7F800000)) {
        return DVSI_FILTER_NONE;
    }
    return DVSI_FILTER_FINITE;
}

float
mbe_low_band_dvsi_gain(float hz) {
    const enum dvsi_filter_class kind = dvsi_filter_class(hz);
    if (kind != DVSI_FILTER_FINITE) {
        return (kind == DVSI_FILTER_INFINITE) ? 1.0f : 0.0f;
    }
    /* x is at most FLT_MAX / 75, so x^6 (under 1e220) stays finite in double. */
    const double x = (double)hz / MBE_LOW_BAND_DVSI_CORNER_HZ;
    const double cube = x * x * x;
    return (float)(cube / sqrt(1.0 + (cube * cube)));
}

/* Above this x the phase is formed so that it does not cancel. */
#define DVSI_PHASE_DIRECT_MAX_X 1024.0

float
mbe_low_band_dvsi_phase(float hz) {
    if (dvsi_filter_class(hz) != DVSI_FILTER_FINITE) {
        return 0.0f;
    }
    const double x = (double)hz / MBE_LOW_BAND_DVSI_CORNER_HZ;
    if (x <= DVSI_PHASE_DIRECT_MAX_X) {
        /* arg s^3 = 3 pi / 2, arg (s + 1) = atan x, arg (s^2 + s + 1) = atan2(x, 1 - x^2).
         * The phase falls with x, toward its high-frequency asymptote 2 / x:
         * here it is at least its value at x 1024, about 0.00195 rad, and the
         * difference's rounding error a few ulps of 3 pi / 2 in double (under
         * 1e-15 rad), far under a float ulp of the result. Every frequency the
         * renderings give (under 4 kHz, x under 54) is formed so. */
        return (float)((1.5 * M_PI) - atan(x) - atan2(x, 1.0 - (x * x)));
    }
    /* The same angle with nothing cancelled: pi / 2 - atan x = atan(1 / x) and
     * pi - atan2(x, 1 - x^2) = atan2(x, x^2 - 1) for x > 1; it tends to 2 / x,
     * so 1.5e-28 rad at 1e30 Hz rather than 0 from the difference above. */
    return (float)(atan(1.0 / x) + atan2(x, (x * x) - 1.0));
}

float
mbe_low_cut_gain(float hz, float corner) {
    if (!(hz < corner)) {
        return 1.0f;
    }
    const double ratio = (double)hz / (double)corner;
    return (float)(ratio * ratio);
}

/* A model's decoded amplitudes and phases, kept while it is rendered. */
struct mbe_render_saved {
    float Ml[57];
    float PHIl[57];
};

/* Give Ml[1..56] and PHIl[1..56] what the output renders: the low-band
 * roll-off's magnitude and phase and the voiced or unvoiced gain. */
static void
mbe_render_prepare(mbe_parms* restrict mp, struct mbe_render_saved* restrict saved, int maxl, float voiced_gain) {
    const float hz_per_harmonic = mp->w0 * (8000.0f / (2.0f * (float)M_PI));
    for (int l = 1; l <= 56; l++) {
        const float hz = hz_per_harmonic * (float)l;
        const int voiced = (l <= mp->L) && (mp->Vl[l] == 1);
        saved->Ml[l] = mp->Ml[l];
        saved->PHIl[l] = mp->PHIl[l];
        mp->Ml[l] *= mbe_low_band_gain(hz) * (voiced ? voiced_gain : MBE_RENDER_UNVOICED_GAIN);
#if !defined(MBE_VF_WIDTH)
        /* Preserve the scalar expression's rounding, including x87 excess
         * precision; only skip phases that no oscillator will read. */
        if (l <= maxl && mp->Vl[l] == 1) {
            mp->PHIl[l] += mbe_low_band_phase(hz);
        }
#endif
    }
#if defined(MBE_VF_WIDTH)
    /* Only voiced oscillators read PHIl; the unvoiced FFT uses Ml. Keep
     * the transcendental calls separate so the amplitude loop can vectorize.
     * Include harmonics fading between different model lengths. */
    for (int l = 1; l <= maxl; l++) {
        if (mp->Vl[l] == 1) {
            mp->PHIl[l] += mbe_low_band_phase(hz_per_harmonic * (float)l);
        }
    }
#endif
}

static void
mbe_render_restore(mbe_parms* mp, const struct mbe_render_saved* saved) {
    for (int l = 1; l <= 56; l++) {
        mp->Ml[l] = saved->Ml[l];
        mp->PHIl[l] = saved->PHIl[l];
    }
}

/*
 * The enhancement's energy sums square the amplitudes: in float they overflow
 * for a peak amplitude above about 2^63 and lose the envelope below 2^-63
 * (decoded models stay within 2^-30..2^49). A model with a peak beyond 2^60 or
 * 2^-60 is enhanced in double precision instead, as is one whose weights float
 * cannot resolve (mbe_spectral_float_unresolved()).
 */
static int
mbe_spectral_needs_wide_range(const mbe_parms* cur_mp) {
    float peak = 0.0f;
    for (int l = 1; l <= cur_mp->L; l++) {
        const float magnitude = fabsf(cur_mp->Ml[l]);
        if (magnitude > peak) {
            peak = magnitude;
        }
    }
    return (peak > 0x1p60f) || (peak > 0.0f && peak < 0x1p-60f);
}

/*
 * Whether float cannot resolve mbe_apply_spectral_weights(): its denominator,
 * w0 (1 - r^2), or the numerator 1 + r^2 - 2 r cos(w0 l) of a weighted
 * harmonic comes too close to 0, where the float cosines' error (about 1e-6)
 * and rounding dominate. Decoded models keep 1 - r^2 above 2.5e-3 and every
 * numerator above 1.5e-3 (energy in harmonic 1 at the lowest pitch).
 */
static int
mbe_spectral_float_unresolved(const mbe_parms* cur_mp, const float cos_tab[57], float Rm0, float Rm1) {
    if (!(Rm0 > 0.0f) || !(cur_mp->w0 > 0.0f)) {
        return 0;
    }
    const float r = Rm1 / Rm0;
    if (!(((float)1 - (r * r)) > 0x1p-14f)) {
        return 1;
    }
    for (int l = 1; l <= cur_mp->L; l++) {
        const float spread = ((float)1 + (r * r)) - ((float)2 * r * cos_tab[l]);
        if ((8 * l) > cur_mp->L && cur_mp->Ml[l] != 0.0f && !(spread > 0x1p-12f)) {
            return 1;
        }
    }
    return 0;
}

/* One weighted amplitude of mbe_apply_spectral_weights(), in double precision.
 * As |r| reaches 1 (all the energy at cos(w0 l) = r = +-1) the numerator and
 * the denominator both vanish and W_l tends to
 * sqrt(M_l / sqrt(Rm0)) (0.96 pi / w0)^(1/4), which is used where even double
 * rounding leaves no positive denominator. */
static double
mbe_wide_weighted_amplitude(const mbe_parms* cur_mp, int l, double cos_w0l, double rm0, double r, int singular) {
    const double w0 = (double)cur_mp->w0;
    const double ml = (double)cur_mp->Ml[l];
    const double spread = singular ? 1.0 : fmax((1.0 + (r * r)) - (2.0 * r * cos_w0l), 0.0);
    const double shape = (0.96 * M_PI) / (singular ? w0 : w0 * (1.0 - (r * r)));
    const double weight = sqrt(ml / sqrt(rm0)) * sqrt(sqrt(shape * spread));
    if ((8 * l) <= cur_mp->L) {
        return ml;
    }
    return ml * ((weight > 1.2) ? 1.2 : ((weight < 0.5) ? 0.5 : weight));
}

/* mbe_spectralAmpEnhanceWithRm0() in double precision, its cosines too (a
 * float cosine of -1 can round past it), amplitudes held to the float range;
 * returns Rm0, held to it too. */
static float
mbe_spectral_enhance_wide_range(mbe_parms* cur_mp, float high_band_gain) {
    const int L = cur_mp->L;
    double rm0 = 0.0, rm1 = 0.0, sum = 0.0;
    double enhanced[57], cosines[57];
    for (int l = 1; l <= L; l++) {
        const double power = (double)cur_mp->Ml[l] * (double)cur_mp->Ml[l];
        cosines[l] = cos((double)cur_mp->w0 * (double)l);
        rm0 += power;
        rm1 += power * cosines[l];
    }
    const int weighted = (rm0 > 0.0) && (cur_mp->w0 > 0.0f);
    const double r = weighted ? rm1 / rm0 : 0.0;
    const int singular = weighted && !((double)cur_mp->w0 * (1.0 - (r * r)) > 0.0);
    for (int l = 1; l <= L; l++) {
        enhanced[l] = (weighted && cur_mp->Ml[l] != 0.0f)
                          ? mbe_wide_weighted_amplitude(cur_mp, l, cosines[l], rm0, r, singular)
                          : (double)cur_mp->Ml[l];
        sum += enhanced[l] * enhanced[l];
    }
    const double gamma = (sum > 0.0) ? sqrt(rm0 / sum) : 1.0;
    const float hz_per_harmonic = cur_mp->w0 * (8000.0f / (2.0f * (float)M_PI));
    for (int l = 1; l <= L; l++) {
        const double value =
            enhanced[l] * gamma * (double)mbe_high_band_gain(hz_per_harmonic * (float)l, high_band_gain);
        cur_mp->Ml[l] = (value > (double)FLT_MAX) ? FLT_MAX : (float)value;
    }
    return (rm0 > (double)FLT_MAX) ? FLT_MAX : (float)rm0;
}

/**
 * @brief Apply spectral amplitude enhancement to the current parameters.
 * @param cur_mp In/out parameter set to enhance.
 *
 * Uses SIMD optimizations for accumulation and scaling loops when available.
 */
float
mbe_spectralAmpEnhanceWithRm0(mbe_parms* cur_mp, float high_band_gain) {

    float Rm0, Rm1;
    float cos_tab[57];

    if (MBE_UNLIKELY(!cur_mp || !mbe_harmonic_count_is_valid(cur_mp->L))) {
        return 0.0f;
    }

    mbe_precompute_harmonic_cosines(cos_tab, cur_mp->L, cur_mp->w0);
    if (MBE_UNLIKELY(mbe_spectral_needs_wide_range(cur_mp))) {
        return mbe_spectral_enhance_wide_range(cur_mp, high_band_gain);
    }
    mbe_accumulate_spectral_energy(cur_mp, cos_tab, &Rm0, &Rm1);
    if (MBE_UNLIKELY(mbe_spectral_float_unresolved(cur_mp, cos_tab, Rm0, Rm1))) {
        return mbe_spectral_enhance_wide_range(cur_mp, high_band_gain);
    }

    mbe_apply_spectral_weights(cur_mp, cos_tab, Rm0, Rm1);

    float sum = mbe_sum_spectral_magnitudes_squared(cur_mp);
    float gamma = (sum == 0.0f) ? 1.0f : sqrtf(Rm0 / sum);
    mbe_scale_spectral_magnitudes(cur_mp, gamma);
    mbe_apply_high_band_gain(cur_mp, high_band_gain);

    return Rm0;
}

void
mbe_spectralAmpEnhance(mbe_parms* cur_mp) {
    (void)mbe_spectralAmpEnhanceWithRm0(cur_mp, MBE_HIGH_BAND_GAIN_IMBE);
}

/* Library output headroom and phase-accumulator policies; independently implemented. */
#define MBE_FLOAT_LIMIT ((32767.0f * 0.95f) / 7.0f)

static float
mbe_clipFloatSample(float sample) {
    uint32_t bits;
    uint32_t limit_bits;
    const float limit = MBE_FLOAT_LIMIT;
    memcpy(&bits, &sample, sizeof(bits));
    memcpy(&limit_bits, &limit, sizeof(limit_bits));
    /* Volatile integer classification remains observable under finite-math assumptions. */
    const volatile uint32_t magnitude = bits & UINT32_C(0x7fffffff);
    if (magnitude > limit_bits && magnitude <= UINT32_C(0x7f800000)) {
        return (bits & UINT32_C(0x80000000)) != 0 ? -limit : limit;
    }
    return sample;
}

void
mbe_clipFloatBuffer(float* samples, int count) {
    if (samples == NULL || count <= 0) {
        return;
    }
    for (int i = 0; i < count; ++i) {
        samples[i] = mbe_clipFloatSample(samples[i]);
    }
}

#ifndef DISABLE_AMBE_TONES
static int
finite_gain(float value) {
    uint32_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT32_C(0x7fffffff)) <= UINT32_C(0x7f7fffff);
}

/* A component of m steps of 2^-19 cycle per sample (mbe_tone_played_steps())
 * as the 32-bit phase step of the tone oscillators: exact. */
static uint32_t
tone_phase_step(unsigned m) {
    return (uint32_t)(m << (32 - MBE_TONE_GRID_BITS));
}

static float
phase_sine(uint32_t phase) {
    const float angle = (float)((double)phase * (2.0 * M_PI / 4294967296.0) - M_PI / 2.0);
    return sinf(angle);
}

static float
mbe_toneAmplitude(double level_db) {
    return (float)(sqrt(2.0) * 32768.0 * pow(10.0, level_db / 20.0) / 7.0);
}

/* A tone whose gain moves linearly from from_gain to gain over the first
 * `ramp` samples (none for 0), each sample limited after the gain is applied,
 * unless the caller fades or mixes the tone and limits the result. The two
 * oscillators advance by step1 and step2 a sample (one component when the two
 * are equal or step2 is 0) from the phases mp carries. The sample loop is the
 * same in both modes: it limits every sample as it computes it (the limiter's
 * classification keeps the loop scalar, so fast-math cannot vectorize it with
 * other sine code) and keeps the unlimited sample beside it, and only after
 * the loop is the unlimited frame chosen. So the limited rendering is computed
 * exactly as before, and limiting the unlimited one gives it bit for bit. */
static void
mbe_renderToneRampf(float* out, mbe_parms* mp, uint32_t step1, uint32_t step2, float from_gain, float gain, int ramp,
                    enum mbe_output_limit output) {
    if (out == NULL) {
        return;
    }
    if (mp == NULL || !finite_gain(gain) || !finite_gain(from_gain) || step1 == 0u) {
        memset(out, 0, 160 * sizeof(*out));
        return;
    }
    const int dual = step2 != 0u && step2 != step1;
    uint32_t phase1;
    memcpy(&phase1, &mp->swn, sizeof(phase1));
    uint32_t phase2 = mp->tonePhase;
    float unlimited[160];
    for (int n = 0; n < 160; ++n) {
        const float g = (n < ramp) ? from_gain + ((gain - from_gain) * ((float)n / (float)ramp)) : gain;
        phase1 += step1;
        float value = g * phase_sine(phase1);
        if (dual) {
            phase2 += step2;
            value += g * phase_sine(phase2);
        }
        out[n] = mbe_clipFloatSample(value);
        unlimited[n] = value;
    }
    if (output == MBE_OUTPUT_UNLIMITED) {
        memcpy(out, unlimited, sizeof(unlimited));
    }
    memcpy(&mp->swn, &phase1, sizeof(phase1));
    mp->tonePhase = phase2;
}

static void
mbe_renderTonef(float* out, mbe_parms* mp, uint32_t step1, uint32_t step2, float gain, enum mbe_output_limit output) {
    mbe_renderToneRampf(out, mp, step1, step2, gain, gain, 0, output);
}

/*
 * The oscillator phase that starts a component of m steps as DVSI's decoders
 * start a tone after a frame with no voiced band: as a fresh pulse train from
 * the frame's start, so that sample n (0..159) of the frame is
 * cos(2 pi F (n - MBE_PULSE_PEAK) + p), F the frequency in cycles per sample
 * and p the phase of DVSI's fixed high-pass at it (mbe_low_band_dvsi_phase(),
 * the third-order filter at 75 Hz of D-STAR's and rate 33's renderings; before
 * the low-band change the second-order one at 100 Hz). Measured on DVSI's
 * rate-33 and D-STAR tone vectors: within one vector the phase at such
 * onsets of one component varies by 0.00013 cycle at the most (between vectors
 * by up to 0.00158), and the measured phase less this one is +0.0003 cycle on
 * average over the 64 DTMF and KNOX components (766 starts), +0.0016 with the
 * second-order filter (docs/dvsi-parity.md, "DVSI's low band"). After a voiced frame DVSI's decoders continue
 * that frame's pulse train into the tone instead (cp2, cp31); the codec paths
 * start every tone after a frame that was not a tone frame this way, the
 * voiced case included, which is a known difference (docs/dvsi-parity.md,
 * "Tones as harmonics of a fundamental"). The renderer plays phase_sine() of the phase
 * after adding the step, -cos of it: the phase before sample 0 is
 * step (-1 - MBE_PULSE_PEAK) + 2^32 (1/2 + p / (2 pi)).
 */
static uint32_t
tone_start_phase(unsigned m) {
    const uint32_t step = tone_phase_step(m);
    const float hz = (float)((double)m * (8000.0 / (double)(1u << MBE_TONE_GRID_BITS)));
    const double turn = 0.5 + ((double)mbe_low_band_dvsi_phase(hz) / (2.0 * M_PI));
    const uint32_t offset = (uint32_t)(uint64_t)floor(((turn - floor(turn)) * 4294967296.0) + 0.5);
    return offset - (step * (uint32_t)(1 + MBE_PULSE_PEAK));
}

static void
mbe_startTone(mbe_parms* mp, unsigned first, unsigned second) {
    const uint32_t phase1 = tone_start_phase(first);
    memcpy(&mp->swn, &phase1, sizeof(phase1));
    mp->tonePhase = (second != 0u && second != first) ? tone_start_phase(second) : 0u;
}

/* The two oscillator steps of an AMBE+2 tone ID (the upper component first, as
 * Table 9's lookup gives the nominal pair); 0 if it has none. */
static int
mbe_tone_steps_ambe2450(int id, uint32_t* step1, uint32_t* step2) {
    unsigned upper, lower;
    const int found = mbe_tone_played_steps(id, &upper, &lower);
    *step1 = tone_phase_step(upper);
    *step2 = tone_phase_step(lower);
    return found;
}

/* The two oscillator steps of a D-STAR tone index, in its lookup's order. */
static int
mbe_tone_steps_dstar(int index, uint32_t* step1, uint32_t* step2) {
    unsigned first, second;
    const int found = mbe_tone_dstar_played_steps(index, &first, &second);
    *step1 = tone_phase_step(first);
    *step2 = tone_phase_step(second);
    return found;
}
#endif

/**
 * @brief Synthesize a tone frame into 160 float samples at 8 kHz.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49) providing tone indices.
 * @param cur_mp   Current parameter set (tone synthesis state).
 */
void
mbe_synthesizeTonef(float* aout_buf, const char* ambe_d, mbe_parms* cur_mp) {
#ifdef DISABLE_AMBE_TONES // generate silence if tones disabled
    (void)ambe_d;
    (void)cur_mp;
    mbe_synthesizeSilencef(aout_buf);
    return;
#else
    int i;

    if (!cur_mp || mbe_validate_bits(ambe_d, 49u) < 0) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }

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
    ID3 = ((u2 & 0x7f) << 1) + ((u3 >> 13) & 0x1);
    ID4 = ((u3 & 0x1fe0) >> 5);

    uint32_t step1, step2;
    (void)ID0;
    (void)ID2;
    (void)ID3; /* parsed for potential validation */
    (void)ID4;

    /* The frequencies DVSI's decoders play for the ID (mbe_tone_played_steps()),
     * continuing the phases cur_mp carries: this call has no frame history. */
    if (!mbe_tone_steps_ambe2450(ID1, &step1, &step2)) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }

    mbe_renderTonef(aout_buf, cur_mp, step1, step2, mbe_toneAmplitude(mbe_tone_ambe2450_level_db(AD)),
                    MBE_OUTPUT_LIMITED);
#endif
}

void
mbe_renderToneAmbe2450f(float* aout_buf, mbe_parms* mp, int id, int ad, enum mbe_output_limit output) {
#ifdef DISABLE_AMBE_TONES
    (void)mp;
    (void)id;
    (void)ad;
    (void)output;
    mbe_synthesizeSilencef(aout_buf);
#else
    uint32_t step1, step2;
    if (!mbe_tone_steps_ambe2450(id, &step1, &step2)) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }
    mbe_renderTonef(aout_buf, mp, step1, step2, mbe_toneAmplitude(mbe_tone_ambe2450_level_db(ad)), output);
#endif
}

void
mbe_startToneAmbe2450(mbe_parms* mp, int id) {
#ifdef DISABLE_AMBE_TONES
    (void)mp;
    (void)id;
#else
    unsigned upper, lower;
    if (mp != NULL && mbe_tone_played_steps(id, &upper, &lower)) {
        mbe_startTone(mp, upper, lower);
    }
#endif
}

/**
 * @brief Synthesize a D-STAR style tone into 160 float samples.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49).
 * @param cur_mp   Current parameter set.
 * @param ID1      Tone index selector.
 */
void
mbe_synthesizeTonefdstar(float* aout_buf, const char* ambe_d, mbe_parms* cur_mp, int ID1) {
#ifdef DISABLE_AMBE_TONES // generate silence if tones disabled
    (void)ambe_d;
    (void)cur_mp;
    (void)ID1;
    mbe_synthesizeSilencef(aout_buf);
    return;
#else
    if (!cur_mp || mbe_validate_bits(ambe_d, 49u) < 0) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }
    mbe_renderToneDstarf(aout_buf, cur_mp, ID1, mbe_tone_dstar_volume(ambe_d));
#endif
}

void
mbe_renderToneDstarf(float* aout_buf, mbe_parms* mp, int index, int volume) {
#ifdef DISABLE_AMBE_TONES
    (void)mp;
    (void)index;
    (void)volume;
    mbe_synthesizeSilencef(aout_buf);
#else
    uint32_t step1, step2;

    if (!mp || !mbe_tone_steps_dstar(index, &step1, &step2)) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }
    mbe_renderTonef(aout_buf, mp, step1, step2, mbe_toneAmplitude(mbe_tone_dstar_level_db(volume)), MBE_OUTPUT_LIMITED);
#endif
}

void
mbe_renderToneDstarFromf(float* aout_buf, mbe_parms* mp, int index, int from_volume, int volume, int ramp) {
#ifdef DISABLE_AMBE_TONES
    (void)mp;
    (void)index;
    (void)from_volume;
    (void)volume;
    (void)ramp;
    mbe_synthesizeSilencef(aout_buf);
#else
    uint32_t step1, step2;

    if (!mp || !mbe_tone_steps_dstar(index, &step1, &step2)) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }
    mbe_renderToneRampf(aout_buf, mp, step1, step2, mbe_toneAmplitude(mbe_tone_dstar_level_db(from_volume)),
                        mbe_toneAmplitude(mbe_tone_dstar_level_db(volume)), (ramp > 0 && ramp <= 160) ? ramp : 0,
                        MBE_OUTPUT_LIMITED);
#endif
}

void
mbe_startToneDstar(mbe_parms* mp, int index) {
#ifdef DISABLE_AMBE_TONES
    (void)mp;
    (void)index;
#else
    unsigned first, second;
    if (mp != NULL && mbe_tone_dstar_played_steps(index, &first, &second)) {
        mbe_startTone(mp, first, second);
    }
#endif
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

/* Regenerated voiced phase (US 5,701,390 Eqs. 7-9): odd 1/m edge-detection
 * kernel over log2 magnitudes. The TIA-102.BABA rendering applies it to the
 * enhanced, smoothed magnitudes it renders (mbe_regenerate_phases()); the
 * codec paths, as DVSI's decoders do, to the frame's decoded log2 amplitudes
 * (mbe_regenerate_phases_log2(); docs/dvsi-parity.md, "Phase regeneration"). */
#define MBE_PHASE_KERNEL_HALF_LEN 19
#define MBE_PHASE_KERNEL_SCALE    0.44f
#define MBE_PHASE_EXTRAP_SLOPE    0.72f
#define MBE_PHASE_MIN_MAGNITUDE   1e-6f

/* IMBE state (any non-AMBE muting threshold) mutes on the error rate. */
static int
mbe_uses_imbe_muting(const mbe_parms* cur_mp) {
    return fabsf(cur_mp->mutingThreshold - MBE_MUTING_THRESHOLD_AMBE) > 1e-6f;
}

/*
 * IMBE mutes here on max repeats or the error rate; D-STAR mutes on max
 * repeats. The caller's codec path picks the noise. The AMBE 3600x2450 process
 * path decides repeats and muting itself (TIA-102.BABA-1 5.6/5.7) and only
 * calls synthesis with repeatCount below MBE_MAX_FRAME_REPEATS, so this check
 * never fires for it.
 */
static int
mbe_should_mute_speech(const mbe_parms* cur_mp) {
    return mbe_isMaxFrameRepeat(cur_mp) || (mbe_uses_imbe_muting(cur_mp) && mbe_requiresMuting(cur_mp));
}

void
mbe_regenerate_phases(const float Ml[57], int L, float phi[57]) {
    const int D = MBE_PHASE_KERNEL_HALF_LEN;
    static const float inv_odd[10] = {
        1.0f, 1.0f / 3, 1.0f / 5, 1.0f / 7, 1.0f / 9, 1.0f / 11, 1.0f / 13, 1.0f / 15, 1.0f / 17, 1.0f / 19,
    };
    float B[56 + 2 * MBE_PHASE_KERNEL_HALF_LEN + 1];
    B[D] = 0.0f;
    for (int l = 1; l <= L; ++l) {
        B[D + l] = log2f(fmaxf(Ml[l], MBE_PHASE_MIN_MAGNITUDE));
    }
    /* Define the tail for all 56 phases, including harmonics fading out. */
    for (int l = L + 1; l <= 56 + D; ++l) {
        B[D + l] = B[D + L] - MBE_PHASE_EXTRAP_SLOPE * (float)(l - L);
    }
    for (int l = 1; l <= D; ++l) {
        B[D - l] = B[D + l];
    }
#if defined(MBE_VF_WIDTH)
    /* Independent harmonics in each lane: retain each phase's tap order
     * instead of horizontally reducing (and reassociating) its ten terms. */
    for (int l = 1; l <= 56; l += MBE_VF_WIDTH) {
        mbe_vf acc = mbe_vf_set(0.0f);
        for (int k = 0; k < 10; ++k) {
            int m = 2 * k + 1;
            mbe_vf delta = mbe_vf_sub(mbe_vf_load(B + D + l + m), mbe_vf_load(B + D + l - m));
            acc = mbe_vf_add(acc, mbe_vf_mul(delta, mbe_vf_set(inv_odd[k])));
        }
        mbe_vf_store(phi + l, mbe_vf_mul(mbe_vf_set(MBE_PHASE_KERNEL_SCALE), acc));
    }
#else
    for (int l = 1; l <= 56; ++l) {
        float acc = 0.0f;
        for (int k = 0; k < 10; ++k) {
            int m = 2 * k + 1;
            acc += (B[D + l + m] - B[D + l - m]) * inv_odd[k];
        }
        phi[l] = MBE_PHASE_KERNEL_SCALE * acc;
    }
#endif
}

/* The kernel of mbe_regenerate_phases() over log2 amplitudes given directly,
 * with its floor (MBE_PHASE_LOG2_FLOOR is log2(MBE_PHASE_MIN_MAGNITUDE)), its
 * tail and its mirror. A loop of its own, so that the TIA-102.BABA rendering's
 * code stays as it is; it adds each phase's taps in the same order. */
void
mbe_regenerate_phases_log2(const float log2_amplitude[57], int L, float phi[57]) {
    const int D = MBE_PHASE_KERNEL_HALF_LEN;
    static const float inv_odd[10] = {
        1.0f, 1.0f / 3, 1.0f / 5, 1.0f / 7, 1.0f / 9, 1.0f / 11, 1.0f / 13, 1.0f / 15, 1.0f / 17, 1.0f / 19,
    };
    float B[56 + 2 * MBE_PHASE_KERNEL_HALF_LEN + 1];
    B[D] = 0.0f;
    for (int l = 1; l <= L; ++l) {
        B[D + l] = fmaxf(log2_amplitude[l], MBE_PHASE_LOG2_FLOOR);
    }
    for (int l = L + 1; l <= 56 + D; ++l) {
        B[D + l] = B[D + L] - MBE_PHASE_EXTRAP_SLOPE * (float)(l - L);
    }
    for (int l = 1; l <= D; ++l) {
        B[D - l] = B[D + l];
    }
    for (int l = 1; l <= 56; ++l) {
        float acc = 0.0f;
        for (int k = 0; k < 10; ++k) {
            int m = 2 * k + 1;
            acc += (B[D + l + m] - B[D + l - m]) * inv_odd[k];
        }
        phi[l] = MBE_PHASE_KERNEL_SCALE * acc;
    }
}

static int
mbe_reconcile_speech_model_lengths(mbe_parms* cur_mp, mbe_parms* prev_mp) {
    if (cur_mp->L > prev_mp->L) {
        int maxl = cur_mp->L;
        for (int l = prev_mp->L + 1; l <= maxl; l++) {
            prev_mp->Ml[l] = 0.0f;
            prev_mp->Vl[l] = 1;
        }
        return maxl;
    }

    int maxl = prev_mp->L;
    for (int l = cur_mp->L + 1; l <= maxl; l++) {
        cur_mp->Ml[l] = 0.0f;
        cur_mp->Vl[l] = 1;
    }
    return maxl;
}

void
mbe_set_frame_mark(mbe_parms* mp, enum mbe_frame_mark mark) {
    mp->PHIl[0] = (float)mark;
}

enum mbe_frame_mark
mbe_frame_mark(const mbe_parms* mp) {
    /* Compare the stored bits: only an exact mark counts, so any other value
     * a caller left there (NaN included, even under fast-math) means none. */
    uint32_t bits;
    memcpy(&bits, &mp->PHIl[0], sizeof(bits));
    switch (bits) {
        case UINT32_C(0x3F800000): return MBE_FRAME_MARK_PITCH_CARRIES;   /* 1.0f */
        case UINT32_C(0x40000000): return MBE_FRAME_MARK_PITCH_ABSENT;    /* 2.0f */
        case UINT32_C(0x40400000): return MBE_FRAME_MARK_PREVIOUS_HELD;   /* 3.0f */
        case UINT32_C(0x40800000): return MBE_FRAME_MARK_TONE;            /* 4.0f */
        case UINT32_C(0x40A00000): return MBE_FRAME_MARK_ONSET_COLUMNS;   /* 5.0f */
        case UINT32_C(0x40C00000): return MBE_FRAME_MARK_MID_COLUMN_0;    /* 6.0f */
        case UINT32_C(0x40E00000): return MBE_FRAME_MARK_MID_COLUMNS_0_1; /* 7.0f */
        case UINT32_C(0x42A00000): return MBE_FRAME_MARK_LATE_ONSET;      /* 80.0f */
        default: return MBE_FRAME_MARK_NONE;
    }
}

float
mbe_pulsed_row_value(const struct mbe_pulsed_row* word) {
    const int value = MBE_PULSED_ROW_WORD_BASE + ((word->held != 0) ? MBE_PULSED_ROW_WORD_HELD : 0)
                      + (MBE_PULSED_ROW_WORD_ROW * (word->row - MBE_PULSED_ROW_FIRST)) + word->b0;
    return (float)value;
}

int
mbe_pulsed_row_read(const float* stored, struct mbe_pulsed_row* word) {
    const struct mbe_pulsed_row last = {MBE_PULSED_ROW_OWN_LAST, MBE_PULSED_ROW_CODES - 1, 1};
    const float lowest = (float)MBE_PULSED_ROW_WORD_BASE;
    const float highest = mbe_pulsed_row_value(&last);
    uint32_t bits, low_bits, high_bits, rebuilt_bits;
    memcpy(&bits, stored, sizeof(bits));
    memcpy(&low_bits, &lowest, sizeof(low_bits));
    memcpy(&high_bits, &highest, sizeof(high_bits));
    /* Positive floats order as their bit patterns do: a negative value, not a
     * number or an infinity lies outside, so the conversion sees 1024..4087. */
    if (bits < low_bits || bits > high_bits) {
        return 0;
    }
    const int value = (int)*stored;
    const float rebuilt = (float)value;
    memcpy(&rebuilt_bits, &rebuilt, sizeof(rebuilt_bits));
    if (rebuilt_bits != bits) {
        return 0;
    }
    const int fields = value - MBE_PULSED_ROW_WORD_BASE;
    const int rest = fields % MBE_PULSED_ROW_WORD_HELD;
    const int row = MBE_PULSED_ROW_FIRST + (rest / MBE_PULSED_ROW_WORD_ROW);
    const int b0 = rest % MBE_PULSED_ROW_WORD_ROW;
    const int held = fields / MBE_PULSED_ROW_WORD_HELD;
    /* A held word is a repeat's, of a row whose own subframe pulses (18 to
     * 25): no writer gives rows 26 to 31 a held form (ambe2450_load_repeat()
     * clears their word), so such a value is no word. */
    if (row > MBE_PULSED_ROW_LAST || b0 >= MBE_PULSED_ROW_CODES || (held && row > MBE_PULSED_ROW_OWN_LAST)) {
        return 0;
    }
    word->row = row;
    word->b0 = b0;
    word->held = held;
    return 1;
}

void
mbe_set_pulsed_row_word(mbe_parms* mp, const struct mbe_pulsed_row* word) {
    mp->PHIl[0] = mbe_pulsed_row_value(word);
}

int
mbe_pulsed_row_word(const mbe_parms* mp, struct mbe_pulsed_row* word) {
    return mbe_pulsed_row_read(&mp->PHIl[0], word);
}

void
mbe_set_voicing_onset_delay(mbe_parms* mp, int late) {
    mbe_set_frame_mark(mp, late ? MBE_FRAME_MARK_LATE_ONSET : MBE_FRAME_MARK_NONE);
}

int
mbe_voicing_onset_delay(const mbe_parms* mp) {
    return (mbe_frame_mark(mp) == MBE_FRAME_MARK_LATE_ONSET) ? MBE_VOICING_ONSET_LATE : 0;
}

/**
 * PSIl is integrated linear phase (pulse position, TIA-102.BABA eq 139). PHIl
 * is the phase at the end of the frame: PSIl plus the regenerated
 * spectral-envelope phase.
 */
static void
mbe_update_speech_phases(mbe_parms* cur_mp, mbe_parms* prev_mp, int N) {
    const float cw0 = cur_mp->w0;
    const float pw0 = prev_mp->w0;
    float phi[57];
    mbe_regenerate_phases(cur_mp->Ml, cur_mp->L, phi); /* L validated by the synthesis core */

    for (int l = 1; l <= 56; l++) {
        float prev_psil_wrapped = fmodf(prev_mp->PSIl[l], MBE_TWO_PI);
        if (prev_psil_wrapped < 0.0f) {
            prev_psil_wrapped += MBE_TWO_PI;
        }
        prev_mp->PSIl[l] = prev_psil_wrapped;

        cur_mp->PSIl[l] = prev_psil_wrapped + ((pw0 + cw0) * ((float)(l * N) / 2.0f));
        cur_mp->PHIl[l] = cur_mp->PSIl[l] + phi[l];
    }
}

#if defined(MBE_VOICED_KERNELS)
static void
mbe_render_voiced_interpolated(float* aout_buf, const mbe_parms* cur_mp, const mbe_parms* prev_mp, int l, int N,
                               float cw0, float pw0, float pw0l) {
    float deltaphil = cur_mp->PHIl[l] - prev_mp->PHIl[l] - (((pw0 + cw0) * (float)(l * N)) / 2.0f);
    float deltawl = (1.0f / (float)N)
                    * (deltaphil - (2.0f * (float)M_PI * floorf((deltaphil + (float)M_PI) / (2.0f * (float)M_PI))));
    mbe_voiced_interpolated(aout_buf, prev_mp->PHIl[l], prev_mp->Ml[l], cur_mp->Ml[l], pw0l + deltawl, cw0 - pw0, l);
}

static void
mbe_render_voiced_windowed(float* aout_buf, const mbe_parms* cur_mp, const mbe_parms* prev_mp, int l, int N, float cw0l,
                           float pw0l, int cur_voiced, int prev_voiced) {
    const struct mbe_voiced_component prev = {prev_mp->PHIl[l], pw0l, 2.0f * prev_mp->Ml[l], prev_voiced};
    const struct mbe_voiced_component cur = {cur_mp->PHIl[l] - cw0l * (float)N, cw0l, 2.0f * cur_mp->Ml[l], cur_voiced};
    mbe_voiced_windowed(aout_buf, Ws, &prev, &cur);
}

#else
static void
mbe_render_voiced_interpolated(float* aout_buf, const mbe_parms* cur_mp, const mbe_parms* prev_mp, int l, int N,
                               float cw0, float pw0, float pw0l) {
    float* Ss = aout_buf;
    float deltaphil = cur_mp->PHIl[l] - prev_mp->PHIl[l] - (((pw0 + cw0) * (float)(l * N)) / 2.0f);
    float deltawl = (1.0f / (float)N)
                    * (deltaphil - (2.0f * (float)M_PI * floorf((deltaphil + (float)M_PI) / (2.0f * (float)M_PI))));

    for (int n = 0; n < N; n++) {
        float thetaln =
            prev_mp->PHIl[l] + ((pw0l + deltawl) * (float)n) + (((cw0 - pw0) * (float)(l * n * n)) / (float)(2 * N));
        float aln = prev_mp->Ml[l] + (((float)n / (float)N) * (cur_mp->Ml[l] - prev_mp->Ml[l]));
        *Ss += 2.0f * aln * cosf(thetaln);
        Ss++;
    }
}

static void
mbe_render_voiced_windowed(float* aout_buf, const mbe_parms* cur_mp, const mbe_parms* prev_mp, int l, int N, float cw0l,
                           float pw0l, int cur_voiced, int prev_voiced) {
    float* Ss = aout_buf;

    if (prev_voiced && cur_voiced) {
        const float gain_prev = 2.0f * prev_mp->Ml[l];
        float sd_prev, cd_prev;
        mbe_sincosf(pw0l, &sd_prev, &cd_prev);
        float s_prev, c_prev;
        mbe_sincosf(prev_mp->PHIl[l], &s_prev, &c_prev);

        const float gain_cur = 2.0f * cur_mp->Ml[l];
        float sd_cur, cd_cur;
        mbe_sincosf(cw0l, &sd_cur, &cd_cur);
        float s_cur, c_cur;
        mbe_sincosf(cur_mp->PHIl[l] - (cw0l * (float)N), &s_cur, &c_cur);

        for (int n = 0; n < N; n += 4) {
            const mbe_voiced_block_component prev_component = {Ws + n + N, gain_prev, &c_prev,
                                                               &s_prev,    sd_prev,   cd_prev};
            const mbe_voiced_block_component cur_component = {Ws + n, gain_cur, &c_cur, &s_cur, sd_cur, cd_cur};
            mbe_add_voiced_dual_block4(Ss, &prev_component, &cur_component);
            Ss += 4;
        }
    } else if (prev_voiced) {
        const float gain_prev = 2.0f * prev_mp->Ml[l];
        float sd_prev, cd_prev;
        mbe_sincosf(pw0l, &sd_prev, &cd_prev);
        float s_prev, c_prev;
        mbe_sincosf(prev_mp->PHIl[l], &s_prev, &c_prev);

        for (int n = 0; n < N; n += 4) {
            mbe_add_voiced_block4(Ss, Ws + n + N, gain_prev, &c_prev, &s_prev, sd_prev, cd_prev);
            Ss += 4;
        }
    } else if (cur_voiced) {
        const float gain_cur = 2.0f * cur_mp->Ml[l];
        float sd_cur, cd_cur;
        mbe_sincosf(cw0l, &sd_cur, &cd_cur);
        float s_cur, c_cur;
        mbe_sincosf(cur_mp->PHIl[l] - (cw0l * (float)N), &s_cur, &c_cur);

        for (int n = 0; n < N; n += 4) {
            mbe_add_voiced_block4(Ss, Ws + n, gain_cur, &c_cur, &s_cur, sd_cur, cd_cur);
            Ss += 4;
        }
    }
}
#endif

/* The voiced part as TIA-102.BABA renders it (the public synthesis; the codec
 * paths render DVSI's subframes, mbe_subframe.h). */
static void
mbe_render_voiced_speech(float* aout_buf, const mbe_parms* cur_mp, const mbe_parms* prev_mp, int maxl, int N) {
    const float cw0 = cur_mp->w0;
    const float pw0 = prev_mp->w0;

    for (int l = 1; l <= maxl; l++) {
        float cw0l = cw0 * (float)l;
        float pw0l = pw0 * (float)l;
        /* A harmonic at or above Nyquist cannot be represented at 8 kHz and
         * would alias back into the band, so it is not rendered. Decoded
         * models always satisfy L * w0 < pi; this guards caller-supplied ones. */
        int cur_voiced = (cur_mp->Vl[l] == 1) && (cw0l < (float)M_PI);
        int prev_voiced = (prev_mp->Vl[l] == 1) && (pw0l < (float)M_PI);

        if (cur_voiced || prev_voiced) {
            /* Phase/amplitude interpolation (Algorithms #134-138) for harmonics
             * 1..7 at steady pitch, the windowed oscillators otherwise. */
            int use_interpolation = (l < 8) && cur_voiced && prev_voiced && (fabsf(cw0 - pw0) < (0.1f * cw0));
            if (use_interpolation) {
                mbe_render_voiced_interpolated(aout_buf, cur_mp, prev_mp, l, N, cw0, pw0, pw0l);
            } else {
                mbe_render_voiced_windowed(aout_buf, cur_mp, prev_mp, l, N, cw0l, pw0l, cur_voiced, prev_voiced);
            }
        }
    }
}

/* How the synthesis core applies adaptive smoothing (Algorithms #111-116). */
enum mbe_speech_smoothing {
    MBE_SMOOTHING_FROM_MODEL,  /* derive the pre-enhancement energy from cur_mp */
    MBE_SMOOTHING_PRE_ENH_RM0, /* the caller captured the pre-enhancement energy */
    MBE_SMOOTHING_NONE,        /* frame repeat: synthesize the model unchanged */
};

/* A frame the core mutes: noise in place of speech (P25 IMBE uses the
 * TIA-102.BABA 7.8 level, [-5, 5] on s(n); D-STAR, ProVoice and direct callers
 * the library comfort-noise level), with the model's progression kept. On the
 * codec paths DVSI's pulse train and noise generator also run through it, so
 * the next audible frame continues from where they stand. */
static void
mbe_render_muted_frame(float* aout_buf, mbe_parms* cur_mp, const mbe_parms* prev_mp, enum mbe_mute_noise mute_noise,
                       enum mbe_voiced_synthesis voiced, enum mbe_synthesis_start start) {
    if (voiced != MBE_VOICED_SYNTHESIS_SPEC) {
        mbe_advance_subframes(cur_mp, prev_mp, start == MBE_SYNTHESIS_FROM_SILENCE, voiced);
    }
    if (mute_noise == MBE_MUTE_NOISE_SPEC) {
        mbe_synthesizeUniformNoisef(aout_buf, MBE_SPEC_MUTE_NOISE_AMPLITUDE * MBE_SPEECH_OUTPUT_GAIN);
    } else {
        mbe_synthesizeComfortNoisef(aout_buf);
    }
}

/* The frame as TIA-102.BABA synthesizes it: s(n) before the output scale. */
static void
mbe_render_spec_frame(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, int N) {
    /* Algorithm #117: generate the unvoiced noise buffer first, preserving
     * per-stream overlap and LCG consumption independently of voiced phase. */
    float noise_buffer[256];
    const int subframe_overlap = mbe_spec_frame_noise(noise_buffer, cur_mp, prev_mp);

    /* eq 128 and 129: Handle different L values between frames */
    int maxl = mbe_reconcile_speech_model_lengths(cur_mp, prev_mp);

    /* Integrate pulse position and regenerate phase from the spectral envelope. */
    mbe_update_speech_phases(cur_mp, prev_mp, N);

    /* Render with the low-band roll-off and the voicing gains; both models
     * keep their decoded amplitudes and phases. A caller may pass one model as
     * both frames: prepare it once. */
    struct mbe_render_saved cur_saved, prev_saved;
    const int shared_model = prev_mp == cur_mp;
    mbe_render_prepare(cur_mp, &cur_saved, maxl, MBE_RENDER_VOICED_GAIN);
    if (!shared_model) {
        mbe_render_prepare(prev_mp, &prev_saved, maxl, MBE_RENDER_VOICED_GAIN);
    }

    mbe_render_voiced_speech(aout_buf, cur_mp, prev_mp, maxl, N);

    /* Unvoiced components from the same per-frame noise buffer (TIA-102.BABA
     * equations 117-126). */
    mbe_fft_plan* plan = mbe_get_fft_plan();
    if (plan) {
        mbe_synthesizeUnvoicedFFTFrom(aout_buf, cur_mp, prev_mp, plan, noise_buffer, subframe_overlap);
    }

    mbe_render_restore(cur_mp, &cur_saved);
    if (!shared_model) {
        mbe_render_restore(prev_mp, &prev_saved);
    }
}

/* The way a frame's rendering starts (mbe_synthesis_start) and what replaces
 * it if the core mutes it. */
struct mbe_speech_options {
    enum mbe_mute_noise mute_noise;
    enum mbe_voiced_synthesis voiced;
    enum mbe_amplitude_limit limit;
    enum mbe_synthesis_start start;
    enum mbe_output_limit output;
    /* A gain on the magnitudes the frame renders, after adaptive smoothing
     * (P25's quiet frames, mbe_synthesizeSpeechScaledf()); NULL for none. */
    const double* magnitude_gain;
};

/* Every magnitude of the model the core renders or advances, times gain: each
 * product one binary64 multiplication, stored (a volatile double, so that no
 * x87 build converts the extended product to float directly, rounding it once
 * where the others round it twice), then rounded to float. On x87 the
 * multiplication is formed in extended precision and rounded to binary64 on
 * its store: the SSE result but where that first rounding lands exactly
 * halfway between two doubles (double rounding), which a product of a float
 * and a double whose exact bits span more than 64 places can do; the IEEE
 * builds agree, an x87 build can differ by one float step where the product
 * lies at a float's rounding edge. */
static void
mbe_scale_rendered_magnitudes(mbe_parms* cur_mp, double gain) {
    for (int l = 1; l <= cur_mp->L; l++) {
        volatile double product = (double)cur_mp->Ml[l] * gain;
        cur_mp->Ml[l] = mbe_round_to_float((float)product);
    }
}

static void
mbe_synthesizeSpeechCore(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, enum mbe_speech_smoothing smoothing,
                         float pre_enh_rm0, const struct mbe_speech_options* options) {

    /* One frame length for the phase arithmetic, the window Ws and the
     * voiced kernels, which process exactly MBE_VOICED_FRAME samples. */
    const int N = MBE_VOICED_FRAME;

    if (MBE_UNLIKELY(!aout_buf)) {
        return;
    }
    if (MBE_UNLIKELY(!cur_mp || !prev_mp || !mbe_harmonic_count_is_valid(cur_mp->L)
                     || !mbe_harmonic_count_is_valid(prev_mp->L))) {
        mbe_synthesizeSilencef(aout_buf);
        return;
    }

    /* TIA-102.BABA equations 111-116: advance smoothing state before the
     * mute decision so a muted frame retains its parameter progression. */
    if (smoothing == MBE_SMOOTHING_PRE_ENH_RM0) {
        mbe_applyAdaptiveSmoothingWithRm0Limit(cur_mp, prev_mp, pre_enh_rm0, options->limit);
    } else if (smoothing == MBE_SMOOTHING_FROM_MODEL) {
        mbe_applyAdaptiveSmoothing(cur_mp, prev_mp);
    }
    /* After smoothing and before the mute decision: the magnitudes smoothing
     * left, scaled, are the model rendered or advanced unheard alike. */
    if (options->magnitude_gain != NULL) {
        mbe_scale_rendered_magnitudes(cur_mp, *options->magnitude_gain);
    }

    /* The shared core checks repeat count and the IMBE error-rate threshold.
     * AMBE codec paths apply their own additional muting rules. */
    if (mbe_should_mute_speech(cur_mp)) {
        mbe_render_muted_frame(aout_buf, cur_mp, prev_mp, options->mute_noise, options->voiced, options->start);
        return;
    }

    memset(aout_buf, 0, N * sizeof(float));
    if (options->voiced != MBE_VOICED_SYNTHESIS_SPEC) {
        /* The codec paths render the frame in DVSI's subframes, each its codec's way. */
        mbe_synthesize_subframes(aout_buf, cur_mp, prev_mp, mbe_get_fft_plan(), MBE_RENDER_VOICED_GAIN_DVSI,
                                 options->start == MBE_SYNTHESIS_FROM_SILENCE, options->voiced);
    } else {
        mbe_render_spec_frame(aout_buf, cur_mp, prev_mp, N);
    }

    /* Convert s(n) to output scale, then limit overloads to 95% headroom,
     * unless the caller mixes another signal in and limits the sum. */
    for (int i = 0; i < N; i++) {
        aout_buf[i] *= MBE_SPEECH_OUTPUT_GAIN;
    }
    if (options->output == MBE_OUTPUT_LIMITED) {
        mbe_clipFloatBuffer(aout_buf, N);
    }
}

void
mbe_synthesizeSpeechScaledf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, float rm0, double magnitude_gain,
                            enum mbe_mute_noise mute_noise, enum mbe_voiced_synthesis voiced,
                            enum mbe_amplitude_limit limit) {
    const struct mbe_speech_options options = {mute_noise,         voiced,         limit, MBE_SYNTHESIS_CONTINUES,
                                               MBE_OUTPUT_LIMITED, &magnitude_gain};
    mbe_synthesizeSpeechCore(aout_buf, cur_mp, prev_mp, MBE_SMOOTHING_PRE_ENH_RM0, rm0, &options);
}

void
mbe_synthesizeSpeechWithPreEnhRm0Startf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, float rm0,
                                        enum mbe_mute_noise mute_noise, enum mbe_voiced_synthesis voiced,
                                        enum mbe_amplitude_limit limit, enum mbe_synthesis_start start,
                                        enum mbe_output_limit output) {
    const struct mbe_speech_options options = {mute_noise, voiced, limit, start, output, NULL};
    mbe_synthesizeSpeechCore(aout_buf, cur_mp, prev_mp, MBE_SMOOTHING_PRE_ENH_RM0, rm0, &options);
}

void
mbe_synthesizeSpeechWithPreEnhRm0f(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp, float rm0,
                                   enum mbe_mute_noise mute_noise, enum mbe_voiced_synthesis voiced,
                                   enum mbe_amplitude_limit limit) {
    mbe_synthesizeSpeechWithPreEnhRm0Startf(aout_buf, cur_mp, prev_mp, rm0, mute_noise, voiced, limit,
                                            MBE_SYNTHESIS_CONTINUES, MBE_OUTPUT_LIMITED);
}

void
mbe_synthesizeRepeatedSpeechStartf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp,
                                   enum mbe_voiced_synthesis voiced, enum mbe_synthesis_start start,
                                   enum mbe_output_limit output) {
    /* Repeat paths, which never mute inside synthesis. */
    const struct mbe_speech_options options = {
        MBE_MUTE_NOISE_SPEC, voiced, MBE_AMPLITUDE_LIMIT_SPEC, start, output, NULL};
    mbe_synthesizeSpeechCore(aout_buf, cur_mp, prev_mp, MBE_SMOOTHING_NONE, 0.0f, &options);
}

void
mbe_synthesizeRepeatedSpeechf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp,
                              enum mbe_voiced_synthesis voiced) {
    mbe_synthesizeRepeatedSpeechStartf(aout_buf, cur_mp, prev_mp, voiced, MBE_SYNTHESIS_CONTINUES, MBE_OUTPUT_LIMITED);
}

void
mbe_synthesizeSpeechf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp) {
    const struct mbe_speech_options options = {MBE_MUTE_NOISE_COMFORT,   MBE_VOICED_SYNTHESIS_SPEC,
                                               MBE_AMPLITUDE_LIMIT_SPEC, MBE_SYNTHESIS_CONTINUES,
                                               MBE_OUTPUT_LIMITED,       NULL};
    mbe_synthesizeSpeechCore(aout_buf, cur_mp, prev_mp, MBE_SMOOTHING_FROM_MODEL, 0.0f, &options);
}

/**
 * @brief Synthesize one speech frame into 160 16-bit samples at 8 kHz.
 * @param aout_buf Output buffer of 160 16-bit samples.
 * @param cur_mp   Current parameter set.
 * @param prev_mp  Previous parameter set.
 */
void
mbe_synthesizeSpeech(short* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp) {
    float float_buf[160];

    if (MBE_UNLIKELY(!aout_buf)) {
        return;
    }
    mbe_synthesizeSpeechf(float_buf, cur_mp, prev_mp);
    mbe_floattoshort(float_buf, aout_buf);
}
