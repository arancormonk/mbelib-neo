// SPDX-License-Identifier: GPL-2.0-or-later
/* Independent implementation of the fixed-size PCM output contract. */
#include "mbe_pcm_convert.h"

#include "mbe_compiler.h"
#include "mbe_math.h"
#include "mbelib-neo/mbelib.h"

#include <math.h>
#include <stddef.h>

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ != 0) || defined(_M_FP_FAST)
#error "mbe_pcm_convert.c requires IEEE floating-point semantics"
#endif

#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#include <immintrin.h>
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h>
#endif

#define PCM_LIMIT (32767.0f * 0.95f)

/*
 * Two conversions share the scale, the limit and the NaN handling:
 * mbe_floattoshort() truncates toward zero, mbe_floattoshort_nearest() rounds
 * to nearest with halves away from zero. Every variant rounds the same way: it
 * truncates the limited value (at most 31128.65, so the integer is exact), takes
 * the fraction value - trunc(value), which is exact in float, and steps one away
 * from zero where the fraction's magnitude is at least a half. Adding 0.5
 * before truncating would not do: 0.49999997f + 0.5f rounds to 1 in float. The
 * vector convert instructions round halves to even or truncate, and the current
 * rounding mode plays no part.
 */
#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2) && defined(__AVX2__)
static __m256
scaled_pcm8(const float* src) {
    __m256 value = _mm256_mul_ps(_mm256_loadu_ps(src), _mm256_set1_ps(7.0f));
    value = _mm256_and_ps(value, _mm256_cmp_ps(value, value, _CMP_ORD_Q));
    value = _mm256_max_ps(value, _mm256_set1_ps(-PCM_LIMIT));
    return _mm256_min_ps(value, _mm256_set1_ps(PCM_LIMIT));
}

static __m256i
nearest_pcm8(const float* src) {
    const __m256 value = scaled_pcm8(src);
    const __m256i truncated = _mm256_cvttps_epi32(value);
    const __m256 fraction = _mm256_sub_ps(value, _mm256_cvtepi32_ps(truncated));
    /* A true comparison is -1 in every bit: subtracting it steps up, adding it steps down. */
    const __m256i up = _mm256_castps_si256(_mm256_cmp_ps(fraction, _mm256_set1_ps(0.5f), _CMP_GE_OQ));
    const __m256i down = _mm256_castps_si256(_mm256_cmp_ps(fraction, _mm256_set1_ps(-0.5f), _CMP_LE_OQ));
    return _mm256_add_epi32(_mm256_sub_epi32(truncated, up), down);
}

static void
store_pcm8(short* dest, __m256i integers) {
    const __m128i low = _mm256_castsi256_si128(integers);
    const __m128i high = _mm256_extracti128_si256(integers, 1);
    _mm_storeu_si128((__m128i*)(void*)dest, _mm_packs_epi32(low, high));
}

static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        store_pcm8(dest + i, _mm256_cvttps_epi32(scaled_pcm8(src + i)));
    }
}

static void
convert_pcm_nearest(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        store_pcm8(dest + i, nearest_pcm8(src + i));
    }
}
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
static __m128
scaled_pcm4(const float* src) {
    __m128 value = _mm_mul_ps(_mm_loadu_ps(src), _mm_set1_ps(7.0f));
    value = _mm_and_ps(value, _mm_cmpord_ps(value, value));
    value = _mm_max_ps(value, _mm_set1_ps(-PCM_LIMIT));
    return _mm_min_ps(value, _mm_set1_ps(PCM_LIMIT));
}

static __m128i
nearest_pcm4(const float* src) {
    const __m128 value = scaled_pcm4(src);
    const __m128i truncated = _mm_cvttps_epi32(value);
    const __m128 fraction = _mm_sub_ps(value, _mm_cvtepi32_ps(truncated));
    /* A true comparison is -1 in every bit: subtracting it steps up, adding it steps down. */
    const __m128i up = _mm_castps_si128(_mm_cmpge_ps(fraction, _mm_set1_ps(0.5f)));
    const __m128i down = _mm_castps_si128(_mm_cmple_ps(fraction, _mm_set1_ps(-0.5f)));
    return _mm_add_epi32(_mm_sub_epi32(truncated, up), down);
}

static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        const __m128i low = _mm_cvttps_epi32(scaled_pcm4(src + i));
        const __m128i high = _mm_cvttps_epi32(scaled_pcm4(src + i + 4));
        _mm_storeu_si128((__m128i*)(void*)(dest + i), _mm_packs_epi32(low, high));
    }
}

static void
convert_pcm_nearest(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        const __m128i low = nearest_pcm4(src + i);
        const __m128i high = nearest_pcm4(src + i + 4);
        _mm_storeu_si128((__m128i*)(void*)(dest + i), _mm_packs_epi32(low, high));
    }
}
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
static float32x4_t
limited_pcm4(const float* src) {
    float32x4_t value = vmulq_n_f32(vld1q_f32(src), 7.0f);
    value = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(value), vceqq_f32(value, value)));
    value = vmaxq_f32(value, vdupq_n_f32(-PCM_LIMIT));
    return vminq_f32(value, vdupq_n_f32(PCM_LIMIT));
}

static int32x4_t
scaled_pcm4(const float* src) {
    return vcvtq_s32_f32(limited_pcm4(src));
}

static int32x4_t
nearest_pcm4(const float* src) {
    const float32x4_t value = limited_pcm4(src);
    const int32x4_t truncated = vcvtq_s32_f32(value);
    const float32x4_t fraction = vsubq_f32(value, vcvtq_f32_s32(truncated));
    /* A true comparison is -1 in every bit: subtracting it steps up, adding it steps down. */
    const int32x4_t up = vreinterpretq_s32_u32(vcgeq_f32(fraction, vdupq_n_f32(0.5f)));
    const int32x4_t down = vreinterpretq_s32_u32(vcleq_f32(fraction, vdupq_n_f32(-0.5f)));
    return vaddq_s32(vsubq_s32(truncated, up), down);
}

static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        const int16x4_t low = vqmovn_s32(scaled_pcm4(src + i));
        const int16x4_t high = vqmovn_s32(scaled_pcm4(src + i + 4));
        vst1q_s16(dest + i, vcombine_s16(low, high));
    }
}

static void
convert_pcm_nearest(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        const int16x4_t low = vqmovn_s32(nearest_pcm4(src + i));
        const int16x4_t high = vqmovn_s32(nearest_pcm4(src + i + 4));
        vst1q_s16(dest + i, vcombine_s16(low, high));
    }
}
#else
static float
limited_sample(float sample) {
    float value = sample * 7.0f;
    if (isnan(value)) {
        value = 0.0f;
    } else if (value > PCM_LIMIT) {
        value = PCM_LIMIT;
    } else if (value < -PCM_LIMIT) {
        value = -PCM_LIMIT;
    }
    return value;
}

static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; ++i) {
        dest[i] = (short)limited_sample(src[i]);
    }
}

/* The limited value is rounded to float first: where float arithmetic keeps
 * excess precision (x87), the product 7 x of an input whose float product is
 * a half can sit just below it (0x3f249249: 4.5f, extended 4.4999999404),
 * and the fraction would then round the wrong way. The truncating
 * conversion above keeps its code as it was. */
static void
convert_pcm_nearest(const float* src, short* dest) {
    for (int i = 0; i < 160; ++i) {
        const float value = mbe_round_to_float(limited_sample(src[i]));
        const int truncated = (int)value;
        const float fraction = value - (float)truncated;
        int rounded = truncated;
        if (fraction >= 0.5f) {
            ++rounded;
        } else if (fraction <= -0.5f) {
            --rounded;
        }
        dest[i] = (short)rounded;
    }
}
#endif

void
mbe_floattoshort(const float* float_buf, short* aout_buf) {
    if (float_buf == NULL || aout_buf == NULL) {
        return;
    }
    convert_pcm(float_buf, aout_buf);
}

void
mbe_floattoshort_nearest(const float* float_buf, short* aout_buf) {
    if (float_buf == NULL || aout_buf == NULL) {
        return;
    }
    convert_pcm_nearest(float_buf, aout_buf);
}
