// SPDX-License-Identifier: GPL-2.0-or-later
/* Independent implementation of the fixed-size PCM output contract. */
#include "mbe_compiler.h"
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

#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2) && defined(__AVX2__)
static __m256
scaled_pcm8(const float* src) {
    __m256 value = _mm256_mul_ps(_mm256_loadu_ps(src), _mm256_set1_ps(7.0f));
    value = _mm256_and_ps(value, _mm256_cmp_ps(value, value, _CMP_ORD_Q));
    value = _mm256_max_ps(value, _mm256_set1_ps(-PCM_LIMIT));
    return _mm256_min_ps(value, _mm256_set1_ps(PCM_LIMIT));
}

static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        const __m256i integers = _mm256_cvttps_epi32(scaled_pcm8(src + i));
        const __m128i low = _mm256_castsi256_si128(integers);
        const __m128i high = _mm256_extracti128_si256(integers, 1);
        _mm_storeu_si128((__m128i*)(void*)(dest + i), _mm_packs_epi32(low, high));
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

static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        const __m128i low = _mm_cvttps_epi32(scaled_pcm4(src + i));
        const __m128i high = _mm_cvttps_epi32(scaled_pcm4(src + i + 4));
        _mm_storeu_si128((__m128i*)(void*)(dest + i), _mm_packs_epi32(low, high));
    }
}
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
static int32x4_t
scaled_pcm4(const float* src) {
    float32x4_t value = vmulq_n_f32(vld1q_f32(src), 7.0f);
    value = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(value), vceqq_f32(value, value)));
    value = vmaxq_f32(value, vdupq_n_f32(-PCM_LIMIT));
    value = vminq_f32(value, vdupq_n_f32(PCM_LIMIT));
    return vcvtq_s32_f32(value);
}

static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; i += 8) {
        const int16x4_t low = vqmovn_s32(scaled_pcm4(src + i));
        const int16x4_t high = vqmovn_s32(scaled_pcm4(src + i + 4));
        vst1q_s16(dest + i, vcombine_s16(low, high));
    }
}
#else
static void
convert_pcm(const float* src, short* dest) {
    for (int i = 0; i < 160; ++i) {
        float value = src[i] * 7.0f;
        if (isnan(value)) {
            value = 0.0f;
        } else if (value > PCM_LIMIT) {
            value = PCM_LIMIT;
        } else if (value < -PCM_LIMIT) {
            value = -PCM_LIMIT;
        }
        dest[i] = (short)value;
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
