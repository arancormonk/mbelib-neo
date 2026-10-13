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
 * @brief Float-to-int16 PCM conversion with scalar and SIMD implementations.
 *
 * Compile this translation unit with IEEE floating-point semantics, even when
 * synthesis uses fast-math, so NaN/Inf classification survives optimization.
 * Under finite-math-only assumptions the compiler may remove the bit-pattern
 * NaN/Inf checks below, so the unit refuses to build that way.
 */

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__) || defined(_M_FP_FAST)
#error "mbe_pcm_convert.c must be compiled with IEEE floating-point semantics (no fast-math)"
#endif

#include <stdint.h>
#include <string.h>
#ifdef MBE_DEBUG
#include <stdio.h>
#endif
#include "mbe_compiler.h"
#if defined(MBELIB_ENABLE_SIMD)
#if defined(MBE_SIMD_TARGET_SSE2)
#include <emmintrin.h>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
#endif
#if defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h>
#endif
#endif

#include "mbelib-neo/mbelib.h"

#if !defined(MBELIB_ENABLE_SIMD) || (!defined(MBE_SIMD_TARGET_SSE2) && !defined(MBE_SIMD_TARGET_NEON))
/**
 * @brief Scale one float sample by the historical gain and clip it for int16.
 *
 * NaN converts to zero; signed infinities clip to the signed limit.
 * @param sample Input float sample (historical mbelib scale, not normalized).
 * @return Scaled and clipped value that fits in an int16.
 */
static inline float
mbe_scaleAndClipFloatSampleForShort(float sample) {
    const float again = 7.0f;
    const float max_amplitude = 32767.0f * 0.95f; /* ~31128.65 */
    uint32_t sample_bits;
    memcpy(&sample_bits, &sample, sizeof(sample_bits));

    const uint32_t abs_bits = sample_bits & 0x7FFFFFFFu;
    if (abs_bits > 0x7F800000u) {
        return 0.0f;
    }
    if (abs_bits == 0x7F800000u) {
        return ((sample_bits & 0x80000000u) != 0u) ? -max_amplitude : max_amplitude;
    }

    float audio = again * sample;
    if (audio > max_amplitude) {
#ifdef MBE_DEBUG
        fprintf(stderr, "audio clip: %f\n", audio);
#endif
        return max_amplitude;
    }
    if (audio < -max_amplitude) {
#ifdef MBE_DEBUG
        fprintf(stderr, "audio clip: %f\n", audio);
#endif
        return -max_amplitude;
    }
    return audio;
}
#endif

#if defined(MBELIB_ENABLE_SIMD)
/**
 * @brief SSE2 specialization for float→int16 conversion.
 */
#if defined(MBE_SIMD_TARGET_SSE2)
/* Classify in integer lanes before any floating-point arithmetic, including
 * for signaling NaNs. Infinities survive and clip to the signed limit. */
#if defined(__AVX2__)
static inline __m256
mbe_pcm_nan_to_zero_avx2(__m256 x) {
    __m256i magnitude = _mm256_and_si256(_mm256_castps_si256(x), _mm256_set1_epi32(0x7fffffff));
    __m256i nan = _mm256_cmpgt_epi32(magnitude, _mm256_set1_epi32(0x7f800000));
    return _mm256_andnot_ps(_mm256_castsi256_ps(nan), x);
}

static void
mbe_floattoshort_avx2(const float* restrict float_buf, short* restrict aout_buf) {
    const __m256 scale = _mm256_set1_ps(7.0f);
    const __m256 hi = _mm256_set1_ps(32767.0f * 0.95f);
    const __m256 lo = _mm256_set1_ps(-32767.0f * 0.95f);
    for (int i = 0; i < 160; i += 16) {
        __m256 a = _mm256_mul_ps(mbe_pcm_nan_to_zero_avx2(_mm256_loadu_ps(float_buf + i)), scale);
        __m256 b = _mm256_mul_ps(mbe_pcm_nan_to_zero_avx2(_mm256_loadu_ps(float_buf + i + 8)), scale);
        a = _mm256_min_ps(_mm256_max_ps(a, lo), hi);
        b = _mm256_min_ps(_mm256_max_ps(b, lo), hi);
        __m256i packed = _mm256_packs_epi32(_mm256_cvttps_epi32(a), _mm256_cvttps_epi32(b));
        /* packs works within 128-bit halves: restore consecutive samples. */
        packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3, 1, 2, 0));
        _mm256_storeu_si256((__m256i*)(aout_buf + i), packed);
    }
}
#else
static inline __m128
mbe_pcm_nan_to_zero_sse2(__m128 x) {
    __m128i magnitude = _mm_and_si128(_mm_castps_si128(x), _mm_set1_epi32(0x7fffffff));
    __m128i nan = _mm_cmpgt_epi32(magnitude, _mm_set1_epi32(0x7f800000));
    return _mm_andnot_ps(_mm_castsi128_ps(nan), x);
}

static void
mbe_floattoshort_sse2(const float* restrict float_buf, short* restrict aout_buf) {
    /* JMBE-compatible soft clipping at 95% of maximum amplitude */
    const __m128 vscale = _mm_set1_ps(7.0f);
    const __m128 vmaxv = _mm_set1_ps(32767.0f * 0.95f);
    const __m128 vminv = _mm_set1_ps(-32767.0f * 0.95f);
    for (int i = 0; i < 160; i += 8) {
        __m128 a = _mm_mul_ps(mbe_pcm_nan_to_zero_sse2(_mm_loadu_ps(float_buf + i)), vscale);
        __m128 b = _mm_mul_ps(mbe_pcm_nan_to_zero_sse2(_mm_loadu_ps(float_buf + i + 4)), vscale);
        a = _mm_min_ps(_mm_max_ps(a, vminv), vmaxv);
        b = _mm_min_ps(_mm_max_ps(b, vminv), vmaxv);
        __m128i ia = _mm_cvttps_epi32(a);
        __m128i ib = _mm_cvttps_epi32(b);
        __m128i packed = _mm_packs_epi32(ia, ib);
        _mm_storeu_si128((__m128i*)(aout_buf + i), packed);
    }
}
#endif
#endif

/**
 * @brief NEON specialization for float→int16 conversion.
 */
#if defined(MBE_SIMD_TARGET_NEON)
static inline float32x4_t
mbe_pcm_nan_to_zero_neon(float32x4_t x) {
    uint32x4_t bits = vreinterpretq_u32_f32(x);
    uint32x4_t magnitude = vandq_u32(bits, vdupq_n_u32(0x7fffffffu));
    uint32x4_t nan = vcgtq_u32(magnitude, vdupq_n_u32(0x7f800000u));
    return vreinterpretq_f32_u32(vbicq_u32(bits, nan));
}

static void
mbe_floattoshort_neon(const float* restrict float_buf, short* restrict aout_buf) {
    /* JMBE-compatible soft clipping at 95% of maximum amplitude */
    const float32x4_t vscale = vdupq_n_f32(7.0f);
    const float32x4_t vmaxv = vdupq_n_f32(32767.0f * 0.95f);
    const float32x4_t vminv = vdupq_n_f32(-32767.0f * 0.95f);
    for (int i = 0; i < 160; i += 8) {
        float32x4_t a = vmulq_f32(mbe_pcm_nan_to_zero_neon(vld1q_f32(float_buf + i)), vscale);
        float32x4_t b = vmulq_f32(mbe_pcm_nan_to_zero_neon(vld1q_f32(float_buf + i + 4)), vscale);
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
#endif /* MBELIB_ENABLE_SIMD */

/**
 * @brief Convert 160 float samples to clipped/scaled 16-bit PCM.
 *
 * Select the compiler target's SIMD specialization, with per-lane NaN/Inf
 * handling. No CPU probing or stronger ISA than the compile target is needed.
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 16-bit samples.
 */
void
mbe_floattoshort(const float* restrict float_buf, short* restrict aout_buf) {
    if (MBE_UNLIKELY(!float_buf || !aout_buf)) {
        return;
    }
#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#if defined(__AVX2__)
    mbe_floattoshort_avx2(float_buf, aout_buf);
#else
    mbe_floattoshort_sse2(float_buf, aout_buf);
#endif
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
    mbe_floattoshort_neon(float_buf, aout_buf);
#else
    for (int i = 0; i < 160; i++) {
        aout_buf[i] = (short)mbe_scaleAndClipFloatSampleForShort(float_buf[i]);
    }
#endif
}
