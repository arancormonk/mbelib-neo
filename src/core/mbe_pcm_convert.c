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
#endif
#if defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h>
#endif
#endif

#include "mbelib-neo/mbelib.h"

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

#if defined(MBELIB_ENABLE_SIMD)
static inline int
mbe_floatSampleIsNonFinite(float sample) {
    uint32_t sample_bits;
    memcpy(&sample_bits, &sample, sizeof(sample_bits));

    return (sample_bits & 0x7FFFFFFFu) >= 0x7F800000u;
}

static int
mbe_floatBufferHasNonFinite(const float* float_buf) {
    for (int i = 0; i < 160; i++) {
        if (mbe_floatSampleIsNonFinite(float_buf[i])) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Portable scalar fallback for float→int16 conversion.
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 int16 samples.
 */
static void
mbe_floattoshort_scalar(const float* restrict float_buf, short* restrict aout_buf) {
    for (int i = 0; i < 160; i++) {
        aout_buf[i] = (short)mbe_scaleAndClipFloatSampleForShort(float_buf[i]);
    }
}

/**
 * @brief SSE2 specialization for float→int16 conversion.
 */
#if defined(MBE_SIMD_TARGET_SSE2)
static void
mbe_floattoshort_sse2(const float* restrict float_buf, short* restrict aout_buf) {
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
#if defined(MBE_SIMD_TARGET_NEON)
static void
mbe_floattoshort_neon(const float* restrict float_buf, short* restrict aout_buf) {
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

typedef void (*mbe_floattoshort_fn)(const float* restrict, short* restrict);
/*
 * Keep dispatch state thread-local so first-use initialization has no cross-thread
 * data races and still amortizes probe cost to one-time per thread.
 */
static MBE_THREAD_LOCAL mbe_floattoshort_fn mbe_floattoshort_impl = NULL; /**< Runtime-selected impl pointer. */

/**
 * @brief Initialize conversion dispatch for the compiled target ISA.
 */
static void
mbe_init_runtime_dispatch(void) {
    if (mbe_floattoshort_impl) {
        return;
    }
    /* Choose implementation first, then publish once. */
    mbe_floattoshort_fn impl = mbe_floattoshort_scalar;

#if defined(MBE_ARCH_AARCH64)
    /* NEON is mandatory on AArch64, including ARM64EC. */
#if defined(MBE_SIMD_TARGET_NEON)
    impl = mbe_floattoshort_neon;
#endif
#elif defined(MBE_ARCH_X86_64)
    /* SSE2 is guaranteed on x86_64 */
#if defined(MBE_SIMD_TARGET_SSE2)
    impl = mbe_floattoshort_sse2;
#endif
#elif defined(MBE_ARCH_X86_32)
    /*
     * 32-bit x86 only compiles the SSE2 specialization when the build targets
     * SSE2 explicitly (for example /arch:SSE2 or -msse2).
     */
#if defined(MBE_SIMD_TARGET_SSE2)
    impl = mbe_floattoshort_sse2;
#endif
#elif defined(MBE_SIMD_TARGET_NEON)
    /* Assume NEON if building with NEON intrinsics */
    impl = mbe_floattoshort_neon;
#endif
    mbe_floattoshort_impl = impl;
}

/**
 * @brief Convert 160 float samples to clipped/scaled 16-bit PCM.
 *
 * Runtime-dispatched conversion with SIMD specializations. Frames containing
 * NaN or infinity always take the scalar path so the documented handling holds.
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 16-bit samples.
 */
void
mbe_floattoshort(const float* restrict float_buf, short* restrict aout_buf) {
    if (MBE_UNLIKELY(!float_buf || !aout_buf)) {
        return;
    }
    if (MBE_UNLIKELY(mbe_floatBufferHasNonFinite(float_buf))) {
        mbe_floattoshort_scalar(float_buf, aout_buf);
        return;
    }
    if (MBE_UNLIKELY(!mbe_floattoshort_impl)) {
        mbe_init_runtime_dispatch();
    }
    mbe_floattoshort_impl(float_buf, aout_buf);
}

#else  /* MBELIB_ENABLE_SIMD not set: keep scalar implementation */
/**
 * @brief Convert 160 float samples to clipped/scaled 16-bit PCM (scalar build).
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 16-bit samples.
 */
void
mbe_floattoshort(const float* restrict float_buf, short* restrict aout_buf) {
    if (MBE_UNLIKELY(!float_buf || !aout_buf)) {
        return;
    }
    for (int i = 0; i < 160; i++) {
        aout_buf[i] = (short)mbe_scaleAndClipFloatSampleForShort(float_buf[i]);
    }
}
#endif /* MBELIB_ENABLE_SIMD */
