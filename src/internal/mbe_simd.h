// SPDX-License-Identifier: GPL-2.0-or-later
/* Private compile-target vector operations; no runtime ISA requirement added. */
#ifndef MBE_INTERNAL_SIMD_H
#define MBE_INTERNAL_SIMD_H
#include "mbe_compiler.h"
#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#include <emmintrin.h>
#if defined(__AVX2__)
#include <immintrin.h>
typedef __m256d mbe_vd;
#define MBE_VD_WIDTH 4
#define mbe_vd_set   _mm256_set1_pd
#define mbe_vd_add   _mm256_add_pd
#define mbe_vd_sub   _mm256_sub_pd
#define mbe_vd_mul   _mm256_mul_pd
#define mbe_vd_load  _mm256_loadu_pd
#define mbe_vd_store _mm256_storeu_pd
#else
typedef __m128d mbe_vd;
#define MBE_VD_WIDTH 2
#define mbe_vd_set   _mm_set1_pd
#define mbe_vd_add   _mm_add_pd
#define mbe_vd_sub   _mm_sub_pd
#define mbe_vd_mul   _mm_mul_pd
#define mbe_vd_load  _mm_loadu_pd
#define mbe_vd_store _mm_storeu_pd
#endif
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON) && defined(MBE_ARCH_AARCH64)
#include <arm_neon.h>
typedef float64x2_t mbe_vd;
#define MBE_VD_WIDTH 2
#define mbe_vd_set   vdupq_n_f64
#define mbe_vd_add   vaddq_f64
#define mbe_vd_sub   vsubq_f64
#define mbe_vd_mul   vmulq_f64
#define mbe_vd_load  vld1q_f64
#define mbe_vd_store vst1q_f64
#else
/* ARM32 has no double-precision NEON. */
typedef double mbe_vd;
#define MBE_VD_WIDTH 1

static inline mbe_vd
mbe_vd_set(double x) {
    return x;
}

static inline mbe_vd
mbe_vd_add(mbe_vd a, mbe_vd b) {
    return a + b;
}

static inline mbe_vd
mbe_vd_sub(mbe_vd a, mbe_vd b) {
    return a - b;
}

static inline mbe_vd
mbe_vd_mul(mbe_vd a, mbe_vd b) {
    return a * b;
}

static inline mbe_vd
mbe_vd_load(const double* x) {
    return *x;
}

static inline void
mbe_vd_store(double* x, mbe_vd a) {
    *x = a;
}
#endif

#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
typedef __m128 mbe_vf;
#define MBE_VF_WIDTH 4
#define mbe_vf_set   _mm_set1_ps
#define mbe_vf_add   _mm_add_ps
#define mbe_vf_sub   _mm_sub_ps
#define mbe_vf_mul   _mm_mul_ps
#define mbe_vf_load  _mm_loadu_ps
#define mbe_vf_store _mm_storeu_ps
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h>
typedef float32x4_t mbe_vf;
#define MBE_VF_WIDTH 4
#define mbe_vf_set   vdupq_n_f32
#define mbe_vf_add   vaddq_f32
#define mbe_vf_sub   vsubq_f32
#define mbe_vf_mul   vmulq_f32
#define mbe_vf_load  vld1q_f32
#define mbe_vf_store vst1q_f32
#endif

static inline double
mbe_vd_sum(mbe_vd v) {
    double lanes[MBE_VD_WIDTH];
    mbe_vd_store(lanes, v);
    double sum = 0.0;
    for (int i = 0; i < MBE_VD_WIDTH; ++i) {
        sum += lanes[i];
    }
    return sum;
}

static inline mbe_vd
mbe_vd_load_float(const float* x) {
#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#if defined(__AVX2__)
    return _mm256_cvtps_pd(_mm_loadu_ps(x));
#else
    return _mm_cvtps_pd(_mm_loadl_pi(_mm_setzero_ps(), (const __m64*)x));
#endif
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON) && defined(MBE_ARCH_AARCH64)
    return vcvt_f64_f32(vld1_f32(x));
#else
    return (double)*x;
#endif
}
#endif
