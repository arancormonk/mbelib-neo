// SPDX-License-Identifier: GPL-2.0-or-later
/* Private fixed-point SIMD for TETRA. Four independent samples per vector; no ISA dispatch. */
#ifndef MBELIB_NEO_INTERNAL_TETRA_SIMD_H
#define MBELIB_NEO_INTERNAL_TETRA_SIMD_H

#include "mbe_compiler.h"

#if defined(MBELIB_ENABLE_SIMD) && (defined(MBE_SIMD_TARGET_SSE2) || defined(MBE_SIMD_TARGET_NEON))
#include <stdint.h>
#include <string.h>
#endif

#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#include <emmintrin.h> // IWYU pragma: export
#include <mmintrin.h>  // IWYU pragma: export
#include <xmmintrin.h> // IWYU pragma: export
#if defined(__AVX2__)
#include <immintrin.h> // IWYU pragma: export
#endif
#define TETRA_SIMD 1
typedef __m128i tetra_v32;
#define tetra_vset _mm_set1_epi32
#define tetra_vadd _mm_add_epi32
#define tetra_vsub _mm_sub_epi32
#define tetra_vshl _mm_slli_epi32
#define tetra_vshr _mm_srai_epi32

static inline tetra_v32
tetra_vload(const int16_t* x) {
    __m128i v = _mm_loadl_epi64((const __m128i*)x);
    return _mm_unpacklo_epi16(v, _mm_cmpgt_epi16(_mm_setzero_si128(), v));
}

static inline tetra_v32
tetra_vmul(const int16_t* x, int16_t c) {
    __m128i v = _mm_loadl_epi64((const __m128i*)x);
    return _mm_madd_epi16(_mm_unpacklo_epi16(v, _mm_setzero_si128()), _mm_set1_epi32((uint16_t)c));
}

static inline tetra_v32
tetra_vmulpair(const int16_t* x, const int16_t* y) {
    __m128i a = _mm_loadl_epi64((const __m128i*)x);
    __m128i b = _mm_loadl_epi64((const __m128i*)y);
    return _mm_unpacklo_epi16(_mm_mullo_epi16(a, b), _mm_mulhi_epi16(a, b));
}

static inline tetra_v32
tetra_vadds(tetra_v32 a, tetra_v32 b) {
    __m128i sum = _mm_add_epi32(a, b);
    __m128i overflow = _mm_srai_epi32(_mm_andnot_si128(_mm_xor_si128(a, b), _mm_xor_si128(a, sum)), 31);
    __m128i limit = _mm_xor_si128(_mm_srai_epi32(a, 31), _mm_set1_epi32(INT32_MAX));
    return _mm_or_si128(_mm_and_si128(overflow, limit), _mm_andnot_si128(overflow, sum));
}

#if defined(__AVX2__)
static inline __m256i
tetra_vadds8(__m256i a, __m256i b) {
    __m256i sum = _mm256_add_epi32(a, b);
    __m256i overflow = _mm256_srai_epi32(_mm256_andnot_si256(_mm256_xor_si256(a, b), _mm256_xor_si256(a, sum)), 31);
    __m256i limit = _mm256_xor_si256(_mm256_srai_epi32(a, 31), _mm256_set1_epi32(INT32_MAX));
    return _mm256_blendv_epi8(sum, limit, overflow);
}
#endif

static inline void
tetra_vstore_sat(int16_t* x, tetra_v32 v) {
    __m128i packed = _mm_packs_epi32(v, v);
    memcpy(x, &packed, 4 * sizeof(*x));
}

static inline void
tetra_vstore_wrap(int16_t* x, tetra_v32 v) {
    tetra_vstore_sat(x, _mm_srai_epi32(_mm_slli_epi32(v, 16), 16));
}

/* Caller proves sum(abs(a)) <= 61439. Every partial sum, including x * 4096,
 * then fits in int32 even when all history samples are -32768. */
static inline int32_t
tetra_dot10(const int16_t* a, const int16_t* x) {
    __m128i p = _mm_madd_epi16(_mm_loadu_si128((const __m128i*)a), _mm_loadu_si128((const __m128i*)x));
    p = _mm_add_epi32(p, _mm_shuffle_epi32(p, _MM_SHUFFLE(1, 0, 3, 2)));
    p = _mm_add_epi32(p, _mm_shuffle_epi32(p, _MM_SHUFFLE(2, 3, 0, 1)));
    return _mm_cvtsi128_si32(p) + (int32_t)a[8] * x[8] + (int32_t)a[9] * x[9];
}

#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h> // IWYU pragma: export
#define TETRA_SIMD  1
typedef int32x4_t tetra_v32;
#define tetra_vset  vdupq_n_s32
#define tetra_vadd  vaddq_s32
#define tetra_vsub  vsubq_s32
#define tetra_vshl  vshlq_n_s32
#define tetra_vshr  vshrq_n_s32
#define tetra_vadds vqaddq_s32

static inline tetra_v32
tetra_vload(const int16_t* x) {
    return vmovl_s16(vld1_s16(x));
}

static inline tetra_v32
tetra_vmul(const int16_t* x, int16_t c) {
    return vmull_n_s16(vld1_s16(x), c);
}

static inline tetra_v32
tetra_vmulpair(const int16_t* x, const int16_t* y) {
    return vmull_s16(vld1_s16(x), vld1_s16(y));
}

static inline void
tetra_vstore_sat(int16_t* x, tetra_v32 v) {
    int16x4_t packed = vqmovn_s32(v);
    memcpy(x, &packed, sizeof(packed));
}

static inline void
tetra_vstore_wrap(int16_t* x, tetra_v32 v) {
    int16x4_t packed = vmovn_s32(v);
    memcpy(x, &packed, sizeof(packed));
}

/* Same coefficient bound as the SSE2 reduction. */
static inline int32_t
tetra_dot10(const int16_t* a, const int16_t* x) {
    int32x4_t p = vmull_s16(vld1_s16(a), vld1_s16(x));
    p = vmlal_s16(p, vld1_s16(a + 4), vld1_s16(x + 4));
    int32x2_t s = vadd_s32(vget_low_s32(p), vget_high_s32(p));
    s = vpadd_s32(s, s);
    return vget_lane_s32(s, 0) + (int32_t)a[8] * x[8] + (int32_t)a[9] * x[9];
}
#endif

#endif /* MBELIB_NEO_INTERNAL_TETRA_SIMD_H */
