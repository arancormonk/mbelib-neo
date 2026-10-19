// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief AVX2 adaptive-codebook interpolation of the TETRA ACELP decoder (tetra_acelp_dsp.c).
 *
 * Eight output lanes, taps in their original order and saturation after every MAC past the safe
 * prefix: the same integers as the four-lane and scalar loops. CMake compiles this unit with AVX2 code
 * generation where the library dispatches at run time (mbe_cpu.h); a target with AVX2 builds it as is.
 * Elsewhere it is empty.
 */

#include "mbe_compiler.h"
#include "tetra_simd.h" // IWYU pragma: keep (TETRA_SIMD and tetra_vadds8)

#if defined(TETRA_SIMD) && defined(__AVX2__) && defined(MBE_SIMD_TARGET_SSE2)
#include <stdint.h>

#include "tetra_acelp.h"

int
tetra_adaptive_avx2(int16_t* exc_cur, int t0, int first, const int16_t* taps) {
    int n = 0;
    if (t0 - (first + 31) < 8) {
        return 0; /* Short lags need the four-lane kernel's earlier stores. */
    }
    for (; n + 8 <= TETRA_SUBFRAME; n += 8) {
        __m256i acc = _mm256_setzero_si256();
        int k;
        for (k = 0; k < TETRA_INTERP_SAFE_TAPS; ++k) {
            __m128i x = _mm_loadu_si128((const __m128i*)(exc_cur + n - t0 + first + k));
            __m256i p = _mm256_mullo_epi32(_mm256_cvtepi16_epi32(x), _mm256_set1_epi32(taps[k]));
            acc = _mm256_add_epi32(acc, p);
        }
        for (; k < 32; ++k) {
            __m128i x = _mm_loadu_si128((const __m128i*)(exc_cur + n - t0 + first + k));
            __m256i p = _mm256_mullo_epi32(_mm256_cvtepi16_epi32(x), _mm256_set1_epi32(taps[k]));
            acc = tetra_vadds8(acc, p);
        }
        acc = tetra_vadds8(acc, acc);
        acc = _mm256_srai_epi32(tetra_vadds8(acc, _mm256_set1_epi32(0x8000)), 16);
        _mm_storeu_si128((__m128i*)(exc_cur + n),
                         _mm_packs_epi32(_mm256_castsi256_si128(acc), _mm256_extracti128_si256(acc, 1)));
    }
    return n;
}
#else
/* ISO C forbids an empty translation unit. */
typedef int tetra_acelp_avx2_unused;
#endif
