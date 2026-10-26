// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Lightweight math helpers for performance-sensitive code paths.
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_MATH_H
#define MBELIB_NEO_INTERNAL_MBE_MATH_H

#include <float.h> // IWYU pragma: keep (FLT_EVAL_METHOD where __FLT_EVAL_METHOD__ is missing)
#include <math.h>  // IWYU pragma: keep (sinf/cosf where no builtin sincosf exists)

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

/* Whether float arithmetic may carry excess precision: x87 (GCC and Clang
 * report __FLT_EVAL_METHOD__ 2 whenever they use it, fast-math included), or
 * an evaluation method the compiler leaves open (MSVC, which lacks the macro,
 * reports FLT_EVAL_METHOD -1 under /fp:fast). */
#if defined(__FLT_EVAL_METHOD__)
#define MBE_FLOAT_EXCESS_PRECISION (__FLT_EVAL_METHOD__ != 0)
#elif defined(FLT_EVAL_METHOD)
#define MBE_FLOAT_EXCESS_PRECISION (FLT_EVAL_METHOD != 0)
#else
#define MBE_FLOAT_EXCESS_PRECISION 0
#endif

/**
 * @brief Compute both sine and cosine of an angle.
 *
 * Uses a combined intrinsic when available for better performance and
 * accuracy; otherwise falls back to separate `sinf` and `cosf` calls.
 *
 * @param x Input angle in radians.
 * @param s Output pointer receiving `sinf(x)`.
 * @param c Output pointer receiving `cosf(x)`.
 */
static inline void
mbe_sincosf(float x, float* s, float* c) {
    /* GCC provides __builtin_sincosf even where __has_builtin cannot say so. */
#if __has_builtin(__builtin_sincosf) || (defined(__GNUC__) && !defined(__clang__))
    __builtin_sincosf(x, s, c);
#else
    *s = sinf(x);
    *c = cosf(x);
#endif
}

/**
 * @brief x rounded to float.
 *
 * Under x87 excess precision (MBE_FLOAT_EXCESS_PRECISION) the compiler may keep
 * a float result's extra bits, even through a store it can see, or round it,
 * case by case, so a comparison against it can fall on either side. Storing
 * through a volatile rounds it as plain float arithmetic does; elsewhere this
 * is the value itself.
 */
static inline float
mbe_round_to_float(float x) {
#if MBE_FLOAT_EXCESS_PRECISION
    volatile float rounded = x;
    return rounded;
#else
    return x;
#endif
}

/**
 * @brief Position (prev_L / L) * l at which spectral-amplitude prediction
 * interpolates the previous frame's log2 amplitudes.
 *
 * Decoders and encoders truncate this position to pick the interpolation taps,
 * and an exact integer position lands on either side of the truncation
 * depending on rounding. Plain float arithmetic rounds each step. Under x87
 * excess precision (MBE_FLOAT_EXCESS_PRECISION) the compiler may keep extra
 * bits or not, case by case, so a decoder and an encoder could pick different
 * taps; there each step is rounded to float explicitly, which gives the same
 * value as plain float arithmetic.
 */
static inline float
mbe_prediction_position(int prev_L, int L, int l) {
    const float ratio = mbe_round_to_float((float)prev_L / (float)L);
    return mbe_round_to_float(ratio * (float)l);
}

#endif /* MBELIB_NEO_INTERNAL_MBE_MATH_H */
