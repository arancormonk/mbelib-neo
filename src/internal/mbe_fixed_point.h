// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Saturating 16/32-bit fixed-point operators for the TETRA ACELP decoder.
 *
 * The usual speech-codec operator set (add, sub, mult, L_mult, L_mac, norm_l, ...). Each operator
 * is written for portable, defined C99 behaviour: sums are formed in 64 bits and saturated, right
 * shifts of negative values are floor divisions, and nothing relies on signed overflow,
 * left-shifting negative values or implementation-defined conversions. Results are therefore the
 * same on every compiler and architecture.
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_FIXED_POINT_H
#define MBELIB_NEO_INTERNAL_MBE_FIXED_POINT_H

#include <stdint.h>

#define MBE_FX_MAX16 32767
#define MBE_FX_MIN16 (-32768)
#define MBE_FX_MAX32 INT32_C(2147483647)
#define MBE_FX_MIN32 (-INT32_C(2147483647) - 1)

/** Saturate to the int16 range. */
static inline int16_t
mbe_fx_sat16(int64_t x) {
    if (x > MBE_FX_MAX16) {
        return MBE_FX_MAX16;
    }
    if (x < MBE_FX_MIN16) {
        return MBE_FX_MIN16;
    }
    return (int16_t)x;
}

/** Saturate to the int32 range. */
static inline int32_t
mbe_fx_sat32(int64_t x) {
    if (x > MBE_FX_MAX32) {
        return MBE_FX_MAX32;
    }
    if (x < MBE_FX_MIN32) {
        return MBE_FX_MIN32;
    }
    return (int32_t)x;
}

/** Floor of x / 2^n for 0 <= n <= 63 (other counts are clamped), without shifting a negative value. */
static inline int64_t
mbe_fx_asr64(int64_t x, int n) {
    if (n <= 0) {
        return x;
    }
    if (n > 63) {
        n = 63;
    }
    if (x >= 0) {
        return x >> n;
    }
    return -((-(x + 1)) >> n) - 1;
}

/** x * 2^n for 0 <= n <= 62 when the result fits; callers saturate. */
static inline int64_t
mbe_fx_lsl64(int64_t x, int n) {
    return x * ((int64_t)1 << n);
}

static inline int16_t
mbe_fx_add(int16_t a, int16_t b) {
    return mbe_fx_sat16((int64_t)a + b);
}

static inline int16_t
mbe_fx_sub(int16_t a, int16_t b) {
    return mbe_fx_sat16((int64_t)a - b);
}

static inline int16_t
mbe_fx_negate(int16_t a) {
    return mbe_fx_sat16(-(int64_t)a);
}

/** Saturating left shift by n (negative n shifts right). */
static inline int16_t
mbe_fx_shl(int16_t a, int n) {
    if (n < 0) {
        return (int16_t)mbe_fx_asr64(a, (n < -15) ? 15 : -n);
    }
    if (n > 16) {
        n = 16;
    }
    return mbe_fx_sat16(mbe_fx_lsl64(a, n));
}

/** Arithmetic right shift by n (negative n shifts left with saturation). */
static inline int16_t
mbe_fx_shr(int16_t a, int n) {
    if (n < 0) {
        return mbe_fx_shl(a, (n < -16) ? 16 : -n);
    }
    return (int16_t)mbe_fx_asr64(a, (n > 15) ? 15 : n);
}

/** (a * b) >> 15, truncating, saturated (only -32768 * -32768 saturates). */
static inline int16_t
mbe_fx_mult(int16_t a, int16_t b) {
    return mbe_fx_sat16(mbe_fx_asr64((int64_t)a * b, 15));
}

/** (a * b + 2^14) >> 15, saturated. */
static inline int16_t
mbe_fx_mult_r(int16_t a, int16_t b) {
    return mbe_fx_sat16(mbe_fx_asr64((int64_t)a * b + 0x4000, 15));
}

/** 2 * a * b, saturated. */
static inline int32_t
mbe_fx_L_mult(int16_t a, int16_t b) {
    return mbe_fx_sat32((int64_t)a * b * 2);
}

/** a * b (always fits in 32 bits). */
static inline int32_t
mbe_fx_L_mult0(int16_t a, int16_t b) {
    return (int32_t)((int64_t)a * b);
}

static inline int32_t
mbe_fx_L_add(int32_t a, int32_t b) {
    return mbe_fx_sat32((int64_t)a + b);
}

static inline int32_t
mbe_fx_L_sub(int32_t a, int32_t b) {
    return mbe_fx_sat32((int64_t)a - b);
}

static inline int32_t
mbe_fx_L_mac(int32_t acc, int16_t a, int16_t b) {
    return mbe_fx_L_add(acc, mbe_fx_L_mult(a, b));
}

static inline int32_t
mbe_fx_L_msu(int32_t acc, int16_t a, int16_t b) {
    return mbe_fx_L_sub(acc, mbe_fx_L_mult(a, b));
}

static inline int32_t
mbe_fx_L_mac0(int32_t acc, int16_t a, int16_t b) {
    return mbe_fx_sat32((int64_t)acc + (int64_t)a * b);
}

static inline int32_t
mbe_fx_L_msu0(int32_t acc, int16_t a, int16_t b) {
    return mbe_fx_sat32((int64_t)acc - (int64_t)a * b);
}

/** Saturating left shift by n (negative n shifts right). */
static inline int32_t
mbe_fx_L_shl(int32_t a, int n) {
    if (n < 0) {
        return (int32_t)mbe_fx_asr64(a, (n < -31) ? 31 : -n);
    }
    if (n > 32) {
        n = 32;
    }
    return mbe_fx_sat32(mbe_fx_lsl64(a, n));
}

/** Arithmetic right shift by n (negative n shifts left with saturation). */
static inline int32_t
mbe_fx_L_shr(int32_t a, int n) {
    if (n < 0) {
        return mbe_fx_L_shl(a, (n < -32) ? 32 : -n);
    }
    return (int32_t)mbe_fx_asr64(a, (n > 31) ? 31 : n);
}

/** Right shift by n rounding half up; n > 31 gives 0, negative n shifts left. */
static inline int32_t
mbe_fx_L_shr_r(int32_t a, int n) {
    if (n > 31) {
        return 0;
    }
    if (n <= 0) {
        return mbe_fx_L_shr(a, n);
    }
    return (int32_t)mbe_fx_asr64((int64_t)a + ((int64_t)1 << (n - 1)), n);
}

/** High 16 bits (floor of a / 65536). */
static inline int16_t
mbe_fx_extract_h(int32_t a) {
    return (int16_t)mbe_fx_asr64(a, 16);
}

/** Low 16 bits, reinterpreted as signed. */
static inline int16_t
mbe_fx_extract_l(int32_t a) {
    int32_t v = (int32_t)((uint32_t)a & 0xFFFFu);
    return (int16_t)(v >= 0x8000 ? v - 0x10000 : v);
}

/** extract_h(L_add(a, 0x8000)). */
static inline int16_t
mbe_fx_round(int32_t a) {
    return mbe_fx_extract_h(mbe_fx_L_add(a, 0x8000));
}

/** a << 16 as a 32-bit value. */
static inline int32_t
mbe_fx_L_deposit_h(int16_t a) {
    return (int32_t)mbe_fx_lsl64(a, 16);
}

/** Left shifts that normalise a (0 for 0, 31 for -1). */
static inline int
mbe_fx_norm_l(int32_t a) {
    int n = 0;
    int64_t v;
    if (a == 0) {
        return 0;
    }
    if (a == -1) {
        return 31;
    }
    v = (a < 0) ? -((int64_t)a + 1) : a;
    while (v < 0x40000000) {
        v *= 2;
        n++;
    }
    return n;
}

/** Split a 32-bit value into hi (Q31 high half) and lo (15 low bits) parts. */
static inline void
mbe_fx_L_extract(int32_t a, int16_t* hi, int16_t* lo) {
    *hi = mbe_fx_extract_h(a);
    *lo = mbe_fx_extract_l(mbe_fx_L_msu(mbe_fx_L_shr(a, 1), *hi, 16384));
}

#endif /* MBELIB_NEO_INTERNAL_MBE_FIXED_POINT_H */
