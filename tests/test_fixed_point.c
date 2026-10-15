// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Saturating fixed-point operators of src/internal/mbe_fixed_point.h at their edges.
 */

#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "mbe_fixed_point.h"

#define MAX16 INT16_MAX
#define MIN16 INT16_MIN
#define MAX32 INT32_MAX
#define MIN32 INT32_MIN

static void
test_saturation(void) {
    assert(mbe_fx_sat16(32768) == MAX16);
    assert(mbe_fx_sat16(-32769) == MIN16);
    assert(mbe_fx_sat16(-32768) == MIN16);
    assert(mbe_fx_sat32((int64_t)MAX32 + 1) == MAX32);
    assert(mbe_fx_sat32((int64_t)MIN32 - 1) == MIN32);

    assert(mbe_fx_add(MAX16, 1) == MAX16);
    assert(mbe_fx_add(MIN16, -1) == MIN16);
    assert(mbe_fx_add(MAX16, MIN16) == -1);
    assert(mbe_fx_sub(MIN16, 1) == MIN16);
    assert(mbe_fx_sub(0, MIN16) == MAX16);
    assert(mbe_fx_negate(MIN16) == MAX16);
    assert(mbe_fx_negate(MAX16) == -MAX16);

    assert(mbe_fx_L_add(MAX32, 1) == MAX32);
    assert(mbe_fx_L_add(MIN32, -1) == MIN32);
    assert(mbe_fx_L_sub(MIN32, 1) == MIN32);
    assert(mbe_fx_L_sub(0, MIN32) == MAX32);
    /* Saturation is applied per operation, so a saturated sum does not cancel later. */
    assert(mbe_fx_L_sub(mbe_fx_L_add(MAX32, MAX32), MAX32) == 0);
}

static void
test_shifts(void) {
    assert(mbe_fx_asr64(-1, 1) == -1);
    assert(mbe_fx_asr64(-3, 1) == -2);
    assert(mbe_fx_asr64(INT64_MIN, 62) == -2);
    assert(mbe_fx_asr64(5, 1) == 2);

    assert(mbe_fx_shl(1, 14) == 16384);
    assert(mbe_fx_shl(1, 15) == MAX16);
    assert(mbe_fx_shl(-1, 15) == MIN16);
    assert(mbe_fx_shl(-2, 15) == MIN16);
    assert(mbe_fx_shl(1, 100) == MAX16);
    assert(mbe_fx_shl(0, 100) == 0);
    assert(mbe_fx_shl(-5, -1) == -3);
    assert(mbe_fx_shl(MIN16, -100) == -1);
    assert(mbe_fx_shr(MAX16, 100) == 0);
    assert(mbe_fx_shr(MIN16, 100) == -1);
    assert(mbe_fx_shr(-5, 1) == -3);
    assert(mbe_fx_shr(1, -20) == MAX16);

    assert(mbe_fx_L_shl(1, 30) == 0x40000000);
    assert(mbe_fx_L_shl(1, 31) == MAX32);
    assert(mbe_fx_L_shl(-1, 31) == MIN32);
    assert(mbe_fx_L_shl(MAX32, 1000) == MAX32);
    assert(mbe_fx_L_shl(MIN32, 1000) == MIN32);
    assert(mbe_fx_L_shl(-7, -1) == -4);
    assert(mbe_fx_L_shr(MIN32, 1000) == -1);
    assert(mbe_fx_L_shr(MAX32, 31) == 0);
    assert(mbe_fx_L_shr(3, -30) == MAX32);

    /* Rounding shift: ties go up, counts beyond the width give 0. */
    assert(mbe_fx_L_shr_r(1, 1) == 1);
    assert(mbe_fx_L_shr_r(-1, 1) == 0);
    assert(mbe_fx_L_shr_r(-3, 1) == -1);
    assert(mbe_fx_L_shr_r(MAX32, 1) == 0x40000000);
    assert(mbe_fx_L_shr_r(MIN32, 31) == -1);
    assert(mbe_fx_L_shr_r(MAX32, 32) == 0);
    assert(mbe_fx_L_shr_r(MIN32, 40) == 0);
    assert(mbe_fx_L_shr_r(5, 0) == 5);
    assert(mbe_fx_L_shr_r(1, -31) == MAX32);
}

static void
test_products(void) {
    assert(mbe_fx_mult(MIN16, MIN16) == MAX16);
    assert(mbe_fx_mult(MIN16, MAX16) == -MAX16);
    assert(mbe_fx_mult(-1, 1) == -1); /* floor of -2^-15 */
    assert(mbe_fx_mult(16384, 16384) == 8192);
    assert(mbe_fx_mult_r(MIN16, MIN16) == MAX16);
    assert(mbe_fx_mult_r(1, 16384) == 1);  /* 0.5 rounds up */
    assert(mbe_fx_mult_r(-1, 16384) == 0); /* -0.5 rounds up too */
    assert(mbe_fx_mult_r(-3, 16384) == -1);

    assert(mbe_fx_L_mult(MIN16, MIN16) == MAX32);
    assert(mbe_fx_L_mult(MAX16, MAX16) == 2 * (int32_t)MAX16 * MAX16);
    assert(mbe_fx_L_mult(MIN16, MAX16) == -2 * (int32_t)MAX16 * 32768);
    assert(mbe_fx_L_mult0(MIN16, MIN16) == (int32_t)1 << 30);
    assert(mbe_fx_L_mult0(MIN16, MAX16) == -(int32_t)32767 * 32768);

    assert(mbe_fx_L_mac(MAX32, 1, 1) == MAX32);
    assert(mbe_fx_L_msu(MIN32, 1, 1) == MIN32);
    assert(mbe_fx_L_mac(0, MIN16, MIN16) == MAX32);
    assert(mbe_fx_L_msu(0, MIN16, MIN16) == -MAX32);
    assert(mbe_fx_L_mac0(MAX32, MIN16, MIN16) == MAX32);
    assert(mbe_fx_L_msu0(MIN32, MIN16, MIN16) == MIN32);
    assert(mbe_fx_L_mac0(-((int32_t)1 << 30), MIN16, MIN16) == 0);
}

static void
test_extract_and_round(void) {
    assert(mbe_fx_extract_h(0x12345678) == 0x1234);
    assert(mbe_fx_extract_h(-1) == -1);
    assert(mbe_fx_extract_h(MIN32) == MIN16);
    assert(mbe_fx_extract_l(0x12348765) == (int16_t)-30875);
    assert(mbe_fx_extract_l(MIN32) == 0);
    assert(mbe_fx_extract_l(MAX32) == -1);
    assert(mbe_fx_round(0x00008000) == 1);
    assert(mbe_fx_round(0x00007FFF) == 0);
    assert(mbe_fx_round(-0x00008000) == 0);
    assert(mbe_fx_round(MAX32) == MAX16);
    assert(mbe_fx_L_deposit_h(MIN16) == MIN32);
    assert(mbe_fx_L_deposit_h(-1) == -65536);
}

static void
test_norm(void) {
    assert(mbe_fx_norm_l(0) == 0);
    assert(mbe_fx_norm_l(-1) == 31);
    assert(mbe_fx_norm_l(1) == 30);
    assert(mbe_fx_norm_l(MAX32) == 0);
    assert(mbe_fx_norm_l(MIN32) == 0);
    assert(mbe_fx_norm_l(0x40000000) == 0);
    assert(mbe_fx_norm_l(0x3FFFFFFF) == 1);
    assert(mbe_fx_norm_l(-0x40000000) == 1);
    assert(mbe_fx_norm_l(-0x40000001) == 0);
    assert(mbe_fx_norm_l(-2) == 30);
}

static void
test_double_precision(void) {
    static const int32_t values[] = {0, 1, -1, 65535, -65536, 0x12345678, -0x12345678, MAX32, MIN32};
    size_t i;
    for (i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        int16_t hi;
        int16_t lo;
        int64_t back;
        mbe_fx_L_extract(values[i], &hi, &lo);
        /* (hi << 16) + (lo << 1) restores the value to its least significant bit. */
        back = mbe_fx_lsl64(hi, 16) + mbe_fx_lsl64(lo, 1);
        assert(lo >= 0);
        assert(back == mbe_fx_lsl64(mbe_fx_asr64(values[i], 1), 1));
    }
    /* 0.5 (Q31) times 0.5 (Q15) is 0.25 (Q31). */
    {
        int16_t hi;
        int16_t lo;
        mbe_fx_L_extract(0x40000000, &hi, &lo);
        assert(mbe_fx_mpy_32_16(hi, lo, 16384) == 0x20000000);
    }
}

int
main(void) {
    test_saturation();
    test_shifts();
    test_products();
    test_extract_and_round();
    test_norm();
    test_double_precision();
    return 0;
}
