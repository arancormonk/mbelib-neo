// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Strict public input validation tests for codec bit arrays.
 */

#include <assert.h>
#include <limits.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

static void
init_params(mbe_parms* cur, mbe_parms* prev, mbe_parms* enh) {
    mbe_initMbeParms(cur, prev, enh);
}

static void
assert_float_buffer_unchanged(const float* out, float value) {
    for (int i = 0; i < 160; ++i) {
        assert(out[i] == value);
    }
}

static void
assert_short_buffer_unchanged(const short* out, short value) {
    for (int i = 0; i < 160; ++i) {
        assert(out[i] == value);
    }
}

static void
test_invalid_hard_frame_rejected_before_decode(void) {
    char frame[4][24] = {{0}};
    char params[49];
    mbe_process_result result;

    memset(params, 7, sizeof(params));
    frame[0][0] = 2;

    assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])frame, params, &result) == MBE_STATUS_INVALID_BITS);
    assert(result.total_errors == 0);
    for (size_t i = 0; i < sizeof(params); ++i) {
        assert(params[i] == 7);
    }
}

static void
test_invalid_soft_frame_rejected_before_decode(void) {
    mbe_soft_bit frame[8][23] = {{{0, 255}}};
    char params[88];
    mbe_process_result result;

    memset(params, 9, sizeof(params));
    frame[0][0].bit = 2;

    assert(mbe_decodeImbe7200x4400SoftFrame((const mbe_soft_bit(*)[23])frame, params, &result)
           == MBE_STATUS_INVALID_BITS);
    assert(result.total_errors == 0);
    for (size_t i = 0; i < sizeof(params); ++i) {
        assert(params[i] == 9);
    }
}

static void
test_invalid_parameter_bits_rejected_before_synthesis(void) {
    char params[49] = {0};
    float out[160];
    mbe_process_result result;
    mbe_parms cur = {0}, prev = {0}, enh = {0};

    for (int i = 0; i < 160; ++i) {
        out[i] = 123.0f;
    }
    init_params(&cur, &prev, &enh);
    mbe_initProcessResult(&result);
    params[8] = 3;

    assert(mbe_processAmbe2400Dataf(out, &result, params, &cur, &prev, &enh, 8) == MBE_STATUS_INVALID_BITS);
    assert_float_buffer_unchanged(out, 123.0f);
}

static int
process_imbe4400_with_result_context(const mbe_process_result* context) {
    char params[88] = {0};
    float out[160] = {0};
    mbe_process_result result = *context;
    mbe_parms cur = {0}, prev = {0}, enh = {0};

    init_params(&cur, &prev, &enh);

    return mbe_processImbe4400Dataf(out, &result, params, &cur, &prev, &enh, 8);
}

static void
expect_invalid_result_context(mbe_process_result result) {
    assert(process_imbe4400_with_result_context(&result) == MBE_STATUS_INVALID_ARGUMENT);
}

static void
expect_valid_result_context(mbe_process_result result) {
    assert(process_imbe4400_with_result_context(&result) >= 0);
}

static void
test_invalid_result_context_rejected(void) {
    mbe_process_result result;

    mbe_initProcessResult(&result);
    result.total_errors = -1;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.total_errors = 185;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.total_errors = INT_MAX;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.c0_errors = INT_MAX;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.protected_errors = 185;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.c4_errors = 185;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.c0_errors = 100;
    result.protected_errors = 85;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.flags = MBE_PROCESS_FLAG_C0_VALID;
    result.c0_errors = 7;
    result.total_errors = 6;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.flags = MBE_PROCESS_FLAG_C4_VALID;
    result.c4_errors = 7;
    result.total_errors = 6;
    expect_invalid_result_context(result);

    mbe_initProcessResult(&result);
    result.c0_errors = 100;
    result.protected_errors = 84;
    expect_valid_result_context(result);

    mbe_initProcessResult(&result);
    result.total_errors = 184;
    expect_valid_result_context(result);
}

static void
test_invalid_harmonic_counts_rejected_by_raw_helpers(void) {
    static const int invalid_l_values[] = {-1, 0, 57, 1000};

    for (size_t i = 0; i < sizeof(invalid_l_values) / sizeof(invalid_l_values[0]); ++i) {
        const int invalid_l = invalid_l_values[i];
        float out[160];
        short out_s[160];
        mbe_parms cur = {0}, prev = {0}, enh = {0};

        init_params(&cur, &prev, &enh);
        for (int j = 0; j < 160; ++j) {
            out[j] = 123.0f;
        }
        cur.L = invalid_l;
        mbe_synthesizeSpeechf(out, &cur, &prev, 8);
        assert_float_buffer_unchanged(out, 0.0f);

        init_params(&cur, &prev, &enh);
        for (int j = 0; j < 160; ++j) {
            out[j] = 123.0f;
        }
        prev.L = invalid_l;
        mbe_synthesizeSpeechf(out, &cur, &prev, 8);
        assert_float_buffer_unchanged(out, 0.0f);

        init_params(&cur, &prev, &enh);
        for (int j = 0; j < 160; ++j) {
            out_s[j] = 123;
        }
        cur.L = invalid_l;
        mbe_synthesizeSpeech(out_s, &cur, &prev, 8);
        assert_short_buffer_unchanged(out_s, 0);

        init_params(&cur, &prev, &enh);
        cur.L = invalid_l;
        cur.Ml[1] = 3.0f;
        cur.Ml[56] = 5.0f;
        cur.gamma = 9.0f;
        mbe_spectralAmpEnhance(&cur);
        assert(cur.L == invalid_l);
        assert(cur.Ml[1] == 3.0f);
        assert(cur.Ml[56] == 5.0f);
        assert(cur.gamma == 9.0f);

        init_params(&cur, &prev, &enh);
        cur.L = invalid_l;
        cur.errorRate = 1.0f;
        cur.errorCountTotal = 10;
        cur.localEnergy = 123.0f;
        cur.amplitudeThreshold = 321;
        mbe_applyAdaptiveSmoothing(&cur, &prev);
        assert(cur.L == invalid_l);
        assert(cur.localEnergy == 123.0f);
        assert(cur.amplitudeThreshold == 321);

        init_params(&cur, &prev, &enh);
        prev.L = invalid_l;
        cur.errorRate = 1.0f;
        cur.errorCountTotal = 10;
        cur.localEnergy = 123.0f;
        cur.amplitudeThreshold = 321;
        mbe_applyAdaptiveSmoothing(&cur, &prev);
        assert(prev.L == invalid_l);
        assert(cur.localEnergy == 123.0f);
        assert(cur.amplitudeThreshold == 321);
    }
}

static void
test_low_level_public_bit_paths_validate(void) {
    char ambe_frame[4][24] = {{0}};
    char imbe7100_frame[7][24] = {{0}};
    char imbe_params[88] = {0};
    char hard_bits[2] = {0, 2};
    mbe_soft_bit soft_bits[2];

    assert(mbe_checkGolayBlock(NULL) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_softBitsFromHard(hard_bits, soft_bits, 2u, 255u) == MBE_STATUS_INVALID_BITS);
    assert(mbe_softBitsFromLlr(NULL, soft_bits, 2u) == MBE_STATUS_INVALID_ARGUMENT);

    ambe_frame[1][1] = 4;
    assert(mbe_eccAmbe3600x2450C0(ambe_frame) == MBE_STATUS_INVALID_BITS);
    assert(mbe_demodulateAmbe3600x2450Data(ambe_frame) == MBE_STATUS_INVALID_BITS);

    imbe7100_frame[2][3] = 5;
    assert(mbe_eccImbe7100x4400Data(imbe7100_frame, imbe_params) == MBE_STATUS_INVALID_BITS);
    assert(mbe_demodulateImbe7100x4400Data(imbe7100_frame) == MBE_STATUS_INVALID_BITS);

    imbe_params[10] = 2;
    assert(mbe_convertImbe7100to7200(imbe_params) == MBE_STATUS_INVALID_BITS);
}

int
main(void) {
    test_invalid_hard_frame_rejected_before_decode();
    test_invalid_soft_frame_rejected_before_decode();
    test_invalid_parameter_bits_rejected_before_synthesis();
    test_invalid_result_context_rejected();
    test_invalid_harmonic_counts_rejected_by_raw_helpers();
    test_low_level_public_bit_paths_validate();
    return 0;
}
