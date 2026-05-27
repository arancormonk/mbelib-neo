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
set_imbe7200_b0(char imbe_d[88], int b0) {
    for (int i = 0; i < 8; ++i) {
        int bit = (b0 >> (7 - i)) & 1;
        if (i < 6) {
            imbe_d[i] = (char)bit;
        } else if (i == 6) {
            imbe_d[85] = (char)bit;
        } else {
            imbe_d[86] = (char)bit;
        }
    }
}

static void
set_valid_ambe_tone_bits(char ambe_d[49]) {
    memset(ambe_d, 0, 49u);
    ambe_d[6] = 1;
    ambe_d[7] = 1;
    ambe_d[8] = 1;
    ambe_d[9] = 1;
    ambe_d[10] = 1;
    ambe_d[11] = 1;
    ambe_d[17] = 1;
    ambe_d[19] = 1;
    ambe_d[44] = 1;
}

static void
fill_float_buffer(float* out, float value) {
    for (int i = 0; i < 160; ++i) {
        out[i] = value;
    }
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
test_imbe_previous_harmonic_count_is_clamped_for_processing(void) {
    static const int invalid_l_values[] = {-1, 0, 57, 1000};

    for (size_t i = 0; i < sizeof(invalid_l_values) / sizeof(invalid_l_values[0]); ++i) {
        char params[88] = {0};
        float out[160];
        mbe_process_result result;
        mbe_parms cur = {0}, prev = {0}, enh = {0};

        set_imbe7200_b0(params, 206);
        init_params(&cur, &prev, &enh);
        mbe_initProcessResult(&result);
        prev.L = invalid_l_values[i];
        prev.Ml[1] = 1.5f;
        prev.log2Ml[1] = 0.25f;
        prev.Ml[56] = 2.5f;
        prev.log2Ml[56] = 0.5f;
        fill_float_buffer(out, 123.0f);

        assert(mbe_processImbe4400Dataf(out, &result, params, &cur, &prev, &enh, 8) >= 0);
        assert(cur.L == 56);
    }
}

static void
test_imbe_previous_harmonic_count_upper_endpoint_is_safe(void) {
    char params[88] = {0};
    float out[160];
    mbe_process_result result;
    mbe_parms cur = {0}, prev = {0}, enh = {0};

    set_imbe7200_b0(params, 206);
    init_params(&cur, &prev, &enh);
    mbe_initProcessResult(&result);
    prev.L = 56;
    prev.Ml[56] = 3.0f;
    prev.log2Ml[56] = 1.0f;
    fill_float_buffer(out, 123.0f);

    assert(mbe_processImbe4400Dataf(out, &result, params, &cur, &prev, &enh, 8) >= 0);
    assert(cur.L == 56);
}

static void
stamp_params(mbe_parms* mp, int marker) {
    mp->w0 = (float)marker + 0.25f;
    mp->L = marker;
    mp->K = marker + 1;
    mp->Vl[1] = marker + 2;
    mp->Ml[1] = (float)marker + 0.5f;
    mp->log2Ml[1] = (float)marker + 0.75f;
    mp->PHIl[1] = (float)marker + 1.0f;
    mp->PSIl[1] = (float)marker + 1.25f;
    mp->gamma = (float)marker + 1.5f;
    mp->un = marker + 3;
    mp->repeat = marker + 4;
    mp->swn = marker + 5;
    mp->localEnergy = (float)marker + 1.75f;
    mp->amplitudeThreshold = marker + 6;
    mp->errorRate = (float)marker + 2.0f;
    mp->errorCountTotal = marker + 7;
    mp->errorCount4 = marker + 8;
    mp->repeatCount = marker + 9;
    mp->mutingThreshold = (float)marker + 2.25f;
    mp->previousUw[0] = (float)marker + 2.5f;
    mp->noiseSeed = (float)marker + 2.75f;
    mp->noiseOverlap[0] = (float)marker + 3.0f;
}

static void
assert_params_stamp(const mbe_parms* mp, int marker) {
    assert(mp->w0 == (float)marker + 0.25f);
    assert(mp->L == marker);
    assert(mp->K == marker + 1);
    assert(mp->Vl[1] == marker + 2);
    assert(mp->Ml[1] == (float)marker + 0.5f);
    assert(mp->log2Ml[1] == (float)marker + 0.75f);
    assert(mp->PHIl[1] == (float)marker + 1.0f);
    assert(mp->PSIl[1] == (float)marker + 1.25f);
    assert(mp->gamma == (float)marker + 1.5f);
    assert(mp->un == marker + 3);
    assert(mp->repeat == marker + 4);
    assert(mp->swn == marker + 5);
    assert(mp->localEnergy == (float)marker + 1.75f);
    assert(mp->amplitudeThreshold == marker + 6);
    assert(mp->errorRate == (float)marker + 2.0f);
    assert(mp->errorCountTotal == marker + 7);
    assert(mp->errorCount4 == marker + 8);
    assert(mp->repeatCount == marker + 9);
    assert(mp->mutingThreshold == (float)marker + 2.25f);
    assert(mp->previousUw[0] == (float)marker + 2.5f);
    assert(mp->noiseSeed == (float)marker + 2.75f);
    assert(mp->noiseOverlap[0] == (float)marker + 3.0f);
}

static void
test_null_parameter_helpers_are_noops(void) {
    mbe_parms cur = {0}, prev = {0}, enh = {0}, src = {0}, dst = {0};

    init_params(&cur, &prev, &enh);
    stamp_params(&cur, 10);
    stamp_params(&prev, 20);
    stamp_params(&enh, 30);

    mbe_initMbeParms(NULL, &prev, &enh);
    assert_params_stamp(&prev, 20);
    assert_params_stamp(&enh, 30);

    mbe_initMbeParms(&cur, NULL, &enh);
    assert_params_stamp(&cur, 10);
    assert_params_stamp(&enh, 30);

    mbe_initMbeParms(&cur, &prev, NULL);
    assert_params_stamp(&cur, 10);
    assert_params_stamp(&prev, 20);
    mbe_initMbeParms(NULL, NULL, NULL);

    stamp_params(&src, 40);
    stamp_params(&dst, 50);
    mbe_moveMbeParms(NULL, &dst);
    assert_params_stamp(&dst, 50);
    mbe_moveMbeParms(&src, NULL);
    mbe_moveMbeParms(NULL, NULL);

    mbe_moveMbeParms(&src, &dst);
    assert_params_stamp(&dst, 40);

    stamp_params(&dst, 60);
    mbe_useLastMbeParms(&dst, NULL);
    assert_params_stamp(&dst, 60);
    mbe_useLastMbeParms(NULL, &src);
    mbe_useLastMbeParms(NULL, NULL);
}

static void
test_tone_helpers_silence_invalid_inputs(void) {
    char ambe_d[49];
    float out[160];
    mbe_parms cur = {0}, prev = {0}, enh = {0};

    init_params(&cur, &prev, &enh);
    set_valid_ambe_tone_bits(ambe_d);

    fill_float_buffer(out, 123.0f);
    mbe_synthesizeTonef(out, NULL, &cur);
    assert_float_buffer_unchanged(out, 0.0f);

    set_valid_ambe_tone_bits(ambe_d);
    ambe_d[0] = 2;
    fill_float_buffer(out, 123.0f);
    mbe_synthesizeTonef(out, ambe_d, &cur);
    assert_float_buffer_unchanged(out, 0.0f);

    set_valid_ambe_tone_bits(ambe_d);
    fill_float_buffer(out, 123.0f);
    mbe_synthesizeTonef(out, ambe_d, NULL);
    assert_float_buffer_unchanged(out, 0.0f);

    fill_float_buffer(out, 123.0f);
    mbe_synthesizeTonefdstar(out, ambe_d, NULL, 5);
    assert_float_buffer_unchanged(out, 0.0f);
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
    test_imbe_previous_harmonic_count_is_clamped_for_processing();
    test_imbe_previous_harmonic_count_upper_endpoint_is_safe();
    test_null_parameter_helpers_are_noops();
    test_tone_helpers_silence_invalid_inputs();
    test_low_level_public_bit_paths_validate();
    return 0;
}
