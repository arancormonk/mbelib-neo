// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief API-level coverage for complete AMBE/IMBE frame processing paths.
 */

#include <assert.h>
#include <math.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

static void
assert_float_pcm_sane(const float* pcm) {
    for (int i = 0; i < 160; ++i) {
        assert(isfinite(pcm[i]));
        assert(pcm[i] > -20000.0f);
        assert(pcm[i] < 20000.0f);
    }
}

static void
assert_status_terminated(const char* status, size_t size) {
    int terminated = 0;
    for (size_t i = 0; i < size; ++i) {
        if (status[i] == '\0') {
            terminated = 1;
            break;
        }
    }
    assert(terminated);
}

static void
assert_result_total(const mbe_process_result* result, int ret) {
    assert(ret == result->total_errors);
    assert(result->total_errors == result->c0_errors + result->protected_errors);
}

static void
init_params(mbe_parms* cur, mbe_parms* prev, mbe_parms* enh) {
    mbe_initMbeParms(cur, prev, enh);
    mbe_setThreadRngSeed(0x12345678u);
}

static void
copy_bits(const char* src, char* dst, size_t count) {
    memcpy(dst, src, count);
}

static void
soft_from_bits(const char* hard, mbe_soft_bit* soft, size_t count) {
    mbe_softBitsFromHard(hard, soft, count, 255u);
}

static void
test_ambe2400_frame_paths(void) {
    char frame[4][24] = {{0}};
    char mutable_frame[4][24];
    char hard_d[49] = {0};
    char soft_d[49] = {0};
    mbe_soft_bit soft_frame[4][24];
    mbe_process_result hard_result;
    mbe_process_result soft_result;
    float out_f[160];
    short out_s[160];
    char status[256];
    int errs = 0;
    int errs2 = 0;
    mbe_parms cur, prev, enh;

    frame[0][0] = 1;
    frame[1][7] = 1;
    frame[2][3] = 1;
    frame[3][11] = 1;
    const char (*const_frame)[24] = (const char (*)[24])frame;
    soft_from_bits(&frame[0][0], &soft_frame[0][0], (size_t)4 * 24u);
    const mbe_soft_bit(*const_soft_frame)[24] = (const mbe_soft_bit(*)[24])soft_frame;

    int hard_ret = mbe_decodeAmbe3600x2400Frame(const_frame, hard_d, &hard_result);
    int soft_ret = mbe_decodeAmbe3600x2400SoftFrame(const_soft_frame, soft_d, &soft_result);
    assert_result_total(&hard_result, hard_ret);
    assert_result_total(&soft_result, soft_ret);
    assert((hard_result.flags & MBE_PROCESS_FLAG_C0_VALID) != 0u);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert((soft_result.flags & MBE_PROCESS_FLAG_C0_VALID) != 0u);

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    mbe_processAmbe2400Dataf(out_f, &hard_result.c0_errors, &hard_result.total_errors, status, hard_d, &cur, &prev,
                             &enh, 8);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    mbe_processAmbe2400Data(out_s, &hard_result.c0_errors, &hard_result.total_errors, status, hard_d, &cur, &prev, &enh,
                            8);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    mbe_process_result data_result = hard_result;
    int ret = mbe_processAmbe2400DatafV2(out_f, &data_result, hard_d, &cur, &prev, &enh, 8);
    assert_result_total(&data_result, ret);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    data_result = hard_result;
    ret = mbe_processAmbe2400DataV2(out_s, &data_result, hard_d, &cur, &prev, &enh, 8);
    assert_result_total(&data_result, ret);

    init_params(&cur, &prev, &enh);
    copy_bits(&frame[0][0], &mutable_frame[0][0], sizeof(frame));
    memset(status, 0, sizeof(status));
    mbe_processAmbe3600x2400Framef(out_f, &errs, &errs2, status, mutable_frame, hard_d, &cur, &prev, &enh, 8);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));
    assert(errs2 >= errs);

    init_params(&cur, &prev, &enh);
    copy_bits(&frame[0][0], &mutable_frame[0][0], sizeof(frame));
    memset(status, 0, sizeof(status));
    mbe_processAmbe3600x2400Frame(out_s, &errs, &errs2, status, mutable_frame, hard_d, &cur, &prev, &enh, 8);
    assert_status_terminated(status, sizeof(status));
    assert(errs2 >= errs);

    init_params(&cur, &prev, &enh);
    ret = mbe_processAmbe3600x2400SoftFramef(out_f, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh, 8);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    ret = mbe_processAmbe3600x2400SoftFrame(out_s, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh, 8);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
}

static void
test_imbe7200_frame_paths(void) {
    char frame[8][23] = {{0}};
    char mutable_frame[8][23];
    char hard_d[88] = {0};
    char soft_d[88] = {0};
    mbe_soft_bit soft_frame[8][23];
    mbe_process_result hard_result;
    mbe_process_result soft_result;
    float out_f[160];
    short out_s[160];
    char status[256];
    int errs = 0;
    int errs2 = 0;
    mbe_parms cur, prev, enh;

    frame[0][22] = 1;
    frame[2][5] = 1;
    frame[5][9] = 1;
    frame[7][2] = 1;
    const char (*const_frame)[23] = (const char (*)[23])frame;
    soft_from_bits(&frame[0][0], &soft_frame[0][0], (size_t)8 * 23u);
    const mbe_soft_bit(*const_soft_frame)[23] = (const mbe_soft_bit(*)[23])soft_frame;

    int hard_ret = mbe_decodeImbe7200x4400Frame(const_frame, hard_d, &hard_result);
    int soft_ret = mbe_decodeImbe7200x4400SoftFrame(const_soft_frame, soft_d, &soft_result);
    assert_result_total(&hard_result, hard_ret);
    assert_result_total(&soft_result, soft_ret);
    assert((hard_result.flags & MBE_PROCESS_FLAG_C0_VALID) != 0u);
    assert((hard_result.flags & MBE_PROCESS_FLAG_C4_VALID) != 0u);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert((soft_result.flags & MBE_PROCESS_FLAG_C4_VALID) != 0u);

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    mbe_processImbe4400Dataf(out_f, &hard_result.c0_errors, &hard_result.total_errors, status, hard_d, &cur, &prev,
                             &enh, 8);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    mbe_processImbe4400Data(out_s, &hard_result.c0_errors, &hard_result.total_errors, status, hard_d, &cur, &prev, &enh,
                            8);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    mbe_process_result data_result = hard_result;
    int ret = mbe_processImbe4400DatafV2(out_f, &data_result, hard_d, &cur, &prev, &enh, 8);
    assert_result_total(&data_result, ret);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    data_result = hard_result;
    ret = mbe_processImbe4400DataV2(out_s, &data_result, hard_d, &cur, &prev, &enh, 8);
    assert_result_total(&data_result, ret);

    init_params(&cur, &prev, &enh);
    copy_bits(&frame[0][0], &mutable_frame[0][0], sizeof(frame));
    memset(status, 0, sizeof(status));
    mbe_processImbe7200x4400Framef(out_f, &errs, &errs2, status, mutable_frame, hard_d, &cur, &prev, &enh, 8);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));
    assert(errs2 >= errs);

    init_params(&cur, &prev, &enh);
    copy_bits(&frame[0][0], &mutable_frame[0][0], sizeof(frame));
    memset(status, 0, sizeof(status));
    mbe_processImbe7200x4400Frame(out_s, &errs, &errs2, status, mutable_frame, hard_d, &cur, &prev, &enh, 8);
    assert_status_terminated(status, sizeof(status));
    assert(errs2 >= errs);

    init_params(&cur, &prev, &enh);
    ret = mbe_processImbe7200x4400SoftFramef(out_f, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh, 8);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    ret = mbe_processImbe7200x4400SoftFrame(out_s, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh, 8);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
}

static void
test_imbe7100_frame_paths(void) {
    char frame[7][24] = {{0}};
    char mutable_frame[7][24];
    char hard_d[88] = {0};
    char soft_d[88] = {0};
    mbe_soft_bit soft_frame[7][24];
    mbe_process_result hard_result;
    mbe_process_result soft_result;
    float out_f[160];
    short out_s[160];
    char status[256];
    int errs = 0;
    int errs2 = 0;
    mbe_parms cur, prev, enh;

    frame[0][18] = 1;
    frame[1][4] = 1;
    frame[4][6] = 1;
    frame[6][2] = 1;
    const char (*const_frame)[24] = (const char (*)[24])frame;
    soft_from_bits(&frame[0][0], &soft_frame[0][0], (size_t)7 * 24u);
    const mbe_soft_bit(*const_soft_frame)[24] = (const mbe_soft_bit(*)[24])soft_frame;

    int hard_ret = mbe_decodeImbe7100x4400Frame(const_frame, hard_d, &hard_result);
    int soft_ret = mbe_decodeImbe7100x4400SoftFrame(const_soft_frame, soft_d, &soft_result);
    assert_result_total(&hard_result, hard_ret);
    assert_result_total(&soft_result, soft_ret);
    assert((hard_result.flags & MBE_PROCESS_FLAG_C0_VALID) != 0u);
    assert((hard_result.flags & MBE_PROCESS_FLAG_C4_VALID) != 0u);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert((soft_result.flags & MBE_PROCESS_FLAG_C4_VALID) != 0u);

    init_params(&cur, &prev, &enh);
    copy_bits(&frame[0][0], &mutable_frame[0][0], sizeof(frame));
    memset(status, 0, sizeof(status));
    mbe_processImbe7100x4400Framef(out_f, &errs, &errs2, status, mutable_frame, hard_d, &cur, &prev, &enh, 8);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));
    assert(errs2 >= errs);

    init_params(&cur, &prev, &enh);
    copy_bits(&frame[0][0], &mutable_frame[0][0], sizeof(frame));
    memset(status, 0, sizeof(status));
    mbe_processImbe7100x4400Frame(out_s, &errs, &errs2, status, mutable_frame, hard_d, &cur, &prev, &enh, 8);
    assert_status_terminated(status, sizeof(status));
    assert(errs2 >= errs);

    init_params(&cur, &prev, &enh);
    int ret = mbe_processImbe7100x4400SoftFramef(out_f, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh, 8);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    ret = mbe_processImbe7100x4400SoftFrame(out_s, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh, 8);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
}

int
main(void) {
    test_ambe2400_frame_paths();
    test_imbe7200_frame_paths();
    test_imbe7100_frame_paths();
    return 0;
}
