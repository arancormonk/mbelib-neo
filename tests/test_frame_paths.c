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
#include <stdint.h>
#include <string.h>

#include "mbe_quality_frames.h"
#include "mbelib-neo/mbelib.h"

static void
assert_float_pcm_sane(const float* pcm) {
    for (int i = 0; i < 160; ++i) {
        assert(isfinite(pcm[i]));
        assert(pcm[i] > -20000.0f);
        assert(pcm[i] < 20000.0f);
    }
}

static int
float_array_equal_exact(const float* a, const float* b, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        uint32_t a_bits, b_bits;
        memcpy(&a_bits, &a[i], sizeof(a_bits));
        memcpy(&b_bits, &b[i], sizeof(b_bits));
        if (a_bits != b_bits) {
            return 0;
        }
    }
    return 1;
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
soft_from_bits(const char* hard, mbe_soft_bit* soft, size_t count) {
    assert(mbe_softBitsFromHard(hard, soft, count, 255u) == 0);
}

static void
test_ambe2400_frame_paths(void) {
    char frame[4][24] = {{0}};
    char hard_d[49] = {0};
    char soft_d[49] = {0};
    mbe_soft_bit soft_frame[4][24];
    mbe_process_result hard_result;
    mbe_process_result soft_result;
    mbe_process_result data_result;
    float out_f[160];
    short out_s[160];
    char status[256];
    int ret;
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
    data_result = hard_result;
    ret = mbe_processAmbe2400Dataf(out_f, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &data_result);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    data_result = hard_result;
    ret = mbe_processAmbe2400Data(out_s, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &data_result);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    data_result = hard_result;
    ret = mbe_processAmbe2400Dataf(out_f, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    data_result = hard_result;
    ret = mbe_processAmbe2400Data(out_s, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    ret = mbe_processAmbe3600x2400Framef(out_f, &hard_result, const_frame, hard_d, &cur, &prev, &enh);
    assert_result_total(&hard_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &hard_result);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    ret = mbe_processAmbe3600x2400Frame(out_s, &hard_result, const_frame, hard_d, &cur, &prev, &enh);
    assert_result_total(&hard_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &hard_result);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    ret = mbe_processAmbe3600x2400SoftFramef(out_f, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    ret = mbe_processAmbe3600x2400SoftFrame(out_s, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
}

static void
test_imbe7200_frame_paths(void) {
    char frame[8][23] = {{0}};
    char hard_d[88] = {0};
    char soft_d[88] = {0};
    mbe_soft_bit soft_frame[8][23];
    mbe_process_result hard_result;
    mbe_process_result soft_result;
    mbe_process_result data_result;
    float out_f[160];
    short out_s[160];
    char status[256];
    int ret;
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
    data_result = hard_result;
    ret = mbe_processImbe4400Dataf(out_f, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &data_result);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    data_result = hard_result;
    ret = mbe_processImbe4400Data(out_s, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &data_result);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    data_result = hard_result;
    ret = mbe_processImbe4400Dataf(out_f, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    data_result = hard_result;
    ret = mbe_processImbe4400Data(out_s, &data_result, hard_d, &cur, &prev, &enh);
    assert_result_total(&data_result, ret);

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    ret = mbe_processImbe7200x4400Framef(out_f, &hard_result, const_frame, hard_d, &cur, &prev, &enh);
    assert_result_total(&hard_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &hard_result);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    ret = mbe_processImbe7200x4400Frame(out_s, &hard_result, const_frame, hard_d, &cur, &prev, &enh);
    assert_result_total(&hard_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &hard_result);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    ret = mbe_processImbe7200x4400SoftFramef(out_f, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    ret = mbe_processImbe7200x4400SoftFrame(out_s, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
}

static void
test_imbe7100_frame_paths(void) {
    char frame[7][24] = {{0}};
    char hard_d[88] = {0};
    char soft_d[88] = {0};
    mbe_soft_bit soft_frame[7][24];
    mbe_process_result hard_result;
    mbe_process_result soft_result;
    float out_f[160];
    short out_s[160];
    char status[256];
    int ret;
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
    memset(status, 0, sizeof(status));
    ret = mbe_processImbe7100x4400Framef(out_f, &hard_result, const_frame, hard_d, &cur, &prev, &enh);
    assert_result_total(&hard_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &hard_result);
    assert_float_pcm_sane(out_f);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    memset(status, 0, sizeof(status));
    ret = mbe_processImbe7100x4400Frame(out_s, &hard_result, const_frame, hard_d, &cur, &prev, &enh);
    assert_result_total(&hard_result, ret);
    mbe_formatProcessResult(status, sizeof(status), &hard_result);
    assert_status_terminated(status, sizeof(status));

    init_params(&cur, &prev, &enh);
    ret = mbe_processImbe7100x4400SoftFramef(out_f, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
    assert_float_pcm_sane(out_f);

    init_params(&cur, &prev, &enh);
    ret = mbe_processImbe7100x4400SoftFrame(out_s, &soft_result, const_soft_frame, soft_d, &cur, &prev, &enh);
    assert_result_total(&soft_result, ret);
    assert((soft_result.flags & MBE_PROCESS_FLAG_SOFT_INPUT) != 0u);
}

enum fixture_mode { IMBE7200, IMBE7100, AMBE2450, AMBE2400 };

static const struct {
    const char* name;
    size_t data_count;
    size_t frame_count;
    int columns;
    int ecc_rows;
} fixture_modes[] = {
    {"imbe7200", 88, 184, 23, 7},
    {"imbe7100", 88, 168, 24, 6},
    {"ambe2450", 49, 96, 24, 2},
    {"ambe2400", 49, 96, 24, 2},
};

static void
tagged_data(char* data, size_t count, unsigned tag) {
    /* Nonperiodic payloads expose misplaced bits in the variable-K permutation. */
    unsigned state = 0x9e3779b9u ^ tag;
    for (size_t i = 0; i < count; ++i) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        data[i] = (char)(state & 1u);
    }
}

static int
decode_fixture(enum fixture_mode mode, const char* frame, char* data, mbe_process_result* result, int soft) {
    mbe_soft_bit reliability[184];
    if (soft) {
        soft_from_bits(frame, reliability, fixture_modes[mode].frame_count);
    }
    switch (mode) {
        case IMBE7200:
            return soft ? mbe_decodeImbe7200x4400SoftFrame((const mbe_soft_bit(*)[23])reliability, data, result)
                        : mbe_decodeImbe7200x4400Frame((const char (*)[23])frame, data, result);
        case IMBE7100:
            return soft ? mbe_decodeImbe7100x4400SoftFrame((const mbe_soft_bit(*)[24])reliability, data, result)
                        : mbe_decodeImbe7100x4400Frame((const char (*)[24])frame, data, result);
        case AMBE2450:
            return soft ? mbe_decodeAmbe3600x2450SoftFrame((const mbe_soft_bit(*)[24])reliability, data, result)
                        : mbe_decodeAmbe3600x2450Frame((const char (*)[24])frame, data, result);
        case AMBE2400:
            return soft ? mbe_decodeAmbe3600x2400SoftFrame((const mbe_soft_bit(*)[24])reliability, data, result)
                        : mbe_decodeAmbe3600x2400Frame((const char (*)[24])frame, data, result);
    }
    assert(0);
    return MBE_STATUS_INVALID_ARGUMENT;
}

static void
assert_same_result(const mbe_process_result* hard, const mbe_process_result* other, int soft) {
    assert(hard->c0_errors == other->c0_errors);
    assert(hard->protected_errors == other->protected_errors);
    assert(hard->c4_errors == other->c4_errors);
    assert(hard->total_errors == other->total_errors);
    assert(other->flags == (hard->flags | (soft ? MBE_PROCESS_FLAG_SOFT_INPUT : 0u)));
}

static void
assert_decoded_fixture(enum fixture_mode mode, const char* frame, const char* expected, int c0, int protected, int c4) {
    char data[88];
    mbe_process_result hard, soft;
    int ret = decode_fixture(mode, frame, data, &hard, 0);
    assert_result_total(&hard, ret);
    assert(ret == c0 + protected);
    assert(hard.c0_errors == c0);
    assert(hard.protected_errors == protected);
    assert(hard.c4_errors == c4);
    assert(hard.flags == (MBE_PROCESS_FLAG_C0_VALID | (mode <= IMBE7100 ? MBE_PROCESS_FLAG_C4_VALID : 0u)));
    assert(memcmp(data, expected, fixture_modes[mode].data_count) == 0);
    ret = decode_fixture(mode, frame, data, &soft, 1);
    assert_result_total(&soft, ret);
    assert_same_result(&hard, &soft, 1);
    assert(memcmp(data, expected, fixture_modes[mode].data_count) == 0);
}

static void
test_fixture_validation(void) {
    char data[89] = {0};
    char frame[184];
    char saved[184];
    memset(frame, 0x5a, sizeof(frame));
    memcpy(saved, frame, sizeof(saved));
    assert(mbe_quality_frame_from_data(NULL, data, 88, frame, sizeof(frame)) == MBE_STATUS_INVALID_ARGUMENT);
    assert(memcmp(frame, saved, sizeof(frame)) == 0);
    assert(mbe_quality_frame_from_data("IMBE7200", data, 88, frame, sizeof(frame)) == MBE_STATUS_INVALID_ARGUMENT);
    assert(memcmp(frame, saved, sizeof(frame)) == 0);
    assert(mbe_quality_frame_from_data("", data, 88, frame, sizeof(frame)) == MBE_STATUS_INVALID_ARGUMENT);
    assert(memcmp(frame, saved, sizeof(frame)) == 0);
    for (int mode = IMBE7200; mode <= AMBE2400; ++mode) {
        const char* codec = fixture_modes[mode].name;
        size_t count = fixture_modes[mode].data_count;
        size_t capacity = fixture_modes[mode].frame_count;
        assert(mbe_quality_frame_from_data(codec, NULL, count, frame, capacity) == MBE_STATUS_INVALID_ARGUMENT);
        assert(memcmp(frame, saved, sizeof(frame)) == 0);
        assert(mbe_quality_frame_from_data(codec, data, count, NULL, capacity) == MBE_STATUS_INVALID_ARGUMENT);
        assert(mbe_quality_frame_from_data(codec, data, 0, frame, capacity) == MBE_STATUS_INVALID_ARGUMENT);
        assert(memcmp(frame, saved, sizeof(frame)) == 0);
        assert(mbe_quality_frame_from_data(codec, data, count - 1, frame, capacity) == MBE_STATUS_INVALID_ARGUMENT);
        assert(memcmp(frame, saved, sizeof(frame)) == 0);
        assert(mbe_quality_frame_from_data(codec, data, count + 1, frame, capacity) == MBE_STATUS_INVALID_ARGUMENT);
        assert(memcmp(frame, saved, sizeof(frame)) == 0);
        assert(mbe_quality_frame_from_data(codec, data, count, frame, capacity - 1) == MBE_STATUS_INVALID_ARGUMENT);
        assert(memcmp(frame, saved, sizeof(frame)) == 0);
        for (size_t i = 0; i < count; ++i) {
            data[i] = (i & 1u) ? (char)-1 : 2;
            assert(mbe_quality_frame_from_data(codec, data, count, frame, capacity) == MBE_STATUS_INVALID_BITS);
            assert(memcmp(frame, saved, sizeof(frame)) == 0);
            data[i] = 0;
        }
    }
}

static void
test_fixture_ecc(void) {
    char data[88], frame[184], damaged[184];
    for (int mode = IMBE7200; mode <= AMBE2400; ++mode) {
        for (unsigned tag = 0; tag < 3; ++tag) {
            tagged_data(data, fixture_modes[mode].data_count, tag);
            int count = mbe_quality_frame_from_data(fixture_modes[mode].name, data, fixture_modes[mode].data_count,
                                                    frame, sizeof(frame));
            assert(count == (int)fixture_modes[mode].frame_count);
            for (int i = 0; i < count; ++i) {
                assert(frame[i] == 0 || frame[i] == 1);
            }
            assert_decoded_fixture((enum fixture_mode)mode, frame, data, 0, 0, 0);
            for (int row = 0; row < fixture_modes[mode].ecc_rows; ++row) {
                int first = 0;
                int end = row < 4 ? 23 : 15;
                if (mode == IMBE7100 && row < 2) {
                    first = 1;
                    end = row == 0 ? 19 : 24;
                } else if (mode >= AMBE2450 && row == 0) {
                    end = 24; /* Include the AMBE C0 overall parity bit. */
                }
                for (int column = first; column < end; ++column) {
                    memcpy(damaged, frame, (size_t)count);
                    damaged[row * fixture_modes[mode].columns + column] ^= 1;
                    /* Golay reports corrected data bits, not parity-bit damage.
                     * AMBE C0 additionally reports its overall parity correction. */
                    int errors = row >= 4 || column - first >= 11;
                    if (mode >= AMBE2450 && row == 0) {
                        errors = 1;
                    }
                    assert_decoded_fixture((enum fixture_mode)mode, damaged, data, row == 0 ? errors : 0,
                                           row != 0 ? errors : 0, mode <= IMBE7100 && row == 4);
                }
            }
        }
    }
}

/*
 * Run each API over the same sequence, not just a cold first frame. The first
 * three frames establish synthesis history; later frames exercise corrected
 * C0/protected (and IMBE C4) context without comparing noisy modes to each other.
 */
#define DEFINE_FRAME_EQUIVALENCE(NAME, MODE, COLS, FRAME_API, DATA_API)                                                \
    static void NAME(void) {                                                                                           \
        char frames[8][184], parameters[8][88], decoded[88];                                                           \
        float expected_f[8][160], pcm[160];                                                                            \
        short expected_s[8][160], pcm_s[160];                                                                          \
        mbe_process_result expected_result[8], result;                                                                 \
        mbe_soft_bit soft[184];                                                                                        \
        mbe_parms cur, prev, enh;                                                                                      \
        for (int n = 0; n < 8; ++n) {                                                                                  \
            tagged_data(parameters[n], fixture_modes[MODE].data_count, (unsigned)n);                                   \
            memset(parameters[n], 0, 6); /* Keep the fundamental in the ordinary speech range. */                      \
            assert(mbe_quality_frame_from_data(fixture_modes[MODE].name, parameters[n],                                \
                                               fixture_modes[MODE].data_count, frames[n], sizeof(frames[n]))           \
                   == (int)fixture_modes[MODE].frame_count);                                                           \
            if (n == 4) {                                                                                              \
                frames[n][(MODE) == IMBE7100 ? 12 : 11] ^= 1;                                                          \
            } else if (n == 5) {                                                                                       \
                frames[n][(COLS) + ((MODE) == IMBE7100 ? 12 : 11)] ^= 1;                                               \
            } else if (n == 6 && (MODE) <= IMBE7100) {                                                                 \
                frames[n][(size_t)4 * (COLS)] ^= 1;                                                                    \
            }                                                                                                          \
        }                                                                                                              \
        for (int path = 0; path < 6; ++path) {                                                                         \
            init_params(&cur, &prev, &enh);                                                                            \
            for (int n = 0; n < 8; ++n) {                                                                              \
                int ret;                                                                                               \
                soft_from_bits(frames[n], soft, fixture_modes[MODE].frame_count);                                      \
                switch (path) {                                                                                        \
                    case 0:                                                                                            \
                        ret = FRAME_API##Framef(pcm, &result, (const char (*)[COLS])frames[n], decoded, &cur, &prev,   \
                                                &enh);                                                                 \
                        break;                                                                                         \
                    case 1:                                                                                            \
                        ret = FRAME_API##Frame(pcm_s, &result, (const char (*)[COLS])frames[n], decoded, &cur, &prev,  \
                                               &enh);                                                                  \
                        break;                                                                                         \
                    case 2:                                                                                            \
                        ret = FRAME_API##SoftFramef(pcm, &result, (const mbe_soft_bit(*)[COLS])soft, decoded, &cur,    \
                                                    &prev, &enh);                                                      \
                        break;                                                                                         \
                    case 3:                                                                                            \
                        ret = FRAME_API##SoftFrame(pcm_s, &result, (const mbe_soft_bit(*)[COLS])soft, decoded, &cur,   \
                                                   &prev, &enh);                                                       \
                        break;                                                                                         \
                    default:                                                                                           \
                        ret = decode_fixture(MODE, frames[n], decoded, &result, 0);                                    \
                        assert(ret >= 0);                                                                              \
                        ret = path == 4 ? DATA_API##Dataf(pcm, &result, decoded, &cur, &prev, &enh)                    \
                                        : DATA_API##Data(pcm_s, &result, decoded, &cur, &prev, &enh);                  \
                        break;                                                                                         \
                }                                                                                                      \
                assert_result_total(&result, ret);                                                                     \
                assert(memcmp(decoded, parameters[n], fixture_modes[MODE].data_count) == 0);                           \
                if ((path & 1) == 0) {                                                                                 \
                    assert_float_pcm_sane(pcm);                                                                        \
                    mbe_floattoshort(pcm, pcm_s);                                                                      \
                }                                                                                                      \
                if (path == 0) {                                                                                       \
                    memcpy(expected_f[n], pcm, sizeof(pcm));                                                           \
                    memcpy(expected_s[n], pcm_s, sizeof(pcm_s));                                                       \
                    expected_result[n] = result;                                                                       \
                } else {                                                                                               \
                    assert_same_result(&expected_result[n], &result, path == 2 || path == 3);                          \
                    assert(memcmp(expected_s[n], pcm_s, sizeof(pcm_s)) == 0);                                          \
                    if ((path & 1) == 0) {                                                                             \
                        assert(float_array_equal_exact(expected_f[n], pcm, 160));                                      \
                    }                                                                                                  \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
    }

DEFINE_FRAME_EQUIVALENCE(test_imbe7200_equivalence, IMBE7200, 23, mbe_processImbe7200x4400, mbe_processImbe4400)
DEFINE_FRAME_EQUIVALENCE(test_imbe7100_equivalence, IMBE7100, 24, mbe_processImbe7100x4400, mbe_processImbe4400)
DEFINE_FRAME_EQUIVALENCE(test_ambe2450_equivalence, AMBE2450, 24, mbe_processAmbe3600x2450, mbe_processAmbe2450)
DEFINE_FRAME_EQUIVALENCE(test_ambe2400_equivalence, AMBE2400, 24, mbe_processAmbe3600x2400, mbe_processAmbe2400)
#undef DEFINE_FRAME_EQUIVALENCE

static void
test_imbe_fundamental_roundtrip(void) {
    static const unsigned b0_indices[] = {0, 1, 2, 3, 4, 5, 85, 86};
    char data[88], frame7200[184], frame7100[184], decoded[88];
    float pcm7200[160], pcm7100[160], warm_pcm[160];
    mbe_parms warm_cur, warm_prev, warm_enh, cur, prev, enh;
    mbe_process_result result7200, result7100;
    char warm_data[88] = {0};

    init_params(&warm_cur, &warm_prev, &warm_enh);
    for (int n = 0; n < 3; ++n) {
        mbe_initProcessResult(&result7200);
        assert(mbe_processImbe4400Dataf(warm_pcm, &result7200, warm_data, &warm_cur, &warm_prev, &warm_enh) == 0);
    }
    /* Exhaustive b0 covers every K transition, including out-of-speech indices.
     * Seven bit planes uniquely tag all payload positions and vary the status
     * bit. Framef below verifies each round trip; soft ECC is covered separately. */
    for (unsigned b0 = 0; b0 < 256; ++b0) {
        for (unsigned tag = 0; tag < 7; ++tag) {
            for (unsigned i = 0; i < sizeof(data); ++i) {
                data[i] = (char)((i >> tag) & 1u);
            }
            for (unsigned i = 0; i < 8; ++i) {
                data[b0_indices[i]] = (char)((b0 >> (7u - i)) & 1u);
            }
            data[87] = (char)(tag & 1u);
            assert(mbe_quality_frame_from_data("imbe7200", data, sizeof(data), frame7200, sizeof(frame7200)) == 184);
            assert(mbe_quality_frame_from_data("imbe7100", data, sizeof(data), frame7100, sizeof(frame7100)) == 168);

            cur = warm_cur;
            prev = warm_prev;
            enh = warm_enh;
            mbe_setThreadRngSeed(0x12345678u);
            assert(mbe_processImbe7200x4400Framef(pcm7200, &result7200, (const char (*)[23])frame7200, decoded, &cur,
                                                  &prev, &enh)
                   == 0);
            assert(memcmp(decoded, data, sizeof(data)) == 0);
            assert(((result7200.flags & MBE_PROCESS_FLAG_REPEAT) != 0u) == (b0 > 207));
            cur = warm_cur;
            prev = warm_prev;
            enh = warm_enh;
            mbe_setThreadRngSeed(0x12345678u);
            assert(mbe_processImbe7100x4400Framef(pcm7100, &result7100, (const char (*)[24])frame7100, decoded, &cur,
                                                  &prev, &enh)
                   == 0);
            assert(memcmp(decoded, data, sizeof(data)) == 0);
            assert_same_result(&result7200, &result7100, 0);
            assert_float_pcm_sane(pcm7200);
            assert_float_pcm_sane(pcm7100);
            assert(float_array_equal_exact(pcm7200, pcm7100, 160));
        }
    }
}

static void
assert_ambe2400_tone_spectrum(const float pcm[20][160], double frequency, int tones_enabled) {
    double energy = 0.0;
    double real = 0.0;
    double imag = 0.0;
    /* The last 16 frames contain an integer number of cycles for IDs 5/6/7. */
    const int count = 16 * 160;
    for (int n = 0; n < count; ++n) {
        double sample = pcm[4 + n / 160][n % 160];
        double phase = 2.0 * 3.14159265358979323846 * frequency * n / 8000.0;
        if (!tones_enabled) {
            assert(sample == 0.0);
        }
        energy += sample * sample;
        real += sample * cos(phase);
        imag += sample * sin(phase);
    }
    if (tones_enabled) {
        assert(energy > 0.0);
        assert(2.0 * (real * real + imag * imag) / (count * energy) > 0.9999);
    }
}

static void
test_ambe2400_accepted_tones(void) {
    static const int low_indices[] = {9, 42, 43, 10, 11};
    static const double frequencies[] = {156.25, 187.5, 218.75};

    static const struct {
        int c0;
        int protected;
        int accepted;
    } cases[] = {
        {0, 0, 1},
        {1, 1, 1}, /* Both error counts immediately below their rejection gates. */
        {2, 0, 0},
        {1, 2, 0},
    };

    char data[49], frame[96], damaged[96], decoded[49];
    float tone[20][160], noise[20][160], pcm[20][160];
    mbe_process_result context, result, frame_result[20];
    mbe_parms cur, prev, enh;

    for (int id = 5; id <= 7; ++id) {
        memset(data, 0, sizeof(data));
        memset(data, 1, 6); /* b0 = 126: bits 0..5 = 1, bit 48 = 0. */
        data[8] = 1;        /* Selector bits 6..8 = 001. */
        for (int bit = 0; bit < 5; ++bit) {
            data[low_indices[bit]] = (char)((id >> (4 - bit)) & 1);
        }
        assert(mbe_quality_frame_from_data("ambe2400", data, sizeof(data), frame, sizeof(frame)) == 96);

        init_params(&cur, &prev, &enh);
        double tone_energy = 0.0;
        for (int n = 0; n < 20; ++n) {
            mbe_synthesizeTonefdstar(tone[n], data, &cur, id);
            mbe_synthesizeComfortNoisef(noise[n]);
            for (int i = 0; i < 160; ++i) {
                tone_energy += (double)tone[n][i] * tone[n][i];
            }
        }
        /* The linked renderer is the oracle in both ordinary and NOTONES builds;
         * NOTONES is a private library definition, not a test compile definition. */
        int tones_enabled = tone_energy > 0.0;

        for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
            memcpy(damaged, frame, sizeof(damaged));
            for (int bit = 0; bit < cases[c].c0; ++bit) {
                damaged[12 + bit] ^= 1; /* Golay counts corrected data bits, not pairs of parity errors. */
            }
            for (int bit = 0; bit < cases[c].protected; ++bit) {
                damaged[24 + 11 + bit] ^= 1;
            }
            int total = cases[c].c0 + cases[c].protected;
            assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])damaged, decoded, &context) == total);
            assert(memcmp(decoded, data, sizeof(data)) == 0);
            assert(context.c0_errors == cases[c].c0);
            assert(context.protected_errors == cases[c].protected);

            /* Framef supplies real ECC context; Dataf must honor that same context. */
            for (int path = 0; path < 2; ++path) {
                init_params(&cur, &prev, &enh);
                for (int n = 0; n < 20; ++n) {
                    int ret;
                    if (path == 0) {
                        ret = mbe_processAmbe3600x2400Framef(pcm[n], &result, (const char (*)[24])damaged, decoded,
                                                             &cur, &prev, &enh);
                    } else {
                        result = context;
                        ret = mbe_processAmbe2400Dataf(pcm[n], &result, decoded, &cur, &prev, &enh);
                    }
                    assert_result_total(&result, ret);
                    assert(ret == total);
                    assert(memcmp(decoded, data, sizeof(data)) == 0);
                    assert(result.flags
                           == (MBE_PROCESS_FLAG_C0_VALID | (cases[c].accepted ? MBE_PROCESS_FLAG_TONE : 0u)));
                    assert_float_pcm_sane(pcm[n]);
                    assert(float_array_equal_exact(pcm[n], cases[c].accepted ? tone[n] : noise[n], 160));
                    if (path == 0) {
                        frame_result[n] = result;
                    } else {
                        assert_same_result(&frame_result[n], &result, 0);
                    }
                }
                if (cases[c].accepted) {
                    assert_ambe2400_tone_spectrum((const float (*)[160])pcm, frequencies[id - 5], tones_enabled);
                }
            }
        }
    }
}

int
main(void) {
    test_ambe2400_frame_paths();
    test_imbe7200_frame_paths();
    test_imbe7100_frame_paths();
    test_fixture_validation();
    test_fixture_ecc();
    test_imbe7200_equivalence();
    test_imbe7100_equivalence();
    test_ambe2450_equivalence();
    test_ambe2400_equivalence();
    test_imbe_fundamental_roundtrip();
    test_ambe2400_accepted_tones();
    return 0;
}
