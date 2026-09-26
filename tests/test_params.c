// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Targeted parameter-derivation tests for IMBE/AMBE paths.
 */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <string.h>

#include "mbe_tone.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

/**
 * @brief Zero-initialize a bit array of length n.
 * @param bits Pointer to bit array (0/1 chars) to clear.
 * @param n    Number of entries to set to zero.
 */
static void
set_bits_zero(char* bits, int n) {
    memset(bits, 0, (size_t)n);
}

static void
init_result_total(mbe_process_result* result, int total_errors) {
    mbe_initProcessResult(result);
    result->total_errors = total_errors;
}

static int
result_has_marker(const mbe_process_result* result, char marker) {
    char status[96];
    mbe_formatProcessResult(status, sizeof(status), result);
    return strchr(status, marker) != NULL;
}

/**
 * @brief Compose IMBE 7200x4400 b0 into the imbe_d vector.
 *
 * Matches the extraction in `mbe_decodeImbe4400Parms`, where b0 is formed
 * from `imbe_d[0..5], imbe_d[85], imbe_d[86]` (MSB-first).
 *
 * @param imbe_d IMBE parameter bit vector (88 entries), modified in-place.
 * @param b0     8-bit value to encode into the proper positions.
 */
static void
set_imbe7200_b0(char imbe_d[88], int b0) {
    for (int i = 0; i < 8; ++i) {
        int bit = (b0 >> (7 - i)) & 1; // MSB-first
        if (i < 6) {
            imbe_d[i] = (char)bit;
        } else if (i == 6) {
            imbe_d[85] = (char)bit;
        } else { // i == 7
            imbe_d[86] = (char)bit;
        }
    }
}

/**
 * @brief Compose AMBE 3600x2450 b0 into the ambe_d vector.
 *
 * Matches the extraction in `mbe_decodeAmbe2450Parms`, mapping
 * `b0 = d0<<6 | d1<<5 | d2<<4 | d3<<3 | d37<<2 | d38<<1 | d39`.
 *
 * @param ambe_d AMBE parameter bit vector (49 entries), modified in-place.
 * @param b0     7-bit value to encode into the proper positions.
 */
static void
set_ambe2450_b0(char ambe_d[49], int b0) {
    ambe_d[0] = (char)((b0 >> 6) & 1);
    ambe_d[1] = (char)((b0 >> 5) & 1);
    ambe_d[2] = (char)((b0 >> 4) & 1);
    ambe_d[3] = (char)((b0 >> 3) & 1);
    ambe_d[37] = (char)((b0 >> 2) & 1);
    ambe_d[38] = (char)((b0 >> 1) & 1);
    ambe_d[39] = (char)((b0 >> 0) & 1);
}

/**
 * @brief Set AMBE 2450 tone signature bits used by tone classification.
 *
 * This sets U0 tone-check bits to 63 and keeps U3 tone-check nibble at 0.
 *
 * @param ambe_d AMBE parameter bit vector (49 entries), modified in-place.
 */
static void
set_ambe2450_tone_signature(char ambe_d[49]) {
    ambe_d[0] = 1;
    ambe_d[1] = 1;
    ambe_d[2] = 1;
    ambe_d[3] = 1;
    ambe_d[4] = 1;
    ambe_d[5] = 1;
    ambe_d[45] = 0;
    ambe_d[46] = 0;
    ambe_d[47] = 0;
    ambe_d[48] = 0;
}

/**
 * @brief Compose AMBE 2450 tone ID1 and the U1 low nibble.
 *
 * Tone ID1 is U1[0..7] (ambe_d[12..19]). U1[8..11] (ambe_d[20..23])
 * belongs to a different field and should not affect tone-ID validity.
 *
 * @param ambe_d      AMBE parameter bit vector (49 entries), modified in-place.
 * @param id1         8-bit tone ID value.
 * @param low_nibble  4-bit U1 low nibble value.
 */
static void
set_ambe2450_tone_id1_and_u1_low_nibble(char ambe_d[49], int id1, int low_nibble) {
    for (int i = 0; i < 8; ++i) {
        ambe_d[12 + i] = (char)((id1 >> (7 - i)) & 1);
    }
    for (int i = 0; i < 4; ++i) {
        ambe_d[20 + i] = (char)((low_nibble >> (3 - i)) & 1);
    }
}

/** Every harmonic 1..L of the model lies below Nyquist (L * w0 < pi). */
static int
harmonics_below_nyquist(const mbe_parms* mp) {
    return (mp->w0 > 0.0f) && (mp->L >= 1) && ((float)mp->L * mp->w0 < (float)M_PI);
}

/**
 * @brief Seed synthetic voiced/unvoiced parameters for synthesis behavior tests.
 *
 * @param cur  Output current parameter set.
 * @param prev Output previous parameter set (copy of cur).
 */
static void
seed_speech_params(mbe_parms* cur, mbe_parms* prev) {
    mbe_parms dummy;
    mbe_initMbeParms(cur, prev, &dummy);
    cur->w0 = 0.10f;
    cur->L = 12;
    for (int l = 1; l <= cur->L; ++l) {
        cur->Vl[l] = (l % 3) ? 1 : 0;
        cur->Ml[l] = 0.03f + (0.001f * (float)l);
        cur->PHIl[l] = 0.0f;
        cur->PSIl[l] = 0.0f;
    }
    *prev = *cur;
}

/**
 * @brief Compare two floats for approximate equality.
 * @param a   First value.
 * @param b   Second value.
 * @param eps Absolute tolerance.
 * @return 1 if |a - b| <= eps, else 0.
 */
static int
approx_equal(float a, float b, float eps) {
    float diff = a - b;
    if (diff < 0) {
        diff = -diff;
    }
    return diff <= eps;
}

/**
 * @brief Check whether two float arrays are exactly element-wise equal.
 * @param a First array.
 * @param b Second array.
 * @param n Number of elements.
 * @return 1 if all elements match exactly, else 0.
 */
static int
float_bits_equal(float a, float b) {
    uint32_t a_bits;
    uint32_t b_bits;
    memcpy(&a_bits, &a, sizeof(a_bits));
    memcpy(&b_bits, &b, sizeof(b_bits));
    return a_bits == b_bits;
}

/**
 * @brief Check whether two float arrays are exactly element-wise equal.
 * @param a First array.
 * @param b Second array.
 * @param n Number of elements.
 * @return 1 if all elements match exactly, else 0.
 */
static int
float_array_equal_exact(const float* a, const float* b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (!float_bits_equal(a[i], b[i])) {
            return 0;
        }
    }
    return 1;
}

static int
float_bits_differ(float a, float b) {
    return !float_bits_equal(a, b);
}

/**
 * @brief Check whether two float arrays differ in at least one element.
 * @param a First array.
 * @param b Second array.
 * @param n Number of elements.
 * @return 1 if any element differs, else 0.
 */
static int
float_array_differs(const float* a, const float* b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (float_bits_differ(a[i], b[i])) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Test entry: IMBE and AMBE parameter derivation spot checks.
 */
int
main(void) {
    // IMBE 7200x4400: verify w0/L/K from b0
    {
        char imbe_d[88];
        mbe_parms cur = {0}, prev = {0};
        mbe_parms dummy;
        mbe_initMbeParms(&cur, &prev, &dummy);

        set_bits_zero(imbe_d, 88);
        // Pick a few valid b0 values well inside range to avoid edge cases
        int b0_values[] = {0, 100, 206};
        for (unsigned i = 0; i < sizeof(b0_values) / sizeof(b0_values[0]); ++i) {
            set_bits_zero(imbe_d, 88);
            set_imbe7200_b0(imbe_d, b0_values[i]);
            int rc = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
            // valid voice frame expected
            assert(rc == 0);
            (void)rc; /* Suppress unused-but-set warning when asserts disabled */
            float w0_expected = (float)((4.0 * M_PI) / ((double)b0_values[i] + 39.5));
            int L_expected = (int)(0.9254 * (int)((M_PI / w0_expected) + 0.25));
            int K_expected = (L_expected < 37) ? ((L_expected + 2) / 3) : 12;
            assert(approx_equal(cur.w0, w0_expected, 1e-6f));
            assert(cur.L == L_expected);
            assert(cur.K == K_expected);
            (void)w0_expected; /* Suppress unused-but-set warning when asserts disabled */
            (void)L_expected;
            (void)K_expected;
        }
    }

    // JMBE parity: phase arrays initialize to zero
    {
        mbe_parms cur = {0}, prev = {0}, prev_enh = {0};
        mbe_initMbeParms(&cur, &prev, &prev_enh);
        for (int l = 0; l <= 56; ++l) {
            assert(approx_equal(prev.PSIl[l], 0.0f, 1e-7f));
            assert(approx_equal(cur.PSIl[l], 0.0f, 1e-7f));
            assert(approx_equal(prev_enh.PSIl[l], 0.0f, 1e-7f));
        }
    }

    // AMBE tone lookup: shared table/formula behavior for synthesis and tone validation
    {
        struct {
            int id;
            int valid;
            float freq1;
            float freq2;
        } cases[] = {
            {4, 0, 0.0f, 0.0f},         {5, 1, 156.25f, 156.25f}, {6, 1, 187.5f, 187.5f},    {7, 1, 218.75f, 218.75f},
            {122, 1, 3812.5f, 3812.5f}, {123, 0, 0.0f, 0.0f},     {128, 1, 1336.0f, 941.0f}, {141, 1, 1633.0f, 941.0f},
            {163, 1, 490.0f, 350.0f},   {255, 0, 0.0f, 0.0f},
        };

        for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
            float freq1 = -1.0f;
            float freq2 = -1.0f;
            int valid = mbe_tone_lookup_freqs(cases[i].id, &freq1, &freq2);
            assert(valid == cases[i].valid);
            assert(mbe_tone_id_is_valid(cases[i].id) == cases[i].valid);
            assert(approx_equal(freq1, cases[i].freq1, 1e-6f));
            assert(approx_equal(freq2, cases[i].freq2, 1e-6f));
        }
    }

    // AMBE+2 3600x2450: verify w0/L mapping from table for a couple of b0 values
    // We avoid including internal headers to prevent duplicate symbol definitions; instead we
    // embed expected reference values for select indices taken from the public table in source.
    {
        char ambe_d[49];
        mbe_parms cur = {0}, prev = {0};
        mbe_parms dummy;
        mbe_initMbeParms(&cur, &prev, &dummy);

        struct {
            int b0;
            float f0; // table value AmbeW0table[b0]
            int L;
        } cases[] = {
            {0, 0.049971f, 9},
            {119, 0.008125f, 56},
        };

        for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
            set_bits_zero(ambe_d, 49);
            set_ambe2450_b0(ambe_d, cases[i].b0);
            int rc = mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev);
            // For these b0 values we expect normal voice decode (rc == 0)
            assert(rc == 0);
            (void)rc; /* Suppress unused-but-set warning when asserts disabled */
            float w0_expected = cases[i].f0 * (float)(2.0 * M_PI);
            assert(approx_equal(cur.w0, w0_expected, 1e-6f));
            assert(cur.L == cases[i].L);
            (void)w0_expected; /* Suppress unused-but-set warning when asserts disabled */
        }
    }

    // AMBE 2450 silence frames (TIA-102.BABA-1 4.1 eqs 1-3): both b0 124 and
    // 125 use w0 = 2*pi/32, L = 14, all bands unvoiced.
    {
        char ambe_d[49];
        mbe_parms cur = {0}, prev = {0};
        mbe_parms dummy;
        const float w0_silence = (float)(2.0 * M_PI / 32.0);

        for (int b0 = 124; b0 <= 125; ++b0) {
            mbe_initMbeParms(&cur, &prev, &dummy);
            set_bits_zero(ambe_d, 49);
            set_ambe2450_b0(ambe_d, b0);
            assert(mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev) == MBE_AMBE2450_FRAME_SILENCE);
            assert(cur.L == 14);
            assert(approx_equal(cur.w0, w0_silence, 1e-6f));
            assert((float)cur.L * cur.w0 < (float)M_PI);
            for (int l = 1; l <= cur.L; ++l) {
                assert(cur.Vl[l] == 0);
            }
        }
        (void)w0_silence;
    }

    // Nyquist: every model the decoders can emit keeps all L harmonics below
    // pi. JMBE's AMBE silence/initial model (w0 = pi^2/16, L = 15) did not.
    {
        char ambe_d[49];
        char imbe_d[88];
        float out[160];
        mbe_parms cur, prev, enh;

        for (int b0 = 0; b0 < 128; ++b0) {
            mbe_initMbeParms(&cur, &prev, &enh);
            set_bits_zero(ambe_d, 49);
            set_ambe2450_b0(ambe_d, b0);
            int rc = mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev);
            if (rc == MBE_AMBE2450_FRAME_VOICE || rc == MBE_AMBE2450_FRAME_SILENCE) {
                assert(harmonics_below_nyquist(&cur));
            }
        }
        /* Initial AMBE model, reached as the repeat model of an erasure at stream start. */
        mbe_initMbeParms(&cur, &prev, &enh);
        set_bits_zero(ambe_d, 49);
        set_ambe2450_b0(ambe_d, 120);
        assert(mbe_processAmbe2450Dataf(out, NULL, ambe_d, &cur, &prev, &enh) >= 0);
        assert(harmonics_below_nyquist(&cur));
        assert(harmonics_below_nyquist(&prev));

        /* D-STAR: every b0 bit pattern (bits 0..5 and 48) that decodes as voice. */
        for (int pattern = 0; pattern < 128; ++pattern) {
            mbe_initMbeParms(&cur, &prev, &enh);
            set_bits_zero(ambe_d, 49);
            for (int i = 0; i < 6; ++i) {
                ambe_d[i] = (char)((pattern >> (6 - i)) & 1);
            }
            ambe_d[48] = (char)(pattern & 1);
            if (mbe_decodeAmbe2400Parms(ambe_d, &cur, &prev) == 0) {
                assert(harmonics_below_nyquist(&cur));
            }
        }
        /* D-STAR state after the standard silence frame (b0 127, tone index 128). */
        mbe_initMbeParms(&cur, &prev, &enh);
        set_bits_zero(ambe_d, 49);
        for (int i = 0; i < 6; ++i) {
            ambe_d[i] = 1;
        }
        ambe_d[48] = 1;
        assert(mbe_processAmbe2400Dataf(out, NULL, ambe_d, &cur, &prev, &enh) >= 0);
        assert(harmonics_below_nyquist(&prev));
        assert(harmonics_below_nyquist(&enh));

        /* IMBE: every voice b0 and the generic initial state. */
        for (int b0 = 0; b0 <= 207; ++b0) {
            mbe_initMbeParms(&cur, &prev, &enh);
            set_bits_zero(imbe_d, 88);
            set_imbe7200_b0(imbe_d, b0);
            if (mbe_decodeImbe4400Parms(imbe_d, &cur, &prev) == 0) {
                assert(harmonics_below_nyquist(&cur));
            }
        }
        mbe_initMbeParms(&cur, &prev, &enh);
        assert(harmonics_below_nyquist(&prev));
    }

    // Shared synthesizer Nyquist guard: voiced harmonics at or above pi cannot
    // be represented at 8 kHz and are skipped instead of aliasing. Only the
    // above-Nyquist harmonics carry amplitude here, so the output is silent.
    {
        float out[160];
        mbe_parms cur, prev;
        seed_speech_params(&cur, &prev);
        cur.w0 = 0.105f;
        cur.L = 36;
        for (int l = 1; l <= 56; ++l) {
            cur.Vl[l] = 1;
            cur.Ml[l] = ((float)l * cur.w0 >= (float)M_PI) ? 1.0f : 0.0f;
        }
        prev = cur;
        mbe_synthesizeSpeechf(out, &cur, &prev);
        for (int i = 0; i < 160; ++i) {
            assert(fabsf(out[i]) < 1e-6f);
        }
    }

    // AMBE 2450 Dataf: absent C0-valid context, repeat decision must depend only on total errors
    {
        char ambe_d[49];
        float out[160];
        mbe_process_result result_a;
        mbe_process_result result_b;
        mbe_parms cur_a = {0}, prev_a = {0}, prev_enh_a = {0};
        mbe_parms cur_b = {0}, prev_b = {0}, prev_enh_b = {0};

        set_bits_zero(ambe_d, 49);
        set_ambe2450_b0(ambe_d, 0);

        mbe_initMbeParms(&cur_a, &prev_a, &prev_enh_a);
        init_result_total(&result_a, 5);
        result_a.c0_errors = 0;
        assert(mbe_processAmbe2450Dataf(out, &result_a, ambe_d, &cur_a, &prev_a, &prev_enh_a) >= 0);

        mbe_initMbeParms(&cur_b, &prev_b, &prev_enh_b);
        init_result_total(&result_b, 5);
        result_b.c0_errors = 4;
        result_b.protected_errors = 1;
        assert(mbe_processAmbe2450Dataf(out, &result_b, ambe_d, &cur_b, &prev_b, &prev_enh_b) >= 0);

        assert(cur_a.repeatCount == cur_b.repeatCount);
        assert(result_has_marker(&result_a, 'R') == result_has_marker(&result_b, 'R'));
    }

    // IMBE 4400 Dataf: absent C0-valid context, repeat decision must depend only on total errors
    {
        char imbe_d[88];
        float out[160];
        mbe_process_result result_a;
        mbe_process_result result_b;
        mbe_parms cur_a = {0}, prev_a = {0}, prev_enh_a = {0};
        mbe_parms cur_b = {0}, prev_b = {0}, prev_enh_b = {0};

        set_bits_zero(imbe_d, 88);
        set_imbe7200_b0(imbe_d, 0);

        mbe_initMbeParms(&cur_a, &prev_a, &prev_enh_a);
        init_result_total(&result_a, 11);
        result_a.c0_errors = 0;
        assert(mbe_processImbe4400Dataf(out, &result_a, imbe_d, &cur_a, &prev_a, &prev_enh_a) >= 0);

        mbe_initMbeParms(&cur_b, &prev_b, &prev_enh_b);
        init_result_total(&result_b, 11);
        result_b.c0_errors = 2;
        result_b.protected_errors = 9;
        assert(mbe_processImbe4400Dataf(out, &result_b, imbe_d, &cur_b, &prev_b, &prev_enh_b) >= 0);

        assert(cur_a.repeatCount == cur_b.repeatCount);
        assert(result_has_marker(&result_a, 'R') == result_has_marker(&result_b, 'R'));
    }

    // Frame repeats keep this frame's error accounting and continue the noise
    // sequence (IMBE and D-STAR). Copying the stale previous parameter set
    // froze the error-rate recursion and replayed the previous frame's noise.
    {
        char imbe_d[88];
        char ambe_d[49];
        float out[160];
        mbe_process_result result;
        mbe_parms cur, prev, enh;

        set_bits_zero(imbe_d, 88);
        mbe_initMbeParms(&cur, &prev, &enh);
        for (int n = 0; n < 2; ++n) {
            init_result_total(&result, 1);
            assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &enh) >= 0);
        }
        float er_before = prev.errorRate;
        float seed_before = enh.noiseSeed;
        init_result_total(&result, 7); /* Dataf fallback: total > 5 repeats */
        assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &enh) >= 0);
        assert(result_has_marker(&result, 'R'));
        assert(fabsf(prev.errorRate - ((0.95f * er_before) + (0.000365f * 7.0f))) < 1e-7f);
        assert(cur.errorCountTotal == 7);
        assert(float_bits_differ(enh.noiseSeed, seed_before));

        set_bits_zero(ambe_d, 49);
        mbe_initMbeParms(&cur, &prev, &enh);
        for (int n = 0; n < 2; ++n) {
            init_result_total(&result, 1);
            assert(mbe_processAmbe2400Dataf(out, &result, ambe_d, &cur, &prev, &enh) >= 0);
        }
        er_before = prev.errorRate;
        seed_before = enh.noiseSeed;
        init_result_total(&result, 4); /* D-STAR: total > 3 repeats */
        assert(mbe_processAmbe2400Dataf(out, &result, ambe_d, &cur, &prev, &enh) >= 0);
        assert(result_has_marker(&result, 'R'));
        assert(fabsf(prev.errorRate - ((0.95f * er_before) + (0.001064f * 4.0f))) < 1e-7f);
        assert(cur.errorCountTotal == 4);
        assert(float_bits_differ(enh.noiseSeed, seed_before));
        (void)er_before;
    }

    // AMBE C0 Golay24 parity behavior: isolated parity-bit error is corrected
    {
        char ambe_fr[4][24] = {{0}};
        int errs = mbe_eccAmbe3600x2450C0(ambe_fr);
        assert(errs == 0);

        ambe_fr[0][0] = 1; /* parity bit only */
        errs = mbe_eccAmbe3600x2450C0(ambe_fr);
        assert(errs == 1);
        assert(ambe_fr[0][0] == 0);
    }

    // Soft AMBE C0 must preserve Golay parity-bit errors for hard-equivalent inputs
    {
        char ambe_fr[4][24] = {{0}};
        char hard_d[49];
        char soft_d[49];
        mbe_soft_bit soft_fr[4][24];
        mbe_process_result hard_result;
        mbe_process_result soft_result;

        ambe_fr[0][1] = 1; /* protected Golay parity bit */
        typedef const char hard_ambe_frame[24];
        hard_ambe_frame* hard_fr = (hard_ambe_frame*)ambe_fr;
        int hard_ret = mbe_decodeAmbe3600x2450Frame(hard_fr, hard_d, &hard_result);
        assert(mbe_softBitsFromHard(&ambe_fr[0][0], &soft_fr[0][0], sizeof(soft_fr) / sizeof(soft_fr[0][0]), 255u)
               == 0);
        typedef const mbe_soft_bit soft_ambe_frame[24];
        soft_ambe_frame* soft_fr_const = (soft_ambe_frame*)soft_fr;
        int soft_ret = mbe_decodeAmbe3600x2450SoftFrame(soft_fr_const, soft_d, &soft_result);

        assert(hard_ret == soft_ret);
        assert(hard_result.c0_errors == 1);
        assert(soft_result.c0_errors == hard_result.c0_errors);
        assert(soft_result.total_errors == hard_result.total_errors);
        assert(memcmp(hard_d, soft_d, sizeof(hard_d)) == 0);
    }

    // AMBE 2450 tone and erasure frames (TIA-102.BABA-1 5.6, 7.3): the repeat
    // criteria apply before tone classification, and an erasure is a frame
    // repeat of the last synthesized frame (not JMBE's W120 model).
    {
        char ambe_d[49];
        float out[160];
        mbe_process_result result;
        mbe_parms cur = {0}, prev = {0}, prev_enh = {0};

        set_bits_zero(ambe_d, 49);
        set_ambe2450_tone_signature(ambe_d);
        set_ambe2450_tone_id1_and_u1_low_nibble(ambe_d, 7, 0x0);

        mbe_initMbeParms(&cur, &prev, &prev_enh);
        init_result_total(&result, 3);
        assert(mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh) >= 0);
        assert(result_has_marker(&result, 'T'));
        assert(!result_has_marker(&result, 'R'));

        /* Without C0 context the Dataf fallback treats total > 3 as a repeat. */
        mbe_initMbeParms(&cur, &prev, &prev_enh);
        init_result_total(&result, 4);
        assert(mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh) >= 0);
        assert(!result_has_marker(&result, 'T'));
        assert(result_has_marker(&result, 'R'));

        set_bits_zero(ambe_d, 49);
        set_ambe2450_b0(ambe_d, 120);
        mbe_initMbeParms(&cur, &prev, &prev_enh);
        init_result_total(&result, 0);
        assert(mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh) >= 0);
        assert(result_has_marker(&result, 'E'));
        assert(result_has_marker(&result, 'R'));
        assert(!result_has_marker(&result, 'M'));
        assert(cur.repeatCount == 1);
        assert(prev.repeatCount == 1);
        /* History stays at the spec initial state (L = 15, gamma = 0). */
        assert(prev.L == 15);
        assert(float_bits_equal(prev.gamma, 0.0f));
        assert(cur.L == 15);
    }

    // AMBE tone ID validity must depend on ID1 only, not U1 low nibble bits
    {
        char ambe_d_a[49];
        char ambe_d_b[49];
        float out[160];
        mbe_process_result result;
        const float custom_w0 = 0.20f;
        const int custom_L = 22;

        set_bits_zero(ambe_d_a, 49);
        set_ambe2450_tone_signature(ambe_d_a);
        set_ambe2450_tone_id1_and_u1_low_nibble(ambe_d_a, 7, 0x0);

        set_bits_zero(ambe_d_b, 49);
        set_ambe2450_tone_signature(ambe_d_b);
        set_ambe2450_tone_id1_and_u1_low_nibble(ambe_d_b, 7, 0xF);

        mbe_parms cur_a = {0}, prev_a = {0}, prev_enh_a = {0};
        mbe_initMbeParms(&cur_a, &prev_a, &prev_enh_a);
        cur_a.w0 = custom_w0;
        cur_a.L = custom_L;
        prev_a = cur_a;
        prev_enh_a = cur_a;
        cur_a.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev_a.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev_enh_a.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev_a.repeatCount = MBE_MAX_FRAME_REPEATS;

        init_result_total(&result, 0);
        assert(mbe_processAmbe2450Dataf(out, &result, ambe_d_a, &cur_a, &prev_a, &prev_enh_a) >= 0);
        assert(result_has_marker(&result, 'T'));
        assert(approx_equal(cur_a.w0, custom_w0, 1e-6f));
        assert(cur_a.L == custom_L);

        mbe_parms cur_b = {0}, prev_b = {0}, prev_enh_b = {0};
        mbe_initMbeParms(&cur_b, &prev_b, &prev_enh_b);
        cur_b.w0 = custom_w0;
        cur_b.L = custom_L;
        prev_b = cur_b;
        prev_enh_b = cur_b;
        cur_b.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev_b.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev_enh_b.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev_b.repeatCount = MBE_MAX_FRAME_REPEATS;

        init_result_total(&result, 0);
        assert(mbe_processAmbe2450Dataf(out, &result, ambe_d_b, &cur_b, &prev_b, &prev_enh_b) >= 0);
        assert(result_has_marker(&result, 'T'));
        assert(approx_equal(cur_b.w0, custom_w0, 1e-6f));
        assert(cur_b.L == custom_L);
    }

    // Shared synthesizer: the AMBE threshold does not error-rate mute (D-STAR
    // relies on this; the AMBE 2450 process path mutes before synthesis per
    // TIA-102.BABA-1 5.7), while the IMBE threshold does.
    {
        float out[160];
        mbe_parms cur = {0}, prev = {0};

        seed_speech_params(&cur, &prev);
        cur.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        cur.errorRate = 1.0f;
        cur.repeatCount = 0;
        float ambe_seed_before = cur.noiseSeed;
        mbe_synthesizeSpeechf(out, &cur, &prev);
        assert(float_bits_differ(cur.noiseSeed, ambe_seed_before));

        seed_speech_params(&cur, &prev);
        cur.mutingThreshold = MBE_MUTING_THRESHOLD_IMBE;
        cur.errorRate = 1.0f;
        cur.repeatCount = 0;
        float imbe_seed_before = cur.noiseSeed;
        float jmbe_noise[160];
        mbe_setThreadRngSeed(0x1234u);
        mbe_synthesizeSpeechf(out, &cur, &prev);
        assert(float_bits_equal(cur.noiseSeed, imbe_seed_before));

        /* The public synthesizer cannot tell the codec, so it keeps JMBE's
         * comfort noise at either threshold (D-STAR: a max-repeat mute). */
        mbe_setThreadRngSeed(0x1234u);
        mbe_synthesizeComfortNoisef(jmbe_noise);
        for (int i = 0; i < 160; ++i) {
            assert(float_bits_equal(out[i], jmbe_noise[i]));
        }
        seed_speech_params(&cur, &prev);
        cur.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        cur.repeatCount = MBE_MAX_FRAME_REPEATS;
        mbe_setThreadRngSeed(0x1234u);
        mbe_synthesizeSpeechf(out, &cur, &prev);
        for (int i = 0; i < 160; ++i) {
            assert(float_bits_equal(out[i], jmbe_noise[i]));
        }
    }

    // Error-rate mutes on the IMBE process paths: P25 (7200x4400 data path)
    // uses the TIA-102.BABA 7.8 noise, uniform in [-5, 5] on s(n); ProVoice
    // (7100x4400, not a TIA-102 codec) keeps JMBE's comfort noise.
    {
        float out[160];
        char imbe_d[88];
        static const char provoice_fr[7][24] = {{0}};
        mbe_parms cur, prev, enh;
        mbe_process_result result;

        set_bits_zero(imbe_d, 88);
        mbe_initMbeParms(&cur, &prev, &enh);
        prev.errorRate = 1.0f;
        mbe_initProcessResult(&result);
        assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &enh) >= 0);
        assert((result.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
        double sumsq = 0.0;
        for (int i = 0; i < 160; ++i) {
            assert(fabsf(out[i]) <= 5.0f);
            sumsq += (double)out[i] * (double)out[i];
        }
        double rms = sqrt(sumsq / 160.0);
        assert(rms > 1.5 && rms < 4.0);
        (void)rms;

        mbe_initMbeParms(&cur, &prev, &enh);
        prev.errorRate = 1.0f;
        mbe_initProcessResult(&result);
        assert(mbe_processImbe7100x4400Framef(out, &result, provoice_fr, imbe_d, &cur, &prev, &enh) >= 0);
        assert((result.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
        float peak = 0.0f;
        for (int i = 0; i < 160; ++i) {
            peak = fmaxf(peak, fabsf(out[i]));
        }
        /* Comfort noise is uniform in about [-14, 14]. */
        assert(peak > 5.0f && peak < 15.0f);
        (void)peak;

        /* Staged ProVoice decoding (hard and soft) carries the ProVoice
         * context into the IMBE 4400 data API and mutes bit-identically. */
        mbe_soft_bit provoice_soft[7][24];
        assert(mbe_softBitsFromHard(&provoice_fr[0][0], &provoice_soft[0][0], sizeof(provoice_fr), 255u) == 0);
        for (int soft = 0; soft < 2; ++soft) {
            float direct[160];
            mbe_process_result direct_result;
            mbe_initMbeParms(&cur, &prev, &enh);
            prev.errorRate = 1.0f;
            mbe_setThreadRngSeed(0x1234u);
            assert(
                (soft ? mbe_processImbe7100x4400SoftFramef(
                            direct, &direct_result, (const mbe_soft_bit(*)[24])provoice_soft, imbe_d, &cur, &prev, &enh)
                      : mbe_processImbe7100x4400Framef(direct, &direct_result, provoice_fr, imbe_d, &cur, &prev, &enh))
                >= 0);

            mbe_initMbeParms(&cur, &prev, &enh);
            prev.errorRate = 1.0f;
            mbe_setThreadRngSeed(0x1234u);
            assert((soft ? mbe_decodeImbe7100x4400SoftFrame((const mbe_soft_bit(*)[24])provoice_soft, imbe_d, &result)
                         : mbe_decodeImbe7100x4400Frame(provoice_fr, imbe_d, &result))
                   >= 0);
            assert((result.flags & MBE_PROCESS_FLAG_PROVOICE) != 0u);
            assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &enh) >= 0);
            assert(result.flags == direct_result.flags);
            assert((result.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
            for (int i = 0; i < 160; ++i) {
                assert(float_bits_equal(out[i], direct[i]));
            }
        }
    }

    // Muted IMBE frames should still advance adaptive smoothing state (JMBE parity)
    {
        float out[160];
        mbe_parms cur = {0}, prev = {0};

        seed_speech_params(&cur, &prev);
        cur.mutingThreshold = MBE_MUTING_THRESHOLD_IMBE;
        cur.errorRate = 1.0f;
        cur.repeatCount = 0;

        float local_before = cur.localEnergy;
        mbe_synthesizeSpeechf(out, &cur, &prev);
        assert(float_bits_differ(cur.localEnergy, local_before));
    }

    // JMBE parity: previous voiced phase is wrapped to [0, 2PI) before advancement
    {
        float out[160];
        mbe_parms cur = {0}, prev = {0};
        seed_speech_params(&cur, &prev);

        const int l_test = 5;
        const float raw_prev_phase = (20.0f * (float)M_PI) + 0.321f;
        prev.PSIl[l_test] = raw_prev_phase;

        mbe_synthesizeSpeechf(out, &cur, &prev);

        float wrapped_prev = fmodf(raw_prev_phase, (float)(2.0 * M_PI));
        if (wrapped_prev < 0.0f) {
            wrapped_prev += (float)(2.0 * M_PI);
        }
        float expected_psil = wrapped_prev + ((prev.w0 + cur.w0) * ((float)(l_test * 160) / 2.0f));

        assert(approx_equal(prev.PSIl[l_test], wrapped_prev, 1e-5f));
        assert(approx_equal(cur.PSIl[l_test], expected_psil, 1e-3f));
    }

    // A negative adaptive threshold attenuates a positive spectrum without inversion.
    {
        mbe_parms cur, prev;
        seed_speech_params(&cur, &prev);
        cur.L = 4;
        cur.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        for (int l = 1; l <= cur.L; ++l) {
            cur.Ml[l] = 10.0f;
            cur.Vl[l] = 1;
        }
        prev = cur;
        cur.errorRate = 0.5f;
        cur.errorCountTotal = 30;
        cur.errorCount4 = 2;
        prev.amplitudeThreshold = 1;

        mbe_parms smoothed = cur;
        mbe_applyAdaptiveSmoothing(&smoothed, &prev);
        float magnitude_sum = 0.0f;
        for (int l = 1; l <= smoothed.L; ++l) {
            assert(smoothed.Ml[l] >= 0.0f);
            magnitude_sum += smoothed.Ml[l];
        }
        assert(magnitude_sum <= 40.0f);
        assert(smoothed.amplitudeThreshold == -2999);

        // The shared synthesizer does not error-rate mute at the AMBE
        // threshold: the previous voiced frame must fade
        // toward zero, bounded by its linearly decreasing amplitude envelope.
        float out[160];
        mbe_synthesizeSpeechf(out, &cur, &prev);
        assert(approx_equal(out[0], 80.0f, 1e-5f));
        for (int n = 0; n < 160; ++n) {
            float envelope = 80.0f * (1.0f - (float)n / 160.0f);
            assert(fabsf(out[n]) <= envelope + 1e-4f);
        }
        for (int l = 1; l <= cur.L; ++l) {
            assert(cur.Ml[l] == 0.0f);
        }

        // Carry the signed threshold forward, then supply a clean positive
        // frame. The existing threshold reset must let speech fade back in.
        mbe_moveMbeParms(&cur, &prev);
        cur.errorRate = 0.0f;
        cur.errorCountTotal = 0;
        cur.errorCount4 = 0;
        for (int l = 1; l <= cur.L; ++l) {
            cur.Ml[l] = 10.0f;
        }
        mbe_synthesizeSpeechf(out, &cur, &prev);
        float recovery_energy = 0.0f;
        for (int n = 0; n < 160; ++n) {
            recovery_energy += out[n] * out[n];
        }
        assert(recovery_energy / 160.0f > 1.0f);
        for (int l = 1; l <= cur.L; ++l) {
            assert(cur.Ml[l] == 10.0f);
        }

        // At the adjacent positive threshold, degraded speech is attenuated,
        // not indiscriminately silenced (Tm = 6000 - 300 * 20 + 1 = 1).
        cur.errorRate = 0.5f;
        cur.errorCountTotal = 20;
        cur.errorCount4 = 2;
        prev.amplitudeThreshold = 1;
        mbe_applyAdaptiveSmoothing(&cur, &prev);
        magnitude_sum = 0.0f;
        for (int l = 1; l <= cur.L; ++l) {
            assert(cur.Ml[l] > 0.0f);
            assert(cur.Ml[l] < 10.0f);
            magnitude_sum += cur.Ml[l];
        }
        assert(approx_equal(magnitude_sum, 1.0f, 1e-6f));
    }

    // Analysis-windowed copies reconstruct every sample across the WOLA join.
    {
        float prev_uw[MBE_FFT_SIZE] = {0};
        float curr_uw[MBE_FFT_SIZE] = {0};
        float out[160] = {0};
        float expected[160];
        for (int n = 0; n < 160; ++n) {
            expected[n] = 0.1f + 0.001f * (float)n + 0.2f * sinf(0.13f * (float)n);
            if (n + 128 < MBE_FFT_SIZE) {
                prev_uw[n + 128] = mbe_synthesisWindow(n) * expected[n];
            }
            if (n >= 32) {
                curr_uw[n - 32] = mbe_synthesisWindow(n - 160) * expected[n];
            }
        }
        mbe_wola_combine(out, prev_uw, curr_uw, 160);
        for (int n = 0; n < 160; ++n) {
            assert(approx_equal(out[n], expected[n], 2e-6f * (1.0f + fabsf(expected[n]))));
        }
    }

    // Thread RNG seeding must drive comfort-noise and unvoiced-noise generators
    {
        float noise_a[160], noise_b[160], noise_c[160];
        float out[160];
        mbe_parms cur = {0}, prev = {0};

        mbe_setThreadRngSeed(0x12345678u);
        mbe_synthesizeComfortNoisef(noise_a);
        mbe_setThreadRngSeed(0x12345678u);
        mbe_synthesizeComfortNoisef(noise_b);
        assert(float_array_equal_exact(noise_a, noise_b, sizeof(noise_a) / sizeof(noise_a[0])) == 1);

        mbe_setThreadRngSeed(0x12340000u);
        mbe_synthesizeComfortNoisef(noise_c);
        assert(float_array_differs(noise_a, noise_c, sizeof(noise_a) / sizeof(noise_a[0])) == 1);

        seed_speech_params(&cur, &prev);
        cur.noiseSeed = -1.0f;
        memset(cur.noiseOverlap, 0, sizeof(cur.noiseOverlap));
        mbe_setThreadRngSeed(0x1234u);
        mbe_synthesizeSpeechf(out, &cur, &prev);
        assert(approx_equal(cur.noiseSeed, (float)(0x1234u % 53125u), 1e-6f));
    }

    // Fully voiced synthesis must not depend on the unvoiced-noise RNG seed.
    {
        float out[2][160];
        mbe_parms cur[2], prev[2], prev_enh;
        for (int i = 0; i < 2; ++i) {
            mbe_initMbeParms(&cur[i], &prev[i], &prev_enh);
            cur[i].w0 = 0.10f;
            cur[i].L = 12;
            for (int l = 1; l <= cur[i].L; ++l) {
                cur[i].Vl[l] = 1;
                cur[i].Ml[l] = 0.05f;
                cur[i].PHIl[l] = cur[i].PSIl[l] = 0.0f;
            }
            prev[i] = cur[i];
            mbe_setThreadRngSeed(i == 0 ? 0x1111u : 0x2222u);
            // Cold-start and second-frame noise heads are zero. The third
            // frame exercises seed-dependent noise at every harmonic index.
            for (int frame = 0; frame < 3; ++frame) {
                mbe_synthesizeSpeechf(out[i], &cur[i], &prev[i]);
                mbe_moveMbeParms(&cur[i], &prev[i]);
            }
        }
        assert(float_array_equal_exact(out[0], out[1], 160));
        assert(float_array_equal_exact(cur[0].PHIl, cur[1].PHIl, 57));
    }

    // Unit flat magnitudes have zero local phase; rising/falling edges lead/lag.
    {
        for (int edge = 0; edge < 3; ++edge) {
            float out[160];
            mbe_parms cur, prev, prev_enh;
            mbe_initMbeParms(&cur, &prev, &prev_enh);
            cur.w0 = 0.10f;
            cur.L = 56;
            for (int l = 1; l <= cur.L; ++l) {
                cur.Vl[l] = 1;
                cur.Ml[l] = ((edge == 1 && l > 28) || (edge == 2 && l <= 28)) ? 16.0f : 1.0f;
                cur.PHIl[l] = cur.PSIl[l] = 0.0f;
            }
            prev = cur;
            mbe_synthesizeSpeechf(out, &cur, &prev);
            if (edge == 0) {
                // Above harmonic 37 the kernel reaches the extrapolated tail.
                for (int l = 1; l <= 37; ++l) {
                    assert(float_bits_equal(cur.PHIl[l], cur.PSIl[l]));
                }
            } else if (edge == 1) {
                assert(cur.PHIl[28] > cur.PSIl[28]);
                assert(cur.PHIl[29] > cur.PSIl[29]);
            } else {
                assert(cur.PHIl[28] < cur.PSIl[28]);
                assert(cur.PHIl[29] < cur.PSIl[29]);
            }
        }
    }

    // AMBE prediction at the maximum previous harmonic must not read past it.
    {
        for (int mode = 0; mode < 2; ++mode) {
            char data[49] = {0};
            if (mode == 0) {
                // b0 = 36: sixteen harmonics, so 56 / L = 3.5 is exact.
                data[1] = data[4] = 1;
            } else {
                set_ambe2450_b0(data, 36);
            }
            float out[160];
            mbe_parms cur, prev, prev_enh;
            mbe_process_result result;
            mbe_initMbeParms(&cur, &prev, &prev_enh);
            int (*process)(float*, mbe_process_result*, const char*, mbe_parms*, mbe_parms*, mbe_parms*) =
                mode == 0 ? mbe_processAmbe2400Dataf : mbe_processAmbe2450Dataf;
            // Establish AMBE defaults before constructing the previous frame.
            mbe_initProcessResult(&result);
            assert(process(out, &result, data, &cur, &prev, &prev_enh) >= 0);
            prev.L = prev_enh.L = 56;
            prev.w0 = prev_enh.w0 = 0.05f;
            mbe_initProcessResult(&result);
            int status = process(out, &result, data, &cur, &prev, &prev_enh);
            assert(status >= 0);
            assert(cur.L < 56);
            for (int l = 1; l <= cur.L; ++l) {
                assert(isfinite(cur.log2Ml[l]));
            }
            for (int i = 0; i < 160; ++i) {
                assert(isfinite(out[i]));
            }
        }
    }

    // IMBE Dataf path must clear stale C4 error state when no C4 context is provided
    {
        char imbe_d[88];
        float out[160];
        mbe_process_result result;
        mbe_parms cur = {0}, prev = {0}, prev_enh = {0};

        mbe_initMbeParms(&cur, &prev, &prev_enh);
        set_bits_zero(imbe_d, 88);
        set_imbe7200_b0(imbe_d, 0);
        init_result_total(&result, 0);

        cur.errorCount4 = 7;
        prev.errorCount4 = 3;

        assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &prev_enh) >= 0);
        assert(cur.errorCount4 == 0);
    }

    // IMBE parameter synthesis must handle partial C0/C4 result metadata independently
    {
        char imbe_d[88];
        float out[160];
        mbe_process_result result;
        mbe_parms cur = {0}, prev = {0}, prev_enh = {0};

        mbe_initMbeParms(&cur, &prev, &prev_enh);
        set_bits_zero(imbe_d, 88);
        set_imbe7200_b0(imbe_d, 0);
        mbe_initProcessResult(&result);

        cur.errorCount4 = 7;
        result.flags = MBE_PROCESS_FLAG_C4_VALID;
        result.c4_errors = 3;
        result.protected_errors = 3;
        result.total_errors = 3;

        assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &prev_enh) >= 0);
        assert(cur.errorCount4 == 3);
    }

    {
        char imbe_d[88];
        float out[160];
        mbe_process_result result;
        mbe_parms cur = {0}, prev = {0}, prev_enh = {0};

        mbe_initMbeParms(&cur, &prev, &prev_enh);
        set_bits_zero(imbe_d, 88);
        set_imbe7200_b0(imbe_d, 0);
        mbe_initProcessResult(&result);

        cur.errorCount4 = 7;
        result.flags = MBE_PROCESS_FLAG_C0_VALID;
        result.c0_errors = 1;
        result.total_errors = 1;

        assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &prev_enh) >= 0);
        assert(cur.errorCount4 == 0);
    }

    // Result formatting must truncate large error counts to the caller-provided buffer
    {
        mbe_process_result result;
        char status[8];
        mbe_initProcessResult(&result);
        result.total_errors = 1000;
        result.flags = MBE_PROCESS_FLAG_REPEAT;
        mbe_formatProcessResult(status, sizeof(status), &result);
        assert(strlen(status) == sizeof(status) - 1u);
        assert(status[sizeof(status) - 1u] == '\0');
    }

    // IMBE repeat headroom parity: prolonged repeats reset to defaults rather than muting indefinitely
    {
        char imbe_d[88];
        float out[160];
        mbe_process_result result;
        mbe_parms cur = {0}, prev = {0}, prev_enh = {0};

        mbe_initMbeParms(&cur, &prev, &prev_enh);
        set_bits_zero(imbe_d, 88);
        set_imbe7200_b0(imbe_d, 0);
        init_result_total(&result, 6);

        prev.repeatCount = 4;
        cur = prev;

        assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &prev_enh) >= 0);

        assert(cur.repeatCount == 0);
        assert(cur.L == 39);
        float w0_default = (float)((4.0 * M_PI) / (134.0 + 39.5));
        assert(approx_equal(cur.w0, w0_default, 1e-5f));
        assert(cur.Vl[1] == 0);
        assert(cur.Ml[1] > 0.0f);
    }

    // IMBE processing must restore IMBE muting semantics after AMBE-style threshold contamination
    {
        char imbe_d[88];
        float out[160];
        mbe_process_result result;
        mbe_parms cur = {0}, prev = {0}, prev_enh = {0};

        mbe_initMbeParms(&cur, &prev, &prev_enh);
        set_bits_zero(imbe_d, 88);
        set_imbe7200_b0(imbe_d, 0);
        init_result_total(&result, 0);

        /* Simulate AMBE state reuse without reinit. */
        cur.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev_enh.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        prev.errorRate = 1.0f;
        cur.repeatCount = 0;

        float noise_seed_before = prev_enh.noiseSeed;
        assert(mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &prev_enh) >= 0);
        assert(float_bits_equal(prev_enh.noiseSeed, noise_seed_before));
    }

    return 0;
}
