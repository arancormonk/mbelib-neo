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

#include "ambe3600x2400_const.h"
#include "ambe3600x2450_const.h"
#include "ambe3600x2450_internal.h"
#include "mbe_adaptive.h"
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
 * Tone ID1 is U1[0..7] (ambe_d[12..19]). U1[8..11] (ambe_d[20..23]) carries a
 * redundant copy of ID(7..4) (TIA-102.BABA-1 Table 10) and does not change the
 * decoded ID.
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

/* Byte identity of two parameter sets, for exact reproducibility. */
static int
parms_identical(const mbe_parms* a, const mbe_parms* b) {
    unsigned char x[sizeof(mbe_parms)];
    unsigned char y[sizeof(mbe_parms)];
    memcpy(x, a, sizeof(x));
    memcpy(y, b, sizeof(y));
    return memcmp(x, y, sizeof(x)) == 0;
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

/* A flat, fully voiced 100 Hz model with harmonics to 3.6 kHz. */
static void
flat_model(mbe_parms* mp) {
    mbe_parms prev, enhanced;
    mbe_initMbeParms(mp, &prev, &enhanced);
    mp->w0 = (float)(2.0 * M_PI * 100.0 / 8000.0);
    mp->L = 36;
    for (int l = 1; l <= mp->L; ++l) {
        mp->Vl[l] = 1;
        mp->Ml[l] = 1.0f;
        mp->log2Ml[l] = 0.0f;
    }
}

/* Shelf shape: 0 dB to 2.5 kHz, rising linearly to the full gain at 3 kHz, flat to 3.6 kHz. */
static void
test_high_band_shelf_shape(void) {
    mbe_parms plain, shelved;
    flat_model(&plain);
    flat_model(&shelved);
    (void)mbe_spectralAmpEnhanceWithRm0(&plain, 1.0f);
    (void)mbe_spectralAmpEnhanceWithRm0(&shelved, MBE_HIGH_BAND_GAIN_AMBE2400);
    for (int l = 1; l <= 36; ++l) {
        float hz = 100.0f * (float)l;
        float t = hz <= 2500.0f ? 0.0f : (hz >= 3000.0f ? 1.0f : (hz - 2500.0f) / 500.0f);
        float want = 1.0f + ((MBE_HIGH_BAND_GAIN_AMBE2400 - 1.0f) * t);
        assert(approx_equal(shelved.Ml[l] / plain.Ml[l], want, 1e-5f));
    }
}

/* The low-band roll-off of DVSI's decoders: a second-order high-pass at 100 Hz. */
static void
test_low_band_gain(void) {
    static const float hz[] = {50.0f, 62.5f, 100.0f, 200.0f, 400.0f, 1000.0f, 3300.0f};
    for (size_t i = 0; i < sizeof(hz) / sizeof(hz[0]); ++i) {
        float want = 1.0f / sqrtf(1.0f + powf(100.0f / hz[i], 4.0f));
        assert(approx_equal(mbe_low_band_gain(hz[i]), want, 1e-5f));
    }
}

/* The roll-off's phase: an analog second-order Butterworth high-pass at
 * 100 Hz, H(s) = s^2 / (s^2 + sqrt(2) s + 1) at s = j hz / 100, leads by
 * 180 degrees toward DC, 90 at the corner and nothing far above it. */
static void
test_low_band_phase(void) {
    static const float hz[] = {20.0f, 62.5f, 100.0f, 200.0f, 400.0f, 1000.0f, 3000.0f};
    for (size_t i = 0; i < sizeof(hz) / sizeof(hz[0]); ++i) {
        double x = hz[i] / 100.0;
        double want = M_PI - atan2(sqrt(2.0) * x, 1.0 - (x * x));
        assert(fabs((double)mbe_low_band_phase(hz[i]) - want) < 1e-5);
    }
    assert(fabsf(mbe_low_band_phase(100.0f) - (float)(M_PI / 2.0)) < 1e-6f);
    assert(mbe_low_band_phase(0.0f) == 0.0f && mbe_low_band_phase(-5.0f) == 0.0f);
}

/* Amplitude of the component at `hz` over a buffer of n samples. */
static double
component_amplitude(const float* x, int n, double hz) {
    double re = 0.0, im = 0.0;
    for (int i = 0; i < n; ++i) {
        re += x[i] * cos(2.0 * M_PI * hz * i / 8000.0);
        im += x[i] * sin(2.0 * M_PI * hz * i / 8000.0);
    }
    return 2.0 * sqrt(re * re + im * im) / n;
}

/* Synthesis renders the roll-off without touching the model: a flat, fully
 * voiced 62.5 Hz series comes out with harmonic 1 at the roll-off's gain
 * relative to harmonic 8 (500 Hz), and the stored amplitudes stay flat. */
static void
test_low_band_rendering(void) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    cur.w0 = (float)(2.0 * M_PI * 62.5 / 8000.0);
    cur.L = 48;
    for (int l = 1; l <= cur.L; ++l) {
        cur.Vl[l] = 1;
        cur.Ml[l] = 1.0f;
    }
    static float out[8 * 160];
    for (int frame = 0; frame < 9; ++frame) {
        prev = cur;
        mbe_synthesizeSpeechf(frame ? out + (size_t)(frame - 1) * 160u : out, &cur, &prev);
    }
    for (int l = 1; l <= cur.L; ++l) {
        assert(cur.Ml[l] == 1.0f);
    }
    double ratio = component_amplitude(out, 8 * 160, 62.5) / component_amplitude(out, 8 * 160, 500.0);
    double want = (double)mbe_low_band_gain(62.5f) / (double)mbe_low_band_gain(500.0f);
    assert(fabs(20.0 * log10(ratio / want)) < 0.2);
}

/* Phase (cosine convention) of the component at `hz` over x[0..n), with the
 * time origin at sample `origin`. */
static double
component_phase(const float* x, int n, double hz, int origin) {
    double re = 0.0, im = 0.0;
    for (int i = 0; i < n; ++i) {
        double t = 2.0 * M_PI * hz * (i - origin) / 8000.0;
        re += x[i] * cos(t);
        im -= x[i] * sin(t);
    }
    return atan2(im, re);
}

static double
wrap_angle(double a) {
    return atan2(sin(a), cos(a));
}

/* Synthesis renders the roll-off's phase too, around the stored model: on a
 * steady, fully voiced 62.5 Hz series the output's harmonic phases at the
 * frame end lead the model's PHIl by the roll-off's phase (125 degrees at
 * 62.5 Hz, 4 at 2 kHz), and PHIl/PSIl stay what the model computed. */
static void
test_low_band_phase_rendering(void) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    cur.w0 = (float)(2.0 * M_PI * 62.5 / 8000.0);
    cur.L = 48;
    for (int l = 1; l <= cur.L; ++l) {
        cur.Vl[l] = 1;
        cur.Ml[l] = 1.0f;
    }
    static float out[8 * 160];
    for (int frame = 0; frame < 9; ++frame) {
        prev = cur;
        const mbe_parms before = prev;
        mbe_synthesizeSpeechf(frame ? out + (size_t)(frame - 1) * 160u : out, &cur, &prev);
        /* The previous model is rendered with the roll-off too, and restored. */
        assert(float_array_equal_exact(prev.Ml, before.Ml, 57));
        assert(float_array_equal_exact(prev.PHIl, before.PHIl, 57));
    }
    /* Two 128-sample periods ending at the last frame's end. */
    const float* tail = out + ((size_t)8 * 160u) - 256u;
    static const int ls[] = {1, 2, 3, 32};
    for (size_t i = 0; i < sizeof(ls) / sizeof(ls[0]); ++i) {
        const double hz = 62.5 * ls[i];
        double measured = component_phase(tail, 256, hz, 256);
        double want = (double)cur.PHIl[ls[i]] + (double)mbe_low_band_phase((float)hz);
        assert(fabs(wrap_angle(measured - want)) < 0.02);
    }
}

/* A caller may pass one model as both frames: it is prepared for rendering
 * once and restored, so its stored amplitudes and phases end as those of a
 * model synthesized from a separate copy of itself. */
static void
test_render_aliased_model(void) {
    mbe_parms model, ignored, enhanced;
    mbe_initMbeParms(&model, &ignored, &enhanced);
    model.w0 = (float)(2.0 * M_PI * 62.5 / 8000.0);
    model.L = 48;
    for (int l = 1; l <= model.L; ++l) {
        model.Vl[l] = (l % 5) != 0; /* voiced, with some unvoiced bands */
        model.Ml[l] = 1.0f / (float)l;
        model.PSIl[l] = 0.1f * (float)l;
        model.PHIl[l] = 0.2f * (float)l;
    }
    mbe_parms aliased = model, cur = model, prev = model;
    float out[160];
    mbe_setThreadRngSeed(0xA11A5u);
    mbe_synthesizeSpeechf(out, &aliased, &aliased);
    mbe_setThreadRngSeed(0xA11A5u);
    mbe_synthesizeSpeechf(out, &cur, &prev);
    assert(float_array_equal_exact(aliased.Ml, cur.Ml, 57));
    assert(float_array_equal_exact(aliased.Ml, model.Ml, 57));
    assert(float_array_equal_exact(aliased.PHIl, cur.PHIl, 57));
    assert(float_array_equal_exact(aliased.PSIl, cur.PSIl, 57));
}

/* Synthesis renders voiced harmonics 0.3 dB up and unvoiced ones 0.4 dB down,
 * shared by every codec: after the decoder fixes, DVSI's decoders play voiced
 * frames 0.32-0.38 dB louder and unvoiced frames 0.32-0.40 dB quieter than
 * this library did on identical bits. A flat voiced harmonic renders at
 * 2 Ml |H| / 7 times the voiced gain; a flat unvoiced model's level, measured
 * before the gains (1.442691668 RMS for this model and seed), times the
 * unvoiced gain. */
static void
test_render_voicing_gains(void) {
    const double gv = pow(10.0, 0.3 / 20.0), gu = pow(10.0, -0.4 / 20.0);
    static float out[16 * 160];
    for (int voiced = 0; voiced < 2; ++voiced) {
        mbe_parms cur, prev, enhanced;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        mbe_setThreadRngSeed(0x12345678u);
        cur.w0 = (float)(2.0 * M_PI * 62.5 / 8000.0);
        cur.L = 48;
        for (int l = 1; l <= cur.L; ++l) {
            cur.Vl[l] = voiced;
            cur.Ml[l] = 1.0f;
        }
        for (int frame = 0; frame < 17; ++frame) {
            prev = cur;
            mbe_synthesizeSpeechf(frame ? out + (size_t)(frame - 1) * 160u : out, &cur, &prev);
        }
        if (voiced) {
            double a = component_amplitude(out, 16 * 160, 500.0);
            assert(fabs(a / (2.0 * mbe_low_band_gain(500.0f) / 7.0) - gv) < 1e-4);
        } else {
            double e = 0.0;
            for (int n = 0; n < 16 * 160; ++n) {
                e += (double)out[n] * out[n];
            }
            assert(fabs(sqrt(e / (16 * 160)) / 1.442691668 - gu) < 1e-4);
        }
    }
}

typedef int (*process_fn)(float*, mbe_process_result*, const char*, mbe_parms*, mbe_parms*, mbe_parms*);

/* A codec's process path enhances with its own shelf: after one voice frame
 * from a fresh stream, the committed (unenhanced) model enhanced with `gain`
 * is the model that was synthesized. */
static void
expect_codec_shelf(process_fn process, const char* data, float gain) {
    mbe_parms cur, prev, enhanced;
    float out[160];
    mbe_initMbeParms(&cur, &prev, &enhanced);
    assert(process(out, NULL, data, &cur, &prev, &enhanced) >= 0);
    mbe_parms model = prev;
    (void)mbe_spectralAmpEnhanceWithRm0(&model, gain);
    assert(model.L == enhanced.L);
    for (int l = 1; l <= model.L; ++l) {
        assert(approx_equal(enhanced.Ml[l], model.Ml[l], 1e-4f * fmaxf(1.0f, fabsf(model.Ml[l]))));
    }
}

/* Each codec's high-band shelf height, matched to DVSI's decoders. */
static void
test_high_band_shelf_per_codec(void) {
    char imbe_d[88] = {0}, ambe_d[49] = {0};
    for (int i = 0; i < 88; ++i) {
        imbe_d[i] = (char)((i * 7 + 3) % 5 == 0);
    }
    for (int i = 0; i < 49; ++i) {
        ambe_d[i] = (char)((i * 5 + 1) % 3 == 0);
    }
    ambe_d[0] = 0; /* keep b0 a voice code */
    /* Heights in dB, each matched to DVSI's decoder for that mode (3-4 kHz
     * band within 0.2 dB on the development vectors). */
    assert(approx_equal(20.0f * log10f(MBE_HIGH_BAND_GAIN_IMBE), 2.0f, 1e-4f));
    assert(approx_equal(20.0f * log10f(MBE_HIGH_BAND_GAIN_AMBE2450), 2.5f, 1e-4f));
    assert(approx_equal(20.0f * log10f(MBE_HIGH_BAND_GAIN_AMBE2400), 2.5f, 1e-4f));
    expect_codec_shelf(mbe_processImbe4400Dataf, imbe_d, MBE_HIGH_BAND_GAIN_IMBE);
    expect_codec_shelf(mbe_processAmbe2450Dataf, ambe_d, MBE_HIGH_BAND_GAIN_AMBE2450);
    expect_codec_shelf(mbe_processAmbe2400Dataf, ambe_d, MBE_HIGH_BAND_GAIN_AMBE2400);
}

/* Direct callers may pass the same model as current and previous frame, and a
 * degenerate pitch: the roll-off must neither compound on the shared model nor
 * overflow (w0 of zero or near zero gives finite, quiet output). */
static void
test_low_band_edge_cases(void) {
    mbe_parms model, unused_prev, unused_enhanced;
    mbe_initMbeParms(&model, &unused_prev, &unused_enhanced);
    model.w0 = (float)(2.0 * M_PI * 62.5 / 8000.0);
    model.L = 48;
    for (int l = 1; l <= model.L; ++l) {
        model.Vl[l] = 1;
        model.Ml[l] = 1.0f;
    }
    float out[160];
    for (int call = 0; call < 3; ++call) {
        mbe_synthesizeSpeechf(out, &model, &model);
    }
    for (int l = 1; l <= model.L; ++l) {
        assert(model.Ml[l] == 1.0f);
    }

    static const float pitches[] = {0.0f, 1e-15f, 1e-6f};
    for (size_t i = 0; i < sizeof(pitches) / sizeof(pitches[0]); ++i) {
        mbe_parms cur, prev, enhanced;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        cur.L = 20;
        for (int l = 1; l <= cur.L; ++l) {
            cur.Vl[l] = 1;
            cur.Ml[l] = 1.0f;
        }
        prev = cur;
        cur.w0 = prev.w0 = pitches[i];
        mbe_synthesizeSpeechf(out, &cur, &prev);
        for (int n = 0; n < 160; ++n) {
            uint32_t bits;
            memcpy(&bits, &out[n], sizeof(bits));
            assert((bits & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000));
            assert(fabsf(out[n]) < 1.0f);
        }
    }
}

/* Rate-33 bit layout of DVSI's AMBE-3000: each of the 49 parameter bits lands
 * in exactly one quantizer value at the listed bit. DVSI sends b4's three low
 * bits at 40-42 and b3's LSB at 43; TIA-102.BABA-1 v1.0.5 Table 8 prints b3(0)
 * first. DVSI's rate-34 frames carry the same 49-bit vector as ours, so only
 * this field assignment differs. */
static void
test_ambe2450_layout(void) {
    static const signed char want[9][10] = {
        {0, 1, 2, 3, 37, 38, 39, -1},
        {4, 5, 6, 7, 35, -1},
        {8, 9, 10, 11, 36, -1},
        {12, 13, 14, 15, 16, 17, 18, 19, 43, -1},
        {20, 21, 22, 23, 40, 41, 42, -1},
        {24, 25, 26, 27, 44, -1},
        {28, 29, 30, 45, -1},
        {31, 32, 33, 46, -1},
        {34, 47, 48, -1},
    };
    int covered = 0;
    for (int i = 0; i < 9; ++i) {
        int width = 0;
        while (want[i][width] >= 0) {
            width++;
        }
        for (int j = 0; j < width; ++j) {
            char ambe_d[49] = {0};
            int b[9];
            ambe_d[(int)want[i][j]] = 1;
            mbe_ambe2450_unpack(ambe_d, b);
            for (int k = 0; k < 9; ++k) {
                assert(b[k] == (k == i ? 1 << (width - 1 - j) : 0));
            }
            int packed_b[9] = {0};
            char packed[49];
            memset(packed, 0, sizeof(packed));
            packed_b[i] = 1 << (width - 1 - j);
            mbe_ambe2450_pack(packed_b, packed);
            assert(memcmp(packed, ambe_d, sizeof(packed)) == 0);
            covered++;
        }
    }
    assert(covered == 49);
}

/* Spectral shape (log2 amplitudes less their mean) of a rate-33 voice frame with
 * a flat prediction history, rebuilt independently from TIA-102.BABA-1 eqs
 * 27-36 and the per-block inverse DCTs: G from the PRBA24 row b3 and the
 * PRBA58 row b4, R by the 8-point inverse DCT, the first two coefficients of
 * each block from R, the rest from the HOC rows b5..b8. */
static void
ambe2450_reference_shape(const int b[9], int L, double shape[57]) {
    double G[9] = {0}, R[9], C[5][18] = {{0}}, T[57] = {0};
    const float (*hoc[4])[4] = {AmbeHOCb5, AmbeHOCb6, AmbeHOCb7, AmbeHOCb8};
    for (int m = 0; m < 3; ++m) {
        G[2 + m] = AmbePRBA24[b[3]][m];
    }
    for (int m = 0; m < 4; ++m) {
        G[5 + m] = AmbePRBA58[b[4]][m];
    }
    for (int i = 1; i <= 8; ++i) {
        R[i] = 0.0;
        for (int m = 1; m <= 8; ++m) {
            R[i] += (m == 1 ? 1.0 : 2.0) * G[m] * cos(M_PI * (m - 1) * (i - 0.5) / 8.0);
        }
    }
    int l = 1;
    for (int i = 1; i <= 4; ++i) {
        const int J = AmbeLmprbl[L][i - 1];
        const int odd = (2 * i) - 1;
        C[i][1] = 0.5 * (R[odd] + R[odd + 1]);
        C[i][2] = (R[odd] - R[odd + 1]) / (2.0 * sqrt(2.0));
        for (int k = 3; k <= J && k <= 6; ++k) {
            C[i][k] = hoc[i - 1][b[4 + i]][k - 3];
        }
        for (int j = 1; j <= J; ++j, ++l) {
            for (int k = 1; k <= J; ++k) {
                T[l] += (k == 1 ? 1.0 : 2.0) * C[i][k] * cos(M_PI * (k - 1) * (j - 0.5) / J);
            }
        }
    }
    double mean = 0.0;
    for (int k = 1; k <= L; ++k) {
        mean += T[k] / L;
    }
    for (int k = 1; k <= L; ++k) {
        shape[k] = T[k] - mean;
    }
}

/* The decoder reads b3 and b4 from DVSI's positions (b3's LSB at 43, b4's low
 * bits at 40-42): frames built from those literal positions, with b3 and b4
 * values whose low bits differ, decode to the shape their PRBA rows give. */
static void
test_ambe2450_prba_reconstruction(void) {
    static const int b3_pos[9] = {12, 13, 14, 15, 16, 17, 18, 19, 43};
    static const int b4_pos[7] = {20, 21, 22, 23, 40, 41, 42};
    static const int cases[][2] = {{301, 72}, {300, 79}, {1, 6}, {510, 1}};
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
        int b[9] = {40, 0, 12, cases[c][0], cases[c][1], 3, 5, 9, 2};
        char ambe_d[49] = {0};
        set_ambe2450_b0(ambe_d, b[0]);
        for (int k = 0; k < 9; ++k) {
            ambe_d[b3_pos[k]] = (char)((b[3] >> (8 - k)) & 1);
        }
        for (int k = 0; k < 7; ++k) {
            ambe_d[b4_pos[k]] = (char)((b[4] >> (6 - k)) & 1);
        }
        static const int hoc_pos[4][5] = {
            {24, 25, 26, 27, 44}, {28, 29, 30, 45, -1}, {31, 32, 33, 46, -1}, {34, 47, 48, -1, -1}};
        static const int hoc_bits[4] = {5, 4, 4, 3};
        for (int i = 0; i < 4; ++i) {
            for (int k = 0; k < hoc_bits[i]; ++k) {
                ambe_d[hoc_pos[i][k]] = (char)((b[5 + i] >> (hoc_bits[i] - 1 - k)) & 1);
            }
        }
        mbe_parms cur, prev, dummy;
        mbe_initMbeParms(&cur, &prev, &dummy);
        assert(mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev) == 0);
        /* A second frame against the first's committed, flat-free history would
         * mix in prediction; decode once more from a flat history instead. */
        mbe_initMbeParms(&cur, &prev, &dummy);
        for (int l = 0; l <= 56; ++l) {
            prev.log2Ml[l] = 0.0f;
        }
        assert(mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev) == 0);
        double want[57], mean = 0.0;
        ambe2450_reference_shape(b, cur.L, want);
        for (int l = 1; l <= cur.L; ++l) {
            mean += cur.log2Ml[l] / cur.L;
        }
        for (int l = 1; l <= cur.L; ++l) {
            assert(fabs((cur.log2Ml[l] - mean) - want[l]) < 1e-4);
        }
    }
}

/* Rate-33 pitch law of DVSI's AMBE-3000 decoder and encoder: log2 f0 =
 * -4.322605 - 0.0220427 b0 and L = floor(0.9254 * 0.5 / f0). TIA-102.BABA-1
 * Annex A tabulates up to 0.2% higher f0 for some codes (b0 70-76); its L
 * column (below) follows the same floor rule from its own f0 and differs only
 * at b0 17, where DVSI renders a 12th harmonic. */
static void
test_ambe2450_pitch_law(void) {
    static const unsigned char annex_a_L[120] = {
        9,  9,  9,  9,  9,  9,  10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13,
        13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18,
        19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, 25, 26, 26, 26, 27,
        27, 28, 28, 29, 29, 30, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 36, 36, 37, 37, 38, 38, 39,
        40, 40, 41, 42, 42, 43, 43, 44, 45, 46, 46, 47, 48, 48, 49, 50, 51, 52, 52, 53, 54, 55, 56, 56};
    for (int b0 = 0; b0 < 120; ++b0) {
        char ambe_d[49];
        mbe_parms cur, prev, dummy;
        mbe_initMbeParms(&cur, &prev, &dummy);
        set_bits_zero(ambe_d, 49);
        set_ambe2450_b0(ambe_d, b0);
        assert(mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev) == 0);
        double f0 = pow(2.0, -4.322605 - (0.0220427 * b0));
        assert(fabs(((double)cur.w0 / (2.0 * M_PI * f0)) - 1.0) < 2e-6);
        assert(cur.L == (b0 == 17 ? 12 : annex_a_L[b0]));
    }
}

/* D-STAR's gain is absolute, as DVSI's AMBE-3000 decodes it: gamma =
 * 2 Dg[b2] - 0.09 whatever the previous frame's gamma (in DVSI's re-
 * quantizations of decoded D-STAR frames into rate 33 the implied gain sits a
 * median 0.037 log2 from this law, 0.35 from rate 33's 0.5 gamma(-1)
 * recursion). Every code, after histories from both ends of the gain range. */
static void
test_ambe2400_absolute_gain(void) {
    static const int b0_bits[7] = {0, 1, 2, 3, 4, 5, 48};
    static const int b2_bits[6] = {6, 7, 8, 9, 42, 43};
    for (int b2 = 0; b2 < 64; ++b2) {
        for (int history = 0; history < 2; ++history) {
            char ambe_d[49] = {0};
            for (int k = 0; k < 7; ++k) {
                ambe_d[b0_bits[k]] = (char)((30 >> (6 - k)) & 1); /* a voice pitch */
            }
            for (int k = 0; k < 6; ++k) {
                ambe_d[b2_bits[k]] = (char)((b2 >> (5 - k)) & 1);
            }
            mbe_parms cur, prev, dummy;
            mbe_initMbeParms(&cur, &prev, &dummy);
            prev.gamma = history ? 10.6f : -2.0f;
            assert(mbe_decodeAmbe2400Parms(ambe_d, &cur, &prev) == 0);
            assert(approx_equal(cur.gamma, (2.0f * AmbePlusDg[b2]) - 0.09f, 1e-5f));
        }
    }
}

/* Rate 33's AMBE+2 gain stays differential (TIA-102.BABA-1 eq 26). */
static void
test_ambe2450_differential_gain(void) {
    static const float prior[2] = {-2.0f, 7.5f};
    for (int b2 = 0; b2 < 32; ++b2) {
        for (int h = 0; h < 2; ++h) {
            char ambe_d[49] = {0};
            int b[9] = {40, 0, b2, 0, 0, 0, 0, 0, 0};
            mbe_ambe2450_pack(b, ambe_d);
            mbe_parms cur, prev, dummy;
            mbe_initMbeParms(&cur, &prev, &dummy);
            prev.gamma = prior[h];
            assert(mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev) == 0);
            assert(approx_equal(cur.gamma, AmbeDg[b2] + (0.5f * prior[h]), 1e-5f));
        }
    }
}

/* One D-STAR voice frame through the staged Data API with a caller-supplied
 * error count (and soft flag); returns the result flags. */
static unsigned
dstar_frame_with_errors(const char ambe_d[49], int errors, int soft, mbe_parms* cur, mbe_parms* prev, mbe_parms* enh,
                        float out[160]) {
    mbe_process_result result;
    init_result_total(&result, errors);
    if (soft) {
        result.flags |= MBE_PROCESS_FLAG_SOFT_INPUT;
    }
    assert(mbe_processAmbe2400Dataf(out, &result, ambe_d, cur, prev, enh) >= 0);
    return result.flags;
}

/* D-STAR mutes hard-decision streams whose error rate (0.95 eps_R + 0.001064
 * eps_T) exceeds 0.042, as DVSI's AMBE-3000 mutes D-STAR at sustained channel
 * error rates near 10% and never at 5%. Three errors a frame (no repeats)
 * take eps_R past 0.042 on the 21st frame, and clean frames bring it back.
 * Soft decisions never mute on the rate. A muted frame is synthesized as usual
 * and only its output replaced, so the hard stream's state (model, phases,
 * overlap and noise) stays byte for byte the soft stream's. */
static void
test_ambe2400_error_rate_mute(void) {
    enum { FRAMES = 52, ERRORED = 40 };

    static mbe_parms hard_state[FRAMES][3];
    static float hard_out[FRAMES][160];
    char ambe_d[49] = {0};
    static const int b0_bits[7] = {0, 1, 2, 3, 4, 5, 48};
    for (int k = 0; k < 7; ++k) {
        ambe_d[b0_bits[k]] = (char)((30 >> (6 - k)) & 1);
    }
    ambe_d[6] = 1; /* a nonzero gain */
    const float comfort = (0.003f * 32767.0f) / 7.0f;
    int muted = 0;
    for (int soft = 0; soft < 2; ++soft) {
        mbe_parms state[3];
        float out[160];
        memset(state, 0, sizeof(state));
        mbe_initMbeParms(&state[0], &state[1], &state[2]);
        mbe_setThreadRngSeed(0x2400u);
        float rate = 0.0f;
        for (int f = 0; f < FRAMES; ++f) {
            const int errors = (f < ERRORED) ? 3 : 0;
            rate = (0.95f * rate) + (0.001064f * (float)errors);
            unsigned flags = dstar_frame_with_errors(ambe_d, errors, soft, &state[0], &state[1], &state[2], out);
            assert((flags & MBE_PROCESS_FLAG_REPEAT) == 0u);
            assert(((flags & MBE_PROCESS_FLAG_MUTE) != 0u) == (!soft && rate > 0.042f));
            assert(f >= ERRORED || ((f >= 20) == (rate > 0.042f)));
            if (!soft) {
                memcpy(hard_state[f], state, sizeof(state));
                memcpy(hard_out[f], out, sizeof(out));
                continue;
            }
            for (int k = 0; k < 3; ++k) {
                assert(parms_identical(&hard_state[f][k], &state[k]));
            }
            if (rate > 0.042f) {
                /* Comfort noise where the soft stream plays the frame's speech. */
                float hard_peak = 0.0f, soft_peak = 0.0f;
                for (int n = 0; n < 160; ++n) {
                    hard_peak = fmaxf(hard_peak, fabsf(hard_out[f][n]));
                    soft_peak = fmaxf(soft_peak, fabsf(out[n]));
                }
                assert(hard_peak <= comfort && soft_peak > 10.0f * comfort);
                muted++;
            }
        }
        assert(rate < 0.042f);
    }
    assert(muted == 25); /* frames 21-40, and five clean ones while eps_R decays */
}

/* The repeat limit's mute resets D-STAR's model but keeps the error rate, so
 * a long run of bad frames keeps accumulating toward the error-rate mute. */
static void
test_ambe2400_repeat_limit_keeps_error_rate(void) {
    char ambe_d[49] = {0};
    ambe_d[6] = 1;
    mbe_parms cur, prev, enh;
    float out[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    double rate = 0.0;
    for (int f = 1; f <= 12; ++f) {
        (void)dstar_frame_with_errors(ambe_d, 5, 0, &cur, &prev, &enh, out); /* repeats: 5 > 3 */
        rate = (0.95 * rate) + (0.001064 * 5.0);
        assert(fabs((double)prev.errorRate - rate) < 1e-6);
    }
}

/**
 * @brief Test entry: IMBE and AMBE parameter derivation spot checks.
 */
int
main(void) {
    test_ambe2450_layout();
    test_ambe2450_prba_reconstruction();
    test_ambe2450_pitch_law();
    test_ambe2400_absolute_gain();
    test_ambe2450_differential_gain();
    test_ambe2400_error_rate_mute();
    test_ambe2400_repeat_limit_keeps_error_rate();
    test_low_band_gain();
    test_low_band_phase();
    test_low_band_rendering();
    test_low_band_phase_rendering();
    test_render_voicing_gains();
    test_render_aliased_model();
    test_low_band_edge_cases();
    test_high_band_shelf_shape();
    test_high_band_shelf_per_codec();

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
            assert(mbe_classifyAmbe2450Frame(ambe_d) == MBE_AMBE2450_FRAME_SILENCE);
            assert(mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev) == MBE_AMBE2450_FRAME_VOICE);
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

        /* Without C0 context the Dataf fallback keeps 2.1's rules: a verified
         * tone with fewer than 6 errors is output, other frames repeat above 3. */
        for (int total = 4; total <= 6; ++total) {
            mbe_initMbeParms(&cur, &prev, &prev_enh);
            init_result_total(&result, total);
            assert(mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh) >= 0);
            assert(result_has_marker(&result, 'T') == (total < 6));
            assert(result_has_marker(&result, 'R') == (total >= 6));
        }
        set_bits_zero(ambe_d, 49);
        set_ambe2450_b0(ambe_d, 10);
        mbe_initMbeParms(&cur, &prev, &prev_enh);
        init_result_total(&result, 4);
        assert(mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh) >= 0);
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
        /* TIA-102.BABA 7.8: [-5, 5] on the 16-bit output scale (float x 7). */
        double sumsq = 0.0;
        for (int i = 0; i < 160; ++i) {
            assert(fabsf(out[i]) * 7.0f <= 5.0f + 1e-4f);
            sumsq += (double)out[i] * (double)out[i] * 49.0;
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
        /* 2 * 10 per harmonic on s(n), rendered with the low-band roll-off's
         * magnitude and phase and the voiced gain (the model's phases are 0). */
        float rendered = 0.0f, first = 0.0f;
        for (int l = 1; l <= 4; ++l) {
            const float hz = prev.w0 * (8000.0f / (2.0f * (float)M_PI)) * (float)l;
            const float amplitude = 20.0f * mbe_low_band_gain(hz) * MBE_RENDER_VOICED_GAIN;
            rendered += amplitude;
            first += amplitude * cosf(mbe_low_band_phase(hz));
        }
        const float peak = rendered * MBE_SPEECH_OUTPUT_GAIN;
        assert(approx_equal(out[0], first * MBE_SPEECH_OUTPUT_GAIN, 1e-4f));
        for (int n = 0; n < 160; ++n) {
            float envelope = peak * (1.0f - (float)n / 160.0f);
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
