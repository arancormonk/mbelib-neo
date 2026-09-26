// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief AMBE 3600x2450 frame-type handling per TIA-102.BABA-1.
 *
 * Covers the Half-Rate Vocoder Addendum rules for silence frames (4.1, 4.3),
 * frame repeats (5.6), muting (5.7), erasures (4.1) and tone frames (7.3),
 * and the triplet
 * semantics they rely on: prev_mp is the last valid voice frame (prediction
 * history), prev_mp_enhanced is the last synthesized frame (repeat source).
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

/* Bit positions of quantizer values b0..b8 in the 49-bit parameter vector. */
static const int k_b0_bits[] = {0, 1, 2, 3, 37, 38, 39};
static const int k_b1_bits[] = {4, 5, 6, 7, 35};
static const int k_b2_bits[] = {8, 9, 10, 11, 36};
static const int k_b3_bits[] = {12, 13, 14, 15, 16, 17, 18, 19, 40};
static const int k_b4_bits[] = {20, 21, 22, 23, 41, 42, 43};
static const int k_b5_bits[] = {24, 25, 26, 27, 44};
static const int k_b6_bits[] = {28, 29, 30, 45};
static const int k_b7_bits[] = {31, 32, 33, 46};
static const int k_b8_bits[] = {34, 47, 48};

struct bfield {
    const int* bits;
    int count;
};

static const struct bfield k_fields[9] = {
    {k_b0_bits, 7}, {k_b1_bits, 5}, {k_b2_bits, 5}, {k_b3_bits, 9}, {k_b4_bits, 7},
    {k_b5_bits, 5}, {k_b6_bits, 4}, {k_b7_bits, 4}, {k_b8_bits, 3},
};

/** Pack quantizer values b0..b8 (MSB first) into a 49-bit parameter vector. */
static void
pack_b(char d[49], const int b[9]) {
    memset(d, 0, 49);
    for (int q = 0; q < 9; ++q) {
        for (int i = 0; i < k_fields[q].count; ++i) {
            int shift = k_fields[q].count - 1 - i;
            d[k_fields[q].bits[i]] = (char)((b[q] >> shift) & 1);
        }
    }
}

/* Loud voice frame (large gain delta), a second voice frame, and an erasure
 * whose u0 top six bits are not 63, so it cannot be classified as a tone. */
static const int k_voice_loud[9] = {40, 3, 31, 100, 40, 9, 5, 7, 2};
static const int k_voice_next[9] = {60, 5, 20, 200, 90, 17, 11, 3, 6};
static const int k_erasure[9] = {121, 0, 0, 0, 0, 0, 0, 0, 0};

/* Standard DVSI AMBE+2 silence vector 0xF801A99F8CE080 and its quantizer values. */
static const unsigned char k_dvsi_silence_bytes[7] = {0xF8, 0x01, 0xA9, 0x9F, 0x8C, 0xE0, 0x80};
static const int k_dvsi_silence[9] = {124, 16, 1, 53, 78, 18, 14, 12, 1};
static const float k_dg_b2_1 = -0.67f; /* Annex D gain level for b2 = 1 */

static void
unpack_bytes(char d[49], const unsigned char bytes[7]) {
    for (int i = 0; i < 49; ++i) {
        d[i] = (char)((bytes[i / 8] >> (7 - (i % 8))) & 1);
    }
}

/** Write the @p n low bits of @p value MSB first at d[pos..pos+n-1]. */
static void
put_bits(char d[49], int pos, int value, int n) {
    for (int i = 0; i < n; ++i) {
        d[pos + i] = (char)((value >> (n - 1 - i)) & 1);
    }
}

/**
 * Tone frame in the TIA-102.BABA-1 7.2 Table 10 layout, at full amplitude
 * (AD = 127): u0 = 63 then AD(6..1); u1 = ID(7..0), ID(7..4); u2 = ID(3..0),
 * ID(7..1); u3 = ID(0), ID(7..0), AD(0), 0000. Parameter bits 0-11 carry u0,
 * 12-23 u1, 24-34 u2 and 35-48 u3.
 */
static void
pack_tone(char d[49], int tone_id) {
    const int ad = 127;
    put_bits(d, 0, 63, 6);
    put_bits(d, 6, ad >> 1, 6);
    put_bits(d, 12, tone_id, 8);
    put_bits(d, 20, tone_id >> 4, 4);
    put_bits(d, 24, tone_id, 4);
    put_bits(d, 28, tone_id >> 1, 7);
    put_bits(d, 35, tone_id, 1);
    put_bits(d, 36, tone_id, 8);
    put_bits(d, 44, ad, 1);
    put_bits(d, 45, 0, 4);
}

static int
float_bits_equal(float a, float b) {
    uint32_t ab;
    uint32_t bb;
    memcpy(&ab, &a, sizeof(ab));
    memcpy(&bb, &b, sizeof(bb));
    return ab == bb;
}

/** Prediction history equality (fields, not bytes: decode writes the eq 44-45 edge values). */
static int
history_equal(const mbe_parms* a, const mbe_parms* b) {
    if (a->L != b->L || !float_bits_equal(a->gamma, b->gamma) || !float_bits_equal(a->w0, b->w0)) {
        return 0;
    }
    for (int l = 1; l <= a->L; ++l) {
        if (a->Vl[l] != b->Vl[l] || !float_bits_equal(a->log2Ml[l], b->log2Ml[l])
            || !float_bits_equal(a->Ml[l], b->Ml[l])) {
            return 0;
        }
    }
    return 1;
}

static int
float_array_bits_equal(const float* a, const float* b, int n) {
    for (int i = 0; i < n; ++i) {
        if (!float_bits_equal(a[i], b[i])) {
            return 0;
        }
    }
    return 1;
}

/** Whole parameter-set equality, compared per field (float bits, no padding bytes). */
static int
parms_equal(const mbe_parms* a, const mbe_parms* b) {
    int ints_equal = a->L == b->L && a->K == b->K && a->swn == b->swn && a->tonePhase == b->tonePhase
                     && a->amplitudeThreshold == b->amplitudeThreshold && a->errorCountTotal == b->errorCountTotal
                     && a->errorCount4 == b->errorCount4 && a->repeatCount == b->repeatCount
                     && memcmp(a->Vl, b->Vl, sizeof(a->Vl)) == 0;
    int floats_equal =
        float_bits_equal(a->w0, b->w0) && float_bits_equal(a->gamma, b->gamma)
        && float_bits_equal(a->localEnergy, b->localEnergy) && float_bits_equal(a->errorRate, b->errorRate)
        && float_bits_equal(a->mutingThreshold, b->mutingThreshold) && float_bits_equal(a->noiseSeed, b->noiseSeed);
    int arrays_equal = float_array_bits_equal(a->Ml, b->Ml, 57) && float_array_bits_equal(a->log2Ml, b->log2Ml, 57)
                       && float_array_bits_equal(a->PHIl, b->PHIl, 57) && float_array_bits_equal(a->PSIl, b->PSIl, 57)
                       && float_array_bits_equal(a->previousUw, b->previousUw, 256)
                       && float_array_bits_equal(a->noiseOverlap, b->noiseOverlap, 96);
    return ints_equal && floats_equal && arrays_equal;
}

/** After a mute, the last synthesized frame is silenced: amplitudes and overlap zeroed, the rest kept. */
static int
enhanced_silenced_from(const mbe_parms* after, const mbe_parms* before) {
    mbe_parms expected = *before;
    memset(expected.Ml, 0, sizeof(expected.Ml));
    memset(expected.previousUw, 0, sizeof(expected.previousUw));
    return parms_equal(after, &expected);
}

struct stream {
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms prev_enh;
    mbe_process_result result;
    float out[160];
};

static void
stream_init(struct stream* s) {
    mbe_setThreadRngSeed(0x5EEDu);
    mbe_initMbeParms(&s->cur, &s->prev, &s->prev_enh);
}

/** Process one frame with C0-valid error context (c0 + protected = total). */
static void
stream_frame(struct stream* s, const char d[49], int c0_errors, int protected_errors) {
    mbe_initProcessResult(&s->result);
    s->result.flags = MBE_PROCESS_FLAG_C0_VALID;
    s->result.c0_errors = c0_errors;
    s->result.protected_errors = protected_errors;
    s->result.total_errors = c0_errors + protected_errors;
    int rc = mbe_processAmbe2450Dataf(s->out, &s->result, d, &s->cur, &s->prev, &s->prev_enh);
    assert(rc >= 0);
    (void)rc;
}

static void
stream_frame_b(struct stream* s, const int b[9], int c0_errors, int protected_errors) {
    char d[49];
    pack_b(d, b);
    stream_frame(s, d, c0_errors, protected_errors);
}

static int
has_flag(const struct stream* s, unsigned flag) {
    return (s->result.flags & flag) != 0u;
}

static float
peak_abs(const float* x) {
    float m = 0.0f;
    for (int i = 0; i < 160; ++i) {
        m = fmaxf(m, fabsf(x[i]));
    }
    return m;
}

static double
rms(const float* x) {
    double acc = 0.0;
    for (int i = 0; i < 160; ++i) {
        acc += (double)x[i] * (double)x[i];
    }
    return sqrt(acc / 160.0);
}

/** Muted output is the 5.7 uniform noise in [-5, 5]. */
static void
assert_mute_noise(const float* x) {
    assert(peak_abs(x) <= 5.0f);
    double r = rms(x);
    assert(r > 1.5 && r < 4.0);
    (void)r;
}

/*
 * 4.1, 5.6, 5.7: erasures are frame repeats of the last synthesized frame;
 * the 4th consecutive invalid frame mutes instead of repeating. The prediction
 * history is never touched, so the next voice frame decodes as if the
 * erasures were absent.
 */
static void
test_erasure_run_repeats_then_mutes(void) {
    struct stream s;
    char erasure[49];
    pack_b(erasure, k_erasure);

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const mbe_parms history = s.prev;

    for (int i = 1; i <= 5; ++i) {
        const mbe_parms enh_before = s.prev_enh;
        stream_frame(&s, erasure, 0, 0);
        const int muted = (i >= 4);
        assert(has_flag(&s, MBE_PROCESS_FLAG_ERASURE));
        assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == muted);
        assert(s.cur.repeatCount == (muted ? MBE_MAX_FRAME_REPEATS : i));
        assert(s.prev.repeatCount == s.cur.repeatCount);
        assert(history_equal(&s.prev, &history));
        if (muted) {
            assert_mute_noise(s.out);
            assert(enhanced_silenced_from(&s.prev_enh, &enh_before));
        } else {
            /* A repeat replays the loud voice model, not comfort noise. */
            assert(peak_abs(s.out) > 5.0f);
            assert(s.cur.L == enh_before.L);
            assert(float_bits_equal(s.cur.w0, enh_before.w0));
        }
    }

    stream_frame_b(&s, k_voice_next, 0, 0);
    assert(s.result.flags == MBE_PROCESS_FLAG_C0_VALID);
    assert(s.cur.repeatCount == 0);

    struct stream ref;
    stream_init(&ref);
    stream_frame_b(&ref, k_voice_loud, 0, 0);
    stream_frame_b(&ref, k_voice_next, 0, 0);
    assert(history_equal(&s.prev, &ref.prev));
}

/*
 * 5.6: a repeat must advance this frame's error-rate recursion. The old path
 * copied the stale previous parameter set over it, freezing errorRate.
 */
static void
test_repeat_advances_error_rate(void) {
    struct stream s;
    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const float before = s.prev.errorRate;

    stream_frame_b(&s, k_voice_next, 4, 3);
    assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_ERASURE));
    const float expected = (0.95f * before) + (0.001064f * 7.0f);
    assert(fabsf(s.prev.errorRate - expected) < 1e-7f);
    assert(float_bits_equal(s.cur.errorRate, s.prev.errorRate));
    assert(s.cur.errorCountTotal == 7);
    (void)expected;
}

/*
 * 5.7: mute when the error rate exceeds 0.096. With 7 errors per frame the
 * recursion crosses the threshold on frame 21. A muted valid voice frame
 * still commits its prediction history (the encoder updated its predictor).
 */
static void
test_error_rate_mute_boundary(void) {
    struct stream s;
    char voice[49];
    pack_b(voice, k_voice_loud);
    stream_init(&s);

    for (int n = 1; n <= 21; ++n) {
        mbe_parms dec_cur = s.cur;
        mbe_parms dec_prev = s.prev;
        const mbe_parms enh_before = s.prev_enh;
        assert(mbe_decodeAmbe2450Parms(voice, &dec_cur, &dec_prev) == 0);

        stream_frame(&s, voice, 1, 6);
        assert(!has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == (n == 21));
        assert(history_equal(&s.prev, &dec_cur));
        if (n == 21) {
            assert_mute_noise(s.out);
            assert(enhanced_silenced_from(&s.prev_enh, &enh_before));
        }
    }
}

/*
 * 7.3: valid tones leave the prediction history untouched (4.3). Not in the
 * spec: the last synthesized frame is silenced, as after a mute.
 */
static void
test_valid_tone_keeps_state(void) {
    struct stream s;
    char tone[49];
    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const mbe_parms history = s.prev;
    const mbe_parms enh_before = s.prev_enh;

    /* The linked renderer is the oracle in both ordinary and NOTONES builds. */
    pack_tone(tone, 7);
    assert(mbe_classifyAmbe2450Frame(tone) == MBE_AMBE2450_FRAME_TONE);
    mbe_parms tone_state = s.cur;
    float expected[160];
    mbe_synthesizeTonef(expected, tone, &tone_state);
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(float_array_bits_equal(s.out, expected, 160));
    assert(history_equal(&s.prev, &history));
    assert(enhanced_silenced_from(&s.prev_enh, &enh_before));

    /* Table 9 / 7.3: ID 255 is a valid zero-amplitude tone. */
    pack_tone(tone, 255);
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(s.cur.repeatCount == 0);
    for (int i = 0; i < 160; ++i) {
        assert(float_bits_equal(s.out[i], 0.0f));
    }
}

/* 7.3: an invalid tone index is an erasure; at the 4th invalid frame it mutes without a reset. */
static void
test_invalid_tone_is_erasure(void) {
    struct stream s;
    char tone[49];
    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const mbe_parms history = s.prev;

    pack_tone(tone, 4);
    mbe_parms dec_cur = s.cur;
    mbe_parms dec_prev = s.prev;
    assert(mbe_classifyAmbe2450Frame(tone) == MBE_AMBE2450_FRAME_ERASURE);
    assert(mbe_decodeAmbe2450Parms(tone, &dec_cur, &dec_prev) == MBE_AMBE2450_FRAME_ERASURE);
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
    assert(has_flag(&s, MBE_PROCESS_FLAG_ERASURE));
    assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    assert(s.cur.repeatCount == 1);

    s.prev.repeatCount = MBE_MAX_FRAME_REPEATS - 1;
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    assert(s.cur.repeatCount == MBE_MAX_FRAME_REPEATS);
    assert_mute_noise(s.out);
    assert(history_equal(&s.prev, &history));
}

/*
 * 7: the u0 tone identifier makes a tone frame whatever b0 would read as. A
 * single tone with ID(6..4) = 5 reads as silence b0 125 (Tables 4 and 10), so
 * a tone frame failing the redundancy check is an erasure (7.3), not decoded
 * as silence.
 */
static void
test_unverified_tone_is_erasure(void) {
    struct stream s;
    char tone[49];
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms enh;

    pack_tone(tone, 80);
    assert(mbe_classifyAmbe2450Frame(tone) == MBE_AMBE2450_FRAME_TONE);
    tone[20] ^= 1; /* u1's copy of ID(7..4) disagrees */
    tone[48] = 1;  /* and the last four bits of u3 are not 0 */
    assert(mbe_classifyAmbe2450Frame(tone) == MBE_AMBE2450_FRAME_ERASURE);
    mbe_initMbeParms(&cur, &prev, &enh);
    assert(mbe_decodeAmbe2450Parms(tone, &cur, &prev) == MBE_AMBE2450_FRAME_ERASURE);

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const mbe_parms history = s.prev;
    stream_frame(&s, tone, 0, 0);
    assert(s.result.flags
           == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT));
    assert(history_equal(&s.prev, &history));
}

/* The frame classifier and the staged decode report the same frame types. */
static void
test_classify_matches_decode(void) {
    char d[49];
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms enh;

    pack_b(d, k_voice_next);
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_VOICE);
    unpack_bytes(d, k_dvsi_silence_bytes);
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_SILENCE);
    pack_b(d, k_erasure);
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_ERASURE);
    mbe_initMbeParms(&cur, &prev, &enh);
    assert(mbe_decodeAmbe2450Parms(d, &cur, &prev) == MBE_AMBE2450_FRAME_ERASURE);
    pack_tone(d, 255);
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_TONE);
    assert(mbe_decodeAmbe2450Parms(d, &cur, &prev) == MBE_AMBE2450_FRAME_TONE);

    d[0] = 2;
    assert(mbe_classifyAmbe2450Frame(d) == MBE_STATUS_INVALID_BITS);
    assert(mbe_classifyAmbe2450Frame(NULL) == MBE_STATUS_INVALID_ARGUMENT);
}

/* 5.6 criteria apply before tone classification; 5.7 muting applies to tones too. */
static void
test_tone_repeat_and_mute_order(void) {
    struct stream s;
    char tone[49];
    pack_tone(tone, 7);

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    stream_frame(&s, tone, 4, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_TONE));

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    s.prev.errorRate = 0.2f;
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    assert_mute_noise(s.out);
}

/* The consecutive-invalid count is caller-owned; clamp it instead of overflowing. */
static void
test_repeat_count_is_clamped(void) {
    struct stream s;
    char erasure[49];
    pack_b(erasure, k_erasure);

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    s.prev.repeatCount = INT_MAX;
    stream_frame(&s, erasure, 0, 0);
    assert(s.cur.repeatCount == MBE_MAX_FRAME_REPEATS);
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    s.prev.repeatCount = -5;
    stream_frame(&s, erasure, 0, 0);
    assert(s.cur.repeatCount == 1);
    assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
}

/* Erasure handling does not depend on a result context. */
static void
test_null_result(void) {
    struct stream s;
    char d[49];
    stream_init(&s);
    pack_b(d, k_voice_loud);
    assert(mbe_processAmbe2450Dataf(s.out, NULL, d, &s.cur, &s.prev, &s.prev_enh) >= 0);
    pack_b(d, k_erasure);
    assert(mbe_processAmbe2450Dataf(s.out, NULL, d, &s.cur, &s.prev, &s.prev_enh) >= 0);
    assert(s.cur.repeatCount == 1);
}

/* The three parameter sets must be distinct; aliasing is rejected before any state changes. */
static void
test_aliased_parms_rejected(void) {
    struct stream s;
    char d[49];
    stream_init(&s);
    pack_b(d, k_voice_loud);
    stream_frame(&s, d, 0, 0);
    const struct stream before = s;
    mbe_parms* const aliases[3][3] = {
        {&s.cur, &s.cur, &s.prev_enh},
        {&s.cur, &s.prev, &s.cur},
        {&s.cur, &s.prev, &s.prev},
    };
    for (int i = 0; i < 3; ++i) {
        assert(mbe_processAmbe2450Dataf(s.out, NULL, d, aliases[i][0], aliases[i][1], aliases[i][2])
               == MBE_STATUS_INVALID_ARGUMENT);
    }
    assert(parms_equal(&s.cur, &before.cur));
    assert(parms_equal(&s.prev, &before.prev));
    assert(parms_equal(&s.prev_enh, &before.prev_enh));
}

/* 4.1 eqs 1-3 on the standard DVSI silence vector: w0 = 2*pi/32, L = 14, unvoiced. */
static void
test_dvsi_silence_vector(void) {
    char d[49];
    char packed[49];
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms enh;

    unpack_bytes(d, k_dvsi_silence_bytes);
    pack_b(packed, k_dvsi_silence);
    assert(memcmp(d, packed, sizeof(d)) == 0);

    mbe_initMbeParms(&cur, &prev, &enh);
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_SILENCE);
    /* As in 2.1, the staged decode reports any decoded model as voice. */
    assert(mbe_decodeAmbe2450Parms(d, &cur, &prev) == MBE_AMBE2450_FRAME_VOICE);
    assert(cur.L == 14);
    assert(fabsf(cur.w0 - (float)(2.0 * M_PI / 32.0)) < 1e-6f);
    const float unvc = 0.2046f / sqrtf(cur.w0);
    for (int l = 1; l <= cur.L; ++l) {
        assert(cur.Vl[l] == 0);
        const float expected = unvc * exp2f(cur.log2Ml[l]);
        assert(fabsf(cur.Ml[l] - expected) <= 1e-5f * expected);
        (void)expected;
    }
    (void)unvc;
}

/*
 * 4.3, eq 26, eqs 43-44: silence frames are decoded against the last voice
 * frame's history and never replace it. A voice frame after a silence run
 * decodes exactly as if the silence frames were absent.
 */
static void
test_silence_freezes_history(void) {
    struct stream s;
    char silence[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const mbe_parms history = s.prev;
    const float silence_gamma = k_dg_b2_1 + (0.5f * history.gamma);

    for (int i = 0; i < 3; ++i) {
        stream_frame(&s, silence, 0, 0);
        assert(s.result.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_SILENCE));
        assert(fabsf(s.cur.gamma - silence_gamma) < 1e-5f);
        assert(history_equal(&s.prev, &history));
        assert(s.prev_enh.L == 14);
        assert(peak_abs(s.out) > 0.0f);
    }
    stream_frame_b(&s, k_voice_next, 0, 0);
    assert(s.result.flags == MBE_PROCESS_FLAG_C0_VALID);

    struct stream ref;
    stream_init(&ref);
    stream_frame_b(&ref, k_voice_loud, 0, 0);
    stream_frame_b(&ref, k_voice_next, 0, 0);
    assert(history_equal(&s.prev, &ref.prev));
    (void)silence_gamma;
}

/* A silence frame right after init leaves the spec initial history (L = 15, gamma = 0). */
static void
test_silence_first_after_init(void) {
    struct stream s;
    char silence[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);

    stream_init(&s);
    stream_frame(&s, silence, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
    assert(s.prev.L == 15);
    assert(float_bits_equal(s.prev.gamma, 0.0f));
    for (int l = 1; l <= s.prev.L; ++l) {
        assert(float_bits_equal(s.prev.log2Ml[l], 0.0f));
    }
    assert(fabsf(s.cur.gamma - k_dg_b2_1) < 1e-6f);
}

/* 5.6: a repeat or erasure after silence repeats the silence model, not the last voice frame. */
static void
test_repeat_after_silence_repeats_silence(void) {
    struct stream s;
    char silence[49];
    char erasure[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);
    pack_b(erasure, k_erasure);

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const mbe_parms history = s.prev;
    stream_frame(&s, silence, 0, 0);

    stream_frame_b(&s, k_voice_next, 4, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
    assert(s.cur.L == 14);
    assert(fabsf(s.cur.w0 - (float)(2.0 * M_PI / 32.0)) < 1e-6f);
    assert(history_equal(&s.prev, &history));

    stream_frame(&s, erasure, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_ERASURE));
    assert(s.cur.L == 14);
    assert(history_equal(&s.prev, &history));
}

/* A result carrying the new SILENCE flag is accepted when passed back in. */
static void
test_silence_flag_round_trip(void) {
    struct stream s;
    char silence[49];
    char with_flag[64];
    char without_flag[64];
    unpack_bytes(silence, k_dvsi_silence_bytes);

    stream_init(&s);
    stream_frame(&s, silence, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
    mbe_formatProcessResult(with_flag, sizeof(with_flag), &s.result);
    mbe_process_result plain = s.result;
    plain.flags &= ~MBE_PROCESS_FLAG_SILENCE;
    mbe_formatProcessResult(without_flag, sizeof(without_flag), &plain);
    assert(strcmp(with_flag, without_flag) == 0);

    assert(mbe_processAmbe2450Dataf(s.out, &s.result, silence, &s.cur, &s.prev, &s.prev_enh) >= 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
}

/*
 * Not in TIA-102.BABA-1: the first frame after a mute fades in instead of
 * overlapping the last frame synthesized before the mute. It starts from
 * exact silence and its head carries far less energy than a steady frame's,
 * whose head is dominated by the previous frame at full amplitude.
 */
static double
head_rms(const float* x) {
    double acc = 0.0;
    for (int i = 0; i < 32; ++i) {
        acc += (double)x[i] * (double)x[i];
    }
    return sqrt(acc / 32.0);
}

/**
 * Head RMS of a steady frame of k_voice_next (unclipped). The gain predictor
 * (eq 26) needs a few frames to settle from gamma = 0, so use the 8th frame.
 */
static double
steady_head_rms(void) {
    struct stream ref;
    stream_init(&ref);
    for (int i = 0; i < 8; ++i) {
        stream_frame_b(&ref, k_voice_next, 0, 0);
    }
    return head_rms(ref.out);
}

static void
test_fade_in_after_repeat_mute(void) {
    struct stream s;
    char erasure[49];
    pack_b(erasure, k_erasure);
    const double steady = steady_head_rms();

    stream_init(&s);
    for (int i = 0; i < 8; ++i) {
        stream_frame_b(&s, k_voice_next, 0, 0);
    }
    for (int i = 0; i < 4; ++i) {
        stream_frame(&s, erasure, 0, 0);
    }
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));

    stream_frame_b(&s, k_voice_next, 0, 0);
    assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    assert(float_bits_equal(s.out[0], 0.0f));
    assert(head_rms(s.out) < 0.5 * steady);
    assert(peak_abs(s.out) > 5.0f);
    (void)steady;
}

static void
test_fade_in_after_error_rate_mute(void) {
    struct stream s;
    char voice[49];
    pack_b(voice, k_voice_next);
    const double steady = steady_head_rms();

    stream_init(&s);
    for (int n = 1; n <= 21; ++n) {
        stream_frame(&s, voice, 1, 6);
    }
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));

    /* The error rate decays below 0.096 on the next clean frame. */
    stream_frame(&s, voice, 0, 0);
    assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    assert(float_bits_equal(s.out[0], 0.0f));
    assert(head_rms(s.out) < 0.5 * steady);
    (void)steady;
}

/*
 * 5.6 step 3: a repeat synthesizes the copied model unchanged. Smoothing it
 * again changed its voicing on this sequence before the fix; this frame's
 * error accounting and noise state must still advance.
 */
static void
test_repeat_preserves_model(void) {
    static const int voice[9] = {60, 6, 20, 200, 90, 17, 11, 3, 6};
    struct stream s;
    stream_init(&s);
    for (int i = 0; i < 8; ++i) {
        stream_frame_b(&s, voice, 0, 0);
    }
    const mbe_parms model = s.prev_enh;

    for (int r = 0; r < 2; ++r) {
        const float seed_before = s.prev_enh.noiseSeed;
        const float er_before = s.prev.errorRate;
        stream_frame_b(&s, voice, 4, 3);
        assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
        assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
        assert(s.cur.L == model.L);
        assert(float_bits_equal(s.cur.w0, model.w0));
        for (int l = 1; l <= model.L; ++l) {
            assert(s.cur.Vl[l] == model.Vl[l]);
            assert(float_bits_equal(s.cur.Ml[l], model.Ml[l]));
        }
        assert(fabsf(s.prev.errorRate - ((0.95f * er_before) + (0.001064f * 7.0f))) < 1e-7f);
        assert(!float_bits_equal(s.prev_enh.noiseSeed, seed_before));
        (void)er_before;
        (void)seed_before;
    }
}

/*
 * Speech after a tone must not overlap the speech from before it (not in
 * TIA-102.BABA-1): voice after the tone fades in, and a repeat right after it
 * replays silence. This holds with or without a mute before the tone, for the
 * zero-amplitude ID 255 and for an audible tone.
 */
static void
test_fade_in_after_tone(void) {
    static const int tone_ids[2] = {255, 7};
    const double steady = steady_head_rms();
    char erasure[49];
    char voice[49];
    char tone[49];
    pack_b(erasure, k_erasure);
    pack_b(voice, k_voice_next);

    for (int k = 0; k < 8; ++k) {
        const int tone_id = tone_ids[k & 1];
        const int repeat_after = (k >> 1) & 1;
        const int mute_before = (k >> 2) & 1;
        struct stream s;
        stream_init(&s);
        for (int i = 0; i < 8; ++i) {
            stream_frame(&s, voice, 0, 0);
        }
        for (int i = 0; mute_before && i < 4; ++i) {
            stream_frame(&s, erasure, 0, 0);
        }
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == mute_before);
        pack_tone(tone, tone_id);
        stream_frame(&s, tone, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
        assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
        if (!repeat_after) {
            stream_frame(&s, voice, 0, 0);
            assert(float_bits_equal(s.out[0], 0.0f));
            assert(head_rms(s.out) < 0.5 * steady);
        } else {
            stream_frame(&s, voice, 4, 0);
            assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
            assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
            assert(peak_abs(s.out) < 1.0f);
        }
    }
    (void)steady;
}

int
main(void) {
    test_repeat_preserves_model();
    test_fade_in_after_tone();
    test_fade_in_after_repeat_mute();
    test_fade_in_after_error_rate_mute();
    test_dvsi_silence_vector();
    test_silence_freezes_history();
    test_silence_first_after_init();
    test_repeat_after_silence_repeats_silence();
    test_silence_flag_round_trip();
    test_erasure_run_repeats_then_mutes();
    test_repeat_advances_error_rate();
    test_error_rate_mute_boundary();
    test_valid_tone_keeps_state();
    test_invalid_tone_is_erasure();
    test_unverified_tone_is_erasure();
    test_classify_matches_decode();
    test_tone_repeat_and_mute_order();
    test_repeat_count_is_clamped();
    test_null_result();
    test_aliased_parms_rejected();
    printf("test_ambe2450_frame_types: OK\n");
    return 0;
}
