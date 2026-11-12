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
 * history), prev_mp_enhanced is the last synthesized frame (repeat source),
 * and, where the spec is silent, what DVSI's decoder does: a mute or a tone
 * silences the output, not the decoder's model, which prev_mp_enhanced keeps.
 * Fades are checked against renderings of the same starting state (twin
 * streams, or the synthesis helpers called directly).
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

#include "ambe_common.h"
#include "mbe_adaptive.h"
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

/* Voice frames at one pitch (b0 60, a period of 50 samples) and one set of
 * amplitude codes: voiced in every 500 Hz column (row 0), voiced only below
 * 500 Hz (row 14), and unvoiced (row 16). */
static const int k_voice_reach[9] = {60, 0, 20, 200, 90, 17, 11, 3, 6};
static const int k_voice_low[9] = {60, 14, 20, 200, 90, 17, 11, 3, 6};
static const int k_voice_unvoiced[9] = {60, 16, 20, 200, 90, 17, 11, 3, 6};
#define K_VOICE_PERIOD 50

/* The decoder's mark of a silence frame, in the otherwise unused PHIl[0]: the
 * previous frame's model stands at mid-frame (stored as 3.0f). */
static const float k_mark_previous_held = 3.0f;

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

/* What the last frame left, in the otherwise unused Vl[0] of the last
 * synthesized frame: a tone sounding (with its AD and ID), and the next frame
 * synthesized starting from silence. */
#define W_TONE_HELD (1 << 16)
#define W_AFRESH    (1 << 17)

/* The word a played full-level tone (AD 127) of this ID leaves. */
static int
tone_word(int tone_id) {
    return W_TONE_HELD | W_AFRESH | (127 << 8) | tone_id;
}

/** After a tone, the last synthesized frame keeps its model and state: only
 * its unvoiced overlap is cleared, and element 0 of Vl holds the word. */
static int
kept_through_tone(const mbe_parms* after, const mbe_parms* before, int word) {
    mbe_parms expected = *before;
    memset(expected.previousUw, 0, sizeof(expected.previousUw));
    expected.Vl[0] = word;
    return parms_equal(after, &expected);
}

/** Gain of the incoming signal at sample n of a tone transition (80 samples, linear). */
static float
fade_in_gain(int n) {
    return (n >= 80) ? 1.0f : (float)n / 80.0f;
}

/** Whether `out` is `incoming` faded in over `outgoing` across the first 80 samples. */
static int
is_crossfade(const float* out, const float* incoming, const float* outgoing) {
    for (int n = 0; n < 160; ++n) {
        const float gain = fade_in_gain(n);
        const double want = ((double)gain * (double)incoming[n]) + ((1.0 - (double)gain) * (double)outgoing[n]);
        if (!(fabs((double)out[n] - want) <= 1e-4 * (1.0 + fabs(want)))) {
            return 0;
        }
    }
    return 1;
}

/** DVSI's unvoiced noise generator (x <- 173 x + 13849 mod 2^16, from 0 for a
 * fresh stream) `frames` frames on from `from`. */
static float
generator_after(float from, int frames) {
    uint32_t state = (from < 0.0f) ? 0u : (uint32_t)from;
    for (int i = 0; i < frames * 160; ++i) {
        state = (173u * state + 13849u) & 0xFFFFu;
    }
    return (float)state;
}

struct stream {
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms prev_enh;
    mbe_process_result result;
    float out[160];
};

/* Every frame starts the thread's comfort-noise generator from this seed, so
 * that a muted frame's noise is mute_noise_reference(). DVSI's unvoiced noise
 * generator, which lives in the stream's state, does not depend on it. */
#define K_NOISE_SEED 0x5EEDu

static void
stream_init(struct stream* s) {
    mbe_setThreadRngSeed(K_NOISE_SEED);
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
    mbe_setThreadRngSeed(K_NOISE_SEED);
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

/** The noise a muted frame plays: DVSI's mute level, uniform in
 * [-MBE_AMBE_MUTE_NOISE_AMPLITUDE, MBE_AMBE_MUTE_NOISE_AMPLITUDE] on the 16-bit
 * output scale (float samples times 7), the level of D-STAR's. */
static void
mute_noise_reference(float noise[160]) {
    mbe_setThreadRngSeed(K_NOISE_SEED);
    mbe_synthesizeUniformNoisef(noise, MBE_AMBE_MUTE_NOISE_AMPLITUDE * MBE_SPEECH_OUTPUT_GAIN);
}

/** Muted output that is the noise alone: bounded by the mute level, at its rms
 * (a uniform noise's is the amplitude over sqrt(3)). */
static void
assert_mute_noise(const float* x) {
    assert(peak_abs(x) * 7.0f <= MBE_AMBE_MUTE_NOISE_AMPLITUDE + 1e-4f);
    double r = rms(x) * 7.0 / MBE_AMBE_MUTE_NOISE_AMPLITUDE;
    assert(r > 0.4 && r < 0.733);
    (void)r;
}

/** The library's float output limit (95% of full scale), applied to a sum. */
static void
limit_sum(float* x, int n) {
    const float limit = (32767.0f * 0.95f) / 7.0f;
    for (int i = 0; i < n; ++i) {
        x[i] = (x[i] > limit) ? limit : ((x[i] < -limit) ? -limit : x[i]);
    }
}

/** A muted frame's second half is the mute noise alone: nothing of the last
 * frame reaches past its first 80 samples (16-bit scale, float times 7). */
static int
second_half_is_noise(const float* out) {
    for (int n = 80; n < 160; ++n) {
        if (!(fabsf(out[n]) * 7.0f <= MBE_AMBE_MUTE_NOISE_AMPLITUDE + 1e-4f)) {
            return 0;
        }
    }
    return 1;
}

static int
close_to(const float* out, const float* want, int from, int to) {
    for (int n = from; n < to; ++n) {
        if (!(fabs((double)out[n] - (double)want[n]) <= 1e-5 * (1.0 + fabs((double)want[n])))) {
            return 0;
        }
    }
    return 1;
}

/**
 * The first muted frame after a heard one, rendered independently of the
 * decoder from the last synthesized frame alone: a repeat of it toward a
 * silent copy of itself (its amplitudes zeroed, no frame mark), which fades it
 * out over the first subframe, as D-STAR's repeat-limit mute renders. `faded`
 * receives the state the rendering leaves (pulse state, phases, noise
 * generator).
 */
static void
fade_reference(float out[160], const mbe_parms* last_synthesized, mbe_parms* faded) {
    mbe_parms last = *last_synthesized;
    *faded = *last_synthesized;
    memset(faded->Ml, 0, sizeof(faded->Ml));
    mbe_set_frame_mark(faded, MBE_FRAME_MARK_NONE);
    faded->repeatCount = 0;
    mbe_synthesizeRepeatedSpeechf(out, faded, &last, MBE_VOICED_SYNTHESIS_DVSI);
}

/** Whether this build plays tones: one built with NOTONES renders them as
 * silence, and so do the references, which render through the library. */
static int
tones_sound(void) {
    char tone[49];
    float out[160];
    mbe_parms phases;
    memset(&phases, 0, sizeof(phases));
    pack_tone(tone, 7);
    mbe_synthesizeTonef(out, tone, &phases);
    return peak_abs(out) > 1000.0f;
}

/** A held tone's tail as a mute or the frame after the tone adds it: the tone
 * going on from the phases `phases` holds, times (1 - n / 80) for n < 80. */
static void
tail_reference(float tail[160], const char tone[49], const mbe_parms* phases) {
    mbe_parms state = *phases;
    mbe_synthesizeTonef(tail, tone, &state);
    for (int n = 0; n < 160; ++n) {
        tail[n] = (n < 80) ? (1.0f - ((float)n / 80.0f)) * tail[n] : 0.0f;
    }
}

/** The model fields of a parameter set (pitch, bands, voicing, amplitudes). */
static int
model_equal(const mbe_parms* a, const mbe_parms* b) {
    if (a->L != b->L || !float_bits_equal(a->w0, b->w0)) {
        return 0;
    }
    for (int l = 1; l <= a->L; ++l) {
        if (a->Vl[l] != b->Vl[l] || !float_bits_equal(a->Ml[l], b->Ml[l])) {
            return 0;
        }
    }
    return 1;
}

static int
overlap_cleared(const mbe_parms* mp) {
    for (int i = 0; i < 256; ++i) {
        if (!float_bits_equal(mp->previousUw[i], 0.0f)) {
            return 0;
        }
    }
    return 1;
}

/**
 * The state a muted frame leaves: the kept model (`model`), its unvoiced
 * overlap cleared, the next frame to start from silence, the noise generator
 * one frame on, nothing in cur_mp's word.
 */
static void
assert_muted_state(const struct stream* s, const mbe_parms* model, float generator_before) {
    assert(has_flag(s, MBE_PROCESS_FLAG_MUTE));
    assert(model_equal(&s->prev_enh, model));
    assert(overlap_cleared(&s->prev_enh));
    assert(s->prev_enh.Vl[0] == W_AFRESH);
    assert(s->cur.Vl[0] == 0);
    assert(float_bits_equal(s->prev_enh.noiseSeed, generator_after(generator_before, 1)));
    assert(float_bits_equal(s->cur.noiseSeed, s->prev_enh.noiseSeed));
}

/*
 * 4.1, 5.6, 5.7: erasures are frame repeats of the last synthesized frame;
 * the 4th consecutive invalid frame mutes instead of repeating. The prediction
 * history is never touched, so the next voice frame decodes as if the
 * erasures were absent. As DVSI's decoder plays it, the 4th plays the last
 * frame fading out over its first 80 samples under the mute noise and keeps
 * that model; the 5th is the noise alone.
 */
static void
test_erasure_run_repeats_then_mutes(void) {
    struct stream s;
    char erasure[49];
    float noise[160];
    pack_b(erasure, k_erasure);
    mute_noise_reference(noise);

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    const mbe_parms history = s.prev;

    for (int i = 1; i <= 5; ++i) {
        const mbe_parms enh_before = s.prev_enh;
        float fade[160];
        mbe_parms faded;
        fade_reference(fade, &enh_before, &faded);
        stream_frame(&s, erasure, 0, 0);
        const int muted = (i >= 4);
        assert(has_flag(&s, MBE_PROCESS_FLAG_ERASURE));
        assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == muted);
        assert(s.cur.repeatCount == (muted ? MBE_MAX_FRAME_REPEATS : i));
        assert(s.prev.repeatCount == s.cur.repeatCount);
        assert(history_equal(&s.prev, &history));
        assert(s.cur.Vl[0] == 0);
        if (i == 4) {
            float want[160];
            for (int n = 0; n < 160; ++n) {
                want[n] = fade[n] + noise[n];
            }
            limit_sum(want, 160);
            assert(peak_abs(fade) > 5.0f);
            assert(close_to(s.out, want, 0, 160));
            assert(second_half_is_noise(s.out));
            /* The repeated model is kept, with the pulses and phases the fade left. */
            assert_muted_state(&s, &enh_before, enh_before.noiseSeed);
            assert(float_array_bits_equal(s.prev_enh.noiseOverlap, faded.noiseOverlap, 96));
            assert(float_array_bits_equal(s.prev_enh.PSIl + 1, faded.PSIl + 1, 56));
            assert(float_array_bits_equal(s.prev_enh.PHIl + 1, faded.PHIl + 1, 56));
        } else if (i == 5) {
            /* A mute right after a mute is the noise alone. */
            assert(float_array_bits_equal(s.out, noise, 160));
            assert_mute_noise(s.out);
            assert_muted_state(&s, &enh_before, enh_before.noiseSeed);
            assert(float_array_bits_equal(s.prev_enh.noiseOverlap, enh_before.noiseOverlap, 96));
        } else {
            /* A repeat replays the loud voice model, not comfort noise. */
            assert(peak_abs(s.out) > 5.0f);
            assert(s.cur.L == enh_before.L);
            assert(float_bits_equal(s.cur.w0, enh_before.w0));
            assert(s.prev_enh.Vl[0] == 0);
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

/**
 * A valid frame `d` as the decoder synthesizes it from the state `before`:
 * decoded against the history, enhanced, smoothed with the frame's error
 * accounting (read from `after`, the stream once it took the frame) and
 * rendered as `start` says over the last synthesized frame. `model` receives
 * the state the synthesis leaves.
 */
static void
frame_reference(float out[160], mbe_parms* model, const struct stream* before, const char d[49],
                const struct stream* after, enum mbe_synthesis_start start) {
    mbe_parms prev = before->prev;
    mbe_parms enh = before->prev_enh;
    *model = before->cur;
    model->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
    model->errorRate = after->cur.errorRate;
    model->errorCountTotal = after->cur.errorCountTotal;
    model->errorCount4 = 0;
    model->repeatCount = 0;
    model->Vl[0] = 0;
    const int kind = mbe_decodeAmbe2450Parms(d, model, &prev);
    assert(kind == MBE_AMBE2450_FRAME_VOICE);
    (void)kind;
    const float rm0 = mbe_spectralAmpEnhanceWithRm0(model, MBE_HIGH_BAND_GAIN_AMBE2450);
    enh.Vl[0] = 0;
    mbe_synthesizeSpeechWithPreEnhRm0Startf(out, model, &enh, rm0, MBE_MUTE_NOISE_SPEC, MBE_VOICED_SYNTHESIS_DVSI,
                                            MBE_AMPLITUDE_LIMIT_NONE, start);
}

/* How the decoder starts the next frame it synthesizes from this state. */
static enum mbe_synthesis_start
start_of(const struct stream* s) {
    return ((s->prev_enh.Vl[0] & W_AFRESH) != 0) ? MBE_SYNTHESIS_FROM_SILENCE : MBE_SYNTHESIS_CONTINUES;
}

/** A valid frame's own model as an unheard synthesis of it leaves it. */
static void
own_model_reference(mbe_parms* model, const struct stream* before, const char d[49], const struct stream* after) {
    float scratch[160];
    frame_reference(scratch, model, before, d, after, start_of(before));
}

/** The kept model repeated from the state `before`, starting as `start` says. */
static void
repeat_reference(float out[160], const struct stream* before, enum mbe_synthesis_start start) {
    mbe_parms repeated = before->prev_enh;
    mbe_parms last = before->prev_enh;
    repeated.Vl[0] = 0;
    last.Vl[0] = 0;
    repeated.repeatCount = 0;
    mbe_synthesizeRepeatedSpeechStartf(out, &repeated, &last, MBE_VOICED_SYNTHESIS_DVSI, start);
}

/*
 * 5.7: mute when the error rate exceeds 0.096. With 7 errors per frame the
 * recursion crosses the threshold on frame 21. A muted valid voice frame
 * still commits its prediction history (the encoder updated its predictor),
 * and, as DVSI's decoder goes on, becomes the kept model, as an unheard
 * synthesis of it leaves it; the frame plays the last one fading out under
 * the mute noise.
 */
static void
test_error_rate_mute_boundary(void) {
    struct stream s;
    char voice[49];
    float noise[160];
    pack_b(voice, k_voice_loud);
    mute_noise_reference(noise);
    stream_init(&s);

    for (int n = 1; n <= 21; ++n) {
        mbe_parms dec_cur = s.cur;
        mbe_parms dec_prev = s.prev;
        const struct stream before = s;
        assert(mbe_decodeAmbe2450Parms(voice, &dec_cur, &dec_prev) == 0);

        stream_frame(&s, voice, 1, 6);
        assert(!has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == (n == 21));
        assert(history_equal(&s.prev, &dec_cur));
        if (n == 21) {
            float fade[160], want[160];
            mbe_parms faded, own;
            fade_reference(fade, &before.prev_enh, &faded);
            for (int i = 0; i < 160; ++i) {
                want[i] = fade[i] + noise[i];
            }
            limit_sum(want, 160);
            assert(close_to(s.out, want, 0, 160));
            assert(second_half_is_noise(s.out));
            own_model_reference(&own, &before, voice, &s);
            assert_muted_state(&s, &own, before.prev_enh.noiseSeed);
            assert(s.cur.repeatCount == 0);
            assert(float_bits_equal(s.prev_enh.localEnergy, own.localEnergy));
            assert(float_array_bits_equal(s.prev_enh.noiseOverlap, own.noiseOverlap, 96));
        }
    }
}

/*
 * 7.3: valid tones leave the prediction history untouched (4.3). Not in the
 * spec: the last synthesized frame is kept, save its unvoiced overlap, and the
 * frame after the tone starts from silence, as after a mute.
 */
static void
test_valid_tone_keeps_state(void) {
    struct stream s;
    char tone[49];
    char erasure[49];
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
    /* What a lost frame would repeat: the last frame, which fades out under
     * the tone's first 80 samples while the tone fades in. */
    struct stream repeat = s;
    pack_b(erasure, k_erasure);
    stream_frame(&repeat, erasure, 0, 0);
    assert(has_flag(&repeat, MBE_PROCESS_FLAG_REPEAT) && peak_abs(repeat.out) > 0.0f);
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(is_crossfade(s.out, expected, repeat.out));
    assert(float_array_bits_equal(s.out + 80, expected + 80, 80));
    assert(history_equal(&s.prev, &history));
    assert(kept_through_tone(&s.prev_enh, &enh_before, tone_word(7)));
    assert(s.cur.Vl[0] == 0);

    /* The tone goes on at once: its second frame is the tone alone. */
    mbe_synthesizeTonef(expected, tone, &tone_state);
    stream_frame(&s, tone, 0, 0);
    assert(float_array_bits_equal(s.out, expected, 160));

    /* The zero-amplitude ID 255 is a lost frame: the tone repeats, at once. */
    mbe_synthesizeTonef(expected, tone, &tone_state);
    pack_tone(tone, 255);
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
    assert(has_flag(&s, MBE_PROCESS_FLAG_ERASURE) && has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(s.cur.repeatCount == 1);
    assert(float_array_bits_equal(s.out, expected, 160));
    assert(s.prev_enh.Vl[0] == tone_word(7));
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
    const mbe_parms enh_before = s.prev_enh;
    float fade[160], noise[160], want[160];
    mbe_parms faded;
    fade_reference(fade, &enh_before, &faded);
    mute_noise_reference(noise);
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    assert(s.cur.repeatCount == MBE_MAX_FRAME_REPEATS);
    /* The repeat it would have played fades out under the mute noise. */
    for (int n = 0; n < 160; ++n) {
        want[n] = fade[n] + noise[n];
    }
    limit_sum(want, 160);
    assert(close_to(s.out, want, 0, 160));
    assert(second_half_is_noise(s.out));
    assert_muted_state(&s, &enh_before, enh_before.noiseSeed);
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
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_ERASURE);
    assert(mbe_decodeAmbe2450Parms(d, &cur, &prev) == MBE_AMBE2450_FRAME_ERASURE);
    pack_tone(d, 7);
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

    /* A tone muted by the error rate: the voice before it fades out under the
     * noise, the tone does not sound (its phases stand), the voice model is
     * kept and the noise generator runs. */
    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    s.prev.errorRate = 0.2f;
    const struct stream before = s;
    float fade[160], noise[160], want[160];
    mbe_parms faded;
    fade_reference(fade, &before.prev_enh, &faded);
    mute_noise_reference(noise);
    stream_frame(&s, tone, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    for (int n = 0; n < 160; ++n) {
        want[n] = fade[n] + noise[n];
    }
    limit_sum(want, 160);
    assert(close_to(s.out, want, 0, 160));
    assert(second_half_is_noise(s.out));
    assert_muted_state(&s, &before.prev_enh, before.prev_enh.noiseSeed);
    assert(s.cur.repeatCount == 0);
    assert(s.cur.swn == before.cur.swn && s.cur.tonePhase == before.cur.tonePhase);
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

/*
 * A silence frame's decoded amplitudes: every band unvoiced, Ml the voice law's
 * unvoiced amplitude 0.2046 / sqrt(w0) * 2^log2Ml a further 0.2046 down
 * (DVSI's decoder; TIA-102.BABA-1 4.1 gives no extra factor).
 */
static int
silence_amplitudes_ok(const mbe_parms* cur) {
    const float unvc = 0.2046f / sqrtf(cur->w0);
    for (int l = 1; l <= cur->L; ++l) {
        const double expected = 0.2046 * (double)unvc * exp2((double)cur->log2Ml[l]);
        if (cur->Vl[l] != 0 || !(fabs((double)cur->Ml[l] - expected) <= 1e-5 * expected)) {
            return 0;
        }
    }
    return 1;
}

/* The standard DVSI silence vector: w0 = 2*pi/32 (4.1 eq 1), and, as DVSI's
 * decoder plays it, 15 unvoiced bands 13.8 dB below the voice law. */
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
    assert(cur.L == 15);
    assert(fabsf(cur.w0 - (float)(2.0 * M_PI / 32.0)) < 1e-6f);
    assert(silence_amplitudes_ok(&cur));
    assert(float_bits_equal(cur.PHIl[0], k_mark_previous_held));
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
        assert(s.prev_enh.L == 15);
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

/* Silence frames b0 124 and 125, from the initial state and after a voice
 * frame: L 15 at w0 = 2*pi/32, unvoiced, 13.8 dB below the voice law, marked;
 * the staged decode reports them as voice and leaves the history as it was. */
static void
test_silence_model_both_codes(void) {
    char d[49];
    unpack_bytes(d, k_dvsi_silence_bytes);
    for (int b0 = 124; b0 <= 125; ++b0) {
        for (int after_voice = 0; after_voice <= 1; ++after_voice) {
            mbe_parms cur, prev, enh;
            mbe_initMbeParms(&cur, &prev, &enh);
            if (after_voice) {
                char voice[49];
                pack_b(voice, k_voice_loud);
                const int voice_kind = mbe_decodeAmbe2450Parms(voice, &cur, &prev);
                assert(voice_kind == MBE_AMBE2450_FRAME_VOICE);
                (void)voice_kind;
                prev = cur;
            }
            const mbe_parms history = prev;
            put_bits(d, 0, b0 >> 3, 4);
            put_bits(d, 37, b0, 3);
            assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_SILENCE);
            const int kind = mbe_decodeAmbe2450Parms(d, &cur, &prev);
            assert(kind == MBE_AMBE2450_FRAME_VOICE);
            (void)kind;
            assert(cur.L == 15);
            assert(float_bits_equal(cur.w0, (float)(2.0 * M_PI / 32.0)));
            assert(silence_amplitudes_ok(&cur));
            assert(float_bits_equal(cur.PHIl[0], k_mark_previous_held));
            assert(history_equal(&prev, &history));
        }
    }
}

static double
span_energy(const float* x, int from, int to) {
    double acc = 0.0;
    for (int n = from; n < to; ++n) {
        acc += (double)x[n] * (double)x[n];
    }
    return acc;
}

/* A settled stream, then `before`, then `silences` silence frames. */
static void
stream_into_silence(struct stream* s, const int before[9], int silences) {
    char silence[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);
    stream_init(s);
    for (int i = 0; i < 8; ++i) {
        stream_frame_b(s, k_voice_next, 0, 0);
    }
    stream_frame_b(s, before, 0, 0);
    for (int i = 0; i < silences; ++i) {
        stream_frame(s, silence, 0, 0);
        assert(has_flag(s, MBE_PROCESS_FLAG_SILENCE));
    }
}

/*
 * As DVSI's decoder plays a silence frame, the previous frame's model stands
 * at mid-frame: the first 10 ms play that model held steady, and the second
 * crossfade from it to the silence model. A twin stream whose next frame is
 * an erasure plays the same model held through the whole frame (a repeat),
 * from the same pulse train and noise generator state:
 *  - samples 0-79 of the silence frame equal the repeat's, whatever the voicing
 *    of the frame before (every column, below 500 Hz only, none), within a
 *    float tolerance (the repeat's mid-frame model passes through the log
 *    domain; the silence frame's is the previous model itself);
 *  - after a fully voiced frame, which leaves no unvoiced band, samples 80-159
 *    are the repeat's faded out linearly, plus the silence model's own noise
 *    rising, far below them.
 */
static void
test_silence_holds_previous_model(void) {
    const int* before[3] = {k_voice_reach, k_voice_low, k_voice_unvoiced};
    char silence[49], erasure[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);
    pack_b(erasure, k_erasure);

    for (int k = 0; k < 3; ++k) {
        struct stream s;
        stream_into_silence(&s, before[k], 0);
        struct stream held = s;
        stream_frame(&s, silence, 0, 0);
        stream_frame(&held, erasure, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
        assert(has_flag(&held, MBE_PROCESS_FLAG_REPEAT) && !has_flag(&held, MBE_PROCESS_FLAG_MUTE));

        assert(span_energy(held.out, 0, 80) > 80.0 * 50.0 * 50.0);
        for (int n = 0; n < 80; ++n) {
            assert(fabs((double)s.out[n] - (double)held.out[n]) <= 1e-3 * (1.0 + fabs((double)held.out[n])));
        }
        if (k == 0) {
            double residual = 0.0;
            for (int n = 80; n < 160; ++n) {
                const double gain = 1.0 - ((double)(n - 80) / 80.0);
                const double r = (double)s.out[n] - (gain * (double)held.out[n]);
                residual += r * r;
            }
            assert(residual < 1e-4 * span_energy(held.out, 80, 160));
            (void)residual;
        }
    }
}

/*
 * The voiced frame after ONE silence frame continues the pulse train that ran
 * through it (the silence frame's mark carries it): it plays otherwise than a
 * twin whose silence frame lost its mark, which starts the train afresh. That
 * holds whether the frame before the silence frame was voiced in every column
 * or only below 500 Hz. After TWO silence frames the train has stopped (the
 * second has no voiced band at either end of either half), and the voiced
 * frame plays exactly as the twin that starts afresh; after an unvoiced frame
 * no train ran at all.
 */
static double
voiced_after_silence(const int before[9], int silences, float* out, float* afresh) {
    struct stream s;
    stream_into_silence(&s, before, silences);
    struct stream twin = s;
    twin.prev_enh.PHIl[0] = 0.0f;
    stream_frame_b(&s, k_voice_next, 0, 0);
    stream_frame_b(&twin, k_voice_next, 0, 0);
    memcpy(out, s.out, sizeof(s.out));
    memcpy(afresh, twin.out, sizeof(twin.out));
    double diff = 0.0;
    for (int n = 0; n < 160; ++n) {
        diff += ((double)out[n] - (double)afresh[n]) * ((double)out[n] - (double)afresh[n]);
    }
    return diff;
}

static void
test_pulse_train_through_one_silence_frame(void) {
    float out[160], afresh[160];
    const int* voiced[2] = {k_voice_reach, k_voice_low};
    for (int k = 0; k < 2; ++k) {
        /* The pulses the silence frames leave (their count in noiseOverlap[0]). */
        struct stream s;
        stream_into_silence(&s, voiced[k], 1);
        assert(s.prev_enh.noiseOverlap[0] >= 1.0f);
        stream_into_silence(&s, voiced[k], 2);
        assert(float_bits_equal(s.prev_enh.noiseOverlap[0], 0.0f));
        const double one = voiced_after_silence(voiced[k], 1, out, afresh);
        assert(one > 0.1 * span_energy(out, 0, 160));
        (void)one;
        (void)voiced_after_silence(voiced[k], 2, out, afresh);
        assert(float_array_bits_equal(out, afresh, 160));
    }
    (void)voiced_after_silence(k_voice_unvoiced, 1, out, afresh);
    assert(float_array_bits_equal(out, afresh, 160));
}

/* Consecutive silence frames: each plays the 15-band model at its level and
 * the history stays the last voice frame's. */
static void
test_consecutive_silence_frames(void) {
    struct stream s;
    stream_into_silence(&s, k_voice_reach, 0);
    const mbe_parms history = s.prev;
    char silence[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);
    double level[4];
    for (int i = 0; i < 4; ++i) {
        stream_frame(&s, silence, 0, 0);
        assert(s.result.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_SILENCE));
        assert(s.prev_enh.L == 15);
        assert(float_bits_equal(s.prev_enh.PHIl[0], k_mark_previous_held));
        assert(history_equal(&s.prev, &history));
        level[i] = rms(s.out);
        /* From the second on, nothing voiced is left to run on. */
        assert(i == 0 || float_bits_equal(s.prev_enh.noiseOverlap[0], 0.0f));
    }
    (void)history;
    /* The first fades the voice frame out; the others hold the silence level. */
    assert(level[1] < 0.05 * level[0]);
    for (int i = 2; i < 4; ++i) {
        assert(level[i] > 0.5 * level[1] && level[i] < 2.0 * level[1]);
    }
    (void)level;
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

/* 5.6: a repeat or erasure after silence repeats the 15-band silence model, not the last voice frame. */
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
    const mbe_parms silent = s.prev_enh;

    stream_frame_b(&s, k_voice_next, 4, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(!has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
    assert(s.cur.L == 15);
    assert(fabsf(s.cur.w0 - (float)(2.0 * M_PI / 32.0)) < 1e-6f);
    assert(history_equal(&s.prev, &history));
    for (int l = 1; l <= 15; ++l) {
        assert(s.cur.Vl[l] == 0);
        assert(float_bits_equal(s.cur.Ml[l], silent.Ml[l]));
    }

    stream_frame(&s, erasure, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_ERASURE));
    assert(s.cur.L == 15);
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
 * The frame after a tone starts from silence (not in TIA-102.BABA-1): it is
 * synthesized from the speech model kept through the tone, as the frame after
 * a mute is, and the tone goes on from its phases, fading out over the first
 * 80 samples, added on top; the sum is limited. This holds with or without a
 * mute before the tone, and after lost frames that repeated the tone (whose
 * phases the tail continues from). The frame after it continues as any does.
 */
static void
test_frame_after_tone(void) {
    char erasure[49];
    char voice[49];
    char tone[49];
    pack_b(erasure, k_erasure);
    pack_b(voice, k_voice_next);
    pack_tone(tone, 7);

    for (int k = 0; k < 4; ++k) {
        const int lost_after = k & 1;
        const int mute_before = (k >> 1) & 1;
        struct stream s, released;
        stream_init(&s);
        for (int i = 0; i < 8; ++i) {
            stream_frame(&s, voice, 0, 0);
        }
        for (int i = 0; mute_before && i < 4; ++i) {
            stream_frame(&s, erasure, 0, 0);
        }
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == mute_before);
        const mbe_parms kept = s.prev_enh;
        stream_frame(&s, tone, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
        assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
        assert(model_equal(&s.prev_enh, &kept));
        for (int i = 0; lost_after && i < 2; ++i) {
            stream_frame(&s, erasure, 0, 0);
            assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT) && s.prev_enh.Vl[0] == tone_word(7));
        }
        float tail[160], want[160];
        mbe_parms model;
        tail_reference(tail, tone, &s.cur);
        const struct stream before = s;
        stream_frame(&s, voice, 0, 0);
        assert(s.result.flags == MBE_PROCESS_FLAG_C0_VALID);
        frame_reference(want, &model, &before, voice, &s, MBE_SYNTHESIS_FROM_SILENCE);
        for (int n = 0; n < 160; ++n) {
            want[n] += tail[n];
        }
        limit_sum(want, 160);
        assert(close_to(s.out, want, 0, 160));
        assert(!tones_sound() || peak_abs(tail) > 1000.0f);
        assert(model_equal(&s.prev_enh, &model));
        assert(s.prev_enh.Vl[0] == 0 && s.cur.Vl[0] == 0);
        assert(s.cur.repeatCount == 0);
        released = s;
        stream_frame(&s, voice, 0, 0);
        frame_reference(want, &model, &released, voice, &s, MBE_SYNTHESIS_CONTINUES);
        assert(close_to(s.out, want, 0, 160));
    }
}

/* ---- Mutes, as DVSI's decoder plays them ---- */

/* A settled stream ending with a frame of `last`. */
static void
stream_settled(struct stream* s, const int last[9]) {
    stream_init(s);
    for (int i = 0; i < 8; ++i) {
        stream_frame_b(s, k_voice_next, 0, 0);
    }
    stream_frame_b(s, last, 0, 0);
}

/* The output a first muted frame plays over the last synthesized frame
 * `last`: its fade plus the mute noise. */
static void
faded_mute_reference(float want[160], const mbe_parms* last, mbe_parms* faded) {
    float noise[160];
    fade_reference(want, last, faded);
    mute_noise_reference(noise);
    for (int n = 0; n < 160; ++n) {
        want[n] += noise[n];
    }
    limit_sum(want, 160);
}

/*
 * The fade, against a twin stream that repeats the last frame instead: the
 * last frame here is voiced in every band and leaves no unvoiced overlap, so
 * its repeat is its harmonics held, and the first muted frame is that repeat
 * times (1 - n / 80) over the first 80 samples, plus the noise, and the noise
 * alone after them; at the repeat limit and at an error-rate mute.
 */
static void
test_mute_fades_the_last_frame(void) {
    char erasure[49];
    float noise[160];
    pack_b(erasure, k_erasure);
    mute_noise_reference(noise);

    for (int by_rate = 0; by_rate < 2; ++by_rate) {
        struct stream s;
        stream_settled(&s, k_voice_reach);
        if (by_rate) {
            s.prev.errorRate = 0.2f;
        } else {
            s.prev.repeatCount = MBE_MAX_FRAME_REPEATS - 1;
        }
        struct stream twin = s;
        twin.prev.errorRate = 0.0f;
        twin.prev.repeatCount = 0;
        stream_frame(&twin, erasure, 0, 0);
        assert(!has_flag(&twin, MBE_PROCESS_FLAG_MUTE));
        stream_frame(&s, erasure, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));

        double residual = 0.0, energy = 0.0;
        for (int n = 0; n < 160; ++n) {
            const double gain = (n < 80) ? 1.0 - ((double)n / 80.0) : 0.0;
            const double r = (double)s.out[n] - (double)noise[n] - (gain * (double)twin.out[n]);
            residual += r * r;
            energy += (n < 80) ? (double)twin.out[n] * (double)twin.out[n] : 0.0;
        }
        assert(energy > 80.0 * 50.0 * 50.0);
        assert(residual < 1e-4 * energy);
        for (int n = 80; n < 160; ++n) {
            assert(fabs((double)s.out[n] - (double)noise[n]) < 1e-4);
        }
        (void)residual;
        (void)energy;
    }
}

/*
 * The first muted frame plays the last synthesized frame fading out, whatever
 * that frame was and whatever mutes: a voice frame muted by the error rate, a
 * silence frame muted by it, an invalid frame muted by it after a repeat, and
 * the 4th invalid frame after three repeats. A valid voice or silence frame
 * becomes the kept model as its unheard synthesis leaves it; an invalid frame
 * keeps the repeated model. Either way the noise generator runs one frame.
 */
static void
test_mute_onsets(void) {
    char voice[49], silence[49], erasure[49];
    pack_b(voice, k_voice_next);
    unpack_bytes(silence, k_dvsi_silence_bytes);
    pack_b(erasure, k_erasure);

    for (int kind = 0; kind < 4; ++kind) {
        struct stream s;
        float want[160];
        mbe_parms faded, own;
        stream_settled(&s, k_voice_loud);
        if (kind >= 2) {
            stream_frame(&s, erasure, 0, 0);
            assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT) && !has_flag(&s, MBE_PROCESS_FLAG_MUTE));
        }
        if (kind == 3) {
            stream_frame(&s, erasure, 0, 0);
            stream_frame(&s, erasure, 0, 0);
        } else {
            s.prev.errorRate = 0.2f;
        }
        const struct stream before = s;
        const char* frame = (kind == 0) ? voice : ((kind == 1) ? silence : erasure);
        faded_mute_reference(want, &before.prev_enh, &faded);
        stream_frame(&s, frame, 0, 0);
        assert(close_to(s.out, want, 0, 160));
        assert(second_half_is_noise(s.out));
        assert(peak_abs(want) > 5.0f);
        assert(has_flag(&s, MBE_PROCESS_FLAG_SILENCE) == (kind == 1));
        assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT) == (kind >= 2));
        if (kind < 2) {
            own_model_reference(&own, &before, frame, &s);
            assert_muted_state(&s, &own, before.prev_enh.noiseSeed);
            assert(s.cur.repeatCount == 0);
        } else {
            assert_muted_state(&s, &before.prev_enh, before.prev_enh.noiseSeed);
            assert(float_array_bits_equal(s.prev_enh.noiseOverlap, faded.noiseOverlap, 96));
            assert(s.cur.repeatCount == ((kind == 3) ? MBE_MAX_FRAME_REPEATS : 2));
        }
        assert(s.prev.repeatCount == s.cur.repeatCount);
        if (kind != 0) {
            assert(history_equal(&s.prev, &before.prev));
        }
    }
}

/*
 * A mute after a silence frame fades the silence model out over the first
 * 80 samples, as after any frame: voice, a silence frame, three erasures
 * repeating it, then the 4th; and voice, a silence frame, then a voice frame
 * the error rate mutes. The faded copy carries no frame mark, so the silence
 * frame's mark (hold the model before to mid-frame) cannot carry the model at
 * full level into the second half: that half is the mute noise alone, and
 * the first is the fade against the reference plus the noise, never above the
 * repeats' peak plus the noise. The same after a V/UV row-15 frame, whose
 * late-onset mark the fade does not depend on either; the kept model keeps
 * its own mark.
 */
static void
test_mute_after_silence_frame(void) {
    char silence[49], erasure[49], voice[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);
    pack_b(erasure, k_erasure);
    pack_b(voice, k_voice_next);
    static const int k_row15[9] = {60, 15, 20, 200, 90, 17, 11, 3, 6};
    for (int k = 0; k < 4; ++k) {
        const int by_rate = k & 1;
        const int row15 = k >> 1;
        struct stream s;
        float want[160];
        mbe_parms faded;
        float repeat_peak = 0.0f;
        stream_settled(&s, k_voice_loud);
        if (row15) {
            stream_frame_b(&s, k_row15, 0, 0);
        } else {
            stream_frame(&s, silence, 0, 0);
            assert(has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
        }
        const float mark = s.prev_enh.PHIl[0];
        if (by_rate) {
            repeat_peak = peak_abs(s.out);
            s.prev.errorRate = 0.2f;
        } else {
            for (int i = 0; i < 3; ++i) {
                stream_frame(&s, erasure, 0, 0);
                assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT) && !has_flag(&s, MBE_PROCESS_FLAG_MUTE));
                repeat_peak = peak_abs(s.out);
            }
        }
        const struct stream before = s;
        faded_mute_reference(want, &before.prev_enh, &faded);
        stream_frame(&s, by_rate ? voice : erasure, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
        assert(close_to(s.out, want, 0, 160));
        assert(second_half_is_noise(s.out));
        assert(peak_abs(s.out) * 7.0f <= (repeat_peak * 7.0f) + MBE_AMBE_MUTE_NOISE_AMPLITUDE + 1e-3f);
        assert(peak_abs(s.out) * 7.0f > 2.0f * MBE_AMBE_MUTE_NOISE_AMPLITUDE);
        /* The kept model keeps its own mark. */
        if (!by_rate) {
            assert(float_bits_equal(s.prev_enh.PHIl[0], mark));
            assert(float_bits_equal(mark, row15 ? 80.0f : k_mark_previous_held));
        }
    }
}

/*
 * A mute while a tone sounds: the tone goes on from its phases and fades out
 * over the first 80 samples, added to the mute noise, which does not fade.
 * The tone is no longer held; the voice frame muted here becomes the kept
 * model, and the noise generator runs one frame.
 */
static void
test_mute_fades_a_held_tone(void) {
    struct stream s;
    char tone[49], voice[49];
    float tail[160], noise[160];
    pack_tone(tone, 7);
    pack_b(voice, k_voice_next);
    stream_settled(&s, k_voice_loud);
    stream_frame(&s, tone, 0, 0);
    assert(s.prev_enh.Vl[0] == tone_word(7));
    s.prev.errorRate = 0.2f;
    const struct stream before = s;
    tail_reference(tail, tone, &before.cur);
    mute_noise_reference(noise);
    stream_frame(&s, voice, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) && !has_flag(&s, MBE_PROCESS_FLAG_TONE));
    for (int n = 0; n < 160; ++n) {
        tail[n] += noise[n];
    }
    limit_sum(tail, 160);
    assert(close_to(s.out, tail, 0, 160));
    assert(!tones_sound() || peak_abs(tail) > 100.0f);
    mbe_parms own;
    own_model_reference(&own, &before, voice, &s);
    assert_muted_state(&s, &own, before.prev_enh.noiseSeed);
}

/* A mute as the very first frame plays the noise alone: nothing was played to
 * fade. The noise generator still runs one frame, from its start. */
static void
test_mute_first_frame(void) {
    char voice[49], erasure[49];
    float noise[160];
    pack_b(voice, k_voice_next);
    pack_b(erasure, k_erasure);
    mute_noise_reference(noise);
    for (int invalid = 0; invalid < 2; ++invalid) {
        struct stream s;
        stream_init(&s);
        /* AMBE state from the start, so that the error state below holds. */
        s.prev.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
        if (invalid) {
            s.prev.repeatCount = MBE_MAX_FRAME_REPEATS - 1;
        } else {
            s.prev.errorRate = 0.2f;
        }
        const mbe_parms enh_before = s.prev_enh;
        stream_frame(&s, invalid ? erasure : voice, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
        assert(float_array_bits_equal(s.out, noise, 160));
        assert(s.prev_enh.Vl[0] == W_AFRESH && s.cur.Vl[0] == 0);
        assert(overlap_cleared(&s.prev_enh));
        assert(float_bits_equal(s.prev_enh.noiseSeed, generator_after(-1.0f, 1)));
        if (invalid) {
            assert(model_equal(&s.prev_enh, &enh_before));
        } else {
            /* The voice frame's own model, which also became the history. */
            assert(s.prev_enh.L == s.prev.L && float_bits_equal(s.prev_enh.w0, s.prev.w0));
            assert(s.prev_enh.L != enh_before.L);
        }
    }
}

/* ---- The frame after a mute or a tone starts from silence ---- */

static double
span_difference(const float* a, const float* b, int from, int to) {
    double acc = 0.0;
    for (int n = from; n < to; ++n) {
        acc += ((double)a[n] - (double)b[n]) * ((double)a[n] - (double)b[n]);
    }
    return acc;
}

/*
 * The frame `d` released after the stream's mute, and the frame after it:
 * the first is synthesized from the kept model, starting from silence (its
 * sample 0 is 0), and consumes the request; the second continues from it as
 * any frame does. Both against renderings from the same state. `released`
 * receives the stream after the first.
 */
static void
assert_release(struct stream* s, const char d[49], struct stream* released) {
    const struct stream before = *s;
    float want[160], next_want[160];
    mbe_parms model, next_model;
    assert(start_of(&before) == MBE_SYNTHESIS_FROM_SILENCE);
    stream_frame(s, d, 0, 0);
    assert(!has_flag(s, MBE_PROCESS_FLAG_MUTE) && !has_flag(s, MBE_PROCESS_FLAG_REPEAT));
    frame_reference(want, &model, &before, d, s, MBE_SYNTHESIS_FROM_SILENCE);
    assert(close_to(s->out, want, 0, 160));
    assert(fabsf(s->out[0]) < 1e-6f);
    assert(span_energy(want, 0, 160) > 160.0);
    assert(model_equal(&s->prev_enh, &model));
    assert(float_array_bits_equal(s->prev_enh.noiseOverlap, model.noiseOverlap, 96));
    assert(s->prev_enh.Vl[0] == 0 && s->cur.Vl[0] == 0);
    assert(float_bits_equal(s->prev_enh.noiseSeed, generator_after(before.prev_enh.noiseSeed, 1)));

    *released = *s;
    stream_frame(s, d, 0, 0);
    frame_reference(next_want, &next_model, released, d, s, MBE_SYNTHESIS_CONTINUES);
    assert(close_to(s->out, next_want, 0, 160));
    assert(s->prev_enh.Vl[0] == 0);
}

/* Frames to release into: voiced in every column, all unvoiced, a silence
 * frame (which, the kept model not having been heard, holds nothing: it rises
 * from silence to its own model), and a V/UV row-15 frame (a late onset after
 * an unvoiced kept model). */
static void
release_frame(char d[49], int kind) {
    static const int k_row15[9] = {60, 15, 20, 200, 90, 17, 11, 3, 6};
    if (kind == 2) {
        unpack_bytes(d, k_dvsi_silence_bytes);
    } else {
        pack_b(d, (kind == 0) ? k_voice_reach : ((kind == 1) ? k_voice_unvoiced : k_row15));
    }
}

/*
 * After the repeat limit's mute: the kept model is the repeated one, and the
 * frame heard again starts from silence with it, whatever the frame is. A
 * second muted frame changes nothing of that.
 */
static void
test_release_after_repeat_mute(void) {
    char erasure[49], d[49];
    pack_b(erasure, k_erasure);
    for (int kind = 0; kind < 4; ++kind) {
        for (int mutes = 1; mutes <= 2; ++mutes) {
            struct stream s, released;
            stream_settled(&s, (kind == 3) ? k_voice_unvoiced : k_voice_loud);
            const mbe_parms repeated = s.prev_enh;
            for (int i = 0; i < 3 + mutes; ++i) {
                stream_frame(&s, erasure, 0, 0);
            }
            assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
            assert(model_equal(&s.prev_enh, &repeated));
            release_frame(d, kind);
            assert_release(&s, d, &released);
        }
    }
}

/*
 * After an error-rate mute the kept model is the last muted valid frame's own
 * model, not the one heard before the mute: the muted frames here differ from
 * the frames before them in pitch and level, and the frame heard again plays
 * the rendering from the muted frame's model, which differs from the one from
 * the model heard before the mute in its first half already. A silence frame
 * plays the same from either: it holds no model that was not heard.
 */
static void
test_release_after_error_rate_mute(void) {
    char muted_voice[49], d[49];
    pack_b(muted_voice, k_voice_loud);
    for (int kind = 0; kind < 4; ++kind) {
        struct stream s, released;
        stream_settled(&s, k_voice_next);
        const mbe_parms heard = s.prev_enh;
        for (int i = 0; i < 2; ++i) {
            mbe_parms own;
            s.prev.errorRate = 0.2f;
            const struct stream before = s;
            stream_frame(&s, muted_voice, 0, 0);
            own_model_reference(&own, &before, muted_voice, &s);
            assert_muted_state(&s, &own, before.prev_enh.noiseSeed);
        }
        assert(s.prev_enh.L != heard.L);
        s.prev.errorRate = 0.0f;
        release_frame(d, kind);
        struct stream other = s;
        other.prev_enh = heard;
        other.prev_enh.noiseSeed = s.prev_enh.noiseSeed;
        other.prev_enh.Vl[0] = s.prev_enh.Vl[0];
        memset(other.prev_enh.previousUw, 0, sizeof(other.prev_enh.previousUw));
        const struct stream at_release = s;
        assert_release(&s, d, &released);

        float from_kept[160], from_heard[160];
        mbe_parms scratch;
        frame_reference(from_kept, &scratch, &at_release, d, &released, MBE_SYNTHESIS_FROM_SILENCE);
        frame_reference(from_heard, &scratch, &other, d, &released, MBE_SYNTHESIS_FROM_SILENCE);
        if (kind == 2) {
            assert(close_to(from_kept, from_heard, 0, 160));
        } else {
            assert(span_difference(from_kept, from_heard, 0, 80) > 0.05 * span_energy(from_kept, 0, 80));
        }
    }
}

/*
 * A silence frame right after a mute or a tone: a silence frame holds the
 * model heard before it to mid-frame, but the kept model was not heard, so
 * from silence the frame rises to its own (quiet) model instead of playing the
 * loud speech from before the mute or the tone. After loud voiced speech and
 * a repeat-limit mute, an error-rate mute or a tone, the silence frame is the
 * from-silence rendering (plus the tone's tail after a tone) and stays within
 * a small multiple of a steady silence frame's peak.
 */
static void
test_silence_frame_after_mute_or_tone(void) {
    char silence[49], erasure[49], voice[49], tone[49];
    unpack_bytes(silence, k_dvsi_silence_bytes);
    pack_b(erasure, k_erasure);
    pack_b(voice, k_voice_loud);
    pack_tone(tone, 7);

    /* A steady silence frame after the same speech: the third in a row. */
    struct stream steady;
    stream_settled(&steady, k_voice_loud);
    const float loud_peak = peak_abs(steady.out);
    for (int i = 0; i < 3; ++i) {
        stream_frame(&steady, silence, 0, 0);
    }
    const float silence_peak = peak_abs(steady.out);
    assert(silence_peak > 0.0f && 30.0f * silence_peak < loud_peak);

    for (int kind = 0; kind < 3; ++kind) {
        struct stream s;
        float want[160], tail[160] = {0};
        mbe_parms model;
        stream_settled(&s, k_voice_loud);
        if (kind == 0) {
            for (int i = 0; i < 4; ++i) {
                stream_frame(&s, erasure, 0, 0);
            }
        } else if (kind == 1) {
            s.prev.errorRate = 0.2f;
            stream_frame(&s, voice, 0, 0);
            s.prev.errorRate = 0.0f;
        } else {
            stream_frame(&s, tone, 0, 0);
            tail_reference(tail, tone, &s.cur);
        }
        assert(has_flag(&s, (kind == 2) ? MBE_PROCESS_FLAG_TONE : MBE_PROCESS_FLAG_MUTE));
        const struct stream before = s;
        stream_frame(&s, silence, 0, 0);
        assert(s.result.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_SILENCE));
        frame_reference(want, &model, &before, silence, &s, MBE_SYNTHESIS_FROM_SILENCE);
        for (int n = 0; n < 160; ++n) {
            assert(fabsf(want[n]) <= 3.0f * silence_peak);
            want[n] += tail[n];
        }
        limit_sum(want, 160);
        assert(close_to(s.out, want, 0, 160));
        for (int n = 0; n < 160; ++n) {
            assert(fabs((double)s.out[n] - (double)tail[n]) <= 3.0 * (double)silence_peak + 1e-3);
        }
        assert(s.prev_enh.Vl[0] == 0 && s.cur.Vl[0] == 0);
    }
}

/*
 * An invalid frame that is not muted while the request stands repeats the
 * kept model from silence and consumes the request: the frame after it
 * continues from the repeat.
 */
static void
test_repeat_after_mute_starts_from_silence(void) {
    char erasure[49], voice[49];
    pack_b(erasure, k_erasure);
    pack_b(voice, k_voice_next);
    struct stream s;
    stream_settled(&s, k_voice_loud);
    s.prev.errorRate = 0.2f;
    stream_frame(&s, voice, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    s.prev.errorRate = 0.0f;
    const struct stream before = s;
    float want[160], next_want[160];
    mbe_parms model;
    repeat_reference(want, &before, MBE_SYNTHESIS_FROM_SILENCE);
    stream_frame(&s, erasure, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_REPEAT) && !has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    assert(close_to(s.out, want, 0, 160));
    assert(fabsf(s.out[0]) < 1e-6f && span_energy(want, 0, 160) > 160.0);
    assert(model_equal(&s.prev_enh, &before.prev_enh));
    assert(s.prev_enh.Vl[0] == 0 && s.cur.Vl[0] == 0);
    assert(s.cur.repeatCount == 1);
    assert(float_bits_equal(s.prev_enh.noiseSeed, generator_after(before.prev_enh.noiseSeed, 1)));
    const struct stream repeated = s;
    stream_frame(&s, voice, 0, 0);
    frame_reference(next_want, &model, &repeated, voice, &s, MBE_SYNTHESIS_CONTINUES);
    assert(close_to(s.out, next_want, 0, 160));
}

/* ---- Lost frames and tones ---- */

/*
 * Table 9 lists ID 255 as a tone of zero amplitude; DVSI's decoder plays it as
 * a lost frame, as it does any unusable index. It classifies and decodes as an
 * erasure, and processes as an erasure flagged as a tone frame, a repeat of
 * the last frame: from parameter bits without C0 context at total errors 3,
 * 4 and 5, and from the hard and soft frame entry points (soft with
 * reliabilities that differ from bit to bit). The repeat criteria still come
 * first: without C0 context the 2.1 heuristic rejects the frame at 6 errors,
 * flagged a repeat alone, and it repeats the same way.
 */
static void
test_tone_255_is_erasure(void) {
    static const unsigned erasure = MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT;
    char d[49], ambe_d[49], fr[4][24];
    mbe_soft_bit soft[4][24];
    mbe_parms cur, prev, enh;
    pack_tone(d, 255);
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_ERASURE);
    mbe_initMbeParms(&cur, &prev, &enh);
    assert(mbe_decodeAmbe2450Parms(d, &cur, &prev) == MBE_AMBE2450_FRAME_ERASURE);

    for (int total = 3; total <= 6; ++total) {
        struct stream s;
        float want[160];
        mbe_process_result result;
        stream_settled(&s, k_voice_loud);
        repeat_reference(want, &s, MBE_SYNTHESIS_CONTINUES);
        mbe_initProcessResult(&result);
        result.total_errors = total;
        const int rc = mbe_processAmbe2450Dataf(s.out, &result, d, &s.cur, &s.prev, &s.prev_enh);
        assert(rc == total);
        (void)rc;
        assert(result.flags == ((total >= 6) ? MBE_PROCESS_FLAG_REPEAT : erasure));
        assert(close_to(s.out, want, 0, 160));
        assert(s.cur.repeatCount == 1 && s.prev_enh.Vl[0] == 0);
    }

    assert(mbe_encodeAmbe3600x2450Frame(d, fr) == 0);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 24; ++j) {
            soft[i][j].bit = (uint8_t)fr[i][j];
            soft[i][j].reliability = (uint8_t)(40 + ((37 * ((i * 24) + j)) % 216));
        }
    }
    for (int soft_input = 0; soft_input < 2; ++soft_input) {
        struct stream s;
        float want[160];
        mbe_process_result result;
        stream_settled(&s, k_voice_loud);
        repeat_reference(want, &s, MBE_SYNTHESIS_CONTINUES);
        const int rc = soft_input ? mbe_processAmbe3600x2450SoftFramef(s.out, &result, (const mbe_soft_bit(*)[24])soft,
                                                                       ambe_d, &s.cur, &s.prev, &s.prev_enh)
                                  : mbe_processAmbe3600x2450Framef(s.out, &result, (const char (*)[24])fr, ambe_d,
                                                                   &s.cur, &s.prev, &s.prev_enh);
        assert(rc == 0);
        (void)rc;
        assert(memcmp(ambe_d, d, sizeof(d)) == 0);
        assert(result.flags == (MBE_PROCESS_FLAG_C0_VALID | erasure | (soft_input ? MBE_PROCESS_FLAG_SOFT_INPUT : 0u)));
        assert(close_to(s.out, want, 0, 160));
    }
}

/* A stream that plays a tone of this ID after settled speech. */
static void
stream_into_tone(struct stream* s, char tone[49], int tone_id) {
    stream_settled(s, k_voice_loud);
    pack_tone(tone, tone_id);
    stream_frame(s, tone, 0, 0);
    assert(has_flag(s, MBE_PROCESS_FLAG_TONE) && !has_flag(s, MBE_PROCESS_FLAG_REPEAT));
}

/*
 * A lost frame repeats the frame last played, a tone included: three lost
 * frames of any kind (an erasure, a tone frame with the zero-amplitude ID 255,
 * a frame the repeat criteria reject) after a tone play exactly what three
 * more tone frames play, from the running tone phases, at full level; each
 * counts toward the repeat limit and advances the error accounting, and the
 * tone, the kept model and the noise generator stay where they were. The
 * 4th fades the tone out under the mute noise, which does not fade, and runs
 * the generator one frame; the next frame starts from silence.
 */
static void
test_lost_frames_repeat_a_tone(void) {
    char tone[49], erasure[49], silent_tone[49], voice[49];
    struct stream s, twin;
    pack_b(erasure, k_erasure);
    pack_tone(silent_tone, 255);
    pack_b(voice, k_voice_next);
    stream_into_tone(&s, tone, 7);
    twin = s;
    const mbe_parms kept = s.prev_enh;
    const float generator = s.prev_enh.noiseSeed;

    for (int i = 1; i <= 3; ++i) {
        const float error_rate = s.prev.errorRate;
        if (i == 1) {
            stream_frame(&s, erasure, 0, 0);
            assert(s.result.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT));
        } else if (i == 2) {
            stream_frame(&s, silent_tone, 0, 0);
            assert(s.result.flags
                   == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT
                       | MBE_PROCESS_FLAG_TONE));
        } else {
            stream_frame(&s, tone, 4, 2);
            assert(s.result.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_REPEAT));
        }
        stream_frame(&twin, tone, 0, 0);
        assert(float_array_bits_equal(s.out, twin.out, 160));
        assert(!tones_sound() || peak_abs(s.out) > 1000.0f);
        assert(s.cur.repeatCount == i && s.prev.repeatCount == i);
        assert(fabsf(s.prev.errorRate - ((0.95f * error_rate) + (0.001064f * (float)s.cur.errorCountTotal))) < 1e-7f);
        assert(s.cur.swn == twin.cur.swn && s.cur.tonePhase == twin.cur.tonePhase);
        assert(kept_through_tone(&s.prev_enh, &kept, tone_word(7)));
        assert(float_bits_equal(s.cur.noiseSeed, generator));
        assert(s.cur.Vl[0] == 0);
        (void)error_rate;
    }

    /* The 4th: the tone's tail under the mute noise. */
    const struct stream before = s;
    float want[160], noise[160];
    tail_reference(want, tone, &before.cur);
    mute_noise_reference(noise);
    for (int n = 0; n < 160; ++n) {
        want[n] += noise[n];
    }
    limit_sum(want, 160);
    stream_frame(&s, erasure, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) && s.cur.repeatCount == MBE_MAX_FRAME_REPEATS);
    assert(close_to(s.out, want, 0, 160));
    assert_muted_state(&s, &kept, generator);

    /* The voice after it starts from silence, with no tail. */
    struct stream released;
    assert_release(&s, voice, &released);
}

/*
 * A valid tone frame after lost frames that repeated the tone goes on at once,
 * with no fade, ends the repeat run and advances the error accounting; a
 * different tone changes at once too.
 * Lost frames inside a tone are not heard: tone, two erasures, tone plays as
 * four tone frames.
 */
static void
test_tone_after_repeated_tone(void) {
    char tone[49], other[49], erasure[49];
    struct stream s, twin;
    pack_b(erasure, k_erasure);
    pack_tone(other, 32);
    stream_into_tone(&s, tone, 7);
    twin = s;
    for (int i = 0; i < 2; ++i) {
        stream_frame(&s, erasure, 0, 0);
        stream_frame(&twin, tone, 0, 0);
        assert(float_array_bits_equal(s.out, twin.out, 160));
    }
    assert(s.cur.repeatCount == 2);
    const float error_rate = s.prev.errorRate;
    stream_frame(&s, tone, 1, 1);
    stream_frame(&twin, tone, 0, 0);
    assert(float_array_bits_equal(s.out, twin.out, 160));
    assert(s.cur.repeatCount == 0 && s.prev.repeatCount == 0);
    assert(fabsf(s.prev.errorRate - ((0.95f * error_rate) + (0.001064f * 2.0f))) < 1e-7f);
    (void)error_rate;
    assert(s.prev_enh.Vl[0] == tone_word(7));

    /* Another tone after a repeated tone. */
    stream_frame(&s, erasure, 0, 0);
    mbe_parms phases = s.cur;
    float want[160];
    mbe_synthesizeTonef(want, other, &phases);
    stream_frame(&s, other, 0, 0);
    assert(float_array_bits_equal(s.out, want, 160));
    assert(s.prev_enh.Vl[0] == tone_word(32) && s.cur.repeatCount == 0);
}

/* The unvoiced noise generator, x <- 173 x + 13849 mod 2^16 from 0, one step
 * per output sample. */
static float
noise_state_after(int frames) {
    uint32_t state = 0u;
    for (int i = 0; i < frames * 160; ++i) {
        state = (173u * state + 13849u) & 0xFFFFu;
    }
    return (float)state;
}

/*
 * DVSI's noise generator runs through every frame that does not play a tone:
 * voice, silence, repeats and mutes each move it 160 samples on (a mute once,
 * whatever it renders), and a tone frame, or a lost frame that repeats a tone,
 * leaves it where it was.
 */
static void
test_noise_generator_through_frame_kinds(void) {
    static const int k_silence[9] = {124, 0, 0, 0, 0, 0, 0, 0, 0};
    struct stream s;
    char d[49];
    int frames = 0;

    stream_init(&s);
    stream_frame_b(&s, k_voice_loud, 0, 0);
    ++frames;
    assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
    stream_frame_b(&s, k_voice_next, 0, 0);
    ++frames;
    assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
    /* Three repeats, then two mutes at the repeat limit. */
    for (int i = 0; i < 5; ++i) {
        stream_frame_b(&s, k_erasure, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == (i >= 3));
        ++frames;
        assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
        assert(float_bits_equal(s.prev_enh.noiseSeed, s.cur.noiseSeed));
    }
    pack_tone(d, 7);
    for (int i = 0; i < 2; ++i) {
        stream_frame(&s, d, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_TONE));
        assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
    }
    /* ID 255 is a lost frame, which repeats the tone. */
    pack_tone(d, 255);
    stream_frame(&s, d, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_TONE) && has_flag(&s, MBE_PROCESS_FLAG_REPEAT));
    assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
    stream_frame_b(&s, k_voice_loud, 0, 0);
    ++frames;
    assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
    pack_b(d, k_silence);
    assert(mbe_classifyAmbe2450Frame(d) == MBE_AMBE2450_FRAME_SILENCE);
    stream_frame(&s, d, 0, 0);
    assert(has_flag(&s, MBE_PROCESS_FLAG_SILENCE));
    ++frames;
    assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
    assert(float_bits_equal(s.prev_enh.noiseSeed, s.cur.noiseSeed));
    /* A voice frame and a tone frame muted by the error rate: the first fades
     * the last frame out (its own model's unheard synthesis moving the state),
     * the second plays the noise alone. One frame each, as any mute. */
    pack_tone(d, 7);
    for (int i = 0; i < 2; ++i) {
        s.prev.errorRate = 0.2f;
        if (i == 0) {
            stream_frame_b(&s, k_voice_loud, 0, 0);
        } else {
            stream_frame(&s, d, 0, 0);
        }
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE));
        ++frames;
        assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
        assert(float_bits_equal(s.prev_enh.noiseSeed, s.cur.noiseSeed));
    }
    /* A tone, three lost frames repeating it (still), then the repeat
     * limit's mute fading it out (one frame). */
    s.prev.errorRate = 0.0f;
    pack_tone(d, 7);
    stream_frame(&s, d, 0, 0);
    assert(!has_flag(&s, MBE_PROCESS_FLAG_MUTE));
    for (int i = 0; i < 4; ++i) {
        stream_frame_b(&s, k_erasure, 0, 0);
        assert(has_flag(&s, MBE_PROCESS_FLAG_MUTE) == (i == 3));
        frames += (i == 3) ? 1 : 0;
        assert(float_bits_equal(s.cur.noiseSeed, noise_state_after(frames)));
        assert(float_bits_equal(s.prev_enh.noiseSeed, s.cur.noiseSeed));
    }
    /* The thread seed does not move it: every stream starts DVSI's generator from 0. */
    assert(frames == 12);
}

int
main(void) {
    test_repeat_preserves_model();
    test_frame_after_tone();
    test_release_after_repeat_mute();
    test_release_after_error_rate_mute();
    test_repeat_after_mute_starts_from_silence();
    test_silence_frame_after_mute_or_tone();
    test_dvsi_silence_vector();
    test_silence_freezes_history();
    test_silence_model_both_codes();
    test_silence_holds_previous_model();
    test_pulse_train_through_one_silence_frame();
    test_consecutive_silence_frames();
    test_silence_first_after_init();
    test_repeat_after_silence_repeats_silence();
    test_silence_flag_round_trip();
    test_erasure_run_repeats_then_mutes();
    test_repeat_advances_error_rate();
    test_error_rate_mute_boundary();
    test_mute_fades_the_last_frame();
    test_mute_onsets();
    test_mute_after_silence_frame();
    test_mute_fades_a_held_tone();
    test_mute_first_frame();
    test_noise_generator_through_frame_kinds();
    test_valid_tone_keeps_state();
    test_tone_255_is_erasure();
    test_lost_frames_repeat_a_tone();
    test_tone_after_repeated_tone();
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
