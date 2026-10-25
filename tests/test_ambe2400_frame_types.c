// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief D-STAR (AMBE 3600x2400) frame classes outside plain voice.
 *
 * DVSI's AMBE-3000 plays b0 125, a b0 126 tone frame whose index carries no
 * tone, and b0 127 frames sample for sample as it plays a lost frame. These
 * tests pin the library's matching rules: such a frame repeats the previous
 * frame exactly as a frame with too many errors does, through every decode
 * API, and counts toward the repeat limit; b0 120..124 (the null frame
 * included) stay voice; a b0 126 tone frame plays at any volume; the all-ones
 * DV frame is a b0 127 frame like any other.
 */

#include <assert.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ambe3600x2400_internal.h"
#include "mbe_tone.h"
#include "mbelib-neo/mbelib.h"

#define NVOICE 6

static void
put_b0(char d[49], int b0) {
    for (int i = 0; i < 6; ++i) {
        d[i] = (char)((b0 >> (6 - i)) & 1);
    }
    d[48] = (char)(b0 & 1);
}

static void
put_b2(char d[49], int b2) {
    static const int bits[6] = {6, 7, 8, 9, 42, 43};
    for (int i = 0; i < 6; ++i) {
        d[bits[i]] = (char)((b2 >> (5 - i)) & 1);
    }
}

static void
pack_tone(char d[49], int b0, int index, int volume) {
    memset(d, 0, 49);
    put_b0(d, b0);
    mbe_ambe2400_set_tone_index(d, index);
    mbe_tone_dstar_set_volume(d, volume);
}

/* Voiced frames with changing pitch, gain and envelope, so a repeat is
 * distinguishable from the frame it replaces. */
static void
make_voice(char voice[NVOICE][49]) {
    for (int f = 0; f < NVOICE; ++f) {
        memset(voice[f], 0, 49);
        put_b0(voice[f], 60 + (7 * f));
        put_b2(voice[f], 40 + (3 * f));
        for (int i = 18; i < 38; ++i) {
            voice[f][i] = (char)(((i * 5) + f) % 3 == 0);
        }
        voice[f][38] = 1; /* b1 >= 8: voiced low bands */
        voice[f][39] = (char)(f & 1);
        voice[f][40] = 1;
        voice[f][41] = 1;
    }
}

static void
air_frame(const char d[49], char fr[4][24]) {
    assert(mbe_encodeAmbe3600x2400Frame(d, fr) == 0);
}

static void
soft_frame(const char fr[4][24], mbe_soft_bit soft[4][24]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 24; ++j) {
            soft[i][j].bit = (uint8_t)fr[i][j];
            soft[i][j].reliability = 255;
        }
    }
}

static void
assert_same_pcm(const float* a, const float* b, size_t n, const char* what) {
    for (size_t i = 0; i < n; ++i) {
        uint32_t x, y;
        memcpy(&x, &a[i], sizeof(x));
        memcpy(&y, &b[i], sizeof(y));
        if (x != y) {
            printf("%s: sample %zu differs (%g vs %g)\n", what, i, (double)a[i], (double)b[i]);
            assert(0);
        }
    }
}

/* Two computed samples, equal up to rounding (optimized builds may fuse the
 * crossfade's multiply-adds differently in the library and here, or
 * reassociate them, which near a zero crossing of a signal in the hundreds
 * leaves a few 1e-5; 1e-4 is under a thousandth of an int16 step). */
static void
assert_close_sample(float got, float want, const char* what, int n) {
    const float tolerance = 1e-5f * fmaxf(10.0f, fabsf(want));
    if (fabsf(got - want) > tolerance) {
        printf("%s: sample %d differs (%g vs %g)\n", what, n, (double)got, (double)want);
        assert(0);
    }
}

/* A clean frame result, or one that reports `errors` protected-bit corrections. */
static void
make_result(mbe_process_result* r, int errors) {
    mbe_initProcessResult(r);
    r->flags = MBE_PROCESS_FLAG_C0_VALID;
    r->protected_errors = errors;
    r->total_errors = errors;
}

/* A frame whose C0 word is uncorrectable, with `errors` protected-bit
 * corrections: the repeat for errors. */
static void
make_uncorrectable_result(mbe_process_result* r, int errors) {
    make_result(r, errors);
    r->flags |= MBE_PROCESS_FLAG_C0_UNCORRECTABLE;
}

/* A frame result with `c0` Golay corrections in C0 and `protected` elsewhere. */
static void
make_result_c0(mbe_process_result* r, int c0, int protected_errors) {
    make_result(r, protected_errors);
    r->c0_errors = c0;
    r->total_errors = c0 + protected_errors;
}

/* Every b0 0..124 is voice; 125 is lost; a b0 126 tone frame plays exactly the
 * indexes D-STAR has a tone for, at any volume; b0 127 is lost. */
static void
test_classes(void) {
    mbe_parms cur, prev, enh;
    char d[49];
    for (int b0 = 0; b0 <= 125; ++b0) {
        memset(d, 0, sizeof(d));
        put_b0(d, b0);
        mbe_initMbeParms(&cur, &prev, &enh);
        int got = mbe_decodeAmbe2400Parms(d, &cur, &prev);
        assert(got == ((b0 == 125) ? MBE_AMBE2400_FRAME_LOST : 0));
    }
    static const int volumes[] = {0, 1, 180, 255};
    for (int index = 0; index < 256; ++index) {
        int playable = (index >= 5 && index <= 122) || (index >= 128 && index <= 147);
        for (size_t v = 0; v < sizeof(volumes) / sizeof(volumes[0]); ++v) {
            pack_tone(d, 126, index, volumes[v]);
            mbe_initMbeParms(&cur, &prev, &enh);
            assert(mbe_decodeAmbe2400Parms(d, &cur, &prev) == (playable ? index : MBE_AMBE2400_FRAME_LOST));
            pack_tone(d, 127, index, volumes[v]);
            mbe_initMbeParms(&cur, &prev, &enh);
            assert(mbe_decodeAmbe2400Parms(d, &cur, &prev) == MBE_AMBE2400_FRAME_LOST);
        }
    }
    puts("classes: b0 0..124 voice, 125 lost, 126 tone or lost by index at any volume, 127 lost");
}

/* The D-STAR null frame decodes cleanly to b0 124, a silence frame decoded as
 * voice: the AMBE+2 silence model (w0 2 pi / 32, L 14), every band unvoiced. */
static void
test_null_frame_is_voice(void) {
    static const unsigned char null_frame[9] = {0x9E, 0x8D, 0x32, 0x88, 0x26, 0x1A, 0x3F, 0x61, 0xE8};
    char fr[4][24], d[49];
    mbe_process_result r;
    mbe_parms cur, prev, enh;
    float out[160];
    assert(mbe_decodeDStarDVData(null_frame, fr) == 0);
    assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])fr, d, &r) == 0);
    mbe_initMbeParms(&cur, &prev, &enh);
    assert(mbe_processAmbe2400Dataf(out, &r, d, &cur, &prev, &enh) == 0);
    assert((r.flags & (MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT)) == 0u);
    assert(cur.L == 14 && fabsf(cur.w0 / (2.0f * (float)M_PI / 32.0f) - 1.0f) < 1e-6f);
    for (int l = 1; l <= cur.L; ++l) {
        assert(cur.Vl[l] == 0);
    }
    puts("null frame: b0 124, a silence frame decoded as voice");
}

/*
 * Decode voice[0..3], `count` copies of `special`, then voice[4..5]. With no
 * `special` the frames are instead voice[3] with an uncorrectable C0 word and
 * `errors` corrections, which takes the repeat path for errors. Writes all PCM
 * and the special frames' flags.
 */
static void
run_stream(const char voice[NVOICE][49], const char* special, int count, int errors, float* pcm, unsigned* flags) {
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    int f = 0;
    mbe_setThreadRngSeed(1);
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 4; ++i, ++f) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(pcm + ((size_t)160 * (size_t)f), &r, voice[i], &cur, &prev, &enh) == 0);
    }
    for (int i = 0; i < count; ++i, ++f) {
        if (special) {
            make_result(&r, 0);
        } else {
            make_uncorrectable_result(&r, errors);
        }
        assert(mbe_processAmbe2400Dataf(pcm + ((size_t)160 * (size_t)f), &r, special ? special : voice[3], &cur, &prev,
                                        &enh)
               >= 0);
        flags[i] = r.flags;
        assert(r.total_errors == (special ? 0 : errors));
    }
    for (int i = 4; i < NVOICE; ++i, ++f) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(pcm + ((size_t)160 * (size_t)f), &r, voice[i], &cur, &prev, &enh) == 0);
    }
}

/* Process voice[0..3], `special` through the parameter-bit (path 0), hard-frame
 * (1) or soft-frame (2) API, then voice[4..5]; writes all PCM, the prediction
 * history before the special frame and the history it commits. */
static void
run_silence_stream(const char voice[NVOICE][49], const char special[49], int path, float* pcm, mbe_parms* before,
                   mbe_parms* history) {
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    int f = 0;
    mbe_setThreadRngSeed(1);
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 4; ++i, ++f) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(pcm + ((size_t)160 * (size_t)f), &r, voice[i], &cur, &prev, &enh) == 0);
    }
    char fr[4][24], d[49];
    mbe_soft_bit soft[4][24];
    air_frame(special, fr);
    soft_frame((const char (*)[24])fr, soft);
    float* out = pcm + ((size_t)160 * (size_t)f++);
    *before = prev;
    if (path == 0) {
        assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])fr, d, &r) == 0);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, d, &cur, &prev, &enh) == 0);
    } else if (path == 1) {
        assert(mbe_processAmbe3600x2400Framef(out, &r, (const char (*)[24])fr, d, &cur, &prev, &enh) == 0);
    } else {
        assert(mbe_processAmbe3600x2400SoftFramef(out, &r, (const mbe_soft_bit(*)[24])soft, d, &cur, &prev, &enh) == 0);
    }
    assert((r.flags & (MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT)) == 0u);
    *history = prev;
    for (int i = 4; i < NVOICE; ++i, ++f) {
        make_result(&r, 0);
        r.flags |= (path == 2) ? MBE_PROCESS_FLAG_SOFT_INPUT : 0u;
        assert(mbe_processAmbe2400Dataf(pcm + ((size_t)160 * (size_t)f), &r, voice[i], &cur, &prev, &enh) == 0);
    }
}

/*
 * b0 124 is a silence frame decoded as voice, as DVSI's AMBE-3000 plays it:
 * the AMBE+2 silence model with every band unvoiced, so all 16 b1 values play
 * sample for sample alike, on every API, through the voice after it; its
 * amplitudes are decoded and committed to the prediction history as a voice
 * frame's (the history after it is the silence model, not the last voice).
 */
static void
test_b0_124_silence_frame(void) {
    char voice[NVOICE][49], silence[49];
    make_voice(voice);
    float first[160 * (NVOICE + 1)], pcm[160 * (NVOICE + 1)];
    mbe_parms before, history, first_history;
    memset(&first_history, 0, sizeof(first_history));
    const size_t n = (size_t)160 * (NVOICE + 1);
    for (int path = 0; path < 3; ++path) {
        for (int b1 = 0; b1 < 16; ++b1) {
            memcpy(silence, voice[3], sizeof(silence));
            put_b0(silence, 124);
            for (int i = 0; i < 4; ++i) {
                silence[38 + i] = (char)((b1 >> (3 - i)) & 1);
            }
            run_silence_stream((const char (*)[49])voice, silence, path, pcm, &before, &history);
            /* The committed history is the silence frame's own decode from the
             * history before it: its amplitudes and gain, not the last voice's. */
            mbe_parms decoded = before, last = before;
            assert(mbe_decodeAmbe2400Parms(silence, &decoded, &last) == 0);
            assert(history.L == 14 && decoded.L == 14);
            for (int l = 1; l <= 14; ++l) {
                assert(history.Vl[l] == 0);
            }
            assert_same_pcm(history.log2Ml + 1, decoded.log2Ml + 1, 14, "committed log2Ml");
            assert_same_pcm(history.Ml + 1, decoded.Ml + 1, 14, "committed Ml");
            assert_same_pcm(&history.gamma, &decoded.gamma, 1, "committed gamma");
            int replaced = 0;
            for (int l = 1; l <= 14; ++l) {
                uint32_t now_bits, before_bits;
                memcpy(&now_bits, &history.log2Ml[l], sizeof(now_bits));
                memcpy(&before_bits, &before.log2Ml[l], sizeof(before_bits));
                replaced |= now_bits != before_bits;
            }
            assert(replaced);
            if (path == 0 && b1 == 0) {
                memcpy(first, pcm, sizeof(first));
                first_history = history;
            }
            assert_same_pcm(pcm, first, n, "b0 124 across b1 and APIs");
        }
    }
    /* The history is the silence frame's own decoded amplitudes, which the
     * voice after it predicts from: it plays otherwise than after a lost frame,
     * which keeps the last voice frame's history. */
    char lost[49];
    memcpy(lost, voice[3], sizeof(lost));
    put_b0(lost, 125);
    unsigned lost_flags[1];
    run_stream((const char (*)[49])voice, lost, 1, 0, pcm, lost_flags);
    assert(first_history.L == 14 && first_history.w0 < 0.2f);
    int differs = 0;
    for (size_t i = (size_t)160 * 5u; i < n; ++i) {
        differs |= fabsf(pcm[i] - first[i]) > 1.0f;
    }
    assert(differs);
    puts("b0 124: a silence frame, every band unvoiced whatever b1, committed as voice");
}

/* A lost frame plays sample for sample as a frame repeated for errors and
 * keeps its own (zero) error count (which the adaptive smoothing of later
 * frames reads, so only the first repeat is compared sample for sample); every
 * kind of lost frame plays identically; the 4th in a row mutes, as on DVSI's
 * decoder. */
static void
test_lost_plays_as_repeat(void) {
    char voice[NVOICE][49];
    make_voice(voice);
    char b125[49], unused[49], old_silence[49], b127_tone[49];
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    pack_tone(unused, 126, 148, 180);
    pack_tone(old_silence, 127, 128, 0);
    pack_tone(b127_tone, 127, 130, 200);
    const char* lost[] = {b125, unused, old_silence, b127_tone};
    for (int count = 1; count <= 5; count += 2) {
        float want[160 * 11], got[160 * 11], first[160 * 11];
        unsigned want_flags[5], got_flags[5];
        const size_t exact = (size_t)160 * 5u; /* through the first lost frame */
        const size_t n = (size_t)160 * (size_t)(NVOICE + count);
        run_stream((const char (*)[49])voice, NULL, count, 4, want, want_flags);
        for (size_t k = 0; k < sizeof(lost) / sizeof(lost[0]); ++k) {
            run_stream((const char (*)[49])voice, lost[k], count, 0, got, got_flags);
            assert_same_pcm(got, want, exact, "lost frame vs error repeat");
            if (k == 0) {
                memcpy(first, got, n * sizeof(float));
            }
            assert_same_pcm(got, first, n, "lost frame kinds");
            for (int i = 0; i < count; ++i) {
                unsigned tone = (k == 0) ? 0u : MBE_PROCESS_FLAG_TONE;
                unsigned mute = want_flags[i] & MBE_PROCESS_FLAG_MUTE;
                assert(
                    got_flags[i]
                    == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT | tone | mute));
                assert(want_flags[i]
                       == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_C0_UNCORRECTABLE | MBE_PROCESS_FLAG_ERASURE
                           | MBE_PROCESS_FLAG_REPEAT | mute));
                assert((i < 3) ? (mute == 0u) : (i == 3) ? (mute != 0u) : 1);
            }
        }
    }
    puts("lost frames: b0 125, unused tone index, b0 127 play as error repeats; the 4th mutes");
}

/* Lost frames and frames repeated for errors share one repeat count. */
static void
test_mixed_repeat_count(void) {
    char voice[NVOICE][49], b125[49];
    make_voice(voice);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    float out[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 4; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
    }
    for (int i = 0; i < 4; ++i) {
        int errored = (i & 1) != 0;
        if (errored) {
            make_uncorrectable_result(&r, 4);
        } else {
            make_result(&r, 0);
        }
        assert(mbe_processAmbe2400Dataf(out, &r, errored ? voice[3] : b125, &cur, &prev, &enh) >= 0);
        assert((r.flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
        assert(((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u) == (i == 3));
    }
    puts("repeat count: lost and errored frames count together");
}

/* Hard and soft air frames reach the same lost-frame result as parameter bits. */
static void
test_frame_apis_agree(void) {
    char voice[NVOICE][49], b125[49], unused[49];
    make_voice(voice);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    pack_tone(unused, 126, 200, 180);
    /* Five lost frames in a row pass the repeat limit, so the mute is covered too. */
    const char* seq[] = {voice[0], voice[1], voice[2], b125, unused, b125, b125, unused, voice[4], voice[5]};
    const size_t n = sizeof(seq) / sizeof(seq[0]);
    mbe_parms c[3], p[3], e[3];
    for (int k = 0; k < 3; ++k) {
        mbe_initMbeParms(&c[k], &p[k], &e[k]);
    }
    for (size_t i = 0; i < n; ++i) {
        char fr[4][24], d[49], dh[49], ds[49];
        mbe_soft_bit soft[4][24];
        mbe_process_result r[3];
        float out[3][160];
        air_frame(seq[i], fr);
        soft_frame((const char (*)[24])fr, soft);
        /* The frame encoder rewrites bit 24 (a parity bit), so decode the bits back. */
        assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])fr, d, &r[0]) == 0);
        /* Mute noise comes from a per-thread RNG: give each path the same draws. */
        make_result(&r[0], 0);
        mbe_setThreadRngSeed((uint32_t)i + 1u);
        assert(mbe_processAmbe2400Dataf(out[0], &r[0], d, &c[0], &p[0], &e[0]) >= 0);
        mbe_setThreadRngSeed((uint32_t)i + 1u);
        assert(mbe_processAmbe3600x2400Framef(out[1], &r[1], (const char (*)[24])fr, dh, &c[1], &p[1], &e[1]) >= 0);
        mbe_setThreadRngSeed((uint32_t)i + 1u);
        assert(
            mbe_processAmbe3600x2400SoftFramef(out[2], &r[2], (const mbe_soft_bit(*)[24])soft, ds, &c[2], &p[2], &e[2])
            >= 0);
        assert(memcmp(dh, d, sizeof(d)) == 0 && memcmp(ds, d, sizeof(d)) == 0);
        assert(r[1].flags == r[0].flags);
        assert(r[2].flags == (r[0].flags | MBE_PROCESS_FLAG_SOFT_INPUT));
        assert_same_pcm(out[1], out[0], 160, "hard frame vs parameter bits");
        assert_same_pcm(out[2], out[0], 160, "soft frame vs parameter bits");
    }
    puts("frame APIs: hard, soft and parameter-bit decoding agree on lost frames");
}

/* A valid b0 126 tone frame at volume 0 is a tone frame, not a reset or a repeat. */
static void
test_zero_volume_tone(void) {
    char voice[NVOICE][49], tone[49];
    make_voice(voice);
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    float out[160];
    for (int index = 5; index <= 147; index += (index == 122) ? 6 : 1) {
        pack_tone(tone, 126, index, 0);
        mbe_initMbeParms(&cur, &prev, &enh);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[0], &cur, &prev, &enh) == 0);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
        assert(r.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE));
        assert(prev.repeatCount == 0);
    }
    puts("zero-volume tones: every playable index is a tone frame");
}

/* A played tone ends a run of repeats: after three lost or errored frames, a
 * tone, then another repeat, the repeat count restarts instead of muting. */
static void
test_tone_ends_repeat_run(void) {
    char voice[NVOICE][49], b125[49];
    make_voice(voice);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    static const int tones[][2] = {{128, 0}, {128, 200}, {30, 204}};
    for (size_t t = 0; t < sizeof(tones) / sizeof(tones[0]); ++t) {
        for (int errored = 0; errored <= 1; ++errored) {
            char tone[49];
            mbe_parms cur, prev, enh;
            mbe_process_result r;
            float out[160];
            pack_tone(tone, 126, tones[t][0], tones[t][1]);
            mbe_initMbeParms(&cur, &prev, &enh);
            make_result(&r, 0);
            assert(mbe_processAmbe2400Dataf(out, &r, voice[0], &cur, &prev, &enh) == 0);
            for (int i = 0; i < 3; ++i) {
                if (errored) {
                    make_uncorrectable_result(&r, 4);
                } else {
                    make_result(&r, 0);
                }
                assert(mbe_processAmbe2400Dataf(out, &r, errored ? voice[1] : b125, &cur, &prev, &enh) >= 0);
                assert((r.flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
            }
            make_result(&r, 0);
            assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
            assert(r.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE));
            assert(prev.repeatCount == 0);
            if (errored) {
                make_uncorrectable_result(&r, 4);
            } else {
                make_result(&r, 0);
            }
            assert(mbe_processAmbe2400Dataf(out, &r, errored ? voice[1] : b125, &cur, &prev, &enh) >= 0);
            assert((r.flags & (MBE_PROCESS_FLAG_REPEAT | MBE_PROCESS_FLAG_MUTE)) == MBE_PROCESS_FLAG_REPEAT);
            assert(prev.repeatCount == 1);
        }
    }
    puts("repeat count: a played tone ends a run of repeats");
}

/*
 * The all-ones DV frame (FF x 9, b0 127 after two corrections) plays as a lost
 * frame, as on DVSI's decoder for runs of 1, 3 and 20: through the hard-frame
 * API its first frame repeats the last frame exactly as a frame repeated for
 * errors does, the 4th mutes, the frames after it play mute noise, and the
 * voice after the run keeps the prediction history.
 */
static void
test_all_ones_frame(void) {
    unsigned char ones[9];
    char fr[4][24], d[49];
    char voice[NVOICE][49];
    mbe_process_result r;
    mbe_parms cur, prev, enh, rcur, rprev, renh;
    float out[160], want[160 * 7];
    unsigned want_flags[1];
    short pcm[160];
    memset(ones, 0xFF, sizeof(ones));
    make_voice(voice);
    assert(mbe_decodeDStarDVData(ones, fr) == 0);
    assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])fr, d, &r) == 2);
    for (int i = 0; i < 6; ++i) {
        assert(d[i] == 1);
    }
    assert(d[48] == 1); /* b0 127 */
    mbe_initMbeParms(&cur, &prev, &enh);
    assert(mbe_decodeAmbe2400Parms(d, &cur, &prev) == MBE_AMBE2400_FRAME_LOST);
    run_stream((const char (*)[49])voice, NULL, 1, 4, want, want_flags);
    static const int runs[] = {1, 3, 20};
    for (size_t k = 0; k < sizeof(runs) / sizeof(runs[0]); ++k) {
        mbe_setThreadRngSeed(1);
        mbe_initMbeParms(&cur, &prev, &enh);
        for (int i = 0; i < 4; ++i) {
            make_result(&r, 0);
            assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
        }
        for (int i = 0; i < runs[k]; ++i) {
            assert(mbe_processAmbe3600x2400Framef(out, &r, (const char (*)[24])fr, d, &cur, &prev, &enh) == 2);
            const unsigned mute = (i >= 3) ? MBE_PROCESS_FLAG_MUTE : 0u;
            assert(r.flags
                   == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT
                       | MBE_PROCESS_FLAG_TONE | mute));
            if (i == 0) {
                assert_same_pcm(out, want + (size_t)160 * 4u, 160, "all-ones frame vs error repeat");
            }
            if (i >= 4) {
                mbe_floattoshort(out, pcm);
                for (int n = 0; n < 160; ++n) {
                    assert(abs(pcm[n]) <= 1);
                }
            }
        }
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[4], &cur, &prev, &enh) == 0);
        /* The same voice with no all-ones frames in between. */
        mbe_initMbeParms(&rcur, &rprev, &renh);
        for (int i = 0; i < 5; ++i) {
            make_result(&r, 0);
            assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &rcur, &rprev, &renh) == 0);
        }
        assert_same_pcm(prev.log2Ml, rprev.log2Ml, sizeof(prev.log2Ml) / sizeof(prev.log2Ml[0]),
                        "log2Ml after the all-ones run");
        assert_same_pcm(&prev.gamma, &rprev.gamma, 1, "gamma after the all-ones run");
    }
    puts("all-ones frame: plays as a lost frame for runs of 1, 3 and 20, prediction kept");
}

/* While the error rate mutes hard-decision voice, a tone frame is muted too
 * (DVSI's decoder stays silent on a channel error that reads as a tone). */
static void
test_tone_muted_by_error_rate(void) {
    char voice[NVOICE][49], tone[49];
    make_voice(voice);
    pack_tone(tone, 126, 76, 212);
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    float out[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 40; ++i) {
        make_result(&r, 3);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i % NVOICE], &cur, &prev, &enh) == 3);
    }
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
    assert(r.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_MUTE));
    float peak = 0.0f;
    for (int i = 0; i < 160; ++i) {
        peak = (out[i] > peak) ? out[i] : ((-out[i] > peak) ? -out[i] : peak);
    }
    assert(peak <= 1.5f / 7.0f); /* mute noise (float is int16 / 7); a tone at volume 212 peaks above 100 */
    /* Soft decisions do not mute on the error rate, so the tone plays. */
    mbe_soft_bit soft[4][24];
    char fr[4][24], d[49];
    air_frame(tone, fr);
    soft_frame((const char (*)[24])fr, soft);
    mbe_processAmbe3600x2400SoftFramef(out, &r, (const mbe_soft_bit(*)[24])soft, d, &cur, &prev, &enh);
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) == 0u);
    /* A muted tone still advances the tone state, as muted voice advances its
     * synthesis: the tone after it continues exactly as after a played tone. */
    mbe_parms mc = cur, mp = prev, me = enh, pc = cur, pp = prev, pe = enh;
    pc.errorRate = 0.0f;
    pp.errorRate = 0.0f;
    float muted[160], played[160], after_muted[160], after_played[160];
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(muted, &r, tone, &mc, &mp, &me) == 0);
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(played, &r, tone, &pc, &pp, &pe) == 0);
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) == 0u);
    assert(mbe_processAmbe3600x2400SoftFramef(after_muted, &r, (const mbe_soft_bit(*)[24])soft, d, &mc, &mp, &me) == 0);
    assert(mbe_processAmbe3600x2400SoftFramef(after_played, &r, (const mbe_soft_bit(*)[24])soft, d, &pc, &pp, &pe)
           == 0);
    assert_same_pcm(after_muted, after_played, 160, "tone after a muted tone vs after a played tone");
    puts("tones: muted with the error rate on hard decisions, tone state still advancing");
}

/*
 * Tone transitions, as on DVSI's decoder: the first tone frame crossfades
 * linearly over its first 80 samples from the last frame repeated (what a lost
 * frame plays there) to the tone, and the first voice frame after it from the
 * tone, which goes on, to the voice (what it plays after that lost frame).
 * Outside the crossfades the output is the plain tone and the plain voice.
 */
static void
test_tone_transitions(void) {
    char voice[NVOICE][49], tone[49], lost[49];
    make_voice(voice);
    pack_tone(tone, 126, 130, 180);
    memcpy(lost, voice[3], sizeof(lost));
    put_b0(lost, 125);
    mbe_parms cur, prev, enh, lcur, lprev, lenh;
    mbe_process_result r;
    float out[160], ref[160], plain[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    mbe_initMbeParms(&lcur, &lprev, &lenh);
    for (int i = 0; i < 4; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &lcur, &lprev, &lenh) == 0);
    }
    /* The tone frame: the lost frame's speech fading out, the tone fading in. */
    mbe_parms tone_state = cur;
    mbe_synthesizeTonefdstar(plain, tone, &tone_state, 130);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
    assert(r.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE));
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(ref, &r, lost, &lcur, &lprev, &lenh) == 0);
    for (int n = 0; n < 160; ++n) {
        const float gain = (n >= 80) ? 1.0f : (float)n / 80.0f;
        const float want = (gain * plain[n]) + ((1.0f - gain) * ref[n]);
        assert_close_sample(out[n], want, "tone frame crossfade", n);
    }
    /* The voice after it: the tone going on and fading out, the voice fading in. */
    mbe_synthesizeTonefdstar(plain, tone, &tone_state, 130);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, voice[4], &cur, &prev, &enh) == 0);
    assert((r.flags & (MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_MUTE)) == 0u);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(ref, &r, voice[4], &lcur, &lprev, &lenh) == 0);
    for (int n = 0; n < 160; ++n) {
        const float gain = (n >= 80) ? 1.0f : (float)n / 80.0f;
        const float want = (gain * ref[n]) + ((1.0f - gain) * plain[n]);
        assert_close_sample(out[n], want, "voice after the tone crossfade", n);
    }
    /* And the voice after that plays as after the lost frame. */
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, voice[5], &cur, &prev, &enh) == 0);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(ref, &r, voice[5], &lcur, &lprev, &lenh) == 0);
    assert_same_pcm(out, ref, 160, "second voice frame after the tone");
    /* After the repeat limit's mute the last frame is silent: the tone fades
     * in with no speech under it. */
    for (int i = 0; i < 5; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, lost, &cur, &prev, &enh) == 0);
    }
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
    tone_state = cur;
    mbe_synthesizeTonefdstar(plain, tone, &tone_state, 130);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
    for (int n = 0; n < 160; ++n) {
        const float gain = (n >= 80) ? 1.0f : (float)n / 80.0f;
        const float want = (gain * plain[n]) + ((1.0f - gain) * 0.0f);
        assert_close_sample(out[n], want, "tone after the mute", n);
    }
    puts("tones: crossfade in and out over the first 80 samples, as on DVSI's decoder");
}

/*
 * Muting is decided per frame from that frame's decisions (hard decisions mute
 * on the error rate, soft ones do not), and a transition looks at whether the
 * last frame was muted, not at how this frame is decoded: hard-decision voice
 * muted by the error rate leaves no speech under a soft-decision tone, and a
 * muted hard-decision tone leaves no tail under soft-decision voice.
 */
static void
test_tone_transitions_after_mute(void) {
    char voice[NVOICE][49], tone[49];
    make_voice(voice);
    pack_tone(tone, 126, 32, 180);
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    float out[160], plain[160], ref[160];
    /* (a) muted hard-decision voice, then a soft-decision tone. */
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 40; ++i) {
        make_result(&r, 3);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i % NVOICE], &cur, &prev, &enh) == 3);
    }
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
    mbe_parms tone_state = cur;
    mbe_synthesizeTonefdstar(plain, tone, &tone_state, 32);
    make_result(&r, 0);
    r.flags |= MBE_PROCESS_FLAG_SOFT_INPUT;
    assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) == 0u);
    for (int n = 0; n < 160; ++n) {
        const float gain = (n >= 80) ? 1.0f : (float)n / 80.0f;
        assert_close_sample(out[n], gain * plain[n], "soft tone after muted hard voice", n);
    }
    /* (b) a muted hard-decision tone, then soft-decision voice. */
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 40; ++i) {
        make_result(&r, 3);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i % NVOICE], &cur, &prev, &enh) == 3);
    }
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
    mbe_parms rc = cur, rp = prev, re = enh;
    rp.Vl[0] = 0; /* the same state with no tone held */
    make_result(&r, 0);
    r.flags |= MBE_PROCESS_FLAG_SOFT_INPUT;
    assert(mbe_processAmbe2400Dataf(out, &r, voice[0], &cur, &prev, &enh) == 0);
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) == 0u);
    make_result(&r, 0);
    r.flags |= MBE_PROCESS_FLAG_SOFT_INPUT;
    assert(mbe_processAmbe2400Dataf(ref, &r, voice[0], &rc, &rp, &re) == 0);
    for (int n = 0; n < 160; ++n) {
        const float gain = (n >= 80) ? 1.0f : (float)n / 80.0f;
        assert_close_sample(out[n], gain * ref[n], "soft voice after a muted hard tone", n);
    }
    puts("tones: a muted last frame leaves no speech under a tone and no tail under voice");
}

/* int16 output of one frame. */
static void
to_int16(const float* pcm, short* out) {
    mbe_floattoshort(pcm, out);
}

/*
 * The repeat limit's mute, as on DVSI's decoder: the 4th lost frame in a row
 * fades the last frame out instead of cutting to noise, the frames after it
 * are DVSI-level mute noise (int16 0 or +-1), and the prediction history is
 * kept, so the voice after the mute decodes the same model as with no lost
 * frames at all.
 */
static void
test_repeat_limit_mute(void) {
    char voice[NVOICE][49], b125[49];
    make_voice(voice);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    mbe_parms cur, prev, enh, rcur, rprev, renh;
    mbe_process_result r;
    float out[160];
    short pcm[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    mbe_initMbeParms(&rcur, &rprev, &renh);
    for (int i = 0; i < 4; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &rcur, &rprev, &renh) == 0);
    }
    for (int i = 0; i < 5; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, b125, &cur, &prev, &enh) == 0);
        assert(((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u) == (i >= 3));
        to_int16(out, pcm);
        double head = 0.0, tail = 0.0;
        int peak = 0;
        for (int n = 0; n < 160; ++n) {
            head += (n < 40) ? (double)pcm[n] * pcm[n] : 0.0;
            tail += (n >= 120) ? (double)pcm[n] * pcm[n] : 0.0;
            peak = (abs(pcm[n]) > peak) ? abs(pcm[n]) : peak;
        }
        if (i == 3) {
            /* Fading out: the previous frame continues at the start, gone by the end. */
            assert(head > 100.0 * 40.0 * 40.0 && tail < 40.0 * 4.0);
        } else if (i == 4) {
            assert(peak <= 1);
        }
    }
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, voice[4], &cur, &prev, &enh) == 0);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, voice[4], &rcur, &rprev, &renh) == 0);
    assert_same_pcm(prev.log2Ml, rprev.log2Ml, sizeof(prev.log2Ml) / sizeof(prev.log2Ml[0]), "log2Ml after the mute");
    assert_same_pcm(&prev.gamma, &rprev.gamma, 1, "gamma after the mute");
    puts("repeat limit: fades out, mutes at DVSI's level, keeps the prediction history");
}

/* Reaching the repeat limit while the error rate mutes does not fade the
 * suppressed speech back in: the output stays DVSI-level noise. */
static void
test_repeat_limit_under_error_mute(void) {
    char voice[NVOICE][49];
    make_voice(voice);
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    float out[160];
    short pcm[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 40; ++i) {
        make_result(&r, 3);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i % NVOICE], &cur, &prev, &enh) == 3);
    }
    assert((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
    for (int i = 0; i < 6; ++i) {
        make_uncorrectable_result(&r, 4);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i % NVOICE], &cur, &prev, &enh) == 4);
        assert((r.flags & (MBE_PROCESS_FLAG_REPEAT | MBE_PROCESS_FLAG_MUTE))
               == (MBE_PROCESS_FLAG_REPEAT | MBE_PROCESS_FLAG_MUTE));
        to_int16(out, pcm);
        for (int n = 0; n < 160; ++n) {
            assert(abs(pcm[n]) <= 1);
        }
    }
    puts("repeat limit under the error-rate mute: stays at the mute level");
}

/* Every frame advances the unvoiced noise generator once, muted frames
 * included, so after voice, 8 lost frames and voice the generator is where
 * 13 voice frames leave it (it used to rewind after the mute). */
static void
test_repeat_limit_noise_continuity(void) {
    char voice[NVOICE][49], b125[49];
    make_voice(voice);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    mbe_parms cur, prev, enh, rcur, rprev, renh;
    mbe_process_result r;
    float out[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    mbe_initMbeParms(&rcur, &rprev, &renh);
    for (int i = 0; i < 13; ++i) {
        const char* d = (i < 4 || i == 12) ? voice[i % NVOICE] : b125;
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, d, &cur, &prev, &enh) == 0);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i % NVOICE], &rcur, &rprev, &renh) == 0);
    }
    assert_same_pcm(&cur.noiseSeed, &rcur.noiseSeed, 1, "noise seed after the mute");
    puts("repeat limit: the noise generator advances on every frame");
}

/* A clean frame at the largest gain and b3 (b2 63, b3 511) held for 50
 * frames plays every frame, as DVSI's AMBE-3000 does (#158). Its prediction
 * builds the amplitudes up to about 2^45; the enhancement's weights used to
 * overflow on the way, and from the 5th copy the output was NaN, played as
 * silence. */
static void
test_loud_held_frame_plays(void) {
    static const unsigned char dv[9] = {0x9e, 0xb0, 0x6f, 0xee, 0xf8, 0x7b, 0x9f, 0xfe, 0x74};
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    char fr[4][24], d[49];
    float out[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int f = 0; f < 50; ++f) {
        assert(mbe_decodeDStarDVData(dv, fr) == 0);
        assert(mbe_processAmbe3600x2400Framef(out, &r, (const char (*)[24])fr, d, &cur, &prev, &enh) == 0);
        assert((r.flags & (MBE_PROCESS_FLAG_REPEAT | MBE_PROCESS_FLAG_MUTE)) == 0u);
        double energy = 0.0;
        for (int n = 0; n < 160; ++n) {
            uint32_t bits;
            memcpy(&bits, &out[n], sizeof(bits));
            assert((bits & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000));
            energy += (double)out[n] * out[n];
        }
        /* At least -30 dBFS (the float output is the int16 scale over 7). */
        assert(energy / 160.0 > pow(32768.0 * 0.0316 / 7.0, 2.0));
    }
    puts("loud held frame: every frame plays");
}

/* One step of the next k-combination of 0..n-1 in idx[]; 0 after the last. */
static int
next_combination(int* idx, int k, int n) {
    int i = k - 1;
    while (i >= 0 && idx[i] == n - k + i) {
        --i;
    }
    if (i < 0) {
        return 0;
    }
    ++idx[i];
    for (int j = i + 1; j < k; ++j) {
        idx[j] = idx[j - 1] + 1;
    }
    return 1;
}

/*
 * The frame decoders flag every C0 error pattern of weight 4, which the
 * extended Golay (24,12) code detects, as uncorrectable, and none of weight 1
 * to 3, the soft decoder from its hard decisions; the error counts stay the
 * corrected data bits (DVSI's AMBE-3000 plays exactly the weight-4 frames as
 * lost, #155).
 */
static void
test_c0_uncorrectable_detection(void) {
    char voice[NVOICE][49], fr[4][24], clean_d[49];
    make_voice(voice);
    air_frame(voice[2], fr);
    mbe_process_result r;
    assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])fr, clean_d, &r) == 0);
    int flagged[5] = {0, 0, 0, 0, 0}, patterns[5] = {0, 0, 0, 0, 0};
    for (int w = 1; w <= 4; ++w) {
        int idx[4] = {0, 1, 2, 3};
        do {
            char bad[4][24], legacy[4][24], d[49];
            memcpy(bad, fr, sizeof(bad));
            for (int i = 0; i < w; ++i) {
                bad[0][idx[i]] ^= 1;
            }
            memcpy(legacy, bad, sizeof(legacy));
            const int c0 = mbe_eccAmbe3600x2400C0(legacy);
            assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])bad, d, &r) >= 0);
            const int is_flagged = (r.flags & MBE_PROCESS_FLAG_C0_UNCORRECTABLE) != 0u;
            assert(is_flagged == (w == 4));
            assert(r.c0_errors == c0);
            if (w < 4) {
                assert(memcmp(d, clean_d, sizeof(d)) == 0);
            }
            mbe_soft_bit soft[4][24];
            soft_frame((const char (*)[24])bad, soft);
            assert(mbe_decodeAmbe3600x2400SoftFrame((const mbe_soft_bit(*)[24])soft, d, &r) >= 0);
            assert(((r.flags & MBE_PROCESS_FLAG_C0_UNCORRECTABLE) != 0u) == is_flagged);
            flagged[w] += is_flagged;
            ++patterns[w];
        } while (next_combination(idx, w, 24));
    }
    assert(patterns[1] == 24 && patterns[2] == 276 && patterns[3] == 2024 && patterns[4] == 10626);
    assert(flagged[4] == 10626);
    puts("C0: every weight-4 error pattern, and none lighter, is uncorrectable");
}

/* Air frames of voice[0..5] with frame 3 replaced by `third`. */
static void
air_stream(const char voice[NVOICE][49], const char third[4][24], char frames[NVOICE][4][24]) {
    for (int f = 0; f < NVOICE; ++f) {
        if (f == 3) {
            memcpy(frames[f], third, sizeof(frames[f]));
        } else {
            air_frame(voice[f], frames[f]);
        }
    }
}

/*
 * A frame with an uncorrectable C0 word plays as a lost frame (sample for
 * sample as a b0 125 frame there), and the frame after it predicts from the
 * frame's own decode, as on DVSI's AMBE-3000, where the frames after it differ
 * from LOST_FRAME's (#155). The staged decode passes the flag on to the
 * parameter-bit API, and fully reliable soft decisions play as hard ones.
 */
static void
test_c0_uncorrectable_plays_lost(void) {
    char voice[NVOICE][49], b125[49], bad[4][24], lost[4][24];
    make_voice(voice);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    air_frame(voice[3], bad);
    for (int i = 0; i < 4; ++i) {
        bad[0][12 + (3 * i)] ^= 1; /* four C0 data bits */
    }
    air_frame(b125, lost);
    char frames[2][NVOICE][4][24];
    air_stream((const char (*)[49])voice, (const char (*)[24])bad, frames[0]);
    air_stream((const char (*)[49])voice, (const char (*)[24])lost, frames[1]);
    float pcm[3][NVOICE][160];
    mbe_parms state[3][3], before, after;
    mbe_process_result r;
    char d[49];
    memset(&before, 0, sizeof(before));
    memset(&after, 0, sizeof(after));
    for (int run = 0; run < 3; ++run) {
        mbe_parms* s3 = state[run];
        mbe_setThreadRngSeed(1);
        mbe_initMbeParms(&s3[0], &s3[1], &s3[2]);
        for (int f = 0; f < NVOICE; ++f) {
            const char (*fr)[24] = (const char (*)[24])frames[run == 1][f];
            if (run == 0 && f == 3) {
                before = s3[1];
            }
            if (run == 2) {
                /* Staged: decode, then process the bits with the decode's result. */
                assert(mbe_decodeAmbe3600x2400Frame(fr, d, &r) >= 0);
                assert(mbe_processAmbe2400Dataf(pcm[run][f], &r, d, &s3[0], &s3[1], &s3[2]) >= 0);
            } else {
                assert(mbe_processAmbe3600x2400Framef(pcm[run][f], &r, fr, d, &s3[0], &s3[1], &s3[2]) >= 0);
            }
            if (run != 1 && f == 3) {
                assert(r.flags
                       == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_C0_UNCORRECTABLE | MBE_PROCESS_FLAG_ERASURE
                           | MBE_PROCESS_FLAG_REPEAT));
                if (run == 0) {
                    after = s3[1];
                }
            }
        }
    }
    assert_same_pcm(pcm[0][0], pcm[2][0], (size_t)160 * NVOICE, "staged decode vs frame API");
    assert_same_pcm(pcm[0][0], pcm[1][0], (size_t)160 * 4, "uncorrectable C0 vs lost frame");
    /* The history it commits is its own decode from the history before it. */
    mbe_parms decoded = before, history = before;
    assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])bad, d, &r) >= 0);
    assert(mbe_decodeAmbe2400Parms(d, &decoded, &history) == 0);
    assert_same_pcm(&after.log2Ml[1], &decoded.log2Ml[1], (size_t)decoded.L, "history after the frame");
    assert(after.L == before.L); /* the repeated model */
    /* The next frame predicts from it, so it is not the frame after a lost frame. */
    mbe_parms next = decoded, next_prev = decoded;
    air_frame(voice[4], bad);
    assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])bad, d, &r) == 0);
    assert(mbe_decodeAmbe2400Parms(d, &next, &next_prev) == 0);
    int differs = 0;
    for (int n = 0; n < 160; ++n) {
        uint32_t x, y;
        memcpy(&x, &pcm[0][4][n], sizeof(x));
        memcpy(&y, &pcm[1][4][n], sizeof(y));
        differs |= x != y;
    }
    assert(differs);
    mbe_parms replay = after, replay_cur = before;
    assert(mbe_decodeAmbe2400Parms(d, &replay_cur, &replay) == 0);
    assert_same_pcm(&replay_cur.log2Ml[1], &next.log2Ml[1], (size_t)next.L, "the next frame's prediction");
    mbe_soft_bit soft[4][24];
    mbe_parms c, p, e;
    float out[NVOICE][160];
    mbe_setThreadRngSeed(1);
    mbe_initMbeParms(&c, &p, &e);
    for (int f = 0; f < NVOICE; ++f) {
        soft_frame((const char (*)[24])frames[0][f], soft);
        assert(mbe_processAmbe3600x2400SoftFramef(out[f], &r, (const mbe_soft_bit(*)[24])soft, d, &c, &p, &e) >= 0);
    }
    assert_same_pcm(out[0], pcm[0][0], (size_t)160 * NVOICE, "soft vs hard frames");
    puts("uncorrectable C0: plays as a lost frame, the next frame predicts from its decode");
}

/* The history an uncorrectable frame leaves outlasts the frames that do not
 * replace it: a lost frame or a tone between it and the next voice frame. */
static void
test_c0_uncorrectable_history_kept(void) {
    char voice[NVOICE][49], b125[49], tone[49], d[49];
    make_voice(voice);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    pack_tone(tone, 126, 32, 180);
    const char* between[] = {b125, tone};
    for (size_t k = 0; k < sizeof(between) / sizeof(between[0]); ++k) {
        mbe_parms cur, prev, enh, decoded, history;
        mbe_process_result r;
        float out[160];
        mbe_initMbeParms(&cur, &prev, &enh);
        for (int i = 0; i < 3; ++i) {
            make_result(&r, 0);
            assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
        }
        decoded = cur;
        history = prev;
        assert(mbe_decodeAmbe2400Parms(voice[3], &decoded, &history) == 0);
        make_uncorrectable_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[3], &cur, &prev, &enh) == 0);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, between[k], &cur, &prev, &enh) == 0);
        mbe_parms want = decoded, want_prev = decoded;
        memcpy(d, voice[4], sizeof(d));
        assert(mbe_decodeAmbe2400Parms(d, &want, &want_prev) == 0);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[4], &cur, &prev, &enh) == 0);
        assert(cur.L == want.L);
        assert_same_pcm(&prev.log2Ml[1], &want.log2Ml[1], (size_t)want.L, "prediction after the frame between");
        assert((prev.Vl[0] & ~0xFFFF) == 0); /* a decoded voice frame drops the held history */
    }
    puts("uncorrectable C0: its history outlasts a lost frame or a tone after it");
}

/* State that holds an uncorrectable frame's history works with the other
 * entry points: a direct parameter decode predicts from it and leaves its own
 * model as the history, and the encoder predicts from it as the decoder does,
 * so decoding its bits from that state gives the model it reports. */
static void
test_c0_uncorrectable_history_entry_points(void) {
    char voice[NVOICE][49];
    make_voice(voice);
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    float out[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 3; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
    }
    make_uncorrectable_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, voice[3], &cur, &prev, &enh) == 0);
    assert((prev.Vl[0] & ~0x3FFFF) != 0);
    /* Direct decoding, as a caller of the staged API would chain it. */
    mbe_parms dcur = cur, dprev = prev, pcur = cur, pprev = prev, penh = enh;
    assert(mbe_decodeAmbe2400Parms(voice[4], &dcur, &dprev) == 0);
    assert((dcur.Vl[0] & ~0x3FFFF) == 0);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(out, &r, voice[4], &pcur, &pprev, &penh) == 0);
    assert(dcur.L == pprev.L);
    assert_same_pcm(&dcur.log2Ml[1], &pprev.log2Ml[1], (size_t)dcur.L, "direct decode after the frame");
    mbe_moveMbeParms(&dcur, &dprev);
    mbe_moveMbeParms(&pprev, &pcur);
    assert(mbe_decodeAmbe2400Parms(voice[5], &dcur, &dprev) == 0);
    assert(mbe_decodeAmbe2400Parms(voice[5], &pcur, &pprev) == 0);
    assert_same_pcm(&dcur.log2Ml[1], &pcur.log2Ml[1], (size_t)dcur.L, "direct decode, the frame after");
    /* The encoder, given the same state. */
    mbe_ambe2400_encoder* enc = mbe_ambe2400EncoderAlloc();
    assert(enc != NULL);
    float pcm[160];
    int voiced = 0;
    for (int f = 0; f < 12; ++f) {
        for (int n = 0; n < 160; ++n) {
            const double t = (double)((f * 160) + n) / 8000.0;
            pcm[n] = (float)(6000.0 * (sin(2.0 * M_PI * 180.0 * t) + (0.5 * sin(2.0 * M_PI * 540.0 * t))));
        }
        char d[49];
        mbe_parms ecur = cur, eprev = prev;
        assert(mbe_encodeAmbe2400Parms(enc, pcm, d, &ecur, &eprev) == 0);
        mbe_parms check = cur, check_prev = prev;
        if (mbe_decodeAmbe2400Parms(d, &check, &check_prev) != 0) {
            continue;
        }
        ++voiced;
        assert(check.L == ecur.L);
        assert_same_pcm(&check.log2Ml[1], &ecur.log2Ml[1], (size_t)check.L, "encoder vs decoder prediction");
        assert((ecur.Vl[0] & ~0x3FFFF) == 0);
    }
    assert(voiced > 0);
    mbe_ambe2400EncoderFree(enc);
    puts("uncorrectable C0: direct decoding and the encoder predict from its history");
}

/*
 * No frame is repeated for its count of corrected errors: DVSI's AMBE-3000
 * plays voice with 3 + 3 corrections in C0 and C1 (#155).
 */
static void
test_corrected_errors_play(void) {
    char voice[NVOICE][49];
    make_voice(voice);
    for (int errors = 4; errors <= 6; ++errors) {
        for (int soft = 0; soft <= 1; ++soft) {
            mbe_parms cur, prev, enh, rcur, rprev;
            mbe_process_result r;
            float out[160];
            mbe_initMbeParms(&cur, &prev, &enh);
            for (int i = 0; i < 3; ++i) {
                make_result(&r, 0);
                assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
            }
            rcur = cur;
            rprev = prev;
            assert(mbe_decodeAmbe2400Parms(voice[3], &rcur, &rprev) == 0);
            make_result_c0(&r, 3, errors - 3);
            r.flags |= soft ? MBE_PROCESS_FLAG_SOFT_INPUT : 0u;
            assert(mbe_processAmbe2400Dataf(out, &r, voice[3], &cur, &prev, &enh) == errors);
            assert((r.flags & MBE_PROCESS_FLAG_REPEAT) == 0u);
            assert_same_pcm(&prev.log2Ml[1], &rcur.log2Ml[1], (size_t)rcur.L, "corrected frame's model");
        }
    }
    puts("corrected errors: a frame with corrections plays, hard or soft");
}

/*
 * A tone frame with any correctable errors plays, hard or soft, as DVSI's
 * AMBE-3000 plays it for 0 + 3, 2 + 0, 1 + 2 and 3 + 3 corrections in C0 and
 * C1 (#155): the tone is sample for sample the clean tone frames' (the speech
 * fading out under the first one is smoothed by its error count); an
 * uncorrectable C0 word makes it a lost frame. The library used to play a tone
 * frame with 2 C0 corrections, or 3 in all, as lost.
 */
static void
test_hard_tone_with_corrections(void) {
    char voice[NVOICE][49], tone[49];
    make_voice(voice);
    pack_tone(tone, 126, 32, 180);
    static const int errors[][3] = {{0, 0, 0}, {0, 3, 0}, {2, 0, 0}, {1, 2, 0}, {3, 3, 0}, {2, 4, 1}, {0, 5, 1}};
    float want[2][160], got[2][160];
    for (size_t k = 0; k < sizeof(errors) / sizeof(errors[0]); ++k) {
        mbe_parms cur, prev, enh;
        mbe_process_result r;
        float out[160];
        mbe_setThreadRngSeed(1);
        mbe_initMbeParms(&cur, &prev, &enh);
        for (int i = 0; i < 4; ++i) {
            make_result(&r, 0);
            assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
        }
        for (int t = 0; t < 2; ++t) {
            const unsigned soft = errors[k][2] ? MBE_PROCESS_FLAG_SOFT_INPUT : 0u;
            make_result_c0(&r, errors[k][0], errors[k][1]);
            r.flags |= soft;
            assert(mbe_processAmbe2400Dataf(k == 0 ? want[t] : got[t], &r, tone, &cur, &prev, &enh) >= 0);
            assert(r.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE | soft));
        }
        if (k > 0) {
            /* The first frame's crossfade can round differently in optimized builds,
             * by a few ulps of the tone's peak. */
            float peak = 0.0f;
            for (int n = 80; n < 160; ++n) {
                peak = fmaxf(peak, fabsf(want[0][n]));
            }
            for (int n = 80; n < 160; ++n) {
                if (fabsf(got[0][n] - want[0][n]) > 1e-5f * peak) {
                    printf("tone frame with corrections: sample %d differs (%g vs %g)\n", n, (double)got[0][n],
                           (double)want[0][n]);
                    assert(0);
                }
            }
            assert_same_pcm(got[1], want[1], 160, "second tone frame with corrections");
        }
        make_uncorrectable_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
        assert((r.flags & (MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT))
               == (MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT));
    }
    puts("tones: a tone frame with corrected errors plays, hard or soft");
}

static void lost_after_tone_run(int unvoiced);

/*
 * Lost frames right after a tone, as on DVSI's AMBE-3000 (#156): the first
 * plays the tone going on and fading out over the first 80 samples, then
 * silence; the next two are silent; the 4th and later are mute noise; the
 * prediction history is the one from before the tone.
 */
static void
test_lost_after_tone(void) {
    for (int unvoiced = 0; unvoiced <= 1; ++unvoiced) {
        lost_after_tone_run(unvoiced);
    }
    puts("lost frames after a tone: the tone fades out, then silence, history kept");
}

static void
lost_after_tone_run(int unvoiced) {
    char voice[NVOICE][49], tone[49], b125[49];
    make_voice(voice);
    for (int f = 0; unvoiced && f < NVOICE; ++f) {
        memset(voice[f] + 38, 0, 4); /* b1 0: every band unvoiced */
    }
    pack_tone(tone, 126, 32, 180);
    memcpy(b125, voice[3], sizeof(b125));
    put_b0(b125, 125);
    mbe_parms cur, prev, enh, tone_state, before;
    mbe_process_result r;
    float out[160], plain[160];
    short pcm[160];
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 4; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
    }
    before = prev;
    tone_state = cur;
    for (int i = 0; i < 2; ++i) {
        mbe_synthesizeTonefdstar(plain, tone, &tone_state, 32);
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
    }
    mbe_synthesizeTonefdstar(plain, tone, &tone_state, 32);
    for (int i = 0; i < 6; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, b125, &cur, &prev, &enh) == 0);
        assert((r.flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
        assert(((r.flags & MBE_PROCESS_FLAG_MUTE) != 0u) == (i >= 3));
        if (i == 0) {
            for (int n = 0; n < 160; ++n) {
                const float gain = (n >= 80) ? 1.0f : (float)n / 80.0f;
                assert_close_sample(out[n], (1.0f - gain) * plain[n], "tone tail", n);
            }
        } else if (i < 3) {
            to_int16(out, pcm); /* digital silence */
            for (int n = 0; n < 160; ++n) {
                assert(pcm[n] == 0);
            }
        } else {
            to_int16(out, pcm);
            for (int n = 0; n < 160; ++n) {
                assert(abs(pcm[n]) <= 1);
            }
        }
    }
    assert_same_pcm(&prev.log2Ml[1], &before.log2Ml[1], (size_t)before.L, "history after the tone");
}

/* The error-rate mute plays the same DVSI-level noise: int16 0 or +-1, about a
 * third of the samples nonzero (DVSI's muted dam_e10_hd: 33%). */
static void
test_error_rate_mute_level(void) {
    char voice[NVOICE][49];
    make_voice(voice);
    mbe_parms cur, prev, enh;
    mbe_process_result r;
    float out[160];
    short pcm[160];
    int nonzero = 0, muted = 0;
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 200; ++i) {
        make_result(&r, 3);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i % NVOICE], &cur, &prev, &enh) == 3);
        if ((r.flags & MBE_PROCESS_FLAG_MUTE) == 0u) {
            continue;
        }
        ++muted;
        to_int16(out, pcm);
        for (int n = 0; n < 160; ++n) {
            assert(abs(pcm[n]) <= 1);
            nonzero += pcm[n] != 0;
        }
    }
    assert(muted > 100);
    double share = (double)nonzero / (160.0 * muted);
    assert(share > 0.30 && share < 0.37);
    puts("error-rate mute: DVSI-level noise");
}

int
main(void) {
    test_classes();
    test_null_frame_is_voice();
    test_b0_124_silence_frame();
    test_lost_plays_as_repeat();
    test_mixed_repeat_count();
    test_frame_apis_agree();
    test_zero_volume_tone();
    test_tone_ends_repeat_run();
    test_tone_muted_by_error_rate();
    test_tone_transitions();
    test_tone_transitions_after_mute();
    test_repeat_limit_mute();
    test_error_rate_mute_level();
    test_repeat_limit_under_error_mute();
    test_repeat_limit_noise_continuity();
    test_all_ones_frame();
    test_loud_held_frame_plays();
    test_c0_uncorrectable_detection();
    test_c0_uncorrectable_plays_lost();
    test_c0_uncorrectable_history_kept();
    test_c0_uncorrectable_history_entry_points();
    test_corrected_errors_play();
    test_hard_tone_with_corrections();
    test_lost_after_tone();
    return 0;
}
