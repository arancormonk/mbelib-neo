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
 * DTX frame keeps its silence handling.
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
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

/* A clean frame result, or one that reports `errors` protected-bit corrections. */
static void
make_result(mbe_process_result* r, int errors) {
    mbe_initProcessResult(r);
    r->flags = MBE_PROCESS_FLAG_C0_VALID;
    r->protected_errors = errors;
    r->total_errors = errors;
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

/* The D-STAR null frame decodes cleanly to b0 124 and is voice. */
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
    assert(cur.L == mbe_ambe2400_harmonic_count(mbe_ambe2400_f0(124)));
    puts("null frame: b0 124, decoded as voice");
}

/*
 * Decode voice[0..3], `count` copies of `special`, then voice[4..5]. With
 * `errors` > 3 the special frames are instead voice[3] reported with that many
 * corrections, which takes the existing repeat path. Writes all PCM and the
 * special frames' flags.
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
        make_result(&r, special ? 0 : errors);
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
        make_result(&r, errored ? 4 : 0);
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
    const char* seq[] = {voice[0], voice[1], voice[2], b125, unused, voice[4], voice[5]};
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
        make_result(&r[0], 0);
        assert(mbe_processAmbe2400Dataf(out[0], &r[0], d, &c[0], &p[0], &e[0]) >= 0);
        assert(mbe_processAmbe3600x2400Framef(out[1], &r[1], (const char (*)[24])fr, dh, &c[1], &p[1], &e[1]) >= 0);
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
                make_result(&r, errored ? 4 : 0);
                assert(mbe_processAmbe2400Dataf(out, &r, errored ? voice[1] : b125, &cur, &prev, &enh) >= 0);
                assert((r.flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
            }
            make_result(&r, 0);
            assert(mbe_processAmbe2400Dataf(out, &r, tone, &cur, &prev, &enh) == 0);
            assert(r.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE));
            assert(prev.repeatCount == 0);
            make_result(&r, errored ? 4 : 0);
            assert(mbe_processAmbe2400Dataf(out, &r, errored ? voice[1] : b125, &cur, &prev, &enh) >= 0);
            assert((r.flags & (MBE_PROCESS_FLAG_REPEAT | MBE_PROCESS_FLAG_MUTE)) == MBE_PROCESS_FLAG_REPEAT);
            assert(prev.repeatCount == 1);
        }
    }
    puts("repeat count: a played tone ends a run of repeats");
}

/* The all-ones DTX frame keeps the silence handling: comfort noise, then a
 * decoder reset, so the next voice frame decodes as from a fresh stream. */
static void
test_dtx_frame(void) {
    unsigned char ones[9];
    char fr[4][24], d[49];
    char voice[NVOICE][49];
    mbe_process_result r;
    mbe_parms cur, prev, enh, fcur, fprev, fenh;
    float out[160], after[160], fresh[160];
    memset(ones, 0xFF, sizeof(ones));
    make_voice(voice);
    assert(mbe_decodeDStarDVData(ones, fr) == 0);
    assert(mbe_decodeAmbe3600x2400Frame((const char (*)[24])fr, d, &r) == 2);
    mbe_initMbeParms(&cur, &prev, &enh);
    assert(mbe_decodeAmbe2400Parms(d, &cur, &prev) == MBE_AMBE2400_FRAME_SILENCE);
    mbe_initMbeParms(&cur, &prev, &enh);
    for (int i = 0; i < 4; ++i) {
        make_result(&r, 0);
        assert(mbe_processAmbe2400Dataf(out, &r, voice[i], &cur, &prev, &enh) == 0);
    }
    assert(mbe_processAmbe3600x2400Framef(out, &r, (const char (*)[24])fr, d, &cur, &prev, &enh) == 2);
    assert(r.flags == (MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_TONE));
    mbe_setThreadRngSeed(2);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(after, &r, voice[4], &cur, &prev, &enh) == 0);
    mbe_initMbeParms(&fcur, &fprev, &fenh);
    mbe_setThreadRngSeed(2);
    make_result(&r, 0);
    assert(mbe_processAmbe2400Dataf(fresh, &r, voice[4], &fcur, &fprev, &fenh) == 0);
    assert_same_pcm(after, fresh, 160, "voice after the DTX frame vs fresh stream");
    puts("DTX frame: silence handling and decoder reset kept");
}

int
main(void) {
    test_classes();
    test_null_frame_is_voice();
    test_lost_plays_as_repeat();
    test_mixed_repeat_count();
    test_frame_apis_agree();
    test_zero_volume_tone();
    test_tone_ends_repeat_run();
    test_dtx_frame();
    return 0;
}
