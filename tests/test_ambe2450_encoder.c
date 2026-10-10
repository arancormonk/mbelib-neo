// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/* AMBE+2 3600x2450 encoder: frame FEC, the quantizers against the decoder,
 * tone frames, and the codec-neutral checks of encoder_test_support.h. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ambe_encoder.h"
#include "encoder_test_support.h"
#include "mbe_speech_analysis.h"
#include "mbe_tone_detect.h"
#include "mbelib-neo/mbelib.h"

static void*
codec_alloc(void) {
    return mbe_ambe2450EncoderAlloc();
}

static void
codec_reset(void* enc) {
    mbe_ambe2450EncoderReset(enc);
}

static void
codec_release(void* enc) {
    mbe_ambe2450EncoderFree(enc);
}

static int
codec_encode(void* enc, const float* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2450Parms(enc, pcm, bits, cur, prev);
}

static int
codec_encode_short(void* enc, const short* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2450ParmsShort(enc, pcm, bits, cur, prev);
}

static int
codec_decode(const char* bits, mbe_parms* cur, mbe_parms* prev) {
    return mbe_decodeAmbe2450Parms(bits, cur, prev);
}

static int
codec_process(float* out, mbe_process_result* result, const char* bits, mbe_parms* cur, mbe_parms* prev,
              mbe_parms* enhanced) {
    return mbe_processAmbe2450Dataf(out, result, bits, cur, prev, enhanced);
}

/* The Table 9 index a tone frame carries (7.2 Table 10, bits 12..19). */
static int
codec_tone_id(const char* bits) {
    int id = 0;
    if (mbe_classifyAmbe2450Frame(bits) != MBE_AMBE2450_FRAME_TONE) {
        return -1;
    }
    for (int i = 12; i < 20; i++) {
        id = (id << 1) | bits[i];
    }
    return id;
}

/* Repeated, this voice frame takes the decoder's log2Ml to 32.3, about the
 * highest a search over repeated frames found. */
static const struct enc_codec ambe2450_codec = {
    "AMBE+2",
    49,
    MBE_AMBE2450_FRAME_TONE,
    codec_alloc,
    codec_reset,
    codec_release,
    codec_encode,
    codec_encode_short,
    codec_decode,
    codec_process,
    "1100100011111110000110100000010101001100010110110",
    codec_tone_id,
};

static void
random_bits(char d[49]) {
    for (int i = 0; i < 49; i++) {
        d[i] = (char)((enc_rnd() >> 12) & 1); /* a high bit: the LCG's low bits repeat quickly */
    }
}

/* All 49 bits survive the frame FEC, bit 24 included. */
static int
test_frame_roundtrip(void) {
    for (int it = 0; it < 4096; it++) {
        char d[49], fr[4][24], out[49];
        random_bits(d);
        if (mbe_encodeAmbe3600x2450Frame(d, fr) != 0
            || mbe_decodeAmbe3600x2450Frame((const char (*)[24])fr, out, NULL) != 0 || memcmp(d, out, 49) != 0) {
            printf("frame FEC round trip failed at %d\n", it);
            return 1;
        }
    }
    puts("AMBE+2 frame FEC: all 49 bits round-trip (4096 frames)");
    return 0;
}

/* The frame may be written over the parameter bits it encodes. */
static int
test_frame_overlap(void) {
    for (int it = 0; it < 64; it++) {
        char d[49], expected[4][24], out[49];

        union {
            char bits[49];
            char frame[4][24];
        } shared;

        random_bits(d);
        memcpy(shared.bits, d, sizeof(shared.bits));
        if (mbe_encodeAmbe3600x2450Frame(d, expected) != 0
            || mbe_encodeAmbe3600x2450Frame(shared.bits, shared.frame) != 0
            || memcmp(expected, shared.frame, sizeof(expected)) != 0
            || mbe_decodeAmbe3600x2450Frame((const char (*)[24])shared.frame, out, NULL) != 0
            || memcmp(d, out, 49) != 0) {
            puts("AMBE+2 frame FEC: overlapping input and output differ");
            return 1;
        }
    }
    puts("AMBE+2 frame FEC: output may overlap the input bits");
    return 0;
}

/* One error in any protected C0 or C1 position is corrected; the decoder
 * counts every error it corrects in a data position. */
static int
test_frame_single_bit_errors(void) {
    for (int it = 0; it < 64; it++) {
        char d[49], fr[4][24], damaged[4][24], out[49];
        random_bits(d);
        if (mbe_encodeAmbe3600x2450Frame(d, fr) != 0) {
            return 1;
        }
        for (int plane = 0; plane < 2; plane++) {
            for (int bit = 0; bit < 24 - plane; bit++) {
                memcpy(damaged, fr, sizeof(fr));
                damaged[plane][bit] ^= 1;
                int errors = mbe_decodeAmbe3600x2450Frame((const char (*)[24])damaged, out, NULL);
                int data = bit >= 12 - plane; /* C0 data at 12..23, C1 data at 11..22 */
                if ((data ? errors != 1 : (errors < 0 || errors > 1)) || memcmp(d, out, 49) != 0) {
                    printf("single-bit error at plane %d bit %d not corrected\n", plane, bit);
                    return 1;
                }
            }
        }
    }
    puts("AMBE+2 frame FEC: every single C0/C1 error corrected");
    return 0;
}

/* b0..b8 of a voice frame (TIA-102.BABA-1 Tables 5-8). */
static int
field(const char d[49], const int* positions, int count) {
    int value = 0;
    for (int i = 0; i < count; i++) {
        value = (value << 1) | d[positions[i]];
    }
    return value;
}

static int
frame_b0(const char d[49]) {
    static const int positions[7] = {0, 1, 2, 3, 37, 38, 39};
    return field(d, positions, 7);
}

static int
frame_b1(const char d[49]) {
    static const int positions[5] = {4, 5, 6, 7, 35};
    return field(d, positions, 5);
}

/*
 * Quantizing a decoded model again gives the same model: random voice frames
 * are decoded, their log2 amplitudes requantized at the same b0 and b1 against
 * the same history, and the result decoded again. This checks the gain,
 * prediction, block DCTs, PRBA and higher-order quantizers and the packing
 * independently of the analysis.
 */
static int
test_requantize(void) {
    int frames = 0;
    double worst = 0.0;
    mbe_parms prev, cur, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int it = 0; it < 3000; it++) {
        char d[49], again[49];
        mbe_parms decoded = cur, history = prev, redecoded = cur, history2 = prev;
        random_bits(d);
        int b0 = frame_b0(d);
        if (b0 >= 120 || frame_b1(d) > 16 || (d[0] & d[1] & d[2] & d[3] & d[4] & d[5])) {
            continue;
        }
        if (mbe_decodeAmbe2450Parms(d, &decoded, &history) != MBE_AMBE2450_FRAME_VOICE) {
            return 1;
        }
        mbe_ambe2450_quantize_amplitudes(b0, frame_b1(d), decoded.log2Ml, &prev, again);
        if (mbe_decodeAmbe2450Parms(again, &redecoded, &history2) != MBE_AMBE2450_FRAME_VOICE
            || redecoded.L != decoded.L || !enc_float_bits_equal(redecoded.gamma, decoded.gamma)) {
            printf("requantize: frame %d changed L or gamma\n", it);
            return 1;
        }
        for (int l = 1; l <= decoded.L; l++) {
            worst = fmax(worst, fabs((double)redecoded.log2Ml[l] - (double)decoded.log2Ml[l]));
            if (redecoded.Vl[l] != decoded.Vl[l]) {
                return 1;
            }
        }
        frames++;
        prev = decoded; /* walk the history through the decoded frames */
    }
    printf("AMBE+2 requantization of %d decoded frames: worst log2 amplitude change %.2g\n", frames, worst);
    return frames < 1000 || worst > 1e-3;
}

/* eq 4 weighs disagreement by squared amplitude: with voicing 111000111 over
 * nine harmonics the squared-error choice differs from a linear one. */
static int
test_vuv_quantizer(void) {
    unsigned char bands[MBE_ANALYSIS_BANDS] = {1, 0, 1};
    const float m[57] = {0, 10, 2, 2, 2, 10, 3, 2, 1, 1};
    int b1 = mbe_ambe2450_quantize_vuv(0, m, bands);
    /* Independent evaluation of eq 4 over the first 17 vectors. */
    int best = -1;
    double best_distance = 0.0;
    for (int n = 0; n < 17; n++) {
        char d[49] = {0};
        mbe_parms cur, prev, enhanced;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int i = 0; i < 4; i++) {
            d[4 + i] = (char)((n >> (4 - i)) & 1);
        }
        d[35] = (char)(n & 1);
        if (mbe_decodeAmbe2450Parms(d, &cur, &prev) != MBE_AMBE2450_FRAME_VOICE || cur.L != 9) {
            return 1;
        }
        double distance = 0.0;
        for (int l = 1; l <= 9; l++) {
            int voiced = bands[(l - 1) / 3];
            distance += (cur.Vl[l] != voiced) ? (double)m[l] * (double)m[l] : 0.0;
        }
        if (best < 0 || distance < best_distance) {
            best = n;
            best_distance = distance;
        }
    }
    printf("AMBE+2 V/UV quantizer: b1 %d, eq 4 by the decoder's vectors %d\n", b1, best);
    return b1 != best;
}

struct tone_case {
    double f1, f2;
    int id;
    const char* name;
};

/* A frame-aligned tone, continuous over frames, at the given phase. */
static void
tone_frame(float pcm[160], int frame, const struct tone_case* t, double amplitude, double phase) {
    for (int i = 0; i < 160; i++) {
        double n = (double)((frame * 160) + i);
        double v = sin((2.0 * M_PI * t->f1 * n / 8000.0) + phase);
        if (t->f2 > 0.0) {
            v += sin((2.0 * M_PI * t->f2 * n / 8000.0) + (2.0 * phase));
        }
        pcm[i] = (float)(amplitude * v);
    }
}

static int
tone_id_of(const char d[49]) {
    int id = 0;
    for (int i = 12; i < 20; i++) {
        id = (id << 1) | d[i];
    }
    return id;
}

/*
 * Encode 25 frames of a tone; every frame after the first two must be a tone
 * frame with the right index and a self-consistent Table 10 layout, cur_mp a
 * copy of prev_mp, and the decoder must play it at about the input level. A
 * call-progress tone fills the span from frame 1, so it is voice until it has
 * filled MBE_TONE_CP_CONFIRM of them.
 */
static int
run_tone(void* enc, const struct tone_case* t, double amplitude, double phase, double* level_db) {
    mbe_parms ec, ep, eh, dc, dp, dh;
    double sum = 0.0;
    int n = 0;
    mbe_ambe2450EncoderReset(enc);
    mbe_initMbeParms(&ec, &ep, &eh);
    mbe_initMbeParms(&dc, &dp, &dh);
    for (int frame = 0; frame < 25; frame++) {
        float pcm[160], out[160];
        char d[49];
        mbe_process_result result;
        tone_frame(pcm, frame, t, amplitude, phase);
        int r = mbe_encodeAmbe2450Parms(enc, pcm, d, &ec, &ep);
        mbe_initProcessResult(&result);
        if (mbe_processAmbe2450Dataf(out, &result, d, &dc, &dp, &dh) < 0) {
            return 1;
        }
        const int first = (t->id >= 160) ? MBE_TONE_CP_CONFIRM : 2;
        if (frame < first && t->id >= 160 && r != MBE_AMBE2450_FRAME_VOICE) {
            printf("  %s at phase %.2f: frame %d returned %d before confirmation\n", t->name, phase, frame, r);
            return 1;
        }
        if (frame >= first) {
            if (r != MBE_AMBE2450_FRAME_TONE || tone_id_of(d) != t->id || !enc_parms_identical(&ec, &ep)
                || mbe_classifyAmbe2450Frame(d) != MBE_AMBE2450_FRAME_TONE || !(result.flags & MBE_PROCESS_FLAG_TONE)) {
                printf("  %s at phase %.2f: frame %d returned %d, id %d\n", t->name, phase, frame, r, tone_id_of(d));
                return 1;
            }
            for (int i = 0; frame >= 5 && i < 160; i++) {
                double x = 7.0 * (double)out[i];
                sum += x * x;
                n++;
            }
        }
        mbe_moveMbeParms(&ec, &ep);
    }
    *level_db = 10.0 * log10((sum / n) + 1e-30);
    return 0;
}

/* DTMF, KNOX, call-progress and single tones at several phases become tone
 * frames with their Table 9 index; the decoder plays each at the input's
 * per-component level within 1.5 dB. */
static int
test_tones(void* enc) {
    static const struct tone_case cases[] = {
        {697, 1209, 129, "DTMF 1"},  {941, 1336, 128, "DTMF 0"},   {941, 1633, 141, "DTMF D"},
        {852, 1477, 137, "DTMF 9"},  {606, 1052, 145, "KNOX 1"},   {820, 1279, 159, "KNOX #"},
        {350, 440, 160, "dial"},     {440, 480, 161, "ringback"},  {480, 620, 162, "busy"},
        {156.25, 0, 5, "156.25 Hz"}, {1000, 0, 32, "1000 Hz"},     {1014, 0, 32, "1014 Hz (off grid)"},
        {3500, 0, 112, "3500 Hz"},   {406.25, 0, 13, "406.25 Hz"}, {468.75, 0, 15, "468.75 Hz"},
    };
    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        for (int p = 0; p < 4; p++) {
            const double amplitude = 3000.0 / 32768.0;
            double level_db;
            if (run_tone(enc, &cases[k], amplitude, 0.7 * p, &level_db) != 0) {
                return 1;
            }
#ifdef MBELIB_TEST_NOTONES
            if (level_db > 0.0) {
                printf("  %s: NOTONES output not silent\n", cases[k].name);
                return 1;
            }
#else
            double components = (cases[k].f2 > 0.0) ? 2.0 : 1.0;
            double input_db = 10.0 * log10(components * 3000.0 * 3000.0 / 2.0);
            if (fabs(level_db - input_db) > 1.5) {
                printf("  %s: decoded at %.1f dB, input %.1f dB\n", cases[k].name, level_db, input_db);
                return 1;
            }
#endif
        }
    }
    puts("AMBE+2 tones: DTMF, KNOX, call progress and single tones sent as tone frames at the input level");
    return 0;
}

/* DTMF with 8 dB of twist or 1.5% off frequency is sent as its digit from the
 * second frame on; 4% off, or one component 20 dB below the other, never
 * becomes a dual-tone frame (the dominant component alone may still be a
 * single tone). */
static int
test_tone_limits(void* enc) {
    const struct {
        double f1, f2, a2;
        int want; /* the digit's index, or -1 for no dual-tone frame */
    } cases[] = {{697, 1209, 0.4, 129},
                 {697 * 1.015, 1209 * 0.985, 1.0, 129},
                 {697 * 1.04, 1209 * 1.04, 1.0, -1},
                 {697, 1209, 0.1, -1}};

    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        mbe_parms ec, ep, eh;
        int expected = 0, dual = 0;
        mbe_ambe2450EncoderReset(enc);
        mbe_initMbeParms(&ec, &ep, &eh);
        for (int frame = 0; frame < 20; frame++) {
            float pcm[160];
            char d[49];
            for (int i = 0; i < 160; i++) {
                double n = (double)((frame * 160) + i);
                pcm[i] = (float)(0.1
                                 * (sin(2.0 * M_PI * cases[k].f1 * n / 8000.0)
                                    + (cases[k].a2 * sin(2.0 * M_PI * cases[k].f2 * n / 8000.0))));
            }
            int r = mbe_encodeAmbe2450Parms(enc, pcm, d, &ec, &ep);
            if (r < 0) {
                return 1;
            }
            int tone = (r == MBE_AMBE2450_FRAME_TONE) ? tone_id_of(d) : -1;
            expected += (frame >= 1 && tone == cases[k].want);
            dual += (tone >= 128);
            mbe_moveMbeParms(&ec, &ep);
        }
        if ((cases[k].want > 0) ? expected != 19 : dual != 0) {
            printf("  DTMF limits case %zu: %d frames as the digit, %d dual-tone frames\n", k, expected, dual);
            return 1;
        }
    }
    puts("AMBE+2 tones: DTMF with 8 dB twist or 1.5% off sent as its digit; 4% off and 20 dB twist never dual");
    return 0;
}

/* Voice, then a tone, then voice again: the voice after the tone predicts from
 * the voice before it, exactly as the decoder does. */
static int
test_voice_tone_voice(void* enc) {
    static const struct tone_case dtmf = {770, 1336, 133, "DTMF 5"};
    mbe_parms ec, ep, eh, dc, dp, dh;
    static short speech[40 * 160];
    enc_speech_like(speech, 40 * 160);
    mbe_ambe2450EncoderReset(enc);
    mbe_initMbeParms(&ec, &ep, &eh);
    mbe_initMbeParms(&dc, &dp, &dh);
    for (int frame = 0; frame < 40; frame++) {
        float pcm[160];
        char d[49];
        if (frame >= 15 && frame < 25) {
            tone_frame(pcm, frame, &dtmf, 0.1, 0.3);
        } else {
            for (int i = 0; i < 160; i++) {
                pcm[i] = (float)speech[(frame * 160) + i] / 32768.0f;
            }
        }
        int r = mbe_encodeAmbe2450Parms(enc, pcm, d, &ec, &ep);
        int kind = mbe_decodeAmbe2450Parms(d, &dc, &dp);
        if (r < 0 || r != kind) {
            return 1;
        }
        if (r == MBE_AMBE2450_FRAME_VOICE) {
            for (int l = 1; l <= dc.L; l++) {
                if (!enc_float_bits_equal(dc.log2Ml[l], ec.log2Ml[l])) {
                    printf("voice-tone-voice: frame %d diverged from the decoder\n", frame);
                    return 1;
                }
            }
            mbe_moveMbeParms(&dc, &dp);
        }
        mbe_moveMbeParms(&ec, &ep);
    }
    puts("AMBE+2 voice-tone-voice: prediction stays with the decoder's");
    return 0;
}

/* A tone frame returns the history as its model, so a history whose model is
 * not finite (w0 or an amplitude) is rejected there too, with cur_mp and the
 * stream untouched. The rejected frame is louder than the one retried, so an
 * analysis it left behind would change the stream. */
static int
test_tone_bad_history(void* enc) {
    static const struct tone_case tone = {1000, 0, 32, "1000 Hz"};
    for (int b = 0; b < 3; b++) {
        char reference[12][49];
        char observed[12][49];
        for (int pass = 0; pass < 2; pass++) {
            mbe_parms cur, prev, enhanced;
            mbe_ambe2450EncoderReset(enc);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            for (int f = 0; f < 12; f++) {
                float pcm[160];
                tone_frame(pcm, f, &tone, 0.1, 0.0);
                if (pass == 1 && f == 6) {
                    mbe_parms bad = prev;
                    const mbe_parms before = cur;
                    char unused[49];
                    float louder[160];
                    tone_frame(louder, f, &tone, 0.11, 0.0);
                    if (b == 0) {
                        bad.w0 = NAN;
                    } else if (b == 1) {
                        bad.Ml[1] = INFINITY;
                    } else {
                        bad.log2Ml[2] = -INFINITY;
                    }
                    if (mbe_encodeAmbe2450Parms(enc, louder, unused, &cur, &bad) != MBE_STATUS_INVALID_ARGUMENT
                        || !enc_parms_identical(&cur, &before)) {
                        printf("tone bad history %d was not rejected cleanly\n", b);
                        return 1;
                    }
                }
                int r = mbe_encodeAmbe2450Parms(enc, pcm, (pass == 0) ? reference[f] : observed[f], &cur, &prev);
                if (r < 0 || (f >= 2 && r != MBE_AMBE2450_FRAME_TONE)) {
                    printf("tone bad history %d: frame %d returned %d\n", b, f, r);
                    return 1;
                }
                mbe_moveMbeParms(&cur, &prev);
            }
        }
        if (memcmp(reference, observed, sizeof(reference)) != 0) {
            printf("tone bad history %d changed the stream\n", b);
            return 1;
        }
    }
    puts("AMBE+2 tones: a history with a non-finite model is rejected without touching state");
    return 0;
}

static int
test_invalid_arguments(void* enc) {
    float pcm[160] = {0};
    short shorts[160] = {0};
    char d[49] = {0}, fr[4][24] = {{0}};
    mbe_parms c, p, h;
    mbe_initMbeParms(&c, &p, &h);
    const int cases[] = {
        mbe_encodeAmbe2450Parms(NULL, pcm, d, &c, &p),
        mbe_encodeAmbe2450ParmsShort(NULL, shorts, d, &c, &p),
        mbe_encodeAmbe2450Parms(enc, NULL, d, &c, &p),
        mbe_encodeAmbe2450Parms(enc, pcm, NULL, &c, &p),
        mbe_encodeAmbe2450Parms(enc, pcm, d, NULL, &p),
        mbe_encodeAmbe2450Parms(enc, pcm, d, &c, NULL),
        mbe_encodeAmbe2450Parms(enc, pcm, d, &p, &p),
        mbe_encodeAmbe2450ParmsShort(enc, NULL, d, &c, &p),
        mbe_encodeAmbe2450ParmsShort(enc, shorts, NULL, &c, &p),
        mbe_encodeAmbe2450ParmsShort(enc, shorts, d, NULL, &p),
        mbe_encodeAmbe2450ParmsShort(enc, shorts, d, &c, NULL),
        mbe_encodeAmbe3600x2450Frame(NULL, fr),
        mbe_encodeAmbe3600x2450Frame(d, NULL),
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        if (cases[i] != MBE_STATUS_INVALID_ARGUMENT) {
            printf("invalid arguments: case %zu returned %d\n", i, cases[i]);
            return 1;
        }
    }
    for (int i = 0; i < 49; i++) {
        d[i] = 2;
        if (mbe_encodeAmbe3600x2450Frame(d, fr) != MBE_STATUS_INVALID_BITS) {
            return 1;
        }
        d[i] = 0;
    }
    puts("AMBE+2 invalid arguments and bits rejected");
    return 0;
}

int
main(void) {
    void* enc = mbe_ambe2450EncoderAlloc();
    if (enc == NULL) {
        return 1;
    }
    int fails = 0;
    fails += test_invalid_arguments(enc);
    fails += test_frame_roundtrip();
    fails += test_frame_single_bit_errors();
    fails += test_frame_overlap();
    fails += test_requantize();
    fails += test_vuv_quantizer();
    fails += enc_run_common(&ambe2450_codec, enc);
    fails += test_tones(enc);
    fails += test_tone_limits(enc);
    fails += test_voice_tone_voice(enc);
    fails += test_tone_bad_history(enc);
    fails += enc_test_call_progress(&ambe2450_codec, enc, 160);
    fails += enc_test_tone_edges(&ambe2450_codec, enc, 137);
    mbe_ambe2450EncoderFree(enc);
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}
