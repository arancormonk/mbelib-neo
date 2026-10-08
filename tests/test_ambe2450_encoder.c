// SPDX-License-Identifier: GPL-2.0-or-later
/* AMBE+2 3600x2450 encoder tests: frame FEC, pitch tracking through the
 * decoder, tone frames, silence, input validation and context replay. */
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "mbelib-neo/mbelib.h"

#define FRAME     160
#define LOOKAHEAD 3 /* frames of encoder delay */

static uint32_t rng = 0x2450u;

static uint32_t
rnd(void) {
    rng = rng * 1664525u + 1013904223u;
    return rng >> 8;
}

/* A harmonic-rich periodic signal (16-bit scale) at a known period. */
static void
harmonic_frame(short* out, long start, double period, double level) {
    for (int i = 0; i < FRAME; i++) {
        double t = (double)(start + i);
        double v = 0.0;
        for (int h = 1; h <= 8; h++) {
            v += level / (double)h * sin(2.0 * M_PI * (double)h * t / period);
        }
        out[i] = (short)lrint(v);
    }
}

static void
tone_frame(short* out, long start, double f1, double f2, double amplitude) {
    for (int i = 0; i < FRAME; i++) {
        double t = (double)(start + i) / 8000.0;
        double v = amplitude * sin(2.0 * M_PI * f1 * t);
        if (f2 > 0.0) {
            v += amplitude * sin(2.0 * M_PI * f2 * t);
        }
        out[i] = (short)lrint(v);
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

static double
rms(const short* x, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        s += (double)x[i] * (double)x[i];
    }
    return sqrt(s / (double)n);
}

static int
test_invalid_arguments(mbe_ambe2450_encoder* enc) {
    int fails = 0;
    float samples[FRAME] = {0};
    short shorts[FRAME] = {0};
    char d[49] = {0};
    char fr[4][24];

    fails += mbe_encodeAmbe2450Parms(NULL, samples, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeAmbe2450Parms(enc, NULL, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeAmbe2450Parms(enc, samples, NULL) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeAmbe2450ParmsShort(NULL, shorts, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeAmbe2450ParmsShort(enc, NULL, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeAmbe3600x2450Frame(NULL, fr) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeAmbe3600x2450Frame(d, NULL) != MBE_STATUS_INVALID_ARGUMENT;
    d[7] = 2;
    fails += mbe_encodeAmbe3600x2450Frame(d, fr) != MBE_STATUS_INVALID_BITS;
    mbe_ambe2450EncoderFree(NULL);
    mbe_ambe2450EncoderReset(NULL);
    mbe_ambe2450EncoderReset(enc);
    printf("invalid arguments: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

static int
test_frame_roundtrip(void) {
    int fails = 0;
    for (int it = 0; it < 4096; it++) {
        char d[49], fr[4][24], out[49];
        mbe_process_result res;
        for (int i = 0; i < 49; i++) {
            d[i] = (char)(rnd() & 1);
        }
        if (mbe_encodeAmbe3600x2450Frame(d, fr) != 0) {
            fails++;
            continue;
        }
        int errs = mbe_decodeAmbe3600x2450Frame((const char (*)[24])fr, out, &res);
        if (errs != 0 || memcmp(d, out, 49) != 0) {
            fails++;
        }
    }
    printf("frame FEC roundtrip: %s (%d fails / 4096)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

/* A single bit error in either Golay codeword (or the C0 parity bit) is
 * corrected. The decoder counts corrected C1 parity bits as no error. */
static int
test_frame_single_bit_errors(void) {
    int fails = 0;
    for (int it = 0; it < 64; it++) {
        char d[49], fr[4][24], out[49];
        mbe_process_result res;
        for (int i = 0; i < 49; i++) {
            d[i] = (char)(rnd() & 1);
        }
        mbe_encodeAmbe3600x2450Frame(d, fr);
        for (int plane = 0; plane < 2; plane++) {
            for (int bit = 0; bit < ((plane == 0) ? 24 : 23); bit++) {
                char bad[4][24];
                memcpy(bad, fr, sizeof(bad));
                bad[plane][bit] ^= 1;
                int errs = mbe_decodeAmbe3600x2450Frame((const char (*)[24])bad, out, &res);
                if (errs < 0 || errs > 1 || memcmp(d, out, 49) != 0) {
                    fails++;
                }
            }
        }
    }
    printf("frame single-bit correction: %s (%d fails)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

/* A steady harmonic signal decodes at its period with a plausible level. */
static int
test_pitch_tracking(mbe_ambe2450_encoder* enc) {
    static const double periods[3] = {25.0, 50.0, 90.0};
    int fails = 0;
    for (int p = 0; p < 3; p++) {
        mbe_parms cur, prev, enhanced;
        int checked = 0;
        mbe_ambe2450EncoderReset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int n = 0; n < 30; n++) {
            short pcm[FRAME];
            char d[49];
            harmonic_frame(pcm, (long)n * FRAME, periods[p], 1500.0);
            if (mbe_encodeAmbe2450ParmsShort(enc, pcm, d) != MBE_AMBE2450_FRAME_VOICE) {
                fails++;
                continue;
            }
            if (mbe_decodeAmbe2450Parms(d, &cur, &prev) != MBE_AMBE2450_FRAME_VOICE) {
                fails++;
                continue;
            }
            mbe_moveMbeParms(&cur, &prev);
            if (n >= LOOKAHEAD + 6) {
                double decoded = 2.0 * M_PI / (double)cur.w0;
                float peak = 0.0f;
                for (int l = 1; l <= cur.L && l <= 8; l++) {
                    peak = (cur.Ml[l] > peak) ? cur.Ml[l] : peak;
                }
                if (fabs((decoded / periods[p]) - 1.0) > 0.05 || peak < 100.0f) {
                    fails++;
                    printf("  period %.0f frame %d: decoded %.2f peak %.1f\n", periods[p], n, decoded, (double)peak);
                }
                checked++;
            }
        }
        fails += checked < 15;
    }
    printf("pitch tracking: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

/* DTMF digits and single tones become tone frames the decoder plays as tones. */
static int
test_tones(mbe_ambe2450_encoder* enc) {
    struct {
        double f1, f2;
        int id;
    } cases[4] = {{770.0, 1336.0, 0x85}, {941.0, 1336.0, 0x80}, {941.0, 1477.0, 0x8F}, {1000.0, 0.0, 32}};

    int fails = 0;
    for (int c = 0; c < 4; c++) {
        mbe_parms cur, prev, enhanced;
        mbe_ambe2450EncoderReset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        double in_rms = 0.0;
        double out_rms = 0.0;
        for (int n = 0; n < 12; n++) {
            short pcm[FRAME];
            short out[FRAME];
            char d[49];
            char fr[4][24];
            char d2[49];
            mbe_process_result res;
            tone_frame(pcm, (long)n * FRAME, cases[c].f1, cases[c].f2, 4000.0);
            int kind = mbe_encodeAmbe2450ParmsShort(enc, pcm, d);
            mbe_encodeAmbe3600x2450Frame(d, fr);
            memset(&res, 0, sizeof(res));
            mbe_processAmbe3600x2450Frame(out, &res, (const char (*)[24])fr, d2, &cur, &prev, &enhanced);
            if (n < LOOKAHEAD) {
                continue;
            }
            if (kind != MBE_AMBE2450_FRAME_TONE || mbe_classifyAmbe2450Frame(d) != MBE_AMBE2450_FRAME_TONE
                || tone_id_of(d) != cases[c].id || (res.flags & MBE_PROCESS_FLAG_TONE) == 0u) {
                fails++;
                printf("  tone %.0f/%.0f frame %d: kind %d id %d flags 0x%x\n", cases[c].f1, cases[c].f2, n, kind,
                       tone_id_of(d), res.flags);
            }
            if (n >= LOOKAHEAD + 2) {
                in_rms += rms(pcm, FRAME);
                out_rms += rms(out, FRAME);
            }
        }
#ifdef MBELIB_TEST_NOTONES
        /* NOTONES decoders play tone frames as silence. */
        if (out_rms > 0.0) {
            fails++;
            printf("  tone %.0f/%.0f: NOTONES output not silent\n", cases[c].f1, cases[c].f2);
        }
        (void)in_rms;
#else
        /* The level field follows the decoder's level law to well within 2 dB. */
        double delta_db = 20.0 * log10(out_rms / in_rms);
        if (fabs(delta_db) > 2.0) {
            fails++;
            printf("  tone %.0f/%.0f: decoded level %+.2f dB\n", cases[c].f1, cases[c].f2, delta_db);
        }
#endif
    }
    printf("tone frames: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

/* Digital silence decodes to (near) silence. */
static int
test_silence(mbe_ambe2450_encoder* enc) {
    mbe_parms cur, prev, enhanced;
    int fails = 0;
    double worst = 0.0;
    mbe_ambe2450EncoderReset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int n = 0; n < 25; n++) {
        short pcm[FRAME] = {0};
        short out[FRAME];
        char d[49];
        char fr[4][24];
        char d2[49];
        mbe_process_result res;
        fails += mbe_encodeAmbe2450ParmsShort(enc, pcm, d) != MBE_AMBE2450_FRAME_VOICE;
        mbe_encodeAmbe3600x2450Frame(d, fr);
        memset(&res, 0, sizeof(res));
        mbe_processAmbe3600x2450Frame(out, &res, (const char (*)[24])fr, d2, &cur, &prev, &enhanced);
        if (n >= 5) {
            double r = rms(out, FRAME);
            worst = (r > worst) ? r : worst;
        }
    }
    fails += worst > 10.0;
    printf("silence: %s (worst frame rms %.2f)\n", fails ? "FAIL" : "ok", worst);
    return fails;
}

/* A varied test stream: a gliding harmonic signal with noise and gaps. */
static void
stream_frame(short* pcm, int n) {
    for (int i = 0; i < FRAME; i++) {
        long t = ((long)n * FRAME) + i;
        double period = 40.0 + (20.0 * sin((double)t / 4000.0));
        double v = 0.0;
        for (int h = 1; h <= 6; h++) {
            v += 2000.0 / (double)h * sin(2.0 * M_PI * (double)h * (double)t / period);
        }
        v += (double)((int)(rnd() % 801u) - 400);
        pcm[i] = (short)(((n / 8) % 4 == 3) ? lrint(v / 50.0) : lrint(v));
    }
}

/* Short and float input agree, two contexts agree, reset replays the stream,
 * and a rejected frame leaves the context unchanged. */
static int
test_contexts(mbe_ambe2450_encoder* enc) {
    mbe_ambe2450_encoder* other = mbe_ambe2450EncoderAlloc();
    int fails = 0;
    static char first[60][49];
    if (other == NULL) {
        return 1;
    }
    mbe_ambe2450EncoderReset(enc);
    for (int pass = 0; pass < 2; pass++) {
        rng = 0x5eed;
        mbe_ambe2450EncoderReset(other);
        for (int n = 0; n < 60; n++) {
            short pcm[FRAME];
            float samples[FRAME];
            char a[49];
            char b[49];
            stream_frame(pcm, n);
            for (int i = 0; i < FRAME; i++) {
                samples[i] = (float)pcm[i] / 32768.0f;
            }
            if (n == 20) {
                float bad[FRAME];
                memcpy(bad, samples, sizeof(bad));
                bad[17] = NAN;
                fails += mbe_encodeAmbe2450Parms(other, bad, b) != MBE_STATUS_INVALID_ARGUMENT;
                bad[17] = 3.0e6f;
                fails += mbe_encodeAmbe2450Parms(other, bad, b) != MBE_STATUS_INVALID_ARGUMENT;
            }
            int ka = mbe_encodeAmbe2450ParmsShort(enc, pcm, a);
            int kb = mbe_encodeAmbe2450Parms(other, samples, b);
            fails += (ka < 0) || (ka != kb) || memcmp(a, b, 49) != 0;
            if (pass == 0) {
                memcpy(first[n], a, 49);
            } else {
                fails += memcmp(first[n], a, 49) != 0;
            }
        }
        mbe_ambe2450EncoderReset(enc);
    }
    mbe_ambe2450EncoderFree(other);
    printf("short/float parity, contexts and reset replay: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

int
main(void) {
    mbe_ambe2450_encoder* enc = mbe_ambe2450EncoderAlloc();
    if (enc == NULL) {
        return 1;
    }
    int fails = 0;
    fails += test_invalid_arguments(enc);
    fails += test_frame_roundtrip();
    fails += test_frame_single_bit_errors();
    fails += test_pitch_tracking(enc);
    fails += test_tones(enc);
    fails += test_silence(enc);
    fails += test_contexts(enc);
    mbe_ambe2450EncoderFree(enc);
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}
