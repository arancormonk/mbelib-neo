// SPDX-License-Identifier: GPL-2.0-or-later
/* IMBE 7200x4400 encoder tests: Hamming and frame FEC, pitch tracking
 * through the decoder, level, silence, input validation and context replay. */
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "mbe_ecc.h"
#include "mbelib-neo/mbelib.h"

#define FRAME     160
#define LOOKAHEAD 3 /* frames of encoder delay */

static uint32_t rng = 0x7200u;

static uint32_t
rnd(void) {
    rng = rng * 1664525u + 1013904223u;
    return rng >> 8;
}

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

static double
rms(const short* x, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        s += (double)x[i] * (double)x[i];
    }
    return sqrt(s / (double)n);
}

static int
test_invalid_arguments(mbe_imbe4400_encoder* enc) {
    int fails = 0;
    float samples[FRAME] = {0};
    short shorts[FRAME] = {0};
    char d[88] = {0};
    char fr[8][23];

    fails += mbe_encodeImbe4400Parms(NULL, samples, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeImbe4400Parms(enc, NULL, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeImbe4400Parms(enc, samples, NULL) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeImbe4400ParmsShort(NULL, shorts, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeImbe4400ParmsShort(enc, NULL, d) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeImbe7200x4400Frame(NULL, fr) != MBE_STATUS_INVALID_ARGUMENT;
    fails += mbe_encodeImbe7200x4400Frame(d, NULL) != MBE_STATUS_INVALID_ARGUMENT;
    d[50] = 2;
    fails += mbe_encodeImbe7200x4400Frame(d, fr) != MBE_STATUS_INVALID_BITS;
    mbe_imbe4400EncoderFree(NULL);
    mbe_imbe4400EncoderReset(NULL);
    mbe_imbe4400EncoderReset(enc);
    printf("invalid arguments: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

/* Every data word encodes to a codeword the decoder accepts unchanged, and
 * every single bit error is corrected. */
static int
test_hamming(void) {
    int fails = 0;
    for (int word = 0; word < 2048; word++) {
        char in[11], cw[15], out[15];
        for (int i = 0; i < 11; i++) {
            in[i] = (char)((word >> (10 - i)) & 1);
        }
        mbe_hamming1511_encode(in, cw);
        for (int flip = -1; flip < 15; flip++) {
            char bad[15];
            memcpy(bad, cw, sizeof(bad));
            if (flip >= 0) {
                bad[flip] ^= 1;
            }
            int errs = mbe_hamming1511(bad, out);
            int ok = errs == ((flip >= 0) ? 1 : 0);
            for (int i = 0; i < 11; i++) {
                ok = ok && (out[14 - i] == in[i]);
            }
            fails += !ok;
        }
    }
    printf("hamming (15,11) encode: %s (%d fails)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

static int
test_frame_roundtrip(void) {
    int fails = 0;
    for (int it = 0; it < 2048; it++) {
        char d[88], fr[8][23], out[88];
        mbe_process_result res;
        for (int i = 0; i < 88; i++) {
            d[i] = (char)(rnd() & 1);
        }
        if (mbe_encodeImbe7200x4400Frame(d, fr) != 0) {
            fails++;
            continue;
        }
        int errs = mbe_decodeImbe7200x4400Frame((const char (*)[23])fr, out, &res);
        fails += (errs != 0) || memcmp(d, out, 88) != 0;
        /* One bit error in each protected vector is corrected; modulation of
         * vectors 1-6 is keyed off the corrected u0. The decoder counts a
         * corrected Golay parity bit as no error. */
        static const int lengths[7] = {23, 23, 23, 23, 15, 15, 15};
        int v = it % 7;
        int bit = (int)(rnd() % (uint32_t)lengths[v]);
        fr[v][bit] ^= 1;
        errs = mbe_decodeImbe7200x4400Frame((const char (*)[23])fr, out, &res);
        fails += (errs < 0) || (errs > 1) || memcmp(d, out, 88) != 0;
    }
    printf("frame FEC roundtrip and correction: %s (%d fails)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

/* A steady harmonic signal decodes at its period (or, for a perfectly
 * periodic signal, its double, which fits equally well). */
static int
test_pitch_tracking(mbe_imbe4400_encoder* enc) {
    static const double periods[3] = {30.0, 60.0, 100.0};
    int fails = 0;
    for (int p = 0; p < 3; p++) {
        mbe_parms cur, prev, enhanced;
        int checked = 0;
        mbe_imbe4400EncoderReset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int n = 0; n < 30; n++) {
            short pcm[FRAME];
            char d[88];
            harmonic_frame(pcm, (long)n * FRAME, periods[p], 1500.0);
            fails += mbe_encodeImbe4400ParmsShort(enc, pcm, d) != 0;
            if (mbe_decodeImbe4400Parms(d, &cur, &prev) != 0) {
                fails++;
                continue;
            }
            mbe_moveMbeParms(&cur, &prev);
            if (n >= LOOKAHEAD + 6) {
                double decoded = 2.0 * M_PI / (double)cur.w0;
                int ok = (fabs((decoded / periods[p]) - 1.0) < 0.03)
                         || (periods[p] * 2.0 <= 123.0 && fabs((decoded / (2.0 * periods[p])) - 1.0) < 0.03);
                if (!ok) {
                    fails++;
                    printf("  period %.0f frame %d: decoded %.2f\n", periods[p], n, decoded);
                }
                checked++;
            }
        }
        fails += checked < 15;
    }
    printf("pitch tracking: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

/* Decoded level of a harmonic signal at `level`, over the frames after the
 * encoder's initial look-ahead period. */
static double
decoded_level(mbe_imbe4400_encoder* enc, double level, double* in_level) {
    mbe_parms cur, prev, enhanced;
    double in_sum = 0.0;
    double out_sum = 0.0;
    mbe_imbe4400EncoderReset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int n = 0; n < 40; n++) {
        short pcm[FRAME];
        short out[FRAME];
        char d[88];
        char fr[8][23];
        char d2[88];
        mbe_process_result res;
        harmonic_frame(pcm, (long)n * FRAME, 57.0, level);
        mbe_encodeImbe4400ParmsShort(enc, pcm, d);
        mbe_encodeImbe7200x4400Frame(d, fr);
        memset(&res, 0, sizeof(res));
        mbe_processImbe7200x4400Frame(out, &res, (const char (*)[23])fr, d2, &cur, &prev, &enhanced);
        if (n >= 10) {
            in_sum += rms(pcm, FRAME);
            out_sum += rms(out, FRAME);
        }
    }
    *in_level = in_sum;
    return out_sum;
}

/* The decoded level follows the input level. */
static int
test_level(mbe_imbe4400_encoder* enc) {
    double in_loud, in_quiet;
    double out_loud = decoded_level(enc, 3000.0, &in_loud);
    double out_quiet = decoded_level(enc, 3000.0 / 8.0, &in_quiet);
    double loud_db = 20.0 * log10(out_loud / in_loud);
    double step_db = 20.0 * log10(out_loud / out_quiet);
    int fails = (fabs(loud_db) > 3.0) || (fabs(step_db - 18.06) > 1.5);
    printf("level: %s (decoded %+.2f dB, 18.06 dB step decodes as %.2f dB)\n", fails ? "FAIL" : "ok", loud_db, step_db);
    return fails;
}

static int
test_silence(mbe_imbe4400_encoder* enc) {
    mbe_parms cur, prev, enhanced;
    int fails = 0;
    double worst = 0.0;
    mbe_imbe4400EncoderReset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int n = 0; n < 25; n++) {
        short pcm[FRAME] = {0};
        short out[FRAME];
        char d[88];
        char fr[8][23];
        char d2[88];
        mbe_process_result res;
        fails += mbe_encodeImbe4400ParmsShort(enc, pcm, d) != 0;
        mbe_encodeImbe7200x4400Frame(d, fr);
        memset(&res, 0, sizeof(res));
        mbe_processImbe7200x4400Frame(out, &res, (const char (*)[23])fr, d2, &cur, &prev, &enhanced);
        if (n >= 5) {
            double r = rms(out, FRAME);
            worst = (r > worst) ? r : worst;
        }
    }
    /* Amplitudes are floored at one PCM step per harmonic, as in OP25's
     * imbe_vocoder: about -68 dBFS. */
    fails += worst > 30.0;
    printf("silence: %s (worst frame rms %.2f)\n", fails ? "FAIL" : "ok", worst);
    return fails;
}

static void
stream_frame(short* pcm, int n) {
    for (int i = 0; i < FRAME; i++) {
        long t = ((long)n * FRAME) + i;
        double period = 45.0 + (25.0 * sin((double)t / 3000.0));
        double v = 0.0;
        for (int h = 1; h <= 6; h++) {
            v += 2000.0 / (double)h * sin(2.0 * M_PI * (double)h * (double)t / period);
        }
        v += (double)((int)(rnd() % 801u) - 400);
        pcm[i] = (short)(((n / 8) % 4 == 3) ? lrint(v / 50.0) : lrint(v));
    }
}

static int
test_contexts(mbe_imbe4400_encoder* enc) {
    mbe_imbe4400_encoder* other = mbe_imbe4400EncoderAlloc();
    int fails = 0;
    static char first[60][88];
    if (other == NULL) {
        return 1;
    }
    mbe_imbe4400EncoderReset(enc);
    for (int pass = 0; pass < 2; pass++) {
        rng = 0x5eed;
        mbe_imbe4400EncoderReset(other);
        for (int n = 0; n < 60; n++) {
            short pcm[FRAME];
            float samples[FRAME];
            char a[88];
            char b[88];
            stream_frame(pcm, n);
            for (int i = 0; i < FRAME; i++) {
                samples[i] = (float)pcm[i] / 32768.0f;
            }
            if (n == 20) {
                float bad[FRAME];
                memcpy(bad, samples, sizeof(bad));
                bad[3] = INFINITY;
                fails += mbe_encodeImbe4400Parms(other, bad, b) != MBE_STATUS_INVALID_ARGUMENT;
            }
            fails += mbe_encodeImbe4400ParmsShort(enc, pcm, a) != 0;
            fails += mbe_encodeImbe4400Parms(other, samples, b) != 0;
            fails += memcmp(a, b, 88) != 0 || a[87] != 0;
            if (pass == 0) {
                memcpy(first[n], a, 88);
            } else {
                fails += memcmp(first[n], a, 88) != 0;
            }
        }
        mbe_imbe4400EncoderReset(enc);
    }
    mbe_imbe4400EncoderFree(other);
    printf("short/float parity, contexts and reset replay: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

int
main(void) {
    mbe_imbe4400_encoder* enc = mbe_imbe4400EncoderAlloc();
    if (enc == NULL) {
        return 1;
    }
    int fails = 0;
    fails += test_invalid_arguments(enc);
    fails += test_hamming();
    fails += test_frame_roundtrip();
    fails += test_pitch_tracking(enc);
    fails += test_level(enc);
    fails += test_silence(enc);
    fails += test_contexts(enc);
    mbe_imbe4400EncoderFree(enc);
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}
