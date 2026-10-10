// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Micro-benchmark for the D-STAR, AMBE+2 and IMBE encoders.
 *
 * Encodes a deterministic speech-like signal (a harmonic series whose pitch
 * glides between 90 and 260 Hz, with a slowly varying spectral tilt and some
 * noise) with each encoder and reports CPU time per 20 ms frame across
 * repeated runs.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FRAMES 1500 /* 30 s per run */
#define RUNS   7

static double
secs(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

static short pcm[FRAMES * 160];

static void
make_signal(void) {
    uint32_t rng = 0x2468ACE1u;
    double phase[40] = {0};
    for (int n = 0; n < FRAMES * 160; ++n) {
        double t = (double)n / 8000.0;
        double f0 = 175.0 + 85.0 * sin(2.0 * M_PI * 0.37 * t);
        double tilt = 0.6 + 0.3 * sin(2.0 * M_PI * 0.11 * t);
        double x = 0.0;
        for (int k = 1; k < 40 && k * f0 < 3800.0; ++k) {
            phase[k] += 2.0 * M_PI * k * f0 / 8000.0;
            x += pow(tilt, k) * sin(phase[k]);
        }
        rng = rng * 1664525u + 1013904223u;
        x = 3000.0 * x + 300.0 * (((double)(rng >> 8) / 16777216.0) - 0.5);
        if (x > 32767.0) {
            x = 32767.0;
        } else if (x < -32768.0) {
            x = -32768.0;
        }
        pcm[n] = (short)x;
    }
}

static int
compare(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

/* The three encoders behind one signature. */
struct bench_codec {
    const char* name;
    int bits;
    void* (*alloc)(void);
    void (*reset)(void* enc);
    void (*release)(void* enc);
    int (*encode)(void* enc, const short* frame, char* bits, mbe_parms* cur, const mbe_parms* prev);
};

static void*
dstar_alloc(void) {
    return mbe_ambe2400EncoderAlloc();
}

static void
dstar_reset(void* enc) {
    mbe_ambe2400EncoderReset(enc);
}

static void
dstar_release(void* enc) {
    mbe_ambe2400EncoderFree(enc);
}

static int
dstar_encode(void* enc, const short* frame, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2400ParmsShort(enc, frame, bits, cur, prev);
}

static void*
ambe2450_alloc(void) {
    return mbe_ambe2450EncoderAlloc();
}

static void
ambe2450_reset(void* enc) {
    mbe_ambe2450EncoderReset(enc);
}

static void
ambe2450_release(void* enc) {
    mbe_ambe2450EncoderFree(enc);
}

static int
ambe2450_encode(void* enc, const short* frame, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2450ParmsShort(enc, frame, bits, cur, prev);
}

static void*
imbe_alloc(void) {
    return mbe_imbe4400EncoderAlloc();
}

static void
imbe_reset(void* enc) {
    mbe_imbe4400EncoderReset(enc);
}

static void
imbe_release(void* enc) {
    mbe_imbe4400EncoderFree(enc);
}

static int
imbe_encode(void* enc, const short* frame, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeImbe4400ParmsShort(enc, frame, bits, cur, prev);
}

static int
bench(const struct bench_codec* codec) {
    void* enc = codec->alloc();
    if (!enc) {
        (void)fprintf(stderr, "%s: encoder allocation failed\n", codec->name);
        return 1;
    }
    double per_frame[RUNS];
    unsigned checksum = 0;
    for (int r = 0; r < RUNS; ++r) {
        mbe_parms cur, prev, enhanced;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        codec->reset(enc);
        char bits[88];
        double t0 = secs();
        for (int f = 0; f < FRAMES; ++f) {
            if (codec->encode(enc, pcm + (size_t)f * 160u, bits, &cur, &prev) < 0) {
                (void)fprintf(stderr, "%s: encode failed at frame %d\n", codec->name, f);
                codec->release(enc);
                return 1;
            }
            mbe_moveMbeParms(&cur, &prev);
            for (int i = 0; i < codec->bits; ++i) {
                checksum = (checksum * 31u) + (unsigned)bits[i];
            }
        }
        per_frame[r] = (secs() - t0) / FRAMES * 1e6;
    }
    codec->release(enc);
    qsort(per_frame, RUNS, sizeof(per_frame[0]), compare);
    printf("encode %-8s: median %.2f us/frame (min %.2f, max %.2f) over %d runs of %d frames, checksum %08x\n",
           codec->name, per_frame[RUNS / 2], per_frame[0], per_frame[RUNS - 1], RUNS, FRAMES, checksum);
    return 0;
}

int
main(void) {
    static const struct bench_codec codecs[] = {
        {"ambe2400", 49, dstar_alloc, dstar_reset, dstar_release, dstar_encode},
        {"ambe2450", 49, ambe2450_alloc, ambe2450_reset, ambe2450_release, ambe2450_encode},
        {"imbe4400", 88, imbe_alloc, imbe_reset, imbe_release, imbe_encode},
    };
    make_signal();
    for (size_t i = 0; i < sizeof(codecs) / sizeof(codecs[0]); ++i) {
        if (bench(&codecs[i]) != 0) {
            return 1;
        }
    }
    return 0;
}
