// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Micro-benchmark for the AMBE 3600x2400 (D-STAR) encoder.
 *
 * Encodes a deterministic speech-like signal (a harmonic series whose pitch
 * glides between 90 and 260 Hz, with a slowly varying spectral tilt and some
 * noise) and reports CPU time per 20 ms frame across repeated runs.
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

int
main(void) {
    make_signal();
    mbe_ambe2400_encoder* enc = mbe_ambe2400EncoderAlloc();
    if (!enc) {
        (void)fprintf(stderr, "encoder allocation failed\n");
        return 1;
    }
    double per_frame[RUNS];
    unsigned checksum = 0;
    for (int r = 0; r < RUNS; ++r) {
        mbe_parms cur, prev, enhanced;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        mbe_ambe2400EncoderReset(enc);
        char bits[49];
        double t0 = secs();
        for (int f = 0; f < FRAMES; ++f) {
            if (mbe_encodeAmbe2400ParmsShort(enc, pcm + (size_t)f * 160u, bits, &cur, &prev) < 0) {
                (void)fprintf(stderr, "encode failed at frame %d\n", f);
                mbe_ambe2400EncoderFree(enc);
                return 1;
            }
            mbe_moveMbeParms(&cur, &prev);
            for (int i = 0; i < 49; ++i) {
                checksum = (checksum * 31u) + (unsigned)bits[i];
            }
        }
        per_frame[r] = (secs() - t0) / FRAMES * 1e6;
    }
    mbe_ambe2400EncoderFree(enc);
    qsort(per_frame, RUNS, sizeof(per_frame[0]), compare);
    printf("encode: median %.2f us/frame (min %.2f, max %.2f) over %d runs of %d frames, checksum %08x\n",
           per_frame[RUNS / 2], per_frame[0], per_frame[RUNS - 1], RUNS, FRAMES, checksum);
    return 0;
}
