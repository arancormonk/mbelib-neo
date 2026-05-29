// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Micro-benchmark focused on the unvoiced oscillator path.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bench_parse.h"
#include "mbelib-neo/mbelib.h"

/**
 * @brief Return elapsed CPU seconds using the C clock.
 */
static double
secs(void) {
    clock_t c = clock();
    return (double)c / (double)CLOCKS_PER_SEC;
}

/**
 * @brief Initialize an all-unvoiced parameter set for benchmarking.
 * @param cur  Output current parameter set.
 * @param prev Output previous parameter set (copy of current).
 * @param L    Number of bands (clamped to [8,56]).
 */
static void
init_unvoiced(mbe_parms* cur, mbe_parms* prev, int L) {
    mbe_parms prev_enh;
    mbe_initMbeParms(cur, prev, &prev_enh);
    if (L < 8) {
        L = 8;
    }
    if (L > 56) {
        L = 56;
    }
    cur->w0 = 0.11f; // moderate pitch
    cur->L = L;
    for (int l = 1; l <= L; ++l) {
        cur->Vl[l] = 0; // all unvoiced
        cur->Ml[l] = 0.03f + 0.002f * (float)(l & 7);
        cur->PHIl[l] = 0.0f;
        cur->PSIl[l] = 0.0f;
    }
    *prev = *cur;
}

/**
 * @brief Benchmark entry: runs repeated unvoiced synthesis and prints timing.
 */
int
main(int argc, char** argv) {
    int iters = 2000;
    int runs = 8;
    int L = 36;
    if (argc > 1 && bench_parse_int_arg(argv[1], 8, 56, &L) < 0) {
        fprintf(stderr, "usage: %s [L] [iters] [runs]\n", argv[0]);
        return 2;
    }
    if (argc > 2 && bench_parse_int_arg(argv[2], 1, INT_MAX, &iters) < 0) {
        fprintf(stderr, "usage: %s [L] [iters] [runs]\n", argv[0]);
        return 2;
    }
    if (argc > 3 && bench_parse_int_arg(argv[3], 1, INT_MAX, &runs) < 0) {
        fprintf(stderr, "usage: %s [L] [iters] [runs]\n", argv[0]);
        return 2;
    }

    float out[160];
    mbe_parms cur, prev;
    mbe_setThreadRngSeed(0xBEEFu);
    init_unvoiced(&cur, &prev, L);

    double total = 0.0, best = 1e9, worst = 0.0;
    for (int r = 0; r < runs; ++r) {
        double t0 = secs();
        for (int i = 0; i < iters; ++i) {
            // Wiggle w0 a bit to avoid degenerate recurrence
            cur.w0 = (i & 1) ? 0.10f : 0.12f;
            mbe_synthesizeSpeechf(out, &cur, &prev);
            mbe_moveMbeParms(&cur, &prev);
        }
        double dt = secs() - t0;
        total += dt;
        if (dt < best) {
            best = dt;
        }
        if (dt > worst) {
            worst = dt;
        }
        printf("run %d: %.6f s\n", r + 1, dt);
    }
    printf("L=%d frames=%d -> avg: %.6f s, best: %.6f s, worst: %.6f s\n", L, iters * runs, total / runs, best, worst);
    return 0;
}
