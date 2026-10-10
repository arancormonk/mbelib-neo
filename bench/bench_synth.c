// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Micro-benchmark for synthesis hot paths.
 *
 * Measures `mbe_synthesizeSpeechf` with a deterministic parameter stream
 * and reports elapsed CPU time across several repeated runs.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mbelib-neo/mbelib.h"

/**
 * @brief Return elapsed CPU seconds using `clock()`.
 */
static double
secs(void) {
    clock_t c = clock();
    return (double)c / (double)CLOCKS_PER_SEC;
}

/**
 * @brief Benchmark entry: runs repeated synthesis loops and prints timing.
 */
int
main(int argc, char** argv) {
    const char* workload = argc > 1 ? argv[1] : "steady";
    const char* voicing = argc > 2 ? argv[2] : "mixed";
    if (argc > 3
        || (strcmp(workload, "steady") != 0 && strcmp(workload, "gradual") != 0 && strcmp(workload, "jump") != 0)
        || (strcmp(voicing, "mixed") != 0 && strcmp(voicing, "voiced") != 0)) {
        fprintf(stderr, "Usage: %s [steady|gradual|jump] [mixed|voiced]\n", argv[0]);
        return 2;
    }
    const int jump = strcmp(workload, "jump") == 0;
    const int gradual = strcmp(workload, "gradual") == 0;
    const int voiced = strcmp(voicing, "voiced") == 0;
    const int iters = 2000; // frames per run
    const int runs = 10;    // repeats
    float out[160];
    mbe_parms cur, prev, prev_enh;

    mbe_setThreadRngSeed(0x123456u);
    mbe_initMbeParms(&cur, &prev, &prev_enh);

    // Create a moderately voiced/unvoiced mix
    cur.w0 = 0.10f;
    cur.L = (int)(0.9254f * 3.14159265f / cur.w0);
    for (int l = 1; l <= cur.L; ++l) {
        cur.Vl[l] = (l % 3) != 0;       // mix of voiced/unvoiced
        cur.Ml[l] = 0.05f + 0.002f * l; // gentle slope
        cur.log2Ml[l] = 0.0f;
        cur.PHIl[l] = (float)l * 0.1f; // varying phases
        cur.PSIl[l] = (float)l * 0.05f;
    }
    prev = cur; // start with same state

    double total = 0.0, best = 1e9, worst = 0.0;
    const mbe_parms initial = cur;
    printf("synth %s %s\n", workload, voicing);
    /* A complete untimed pass also allocates the thread-local FFT plan. */
    for (int r = -1; r < runs; ++r) {
        cur = initial;
        prev = initial;
        mbe_setThreadRngSeed(0x123456u);
        double t0 = secs();
        for (int i = 0; i < iters; ++i) {
            if (jump) {
                cur.w0 = (i & 1) ? 0.09f : 0.11f;
            } else if (gradual) {
                cur.w0 = 0.10f + 0.01f * sinf((float)i * 0.02f);
            }
            cur.L = (int)(0.9254f * 3.14159265f / cur.w0);
            for (int l = 1; l <= cur.L; ++l) {
                cur.Vl[l] = voiced || ((i + l) % 5) != 0;
                cur.Ml[l] = 0.04f + 0.003f * (float)((i + l) % 7);
            }
            mbe_synthesizeSpeechf(out, &cur, &prev);
            mbe_moveMbeParms(&cur, &prev);
        }
        double dt = secs() - t0;
        if (r < 0) {
            continue;
        }
        total += dt;
        if (dt < best) {
            best = dt;
        }
        if (dt > worst) {
            worst = dt;
        }
        printf("run %d: %.6f s\n", r + 1, dt);
    }
    printf("avg: %.6f s, best: %.6f s, worst: %.6f s (frames=%d)\n", total / runs, best, worst, iters * runs);
    return 0;
}
