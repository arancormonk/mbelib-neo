// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Exact parity checks for mbe_floattoshort().
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

#define TEST_FRAME_LEN 160

static short
reference_floattoshort(float sample) {
    const float again = 7.0f;
    const float max_amplitude = 32767.0f * 0.95f;
    float audio = again * sample;
    if (audio > max_amplitude) {
        audio = max_amplitude;
    } else if (audio < -max_amplitude) {
        audio = -max_amplitude;
    }
    return (short)audio;
}

static void
fill_test_input(float* input, uint32_t seed) {
    uint32_t state = seed;
    const float clip_point = (32767.0f * 0.95f) / 7.0f;

    for (int i = 0; i < TEST_FRAME_LEN; ++i) {
        state = (state * 1664525u) + 1013904223u;
        int32_t v = (int32_t)(state >> 8) - 0x007FFFFF;
        input[i] = (float)v / 65536.0f;
    }

    input[0] = 0.0f;
    input[1] = clip_point;
    input[2] = clip_point + (1.0f / 32768.0f);
    input[3] = clip_point - (1.0f / 32768.0f);
    input[4] = -clip_point;
    input[5] = -clip_point - (1.0f / 32768.0f);
    input[6] = -clip_point + (1.0f / 32768.0f);
    input[7] = 1.0f / 7.0f;
    input[8] = -1.0f / 7.0f;
}

int
main(void) {
    static const uint32_t seeds[] = {
        0x00000001u,
        0x12345678u,
        0x00C0FFEEu,
        0xFFFFFFFFu,
    };

    float input[TEST_FRAME_LEN];
    short out[TEST_FRAME_LEN];
    short out2[TEST_FRAME_LEN];

    for (size_t s = 0; s < (sizeof(seeds) / sizeof(seeds[0])); ++s) {
        fill_test_input(input, seeds[s]);

        mbe_floattoshort(input, out);
        mbe_floattoshort(input, out2);
        if (memcmp(out, out2, sizeof(out)) != 0) {
            fprintf(stderr, "determinism failure for seed=0x%08X\n", (unsigned)seeds[s]);
            return 1;
        }

        for (int i = 0; i < TEST_FRAME_LEN; ++i) {
            short expected = reference_floattoshort(input[i]);
            if (out[i] != expected) {
                fprintf(stderr, "parity failure seed=0x%08X idx=%d in=%f got=%d exp=%d\n", (unsigned)seeds[s], i,
                        input[i], (int)out[i], (int)expected);
                return 1;
            }
        }
    }

    return 0;
}
