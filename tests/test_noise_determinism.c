// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Sample-exact warm replay and interleaving of nonmuted FFT synthesis.
 *
 * Unvoiced excitation and overlap belong to mbe_parms. Comfort noise instead
 * uses the documented thread-local generator and is intentionally not tested
 * as independently replayable stream state here.
 */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "mbelib-neo/mbelib.h"

typedef struct {
    mbe_parms cur;
    mbe_parms prev;
} synthesis_state;

static int
float_array_equal_exact(const float* a, const float* b, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        uint32_t a_bits, b_bits;
        memcpy(&a_bits, &a[i], sizeof(a_bits));
        memcpy(&b_bits, &b[i], sizeof(b_bits));
        if (a_bits != b_bits) {
            return 0;
        }
    }
    return 1;
}

static void
fill_params(synthesis_state* state, int mixed, int stream) {
    mbe_parms enhanced;
    mbe_initMbeParms(&state->cur, &state->prev, &enhanced);
    state->cur.w0 = stream ? 0.125f : 0.10f;
    state->cur.L = stream ? 20 : 24;
    for (int l = 1; l <= state->cur.L; ++l) {
        state->cur.Vl[l] = mixed && ((l + stream) % 3 == 0);
        state->cur.Ml[l] = 0.04f + 0.001f * (float)(l + 3 * stream);
    }
    state->prev = state->cur;
}

static double
advance(synthesis_state* state, float pcm[160]) {
    mbe_synthesizeSpeechf(pcm, &state->cur, &state->prev);
    double energy = 0;
    for (int i = 0; i < 160; ++i) {
        assert(isfinite(pcm[i]));
        energy += (double)pcm[i] * pcm[i];
    }
    mbe_moveMbeParms(&state->cur, &state->prev);
    return energy;
}

static void
test_warm_streams(int mixed) {
    synthesis_state a, b;
    float discarded[160], expected[12][160], actual[160];
    fill_params(&a, mixed, 0);
    fill_params(&b, !mixed, 1);
    for (int n = 0; n < 3; ++n) {
        advance(&a, discarded);
        advance(&b, discarded);
    }

    synthesis_state replay = a, interleaved = a;
    /* Warm replay must exercise nonzero excitation, not equal muted buffers. */
    for (int n = 0; n < 12; ++n) {
        assert(advance(&a, expected[n]) > 0);
    }
    assert(!float_array_equal_exact(expected[0], expected[1], 160));
    for (int n = 0; n < 12; ++n) {
        advance(&replay, actual);
        assert(float_array_equal_exact(expected[n], actual, 160));
    }
    for (int n = 0; n < 12; ++n) {
        assert(advance(&b, discarded) > 0);
        advance(&interleaved, actual);
        assert(float_array_equal_exact(expected[n], actual, 160));
    }
}

int
main(void) {
    test_warm_streams(0);
    test_warm_streams(1);
    return 0;
}
