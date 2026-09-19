// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief Simple additive AMBE 2400 synthesizer.
 *
 * Renders decoded AMBE 2400 parameters (w0, L, Vl, Ml) directly into 160
 * float samples per 20 ms frame. Voiced bands are phase-continuous
 * sinusoids; unvoiced bands are noise shaped by their amplitude. A slow
 * output AGC normalizes the level to a comfortable listening level.
 *
 * This is a compact, self-contained reference synthesizer intended for
 * applications that want speech directly from the decoded parameters
 * without the full parametric-synthesis pipeline.
 */

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

#define MBE_SYNTH_SAMPLES    160
#define MBE_SYNTH_TARGET_RMS 0.10f

/**
 * @brief Reset synthesizer state (call before the first frame).
 */
void
mbe_synthInit(mbe_synth_state* st) {
    memset(st, 0, sizeof(*st));
    st->rng = 0x12345u;
    st->gain = 0.1f;
}

/**
 * @brief Synthesize one 20 ms frame of AMBE 2400 speech from decoded
 *        parameters into 160 float samples.
 */
void
mbe_synthFrame(mbe_synth_state* st, const mbe_parms* cur, float* out160) {
    float w0 = cur->w0;
    int   L = cur->L;
    float Ml[57];
    float Vl[57];
    float s;

    if (L < 1) {
        L = 1;
    }
    if (L > 56) {
        L = 56;
    }

    /* Interpolate amplitudes toward the current frame to avoid pops. */
    for (int l = 1; l <= 56; l++) {
        Ml[l] = st->inited ? (0.5f * (cur->Ml[l] + st->prev_Ml[l])) : cur->Ml[l];
        Vl[l] = (float)cur->Vl[l];
    }

    for (int n = 0; n < MBE_SYNTH_SAMPLES; n++) {
        s = 0.0f;
        for (int l = 1; l <= L; l++) {
            if (Vl[l] > 0.5f) {
                s += Ml[l] * sinf(st->phase[l] + (float)l * w0 * (float)n);
            } else {
                /* unvoiced band: fast noise scaled by the amplitude */
                st->rng = st->rng * 1103515245u + 12345u;
                float r = (float)((st->rng >> 16) & 0x7fffu) / 16384.0f - 1.0f;
                s += Ml[l] * r;
            }
        }
        out160[n] = s;
    }

    for (int l = 1; l <= 56; l++) {
        st->phase[l] = fmodf(st->phase[l] + (float)l * w0 * (float)MBE_SYNTH_SAMPLES, (float)(2.0 * M_PI));
    }
    memcpy(st->prev_Ml, cur->Ml, sizeof(st->prev_Ml));
    memcpy(st->prev_Vl, cur->Vl, sizeof(st->prev_Vl));
    st->prev_L = L;
    st->inited = 1;

    /* Slow peak-based normalization to a comfortable level (robust to the
     * frame's harmonic count and crest factor). */
    {
        float pk = 0.0f;
        for (int n = 0; n < MBE_SYNTH_SAMPLES; n++) {
            float a = fabsf(out160[n]);
            if (a > pk) {
                pk = a;
            }
        }
        if (pk > 1e-4f) {
            float want = 0.75f / pk;
            st->gain = (0.5f * st->gain) + (0.5f * want);
            if (st->gain > 100.0f) {
                st->gain = 100.0f;
            }
            for (int n = 0; n < MBE_SYNTH_SAMPLES; n++) {
                out160[n] *= st->gain;
                if (out160[n] > 0.85f) out160[n] = 0.85f;
                else if (out160[n] < -0.85f) out160[n] = -0.85f;
            }
        }
    }
}
