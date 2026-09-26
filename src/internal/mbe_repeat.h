// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Frame-repeat helpers shared by the IMBE and AMBE process paths.
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_REPEAT_H
#define MBELIB_NEO_INTERNAL_MBE_REPEAT_H

#include <stdint.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

/**
 * State a repeated frame keeps from the current frame instead of the model it
 * repeats: this frame's error accounting (the error-rate recursion and mute
 * tests advance on every frame), the tone oscillator phases, and the
 * unvoiced-noise generator state, which must continue rather than replay the
 * previous frame's noise.
 */
struct mbe_repeat_carry {
    float errorRate;
    int errorCountTotal;
    int errorCount4;
    int repeatCount;
    float mutingThreshold;
    int swn;
    uint32_t tonePhase;
    float noiseSeed;
    float noiseOverlap[96];
};

static inline void
mbe_repeat_save_carry(struct mbe_repeat_carry* c, const mbe_parms* mp) {
    c->errorRate = mp->errorRate;
    c->errorCountTotal = mp->errorCountTotal;
    c->errorCount4 = mp->errorCount4;
    c->repeatCount = mp->repeatCount;
    c->mutingThreshold = mp->mutingThreshold;
    c->swn = mp->swn;
    c->tonePhase = mp->tonePhase;
    c->noiseSeed = mp->noiseSeed;
    memcpy(c->noiseOverlap, mp->noiseOverlap, sizeof(c->noiseOverlap));
}

static inline void
mbe_repeat_apply_carry(mbe_parms* mp, const struct mbe_repeat_carry* c) {
    mp->errorRate = c->errorRate;
    mp->errorCountTotal = c->errorCountTotal;
    mp->errorCount4 = c->errorCount4;
    mp->repeatCount = c->repeatCount;
    mp->mutingThreshold = c->mutingThreshold;
    mp->swn = c->swn;
    mp->tonePhase = c->tonePhase;
    mp->noiseSeed = c->noiseSeed;
    memcpy(mp->noiseOverlap, c->noiseOverlap, sizeof(mp->noiseOverlap));
}

/**
 * @brief Replace cur_mp's model with model_src for a frame repeat.
 *
 * Everything in struct mbe_repeat_carry stays with cur_mp. Copying the
 * previous parameter set wholesale would freeze the error-rate recursion
 * across a repeat run and replay the previous frame's noise buffer.
 *
 * @param cur_mp    Current frame parameters, rewritten in place.
 * @param model_src Parameter set whose model is repeated (distinct from cur_mp).
 */
static inline void
mbe_repeat_load_model(mbe_parms* cur_mp, const mbe_parms* model_src) {
    struct mbe_repeat_carry c;
    mbe_repeat_save_carry(&c, cur_mp);
    *cur_mp = *model_src;
    mbe_repeat_apply_carry(cur_mp, &c);
}

/** Next consecutive-repeat count, clamped so caller-owned state cannot overflow. */
static inline int
mbe_repeat_next_count(int prev_count) {
    if (prev_count < 0) {
        return 1;
    }
    return (prev_count >= MBE_MAX_FRAME_REPEATS) ? MBE_MAX_FRAME_REPEATS : prev_count + 1;
}

#endif /* MBELIB_NEO_INTERNAL_MBE_REPEAT_H */
