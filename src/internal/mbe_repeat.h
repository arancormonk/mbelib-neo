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

#include <string.h>

#include "mbelib-neo/mbelib.h"

/**
 * Per-frame error accounting. It advances on every frame, whatever the frame
 * type, so the error-rate recursion and the mute tests see each frame.
 */
static inline void
mbe_copy_error_state(mbe_parms* dst, const mbe_parms* src) {
    dst->errorRate = src->errorRate;
    dst->errorCountTotal = src->errorCountTotal;
    dst->errorCount4 = src->errorCount4;
    dst->repeatCount = src->repeatCount;
    dst->mutingThreshold = src->mutingThreshold;
}

/**
 * Tone oscillator phases and unvoiced-noise generator state, which must
 * continue rather than replay the previous frame's.
 */
static inline void
mbe_copy_generator_state(mbe_parms* dst, const mbe_parms* src) {
    dst->swn = src->swn;
    dst->tonePhase = src->tonePhase;
    dst->noiseSeed = src->noiseSeed;
    memcpy(dst->noiseOverlap, src->noiseOverlap, sizeof(dst->noiseOverlap));
}

/**
 * @brief Replace cur_mp's model with model_src for a frame repeat.
 *
 * cur_mp keeps its error state and generator state. Copying the previous
 * parameter set wholesale would freeze the error-rate recursion across a
 * repeat run and replay the previous frame's noise buffer.
 *
 * @param cur_mp    Current frame parameters, rewritten in place.
 * @param model_src Parameter set whose model is repeated.
 */
static inline void
mbe_repeat_load_model(mbe_parms* cur_mp, const mbe_parms* model_src) {
    mbe_parms model = *model_src;
    mbe_copy_error_state(&model, cur_mp);
    mbe_copy_generator_state(&model, cur_mp);
    *cur_mp = model;
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
