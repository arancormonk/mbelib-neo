// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_INTERNAL_VOICED_H
#define MBE_INTERNAL_VOICED_H

#include "mbe_simd.h" // IWYU pragma: keep (selects MBE_VOICED_KERNELS)

/* Samples per synthesized frame. The synthesis window Ws spans two frames. */
#define MBE_VOICED_FRAME 160

/* Private one-frame vector voiced kernels, only on targets with mbe_vf
 * vectors; scalar targets keep the historical loops in mbelib.c. They do not
 * modify phase/RNG history. */
#if defined(MBE_VF_WIDTH)
#define MBE_VOICED_KERNELS 1

struct mbe_voiced_component {
    float phase, omega, gain;
    int active;
};

void mbe_voiced_windowed(float out[MBE_VOICED_FRAME], const float window[2 * MBE_VOICED_FRAME],
                         const struct mbe_voiced_component* prev, const struct mbe_voiced_component* cur);
void mbe_voiced_interpolated(float out[MBE_VOICED_FRAME], float phase, float prev_amp, float cur_amp, float frequency,
                             float pitch_delta, int harmonic);
#endif
#endif
