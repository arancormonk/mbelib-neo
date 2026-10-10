// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_INTERNAL_VOICED_H
#define MBE_INTERNAL_VOICED_H

/* Private 160-sample voiced kernels. They do not modify phase/RNG history. */
struct mbe_voiced_component {
    float phase, omega, gain;
    int active;
};

void mbe_voiced_windowed(float out[160], const float window[320], const struct mbe_voiced_component* prev,
                         const struct mbe_voiced_component* cur);
void mbe_voiced_interpolated(float out[160], float phase, float prev_amp, float cur_amp, float frequency,
                             float pitch_delta, int harmonic);
#endif
