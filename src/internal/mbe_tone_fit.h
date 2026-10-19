// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_INTERNAL_TONE_FIT_H
#define MBE_INTERNAL_TONE_FIT_H
#include "mbe_cpu.h" // IWYU pragma: keep (MBE_AVX2_DISPATCH declares the dispatched copies)
/* Sample rate of the tone detector's input, in Hz. */
#define MBE_TONE_FS 8000.0

/* Least-squares projection onto count (1 or 2) sinusoids at hz. The private
 * caller supplies 1..160 samples at MBE_TONE_FS and distinct frequencies.
 * Returns explained energy, or zero for a singular fit; writes amplitudes on
 * success. */
double mbe_tone_fit(const float* samples, int n, const double hz[2], int count, double amplitude[2]);

/* The baseline and AVX2 copies behind the run-time dispatch (mbe_cpu.h). */
#if defined(MBE_AVX2_DISPATCH) || defined(MBE_KERNEL_AVX2)
double mbe_tone_fit_sse2(const float* samples, int n, const double hz[2], int count, double amplitude[2]);
double mbe_tone_fit_avx2(const float* samples, int n, const double hz[2], int count, double amplitude[2]);
#endif
#endif
