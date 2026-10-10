// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_INTERNAL_TONE_FIT_H
#define MBE_INTERNAL_TONE_FIT_H
/* Sample rate of the tone detector's input, in Hz. */
#define MBE_TONE_FS 8000.0

/* Least-squares projection onto count (1 or 2) sinusoids at hz. The private
 * caller supplies 1..160 samples at MBE_TONE_FS and distinct frequencies.
 * Returns explained energy, or zero for a singular fit; writes amplitudes on
 * success. */
double mbe_tone_fit(const float* samples, int n, const double hz[2], int count, double amplitude[2]);
#endif
