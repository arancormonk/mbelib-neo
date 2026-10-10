// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_INTERNAL_ANALYSIS_KERNELS_H
#define MBE_INTERNAL_ANALYSIS_KERNELS_H
/* Private reductions over n contiguous bins (0..16), with no padding required.
 * sums = { dot(re,w), dot(im,w), dot(w,w), dot(re,re)+dot(im,im) }. */
void mbe_analysis_band_sums(const float* re, const float* im, const float* w, int n, float sums[4]);
float mbe_analysis_residual(const float* re, const float* im, const float* w, int n, float ar, float ai);
#endif
