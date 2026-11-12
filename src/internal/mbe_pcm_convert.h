// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief 16-bit PCM output of the AMBE decode paths.
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_PCM_CONVERT_H
#define MBELIB_NEO_INTERNAL_MBE_PCM_CONVERT_H

/**
 * @brief Convert 160 float samples to 16-bit PCM, rounding to nearest.
 *
 * The `short` entry points of the D-STAR (AMBE 3600x2400) and AMBE 3600x2450
 * paths convert with this. It applies the same 7x scale, the same limit
 * (+-32767 * 0.95) and the same NaN handling (NaN becomes 0) as
 * mbe_floattoshort(), then rounds to the nearest integer, halves away from
 * zero, where mbe_floattoshort() truncates toward zero. DVSI's D-STAR and
 * AMBE 3600x2450 decoders round to nearest (the share of samples identical to
 * DVSI's decodes of its test vectors rises with it); its P25 decoder does not,
 * so the IMBE paths keep truncating. DVSI's vectors separate rounding to
 * nearest from truncation, not one tie rule from another: halves away from
 * zero is this library's choice. Every compiled variant (scalar, SSE2, AVX2,
 * NEON) gives the same output for every input. Does nothing when either
 * pointer is NULL.
 *
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 16-bit samples.
 */
void mbe_floattoshort_nearest(const float* float_buf, short* aout_buf);

#endif /* MBELIB_NEO_INTERNAL_MBE_PCM_CONVERT_H */
