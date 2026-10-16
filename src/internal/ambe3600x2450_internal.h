// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief AMBE+2 3600x2450 (rate 33) conventions shared by the decoder and the
 *        encoder (private API).
 */
#ifndef MBELIB_NEO_INTERNAL_AMBE3600X2450_INTERNAL_H
#define MBELIB_NEO_INTERNAL_AMBE3600X2450_INTERNAL_H

#include <math.h>

/*
 * Pitch law of DVSI's AMBE-3000 for rate 33, in its decoder and its encoder:
 * log2 f0 = -4.322605 - 0.0220427 b0 (f0 in cycles per sample) for the voice
 * codes 0..119. TIA-102.BABA-1 Annex A tabulates f0 up to 0.2% higher for some
 * codes (b0 70-76, with oversized steps at 30-31 and 76-77); DVSI's decoded
 * harmonics sit on this law (in-frame fit residual rms 1.2e-4 in log2 f0) and
 * its re-quantizations of a decoded pitch fall on this law's decision points.
 */
#define MBE_AMBE2450_LOG2_F0_OFFSET (-4.322605f)
#define MBE_AMBE2450_LOG2_F0_SLOPE  0.0220427f
#define MBE_AMBE2450_PITCH_CODES    120

/* Fundamental (cycles per sample) of voice code b0. */
static inline float
mbe_ambe2450_f0(int b0) {
    return exp2f(MBE_AMBE2450_LOG2_F0_OFFSET - (MBE_AMBE2450_LOG2_F0_SLOPE * (float)b0));
}

/*
 * Harmonic count of voice code b0: every harmonic at or below 0.9254 of
 * Nyquist, the rule Annex A's L column follows from its own f0. On this law it
 * differs from that column only at b0 17 (12, not 11), where DVSI's output
 * matches 12. The nearest integer is at least 0.004 away for every code.
 */
static inline int
mbe_ambe2450_harmonic_count(int b0) {
    int L = (int)(0.9254f * 0.5f / mbe_ambe2450_f0(b0));
    return (L < 9) ? 9 : ((L > 56) ? 56 : L);
}

/* Voice code nearest to f0 (cycles per sample) in log2, as DVSI's encoder
 * quantizes; clamped to 0..119. */
int mbe_ambe2450_pitch_code(float f0);

/*
 * The 49 parameter bits as b0..b8: unpack reads them as the decoder does and
 * pack writes them as the encoder does, both from one table of bit positions
 * (most significant bit first).
 */
void mbe_ambe2450_unpack(const char ambe_d[49], int b[9]);
void mbe_ambe2450_pack(const int b[9], char ambe_d[49]);

#endif /* MBELIB_NEO_INTERNAL_AMBE3600X2450_INTERNAL_H */
