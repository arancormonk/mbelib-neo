// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief DVSI's low band: the fixed high-pass and the adaptive low cut of the
 * D-STAR, P25 and rate-33 renderings (docs/dvsi-parity.md, "DVSI's low band").
 *
 * DVSI's decoders, measured on DVSI's test vectors, play the low band in two
 * parts. A fixed high-pass, a third-order Butterworth at
 * MBE_LOW_BAND_DVSI_CORNER_HZ, in magnitude and in phase. And an adaptive low
 * cut: a tracker that starts each stream at MBE_LOW_CUT_START_HZ moves toward
 * the fundamental of every frame the synthesis renders with every band voiced
 * (T <- T + MBE_LOW_CUT_RATE (f0 - T)), but for a P25 tone model and a silence
 * frame (mbe_subframe.c, low_cut_moves()), and stands still through every
 * other frame; the corner a frame applies is min(T, MBE_LOW_CUT_CEILING_HZ); a voiced
 * harmonic at f under the corner is scaled by (f / corner)^2 (the power, 2, is
 * a fit too), its phase untouched, one at or above it not at all; and the bins
 * of the unvoiced transform whose centre frequency lies under the corner are
 * zero. What DVSI's decoders compute is not known: the filter and the tracker
 * are the simplest forms found to fit. The corner of the filter, the tracker's
 * first value, its rate, the ceiling and the cut's power are fits to DVSI's
 * decodes of the development and validation vectors and of its synthetic
 * pulse-train vectors, not values from a document; the tracker's range
 * (MBE_LOW_CUT_TRACK_MIN_HZ, MBE_LOW_CUT_TRACK_MAX_HZ: the clamp, and with the
 * ceiling the range a stored tracker or corner is read in) is the library's
 * choice, no fit.
 *
 * ProVoice (no DVSI material) keeps the second-order filter at 100 Hz
 * (mbe_low_band_gain(), mbe_adaptive.h) and no cut; so does the public
 * TIA-102.BABA rendering of mbe_synthesizeSpeechf().
 */

#ifndef MBEINT_MBE_LOW_BAND_H
#define MBEINT_MBE_LOW_BAND_H

/* The fixed high-pass's corner and order: a third-order Butterworth at 75 Hz.
 * Taken from the phase DVSI's harmonics carry under 300 Hz against ours, on
 * frames the cut does not reach; how finely the corner is known is in
 * docs/dvsi-parity.md ("DVSI's low band", the fixed filter). */
#define MBE_LOW_BAND_DVSI_CORNER_HZ 75.0

/* The tracker's value at the start of a stream, in Hz. Fitted on DVSI's
 * pulse-train vector ptrain5 (no frame fully voiced, so the corner never moves
 * from it): 105 and 115 Hz fit worse. */
#define MBE_LOW_CUT_START_HZ        110.0f

/* The tracker's step toward a fully voiced frame's fundamental, once a frame.
 * 0.042 to 0.044 fit alike, and 1/24; 0.040 and 0.046 worse. */
#define MBE_LOW_CUT_RATE            0.043

/* The ceiling of the corner (not of the tracker, which may stand above it):
 * 196 to 202 Hz fit alike. */
#define MBE_LOW_CUT_CEILING_HZ      200.0f

/* The range a stored tracker is read in, Hz: every pitch code of the three
 * codecs decodes to a fundamental inside it, the pulsed rows' band grid
 * included (test_low_cut_range() in test_subframe_synthesis.c), so a tracker that starts at
 * MBE_LOW_CUT_START_HZ and moves only toward such fundamentals stays inside it.
 * A caller's model at a pitch outside it moves the tracker to the bound at
 * most, and a stored value outside it (or any value no frame stores) reads as
 * a fresh stream's. Not a fit: bounds with room around the codecs' pitches. */
#define MBE_LOW_CUT_TRACK_MIN_HZ    50.0f
#define MBE_LOW_CUT_TRACK_MAX_HZ    500.0f

/**
 * @brief Gain of DVSI's fixed high-pass at a frequency: x^3 / sqrt(1 + x^6),
 *        x = hz / MBE_LOW_BAND_DVSI_CORNER_HZ, formed in double and rounded to
 *        float once, for a positive normal finite hz; by the bit pattern, 0 for
 *        zero, a subnormal, a negative value or a NaN, and 1 for +infinity.
 */
float mbe_low_band_dvsi_gain(float hz);

/**
 * @brief Phase (radians) of the same high-pass: the argument of
 *        s^3 / ((s + 1)(s^2 + s + 1)) at s = j x, 3 pi / 2 toward DC, 3 pi / 4
 *        at the corner, about 2 / x far above (formed there so that it does
 *        not cancel); in double with the double-precision functions and
 *        rounded to float once, for a positive normal finite hz; by the bit
 *        pattern, 0 for zero, a subnormal, a negative value, a NaN and
 *        +infinity.
 */
float mbe_low_band_dvsi_phase(float hz);

/**
 * @brief The cut's scale on a voiced harmonic at @p hz under a corner:
 *        (hz / corner)^2 when hz < corner, formed in double and rounded to float
 *        once, and exactly 1 when hz is at or above it. @p corner is a corner
 *        the renderer read (50 to 200 Hz); @p hz the harmonic's frequency as
 *        the renderer forms it (one float, also the filter's argument), for a
 *        harmonic under Nyquist of a positive normal pitch, so a positive
 *        normal value of at most about 4 kHz whatever a caller's state holds.
 */
float mbe_low_cut_gain(float hz, float corner);

#endif /* MBEINT_MBE_LOW_BAND_H */
