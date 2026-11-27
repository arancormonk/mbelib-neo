// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief The reserved pitch codes of the pulsed voicing state (private API).
 *
 * DVSI's P25 and D-STAR encoders give a voice frame with no voiced band
 * either one default pitch code (P25 b0 25, D-STAR b0 34) or one of 56 others
 * (P25 b0 128..183, D-STAR b0 64..119). Each of the 56 names, for the frame's
 * own subframe and for its mid-frame subframe, which of eight 500 Hz channels
 * are pulsed: a band there plays no noise and, instead, the response to a
 * single pulse (docs/dvsi-parity.md, "D-STAR and P25: the pulsed state"). The
 * patterns are Table 2 of US 7,970,606, read; that D-STAR's codes carry the
 * same patterns at b0 - 64, and that DVSI's decoders play them, were measured
 * on DVSI's decodes of its own encodings. The synthesizer (mbe_subframe.c) holds the patterns and
 * renders the state; it asks the codec whose rendering it serves which code,
 * if any, a model's pitch is, through the functions below (and, for P25's
 * validity counter, which reads a code only while it stands at 0 or more,
 * which pitch code, if any, a pitch is: US 7,970,606's table lists b0 25 too,
 * and only a model at some code's pitch moves the counter). Each codec answers
 * with its decoder's own pitch function, the one body (not inlined) the
 * decode uses, so that the decode and the test cannot round apart in any
 * build.
 */
#ifndef MBEINT_MBE_PULSED_H
#define MBEINT_MBE_PULSED_H

/* The reserved codes, and the index of a code among them (0..55). */
#define MBE_PULSED_CODES       56
#define MBE_PULSED_P25_FIRST   128
#define MBE_PULSED_DSTAR_FIRST 64
/* P25's default code for a frame with no voiced band, which the validity
 * counter counts as a code of the table (mbe_subframe.c). */
#define MBE_PULSED_P25_DEFAULT 25

/**
 * @brief The reserved P25 pitch code whose pitch @p w0 is, as an index.
 *
 * The pitch P25's decoder stores for b0 128..183 (imbe_code_pitch(), DVSI's
 * pitch: the nominal float times 1 + g, rounded once), compared by bit
 * pattern. ProVoice, which keeps the nominal pitch, has no reserved code.
 * @param w0 A pitch, positive, finite and normal (the caller checks).
 * @return b0 - MBE_PULSED_P25_FIRST (0..55), or -1 for any other pitch.
 */
int mbe_imbe4400_pulsed_code(float w0);

/**
 * @brief The P25 pitch code whose pitch @p w0 is.
 *
 * The pitch P25's decoder stores for b0 0..207 (imbe_code_pitch(), DVSI's
 * pitch), compared by bit pattern. P25's validity counter (mbe_subframe.c)
 * moves only for a model at such a pitch: a code of US 7,970,606's Table 2
 * (b0 25, listed with both subframes unvoiced, and the 56) raises it, any
 * other code lowers it, and a pitch no code decodes to leaves it.
 * @param w0 A pitch, positive, finite and normal (the caller checks).
 * @return b0 (0..207), or -1 for a pitch no code decodes to.
 */
int mbe_imbe4400_pitch_code(float w0);

/**
 * @brief The reserved D-STAR pitch code whose pitch @p w0 is, as an index.
 *
 * The pitch D-STAR's decoder stores for b0 64..119 (the pitch law,
 * mbe_ambe2400_f0(), times 2 pi, as the decode forms it), compared by bit
 * pattern.
 * @param w0 A pitch, positive, finite and normal (the caller checks).
 * @return b0 - MBE_PULSED_DSTAR_FIRST (0..55), or -1 for any other pitch.
 */
int mbe_ambe2400_pulsed_code(float w0);

/**
 * @brief Whether a model is the rate-33 band-grid model of pitch code @p b0.
 *
 * Rate 33's form of the pulsed state lives in V/UV rows 18 to 31, whose
 * frames the decoder plays on the band grid: L = mbe_ambe2450_harmonic_count(b0)
 * bands and w0 = MBE_AMBE2450_GRID_SPAN / L (ambe3600x2450.c). A pulsed-row
 * word (mbe_adaptive.h) counts only on such a model (mbe_subframe.c): the
 * harmonic count equal and the pitch equal by bit pattern to what the decoder's
 * own grid function (ambe2450_grid_model(), the one body, not inlined, the decode
 * stores the model with) gives for @p b0.
 * @param w0 The model's pitch.
 * @param L  The model's harmonic count.
 * @param b0 A pitch code, 0..119 (the caller checks).
 * @return 1 for the grid model of @p b0, 0 otherwise.
 */
int mbe_ambe2450_pulsed_row_model(float w0, int L, int b0);

#endif /* MBEINT_MBE_PULSED_H */
