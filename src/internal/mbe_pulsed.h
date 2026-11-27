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
 * on DVSI's decodes. The synthesizer (mbe_subframe.c) holds the patterns and
 * renders the state; it asks the codec whose rendering it serves which code,
 * if any, a model's pitch is, through the functions below. Each codec answers
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
 * @brief The reserved D-STAR pitch code whose pitch @p w0 is, as an index.
 *
 * The pitch D-STAR's decoder stores for b0 64..119 (the pitch law,
 * mbe_ambe2400_f0(), times 2 pi, as the decode forms it), compared by bit
 * pattern.
 * @param w0 A pitch, positive, finite and normal (the caller checks).
 * @return b0 - MBE_PULSED_DSTAR_FIRST (0..55), or -1 for any other pitch.
 */
int mbe_ambe2400_pulsed_code(float w0);

#endif /* MBEINT_MBE_PULSED_H */
