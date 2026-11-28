// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_NEO_INTERNAL_IMBE4400_INTERNAL_H
#define MBELIB_NEO_INTERNAL_IMBE4400_INTERNAL_H

#include "mbelib-neo/mbelib.h"

int mbe_processImbe4400Dataf_internal(float* aout_buf, mbe_process_result* result, const char imbe_d[88],
                                      mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/*
 * The IMBE encoder's quantization (imbe7200x4400_enc.c), exposed for the
 * tests: b0 (0..207), the band decisions (eq 49) and log2 amplitudes logm
 * indexed by harmonic 1..L of b0's harmonic count, quantized against prev_mp
 * (TIA-102.BABA 6.3) and packed with synchronization bit sync.
 */
void mbe_imbe4400_quantize_amplitudes(int b0, const unsigned char bands[12], const float logm[57],
                                      const mbe_parms* prev_mp, int sync, char imbe_d[88]);

/*
 * The P25 tone table (imbe7200x4400.c, imbe_tone_rows and imbe_tone_first_row),
 * exposed for the tests that pin it. mbe_imbe4400_tone_row() gives row i's
 * pitch code, voicing word (band 1 the most significant bit) and harmonics
 * and returns 1 for a row, 0 past the last (any field pointer may be NULL).
 * mbe_imbe4400_tone_first_row() gives the index entry of pitch code b0: the
 * first row of b0 for 0..197, the number of rows for 198, -1 otherwise.
 */
int mbe_imbe4400_tone_row(int i, int* b0, int* b1, int* l1, int* l2);
int mbe_imbe4400_tone_first_row(int b0);

/*
 * P25's quiet-frame gain (imbe7200x4400.c; docs/dvsi-parity.md, "P25: quiet
 * frames"), exposed for the tests: g = 2^clamp(K (x - X0), FLOOR, 0), x the
 * mean of mp's log2Ml[1..L]. Returns 1 and stores g in *gain where g < 1, and
 * returns 0 and stores 1 otherwise (x at or above X0, L outside 1..56, or a
 * log2Ml in 1..L beyond 2^20 in magnitude, a NaN or an infinity among them). P25's data path multiplies every
 * magnitude it renders by g, after adaptive smoothing (not on a tone model,
 * not in ProVoice).
 */
int mbe_imbe4400_quiet_gain(const mbe_parms* mp, double* gain);

#endif /* MBELIB_NEO_INTERNAL_IMBE4400_INTERNAL_H */
