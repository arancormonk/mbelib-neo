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

#endif /* MBELIB_NEO_INTERNAL_IMBE4400_INTERNAL_H */
