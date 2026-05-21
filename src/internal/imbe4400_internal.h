// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_NEO_INTERNAL_IMBE4400_INTERNAL_H
#define MBELIB_NEO_INTERNAL_IMBE4400_INTERNAL_H

#include "mbelib-neo/mbelib.h"

void mbe_processImbe4400Dataf_withC0(float* aout_buf, const int* errs2, char* err_str, const char imbe_d[88],
                                     mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced, int uvquality,
                                     int c0_errors, int c0_errors_valid, int c4_errors_valid);

#endif /* MBELIB_NEO_INTERNAL_IMBE4400_INTERNAL_H */
