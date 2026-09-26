// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_VALIDATION_H
#define MBELIB_NEO_INTERNAL_MBE_VALIDATION_H

#define MBE_MIN_HARMONIC_BANDS 1
#define MBE_MAX_HARMONIC_BANDS 56
#define MBE_MAX_FRAME_BITS     184

static inline int
mbe_harmonic_count_is_valid(int L) {
    return L >= MBE_MIN_HARMONIC_BANDS && L <= MBE_MAX_HARMONIC_BANDS;
}

static inline int
mbe_clamp_harmonic_count(int L) {
    if (L < MBE_MIN_HARMONIC_BANDS) {
        return MBE_MIN_HARMONIC_BANDS;
    }
    if (L > MBE_MAX_HARMONIC_BANDS) {
        return MBE_MAX_HARMONIC_BANDS;
    }
    return L;
}

static inline int
mbe_error_count_is_valid(int count) {
    return count >= 0 && count <= MBE_MAX_FRAME_BITS;
}

/** Nonzero when the caller-owned cur/prev/prev_enhanced parameter sets are non-null, distinct objects. */
static inline int
mbe_parms_triplet_is_valid(const void* cur_mp, const void* prev_mp, const void* prev_mp_enhanced) {
    return cur_mp && prev_mp && prev_mp_enhanced && (cur_mp != prev_mp) && (cur_mp != prev_mp_enhanced)
           && (prev_mp != prev_mp_enhanced);
}

#endif /* MBELIB_NEO_INTERNAL_MBE_VALIDATION_H */
