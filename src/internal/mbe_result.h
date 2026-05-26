// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_RESULT_H
#define MBELIB_NEO_INTERNAL_MBE_RESULT_H

#include <stddef.h>
#include "mbelib-neo/mbelib.h"

#define MBE_RESULT_CONTEXT_FLAGS (MBE_PROCESS_FLAG_SOFT_INPUT | MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_C4_VALID)

static inline int
mbe_validate_bits(const char* bits, size_t count) {
    if (!bits) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (size_t i = 0u; i < count; ++i) {
        if (bits[i] != 0 && bits[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }
    return 0;
}

static inline int
mbe_validate_soft_bits(const mbe_soft_bit* bits, size_t count) {
    if (!bits) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (size_t i = 0u; i < count; ++i) {
        if (bits[i].bit > 1u) {
            return MBE_STATUS_INVALID_BITS;
        }
    }
    return 0;
}

static inline int
mbe_result_resolve_total_errors(const mbe_process_result* result, int* total_errors) {
    int total = 0;

    if (!total_errors) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    if (result) {
        if (result->c0_errors < 0 || result->protected_errors < 0 || result->c4_errors < 0
            || result->total_errors < 0) {
            return MBE_STATUS_INVALID_ARGUMENT;
        }
        total = result->total_errors;
        if (total == 0 && (result->c0_errors != 0 || result->protected_errors != 0)) {
            total = result->c0_errors + result->protected_errors;
        }
        if (((result->flags & MBE_PROCESS_FLAG_C0_VALID) != 0u) && total < result->c0_errors) {
            return MBE_STATUS_INVALID_ARGUMENT;
        }
    }

    *total_errors = total;
    return 0;
}

static inline void
mbe_result_prepare_synthesis(mbe_process_result* result, int total_errors) {
    if (!result) {
        return;
    }

    const unsigned context_flags = result->flags & MBE_RESULT_CONTEXT_FLAGS;
    const int c0_errors = ((context_flags & MBE_PROCESS_FLAG_C0_VALID) != 0u) ? result->c0_errors : 0;
    const int c4_errors = ((context_flags & MBE_PROCESS_FLAG_C4_VALID) != 0u) ? result->c4_errors : 0;

    result->flags = context_flags;
    result->c0_errors = c0_errors;
    result->c4_errors = c4_errors;
    result->total_errors = total_errors;
    result->protected_errors = total_errors - c0_errors;
}

static inline void
mbe_result_set_flag(mbe_process_result* result, unsigned flag) {
    if (result) {
        result->flags |= flag;
    }
}

#endif /* MBELIB_NEO_INTERNAL_MBE_RESULT_H */
