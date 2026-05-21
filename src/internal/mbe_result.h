// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_RESULT_H
#define MBELIB_NEO_INTERNAL_MBE_RESULT_H

#include <stddef.h>
#include "mbelib-neo/mbelib.h"

#define MBE_RESULT_STATUS_SIZE       384u
#define MBE_RESULT_STATUS_SUFFIX_MAX 2u

static inline int
mbe_result_clamp_status_errors(int errors, size_t status_size) {
    if (errors <= 0) {
        return 0;
    }
    if (status_size <= (MBE_RESULT_STATUS_SUFFIX_MAX + 1u)) {
        return 0;
    }

    const size_t max_errors = status_size - MBE_RESULT_STATUS_SUFFIX_MAX - 1u;
    if ((size_t)errors > max_errors) {
        return (int)max_errors;
    }
    return errors;
}

void mbe_result_apply_status(mbe_process_result* result, const char* status);

#endif /* MBELIB_NEO_INTERNAL_MBE_RESULT_H */
