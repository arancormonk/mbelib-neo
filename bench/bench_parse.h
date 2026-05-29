// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_NEO_BENCH_PARSE_H
#define MBELIB_NEO_BENCH_PARSE_H

#include <errno.h>
#include <limits.h>
#include <stdlib.h>

static inline int
bench_parse_int_arg(const char* text, int min_value, int max_value, int* out_value) {
    char* end = NULL;
    long value;

    if (!text || !out_value || min_value > max_value) {
        return -1;
    }

    errno = 0;
    value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value < (long)min_value || value > (long)max_value) {
        return -1;
    }

    *out_value = (int)value;
    return 0;
}

#endif /* MBELIB_NEO_BENCH_PARSE_H */
