// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef INDEPENDENT_TEST_COMMON_H
#define INDEPENDENT_TEST_COMMON_H
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "mbelib-neo/mbelib.h"

#define TEST_PI 3.14159265358979323846264338327950288

/* Exact reproducibility and unchanged-state checks compare object bytes. */
static inline int
compare_bytes(const void* left, const void* right, size_t count) {
    const unsigned char* a = left;
    const unsigned char* b = right;
    for (size_t i = 0; i < count; ++i) {
        if (a[i] != b[i]) {
            return 1;
        }
    }
    return 0;
}

static inline uint32_t
float_bits(float value) {
    uint32_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static inline void
assert_close(double actual, double expected, double absolute, double relative) {
    if (!(fabs(actual - expected) <= absolute + relative * fabs(expected))) {
        fprintf(stderr, "actual %.17g, expected %.17g, difference %.17g\n", actual, expected, actual - expected);
        assert(0);
    }
}

static inline void
new_model(mbe_parms* mp) {
    mbe_initMbeParms(mp, mp, mp);
}

static inline double
window_equation(int n) {
    const double distance = fabs((double)n);
    if (distance <= 55.0) {
        return 1.0;
    }
    return distance < 105.0 ? (105.0 - distance) / 50.0 : 0.0;
}
#endif
