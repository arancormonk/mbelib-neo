// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief The AMBE paths' 16-bit conversion, mbe_floattoshort_nearest(),
 *        against an independent roundf() reference.
 *
 * The library compiles one variant of the conversion (scalar, SSE2, AVX2 or
 * NEON, chosen at compile time), and this test checks whichever one it holds:
 * exact halves of the scaled value and the floats next to each half on both
 * sides, both signs, small and large magnitudes, the limit, values beyond it,
 * infinities, NaN, signed zero and subnormals, a strided sweep over every float
 * bit pattern, each value in every lane position, at unaligned addresses, with
 * canaries around the frame. mbe_floattoshort() must still truncate on the
 * same inputs.
 */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbe_pcm_convert.h"
#include "mbelib-neo/mbelib.h"

#define FRAME     160
#define LIMIT     (32767.0f * 0.95f)
#define MAX_INPUT 600000
#define CANARY    12345

static float inputs[MAX_INPUT];
static size_t input_count;

/* NaN, read from the bit pattern: a float compared with itself is folded
 * away where the compiler assumes finite math. */
static int
is_nan_bits(float x) {
    uint32_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return (bits & UINT32_C(0x7FFFFFFF)) > UINT32_C(0x7F800000);
}

/* The scaled, limited sample both conversions start from. */
static float
reference_limited(float sample) {
    float value = sample * 7.0f;
    if (is_nan_bits(value)) {
        return 0.0f;
    }
    if (value > LIMIT) {
        return LIMIT;
    }
    if (value < -LIMIT) {
        return -LIMIT;
    }
    return value;
}

/* roundf() rounds halves away from zero whatever the rounding mode. */
static short
reference_nearest(float sample) {
    return (short)roundf(reference_limited(sample));
}

static short
reference_truncate(float sample) {
    return (short)reference_limited(sample);
}

static void
add_input(float value) {
    assert(input_count < MAX_INPUT);
    inputs[input_count] = value;
    ++input_count;
}

static void
add_bits(uint32_t bits) {
    float value;
    memcpy(&value, &bits, sizeof(value));
    add_input(value);
}

/* An input whose scaled value is exactly `target`, searched among the floats
 * next to target / 7. Returns 0 when no such input exists nearby. */
static int
input_scaling_to(float target, float* found) {
    float x = target / 7.0f;
    for (int step = 0; step < 16; ++step) {
        x = nextafterf(x, -INFINITY);
    }
    for (int step = 0; step < 33; ++step) {
        const float scaled = x * 7.0f;
        uint32_t a, b;
        memcpy(&a, &scaled, sizeof(a));
        memcpy(&b, &target, sizeof(b));
        if (a == b) {
            *found = x;
            return 1;
        }
        x = nextafterf(x, INFINITY);
    }
    return 0;
}

struct half_coverage {
    size_t exact, below, above;
};

/* The inputs around the half h (both signs): the floats within four steps of
 * h / 7 and, where they exist, inputs scaling exactly to h and to the floats
 * just below and just above it. */
static void
add_half(float h, struct half_coverage* coverage) {
    for (int sign = 0; sign < 2; ++sign) {
        const float half = sign ? -h : h;
        const float targets[3] = {nextafterf(half, -INFINITY), half, nextafterf(half, INFINITY)};
        float x = half / 7.0f;
        for (int step = 0; step < 4; ++step) {
            x = nextafterf(x, -INFINITY);
        }
        for (int step = 0; step < 9; ++step) {
            add_input(x);
            x = nextafterf(x, INFINITY);
        }
        for (int t = 0; t < 3; ++t) {
            float found;
            if (input_scaling_to(targets[t], &found)) {
                add_input(found);
                coverage->below += (t == 0);
                coverage->exact += (t == 1);
                coverage->above += (t == 2);
            }
        }
    }
}

static void
build_inputs(void) {
    struct half_coverage coverage = {0, 0, 0};
    input_count = 0;
    /* Every half up to 1024.5, then halves spread up to the limit. */
    for (int k = 0; k <= 1024; ++k) {
        add_half((float)k + 0.5f, &coverage);
    }
    for (int k = 1025; k < 31140; k += 37) {
        add_half((float)k + 0.5f, &coverage);
    }
    for (int k = 31100; k < 31140; ++k) {
        add_half((float)k + 0.5f, &coverage);
    }
    printf("halves: %zu exact, %zu just below, %zu just above\n", coverage.exact, coverage.below, coverage.above);
    assert(coverage.exact > 1000 && coverage.below > 1000 && coverage.above > 1000);

    /* Around the limit, beyond it and the special values. */
    float x = LIMIT / 7.0f;
    for (int step = 0; step < 8; ++step) {
        x = nextafterf(x, -INFINITY);
    }
    for (int step = 0; step < 17; ++step) {
        add_input(x);
        add_input(-x);
        x = nextafterf(x, INFINITY);
    }
    static const float specials[] = {0.0f, 1.0f, 1000.0f, 4447.0f, 4448.0f, 1e6f, 1e30f, FLT_MAX, FLT_MIN, INFINITY};
    for (size_t i = 0; i < sizeof(specials) / sizeof(specials[0]); ++i) {
        add_input(specials[i]);
        add_input(-specials[i]);
    }
    static const uint32_t special_bits[] = {0x00000001u, 0x80000001u, 0x007fffffu, 0x807fffffu, 0x7fc00000u,
                                            0xffc00001u, 0x7f800001u, 0xff800001u, 0x7fffffffu, 0xffffffffu};
    for (size_t i = 0; i < sizeof(special_bits) / sizeof(special_bits[0]); ++i) {
        add_bits(special_bits[i]);
    }
    /* Every 16411th float bit pattern, finite or not. */
    for (uint64_t bits = 0; bits <= UINT32_MAX; bits += 16411u) {
        add_bits((uint32_t)bits);
    }
}

static int
check_frame(const float* in, const short* nearest, const short* truncated, int offset, size_t base) {
    for (int i = 0; i < FRAME; ++i) {
        const short want_nearest = reference_nearest(in[i]);
        const short want_truncated = reference_truncate(in[i]);
        if (nearest[i] != want_nearest || truncated[i] != want_truncated) {
            uint32_t bits;
            memcpy(&bits, &in[i], sizeof(bits));
            fprintf(stderr,
                    "offset %d frame at %zu lane %d input 0x%08x: nearest %d (want %d), truncated %d (want %d)\n",
                    offset, base, i, (unsigned)bits, nearest[i], want_nearest, truncated[i], want_truncated);
            return 1;
        }
    }
    return 0;
}

static void
fill_canaries(short* buffer, int count) {
    for (int i = 0; i < count; ++i) {
        buffer[i] = CANARY;
    }
}

static int
canaries_intact(const short* out) {
    return out[-1] == CANARY && out[FRAME] == CANARY;
}

/* Every input in every lane position: the frames start `offset` inputs early,
 * so input j sits in lane (j + offset) of its frame, and the buffers sit at
 * odd addresses. */
static int
check_all_lanes(void) {
    for (int offset = 0; offset < 8; ++offset) {
        float storage[FRAME + 8];
        short near_storage[FRAME + 16];
        short trunc_storage[FRAME + 16];
        float* in = storage + offset;
        short* nearest = near_storage + offset + 1;
        short* truncated = trunc_storage + offset + 1;
        for (size_t base = 0; base < input_count + (size_t)offset; base += FRAME) {
            for (int i = 0; i < FRAME; ++i) {
                const size_t index = base + (size_t)i;
                in[i] = (index >= (size_t)offset && index - (size_t)offset < input_count)
                            ? inputs[index - (size_t)offset]
                            : 0.0f;
            }
            fill_canaries(near_storage, FRAME + 16);
            fill_canaries(trunc_storage, FRAME + 16);
            mbe_floattoshort_nearest(in, nearest);
            mbe_floattoshort(in, truncated);
            if (!canaries_intact(nearest) || !canaries_intact(truncated)) {
                fprintf(stderr, "conversion wrote outside the frame (offset %d)\n", offset);
                return 1;
            }
            if (check_frame(in, nearest, truncated, offset, base)) {
                return 1;
            }
        }
    }
    return 0;
}

/* One sample through both conversions, in lane `lane` of a frame of zeros. */
static void
convert_one(float sample, int lane, short* nearest, short* truncated) {
    float in[FRAME];
    short out_nearest[FRAME];
    short out_truncated[FRAME];
    for (int i = 0; i < FRAME; ++i) {
        in[i] = 0.0f;
    }
    in[lane] = sample;
    mbe_floattoshort_nearest(in, out_nearest);
    mbe_floattoshort(in, out_truncated);
    *nearest = out_nearest[lane];
    *truncated = out_truncated[lane];
}

/* Named cases, each in every lane of an 8-sample block. */
static void
check_named_cases(void) {
    static const struct {
        float scaled;
        short nearest, truncated;
    } cases[] = {
        {0.5f, 1, 0},
        {-0.5f, -1, 0},
        {1.5f, 2, 1},
        {-1.5f, -2, -1},
        {2.5f, 3, 2}, /* away from zero, not to even */
        {-2.5f, -3, -2},
        {0.49999997f, 0, 0},
        {-0.49999997f, 0, 0},
        {0.50000006f, 1, 0},
        {1.4999999f, 1, 1},
        {1000.5f, 1001, 1000},
        {-1000.5f, -1001, -1000},
        {31127.5f, 31128, 31127},
    };

    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
        float x;
        const int found = input_scaling_to(cases[c].scaled, &x);
        assert(found);
        for (int lane = 0; lane < 8; ++lane) {
            short nearest, truncated;
            convert_one(x, 152 + lane, &nearest, &truncated);
            assert(nearest == cases[c].nearest);
            assert(truncated == cases[c].truncated);
        }
        (void)found;
    }
    /* The limit 31128.65 rounds to 31129 and truncates to 31128. */
    static const float beyond[] = {INFINITY, 1e30f, 4448.0f, LIMIT / 7.0f};
    for (size_t i = 0; i < sizeof(beyond) / sizeof(beyond[0]); ++i) {
        for (int lane = 0; lane < 8; ++lane) {
            short nearest, truncated;
            convert_one(beyond[i], lane, &nearest, &truncated);
            assert(nearest == 31129 && truncated == 31128);
            convert_one(-beyond[i], lane, &nearest, &truncated);
            assert(nearest == -31129 && truncated == -31128);
        }
    }
    for (int lane = 0; lane < 8; ++lane) {
        short nearest, truncated;
        convert_one(NAN, lane, &nearest, &truncated);
        assert(nearest == 0 && truncated == 0);
        convert_one(-NAN, lane, &nearest, &truncated);
        assert(nearest == 0 && truncated == 0);
    }
    /* NULL buffers are ignored. */
    const float in[FRAME] = {0.0f};
    short out[FRAME];
    mbe_floattoshort_nearest(NULL, out);
    mbe_floattoshort_nearest(in, NULL);
}

int
main(void) {
    check_named_cases();
    build_inputs();
    printf("%zu inputs\n", input_count);
    if (check_all_lanes()) {
        return 1;
    }
    puts("mbe_floattoshort_nearest(): equals roundf() of the limited 7x value; mbe_floattoshort() truncates");
    return 0;
}
