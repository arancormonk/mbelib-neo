// SPDX-License-Identifier: GPL-2.0-or-later
/* Independent, ordered arithmetic oracles for the optimized TETRA stages. */
#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "tetra_acelp.h"
#include "tetra_acelp_constants.h"

static uint32_t seed = 0x752A1CE3u;

static int16_t
sample(void) {
    seed = seed * 1664525u + 1013904223u;
    int32_t v = (int32_t)(seed >> 16);
    return (int16_t)(v >= 32768 ? v - 65536 : v);
}

static int32_t
sat32(int64_t v) {
    return v > INT32_MAX ? INT32_MAX : (v < INT32_MIN ? INT32_MIN : (int32_t)v);
}

static int16_t
sat16(int64_t v) {
    return v > INT16_MAX ? INT16_MAX : (v < INT16_MIN ? INT16_MIN : (int16_t)v);
}

static int64_t
floor_shift(int64_t v, int bits) {
    int64_t divisor = INT64_C(1) << bits;
    return v / divisor - (v < 0 && v % divisor != 0);
}

static int16_t
wrap16(int64_t v) {
    int32_t low = (int32_t)((uint64_t)v & 65535u);
    return (int16_t)(low >= 32768 ? low - 65536 : low);
}

static void
syn_reference(const int16_t* a, const int16_t* x, int16_t* y, int len, int16_t* mem, int update) {
    int16_t past[TETRA_M + TETRA_FRAME];
    memcpy(past, mem, TETRA_M * sizeof(*mem));
    for (int n = 0; n < len; ++n) {
        int32_t acc = (int32_t)x[n] * 4096;
        for (int i = 1; i <= TETRA_M; ++i) {
            acc = sat32((int64_t)acc - (int64_t)a[i] * past[TETRA_M + n - i]);
        }
        y[n] = (int16_t)floor_shift(sat32((int64_t)sat32((int64_t)acc * 16) + 32768), 16);
        past[TETRA_M + n] = y[n];
    }
    if (update) {
        memcpy(mem, past + len, TETRA_M * sizeof(*mem));
    }
}

static void
test_synthesis(void) {
    static const int lengths[] = {0, 1, 2, 3, 4, 9, 10, 59, 60, 239, 240};
    for (int trial = 0; trial < 256; ++trial) {
        int16_t a_storage[TETRA_M + 2];
        int16_t* a = a_storage + 1;
        int16_t input[TETRA_FRAME + 2], actual[TETRA_FRAME + 2], expect[TETRA_FRAME + 2];
        int16_t mem[TETRA_M], ref_mem[TETRA_M];
        a[0] = 4096;
        for (int i = 1; i <= TETRA_M; ++i) {
            a[i] = (int16_t)(sample() / (trial % 2 ? 8 : 1));
        }
        if (trial < 4) {
            /* Just below/at/above the SIMD bound, including an INT16_MIN tap. */
            memset(a + 1, 0, TETRA_M * sizeof(*a));
            a[1] = INT16_MIN;
            a[2] = (int16_t)(28670 + trial);
        }
        for (int n = 0; n < TETRA_FRAME + 2; ++n) {
            input[n] = trial % 3 ? sample() : INT16_MIN;
        }
        for (size_t l = 0; l < sizeof(lengths) / sizeof(lengths[0]); ++l) {
            int len = lengths[l];
            for (int i = 0; i < TETRA_M; ++i) {
                mem[i] = ref_mem[i] = sample();
            }
            memcpy(actual, input, sizeof(input));
            memcpy(expect, input, sizeof(input));
            syn_reference(a, input + 1, expect + 1, len, ref_mem, trial % 2);
            /* Both separate and aliased input/output, with unaligned buffers and canaries. */
            tetra_syn_filt(a, trial % 3 ? input + 1 : actual + 1, actual + 1, len, mem, trial % 2);
            assert(memcmp(actual, expect, sizeof(actual)) == 0);
            assert(memcmp(mem, ref_mem, sizeof(mem)) == 0);
        }
    }
}

static void
adaptive_reference(int16_t* cur, int t0, int frac) {
    const int16_t* taps = frac > 0 ? tetra_interp_p13 : tetra_interp_m13;
    int first = frac > 0 ? -16 : -15;
    for (int n = 0; n < TETRA_SUBFRAME; ++n) {
        if (frac == 0) {
            cur[n] = cur[n - t0];
        } else {
            int32_t sum = 0;
            for (int k = 0; k < 32; ++k) {
                sum = sat32((int64_t)sum + (int64_t)cur[n - t0 + first + k] * taps[k]);
            }
            cur[n] = (int16_t)floor_shift(sat32((int64_t)sat32((int64_t)sum * 2) + 32768), 16);
        }
    }
}

static void
test_adaptive(void) {
    int bound_p = 0, bound_m = 0;
    for (int k = 0; k < 18; ++k) {
        bound_p += tetra_interp_p13[k] < 0 ? -tetra_interp_p13[k] : tetra_interp_p13[k];
        bound_m += tetra_interp_m13[k] < 0 ? -tetra_interp_m13[k] : tetra_interp_m13[k];
    }
    assert(bound_p <= 65535 && bound_m <= 65535);
    for (int lag = 19; lag <= 144; ++lag) {
        for (int frac = -1; frac <= 1; ++frac) {
            for (int trial = 0; trial < 8; ++trial) {
                int16_t actual[TETRA_EXC_HIST + TETRA_SUBFRAME + 2];
                int16_t expect[TETRA_EXC_HIST + TETRA_SUBFRAME + 2];
                int16_t* cur = actual + 1 + TETRA_EXC_HIST;
                for (size_t i = 0; i < sizeof(actual) / sizeof(actual[0]); ++i) {
                    actual[i] = trial % 2 ? sample() : INT16_MIN;
                }
                if (trial >= 4 && frac != 0) {
                    /* Align signs with the taps to force intermediate saturation, then
                     * reverse the last taps so reassociation would change the result. */
                    const int16_t* taps = frac > 0 ? tetra_interp_p13 : tetra_interp_m13;
                    int first = frac > 0 ? -16 : -15;
                    for (int k = 0; k < 32; ++k) {
                        int positive = (taps[k] > 0) ^ (trial & 1) ^ (k >= 22);
                        cur[-lag + first + k] = positive ? INT16_MAX : INT16_MIN;
                    }
                }
                memcpy(expect, actual, sizeof(actual));
                adaptive_reference(expect + 1 + TETRA_EXC_HIST, lag, frac);
                tetra_adaptive_vector(cur, lag, frac);
                assert(memcmp(actual, expect, sizeof(actual)) == 0);
            }
        }
    }
}

static void
code_reference(const int16_t* a, int lag, const tetra_subframe_params* sp, int16_t* code) {
    int16_t ap1[TETRA_M + 1], ap2[TETRA_M + 1];
    int16_t x[TETRA_SUBFRAME] = {0}, mem[TETRA_M] = {0}, f[TETRA_SUBFRAME];
    ap1[0] = ap2[0] = a[0];
    for (int i = 1; i <= TETRA_M; ++i) {
        ap1[i] = sat16(floor_shift((int64_t)a[i] * tetra_gamma1_pow[i - 1] + 16384, 15));
        ap2[i] = sat16(floor_shift((int64_t)a[i] * tetra_gamma2_pow[i - 1] + 16384, 15));
    }
    memcpy(x, ap1, sizeof(ap1));
    syn_reference(ap2, x, f, TETRA_SUBFRAME, mem, 0);
    for (int n = lag; n < TETRA_SUBFRAME; ++n) {
        f[n] = sat16(f[n] + floor_shift((int64_t)f[n - lag] * 26216, 15));
    }
    int pos[4] = {2 * sp->pulse1 + sp->shift, 2 + 8 * sp->pulse2 + sp->shift, 4 + 8 * sp->pulse3 + sp->shift,
                  6 + 8 * sp->pulse4 + sp->shift};
    static const int weights[4] = {23168, -16384, 16384, -16384};
    for (int n = 0; n < TETRA_SUBFRAME; ++n) {
        int64_t sum = 0;
        for (int k = 0; k < 4; ++k) {
            if (n >= pos[k]) {
                sum += (int64_t)weights[k] * f[n - pos[k]];
            }
        }
        code[n] = wrap16(floor_shift(sp->sign ? -sum : sum, 14));
    }
}

static void
test_codevector(void) {
    for (int trial = 0; trial < 1024; ++trial) {
        int16_t a[TETRA_M + 1], actual[TETRA_SUBFRAME + 2], expect[TETRA_SUBFRAME + 2];
        tetra_subframe_params sp = {0};
        a[0] = 4096;
        for (int i = 1; i <= TETRA_M; ++i) {
            a[i] = (int16_t)(sample() / (trial % 2 ? 8 : 1));
        }
        sp.pulse1 = trial % 32;
        sp.pulse2 = trial / 2 % 8;
        sp.pulse3 = trial / 4 % 8;
        sp.pulse4 = trial / 8 % 8;
        sp.shift = trial / 32 % 2;
        sp.sign = trial / 64 % 2;
        int lag = 19 + trial % 126;
        memset(actual, 0x5A, sizeof(actual));
        memcpy(expect, actual, sizeof(actual));
        code_reference(a, lag, &sp, expect + 1);
        tetra_codevector(a, lag, &sp, actual + 1);
        assert(memcmp(actual, expect, sizeof(actual)) == 0);
    }
}

static void
test_energy_and_excitation(void) {
    for (int trial = 0; trial < 256; ++trial) {
        int16_t v[TETRA_SUBFRAME + 2], c[TETRA_SUBFRAME + 2], expect[TETRA_SUBFRAME + 2];
        for (int n = 0; n < TETRA_SUBFRAME + 2; ++n) {
            v[n] = trial % 3 ? sample() : INT16_MIN;
            if (trial % 4 == 1) {
                v[n] = (int16_t)(v[n] / 32); /* Multi-block energy below the saturation limit. */
            }
            c[n] = trial % 2 ? sample() : INT16_MIN;
        }
        for (int len = 0; len <= TETRA_SUBFRAME; ++len) {
            int64_t sum = trial % 2;
            for (int n = 0; n < len; ++n) {
                sum += (int64_t)v[n + 1] * v[n + 1];
            }
            assert(tetra_energy(v + 1, len, trial % 2) == sat32(sum));
        }
        int16_t gp = trial % 2 ? (int16_t)((uint16_t)sample() % 4916u) : 4915;
        int16_t gc = trial % 2 ? sample() : INT16_MIN;
        memcpy(expect, v, sizeof(v));
        for (int n = 1; n <= TETRA_SUBFRAME; ++n) {
            expect[n] = wrap16(floor_shift((int64_t)gp * v[n] + (int64_t)gc * c[n] + 2048, 12));
        }
        tetra_mix_excitation(gp, gc, v + 1, c + 1);
        assert(memcmp(v, expect, sizeof(v)) == 0);
    }
}

static void
test_interpolation(void) {
    for (int trial = 0; trial < 256; ++trial) {
        int16_t old[TETRA_M + 1], next[TETRA_M + 1], actual[TETRA_SUBFRAMES][TETRA_M];
        for (int i = 0; i < TETRA_M + 1; ++i) {
            old[i] = trial % 2 ? sample() : INT16_MIN;
            next[i] = trial % 3 ? sample() : INT16_MIN;
        }
        tetra_lsp_interpolate(old + 1, next + 1, actual);
        for (int k = 0; k < TETRA_SUBFRAMES; ++k) {
            for (int i = 0; i < TETRA_M; ++i) {
                int64_t sum = (int64_t)(3 - k) * old[i + 1] + (int64_t)(k + 1) * next[i + 1];
                assert(actual[k][i] == floor_shift(sum, 2));
            }
        }
    }
}

int
main(void) {
    test_synthesis();
    test_adaptive();
    test_codevector();
    test_energy_and_excitation();
    test_interpolation();
    return 0;
}
