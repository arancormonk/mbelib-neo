// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Run-time AVX2 dispatch: CPU probe, ISA selection and agreement of the kernel copies.
 *
 * The voiced oscillators run at a different vector width per ISA, so their output differs in the last bits but
 * stays within tolerance; the soft ECC search and the TETRA adaptive codebook are exact integer kernels and must
 * agree bit for bit.
 */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbe_compiler.h"
#include "mbe_cpu.h"
#include "mbe_tone_fit.h"
#include "mbe_voiced.h" // IWYU pragma: keep (MBE_VOICED_KERNELS selects the voiced checks)
#include "mbelib-neo/mbelib.h"
#include "tetra_acelp.h"

static uint32_t rng_state = 0x2c9e4f17u;

static uint32_t
next_random(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

/* Probe against the compiler's own CPU-feature runtime where it has one. */
static void
test_probe(void) {
    const int has = mbe_cpu_has_avx2();
    assert(has == 0 || has == 1);
#if defined(MBE_ARCH_X86_64) && !defined(MBE_ARCH_ARM64EC) && defined(__GNUC__)
    __builtin_cpu_init();
    const int builtin = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma") && __builtin_cpu_supports("bmi")
                        && __builtin_cpu_supports("bmi2") && __builtin_cpu_supports("avx")
                        && __builtin_cpu_supports("popcnt") && __builtin_cpu_supports("sse4.2");
    /* Every CPU with AVX2, FMA and BMI2 also has the rest of x86-64-v3 (LZCNT, MOVBE, F16C). */
    printf("probe: AVX2 set %d, compiler runtime %d\n", has, builtin);
    assert(has == builtin);
#elif !defined(MBE_ARCH_X86_64) || defined(MBE_ARCH_ARM64EC)
    assert(has == 0);
#else
    printf("probe: AVX2 set %d\n", has);
#endif
}

static void
test_selection(void) {
    assert(mbe_cpu_select_isa(MBE_CPU_ISA_AUTO) == 0);
    const int automatic = mbe_cpu_isa();
    assert(mbe_cpu_select_isa(-1) == -1);
    assert(mbe_cpu_select_isa(99) == -1);
    if (mbe_cpu_dispatch_enabled()) {
        assert(automatic == (mbe_cpu_has_avx2() ? MBE_CPU_ISA_AVX2 : MBE_CPU_ISA_SSE2));
        assert(mbe_cpu_select_isa(MBE_CPU_ISA_SSE2) == 0 && mbe_cpu_isa() == MBE_CPU_ISA_SSE2);
        assert(mbe_cpu_select_isa(MBE_CPU_ISA_AVX2) == (mbe_cpu_has_avx2() ? 0 : -1));
        assert(mbe_cpu_select_isa(MBE_CPU_ISA_SCALAR) == -1);
        assert(mbe_cpu_select_isa(MBE_CPU_ISA_NEON) == -1);
    } else {
        /* Without dispatch the build runs its compile-time ISA and accepts only that. */
        const int others[] = {MBE_CPU_ISA_SCALAR, MBE_CPU_ISA_SSE2, MBE_CPU_ISA_AVX2, MBE_CPU_ISA_NEON};
        for (unsigned i = 0; i < sizeof(others) / sizeof(others[0]); ++i) {
            assert(mbe_cpu_select_isa(others[i]) == (others[i] == automatic ? 0 : -1));
            assert(mbe_cpu_isa() == automatic);
        }
    }
    assert(mbe_cpu_select_isa(MBE_CPU_ISA_AUTO) == 0 && mbe_cpu_isa() == automatic);
}

/* Both ISAs selectable: the dispatch can be compared against itself. */
static int
both_isas(void) {
    return mbe_cpu_dispatch_enabled() && mbe_cpu_has_avx2();
}

#if defined(MBE_VOICED_KERNELS)
struct voiced_run {
    float windowed[MBE_VOICED_FRAME], interpolated[MBE_VOICED_FRAME];
};

static void
voiced_run(int isa, int trial, struct voiced_run* r) {
    float window[2 * MBE_VOICED_FRAME];
    for (int n = 0; n < 2 * MBE_VOICED_FRAME; ++n) {
        window[n] = n < MBE_VOICED_FRAME ? (float)n / 160.0f : (float)(320 - n) / 160.0f;
    }
    struct mbe_voiced_component prev = {(float)trial * 1.37f, 0.03f + 0.02f * (float)trial, 0.9f, 1};
    struct mbe_voiced_component cur = {(float)trial * -0.71f, 0.05f + 0.019f * (float)trial, 1.3f, 1};
    memset(r, 0, sizeof(*r));
    assert(mbe_cpu_select_isa(isa) == 0);
    mbe_voiced_windowed(r->windowed, window, &prev, &cur);
    mbe_voiced_interpolated(r->interpolated, (float)trial * 2.3f, 0.4f, 1.1f, 0.02f * (float)(trial + 1), 0.004f,
                            trial % 7 + 1);
}

static void
test_voiced_across_isas(void) {
    int differs = 0;
    for (int trial = 0; trial < 64; ++trial) {
        struct voiced_run sse2, avx2;
        voiced_run(MBE_CPU_ISA_SSE2, trial, &sse2);
        voiced_run(MBE_CPU_ISA_AVX2, trial, &avx2);
        for (int n = 0; n < MBE_VOICED_FRAME; ++n) {
            assert(fabsf(sse2.windowed[n] - avx2.windowed[n]) <= 1e-4f);
            assert(fabsf(sse2.interpolated[n] - avx2.interpolated[n]) <= 5e-4f);
        }
        differs |= memcmp(&sse2, &avx2, sizeof(sse2)) != 0;
    }
    /* Equal bits for every trial would mean the AVX2 copy never ran. */
    assert(differs);
}
#endif

static void
test_tone_fit_across_isas(void) {
    const double pairs[][2] = {{350, 440}, {697, 1209}, {941, 1633}};
    for (unsigned p = 0; p < sizeof(pairs) / sizeof(pairs[0]); ++p) {
        float samples[160];
        for (int i = 0; i < 160; ++i) {
            samples[i] = (float)(800.0 * sin(2.0 * 3.14159265358979 * pairs[p][0] * i / 8000.0)
                                 + 300.0 * cos(2.0 * 3.14159265358979 * pairs[p][1] * i / 8000.0));
        }
        for (int count = 1; count <= 2; ++count) {
            for (int n = 100; n <= 160; n += 7) {
                double amp_sse2[2] = {0}, amp_avx2[2] = {0};
                assert(mbe_cpu_select_isa(MBE_CPU_ISA_SSE2) == 0);
                const double e_sse2 = mbe_tone_fit(samples, n, pairs[p], count, amp_sse2);
                assert(mbe_cpu_select_isa(MBE_CPU_ISA_AVX2) == 0);
                const double e_avx2 = mbe_tone_fit(samples, n, pairs[p], count, amp_avx2);
                assert(fabs(e_sse2 - e_avx2) <= 1e-9 * fmax(1.0, fabs(e_sse2)));
                for (int k = 0; k < count; ++k) {
                    assert(fabs(amp_sse2[k] - amp_avx2[k]) <= 1e-9 * fmax(1.0, amp_sse2[k]));
                }
            }
        }
    }
}

typedef int (*soft_decoder)(const mbe_soft_bit* in, char* out);

static void
test_soft_ecc_across_isas(void) {
    const soft_decoder decoders[] = {mbe_golay2312Soft, mbe_hamming1511Soft, mbe_7100x4400hamming1511Soft};
    const int widths[] = {23, 15, 15};
    for (unsigned d = 0; d < sizeof(decoders) / sizeof(decoders[0]); ++d) {
        for (int trial = 0; trial < 4000; ++trial) {
            mbe_soft_bit in[23];
            for (int i = 0; i < widths[d]; ++i) {
                const uint32_t r = next_random();
                in[i].bit = (uint8_t)(r & 1u);
                /* Many equal weights, to exercise the tie-breaks, and the full range. */
                in[i].reliability = (uint8_t)((trial & 1) ? (r >> 8) & 0xffu : ((r >> 8) & 3u) * 64u);
            }
            char out_sse2[23] = {0}, out_avx2[23] = {0};
            assert(mbe_cpu_select_isa(MBE_CPU_ISA_SSE2) == 0);
            const int r_sse2 = decoders[d](in, out_sse2);
            assert(mbe_cpu_select_isa(MBE_CPU_ISA_AVX2) == 0);
            const int r_avx2 = decoders[d](in, out_avx2);
            assert(r_sse2 == r_avx2 && memcmp(out_sse2, out_avx2, (size_t)widths[d]) == 0);
        }
    }
}

static void
test_tetra_adaptive_across_isas(void) {
    for (int lag = TETRA_PIT_MIN - 1; lag <= TETRA_PIT_MAX + 1; ++lag) {
        for (int frac = -1; frac <= 1; ++frac) {
            for (int trial = 0; trial < 6; ++trial) {
                int16_t history[TETRA_EXC_HIST + TETRA_SUBFRAME];
                for (size_t i = 0; i < sizeof(history) / sizeof(history[0]); ++i) {
                    const uint32_t r = next_random();
                    history[i] = trial < 2 ? (int16_t)((r & 1u) ? INT16_MAX : INT16_MIN) : (int16_t)(r >> 16);
                }
                int16_t sse2[TETRA_EXC_HIST + TETRA_SUBFRAME], avx2[TETRA_EXC_HIST + TETRA_SUBFRAME];
                memcpy(sse2, history, sizeof(history));
                memcpy(avx2, history, sizeof(history));
                assert(mbe_cpu_select_isa(MBE_CPU_ISA_SSE2) == 0);
                tetra_adaptive_vector(sse2 + TETRA_EXC_HIST, lag, frac);
                assert(mbe_cpu_select_isa(MBE_CPU_ISA_AVX2) == 0);
                tetra_adaptive_vector(avx2 + TETRA_EXC_HIST, lag, frac);
                assert(memcmp(sse2, avx2, sizeof(sse2)) == 0);
            }
        }
    }
}

int
main(void) {
    test_probe();
    test_selection();
    printf("dispatch %s, AVX2 set %s, kernel ISA %d\n", mbe_cpu_dispatch_enabled() ? "on" : "off",
           mbe_cpu_has_avx2() ? "yes" : "no", mbe_cpu_isa());
    if (both_isas()) {
#if defined(MBE_VOICED_KERNELS)
        test_voiced_across_isas();
#endif
        test_tone_fit_across_isas();
        test_soft_ecc_across_isas();
        test_tetra_adaptive_across_isas();
        assert(mbe_cpu_select_isa(MBE_CPU_ISA_AUTO) == 0);
        puts("SSE2 and AVX2 kernels agree");
    } else {
        puts("one kernel ISA in this build or on this CPU: cross-ISA checks not applicable");
    }
    return 0;
}
