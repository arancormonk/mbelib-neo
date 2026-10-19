// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief CPU probe and kernel-ISA selection for the run-time AVX2 dispatch (see mbe_cpu.h).
 *
 * Built with the library's baseline flags only.
 */

#include "mbe_cpu.h"
#include "mbe_simd.h" // IWYU pragma: keep (MBE_VF_WIDTH gives the compile-time kernel ISA)

#if defined(MBE_ARCH_X86_64) && !defined(MBE_ARCH_ARM64EC)                                                             \
    && (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__))
#define MBE_CPU_X86_PROBE 1
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif

static void
cpu_id(unsigned leaf, unsigned subleaf, unsigned r[4]) {
    unsigned a = 0, b = 0, c = 0, d = 0;
#if defined(_MSC_VER)
    int info[4];
    __cpuidex(info, (int)leaf, (int)subleaf);
    a = (unsigned)info[0];
    b = (unsigned)info[1];
    c = (unsigned)info[2];
    d = (unsigned)info[3];
#else
    __cpuid_count(leaf, subleaf, a, b, c, d);
#endif
    r[0] = a;
    r[1] = b;
    r[2] = c;
    r[3] = d;
}

/* XCR0. GCC and Clang (clang-cl too) expose XGETBV only to code built for
 * XSAVE, so this wrapper alone is; it runs only after CPUID reports OSXSAVE,
 * and GCC and Clang do not inline it into the baseline callers. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((target("xsave"))) static unsigned long long
cpu_xgetbv0(void) {
    return __builtin_ia32_xgetbv(0);
}
#else
static unsigned long long
cpu_xgetbv0(void) {
    return _xgetbv(0);
}
#endif

static int
cpu_bits(unsigned reg, unsigned mask) {
    return (reg & mask) == mask;
}

/* The x86-64-v3 set: compilers given /arch:AVX2 (clang-cl) or a v3 target
 * may emit any of it, not only AVX2 and FMA. */
static int
cpu_probe_avx2(void) {
    unsigned r[4];
    cpu_id(0u, 0u, r);
    const unsigned max_leaf = r[0];
    if (max_leaf < 7u) {
        return 0;
    }
    cpu_id(1u, 0u, r);
    /* SSE3, SSSE3, FMA, CX16, SSE4.1, SSE4.2, MOVBE, POPCNT, OSXSAVE, AVX, F16C */
    const unsigned leaf1_ecx = (1u << 0) | (1u << 9) | (1u << 12) | (1u << 13) | (1u << 19) | (1u << 20) | (1u << 22)
                               | (1u << 23) | (1u << 27) | (1u << 28) | (1u << 29);
    if (!cpu_bits(r[2], leaf1_ecx)) {
        return 0;
    }
    /* The OS saves XMM and YMM state. */
    if ((cpu_xgetbv0() & 0x6u) != 0x6u) {
        return 0;
    }
    cpu_id(7u, 0u, r);
    /* BMI1, AVX2, BMI2 */
    if (!cpu_bits(r[1], (1u << 3) | (1u << 5) | (1u << 8))) {
        return 0;
    }
    cpu_id(0x80000000u, 0u, r);
    if (r[0] < 0x80000001u) {
        return 0;
    }
    cpu_id(0x80000001u, 0u, r);
    /* LAHF/SAHF, LZCNT */
    return cpu_bits(r[2], (1u << 0) | (1u << 5));
}
#endif

int
mbe_cpu_has_avx2(void) {
#if defined(MBE_CPU_X86_PROBE)
    return cpu_probe_avx2();
#else
    return 0;
#endif
}

int
mbe_cpu_dispatch_enabled(void) {
#if defined(MBE_AVX2_DISPATCH)
    return 1;
#else
    return 0;
#endif
}

#if defined(MBE_AVX2_DISPATCH)
int mbe_cpu_active_isa = MBE_CPU_ISA_AUTO;

MBE_NOINLINE int
mbe_cpu_resolve(void) {
    const int isa = cpu_probe_avx2() ? MBE_CPU_ISA_AVX2 : MBE_CPU_ISA_SSE2;
    MBE_STORE_RELAXED(&mbe_cpu_active_isa, isa);
    return isa;
}

int
mbe_cpu_isa(void) {
    const int isa = MBE_LOAD_RELAXED(&mbe_cpu_active_isa);
    return isa != MBE_CPU_ISA_AUTO ? isa : mbe_cpu_resolve();
}

int
mbe_cpu_select_isa(int isa) {
    if (isa == MBE_CPU_ISA_AUTO) {
        (void)mbe_cpu_resolve();
        return 0;
    }
    if (isa == MBE_CPU_ISA_SSE2 || (isa == MBE_CPU_ISA_AVX2 && cpu_probe_avx2())) {
        MBE_STORE_RELAXED(&mbe_cpu_active_isa, isa);
        return 0;
    }
    return -1;
}
#else
/* The width mbe_simd.h selects for this unit, which shares the library's flags. */
static int
compiled_isa(void) {
#if defined(MBE_VF_WIDTH) && defined(MBE_SIMD_TARGET_SSE2)
    return MBE_VF_WIDTH == 8 ? MBE_CPU_ISA_AVX2 : MBE_CPU_ISA_SSE2;
#elif defined(MBE_VF_WIDTH) && defined(MBE_SIMD_TARGET_NEON)
    return MBE_CPU_ISA_NEON;
#else
    return MBE_CPU_ISA_SCALAR;
#endif
}

int
mbe_cpu_isa(void) {
    return compiled_isa();
}

int
mbe_cpu_select_isa(int isa) {
    return (isa == MBE_CPU_ISA_AUTO || isa == compiled_isa()) ? 0 : -1;
}
#endif
