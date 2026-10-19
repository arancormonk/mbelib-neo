// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Private run-time selection of the AVX2/FMA kernels on x86-64.
 *
 * The kernels whose output depends on the vector width (voiced oscillators,
 * tone projection), and the exact integer kernels that gain from AVX2 (soft
 * ECC search, TETRA adaptive codebook), are built twice when the library
 * targets baseline x86-64: once for the target and once for AVX2/FMA. A
 * baseline caller then picks one per call from a cached CPU probe.
 *
 * `MBE_AVX2_DISPATCH` is defined in the baseline translation units of such a
 * build: SIMD x86-64 (not ARM64EC), CMake built the AVX2 copies
 * (`MBELIB_AVX2_KERNELS`), and the target itself lacks AVX2. A target with
 * AVX2 (`-march=native` on an AVX2 host, `-march=x86-64-v3`, `/arch:AVX2`)
 * calls its kernels directly, with no check. The AVX2 copies define
 * `MBE_KERNEL_AVX2` before including the baseline source.
 *
 * `MBE_KFN(name)` names a kernel in the source that is built twice: `name_avx2`
 * in the AVX2 copy, `name_sse2` in a dispatching baseline unit, else `name`.
 * The dispatcher, defined in the baseline unit under `MBE_AVX2_DISPATCH`,
 * keeps the plain name, so callers do not change.
 *
 * An AVX2 translation unit defines only static functions and `*_avx2`
 * externals, so no inline or shared helper compiled for AVX2 can be picked by
 * a baseline caller (test_avx2_symbols checks this on Linux).
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_CPU_H
#define MBELIB_NEO_INTERNAL_MBE_CPU_H

#include "mbe_compiler.h"

/** Kernel ISA: the vector width the width-dependent kernels run at. */
#define MBE_CPU_ISA_AUTO   0 /* mbe_cpu_select_isa(): back to the CPU probe */
#define MBE_CPU_ISA_SCALAR 1
#define MBE_CPU_ISA_SSE2   2
#define MBE_CPU_ISA_AVX2   3
#define MBE_CPU_ISA_NEON   4

#if defined(MBELIB_ENABLE_SIMD) && defined(MBELIB_AVX2_KERNELS) && defined(MBE_SIMD_TARGET_SSE2)                       \
    && defined(MBE_ARCH_X86_64) && !defined(MBE_ARCH_ARM64EC) && !defined(MBE_ARCH_X86_32) && !defined(__AVX2__)       \
    && !defined(MBE_KERNEL_AVX2)
#define MBE_AVX2_DISPATCH 1
#endif

#if defined(MBE_KERNEL_AVX2)
#define MBE_KFN(name) name##_avx2
#elif defined(MBE_AVX2_DISPATCH)
#define MBE_KFN(name) name##_sse2
#else
#define MBE_KFN(name) name
#endif

/** The ISA the dispatched kernels use now; a build without dispatch reports its compile-time ISA. */
int mbe_cpu_isa(void);
/** Nonzero when this x86-64 CPU and OS support the whole x86-64-v3 set (AVX2, FMA, BMI1/2, LZCNT, MOVBE, F16C). */
int mbe_cpu_has_avx2(void);
/** Nonzero when the library was built with run-time AVX2 dispatch. */
int mbe_cpu_dispatch_enabled(void);
/** Test-only: force a kernel ISA (MBE_CPU_ISA_AUTO restores the probe). Returns 0, or -1 when this build or CPU cannot
 * run it. Not thread-safe against concurrent decoding. */
int mbe_cpu_select_isa(int isa);

#if defined(MBE_AVX2_DISPATCH)
/* MBE_CPU_ISA_AUTO until the first dispatched call resolves it. */
extern int mbe_cpu_active_isa;
/** Probe the CPU and publish the result; idempotent, so racing first calls agree. */
int mbe_cpu_resolve(void);

static inline int
mbe_cpu_use_avx2(void) {
    int isa = MBE_LOAD_RELAXED(&mbe_cpu_active_isa);
    if (MBE_UNLIKELY(isa == MBE_CPU_ISA_AUTO)) {
        isa = mbe_cpu_resolve();
    }
    return isa == MBE_CPU_ISA_AVX2;
}
#endif

#endif /* MBELIB_NEO_INTERNAL_MBE_CPU_H */
