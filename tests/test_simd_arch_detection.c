// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Compile-time regression checks for SIMD architecture detection.
 */

#if defined(MBE_TEST_SCENARIO_X64)
/*
 * Exercise the normalized x64 path without redefining compiler-provided
 * target macros such as `_M_X64`. SIMD must be inferred from the normalized
 * architecture classification, not hard-coded by the test.
 */
#define MBE_TEST_OVERRIDE_ARCH_X86_64 1
#elif defined(MBE_TEST_SCENARIO_ARM64EC)
#define MBE_TEST_OVERRIDE_ARCH_ARM64EC 1
#elif defined(MBE_TEST_SCENARIO_X86_32_MSVC_SSE2)
/*
 * Exercise 32-bit x86 SIMD inference via the MSVC `_M_IX86_FP` contract while
 * forcing the normalized architecture classification away from the host.
 */
#define MBE_TEST_OVERRIDE_ARCH_X86_32 1
#endif

/* A SIMD build whose CMake configuration built the AVX2 kernel copies: only
 * baseline x86-64 dispatches at run time (mbe_cpu.h). */
#define MBELIB_ENABLE_SIMD  1
#define MBELIB_AVX2_KERNELS 1
#include "mbe_compiler.h"
#include "mbe_cpu.h"

int
main(void) {
#if defined(MBE_TEST_SCENARIO_X64)
#if !defined(MBE_ARCH_X86_64)
#error "x64 scenario must classify as x86_64"
#endif
#if defined(MBE_ARCH_ARM64EC)
#error "x64 scenario must not classify as ARM64EC"
#endif
#if defined(MBE_ARCH_AARCH64)
#error "x64 scenario must not classify as AArch64"
#endif
#if !defined(MBE_SIMD_TARGET_SSE2)
#error "x64 scenario must enable SSE2 intrinsics"
#endif
#if defined(MBE_SIMD_TARGET_NEON)
#error "x64 scenario must not enable NEON intrinsics"
#endif
#if defined(__AVX2__) && defined(MBE_AVX2_DISPATCH)
#error "an x64 target with AVX2 must call its kernels directly"
#endif
#if !defined(__AVX2__) && !defined(MBE_AVX2_DISPATCH)
#error "a baseline x64 target must dispatch the AVX2 kernels at run time"
#endif
#elif defined(MBE_TEST_SCENARIO_ARM64EC)
#if !defined(MBE_ARCH_ARM64EC)
#error "ARM64EC scenario must classify as ARM64EC"
#endif
#if defined(MBE_ARCH_X86_64)
#error "ARM64EC scenario must not classify as x86_64"
#endif
#if !defined(MBE_ARCH_AARCH64)
#error "ARM64EC scenario must classify as AArch64"
#endif
#if !defined(MBE_SIMD_TARGET_NEON)
#error "ARM64EC scenario must enable NEON intrinsics"
#endif
#if defined(MBE_SIMD_TARGET_SSE2)
#error "ARM64EC scenario must not select SSE2 over NEON"
#endif
#if defined(MBE_AVX2_DISPATCH)
#error "ARM64EC scenario must not dispatch AVX2 kernels"
#endif
#elif defined(MBE_TEST_SCENARIO_X86_32_MSVC_SSE2)
#if !defined(MBE_ARCH_X86_32)
#error "x86_32 MSVC SSE2 scenario must classify as x86_32"
#endif
#if defined(MBE_ARCH_X86_64)
#error "x86_32 MSVC SSE2 scenario must not classify as x86_64"
#endif
#if defined(MBE_ARCH_AARCH64)
#error "x86_32 MSVC SSE2 scenario must not classify as AArch64"
#endif
#if defined(MBE_ARCH_ARM64EC)
#error "x86_32 MSVC SSE2 scenario must not classify as ARM64EC"
#endif
#if !defined(MBE_SIMD_TARGET_SSE2)
#error "x86_32 MSVC SSE2 scenario must enable SSE2 intrinsics"
#endif
#if defined(MBE_SIMD_TARGET_NEON)
#error "x86_32 MSVC SSE2 scenario must not enable NEON intrinsics"
#endif
#if defined(MBE_AVX2_DISPATCH)
#error "x86_32 scenario must not dispatch AVX2 kernels"
#endif
#else
#error "test_simd_arch_detection requires an explicit scenario"
#endif

    return 0;
}
