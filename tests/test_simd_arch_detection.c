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
#elif defined(MBE_TEST_SCENARIO_ARM32_NEON)
/*
 * Exercise 32-bit ARM whose compiler target enables NEON. `__ARM_NEON` is
 * supplied by the build the same way `_M_IX86_FP` is above, so this drives the
 * real (non-override) SIMD inference path rather than a synthetic shortcut.
 */
#define MBE_TEST_OVERRIDE_ARCH_ARM_32 1
#elif defined(MBE_TEST_SCENARIO_ARM32_NO_NEON)
/*
 * The same target without NEON. This is the configuration that silently
 * compiled the NEON paths out before MBELIB_ENABLE_SIMD learned to request
 * -mfpu=neon on 32-bit ARM; the scenario keeps that fallback explicit.
 */
#define MBE_TEST_OVERRIDE_ARCH_ARM_32 1
#endif

#include "mbe_compiler.h"

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
#elif defined(MBE_TEST_SCENARIO_ARM32_NEON)
#if !defined(MBE_ARCH_ARM_32)
#error "ARM32 NEON scenario must classify as 32-bit ARM"
#endif
#if defined(MBE_ARCH_AARCH64)
#error "ARM32 NEON scenario must not classify as AArch64"
#endif
#if defined(MBE_ARCH_X86_64) || defined(MBE_ARCH_X86_32)
#error "ARM32 NEON scenario must not classify as x86"
#endif
#if !defined(MBE_SIMD_TARGET_NEON)
#error "ARM32 NEON scenario must enable NEON intrinsics"
#endif
#if defined(MBE_SIMD_TARGET_SSE2)
#error "ARM32 NEON scenario must not enable SSE2 intrinsics"
#endif
#elif defined(MBE_TEST_SCENARIO_ARM32_NO_NEON)
#if !defined(MBE_ARCH_ARM_32)
#error "ARM32 scalar scenario must classify as 32-bit ARM"
#endif
#if defined(MBE_SIMD_TARGET_NEON)
#error "ARM32 scalar scenario must not enable NEON intrinsics"
#endif
#if defined(MBE_SIMD_TARGET_SSE2)
#error "ARM32 scalar scenario must not enable SSE2 intrinsics"
#endif
#else
#error "test_simd_arch_detection requires an explicit scenario"
#endif

    return 0;
}
