// SPDX-License-Identifier: GPL-2.0-or-later
/** @file AVX2/FMA build of mbe_tone_fit.c for the run-time dispatch (mbe_cpu.h).
 *
 * CMake compiles this unit with AVX2/FMA code generation and defines
 * MBELIB_AVX2_KERNELS only for configurations that dispatch; elsewhere it is
 * empty. */
#if defined(MBELIB_AVX2_KERNELS)
#if !defined(__AVX2__)
#error "mbe_tone_fit_avx2.c must be compiled with AVX2 code generation"
#endif
#define MBE_KERNEL_AVX2 1
#include "mbe_tone_fit.c" // NOLINT(bugprone-suspicious-include): the same kernels, built for AVX2
#else
/* ISO C forbids an empty translation unit. */
typedef int mbe_tone_fit_avx2_unused;
#endif
