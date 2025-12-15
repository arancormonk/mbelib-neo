// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Compiler-specific macros for optimization hints and portability.
 */

#ifndef MBELIB_NEO_INTERNAL_MBE_COMPILER_H
#define MBELIB_NEO_INTERNAL_MBE_COMPILER_H

/**
 * @brief Branch prediction hint for likely conditions.
 *
 * Use in hot loops where the condition is almost always true.
 * Example: if (MBE_LIKELY(ptr != NULL))
 */
#if defined(__GNUC__) || defined(__clang__)
#define MBE_LIKELY(x)   __builtin_expect(!!(x), 1)
#define MBE_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define MBE_LIKELY(x)   (x)
#define MBE_UNLIKELY(x) (x)
#endif

/**
 * @brief Thread-local storage specifier.
 *
 * Portable macro for thread-local variables across compilers.
 */
#if defined(_MSC_VER)
#define MBE_THREAD_LOCAL __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
#define MBE_THREAD_LOCAL __thread
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define MBE_THREAD_LOCAL _Thread_local
#else
/* Fallback for non-C11 compilers: most support __thread as an extension */
#define MBE_THREAD_LOCAL __thread
#endif

#endif /* MBELIB_NEO_INTERNAL_MBE_COMPILER_H */
