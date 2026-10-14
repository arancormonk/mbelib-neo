// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#ifndef MBELIB_TESTS_ENABLE_ASSERTIONS_H
#define MBELIB_TESTS_ENABLE_ASSERTIONS_H

/* Force-included by MSVC before any test header can define assert(). */
#ifdef NDEBUG
#undef NDEBUG
#endif

#endif
