// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#include <assert.h>

#ifdef NDEBUG
#error "Test assertions must remain enabled in every build configuration"
#endif

int
main(void) {
    int checks = 0;
    // NOLINTNEXTLINE(bugprone-assert-side-effect): prove assert() evaluates its operand.
    assert(++checks == 1);
    return checks != 1;
}
