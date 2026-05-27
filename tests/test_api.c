// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Basic API smoke test for version reporting.
 */
#include <assert.h>
#include <string.h>
#include "mbelib-neo/mbelib.h"
#include "mbelib-neo/version.h"

/**
 * @brief Test entry: verifies version helpers match MBELIB_VERSION.
 */
int
main(void) {
    assert(strcmp(mbe_versionString(), MBELIB_VERSION) == 0);

    mbe_process_result result;
    mbe_initProcessResult(&result);
    result.total_errors = 3;
    result.flags = MBE_PROCESS_FLAG_REPEAT | MBE_PROCESS_FLAG_MUTE;
    char status[16] = {0};
    mbe_formatProcessResult(status, sizeof(status), &result);
    assert(strcmp(status, "===RM") == 0);

    mbe_soft_bit hard = mbe_softBitFromHard(2, 77);
    assert(hard.bit == 1);
    assert(hard.reliability == 77);
    mbe_soft_bit llr = mbe_softBitFromLlr(-123);
    assert(llr.bit == 0);
    assert(llr.reliability == 123);
    return 0;
}
