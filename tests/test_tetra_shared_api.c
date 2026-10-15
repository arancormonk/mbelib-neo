// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief The TETRA entry points, called through the shared library on every platform.
 *
 * test_api links the static library on Windows; this test makes sure the TETRA functions are
 * exported from the DLL as well.
 */

#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

int
main(void) {
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    uint8_t bits[MBE_TETRA_FRAME_BITS];
    int16_t first[MBE_TETRA_FRAME_SAMPLES];
    int16_t again[MBE_TETRA_FRAME_SAMPLES];
    int i;
    assert(dec);
    for (i = 0; i < MBE_TETRA_FRAME_BITS; i++) {
        bits[i] = (uint8_t)((i * 5 + 1) % 3 == 0);
    }
    assert(mbe_tetraDecodeFrame(dec, bits, 0, first) == 0);
    assert(mbe_tetraDecodeFrame(dec, NULL, 1, again) == 0);
    bits[3] = 2;
    assert(mbe_tetraDecodeFrame(dec, bits, 0, again) == MBE_STATUS_INVALID_BITS);
    bits[3] = 0;
    mbe_tetraDecoderReset(dec);
    assert(mbe_tetraDecodeFrame(dec, bits, 0, again) == 0);
    assert(memcmp(first, again, sizeof(first)) == 0);
    mbe_tetraDecoderFree(dec);
    return 0;
}
