// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 *
 * Copyright (C) 2010 mbelib Author
 * GPG Key ID: 0xEA5EFE2C (9E7A 5527 9CDC EBF7 BF1B  D772 4F98 E863 EA5E FE2C)
 *
 * Portions were originally under the ISC license; this mbelib-neo
 * distribution is provided under GPL-2.0-or-later. See LICENSE for details.
 */

/** @file D-STAR DV air-interface mapping for the 72-bit AMBE frame. */

#include <string.h>

#include "mbelib-neo/mbelib.h"

/*
 * D-STAR DV frame interleave tables (air bit -> plane/index), from the
 * D-STAR specification. Air bit i (0..71, after the 24-bit sync word)
 * carries ambe_fr[dW[i]][dX[i]].
 */
static const int dstar_dW[72] = {
    0, 0, 3, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, 1, 3, 2, 1, 1, 0, 0, 3, 2, 0, 0, 3, 2, 1, 1, 0, 0, 1, 1, 0, 0,
    3, 2, 1, 1, 3, 2, 1, 1, 0, 0, 3, 2, 0, 0, 3, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, 1, 3, 3, 2, 1, 0, 0, 3, 3,
};

static const int dstar_dX[72] = {
    10, 22, 11, 9, 10, 22, 11, 23, 8, 20, 9, 21, 10, 8, 9, 21, 8, 6,  7,  19, 8, 20, 9, 7,
    6,  18, 7,  5, 6,  18, 7,  19, 4, 16, 5, 17, 6,  4, 5, 17, 4, 2,  3,  15, 4, 16, 5, 3,
    2,  14, 3,  1, 2,  14, 3,  15, 0, 12, 1, 13, 2,  0, 1, 13, 0, 12, 10, 11, 0, 12, 1, 13,
};

/**
 * @brief Serialize a 72-bit AMBE 3600x2400 frame into the 9 data bytes
 *        of a D-STAR DV frame (the 24-bit sync word is NOT included).
 *
 * Bytes are packed in air order (LSB first within each byte, matching
 * the GMSK modulator). This is the exact inverse of
 * mbe_decodeDStarDVData().
 *
 * @param ambe_fr  Input frame as 4x24 bitplanes.
 * @param bytes9   Output 9 bytes (72 bits).
 * @return 0 on success.
 */
int
mbe_encodeDStarDVData(const char ambe_fr[4][24], unsigned char bytes9[9]) {
    if (ambe_fr == NULL || bytes9 == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    memset(bytes9, 0, 9);
    for (int i = 0; i < 72; i++) {
        if (ambe_fr[dstar_dW[i]][dstar_dX[i]] & 1) {
            bytes9[i >> 3] |= (unsigned char)(1u << (i & 7));
        }
    }
    return 0;
}

/**
 * @brief Extract a 72-bit AMBE 3600x2400 frame from the 9 data bytes of
 *        a D-STAR DV frame.
 *
 * @see mbe_encodeDStarDVData for the byte/bit convention.
 *
 * @param bytes9  Input 9 bytes (72 bits).
 * @param ambe_fr Output frame as 4x24 bitplanes.
 * @return 0 on success.
 */
int
mbe_decodeDStarDVData(const unsigned char bytes9[9], char ambe_fr[4][24]) {
    if (bytes9 == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    memset(ambe_fr, 0, 4 * sizeof(ambe_fr[0]));
    for (int i = 0; i < 72; i++) {
        ambe_fr[dstar_dW[i]][dstar_dX[i]] = (char)((bytes9[i >> 3] >> (i & 7)) & 1u);
    }
    return 0;
}
