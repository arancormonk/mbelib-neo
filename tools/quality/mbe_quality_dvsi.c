// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Interleave schedules below are from DSD (include/p25p1_const.h and
 * include/dmr_const.h):
 *
 * Copyright (C) 2010 DSD Author
 * GPG Key ID: 0x3F1D7FD0 (74EF 430D F7F2 0A48 FCE6  F630 FAA2 635D 3F1D 7FD0)
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND ISC DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS.  IN NO EVENT SHALL ISC BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 *
 * The rest of this file is distributed under GPL-2.0-or-later.
 */
/**
 * @file
 * @brief Import DVSI hard-decision test-vector frames as rectangular frames.
 *
 * The bits are only deinterleaved: channel errors and FEC parity are kept, so
 * error-variant vectors still exercise the decoder's error handling.
 */
#include "mbe_quality_dvsi.h"

#include <string.h>

#include "mbelib-neo/mbelib.h"

/* Interleaved bit 2i lands at [row_high[i]][col_high[i]], bit 2i+1 at
 * [row_low[i]][col_low[i]]. */
static const unsigned char p25_row_high[72] = {
    0, 2, 4, 1, 3, 5, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6,
    0, 2, 5, 1, 3, 6, 0, 2, 5, 1, 3, 6, 0, 2, 5, 1, 3, 7, 0, 2, 5, 1, 3, 7, 0, 2, 5, 1, 4, 7, 0, 3, 5, 2, 4, 7,
};
static const unsigned char p25_col_high[72] = {
    22, 20, 10, 20, 18, 0, 20, 18, 8, 18, 16, 13, 18, 16, 6,  16, 14, 11, 16, 14, 4,  14, 12, 9,
    14, 12, 2,  12, 10, 7, 12, 10, 0, 10, 8,  5,  10, 8,  13, 8,  6,  3,  8,  6,  11, 6,  4,  1,
    6,  4,  9,  4,  2,  6, 4,  2,  7, 2,  0,  4,  2,  0,  5,  0,  13, 2,  0,  21, 3,  21, 11, 0,
};
static const unsigned char p25_row_low[72] = {
    1, 3, 5, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 5,
    1, 3, 6, 0, 2, 5, 1, 3, 6, 0, 2, 5, 1, 3, 6, 0, 2, 5, 1, 3, 7, 0, 2, 5, 1, 4, 7, 0, 3, 5, 2, 4, 7, 1, 3, 5,
};
static const unsigned char p25_col_low[72] = {
    21, 19, 1, 21, 19, 9, 19, 17, 14, 19, 17, 7,  17, 15, 12, 17, 15, 5,  15, 13, 10, 15, 13, 3,
    13, 11, 8, 13, 11, 1, 11, 9,  6,  11, 9,  14, 9,  7,  4,  9,  7,  12, 7,  5,  2,  7,  5,  10,
    5,  3,  0, 5,  3,  8, 3,  1,  5,  3,  1,  6,  1,  14, 3,  1,  22, 4,  22, 12, 1,  22, 20, 2,
};
static const unsigned char dmr_row_high[36] = {
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
};
static const unsigned char dmr_col_high[36] = {
    23, 10, 22, 9, 21, 8,  20, 7, 19, 6, 18, 5, 17, 4, 16, 3, 15, 2,
    14, 1,  13, 0, 12, 10, 11, 9, 10, 8, 9,  7, 8,  6, 7,  5, 6,  4,
};
static const unsigned char dmr_row_low[36] = {
    0, 2, 0, 2, 0, 2, 0, 2, 0, 3, 0, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3,
};
static const unsigned char dmr_col_low[36] = {
    5,  3, 4,  2, 3,  1, 2,  0, 1,  13, 0,  12, 22, 11, 21, 10, 20, 9,
    19, 8, 18, 7, 17, 6, 16, 5, 15, 4,  14, 3,  13, 2,  12, 1,  11, 0,
};

static void
deinterleave(const unsigned char* bytes, int symbols, const unsigned char* row_high, const unsigned char* col_high,
             const unsigned char* row_low, const unsigned char* col_low, char* frame, int width) {
    for (int i = 0; i < symbols; ++i) {
        int high = 2 * i;
        int low = high + 1;
        frame[(row_high[i] * width) + col_high[i]] = (char)((bytes[high >> 3] >> (7 - (high & 7))) & 1u);
        frame[(row_low[i] * width) + col_low[i]] = (char)((bytes[low >> 3] >> (7 - (low & 7))) & 1u);
    }
}

int
mbe_quality_dvsi_frame_bytes(const char* codec) {
    if (!codec) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    if (strcmp(codec, "imbe7200") == 0) {
        return 18;
    }
    if (strcmp(codec, "ambe2400") == 0 || strcmp(codec, "ambe2450") == 0) {
        return 9;
    }
    return MBE_STATUS_INVALID_ARGUMENT;
}

int
mbe_quality_frame_from_dvsi(const char* codec, const unsigned char* bytes, size_t count, char* frame, size_t capacity) {
    int frame_bytes = mbe_quality_dvsi_frame_bytes(codec);
    if (frame_bytes < 0 || !bytes || !frame || count != (size_t)frame_bytes) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    union {
        char imbe7200[8][23];
        char ambe[4][24];
    } packed;

    memset(&packed, 0, sizeof(packed));
    int frame_count;
    if (strcmp(codec, "imbe7200") == 0) {
        frame_count = 184;
        deinterleave(bytes, 72, p25_row_high, p25_col_high, p25_row_low, p25_col_low, &packed.imbe7200[0][0], 23);
    } else if (strcmp(codec, "ambe2450") == 0) {
        frame_count = 96;
        deinterleave(bytes, 36, dmr_row_high, dmr_col_high, dmr_row_low, dmr_col_low, &packed.ambe[0][0], 24);
    } else {
        frame_count = 96;
        /* D-STAR vectors use the DV data byte order (LSB first). */
        int ret = mbe_decodeDStarDVData(bytes, packed.ambe);
        if (ret < 0) {
            return ret;
        }
    }
    if (capacity < (size_t)frame_count) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    memcpy(frame, &packed, (size_t)frame_count);
    return frame_count;
}
