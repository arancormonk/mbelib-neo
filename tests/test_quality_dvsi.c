// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief DVSI test-vector import: deinterleave structure and bit preservation.
 */

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "mbe_quality_dvsi.h"
#include "mbe_quality_frames.h"
#include "mbelib-neo/mbelib.h"

struct layout {
    const char* codec;
    int bytes;
    int width;
    int rows;
    int lsb_first;
    int row_bits[8];
};

static const struct layout layouts[] = {
    {"imbe7200", 18, 23, 8, 0, {23, 23, 23, 23, 15, 15, 15, 7}},
    {"ambe2450", 9, 24, 4, 0, {24, 23, 11, 14}},
    {"ambe2400", 9, 24, 4, 1, {24, 23, 11, 14}},
};

static void
set_input_bit(const struct layout* l, unsigned char* bytes, int bit) {
    int shift = l->lsb_first ? (bit & 7) : (7 - (bit & 7));
    bytes[bit >> 3] ^= (unsigned char)(1u << shift);
}

/* Map each input bit to the single frame cell it lands in. */
static void
one_hot_map(const struct layout* l, int cell_of[144]) {
    for (int bit = 0; bit < l->bytes * 8; ++bit) {
        unsigned char bytes[18] = {0};
        char frame[184];
        set_input_bit(l, bytes, bit);
        int count = mbe_quality_frame_from_dvsi(l->codec, bytes, (size_t)l->bytes, frame, sizeof(frame));
        assert(count == l->rows * l->width);
        int found = -1;
        for (int c = 0; c < count; ++c) {
            assert(frame[c] == 0 || frame[c] == 1);
            if (frame[c]) {
                assert(found < 0);
                found = c;
            }
        }
        assert(found >= 0);
        cell_of[bit] = found;
    }
}

static void
test_bijection(const struct layout* l) {
    int cell_of[144];
    char used[184] = {0};
    one_hot_map(l, cell_of);
    for (int bit = 0; bit < l->bytes * 8; ++bit) {
        int row = cell_of[bit] / l->width;
        int col = cell_of[bit] % l->width;
        assert(!used[cell_of[bit]]);
        used[cell_of[bit]] = 1;
        assert(col < l->row_bits[row]);
    }
    int total = 0;
    for (int r = 0; r < l->rows; ++r) {
        total += l->row_bits[r];
    }
    assert(total == l->bytes * 8);
}

/* TIA-102.BABA 7.5: bits of one code vector are at least 3 dibits apart. */
static void
test_p25_symbol_separation(void) {
    int cell_of[144];
    one_hot_map(&layouts[0], cell_of);
    for (int a = 0; a < 144; ++a) {
        for (int b = a + 1; b < 144; ++b) {
            if (cell_of[a] / 23 == cell_of[b] / 23) {
                assert((b / 2) - (a / 2) >= 3);
            }
        }
    }
}

static void
interleave(const struct layout* l, const char* frame, unsigned char* bytes) {
    int cell_of[144];
    one_hot_map(l, cell_of);
    memset(bytes, 0, (size_t)l->bytes);
    for (int bit = 0; bit < l->bytes * 8; ++bit) {
        if (frame[cell_of[bit]]) {
            set_input_bit(l, bytes, bit);
        }
    }
}

static int
process_errors(const char* codec, char* frame) {
    float pcm[160];
    mbe_process_result result;
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    char data[88];
    int status;
    if (strcmp(codec, "imbe7200") == 0) {
        status = mbe_processImbe7200x4400Framef(pcm, &result, (const char (*)[23])frame, data, &cur, &prev, &enhanced);
    } else if (strcmp(codec, "ambe2450") == 0) {
        status = mbe_processAmbe3600x2450Framef(pcm, &result, (const char (*)[24])frame, data, &cur, &prev, &enhanced);
    } else {
        status = mbe_processAmbe3600x2400Framef(pcm, &result, (const char (*)[24])frame, data, &cur, &prev, &enhanced);
    }
    assert(status >= 0);
    return result.total_errors;
}

/* Fixture frames survive interleave+import bit-exactly and decode cleanly. */
static void
test_round_trip(const struct layout* l) {
    unsigned int seed = 12345u;
    for (int trial = 0; trial < 32; ++trial) {
        char data[88];
        char frame[184];
        char imported[184];
        unsigned char bytes[18];
        size_t data_bits = l->width == 23 ? 88u : 49u;
        for (size_t i = 0; i < data_bits; ++i) {
            seed = (seed * 1103515245u) + 12345u;
            data[i] = (char)((seed >> 16) & 1u);
        }
        int count = mbe_quality_frame_from_data(l->codec, data, data_bits, frame, sizeof(frame));
        assert(count == l->rows * l->width);
        interleave(l, frame, bytes);
        assert(mbe_quality_frame_from_dvsi(l->codec, bytes, (size_t)l->bytes, imported, sizeof(imported)) == count);
        assert(memcmp(frame, imported, (size_t)count) == 0);
        assert(process_errors(l->codec, imported) == 0);

        /* A flipped channel bit is carried through, not corrected. */
        set_input_bit(l, bytes, trial % (l->bytes * 8));
        char flipped[184];
        assert(mbe_quality_frame_from_dvsi(l->codec, bytes, (size_t)l->bytes, flipped, sizeof(flipped)) == count);
        int differing = 0;
        for (int c = 0; c < count; ++c) {
            differing += frame[c] != flipped[c];
        }
        assert(differing == 1);
    }
}

static void
test_arguments(void) {
    unsigned char bytes[18] = {0};
    char frame[184];
    assert(mbe_quality_dvsi_frame_bytes("imbe7100") == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_quality_dvsi_frame_bytes(NULL) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_quality_frame_from_dvsi("imbe7200", bytes, 17, frame, sizeof(frame)) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_quality_frame_from_dvsi("ambe2400", NULL, 9, frame, sizeof(frame)) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_quality_frame_from_dvsi("ambe2400", bytes, 9, NULL, sizeof(frame)) == MBE_STATUS_INVALID_ARGUMENT);
    memset(frame, 7, sizeof(frame));
    assert(mbe_quality_frame_from_dvsi("imbe7200", bytes, 18, frame, 183) == MBE_STATUS_INVALID_ARGUMENT);
    assert(frame[0] == 7);
}

int
main(void) {
    for (size_t i = 0; i < sizeof(layouts) / sizeof(layouts[0]); ++i) {
        test_bijection(&layouts[i]);
        test_round_trip(&layouts[i]);
    }
    test_p25_symbol_separation();
    test_arguments();
    puts("DVSI import: ok");
    return 0;
}
