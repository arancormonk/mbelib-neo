// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Exact weighted soft ECC over immutable, enumeration-ordered codebooks. */
#include <stddef.h>
#include <stdint.h>
#include "mbe_result.h"
#include "mbelib-neo/mbelib.h"

#include "mbe_ecc_codebooks.h"

#include "mbe_ecc_soft_search.h"

#define COUNT(array) ((unsigned)(sizeof(array) / sizeof((array)[0])))
static const struct mbe_soft_code golay_code = {soft_golay, COUNT(soft_golay), soft_golay_data, COUNT(soft_golay_data),
                                                23,         0x7ff800u};
static const struct mbe_soft_code hamming_code = {
    soft_hamming, COUNT(soft_hamming), soft_hamming_data, COUNT(soft_hamming_data), 15, 0x7fffu};
static const struct mbe_soft_code provoice_code = {
    soft_provoice, COUNT(soft_provoice), soft_provoice_data, COUNT(soft_provoice_data), 15, 0x7fffu};

static void
soft_input_init(struct mbe_soft_input* s, const mbe_soft_bit* in, const char* hard, int width) {
    s->received = 0;
    s->hard = 0;
    for (int i = 0; i < width; ++i) {
        s->received |= (uint32_t)in[i].bit << i;
        s->hard |= (uint32_t)hard[i] << i;
        s->weights[i] = in[i].reliability;
    }
}

static int
soft_decode(const mbe_soft_bit* in, char* out, const struct mbe_soft_code* code) {
    if (!out) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    const int width = code->width;
    int status = mbe_validate_soft_bits(in, (size_t)width);
    if (status < 0) {
        return status;
    }
    char hard_in[23], hard_out[23];
    for (int i = 0; i < width; ++i) {
        hard_in[i] = (char)in[i].bit;
    }
    if (code == &golay_code) {
        (void)mbe_golay2312(hard_in, hard_out);
    } else if (code == &provoice_code) {
        (void)mbe_7100x4400hamming1511(hard_in, hard_out);
    } else {
        (void)mbe_hamming1511(hard_in, hard_out);
    }
    struct mbe_soft_input s;
    soft_input_init(&s, in, hard_out, width);
    uint32_t best = mbe_ecc_soft_search(&s, code);
    int diffs = (int)mbe_soft_bit_count((best ^ s.received) & code->data_mask);
    /* Golay corrects/reports data only; preserve the received parity bits. */
    best = (best & code->data_mask) | (s.received & ~code->data_mask);
    for (int i = 0; i < width; ++i) {
        out[i] = (char)((best >> i) & 1u);
    }
    return diffs;
}

int
mbe_golay2312Soft(const mbe_soft_bit* in, char* out) {
    return soft_decode(in, out, &golay_code);
}

int
mbe_hamming1511Soft(const mbe_soft_bit* in, char* out) {
    return soft_decode(in, out, &hamming_code);
}

int
mbe_7100x4400hamming1511Soft(const mbe_soft_bit* in, char* out) {
    return soft_decode(in, out, &provoice_code);
}
