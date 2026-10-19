// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Private exact weighted soft-decision codebook search, shared by ecc_soft.c and its search units. */
#ifndef MBELIB_NEO_INTERNAL_MBE_ECC_SOFT_SEARCH_H
#define MBELIB_NEO_INTERNAL_MBE_ECC_SOFT_SEARCH_H

#include <stdint.h>

#include "mbe_cpu.h" // IWYU pragma: keep (MBE_AVX2_DISPATCH declares the dispatched copies)

struct mbe_soft_code {
    const uint32_t* words; /* codewords in enumeration order */
    unsigned count;        /* 2^index_bits */
    const uint8_t* data;   /* codeword bit of each enumeration index bit */
    unsigned index_bits;
    int width;
    uint32_t data_mask; /* bits a decode reports; Golay keeps the received parity */
};

struct mbe_soft_input {
    uint32_t received, hard;
    uint8_t weights[23];
};

static inline unsigned
mbe_soft_bit_count(uint32_t x) {
    x -= (x >> 1) & 0x55555555u;
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    x = (x + (x >> 4)) & 0x0f0f0f0fu;
    return (x * 0x01010101u) >> 24;
}

/** The least-cost codeword for s, with the tie-breaks documented in ecc_soft_search.c. Do not reorder the
 * codebooks. */
uint32_t mbe_ecc_soft_search(const struct mbe_soft_input* s, const struct mbe_soft_code* code);

/* The baseline and AVX2 copies behind the run-time dispatch (mbe_cpu.h). */
#if defined(MBE_AVX2_DISPATCH) || defined(MBE_KERNEL_AVX2)
uint32_t mbe_ecc_soft_search_sse2(const struct mbe_soft_input* s, const struct mbe_soft_code* code);
uint32_t mbe_ecc_soft_search_avx2(const struct mbe_soft_input* s, const struct mbe_soft_code* code);
#endif

#endif /* MBELIB_NEO_INTERNAL_MBE_ECC_SOFT_SEARCH_H */
