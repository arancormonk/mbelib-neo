// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Exact weighted soft ECC over immutable, enumeration-ordered codebooks. */
#include <stddef.h>
#include <stdint.h>
#include "mbe_result.h"
#include "mbelib-neo/mbelib.h"

#include "mbe_ecc_codebooks.h"

/* The search visits codewords in groups of 64 that share their high data bits. */
#define GROUP_BITS 6u
#define GROUP_SIZE (1u << GROUP_BITS)

struct soft_code {
    const uint32_t* words; /* codewords in enumeration order */
    unsigned count;        /* 2^index_bits */
    const uint8_t* data;   /* codeword bit of each enumeration index bit */
    unsigned index_bits;
    int width;
    uint32_t data_mask; /* bits a decode reports; Golay keeps the received parity */
};

#define COUNT(array) ((unsigned)(sizeof(array) / sizeof((array)[0])))
static const struct soft_code golay_code = {soft_golay, COUNT(soft_golay), soft_golay_data, COUNT(soft_golay_data),
                                            23,         0x7ff800u};
static const struct soft_code hamming_code = {
    soft_hamming, COUNT(soft_hamming), soft_hamming_data, COUNT(soft_hamming_data), 15, 0x7fffu};
static const struct soft_code provoice_code = {
    soft_provoice, COUNT(soft_provoice), soft_provoice_data, COUNT(soft_provoice_data), 15, 0x7fffu};

struct soft_input {
    uint32_t received, hard;
    uint8_t weights[23];
};

/* A codeword's cost, the summed reliability of the bits where it differs from
 * the received word, splits by bits:
 *   cost = high[index >> GROUP_BITS] + low[index % GROUP_SIZE] + parity[word & parity_mask]
 * low and high cover the data bits through the enumeration index; parity
 * covers codeword bits below the highest parity bit, with data bits weighted
 * zero. At most 23 * 255 = 5865, so 16-bit entries hold every sum. */
struct soft_tables {
    uint16_t low[GROUP_SIZE], high[GROUP_SIZE], parity[2048];
    uint32_t parity_mask;
};

static unsigned
bit_count(uint32_t x) {
    x -= (x >> 1) & 0x55555555u;
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    x = (x + (x >> 4)) & 0x0f0f0f0fu;
    return (x * 0x01010101u) >> 24;
}

static void
soft_input_init(struct soft_input* s, const mbe_soft_bit* in, const char* hard, int width) {
    s->received = 0;
    s->hard = 0;
    for (int i = 0; i < width; ++i) {
        s->received |= (uint32_t)in[i].bit << i;
        s->hard |= (uint32_t)hard[i] << i;
        s->weights[i] = in[i].reliability;
    }
}

/* t[v] is the summed weight of the table bits b where bit b of v differs from
 * bit b of received, for every v below 2^bits. Doubling the table per bit adds
 * that bit's weight, or takes it away when the received bit is set (the
 * arithmetic wraps modulo 2^16, but every true sum fits). */
static void
subset_costs(uint16_t* t, unsigned bits, const uint8_t* weight, uint32_t received) {
    unsigned sum = 0;
    for (unsigned b = 0; b < bits; ++b) {
        sum += ((received >> b) & 1u) * weight[b];
    }
    t[0] = (uint16_t)sum;
    for (unsigned b = 0; b < bits; ++b) {
        const unsigned half = 1u << b;
        const uint16_t delta = (uint16_t)(((received >> b) & 1u) ? 0u - weight[b] : weight[b]);
        for (unsigned j = 0; j < half; ++j) {
            t[half + j] = (uint16_t)(t[j] + delta);
        }
    }
}

static void
soft_tables_init(struct soft_tables* t, const struct soft_input* s, const struct soft_code* code) {
    /* Each enumeration index bit's weight and received value, in index order. */
    uint8_t weight[23] = {0};
    uint32_t received = 0, data_bits = 0;
    for (unsigned b = 0; b < code->index_bits; ++b) {
        unsigned bit = code->data[b];
        data_bits |= 1u << bit;
        weight[b] = s->weights[bit];
        received |= ((s->received >> bit) & 1u) << b;
    }
    subset_costs(t->low, GROUP_BITS, weight, received);
    subset_costs(t->high, code->index_bits - GROUP_BITS, weight + GROUP_BITS, received >> GROUP_BITS);

    uint8_t parity_weight[23] = {0};
    unsigned parity_bits = 0;
    for (unsigned bit = 0; bit < (unsigned)code->width; ++bit) {
        if (!((data_bits >> bit) & 1u)) {
            parity_weight[bit] = s->weights[bit];
            parity_bits = bit + 1u;
        }
    }
    subset_costs(t->parity, parity_bits, parity_weight, s->received);
    t->parity_mask = (1u << parity_bits) - 1u;
}

/* Equal cost prefers the hard decision, then fewer protected-bit changes,
 * then the first enumerated candidate. Do not reorder the codebooks. */
static uint32_t
soft_search(const struct soft_input* s, const struct soft_code* code) {
    struct soft_tables t;
    soft_tables_init(&t, s, code);
    /* A codeword costing more than another codeword can never win, so the
     * codeword with the hard decision's data bits bounds the search. Skipping
     * only costlier codewords leaves the result exactly that of a full scan. */
    unsigned seed = 0;
    for (unsigned b = 0; b < code->index_bits; ++b) {
        seed |= ((s->hard >> code->data[b]) & 1u) << b;
    }
    unsigned bound =
        t.high[seed >> GROUP_BITS] + t.low[seed % GROUP_SIZE] + t.parity[code->words[seed] & t.parity_mask];
    unsigned best_score = 5866u, best_rank = 64u;
    uint32_t best = 0;
    for (unsigned first = 0; first < code->count; first += GROUP_SIZE) {
        const unsigned base = t.high[first >> GROUP_BITS];
        if (base > bound) {
            continue; /* every codeword in the group costs more */
        }
        const uint32_t* row = code->words + first;
        for (unsigned k = 0; k < GROUP_SIZE; ++k) {
            const unsigned score = base + t.low[k] + t.parity[row[k] & t.parity_mask];
            if (score > bound) {
                continue;
            }
            const uint32_t word = row[k];
            unsigned diffs = bit_count((word ^ s->received) & code->data_mask);
            unsigned rank = (((word ^ s->hard) & code->data_mask) != 0u ? 32u : 0u) + diffs;
            if (score < best_score || rank < best_rank) {
                best = word;
                best_score = score;
                best_rank = rank;
                bound = score;
            }
        }
    }
    return best;
}

static int
soft_decode(const mbe_soft_bit* in, char* out, const struct soft_code* code) {
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
    struct soft_input s;
    soft_input_init(&s, in, hard_out, width);
    uint32_t best = soft_search(&s, code);
    int diffs = (int)bit_count((best ^ s.received) & code->data_mask);
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
