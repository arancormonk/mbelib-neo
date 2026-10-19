// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Exact weighted soft-decision search over an enumeration-ordered codebook.
 *
 * ecc_soft_search_avx2.c builds this file a second time for AVX2 (mbe_cpu.h);
 * both copies return the same codeword. */
#include <stdint.h>

#include "mbe_compiler.h"
#include "mbe_cpu.h"
#include "mbe_ecc_soft_search.h"

#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#if defined(__AVX2__)
#include <immintrin.h>
#define SOFT_SCORE_LANES 8
#else
#include <emmintrin.h>
#define SOFT_SCORE_LANES 4
#endif
#endif

/* The search visits codewords in groups of 64 that share their high data bits. */
#define GROUP_BITS 6u
#define GROUP_SIZE (1u << GROUP_BITS)

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
soft_tables_init(struct soft_tables* t, const struct mbe_soft_input* s, const struct mbe_soft_code* code) {
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

struct soft_match {
    unsigned bound, score, rank;
    uint32_t word;
};

static void
soft_consider(const struct mbe_soft_input* s, const struct mbe_soft_code* code, uint32_t word, unsigned score,
              struct soft_match* best) {
    if (score > best->bound) {
        return;
    }
    unsigned diffs = mbe_soft_bit_count((word ^ s->received) & code->data_mask);
    unsigned rank = (((word ^ s->hard) & code->data_mask) != 0u ? 32u : 0u) + diffs;
    if (score < best->score || rank < best->rank) {
        best->word = word;
        best->score = score;
        best->rank = rank;
        best->bound = score;
    }
}

#if defined(SOFT_SCORE_LANES)
#if defined(__AVX2__)
/* Gather pairs of 16-bit parity costs, then select the wanted half. This
 * reads only complete pairs within parity[2048], including its last entry. */
static unsigned
soft_scores(const struct soft_tables* t, const uint32_t* words, unsigned k, unsigned limit,
            uint32_t scores[SOFT_SCORE_LANES]) {
    __m256i index = _mm256_and_si256(_mm256_loadu_si256((const __m256i*)words), _mm256_set1_epi32((int)t->parity_mask));
    __m256i pairs = _mm256_i32gather_epi32((const int*)(const void*)t->parity, _mm256_srli_epi32(index, 1), 4);
    __m256i shift = _mm256_slli_epi32(_mm256_and_si256(index, _mm256_set1_epi32(1)), 4);
    __m256i parity = _mm256_and_si256(_mm256_srlv_epi32(pairs, shift), _mm256_set1_epi32(0xffff));
    __m256i low = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i*)(t->low + k)));
    __m256i sum = _mm256_add_epi32(low, parity);
    _mm256_storeu_si256((__m256i*)scores, sum);
    __m256i above = _mm256_cmpgt_epi32(sum, _mm256_set1_epi32((int)limit));
    return (unsigned)~_mm256_movemask_ps(_mm256_castsi256_ps(above)) & 255u;
}
#else
static unsigned
soft_scores(const struct soft_tables* t, const uint32_t* words, unsigned k, unsigned limit,
            uint32_t scores[SOFT_SCORE_LANES]) {
    __m128i parity = _mm_setr_epi32(t->parity[words[0] & t->parity_mask], t->parity[words[1] & t->parity_mask],
                                    t->parity[words[2] & t->parity_mask], t->parity[words[3] & t->parity_mask]);
    __m128i low = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i*)(t->low + k)), _mm_setzero_si128());
    __m128i sum = _mm_add_epi32(low, parity);
    _mm_storeu_si128((__m128i*)scores, sum);
    __m128i above = _mm_cmpgt_epi32(sum, _mm_set1_epi32((int)limit));
    return (unsigned)~_mm_movemask_ps(_mm_castsi128_ps(above)) & 15u;
}
#endif
#endif

static void
soft_search_group(const struct mbe_soft_input* s, const struct mbe_soft_code* code, const struct soft_tables* t,
                  unsigned first, struct soft_match* best) {
    const unsigned base = t->high[first >> GROUP_BITS];
    if (base > best->bound) {
        return;
    }
    const uint32_t* row = code->words + first;
#if defined(SOFT_SCORE_LANES)
    for (unsigned k = 0; k < GROUP_SIZE; k += SOFT_SCORE_LANES) {
        uint32_t scores[SOFT_SCORE_LANES];
        unsigned mask = soft_scores(t, row + k, k, best->bound - base, scores);
        if (mask == 0u) {
            continue;
        }
        for (unsigned lane = 0; lane < SOFT_SCORE_LANES; ++lane) {
            if ((mask >> lane) & 1u) {
                soft_consider(s, code, row[k + lane], base + scores[lane], best);
            }
        }
    }
#else
    for (unsigned k = 0; k < GROUP_SIZE; ++k) {
        unsigned score = base + t->low[k] + t->parity[row[k] & t->parity_mask];
        soft_consider(s, code, row[k], score, best);
    }
#endif
}

/* Equal cost prefers the hard decision, then fewer protected-bit changes,
 * then the first enumerated candidate. Do not reorder the codebooks. */
uint32_t
MBE_KFN(mbe_ecc_soft_search)(const struct mbe_soft_input* s, const struct mbe_soft_code* code) {
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
    struct soft_match best = {bound, 5866u, 64u, 0};
    for (unsigned first = 0; first < code->count; first += GROUP_SIZE) {
        soft_search_group(s, code, &t, first, &best);
    }
    return best.word;
}

#if defined(MBE_AVX2_DISPATCH)
uint32_t
mbe_ecc_soft_search(const struct mbe_soft_input* s, const struct mbe_soft_code* code) {
    if (mbe_cpu_use_avx2()) {
        return mbe_ecc_soft_search_avx2(s, code);
    }
    return mbe_ecc_soft_search_sse2(s, code);
}
#endif
