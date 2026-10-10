// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Exact weighted soft ECC over immutable, enumeration-ordered codebooks. */
#include <stdint.h>
#include <string.h>
#include "mbe_compiler.h"
#include "mbe_result.h"
#include "mbelib-neo/mbelib.h"

#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#include <emmintrin.h>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
#define MBE_ECC_VECTOR 1
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h>
#define MBE_ECC_VECTOR 1
#endif

#if defined(MBE_ECC_VECTOR)
typedef uint8_t mbe_ecc_word23[32];
typedef uint8_t mbe_ecc_word15[16];
#define MBE_ECC_MASK(v, b) ((uint8_t)(0u - (((v) >> (b)) & 1u)))
#define MBE_ECC_LOW(v)                                                                                                 \
    MBE_ECC_MASK(v, 0), MBE_ECC_MASK(v, 1), MBE_ECC_MASK(v, 2), MBE_ECC_MASK(v, 3), MBE_ECC_MASK(v, 4),                \
        MBE_ECC_MASK(v, 5), MBE_ECC_MASK(v, 6), MBE_ECC_MASK(v, 7), MBE_ECC_MASK(v, 8), MBE_ECC_MASK(v, 9),            \
        MBE_ECC_MASK(v, 10), MBE_ECC_MASK(v, 11), MBE_ECC_MASK(v, 12), MBE_ECC_MASK(v, 13), MBE_ECC_MASK(v, 14),       \
        MBE_ECC_MASK(v, 15)
#define MBE_ECC_WORD15(v) {MBE_ECC_LOW(v)}
#define MBE_ECC_WORD23(v)                                                                                              \
    {MBE_ECC_LOW(v),                                                                                                   \
     MBE_ECC_MASK(v, 16),                                                                                              \
     MBE_ECC_MASK(v, 17),                                                                                              \
     MBE_ECC_MASK(v, 18),                                                                                              \
     MBE_ECC_MASK(v, 19),                                                                                              \
     MBE_ECC_MASK(v, 20),                                                                                              \
     MBE_ECC_MASK(v, 21),                                                                                              \
     MBE_ECC_MASK(v, 22),                                                                                              \
     0,                                                                                                                \
     0,                                                                                                                \
     0,                                                                                                                \
     0,                                                                                                                \
     0,                                                                                                                \
     0,                                                                                                                \
     0,                                                                                                                \
     0,                                                                                                                \
     0}
#else
typedef uint32_t mbe_ecc_word23;
typedef uint32_t mbe_ecc_word15;
#define MBE_ECC_WORD15(v) (v)
#define MBE_ECC_WORD23(v) (v)
#endif
#include "mbe_ecc_codebooks.h"
#undef MBE_ECC_WORD15
#undef MBE_ECC_WORD23

struct soft_input {
    uint32_t received, hard, data_mask;
#if defined(MBE_ECC_VECTOR)
    uint8_t bits[32], weights[32];
#else
    unsigned costs[3][256];
#endif
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
    memset(s, 0, sizeof(*s));
    s->data_mask = width == 23 ? 0x7ff800u : 0x7fffu;
    for (int i = 0; i < width; ++i) {
        s->received |= (uint32_t)in[i].bit << i;
        s->hard |= (uint32_t)hard[i] << i;
#if defined(MBE_ECC_VECTOR)
        s->bits[i] = (uint8_t)(0u - (unsigned)in[i].bit);
        s->weights[i] = in[i].reliability;
#else
        /* Subset-sum costs for each byte of candidate XOR received. */
        int byte = i / 8;
        unsigned bit = 1u << (i % 8);
        for (unsigned j = 0; j < bit; ++j) {
            s->costs[byte][j + bit] = s->costs[byte][j] + in[i].reliability;
        }
#endif
    }
}

#if defined(MBE_ECC_VECTOR)
/* At most 23 * 255 = 5865: neither byte SAD nor its widened sums overflow. */
static unsigned
candidate_cost(const struct soft_input* s, const uint8_t* row, int width) {
#if defined(MBE_SIMD_TARGET_SSE2)
#if defined(__AVX2__)
    if (width == 23) {
        __m256i bits = _mm256_loadu_si256((const __m256i*)s->bits);
        __m256i weights = _mm256_loadu_si256((const __m256i*)s->weights);
        __m256i v = _mm256_loadu_si256((const __m256i*)row);
        __m256i sad = _mm256_sad_epu8(_mm256_and_si256(_mm256_xor_si256(v, bits), weights), _mm256_setzero_si256());
        __m128i sum = _mm_add_epi64(_mm256_castsi256_si128(sad), _mm256_extracti128_si256(sad, 1));
        return (unsigned)_mm_cvtsi128_si32(_mm_add_epi64(sum, _mm_srli_si128(sum, 8)));
    }
#endif
    __m128i sum = _mm_setzero_si128();
    for (int i = 0; i < width; i += 16) {
        __m128i v = _mm_loadu_si128((const __m128i*)(row + i));
        __m128i bits = _mm_loadu_si128((const __m128i*)(s->bits + i));
        __m128i weights = _mm_loadu_si128((const __m128i*)(s->weights + i));
        sum = _mm_add_epi64(sum, _mm_sad_epu8(_mm_and_si128(_mm_xor_si128(v, bits), weights), _mm_setzero_si128()));
    }
    return (unsigned)_mm_cvtsi128_si32(_mm_add_epi64(sum, _mm_srli_si128(sum, 8)));
#else
    uint16x8_t sum = vdupq_n_u16(0);
    for (int i = 0; i < width; i += 16) {
        uint8x16_t diff = veorq_u8(vld1q_u8(row + i), vld1q_u8(s->bits + i));
        sum = vaddq_u16(sum, vpaddlq_u8(vandq_u8(diff, vld1q_u8(s->weights + i))));
    }
    uint64x2_t total = vpaddlq_u32(vpaddlq_u16(sum));
    return (unsigned)(vgetq_lane_u64(total, 0) + vgetq_lane_u64(total, 1));
#endif
}

static uint32_t
candidate_word(const uint8_t* row, int width) {
#if defined(MBE_SIMD_TARGET_SSE2)
    uint32_t word = (uint32_t)_mm_movemask_epi8(_mm_loadu_si128((const __m128i*)row));
    if (width == 23) {
        word |= (uint32_t)_mm_movemask_epi8(_mm_loadu_si128((const __m128i*)(row + 16))) << 16;
    }
    return word;
#else
    uint32_t word = 0;
    for (int i = 0; i < width; ++i) {
        word |= (uint32_t)(row[i] & 1u) << i;
    }
    return word;
#endif
}
#else
static uint32_t
candidate_word(const uint8_t* row, int width) {
    uint32_t word;
    (void)width;
    memcpy(&word, row, sizeof(word));
    return word;
}

static unsigned
candidate_cost(const struct soft_input* s, const uint8_t* row, int width) {
    uint32_t diff = candidate_word(row, width) ^ s->received;
    return s->costs[0][diff & 255u] + s->costs[1][(diff >> 8) & 255u] + s->costs[2][diff >> 16];
}
#endif

/* Equal cost prefers the hard decision, then fewer protected-bit changes,
 * then the first enumerated candidate. Do not reorder the codebooks. */
static uint32_t
soft_search(const struct soft_input* s, const void* words, int count, int width, size_t stride) {
    unsigned best_score = 5866u, best_rank = 64u;
    uint32_t best = 0;
    const uint8_t* row = words;
    for (int i = 0; i < count; ++i, row += stride) {
        unsigned score = candidate_cost(s, row, width);
        if (score > best_score) {
            continue;
        }
        uint32_t word = candidate_word(row, width);
        unsigned diffs = bit_count((word ^ s->received) & s->data_mask);
        unsigned rank = (((word ^ s->hard) & s->data_mask) != 0u ? 32u : 0u) + diffs;
        if (score < best_score || rank < best_rank) {
            best = word;
            best_score = score;
            best_rank = rank;
        }
    }
    return best;
}

static int
soft_decode(const mbe_soft_bit* in, char* out, int width, int provoice) {
    if (!out) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    int status = mbe_validate_soft_bits(in, (size_t)width);
    if (status < 0) {
        return status;
    }
    char hard_in[23], hard_out[23];
    for (int i = 0; i < width; ++i) {
        hard_in[i] = (char)in[i].bit;
    }
    if (width == 23) {
        (void)mbe_golay2312(hard_in, hard_out);
    } else if (provoice) {
        (void)mbe_7100x4400hamming1511(hard_in, hard_out);
    } else {
        (void)mbe_hamming1511(hard_in, hard_out);
    }
    struct soft_input s;
    soft_input_init(&s, in, hard_out, width);
    uint32_t best;
    if (width == 23) {
        best = soft_search(&s, soft_golay, 4096, width, sizeof(soft_golay[0]));
    } else {
        best = soft_search(&s, provoice ? soft_provoice : soft_hamming, 2048, width, sizeof(soft_hamming[0]));
    }
    int diffs = (int)bit_count((best ^ s.received) & s.data_mask);
    /* Golay corrects/reports data only; preserve the received parity bits. */
    best = (best & s.data_mask) | (s.received & ~s.data_mask);
    for (int i = 0; i < width; ++i) {
        out[i] = (char)((best >> i) & 1u);
    }
    return diffs;
}

int
mbe_golay2312Soft(const mbe_soft_bit* in, char* out) {
    return soft_decode(in, out, 23, 0);
}

int
mbe_hamming1511Soft(const mbe_soft_bit* in, char* out) {
    return soft_decode(in, out, 15, 0);
}

int
mbe_7100x4400hamming1511Soft(const mbe_soft_bit* in, char* out) {
    return soft_decode(in, out, 15, 1);
}
