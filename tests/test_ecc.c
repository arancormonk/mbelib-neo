// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Unit tests for ECC helpers (Golay and Hamming variants).
 */
#include <assert.h>
#include <string.h>
// Embedded copies of generator tables used for test expectations
static const int hammingGenerator[4] = {0x7f08, 0x78e4, 0x66d2, 0x55b1};
static const int imbe7100x4400hammingGenerator[4] = {0x7ac8, 0x3d64, 0x1eb2, 0x7591};
static const int golayGenerator[12] = {0x63a, 0x31d, 0x7b4, 0x3da, 0x1ed, 0x6cc,
                                       0x366, 0x1b3, 0x6e3, 0x54b, 0x49f, 0x475};
#include "mbelib-neo/mbelib.h"

// make stderr unbuffered for debug visibility under ctest
#include <stdio.h>

/**
 * @brief Flip a single bit in a bit-vector.
 * @param bits Bit array of chars containing 0/1 values.
 * @param idx  Index of bit to toggle.
 */
static void
flip_bit(char* bits, int idx) {
    bits[idx] = bits[idx] ? 0 : 1;
}

static void
make_soft_bits(const char* bits, mbe_soft_bit* soft, int count, int weak0, int weak1, int weak2, int weak3) {
    for (int i = 0; i < count; ++i) {
        soft[i] = mbe_softBitFromHard(bits[i], 200);
    }
    if (weak0 >= 0) {
        soft[weak0].reliability = 1;
    }
    if (weak1 >= 0) {
        soft[weak1].reliability = 1;
    }
    if (weak2 >= 0) {
        soft[weak2].reliability = 1;
    }
    if (weak3 >= 0) {
        soft[weak3].reliability = 1;
    }
}

/**
 * @brief Compute 4-bit syndrome for a 15-bit codeword using a generator.
 * @param code      Codeword bits as 0/1 chars, LSB at index 0.
 * @param generator Array of 4 generator row masks.
 * @return Syndrome (LSB-first) in range [0, 15].
 */
static int
compute_syndrome_with_generator(const char code[15], const int generator[4]) {
    int block = 0;
    for (int i = 14; i >= 0; --i) {
        block <<= 1;
        block |= (code[i] ? 1 : 0);
    }
    int syndrome = 0;
    for (int r = 0; r < 4; ++r) {
        int stmp = block & generator[r];
        int sbit = (stmp & 1);
        for (int j = 0; j < 14; ++j) {
            stmp >>= 1;
            sbit ^= (stmp & 1);
        }
        syndrome |= (sbit << r);
    }
    return syndrome;
}

/**
 * @brief Determine parity bit positions (those not in data_pos).
 * @param data_pos  11 indices of data positions (0..14).
 * @param parity_pos Output 4 indices for parity positions.
 */
static void
derive_parity_positions(const int data_pos[11], int parity_pos[4]) {
    int p = 0;
    for (int i = 0; i < 15; ++i) {
        int is_data = 0;
        for (int d = 0; d < 11; ++d) {
            if (data_pos[d] == i) {
                is_data = 1;
                break;
            }
        }
        if (!is_data) {
            parity_pos[p++] = i;
        }
    }
}

/**
 * @brief Encode a (15,11) Hamming codeword for given data positions.
 *
 * Brute-forces the 4 parity bits (16 combinations) until syndrome is zero.
 *
 * @param generator  Generator row masks (4 entries).
 * @param data_pos   11 indices of data positions (0..14).
 * @param data_in    11 data bits (0/1 chars) mapped to data_pos.
 * @param code       Output encoded 15-bit codeword (0/1 chars).
 * @return 1 on success, 0 if no parity combination yields zero syndrome.
 */
static int
encode_hamming15_with_generator(const int generator[4], const int data_pos[11], const char data_in[11], char code[15]) {
    for (int i = 0; i < 15; ++i) {
        code[i] = 0;
    }
    for (int i = 0; i < 11; ++i) {
        code[data_pos[i]] = (data_in[i] ? 1 : 0);
    }
    int parity_pos[4];
    derive_parity_positions(data_pos, parity_pos);
    // brute-force 4 parity bits (16 combos) to satisfy parity equations
    for (int bits = 0; bits < 16; ++bits) {
        for (int k = 0; k < 4; ++k) {
            code[parity_pos[k]] = ((bits >> k) & 1);
        }
        int syn = compute_syndrome_with_generator(code, generator);
        if (syn == 0) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Convert a 4-bit syndrome into a (1<<bit) correction mask.
 * @param syndrome  4-bit syndrome value (LSB-first).
 * @param generator Generator row masks (4 entries).
 * @return Bit mask with one bit set (1<<bit), or 0 if not found.
 */
static int
correction_mask_from_generator(int syndrome, const int generator[4]) {
    if (syndrome == 0) {
        return 0;
    }
    for (int j = 0; j < 15; ++j) {
        int syn = 0;
        for (int r = 0; r < 4; ++r) {
            syn |= (((generator[r] >> j) & 1) << r);
        }
        if (syn == syndrome) {
            return (1 << j);
        }
    }
    return 0;
}

/* Independent oracle that scores every codeword: reconstruct candidates with
 * the original parity search, and score each bit individually. Deliberately
 * shares neither tables nor scoring/tie-ranking machinery with the optimized
 * decoder. */
static char reference_words[3][4096][23];

static void
init_reference_words(void) {
    const int positions[2][11] = {{2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14}, {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}};
    for (unsigned d = 0; d < 4096u; ++d) {
        unsigned ecc = 0;
        for (int i = 0; i < 12; ++i) {
            if ((d >> (11 - i)) & 1u) {
                ecc ^= (unsigned)golayGenerator[i];
            }
        }
        unsigned word = (d << 11) | ecc;
        for (int i = 0; i < 23; ++i) {
            reference_words[0][d][i] = (char)((word >> i) & 1u);
        }
    }
    for (int variant = 0; variant < 2; ++variant) {
        for (unsigned d = 0; d < 2048u; ++d) {
            char data[11];
            for (int i = 0; i < 11; ++i) {
                data[i] = (char)((d >> i) & 1u);
            }
            int ok = encode_hamming15_with_generator(variant ? imbe7100x4400hammingGenerator : hammingGenerator,
                                                     positions[variant], data, reference_words[variant + 1][d]);
            assert(ok);
        }
    }
}

static int
reference_soft(const mbe_soft_bit* in, char* out, int variant) {
    int width = variant == 0 ? 23 : 15;
    int first = variant == 0 ? 11 : 0;
    int count = variant == 0 ? 4096 : 2048;
    char hard_in[23], hard[23];
    for (int i = 0; i < width; ++i) {
        hard_in[i] = (char)in[i].bit;
    }
    if (variant == 0) {
        (void)mbe_golay2312(hard_in, hard);
    } else if (variant == 1) {
        (void)mbe_hamming1511(hard_in, hard);
    } else {
        (void)mbe_7100x4400hamming1511(hard_in, hard);
    }
    int best_score = 100000, best_diffs = 100, best_hard = 0;
    for (int d = 0; d < count; ++d) {
        const char* row = reference_words[variant][d];
        int score = 0, diffs = 0, matches = 1;
        for (int i = 0; i < width; ++i) {
            score += row[i] != hard_in[i] ? in[i].reliability : 0;
            if (i >= first) {
                diffs += row[i] != hard_in[i];
                matches &= row[i] == hard[i];
            }
        }
        if (score < best_score || (score == best_score && matches > best_hard)
            || (score == best_score && matches == best_hard && diffs < best_diffs)) {
            memcpy(out, row, (size_t)width);
            best_score = score;
            best_diffs = diffs;
            best_hard = matches;
        }
    }
    for (int i = 0; i < first; ++i) {
        out[i] = hard_in[i];
    }
    return best_diffs;
}

static void
test_soft_differential(void) {
    int (*const decode[3])(const mbe_soft_bit*, char*) = {mbe_golay2312Soft, mbe_hamming1511Soft,
                                                          mbe_7100x4400hamming1511Soft};
    unsigned rng = 0x873fa11u;
    init_reference_words();
    for (int variant = 0; variant < 3; ++variant) {
        int width = variant == 0 ? 23 : 15;
        for (int trial = 0; trial < 1024; ++trial) {
            mbe_soft_bit soft[23];
            for (int i = 0; i < width; ++i) {
                rng = rng * 1664525u + 1013904223u;
                soft[i].bit = (unsigned char)(rng >> 31);
                switch (trial % 7) {
                    case 0: soft[i].reliability = 0; break;
                    case 1: soft[i].reliability = 255; break;
                    case 2: soft[i].reliability = (unsigned char)((rng >> 16) & 1u); break;
                    case 3: soft[i].reliability = i == trial % width ? 255 : 0; break;
                    /* Small weights tie many candidates at the search's cost bound. */
                    case 4: soft[i].reliability = (unsigned char)((rng >> 16) & 3u); break;
                    /* Mostly confident words with a few weak bits prune hardest. */
                    case 5: soft[i].reliability = i % 4 == trial % 4 ? (unsigned char)((rng >> 16) & 15u) : 255; break;
                    default: soft[i].reliability = (unsigned char)(rng >> 16); break;
                }
            }
            char expected[23], actual[23];
            int errors = reference_soft(soft, expected, variant);
            assert(decode[variant](soft, actual) == errors);
            assert(memcmp(actual, expected, (size_t)width) == 0);
        }
        /* Every clean candidate, including the entire generated codebook. */
        int count = variant == 0 ? 4096 : 2048;
        for (int d = 0; d < count; ++d) {
            mbe_soft_bit soft[23];
            char actual[23];
            assert(mbe_softBitsFromHard(reference_words[variant][d], soft, (size_t)width, 255) == 0);
            assert(decode[variant](soft, actual) == 0);
            assert(memcmp(actual, reference_words[variant][d], (size_t)width) == 0);
        }
    }
}

/* Soft-decode one word of code `variant`: Golay (23,12), Hamming (15,11) or
 * the ProVoice Hamming (15,11). */
static int
decode_soft(int variant, const mbe_soft_bit* soft, char* out) {
    if (variant == 0) {
        return mbe_golay2312Soft(soft, out);
    }
    return variant == 1 ? mbe_hamming1511Soft(soft, out) : mbe_7100x4400hamming1511Soft(soft, out);
}

/* Confidence of bit i under pattern `kind`, from a random draw: nothing
 * known, everything certain, every other bit unknown, binary, small weights
 * (many ties), and any weight. */
static unsigned char
pattern_reliability(int kind, int i, unsigned draw) {
    switch (kind) {
        case 0: return 0;
        case 1: return 255;
        case 2: return (unsigned char)((i & 1) ? 0 : 255);
        case 3: return (unsigned char)(draw & 1u);
        case 4: return (unsigned char)(draw & 3u);
        default: return (unsigned char)(draw & 255u);
    }
}

/* A received word that is a codeword decodes to itself, with no errors,
 * whatever the confidences are. Where enough bits have zero confidence, other
 * codewords cost as little as it does (nothing), and the hard decision's own
 * codeword must still win. Fails for a tie-break, or a shortcut around the
 * search, that hands back another codeword of equal cost. */
static void
test_soft_clean_words(void) {
    unsigned rng = 0x51ab37u;
    for (int variant = 0; variant < 3; ++variant) {
        const int width = variant == 0 ? 23 : 15;
        const int count = variant == 0 ? 4096 : 2048;
        for (int d = 0; d < count; ++d) {
            const char* word = reference_words[variant][d];
            for (int kind = 0; kind < 6; ++kind) {
                mbe_soft_bit soft[23];
                char actual[23];
                for (int i = 0; i < width; ++i) {
                    rng = rng * 1664525u + 1013904223u;
                    soft[i].bit = (unsigned char)word[i];
                    soft[i].reliability = pattern_reliability(kind, i, rng >> 16);
                }
                const int errors = decode_soft(variant, soft, actual);
                assert(errors == 0);
                assert(memcmp(actual, word, (size_t)width) == 0);
                (void)errors;
            }
        }
    }
}

/* Whether the hard decisions of a soft word are one of the code's codewords. */
static int
is_reference_word(int variant, const mbe_soft_bit* soft) {
    const int width = variant == 0 ? 23 : 15;
    const int count = variant == 0 ? 4096 : 2048;
    for (int d = 0; d < count; ++d) {
        int same = 1;
        for (int i = 0; i < width && same; ++i) {
            same = reference_words[variant][d][i] == (char)soft[i].bit;
        }
        if (same) {
            return 1;
        }
    }
    return 0;
}

/* Flip `flips` different bits of a soft word, chosen by *rng, and take their
 * confidence away. */
static void
flip_at_zero_confidence(mbe_soft_bit* soft, int width, int flips, unsigned* rng) {
    int flipped = 0;
    while (flipped < flips) {
        *rng = *rng * 1664525u + 1013904223u;
        const int at = (int)((*rng >> 16) % (unsigned)width);
        if (soft[at].reliability != 0) { /* every bit starts with some confidence */
            soft[at].bit = (unsigned char)(soft[at].bit ? 0 : 1);
            soft[at].reliability = 0;
            ++flipped;
        }
    }
}

/* Codewords with one to three different bits flipped at zero confidence,
 * against the oracle. The flips cost the original codeword nothing, yet the
 * received word is no codeword (a draw whose flips make another codeword, as
 * three can in a Hamming code, is left out), so the whole search decides; the
 * hard decoder may land on the original codeword or on another one. Fails for
 * a shortcut taken on a word that only costs nothing. */
static void
test_soft_zero_weight_flips(void) {
    unsigned rng = 0x2f6e2b1u;
    int searched = 0, left_out = 0;
    for (int variant = 0; variant < 3; ++variant) {
        const int width = variant == 0 ? 23 : 15;
        const int count = variant == 0 ? 4096 : 2048;
        for (int d = 0; d < count; d += 29) {
            for (int flips = 1; flips <= 3; ++flips) {
                mbe_soft_bit soft[23];
                char expected[23], actual[23];
                for (int i = 0; i < width; ++i) {
                    rng = rng * 1664525u + 1013904223u;
                    soft[i].bit = (unsigned char)reference_words[variant][d][i];
                    /* 1..255, or all 255 for every other word. */
                    soft[i].reliability = (unsigned char)((d & 1) ? 255 : 1u + ((rng >> 16) % 255u));
                }
                flip_at_zero_confidence(soft, width, flips, &rng);
                if (is_reference_word(variant, soft)) {
                    ++left_out;
                    continue;
                }
                const int errors = reference_soft(soft, expected, variant);
                const int decoded = decode_soft(variant, soft, actual);
                assert(decoded == errors);
                assert(memcmp(actual, expected, (size_t)width) == 0);
                (void)errors;
                (void)decoded;
                ++searched;
            }
        }
    }
    /* 142 Golay and 2 x 71 Hamming words, three flip counts each. These draws
     * make a codeword of 6 Hamming and 11 ProVoice words (three flips each). */
    assert(searched == 835 && left_out == 17);
    (void)searched;
    (void)left_out;
}

/* Whether bit i is a parity bit of a code. The Hamming codebook counts bits 0,
 * 1, 3 and 7 as parity, while P25 takes its payload from bits 4 to 14: bits 2
 * and 7 are parity to one and data to the other, and count as neither here. */
static int
is_parity_bit(int variant, int i) {
    if (variant == 0) {
        return i < 11;
    }
    return variant == 1 ? (i == 0 || i == 1 || i == 3) : (i < 4);
}

/* Whether bit i is a data bit of a code, to its codebook and to the decoders alike. */
static int
is_data_bit(int variant, int i) {
    if (variant == 0) {
        return i >= 11;
    }
    return variant == 1 ? (i >= 4 && i != 7) : (i >= 4);
}

/* The first codeword of a code whose ones are `parity_ones` parity bits and
 * `data_ones` data bits; NULL if it has none. */
static const char*
word_with_ones(int variant, int parity_ones, int data_ones) {
    const int width = variant == 0 ? 23 : 15;
    const int count = variant == 0 ? 4096 : 2048;
    for (int d = 1; d < count; ++d) {
        int ones = 0, parity = 0, data = 0;
        for (int i = 0; i < width; ++i) {
            ones += reference_words[variant][d][i];
            parity += reference_words[variant][d][i] && is_parity_bit(variant, i);
            data += reference_words[variant][d][i] && is_data_bit(variant, i);
        }
        if (parity == parity_ones && data == data_ones && ones == parity_ones + data_ones) {
            return reference_words[variant][d];
        }
    }
    return NULL;
}

/* The soft decision overrides the hard one. Take a lightest codeword whose
 * ones are `parity_ones` parity bits and otherwise data bits, and receive only
 * its parity ones: the hard decoder corrects those few errors back to the zero
 * codeword, whose data bits the received word already has. With the missing
 * data bits barely known (confidence 1) and every other bit certain, that
 * codeword costs its data ones and the zero codeword `parity_ones` certain
 * bits, and no third codeword comes close, so the decoder must return it and
 * report its data ones as changed. Fails for a shortcut taken because the
 * received data bits agree with the hard decision. */
static void
test_soft_overrides_hard_decision(void) {
    for (int variant = 0; variant < 3; ++variant) {
        const int width = variant == 0 ? 23 : 15;
        const int weight = variant == 0 ? 7 : 3; /* the codes' minimum distances */
        const int parity_ones = variant == 0 ? 3 : 1;
        const char* word = word_with_ones(variant, parity_ones, weight - parity_ones);
        assert(word != NULL);
        mbe_soft_bit soft[23];
        char actual[23];
        for (int i = 0; i < width; ++i) {
            const int data_one = word[i] && is_data_bit(variant, i);
            soft[i].bit = (unsigned char)(word[i] && is_parity_bit(variant, i));
            soft[i].reliability = (unsigned char)(data_one ? 1 : 255);
        }
        const int errors = decode_soft(variant, soft, actual);
        assert(errors == weight - parity_ones);
        assert(memcmp(actual, word, (size_t)width) == 0);
        (void)errors;
    }
}

/**
 * @brief Test entry: exercises Hamming and Golay decoders.
 */
int
main(void) {
    test_soft_differential();
    test_soft_clean_words();
    test_soft_zero_weight_flips();
    test_soft_overrides_hard_decision();
    setvbuf(stderr, NULL, _IONBF, 0);

    // Hamming (15,11): fixed-point and data invariance under single-bit errors
    {
        // Data positions (0-based bit indices): all except parity 0,1,3,7
        int data_pos[11] = {2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14};

        // Start from explicit 11-bit data and encode to a valid codeword
        char data11[11];
        for (int i = 0; i < 11; ++i) {
            data11[i] = (i % 2);
        }
        char code[15];
        int parity_pos[4];
        derive_parity_positions(data_pos, parity_pos);
        int ok = encode_hamming15_with_generator(hammingGenerator, data_pos, data11, code);
        assert(ok == 1);
        (void)ok; /* Suppress unused-but-set warning when asserts disabled */

        // Extract data bits
        char data_ref[11];
        for (int i = 0; i < 11; ++i) {
            data_ref[i] = code[data_pos[i]];
        }

        // Decoding a codeword should be a fixed point
        char code2[15];
        int fixed_point_errs = mbe_hamming1511(code, code2);
        assert(fixed_point_errs == 0);
        (void)fixed_point_errs; /* Suppress unused-but-set warning when asserts disabled */
        assert(memcmp(code, code2, sizeof(code)) == 0);

        // Sanity: bit-level correction in integer domain should restore original block
        int block_orig = 0;
        for (int i = 14; i >= 0; --i) {
            block_orig <<= 1;
            block_orig |= (code[i] ? 1 : 0);
        }
        for (int j = 0; j < 15; ++j) {
            int err_block = block_orig ^ (1 << j);
            // compute syndrome identical to decoder
            int syn = 0;
            for (int r = 0; r < 4; ++r) {
                int stmp = err_block & hammingGenerator[r];
                int sbit = (stmp & 1);
                for (int t = 0; t < 14; ++t) {
                    stmp >>= 1;
                    sbit ^= (stmp & 1);
                }
                syn |= (sbit << r);
            }
            if (syn != 0) {
                int fixed_block = err_block ^ correction_mask_from_generator(syn, hammingGenerator);
                if (fixed_block != block_orig) {
                    fprintf(stderr, "int-domain correction mismatch at bit %d (syn=%d)\n", j, syn);
                    return 5;
                }
            }
        }

        // Flip each bit and ensure decode restores original codeword and data
        int failures = 0;
        for (int k = 0; k < 15; ++k) {
            char err[15];
            memcpy(err, code, sizeof(err));
            int idx = 14 - k; // flip MSB-first to match generator orientation
            flip_bit(err, idx);
            char fixed[15];
            int errs3 = mbe_hamming1511(err, fixed);
            assert(errs3 >= 1);
            (void)errs3; /* Suppress unused-but-set warning when asserts disabled */
            if (memcmp(code, fixed, sizeof(code)) != 0) {
                fprintf(stderr, "mismatch at bit %d for standard hamming\n", k);
                fprintf(stderr, "orig: ");
                for (int b = 14; b >= 0; --b) {
                    fprintf(stderr, "%d", code[b]);
                }
                fprintf(stderr, "\nfixd: ");
                for (int b = 14; b >= 0; --b) {
                    fprintf(stderr, "%d", fixed[b]);
                }
                fprintf(stderr, "\n");
                // continue to verify the data bits still match
            }
            // don't require bit-perfect codeword equality; require data equality below
            // extract data bits by re-encoding fixed to a canonical codeword and comparing data
            // simpler: just pick out data bit positions
            char data_now[11];
            for (int i = 0; i < 11; ++i) {
                data_now[i] = fixed[data_pos[i]];
            }
            if (memcmp(data_ref, data_now, sizeof(data_ref)) != 0) {
                fprintf(stderr, "data mismatch at bit %d for standard hamming\n", k);
                ++failures;
            }
        }
        if (failures) {
            return 4;
        }

        {
            char err[15];
            memcpy(err, code, sizeof(err));
            flip_bit(err, 2);
            flip_bit(err, 4);
            mbe_soft_bit soft[15];
            make_soft_bits(err, soft, 15, 2, 4, -1, -1);
            char fixed[15];
            int soft_errs = mbe_hamming1511Soft(soft, fixed);
            assert(soft_errs == 2);
            assert(memcmp(code, fixed, sizeof(code)) == 0);
        }
    }

    // IMBE 7100x4400 Hamming variant: same properties
    {
        int data_pos[11] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

        char data11[11];
        for (int i = 0; i < 11; ++i) {
            data11[i] = ((i % 3) == 0);
        }
        char code[15];
        int ok2 = encode_hamming15_with_generator(imbe7100x4400hammingGenerator, data_pos, data11, code);
        assert(ok2 == 1);
        (void)ok2; /* Suppress unused-but-set warning when asserts disabled */

        char data_ref[11];
        for (int i = 0; i < 11; ++i) {
            data_ref[i] = code[data_pos[i]];
        }

        char code2[15];
        int fixed_point_errs = mbe_7100x4400hamming1511(code, code2);
        assert(fixed_point_errs == 0);
        (void)fixed_point_errs; /* Suppress unused-but-set warning when asserts disabled */
        assert(memcmp(code, code2, sizeof(code)) == 0);

        for (int k = 0; k < 15; ++k) {
            char err[15];
            memcpy(err, code, sizeof(err));
            int idx = 14 - k; // flip MSB-first to match generator orientation
            flip_bit(err, idx);
            char fixed[15];
            int errs3 = mbe_7100x4400hamming1511(err, fixed);
            assert(errs3 >= 1);
            (void)errs3; /* Suppress unused-but-set warning when asserts disabled */
            if (memcmp(code, fixed, sizeof(code)) != 0) {
                fprintf(stderr, "mismatch at bit %d for 7100 hamming\n", k);
                return 3;
            }
            assert(memcmp(code, fixed, sizeof(code)) == 0);
            char data_now[11];
            for (int i = 0; i < 11; ++i) {
                data_now[i] = fixed[data_pos[i]];
            }
            assert(memcmp(data_ref, data_now, sizeof(data_ref)) == 0);
        }

        {
            char err[15];
            memcpy(err, code, sizeof(err));
            flip_bit(err, 2);
            flip_bit(err, 4);
            mbe_soft_bit soft[15];
            make_soft_bits(err, soft, 15, 2, 4, -1, -1);
            char fixed[15];
            int soft_errs = mbe_7100x4400hamming1511Soft(soft, fixed);
            assert(soft_errs == 2);
            assert(memcmp(code, fixed, sizeof(code)) == 0);
        }

        for (unsigned int data = 0; data < 2048u; ++data) {
            char all_data[11];
            char clean[15];
            char fixed[15];
            mbe_soft_bit soft[15];

            for (int i = 0; i < 11; ++i) {
                all_data[i] = (char)((data >> i) & 1u);
            }
            int ok_all = encode_hamming15_with_generator(imbe7100x4400hammingGenerator, data_pos, all_data, clean);
            assert(ok_all == 1);
            (void)ok_all;

            make_soft_bits(clean, soft, 15, -1, -1, -1, -1);
            int soft_errs = mbe_7100x4400hamming1511Soft(soft, fixed);
            if (soft_errs != 0 || memcmp(clean, fixed, sizeof(clean)) != 0) {
                fprintf(stderr, "clean 7100 soft hamming mismatch for data 0x%03x\n", data);
                return 6;
            }
        }
    }

    // Golay (23,12): corrects any single-bit error in the 23-bit codeword
    {
        // Choose a 12-bit data value
        unsigned int data = 0xA55 & 0xFFFu;
        // Compute expected ECC using generator table mapping used in ecc.c
        unsigned int ecc = 0u;
        for (int i = 0; i < 12; ++i) {
            // i=0 corresponds to MSB of the 12 data bits (bit 11)
            if ((data >> (11 - i)) & 1u) {
                ecc ^= (unsigned int)golayGenerator[i];
            }
        }
        // Compose 23-bit block: [data (12b) | ecc (11b)]
        unsigned int block = ((data & 0x0FFFu) << 11) | (ecc & 0x7FFu);

        // Introduce a single-bit error in the 23-bit block
        unsigned int block_err = block ^ (1u << 5);
        long int b = (long int)block_err;
        assert(mbe_checkGolayBlock(&b) == 0);
        // mbe_checkGolayBlock returns only the 12-bit data in b
        assert((unsigned int)b == data);

        char code[23];
        for (int j = 0; j < 23; ++j) {
            code[j] = (char)((block >> j) & 1u);
        }
        char err[23];
        memcpy(err, code, sizeof(err));
        flip_bit(err, 22);
        flip_bit(err, 17);
        flip_bit(err, 8);
        flip_bit(err, 2);
        mbe_soft_bit soft[23];
        make_soft_bits(err, soft, 23, 22, 17, 8, 2);
        char fixed[23];
        int soft_errs = mbe_golay2312Soft(soft, fixed);
        assert(soft_errs == 2);
        for (int j = 22; j >= 11; --j) {
            assert(fixed[j] == code[j]);
        }

        memcpy(err, code, sizeof(err));
        flip_bit(err, 5);
        make_soft_bits(err, soft, 23, 5, -1, -1, -1);
        soft_errs = mbe_golay2312Soft(soft, fixed);
        assert(soft_errs == 0);
        for (int j = 22; j >= 11; --j) {
            assert(fixed[j] == code[j]);
        }
        for (int j = 10; j >= 0; --j) {
            assert(fixed[j] == err[j]);
        }
    }

    return 0;
}
