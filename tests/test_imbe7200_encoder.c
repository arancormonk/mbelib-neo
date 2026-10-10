// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/* IMBE 7200x4400 encoder: Hamming and frame FEC, the synchronization bit, the
 * quantizers against the decoder, and the codec-neutral checks of
 * encoder_test_support.h. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "encoder_test_support.h"
#include "imbe4400_internal.h"
#include "mbe_ecc.h"
#include "mbelib-neo/mbelib.h"

static void*
codec_alloc(void) {
    return mbe_imbe4400EncoderAlloc();
}

static void
codec_reset(void* enc) {
    mbe_imbe4400EncoderReset(enc);
}

static void
codec_release(void* enc) {
    mbe_imbe4400EncoderFree(enc);
}

static int
codec_encode(void* enc, const float* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeImbe4400Parms(enc, pcm, bits, cur, prev);
}

static int
codec_encode_short(void* enc, const short* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeImbe4400ParmsShort(enc, pcm, bits, cur, prev);
}

static int
codec_decode(const char* bits, mbe_parms* cur, mbe_parms* prev) {
    return mbe_decodeImbe4400Parms(bits, cur, prev);
}

static int
codec_process(float* out, mbe_process_result* result, const char* bits, mbe_parms* cur, mbe_parms* prev,
              mbe_parms* enhanced) {
    return mbe_processImbe4400Dataf(out, result, bits, cur, prev, enhanced);
}

/* Repeated, this voice frame takes the decoder's log2Ml to 34.8, about the
 * highest a search over repeated frames found. */
static const struct enc_codec imbe_codec = {
    "IMBE",
    88,
    -1,
    codec_alloc,
    codec_reset,
    codec_release,
    codec_encode,
    codec_encode_short,
    codec_decode,
    codec_process,
    "0100101111111111111111111100111111111110010000111101001001111111111111110101100001001011",
    NULL,
};

/* Every data word encodes to a codeword the decoder accepts unchanged, and
 * every single error is corrected. */
static int
test_hamming(void) {
    for (int word = 0; word < 2048; word++) {
        char in[11], cw[15], damaged[15], out[15];
        for (int i = 0; i < 11; i++) {
            in[i] = (char)((word >> (10 - i)) & 1);
        }
        mbe_hamming1511_encode(in, cw);
        if (mbe_hamming1511(cw, out) != 0 || memcmp(cw, out, 15) != 0) {
            printf("Hamming: word %d is not a codeword\n", word);
            return 1;
        }
        for (int i = 0; i < 11; i++) {
            if (out[14 - i] != in[i]) {
                return 1;
            }
        }
        for (int bit = 0; bit < 15; bit++) {
            memcpy(damaged, cw, sizeof(cw));
            damaged[bit] ^= 1;
            if (mbe_hamming1511(damaged, out) != 1 || memcmp(cw, out, 15) != 0) {
                printf("Hamming: word %d bit %d not corrected\n", word, bit);
                return 1;
            }
        }
    }
    puts("IMBE (15,11) Hamming: all 2048 words, every single error corrected");
    return 0;
}

static void
random_bits(char d[88]) {
    for (int i = 0; i < 88; i++) {
        d[i] = (char)((enc_rnd() >> 12) & 1); /* a high bit: the LCG's low bits repeat quickly */
    }
}

/* All 88 bits survive the frame FEC; cells past each vector are 0; one error
 * anywhere in a Golay or Hamming vector is corrected (the Golay decoder counts
 * the errors it corrects in data positions, the Hamming decoder any). */
static int
test_frame_fec(void) {
    static const int lengths[8] = {23, 23, 23, 23, 15, 15, 15, 7};
    for (int it = 0; it < 512; it++) {
        char d[88], fr[8][23], damaged[8][23], out[88];
        random_bits(d);
        if (mbe_encodeImbe7200x4400Frame(d, fr) != 0
            || mbe_decodeImbe7200x4400Frame((const char (*)[23])fr, out, NULL) != 0 || memcmp(d, out, 88) != 0) {
            printf("IMBE frame FEC: round trip failed at %d\n", it);
            return 1;
        }
        for (int v = 0; v < 8; v++) {
            for (int j = lengths[v]; j < 23; j++) {
                if (fr[v][j] != 0) {
                    return 1;
                }
            }
        }
        for (int v = 0; v < 7 && it < 64; v++) {
            for (int j = 0; j < lengths[v]; j++) {
                memcpy(damaged, fr, sizeof(fr));
                damaged[v][j] ^= 1;
                int counted = (v >= 4) || (j >= 11); /* Golay counts data errors, Hamming any */
                if (mbe_decodeImbe7200x4400Frame((const char (*)[23])damaged, out, NULL) != counted
                    || memcmp(d, out, 88) != 0) {
                    printf("IMBE frame FEC: error in vector %d bit %d not corrected\n", v, j);
                    return 1;
                }
            }
        }
    }
    puts("IMBE frame FEC: all 88 bits round-trip, every single error in u0..u6 corrected");
    return 0;
}

/* The frame may be written over the parameter bits it encodes. */
static int
test_frame_overlap(void) {
    for (int it = 0; it < 64; it++) {
        char d[88], expected[8][23], out[88];

        union {
            char bits[88];
            char frame[8][23];
        } shared;

        random_bits(d);
        memcpy(shared.bits, d, sizeof(shared.bits));
        if (mbe_encodeImbe7200x4400Frame(d, expected) != 0
            || mbe_encodeImbe7200x4400Frame(shared.bits, shared.frame) != 0
            || memcmp(expected, shared.frame, sizeof(expected)) != 0
            || mbe_decodeImbe7200x4400Frame((const char (*)[23])shared.frame, out, NULL) != 0
            || memcmp(d, out, 88) != 0) {
            puts("IMBE frame FEC: overlapping input and output differ");
            return 1;
        }
    }
    puts("IMBE frame FEC: output may overlap the input bits");
    return 0;
}

/* TIA-102.BABA 6.5: the synchronization bit is 0 after allocation and reset,
 * alternates with every encoded frame, and a rejected frame does not advance it. */
static int
test_sync_bit(void* enc) {
    mbe_parms cur, prev, enhanced;
    float pcm[160] = {0}, bad[160] = {0};
    char d[88];
    bad[3] = NAN;
    for (int pass = 0; pass < 2; pass++) {
        mbe_imbe4400EncoderReset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int frame = 0; frame < 9; frame++) {
            if (frame == 4 && mbe_encodeImbe4400Parms(enc, bad, d, &cur, &prev) != MBE_STATUS_INVALID_ARGUMENT) {
                return 1;
            }
            if (mbe_encodeImbe4400Parms(enc, pcm, d, &cur, &prev) != 0 || d[87] != (frame & 1)) {
                printf("IMBE sync bit: frame %d sent %d\n", frame, d[87]);
                return 1;
            }
            mbe_moveMbeParms(&cur, &prev);
        }
    }
    puts("IMBE synchronization bit: 0, 1, 0, ... from allocation and reset");
    return 0;
}

static int
frame_b0(const char d[88]) {
    static const int positions[8] = {0, 1, 2, 3, 4, 5, 85, 86};
    int b0 = 0;
    for (int i = 0; i < 8; i++) {
        b0 = (b0 << 1) | d[positions[i]];
    }
    return b0;
}

/*
 * Quantizing a decoded model again gives the same bits: random voice frames
 * are decoded and their log2 amplitudes and voicing requantized against the
 * same history. Every uniform quantizer's cell centre and the gain table's
 * levels map back to their own codes, so all 88 bits must match; this checks
 * the residual, the block and gain DCTs, the quantizers and the bit layout for
 * every reachable harmonic count independently of the analysis.
 */
static int
test_requantize(void) {
    int frames = 0;
    int counts[57] = {0};
    mbe_parms prev, cur, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int it = 0; it < 20000; it++) {
        char d[88], again[88];
        unsigned char bands[12] = {0};
        mbe_parms decoded = cur, history = prev;
        random_bits(d);
        if (frame_b0(d) > 207 || mbe_decodeImbe4400Parms(d, &decoded, &history) != 0) {
            continue;
        }
        for (int l = 1; l <= decoded.L; l++) {
            bands[((l <= 36) ? (l + 2) / 3 : 12) - 1] = (unsigned char)decoded.Vl[l];
        }
        mbe_imbe4400_quantize_amplitudes(frame_b0(d), bands, decoded.log2Ml, &prev, d[87], again);
        if (memcmp(d, again, 88) != 0) {
            int first = 0;
            while (d[first] == again[first]) {
                first++;
            }
            printf("IMBE requantize: frame %d (L %d) differs first at bit %d\n", it, decoded.L, first);
            return 1;
        }
        counts[decoded.L]++;
        frames++;
        prev = decoded;
    }
    for (int L = 9; L <= 56; L++) {
        if (counts[L] == 0) {
            printf("IMBE requantize: no frame with L %d\n", L);
            return 1;
        }
    }
    printf("IMBE requantization of %d decoded frames (L 9..56): bits identical\n", frames);
    return 0;
}

static int
test_invalid_arguments(void* enc) {
    float pcm[160] = {0};
    short shorts[160] = {0};
    char d[88] = {0}, fr[8][23] = {{0}};
    mbe_parms c, p, h;
    mbe_initMbeParms(&c, &p, &h);
    const int cases[] = {
        mbe_encodeImbe4400Parms(NULL, pcm, d, &c, &p),
        mbe_encodeImbe4400ParmsShort(NULL, shorts, d, &c, &p),
        mbe_encodeImbe4400Parms(enc, NULL, d, &c, &p),
        mbe_encodeImbe4400Parms(enc, pcm, NULL, &c, &p),
        mbe_encodeImbe4400Parms(enc, pcm, d, NULL, &p),
        mbe_encodeImbe4400Parms(enc, pcm, d, &c, NULL),
        mbe_encodeImbe4400Parms(enc, pcm, d, &p, &p),
        mbe_encodeImbe4400ParmsShort(enc, NULL, d, &c, &p),
        mbe_encodeImbe4400ParmsShort(enc, shorts, NULL, &c, &p),
        mbe_encodeImbe4400ParmsShort(enc, shorts, d, NULL, &p),
        mbe_encodeImbe4400ParmsShort(enc, shorts, d, &c, NULL),
        mbe_encodeImbe7200x4400Frame(NULL, fr),
        mbe_encodeImbe7200x4400Frame(d, NULL),
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        if (cases[i] != MBE_STATUS_INVALID_ARGUMENT) {
            printf("invalid arguments: case %zu returned %d\n", i, cases[i]);
            return 1;
        }
    }
    for (int i = 0; i < 88; i++) {
        d[i] = 2;
        if (mbe_encodeImbe7200x4400Frame(d, fr) != MBE_STATUS_INVALID_BITS) {
            return 1;
        }
        d[i] = 0;
    }
    puts("IMBE invalid arguments and bits rejected");
    return 0;
}

int
main(void) {
    void* enc = mbe_imbe4400EncoderAlloc();
    if (enc == NULL) {
        return 1;
    }
    int fails = 0;
    fails += test_invalid_arguments(enc);
    fails += test_hamming();
    fails += test_frame_fec();
    fails += test_frame_overlap();
    fails += test_sync_bit(enc);
    fails += test_requantize();
    fails += enc_run_common(&imbe_codec, enc);
    mbe_imbe4400EncoderFree(enc);
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}
