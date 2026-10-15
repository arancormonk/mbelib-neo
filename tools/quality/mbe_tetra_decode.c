// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Decode a TETRA speech frame file (.CHD layout) to raw 16-bit PCM.
 *
 * Input: records of 138 little-endian 16-bit words, the bad-frame indicator followed by B1..B137,
 * as in EN 300 395-2 clause 8 and the EN 300 395-4 reference files (TE1_46.CHD). A nonzero
 * indicator marks a bad frame, whose payload words are not inspected; on a good frame every
 * payload word must be 0 or 1. Output: 240 little-endian 16-bit samples per frame at 8 kHz.
 *
 *     mbe_tetra_decode INPUT.chd OUTPUT.pcm
 *
 * The whole input is checked and decoded before the output is opened, so a malformed input
 * leaves an existing output untouched.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "mbe_quality_fs.h"
#include "mbelib-neo/mbelib.h"

#define RECORD_WORDS (1 + MBE_TETRA_FRAME_BITS)
#define RECORD_BYTES (2 * RECORD_WORDS)

static int
read_all(const char* path, unsigned char** data, size_t* size) {
    FILE* in = fopen(path, "rb");
    unsigned char* buf = NULL;
    size_t cap = 0;
    size_t len = 0;
    if (!in) {
        fprintf(stderr, "Cannot open input: %s\n", path);
        return -1;
    }
    for (;;) {
        size_t got;
        if (len == cap) {
            size_t ncap = cap ? cap * 2 : 65536;
            unsigned char* nbuf = (cap > SIZE_MAX / 2) ? NULL : (unsigned char*)realloc(buf, ncap);
            if (!nbuf) {
                fprintf(stderr, "Out of memory reading %s\n", path);
                free(buf);
                fclose(in);
                return -1;
            }
            buf = nbuf;
            cap = ncap;
        }
        got = fread(buf + len, 1, cap - len, in);
        len += got;
        if (got == 0) {
            break;
        }
    }
    if (ferror(in)) {
        fprintf(stderr, "Error reading %s\n", path);
        free(buf);
        fclose(in);
        return -1;
    }
    fclose(in);
    *data = buf;
    *size = len;
    return 0;
}

static int
word_at(const unsigned char* rec, size_t i) {
    unsigned v = (unsigned)rec[2 * i] | ((unsigned)rec[2 * i + 1] << 8);
    return (v >= 0x8000u) ? (int)v - 0x10000 : (int)v;
}

/** Decode every record into pcm (frames * 240 samples). Returns 0 or -1. */
static int
decode_all(const unsigned char* data, size_t frames, int16_t* pcm) {
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    size_t f;
    if (!dec) {
        fprintf(stderr, "Cannot allocate the decoder\n");
        return -1;
    }
    for (f = 0; f < frames; f++) {
        const unsigned char* rec = data + f * RECORD_BYTES;
        uint8_t bits[MBE_TETRA_FRAME_BITS];
        int bad = word_at(rec, 0) != 0;
        int i;
        for (i = 0; i < MBE_TETRA_FRAME_BITS && !bad; i++) {
            int w = word_at(rec, (size_t)i + 1);
            if (w != 0 && w != 1) {
                fprintf(stderr, "Frame %zu: bit B%d is %d, not 0 or 1\n", f + 1, i + 1, w);
                mbe_tetraDecoderFree(dec);
                return -1;
            }
            bits[i] = (uint8_t)w;
        }
        if (mbe_tetraDecodeFrame(dec, bad ? NULL : bits, bad, &pcm[f * MBE_TETRA_FRAME_SAMPLES]) != 0) {
            fprintf(stderr, "Frame %zu: decoding failed\n", f + 1);
            mbe_tetraDecoderFree(dec);
            return -1;
        }
    }
    mbe_tetraDecoderFree(dec);
    return 0;
}

static int
write_pcm(const char* path, const int16_t* pcm, size_t samples) {
    FILE* out = mbe_quality_open_output(path);
    size_t i;
    int ok = 1;
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", path);
        return -1;
    }
    for (i = 0; i < samples && ok; i++) {
        unsigned v = (unsigned)(uint16_t)pcm[i];
        unsigned char b[2] = {(unsigned char)(v & 0xFFu), (unsigned char)(v >> 8)};
        ok = fwrite(b, 1, 2, out) == 2;
    }
    if (fclose(out) != 0 || !ok) {
        fprintf(stderr, "Error writing %s\n", path);
        return -1;
    }
    return 0;
}

int
main(int argc, char** argv) {
    unsigned char* data = NULL;
    size_t size = 0;
    size_t frames;
    int16_t* pcm;
    int ret;
    if (argc != 3) {
        fprintf(stderr, "usage: %s INPUT.chd OUTPUT.pcm\n", argv[0]);
        return 2;
    }
    if (mbe_quality_paths_equal(argv[1], argv[2]) || mbe_quality_same_file(argv[1], argv[2])) {
        fprintf(stderr, "Input and output must be different files.\n");
        return 2;
    }
    if (read_all(argv[1], &data, &size) != 0) {
        return 2;
    }
    if (size % RECORD_BYTES != 0) {
        fprintf(stderr, "%s: %zu bytes is not a whole number of %d-byte records\n", argv[1], size, RECORD_BYTES);
        free(data);
        return 2;
    }
    frames = size / RECORD_BYTES;
    pcm = (int16_t*)calloc(frames ? frames * MBE_TETRA_FRAME_SAMPLES : 1, sizeof(int16_t));
    if (!pcm) {
        fprintf(stderr, "Out of memory\n");
        free(data);
        return 2;
    }
    ret = decode_all(data, frames, pcm) == 0 && write_pcm(argv[2], pcm, frames * MBE_TETRA_FRAME_SAMPLES) == 0 ? 0 : 2;
    free(pcm);
    free(data);
    return ret;
}
