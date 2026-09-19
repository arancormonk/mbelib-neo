// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief Decode the private .dstar example container to an 8 kHz mono PCM WAV.
 *
 * Each 20 ms frame is a sync word (0x55 0x2D 0x16) plus 9 AMBE payload bytes.
 * This is not the D-STAR air-interface stream, which has a radio header and
 * slow-data/sync framing every 21 voice frames.
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <mbelib-neo/mbelib.h>

#include "example_file.h"

#define FRAME_SAMPLES  160
#define DV_FRAME_BYTES 12

static void
write_le32(unsigned char* bytes, uint32_t value) {
    for (int i = 0; i < 4; i++) {
        bytes[i] = (unsigned char)(value >> (8 * i));
    }
}

static int
write_wav_header(FILE* fp, uint32_t data_size) {
    unsigned char header[44] = {'R', 'I', 'F', 'F', 0,  0, 0,   0,   'W', 'A',  'V',  'E', 'f', 'm',  't',
                                ' ', 16,  0,   0,   0,  1, 0,   1,   0,   0x40, 0x1f, 0,   0,   0x80, 0x3e,
                                0,   0,   2,   0,   16, 0, 'd', 'a', 't', 'a',  0,    0,   0,   0};
    write_le32(header + 4, 36 + data_size);
    write_le32(header + 40, data_size);
    return fwrite(header, 1, sizeof(header), fp) == sizeof(header) ? 0 : -1;
}

static int
write_samples(FILE* fp, short pcm[FRAME_SAMPLES], const char bits[49]) {
    /* The example emits true silence for the codec's comfort-noise frame. */
    if (bits[0] && bits[1] && bits[2] && bits[3] && bits[4] && bits[5] && bits[48]) {
        memset(pcm, 0, FRAME_SAMPLES * sizeof(short));
    }
    unsigned char bytes[FRAME_SAMPLES * 2];
    for (int i = 0; i < FRAME_SAMPLES; i++) {
        int value = pcm[i];
        if (value > 31128) {
            value = 31128;
        }
        if (value < -31128) {
            value = -31128;
        }
        uint16_t sample = (uint16_t)value;
        bytes[2 * (size_t)i] = (unsigned char)sample;
        bytes[2 * i + 1] = (unsigned char)(sample >> 8);
    }
    return fwrite(bytes, 1, sizeof(bytes), fp) == sizeof(bytes) ? 0 : -1;
}

static int
decode_frames(FILE* fin, FILE* fout, uint32_t* data_size) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    unsigned char dv[DV_FRAME_BYTES];
    uint32_t frames = 0, sync_errors = 0;
    size_t n;
    while ((n = fread(dv, 1, sizeof(dv), fin)) == sizeof(dv)) {
        short pcm[FRAME_SAMPLES];
        char frame[4][24], bits[49];
        if (dv[0] != 0x55 || dv[1] != 0x2D || dv[2] != 0x16) {
            sync_errors++;
        }
        if (mbe_decodeDStarDVData(dv + 3, frame) < 0
            || mbe_processAmbe3600x2400Frame(pcm, NULL, (const char (*)[24])frame, bits, &cur, &prev, &enhanced) < 0) {
            (void)fprintf(stderr, "decode error\n");
            return -1;
        }
        if (*data_size > UINT32_MAX - 36u - FRAME_SAMPLES * 2u) {
            (void)fprintf(stderr, "WAV size limit exceeded\n");
            return -1;
        }
        if (write_samples(fout, pcm, bits) < 0) {
            (void)fprintf(stderr, "write error: %s\n", strerror(errno));
            return -1;
        }
        *data_size += FRAME_SAMPLES * 2u;
        frames++;
    }
    if (ferror(fin) != 0 || n != 0) {
        (void)fprintf(stderr, "read error or trailing partial .dstar frame\n");
        return -1;
    }
    if (sync_errors > 0) {
        (void)fprintf(stderr, "warning: %u of %u frames had sync word mismatches\n", sync_errors, frames);
    }
    (void)fprintf(stderr, "decoded %u frames\n", frames);
    return 0;
}

int
main(int argc, char** argv) {
    if (argc != 3) {
        (void)fprintf(stderr, "usage: %s input.dstar output.wav\n", argv[0]);
        return 1;
    }
    if (!example_paths_differ(argv[1], argv[2])) {
        return 1;
    }
    FILE* fin = example_open_file(argv[1], "rb");
    if (fin == NULL) {
        return 1;
    }
    FILE* fout = example_open_file(argv[2], "wb");
    int ret = 1;
    uint32_t data_size = 0;
    if (fout == NULL) {
        goto done;
    }
    if (write_wav_header(fout, 0) < 0 || decode_frames(fin, fout, &data_size) < 0) {
        goto done;
    }
    if (fseek(fout, 0, SEEK_SET) != 0 || write_wav_header(fout, data_size) < 0) {
        (void)fprintf(stderr, "WAV header write error: %s\n", strerror(errno));
        goto done;
    }
    ret = 0;
done:
    if (example_close_file(fin) < 0) {
        ret = 1;
    }
    if (example_close_file(fout) < 0) {
        ret = 1;
    }
    return ret;
}
