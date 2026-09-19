// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief Decode the private .dstar example container to an 8 kHz mono PCM WAV.
 *
 * Usage: dstar_decode < input.dstar > output.wav
 *
 * Each 20 ms frame is a sync word (0x55 0x2D 0x16) plus 9 AMBE payload bytes.
 * This is not the D-STAR air-interface stream, which has a radio header and
 * slow-data/sync framing every 21 voice frames.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mbelib-neo/mbelib.h>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#define FRAME_SAMPLES     160
#define DV_FRAME_BYTES    12
#define MAX_WAV_DATA_SIZE (UINT32_MAX - 36u)

struct pcm_buffer {
    unsigned char* data;
    size_t size;
    size_t capacity;
};

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

static void
store_samples(unsigned char* bytes, short pcm[FRAME_SAMPLES], const char bits[49]) {
    /* Only b0 == 127 with tone index 128 (all eight tone bits zero) is silence. */
    int tone_bits = bits[6] | bits[7] | bits[8] | bits[9] | bits[10] | bits[11] | bits[42] | bits[43];
    if (bits[0] && bits[1] && bits[2] && bits[3] && bits[4] && bits[5] && bits[48] && tone_bits == 0) {
        memset(pcm, 0, FRAME_SAMPLES * sizeof(short));
    }
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
        bytes[2 * (size_t)i + 1] = (unsigned char)(sample >> 8);
    }
}

static int
append_samples(struct pcm_buffer* buffer, short pcm[FRAME_SAMPLES], const char bits[49]) {
    if (buffer->size > MAX_WAV_DATA_SIZE - FRAME_SAMPLES * 2u) {
        (void)fprintf(stderr, "WAV size limit exceeded\n");
        return -1;
    }
    size_t needed = buffer->size + (size_t)FRAME_SAMPLES * 2;
    if (needed > buffer->capacity) {
        size_t capacity = buffer->capacity;
        if (capacity > MAX_WAV_DATA_SIZE / 2u) {
            capacity = MAX_WAV_DATA_SIZE;
        } else {
            capacity *= 2;
        }
        unsigned char* data = (unsigned char*)realloc(buffer->data, capacity);
        if (data == NULL) {
            (void)fprintf(stderr, "out of memory buffering decoded PCM\n");
            return -1;
        }
        buffer->data = data;
        buffer->capacity = capacity;
    }
    store_samples(buffer->data + buffer->size, pcm, bits);
    buffer->size = needed;
    return 0;
}

static int
decode_frames(FILE* fin, struct pcm_buffer* buffer) {
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
        if (append_samples(buffer, pcm, bits) < 0) {
            return -1;
        }
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
    if (argc != 1) {
        (void)fprintf(stderr, "usage: %s < input.dstar > output.wav\n", argv[0]);
        return 1;
    }
#ifdef _WIN32
    if (_setmode(_fileno(stdin), _O_BINARY) == -1 || _setmode(_fileno(stdout), _O_BINARY) == -1) {
        perror("cannot set binary stdin/stdout");
        return 1;
    }
#endif
    struct pcm_buffer buffer = {NULL, 0, (size_t)FRAME_SAMPLES * 2};
    buffer.data = (unsigned char*)malloc(buffer.capacity);
    if (buffer.data == NULL) {
        (void)fprintf(stderr, "out of memory buffering decoded PCM\n");
        return 1;
    }
    int ret = 1;
    if (decode_frames(stdin, &buffer) == 0) {
        if (write_wav_header(stdout, (uint32_t)buffer.size) < 0
            || fwrite(buffer.data, 1, buffer.size, stdout) != buffer.size) {
            perror("WAV write error");
        } else {
            ret = 0;
        }
    }
    free(buffer.data);
    if (fflush(stdout) != 0 || ferror(stdout) != 0) {
        perror("write error");
        ret = 1;
    }
    return ret;
}
