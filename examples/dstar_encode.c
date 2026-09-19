// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief Encode an 8 kHz mono 16-bit PCM WAV into the private .dstar example container.
 *
 * Usage: dstar_encode < input.wav > output.dstar
 *
 * Each 20 ms frame is a sync word (0x55 0x2D 0x16) plus 9 AMBE payload bytes.
 * This is not the D-STAR air-interface stream, which has a radio header and
 * slow-data/sync framing every 21 voice frames.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <mbelib-neo/mbelib.h>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#define FRAME_SAMPLES  160
#define DV_FRAME_BYTES 12

static uint16_t
read_le16(const unsigned char* p) {
    return (uint16_t)((uint16_t)p[0] | ((uint16_t)p[1] << 8));
}

static uint32_t
read_le32(const unsigned char* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static int
skip_bytes(FILE* fp, uint32_t count) {
    unsigned char buf[256];
    while (count > 0) {
        size_t n = count < sizeof(buf) ? (size_t)count : sizeof(buf);
        if (fread(buf, 1, n, fp) != n) {
            return -1;
        }
        count -= (uint32_t)n;
    }
    return 0;
}

static int
read_wav_format(FILE* fp, uint32_t size) {
    unsigned char fmt[16];
    if (size < sizeof(fmt) || fread(fmt, 1, sizeof(fmt), fp) != sizeof(fmt)) {
        return -1;
    }
    if (read_le16(fmt) != 1 || read_le16(fmt + 2) != 1 || read_le32(fmt + 4) != 8000 || read_le32(fmt + 8) != 16000
        || read_le16(fmt + 12) != 2 || read_le16(fmt + 14) != 16) {
        (void)fprintf(stderr, "expected uncompressed 8 kHz mono 16-bit PCM WAV\n");
        return -1;
    }
    return skip_bytes(fp, size - 16 + (size & 1u));
}

/* Account for the entire padded chunk within the declared RIFF boundary. */
static int
read_chunk_header(FILE* fp, uint32_t* remaining, unsigned char chunk[8], uint32_t* size) {
    if (*remaining < 8 || fread(chunk, 1, 8, fp) != 8) {
        return -1;
    }
    *remaining -= 8;
    *size = read_le32(chunk + 4);
    if ((uint64_t)*size + (*size & 1u) > *remaining) {
        return -1;
    }
    *remaining -= *size + (*size & 1u);
    return 0;
}

/* Return 1 at PCM data, 0 after another chunk, or -1 on invalid/truncated input. */
static int
read_wav_chunk(FILE* fp, uint32_t* remaining, bool* have_fmt, uint32_t* data_size) {
    unsigned char chunk[8];
    uint32_t size;
    if (read_chunk_header(fp, remaining, chunk, &size) < 0) {
        return -1;
    }
    if (memcmp(chunk, "data", 4) == 0) {
        if (!*have_fmt || (size & 1u) != 0) {
            return -1;
        }
        *data_size = size;
        return 1;
    }
    if (memcmp(chunk, "fmt ", 4) == 0) {
        if (read_wav_format(fp, size) < 0) {
            return -1;
        }
        *have_fmt = true;
    } else if (skip_bytes(fp, size + (size & 1u)) < 0) {
        return -1;
    }
    return 0;
}

static int
read_wav_pcm(FILE* fp, uint32_t* data_size, uint32_t* trailing_size) {
    unsigned char hdr[12];
    if (fread(hdr, 1, sizeof(hdr), fp) != sizeof(hdr)) {
        return -1;
    }
    uint32_t remaining = read_le32(hdr + 4);
    if (memcmp(hdr, "RIFF", 4) != 0 || memcmp(hdr + 8, "WAVE", 4) != 0 || remaining < 4) {
        return -1;
    }
    remaining -= 4;
    bool have_fmt = false;
    while (remaining >= 8) {
        int result = read_wav_chunk(fp, &remaining, &have_fmt, data_size);
        if (result != 0) {
            *trailing_size = remaining;
            return result > 0 ? 0 : -1;
        }
    }
    return -1;
}

static int
read_wav_tail(FILE* fp, uint32_t remaining) {
    while (remaining > 0) {
        unsigned char chunk[8];
        uint32_t size;
        if (read_chunk_header(fp, &remaining, chunk, &size) < 0 || skip_bytes(fp, size + (size & 1u)) < 0) {
            return -1;
        }
    }
    return 0;
}

static int
read_samples(FILE* fp, short pcm[FRAME_SAMPLES], uint32_t* remaining) {
    unsigned char bytes[FRAME_SAMPLES * 2];
    size_t count = *remaining < sizeof(bytes) ? (size_t)*remaining : sizeof(bytes);
    if (fread(bytes, 1, count, fp) != count) {
        return -1;
    }
    *remaining -= (uint32_t)count;
    for (size_t i = 0; i < count / 2; i++) {
        int value = read_le16(bytes + 2 * i);
        pcm[i] = (short)(value >= 32768 ? value - 65536 : value);
    }
    return (int)(count / 2);
}

static int
write_frame(mbe_ambe2400_encoder* enc, FILE* fp, const short pcm[FRAME_SAMPLES], mbe_parms* cur, mbe_parms* prev) {
    char bits[49], frame[4][24];
    unsigned char dv[DV_FRAME_BYTES] = {0x55, 0x2D, 0x16};
    if (mbe_encodeAmbe2400ParmsShort(enc, pcm, bits, cur, prev) < 0 || mbe_encodeAmbe3600x2400Frame(bits, frame) < 0
        || mbe_encodeDStarDVData((const char (*)[24])frame, dv + 3) < 0) {
        (void)fprintf(stderr, "encode error\n");
        return -1;
    }
    if (fwrite(dv, 1, sizeof(dv), fp) != sizeof(dv)) {
        perror("write error");
        return -1;
    }
    mbe_moveMbeParms(cur, prev);
    return 0;
}

static int
encode_wav(FILE* fin, FILE* fout, uint32_t data_size) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    mbe_ambe2400_encoder* enc = mbe_ambe2400EncoderAlloc();
    if (enc == NULL) {
        (void)fprintf(stderr, "cannot allocate encoder context\n");
        return 1;
    }
    int ret = 0;
    /* Emit exactly one zero frame after the last full input frame. */
    int flushed = 0;
    while (!flushed) {
        short pcm[FRAME_SAMPLES] = {0};
        if (data_size > 0) {
            int n = read_samples(fin, pcm, &data_size);
            if (n < 0) {
                (void)fprintf(stderr, "truncated WAV data or read error\n");
                ret = 1;
                break;
            }
            if (n < FRAME_SAMPLES) {
                (void)fprintf(stderr, "warning: %d trailing samples ignored (not a full frame)\n", n);
                continue;
            }
        } else {
            flushed = 1;
        }
        if (write_frame(enc, fout, pcm, &cur, &prev) < 0) {
            ret = 1;
            break;
        }
    }
    mbe_ambe2400EncoderFree(enc);
    return ret;
}

int
main(int argc, char** argv) {
    if (argc != 1) {
        (void)fprintf(stderr, "usage: %s < input.wav > output.dstar\n", argv[0]);
        return 1;
    }
#ifdef _WIN32
    if (_setmode(_fileno(stdin), _O_BINARY) == -1 || _setmode(_fileno(stdout), _O_BINARY) == -1) {
        perror("cannot set binary stdin/stdout");
        return 1;
    }
#endif
    uint32_t data_size, trailing_size;
    if (read_wav_pcm(stdin, &data_size, &trailing_size) < 0) {
        (void)fprintf(stderr, "invalid or truncated WAV\n");
        return 1;
    }
    int ret = encode_wav(stdin, stdout, data_size);
    if (ret == 0 && read_wav_tail(stdin, trailing_size) < 0) {
        (void)fprintf(stderr, "invalid or truncated trailing WAV chunk\n");
        ret = 1;
    }
    if (ferror(stdin) != 0) {
        perror("read error");
        ret = 1;
    }
    if (fflush(stdout) != 0 || ferror(stdout) != 0) {
        perror("write error");
        ret = 1;
    }
    return ret;
}
