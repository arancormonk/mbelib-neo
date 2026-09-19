// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief Encode an 8 kHz mono 16-bit PCM WAV into the private .dstar example container.
 *
 * Each 20 ms frame is a sync word (0x55 0x2D 0x16) plus 9 AMBE payload bytes.
 * This is not the D-STAR air-interface stream, which has a radio header and
 * slow-data/sync framing every 21 voice frames.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mbelib-neo/mbelib.h>

#include "example_file.h"

/* Optional in-speech noise reduction; DSTAR_DENOISE=0 disables it. */
#ifdef HAVE_SPECBLEACH
#include <specbleach_denoiser.h>
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

static int
read_wav_pcm(FILE* fp, uint32_t* data_size) {
    unsigned char hdr[12];
    if (fseek(fp, 0, SEEK_END) != 0) {
        return -1;
    }
    long file_size = ftell(fp);
    if (file_size < 12 || fseek(fp, 0, SEEK_SET) != 0 || fread(hdr, 1, sizeof(hdr), fp) != sizeof(hdr)) {
        return -1;
    }
    uint32_t remaining = read_le32(hdr + 4);
    if (memcmp(hdr, "RIFF", 4) != 0 || memcmp(hdr + 8, "WAVE", 4) != 0 || remaining < 4
        || (uint64_t)remaining + 8 > (uint64_t)file_size) {
        return -1;
    }
    remaining -= 4;
    bool have_fmt = false;
    while (remaining >= 8) {
        unsigned char chunk[8];
        if (fread(chunk, 1, sizeof(chunk), fp) != sizeof(chunk)) {
            return -1;
        }
        remaining -= 8;
        uint32_t size = read_le32(chunk + 4);
        if ((uint64_t)size + (size & 1u) > remaining) {
            return -1;
        }
        if (memcmp(chunk, "data", 4) == 0) {
            if (!have_fmt || (size & 1u) != 0) {
                return -1;
            }
            *data_size = size;
            return 0;
        }
        if (memcmp(chunk, "fmt ", 4) == 0) {
            if (read_wav_format(fp, size) < 0) {
                return -1;
            }
            have_fmt = true;
        } else if (skip_bytes(fp, size + (size & 1u)) < 0) {
            return -1;
        }
        remaining -= size + (size & 1u);
    }
    return -1;
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

#ifdef HAVE_SPECBLEACH
struct denoiser {
    void* handle;
    uint32_t latency_remaining;
};

static struct denoiser
init_denoiser(void) {
    struct denoiser nr = {0};
    const char* env = getenv("DSTAR_DENOISE");
    if (env != NULL && strcmp(env, "0") == 0) {
        return nr;
    }
    nr.handle = specbleach_initialize(8000, 20.0f);
    if (nr.handle != NULL) {
        SpectralBleachDenoiserParameters p = {0};
        p.reduction_amount = 12.0f;
        p.smoothing_factor = 30.0f;
        p.whitening_factor = 15.0f;
        p.adaptive_noise = 1;
        p.noise_estimation_method = 2; /* Martin Minimum Statistics */
        p.masking_depth = 0.5f;
        p.suppression_strength = 0.6f;
        specbleach_load_parameters(nr.handle, p);
        nr.latency_remaining = specbleach_get_latency(nr.handle);
    }
    return nr;
}

static bool
denoise_frame(struct denoiser* nr, short pcm[FRAME_SAMPLES]) {
    if (nr->handle == NULL) {
        return true;
    }
    float input[FRAME_SAMPLES], output[FRAME_SAMPLES];
    for (int i = 0; i < FRAME_SAMPLES; i++) {
        input[i] = (float)pcm[i] / 32768.0f;
    }
    specbleach_process(nr->handle, FRAME_SAMPLES, input, output);
    if (nr->latency_remaining >= FRAME_SAMPLES) {
        nr->latency_remaining -= FRAME_SAMPLES;
        return false;
    }
    nr->latency_remaining = 0;
    for (int i = 0; i < FRAME_SAMPLES; i++) {
        float sample = output[i];
        if (sample > 1.0f) {
            sample = 1.0f;
        }
        if (sample < -1.0f) {
            sample = -1.0f;
        }
        pcm[i] = (short)(sample * 32767.0f);
    }
    return true;
}
#endif

static int
write_frame(FILE* fp, const short pcm[FRAME_SAMPLES], mbe_parms* cur, mbe_parms* prev) {
    char bits[49], frame[4][24];
    unsigned char dv[DV_FRAME_BYTES] = {0x55, 0x2D, 0x16};
    if (mbe_encodeAmbe2400ParmsShort(pcm, bits, cur, prev) < 0 || mbe_encodeAmbe3600x2400Frame(bits, frame) < 0
        || mbe_encodeDStarDVData((const char (*)[24])frame, dv + 3) < 0) {
        (void)fprintf(stderr, "encode error\n");
        return -1;
    }
    if (fwrite(dv, 1, sizeof(dv), fp) != sizeof(dv)) {
        (void)fprintf(stderr, "write error: %s\n", strerror(errno));
        return -1;
    }
    mbe_moveMbeParms(cur, prev);
    return 0;
}

static int
encode_wav(FILE* fin, FILE* fout, uint32_t data_size) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
#ifdef HAVE_SPECBLEACH
    struct denoiser nr = init_denoiser();
#endif
    int ret = 0;
    while (data_size > 0) {
        short pcm[FRAME_SAMPLES];
        int n = read_samples(fin, pcm, &data_size);
        if (n < 0) {
            (void)fprintf(stderr, "truncated WAV data or read error\n");
            ret = 1;
            break;
        }
        if (n < FRAME_SAMPLES) {
            (void)fprintf(stderr, "warning: %d trailing samples ignored (not a full frame)\n", n);
            break;
        }
#ifdef HAVE_SPECBLEACH
        if (!denoise_frame(&nr, pcm)) {
            continue;
        }
#endif
        if (write_frame(fout, pcm, &cur, &prev) < 0) {
            ret = 1;
            break;
        }
    }
#ifdef HAVE_SPECBLEACH
    if (nr.handle != NULL) {
        specbleach_free(nr.handle);
    }
#endif
    return ret;
}

int
main(int argc, char** argv) {
    if (argc != 3) {
        (void)fprintf(stderr, "usage: %s input.wav output.dstar\n", argv[0]);
        return 1;
    }
    if (!example_paths_differ(argv[1], argv[2])) {
        return 1;
    }
    FILE* fin = example_open_file(argv[1], "rb");
    if (fin == NULL) {
        return 1;
    }
    FILE* fout = NULL;
    int ret = 1;
    uint32_t data_size;
    if (read_wav_pcm(fin, &data_size) < 0) {
        (void)fprintf(stderr, "invalid or truncated WAV\n");
        goto done;
    }
    fout = example_open_file(argv[2], "wb");
    if (fout != NULL) {
        ret = encode_wav(fin, fout, data_size);
    }
done:
    if (example_close_file(fin) < 0) {
        ret = 1;
    }
    if (example_close_file(fout) < 0) {
        ret = 1;
    }
    return ret;
}
