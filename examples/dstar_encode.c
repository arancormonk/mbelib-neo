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

/* Optional in-speech noise reduction; DSTAR_DENOISE=0 disables it. */
#ifdef HAVE_SPECBLEACH
#include <stdlib.h>

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

#ifdef HAVE_SPECBLEACH
struct denoiser {
    void* handle;
    uint32_t latency_remaining;
    size_t pending_samples;
    short pending_pcm[FRAME_SAMPLES];
};

static int
denoiser_disabled(void) {
#ifdef _MSC_VER
    char env[2];
    size_t length;
    if (getenv_s(&length, NULL, 0, "DSTAR_DENOISE") != 0 || length != sizeof(env)) {
        return 0;
    }
    return getenv_s(&length, env, sizeof(env), "DSTAR_DENOISE") == 0 && strcmp(env, "0") == 0;
#else
    const char* env = getenv("DSTAR_DENOISE");
    return env != NULL && strcmp(env, "0") == 0;
#endif
}

static struct denoiser
init_denoiser(void) {
    struct denoiser nr = {0};
    if (denoiser_disabled()) {
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
    bool ready = false;
    for (int i = 0; i < FRAME_SAMPLES; i++) {
        /* Drop the exact delay, retaining partial output until a full frame
         * is ready. Latency need not be a multiple of FRAME_SAMPLES. */
        if (nr->latency_remaining > 0) {
            nr->latency_remaining--;
            continue;
        }
        float sample = output[i];
        if (sample > 1.0f) {
            sample = 1.0f;
        }
        if (sample < -1.0f) {
            sample = -1.0f;
        }
        nr->pending_pcm[nr->pending_samples++] = (short)(sample * 32767.0f);
        if (nr->pending_samples == FRAME_SAMPLES) {
            memcpy(pcm, nr->pending_pcm, sizeof(nr->pending_pcm));
            nr->pending_samples = 0;
            ready = true;
        }
    }
    return ready;
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
#ifdef HAVE_SPECBLEACH
    struct denoiser nr = init_denoiser();
#endif
    int ret = 0;
    uint32_t frames_read = 0, frames_emitted = 0;
    /* Emit one flush frame beyond the full input frames, feeding extra zeros
     * through the same path until any denoiser latency has been drained. */
    while (data_size > 0 || frames_emitted < frames_read + 1) {
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
            frames_read++;
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
        frames_emitted++;
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
