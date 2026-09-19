// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief D-STAR DV audio encoder sample app.
 *
 * Reads an 8 kHz mono 16-bit WAV and writes a .dstar file made of
 * 12-byte D-STAR DV frames (24-bit sync word 0x55 0x2D 0x16 + 72-bit
 * AMBE+FEC data), one per 20 ms, in air order.
 */

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mbelib-neo/mbelib.h>

/* Optional in-speech noise reduction front-end (libspecbleach). Disable
 * with DSTAR_DENOISE=0 if the input is already clean. */
#ifdef HAVE_SPECBLEACH
#include <specbleach_denoiser.h>
#endif

#define FRAME_SAMPLES 160
#define DV_FRAME_BYTES 12
#define DV_SYNC_0      0x55
#define DV_SYNC_1      0x2D
#define DV_SYNC_2      0x16

struct wav_header {
    uint32_t riff_tag;
    uint32_t riff_size;
    uint32_t wave_tag;
};

static int
read_samples(FILE* fp, short* buf, int n) {
    return (int)fread(buf, sizeof(short), (size_t)n, fp);
}

static int
read_wav_pcm(FILE* fp, int* rate, int* channels, uint32_t* data_size) {
    uint8_t tag[4];
    uint32_t size;
    uint16_t fmt_tag;
    uint16_t bits;
    struct wav_header hdr;

    if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        return -1;
    }
    if (hdr.riff_tag != 0x46464952U || hdr.wave_tag != 0x45564157U) { /* RIFF/WAVE */
        return -1;
    }

    for (;;) {
        if (fread(tag, 4, 1, fp) != 1) {
            return -1;
        }
        if (fread(&size, 4, 1, fp) != 1) {
            return -1;
        }

        if (memcmp(tag, "fmt ", 4) == 0) {
            if (fread(&fmt_tag, 2, 1, fp) != 1) {
                return -1;
            }
            if (fmt_tag != 1) {
                fprintf(stderr, "only uncompressed PCM WAV is supported\n");
                return -1;
            }
            if (fread(channels, 2, 1, fp) != 1) {
                return -1;
            }
            if (fread(rate, 4, 1, fp) != 1) {
                return -1;
            }
            fseek(fp, 6, SEEK_CUR); /* byte rate + block align */
            if (fread(&bits, 2, 1, fp) != 1) {
                return -1;
            }
            if (*rate != 8000 || *channels != 1 || bits != 16) {
                fprintf(stderr, "expected 8 kHz mono 16-bit WAV, got %d Hz %d ch %d bit\n", *rate, *channels, bits);
                return -1;
            }
            if (size > 16) {
                fseek(fp, (long)(size - 16), SEEK_CUR);
            }
        } else if (memcmp(tag, "data", 4) == 0) {
            *data_size = size;
            return 0;
        } else {
            if (size & 1) {
                size++;
            }
            fseek(fp, (long)size, SEEK_CUR);
        }
    }
}

int
main(int argc, char** argv) {
    FILE* fin = NULL;
    FILE* fout = NULL;
    short pcm[FRAME_SAMPLES];
    mbe_parms cur_mp;
    mbe_parms prev_mp;
    mbe_parms prev_mp_enhanced;
    uint32_t data_size;
    int rate;
    int channels;
    int ret = 1;
    bool denoise = (getenv("DSTAR_DENOISE") == NULL) || (strcmp(getenv("DSTAR_DENOISE"), "0") != 0);
    void* nr = NULL;
    uint32_t nr_latency = 0U;

    if (argc != 3) {
        fprintf(stderr, "usage: %s input.wav output.dstar\n", argv[0]);
        return 1;
    }

    fin = fopen(argv[1], "rb");
    if (fin == NULL) {
        fprintf(stderr, "cannot open %s: %s\n", argv[1], strerror(errno));
        return 1;
    }
    fout = fopen(argv[2], "wb");
    if (fout == NULL) {
        fprintf(stderr, "cannot open %s: %s\n", argv[2], strerror(errno));
        fclose(fin);
        return 1;
    }

    if (read_wav_pcm(fin, &rate, &channels, &data_size) < 0) {
        goto done;
    }
    (void)rate;
    (void)channels;

#ifdef HAVE_SPECBLEACH
    if (denoise) {
        nr = specbleach_initialize(8000, 20.0f);
        if (nr != NULL) {
            SpectralBleachDenoiserParameters p;
            memset(&p, 0, sizeof(p));
            p.learn_noise = 0;
            p.residual_listen = false;
            p.reduction_amount = 12.0f;
            p.smoothing_factor = 30.0f;
            p.whitening_factor = 15.0f;
            p.adaptive_noise = 1;
            p.noise_estimation_method = 2;   /* Martin Minimum Statistics */
            p.masking_depth = 0.5f;
            p.suppression_strength = 0.6f;
            p.aggressiveness = 0.0f;
            p.tonal_reduction = 0.0f;
            specbleach_load_parameters(nr, p);
            nr_latency = specbleach_get_latency(nr);
        } else {
            denoise = false;
        }
    }
#endif

    mbe_initMbeParms(&cur_mp, &prev_mp, &prev_mp_enhanced);

    {
        int n;
        /* Trim the denoiser latency off the front of the stream. */
        uint32_t denoised_total = 0U;
        while ((n = read_samples(fin, pcm, FRAME_SAMPLES)) == FRAME_SAMPLES) {
            char ambe_d[49];
            char frame_buf[4 * 24];
            unsigned char dv[DV_FRAME_BYTES];

            if (denoise) {
#ifdef HAVE_SPECBLEACH
                float finf[FRAME_SAMPLES];
                float foutf[FRAME_SAMPLES];
                for (int i = 0; i < FRAME_SAMPLES; i++)
                    finf[i] = (float)pcm[i] / 32768.0f;
                specbleach_process(nr, FRAME_SAMPLES, finf, foutf);
                denoised_total += FRAME_SAMPLES;
                if (denoised_total <= nr_latency) {
                    continue;   /* still inside the filter warm-up */
                }
                for (int i = 0; i < FRAME_SAMPLES; i++) {
                    float s = foutf[i];
                    if (s > 1.0f) s = 1.0f;
                    if (s < -1.0f) s = -1.0f;
                    pcm[i] = (short)(s * 32767.0f);
                }
#endif
            }

            mbe_encodeAmbe2400ParmsShort(pcm, ambe_d, &cur_mp, &prev_mp);
            mbe_encodeAmbe3600x2400Frame(ambe_d, (char(*)[24])frame_buf);

            dv[0] = DV_SYNC_0;
            dv[1] = DV_SYNC_1;
            dv[2] = DV_SYNC_2;
            mbe_encodeDStarDVData((const char(*)[24])frame_buf, dv + 3);

            if (fwrite(dv, 1, DV_FRAME_BYTES, fout) != DV_FRAME_BYTES) {
                fprintf(stderr, "write error: %s\n", strerror(errno));
                goto done;
            }

            mbe_moveMbeParms(&cur_mp, &prev_mp);
        }
        if (n > 0) {
            fprintf(stderr, "warning: %d trailing samples ignored (not a full frame)\n", n);
        }
    }

    ret = 0;

done:
#ifdef HAVE_SPECBLEACH
    if (nr != NULL)
        specbleach_free(nr);
#endif
    fclose(fin);
    fclose(fout);
    return ret;
}
