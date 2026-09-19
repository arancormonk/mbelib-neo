// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief D-STAR DV audio decoder sample app.
 *
 * Reads a .dstar file of 12-byte D-STAR DV frames (sync word + 72-bit
 * AMBE+FEC data) as produced by dstar_encode, and writes an 8 kHz mono
 * 16-bit WAV.
 */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mbelib-neo/mbelib.h>

#define FRAME_SAMPLES 160
#define DV_FRAME_BYTES 12
#define DV_SYNC_0      0x55
#define DV_SYNC_1      0x2D
#define DV_SYNC_2      0x16

static int
write_wav_header(FILE* fp, uint32_t data_size) {
    uint32_t u32;
    uint16_t u16;

    u32 = 0x46464952U; /* RIFF */
    fwrite(&u32, 4, 1, fp);
    u32 = 36 + data_size;
    fwrite(&u32, 4, 1, fp);
    u32 = 0x45564157U; /* WAVE */
    fwrite(&u32, 4, 1, fp);

    u32 = 0x20746d66U; /* fmt  */
    fwrite(&u32, 4, 1, fp);
    u32 = 16;
    fwrite(&u32, 4, 1, fp);
    u16 = 1;
    fwrite(&u16, 2, 1, fp);
    u16 = 1;
    fwrite(&u16, 2, 1, fp);
    u32 = 8000;
    fwrite(&u32, 4, 1, fp);
    u32 = 16000;
    fwrite(&u32, 4, 1, fp);
    u16 = 2;
    fwrite(&u16, 2, 1, fp);
    u16 = 16;
    fwrite(&u16, 2, 1, fp);

    u32 = 0x61746164U; /* data */
    fwrite(&u32, 4, 1, fp);
    u32 = data_size;
    fwrite(&u32, 4, 1, fp);
    return 0;
}

int
main(int argc, char** argv) {
    FILE* fin = NULL;
    FILE* fout = NULL;
    unsigned char dv[DV_FRAME_BYTES];
    mbe_parms cur_mp;
    mbe_parms prev_mp;
    mbe_parms prev_mp_enhanced;
    int ret = 1;
    long frames = 0;
    long sync_errors = 0;

    if (argc != 3) {
        fprintf(stderr, "usage: %s input.dstar output.wav\n", argv[0]);
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

    /* Placeholder header; patched before closing. */
    write_wav_header(fout, 0);

    mbe_initMbeParms(&cur_mp, &prev_mp, &prev_mp_enhanced);

    while (fread(dv, 1, DV_FRAME_BYTES, fin) == DV_FRAME_BYTES) {
        short pcm[FRAME_SAMPLES];
        char frame_buf[4 * 24];
        char ambe_d[49];

        if (dv[0] != DV_SYNC_0 || dv[1] != DV_SYNC_1 || dv[2] != DV_SYNC_2) {
            sync_errors++;
        }

        mbe_decodeDStarDVData(dv + 3, (char(*)[24])frame_buf);
        mbe_processAmbe3600x2400Frame(pcm, NULL, (const char(*)[24])frame_buf, ambe_d,
                                      &cur_mp, &prev_mp, &prev_mp_enhanced);

        /* A silence frame (b0 == 127) would otherwise be rendered as loud
         * comfort noise; emit true silence instead. */
        if (ambe_d[0] && ambe_d[1] && ambe_d[2] && ambe_d[3] && ambe_d[4] && ambe_d[5] && ambe_d[48])
            memset(pcm, 0, sizeof(pcm));

        for (int i = 0; i < FRAME_SAMPLES; i++) {
            int v = (int)pcm[i];
            if (v > 31128)  v = 31128;
            if (v < -31128) v = -31128;
            pcm[i] = (short)v;
        }

        if (fwrite(pcm, sizeof(short), FRAME_SAMPLES, fout) != FRAME_SAMPLES) {
            fprintf(stderr, "write error: %s\n", strerror(errno));
            goto done;
        }

        mbe_moveMbeParms(&cur_mp, &prev_mp);
        frames++;
    }

    if (!feof(fin)) {
        fprintf(stderr, "read error: %s\n", strerror(errno));
        goto done;
    }

    ret = 0;

done:
    if (ret == 0) {
        long bytes = frames * FRAME_SAMPLES * (long)sizeof(short);
        fseek(fout, 0, SEEK_SET);
        write_wav_header(fout, (uint32_t)bytes);
        if (sync_errors > 0) {
            fprintf(stderr, "warning: %ld of %ld frames had sync word mismatches\n", sync_errors, frames);
        }
        fprintf(stderr, "decoded %ld frames (%ld ms)\n", frames, frames * 20);
    }

    fclose(fin);
    fclose(fout);
    return ret;
}
