// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 *
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Spectral amplitude quantization and FEC shared by the AMBE encoders.
 *
 * Every step inverts the decoders' dequantization with their own tables and
 * float evaluation order (mbe_decodeAmbe2400Parms, mbe_decodeAmbe2450Parms),
 * so the residual is formed against exactly the prediction the decoder makes.
 */

#include "ambe_encoder.h"

#include <math.h>
#include <string.h>

#include "ambe3600x2400_internal.h"
#include "mbe_ecc.h"
#include "mbe_encoder.h"
#include "mbe_validation.h"

void
ambe_enc_mean_amplitude(struct ambe_enc_frame* q) {
    q->mean_a = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        q->mean_a += q->a[l];
    }
    q->mean_a /= (float)q->L;
}

/* The decoder sets the mean log2 magnitude to the transmitted gain, which
 * cannot go below the first gain step; see mbe_encoder_fit_floor(). */
static void
ambe_enc_fit_floor_gain(struct ambe_enc_frame* q) {
    const float floor_mean = q->gamma_q - (0.5f * log2f((float)q->L));
    q->mean_a = mbe_encoder_fit_floor(q->a, q->L, q->mean_a, floor_mean);
}

void
ambe_enc_quantize_gain(struct ambe_enc_frame* q, const struct ambe_enc_tables* t, const mbe_parms* prev_mp) {
    float gamma_raw = q->mean_a + (0.5f * log2f((float)q->L));
    float target = gamma_raw - (0.5f * prev_mp->gamma);
    int best = 0;
    float best_err = 1e30f;

    for (int c = 0; c < t->dg_count; c++) {
        float err = fabsf(t->dg[c] - target);
        if (err < best_err) {
            best_err = err;
            best = c;
        }
    }
    q->b[2] = best;
    q->gamma_q = t->dg[q->b[2]] + (0.5f * prev_mp->gamma);
    if (t->fit_floor_gain && target < t->dg[0]) {
        ambe_enc_fit_floor_gain(q);
    }
}

void
ambe_enc_prediction(struct ambe_enc_frame* q, const struct ambe_enc_tables* t, const mbe_parms* prev_mp) {
    int prev_L = mbe_clamp_harmonic_count(prev_mp->L);
    q->L = mbe_clamp_harmonic_count(q->L);
    float prev_log2Ml[57];
    memcpy(prev_log2Ml, prev_mp->log2Ml, sizeof(prev_log2Ml));
    /* The decoder's edge values (TIA-102.BABA-1 eqs 14-15). */
    for (int l = prev_L + 1; l <= q->L; l++) {
        prev_log2Ml[l] = prev_log2Ml[prev_L];
    }
    prev_log2Ml[0] = prev_log2Ml[1];

    q->mean_p = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        float flokl = ambe2400_prediction_position(prev_L, q->L, l);
        int intkl = (int)flokl;
        float deltal = flokl - (float)intkl;
        int upper = intkl + 1;
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }
        q->p[l] = ambe2400_interpolate_prediction(deltal, prev_log2Ml[intkl], prev_log2Ml[upper]);
        q->mean_p += q->p[l];
    }
    q->mean_p /= (float)q->L;

    /* eq 13 less its mean, which the gain carries. */
    for (int l = 1; l <= q->L; l++) {
        q->Tl[l] = q->a[l] - q->mean_a - (t->rho * (q->p[l] - q->mean_p));
    }
}

/* eq 16: the DCT of each of the four blocks. */
static void
ambe_enc_block_dct(struct ambe_enc_frame* q, const struct ambe_enc_tables* t) {
    q->Ji[1] = t->lmprbl[q->L][0];
    q->Ji[2] = t->lmprbl[q->L][1];
    q->Ji[3] = t->lmprbl[q->L][2];
    q->Ji[4] = t->lmprbl[q->L][3];

    int l = 1;
    for (int blk = 1; blk <= 4; blk++) {
        int ji = q->Ji[blk];
        for (int k = 1; k <= ji; k++) {
            float sum = 0.0f;
            for (int j = 1; j <= ji; j++) {
                sum += q->Tl[l + j - 1] * q->cache->idct_cos[ji][j][k];
            }
            /* Decoder IDCT is Tl[j] = sum_k a_k Cik[k] cos(theta_kj),
             * a_1=1, a_k=2. Exact inverse: Cik[k] = (1/ji) sum_j Tl[j] cos. */
            q->Cik[blk][k] = (1.0f / (float)ji) * sum;
        }
        l += ji;
    }
}

/* eqs 17-25: the PRBA vector and its 8-point DCT (Gm[1], the mean, is not sent). */
static void
ambe_enc_prba_dct(struct ambe_enc_frame* q) {
    float Ri[9];
    const float sqrt2 = 1.41421356237f;
    Ri[1] = q->Cik[1][1] + (sqrt2 * q->Cik[1][2]);
    Ri[2] = q->Cik[1][1] - (sqrt2 * q->Cik[1][2]);
    Ri[3] = q->Cik[2][1] + (sqrt2 * q->Cik[2][2]);
    Ri[4] = q->Cik[2][1] - (sqrt2 * q->Cik[2][2]);
    Ri[5] = q->Cik[3][1] + (sqrt2 * q->Cik[3][2]);
    Ri[6] = q->Cik[3][1] - (sqrt2 * q->Cik[3][2]);
    Ri[7] = q->Cik[4][1] + (sqrt2 * q->Cik[4][2]);
    Ri[8] = q->Cik[4][1] - (sqrt2 * q->Cik[4][2]);

    for (int m = 2; m <= 8; m++) {
        float sum = 0.0f;
        for (int i = 1; i <= 8; i++) {
            sum += Ri[i] * q->cache->ri_cos[m][i];
        }
        q->Gm[m] = sum / 8.0f;
    }
}

/* Index of the nearest of `rows` vectors (every `step`th) in their first
 * `used` elements, from element `first` of the target; ties keep the lower. */
static int
ambe_enc_nearest(const float* table, int width, int rows, int step, const float* target, int used) {
    int best = 0;
    float best_err = 1e30f;
    for (int c = 0; c < rows; c += step) {
        float err = 0.0f;
        for (int m = 0; m < used; m++) {
            float d = target[m] - table[(c * width) + m];
            err += d * d;
        }
        if (err < best_err) {
            best_err = err;
            best = c;
        }
    }
    return best;
}

void
ambe_enc_quantize_residual(struct ambe_enc_frame* q, const struct ambe_enc_tables* t) {
    ambe_enc_block_dct(q, t);
    ambe_enc_prba_dct(q);
    /* 4.3.1: b3 from Gm[2..4], b4 from Gm[5..8]. */
    q->b[3] = ambe_enc_nearest(&t->prba24[0][0], 3, 512, 1, &q->Gm[2], 3);
    q->b[4] = ambe_enc_nearest(&t->prba58[0][0], 4, 128, 1, &q->Gm[5], 4);
    /* 4.3.2: the HOC vector Cik[3..min(Ji, 6)]; an empty one quantizes to 0. */
    for (int blk = 0; blk < 4; blk++) {
        int ji = q->Ji[blk + 1];
        int used = ((ji < 6) ? ji : 6) - 2;
        q->b[5 + blk] = (used > 0) ? ambe_enc_nearest(&t->hoc[blk][0][0], 4, t->hoc_rows[blk], t->hoc_step[blk],
                                                      &q->Cik[blk + 1][3], used)
                                   : 0;
    }
}

static void
ambe_enc_c0(const char* ambe_d, char ambe_fr[4][24]) {
    char cw[23], data12[12];
    /* C0: 12 data bits (ambe_d[0..11]) -> (23,12) Golay + even parity */
    for (int i = 0; i < 12; i++) {
        data12[i] = ambe_d[i];
    }
    mbe_golay2312_encode(data12, cw);
    for (int j = 0; j < 23; j++) {
        ambe_fr[0][j + 1] = cw[j];
    }
    int ones = 0;
    for (int j = 0; j < 24; j++) {
        ones += (ambe_fr[0][j] & 1);
    }
    ambe_fr[0][0] = (char)(ones & 1);
}

static void
ambe_enc_pr(const char ambe_fr[4][24], unsigned short pr[25]) {
    unsigned short foo = 0;
    /* PR scramble sequence from the C0 data word */
    for (int i = 23; i >= 12; i--) {
        foo <<= 1;
        foo |= (unsigned short)(ambe_fr[0][i] & 1);
    }
    pr[0] = (unsigned short)(16 * foo);
    for (int i = 1; i < 25; i++) {
        pr[i] = (unsigned short)((173 * pr[i - 1]) + 13849 - (65536 * (((173 * pr[i - 1]) + 13849) / 65536)));
    }
}

/* C1: 12 data bits (ambe_d[12..23]) -> Golay, XOR p bits 23..1. Returns the
 * codeword's even parity scrambled with p bit 0. */
static int
ambe_enc_c1(const char* ambe_d, char ambe_fr[4][24], const unsigned short pr[25]) {
    char cw[23], data12[12];
    for (int i = 0; i < 12; i++) {
        data12[i] = ambe_d[12 + i];
    }
    mbe_golay2312_encode(data12, cw);
    for (int j = 0; j < 23; j++) {
        /* fr[1][j] = b word bit (j+1), scrambled with p bit (j+1) = pr[24-(j+1)] >> 15 */
        int prbit = pr[23 - j] / 32768;
        ambe_fr[1][j] = (char)((cw[j] & 1) ^ prbit);
    }

    int ones = 0;
    for (int j = 0; j < 23; j++) {
        ones += (cw[j] & 1);
    }
    return ((ones & 1) ^ (pr[24] / 32768)) & 1;
}

static void
ambe_enc_raw(const char* ambe_d, char ambe_fr[4][24]) {
    /* C2/C3: raw bits (decoder reads them MSB-first) */
    for (int j = 0; j < 11; j++) {
        ambe_fr[2][j] = ambe_d[24 + (10 - j)];
    }
    for (int j = 0; j < 14; j++) {
        ambe_fr[3][j] = ambe_d[35 + (13 - j)];
    }
}

int
ambe_enc_fec(const char ambe_d[49], char ambe_fr[4][24]) {
    unsigned short pr[25];
    char d[49];
    memcpy(d, ambe_d, sizeof(d)); /* the caller's bits may share the frame's memory */
    memset(ambe_fr, 0, 4 * sizeof(ambe_fr[0]));
    ambe_enc_c0(d, ambe_fr);
    ambe_enc_pr((const char (*)[24])ambe_fr, pr);
    ambe_enc_raw(d, ambe_fr);
    return ambe_enc_c1(d, ambe_fr, pr);
}
