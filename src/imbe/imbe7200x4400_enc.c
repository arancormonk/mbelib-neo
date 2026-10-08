// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 *
 * Based on Bruce Perens' ham_digital_modes
 * (https://github.com/BrucePerens/hams_open).
 */

/**
 * @file
 * @brief IMBE 7200x4400 (P25 Phase 1 full rate) speech encoder.
 *
 * C port of ham_digital_modes' float TIA-102.BABA encoder
 * (tia_102_baba/encoder.rs, mod.rs, prediction.rs, quantize.rs and
 * parameter_encoding.rs):
 *
 *  - Analysis (mbe_frame_analysis.c): pitch tracking and refinement, V/UV
 *    bands and spectral amplitudes at the refined fundamental (5.1-5.3).
 *  - b0 (eq 45), b1 from the band decisions, then the amplitude residual of
 *    eq 54 against the decoder's previous amplitudes, split into the Annex J
 *    blocks, block DCT (eq 60), gain vector DCT (eq 61) and the scalar
 *    quantizers of eq 62-64 (Annex E-G tables).
 *  - Bits are laid out with this library's bit order table, so
 *    mbe_decodeImbe4400Parms() reads them back exactly, and the prediction
 *    history is advanced by decoding each emitted frame with it.
 *
 * Deviation from the reference: the prediction follows this library's
 * decoder, which reads log2 M_0(-1) as log2 M_1(-1) where eq 56 has 0.
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "imbe7200x4400_const.h"
#include "mbe_ecc.h"
#include "mbe_frame_analysis.h"
#include "mbelib-neo/mbelib.h"

#define IMBE4400_ENC_SAMPLES   160
#define IMBE4400_ENC_PCM_SCALE 32768.0 /* analysis runs on the 16-bit scale */

struct mbe_imbe4400_encoder {
    struct mbe_fa_tables tables;
    struct mbe_fa_state analysis;
    struct mbe_fa_voicing_state voicing;
    struct mbe_fa_result frame;
    mbe_parms prev; /* the decoder's prediction history */
};

void
mbe_imbe4400EncoderReset(mbe_imbe4400_encoder* enc) {
    mbe_parms cur;
    mbe_parms prev_enhanced;

    if (enc == NULL) {
        return;
    }
    mbe_fa_reset(&enc->analysis);
    mbe_fa_voicing_reset(&enc->voicing);
    mbe_initMbeParms(&cur, &enc->prev, &prev_enhanced);
}

mbe_imbe4400_encoder*
mbe_imbe4400EncoderAlloc(void) {
    mbe_imbe4400_encoder* enc = calloc(1, sizeof(*enc));
    if (enc == NULL) {
        return NULL;
    }
    mbe_fa_init_tables(&enc->tables);
    mbe_imbe4400EncoderReset(enc);
    return enc;
}

void
mbe_imbe4400EncoderFree(mbe_imbe4400_encoder* enc) {
    free(enc);
}

/* Frames emitted before the first analysed frame are encoded as silent
 * input: a zero spectrum, the standard's initial pitch period of 100 samples,
 * and a pitch error of 1 (no periodicity). */
static void
imbe4400_enc_quiet_frame(struct mbe_fa_result* frame) {
    memset(frame, 0, sizeof(*frame));
    frame->initial_pitch = 100.0;
    frame->omega0_hat = 2.0 * M_PI / 100.0;
    frame->initial_pitch_error = 1.0;
}

/* Eq 55, in the decoder's float arithmetic. */
static float
imbe4400_enc_rho(int L) {
    if (L <= 15) {
        return 0.4f;
    }
    if (L <= 24) {
        return (0.03f * (float)L) - 0.05f;
    }
    return 0.7f;
}

/* The decoder's previous log2Ml: index 0 is read as index 1, and harmonics
 * past the previous frame's last hold its value (eq 57). */
static double
imbe4400_enc_prev_log2(const mbe_parms* prev, int prev_l, int j) {
    if (j < 1) {
        j = 1;
    }
    if (j > prev_l) {
        j = prev_l;
    }
    if (j > 56) {
        j = 56;
    }
    return (double)prev->log2Ml[j];
}

/* Eq 62-64: uniform quantizer of `bits` bits and step `step`, saturating. */
static int
imbe4400_enc_uniform(double value, int bits, double step) {
    int64_t half = (int64_t)1 << (bits - 1);
    int64_t idx = (int64_t)floor(value / step);
    if (idx < -half) {
        return 0;
    }
    if (idx >= half) {
        return (int)(((int64_t)1 << bits) - 1);
    }
    return (int)(idx + half);
}

/* Eq 52-54: the residual of log2 m[1..L] against the decoder's prediction. */
static void
imbe4400_enc_residual(int L, const double* m, const mbe_parms* prev, double residual[57]) {
    double interp[57];
    double sum77 = 0.0;
    double rho = (double)imbe4400_enc_rho(L);
    int prev_l = prev->L;
    if (prev_l < 1) {
        prev_l = 1;
    }
    if (prev_l > 56) {
        prev_l = 56;
    }
    for (int l = 1; l <= L; l++) {
        float flokl = ((float)prev_l / (float)L) * (float)l;
        int intkl = (int)flokl;
        double delta = (double)(flokl - (float)intkl);
        interp[l] = ((1.0 - delta) * imbe4400_enc_prev_log2(prev, prev_l, intkl))
                    + (delta * imbe4400_enc_prev_log2(prev, prev_l, intkl + 1));
        sum77 += interp[l];
    }
    sum77 *= rho / (double)L;
    for (int l = 1; l <= L; l++) {
        residual[l] = log2(m[l]) - (rho * interp[l]) + sum77;
    }
}

/* Split the residual into the Annex J blocks and apply the DCT to each
 * (eq 58-60); the gain vector is the DCT of the six block means (eq 61). */
static void
imbe4400_enc_transform(int L, const double residual[57], double dct[6][10], double g[7]) {
    const int L9 = L - 9;
    int offset = 1;
    for (int i = 0; i < 6; i++) {
        int j_len = ImbeJi[L9][i];
        for (int k = 1; k <= j_len; k++) {
            double sum = 0.0;
            for (int j = 1; j <= j_len; j++) {
                sum += residual[offset + j - 1] * cos(M_PI * ((double)k - 1.0) * ((double)j - 0.5) / (double)j_len);
            }
            dct[i][k - 1] = sum / (double)j_len;
        }
        offset += j_len;
    }
    for (int mm = 0; mm < 6; mm++) {
        double sum = 0.0;
        for (int i = 0; i < 6; i++) {
            sum += dct[i][0] * cos(M_PI * (double)mm * ((double)(i + 1) - 0.5) / 6.0);
        }
        g[mm + 1] = sum / 6.0;
    }
}

/* Eq 62-64: b2 (Annex E), b3..b7 (Annex F) and b8..b(L+1), the higher-order
 * coefficients block by block (Annex G). */
static void
imbe4400_enc_quantize_transform(int L, const double dct[6][10], const double g[7], int bvals[58]) {
    const int L9 = L - 9;
    int b2 = 0;
    for (int i = 1; i < 64; i++) {
        if (fabs((double)B2[i] - g[1]) < fabs((double)B2[b2] - g[1])) {
            b2 = i;
        }
    }
    bvals[2] = b2;
    for (int mm = 2; mm <= 6; mm++) {
        int bits = (int)ba[L9][mm - 2][0];
        double step = (double)ba[L9][mm - 2][1];
        bvals[mm + 1] = imbe4400_enc_uniform(g[mm], bits, step);
    }
    int pos = 8;
    for (int i = 0; i < 6; i++) {
        for (int k = 2; k <= ImbeJi[L9][i]; k++) {
            int bits = hoba[L9][pos - 8];
            double step = (bits > 0) ? (double)quantstep[bits - 1] * (double)standdev[k - 2] : 0.0;
            bvals[pos] = (bits > 0) ? imbe4400_enc_uniform(dct[i][k - 1], bits, step) : 0;
            pos++;
        }
    }
}

/* Quantize L amplitudes m[1..L] into b2..b(L+1) of bvals, against the
 * decoder's previous frame. */
static void
imbe4400_enc_quantize_amplitudes(int L, const double* m, const mbe_parms* prev, int bvals[58]) {
    double residual[57] = {0};
    double dct[6][10] = {{0}};
    double g[7] = {0};
    if (L < 9 || L > 56) {
        return;
    }
    imbe4400_enc_residual(L, m, prev, residual);
    imbe4400_enc_transform(L, residual, dct, g);
    imbe4400_enc_quantize_transform(L, (const double (*)[10])dct, g, bvals);
}

/* Lay b0..b(L+1) out as mbe_decodeImbe4400Parms() reads them. */
static void
imbe4400_enc_pack(int L, const int bvals[58], char imbe_d[88]) {
    const int L9 = L - 9;
    static const int b0_indices[8] = {0, 1, 2, 3, 4, 5, 85, 86};

    memset(imbe_d, 0, 88);
    if (L9 < 0 || L9 > 47) {
        return;
    }
    for (int i = 0; i < 8; i++) {
        imbe_d[b0_indices[i]] = (char)((bvals[0] >> (7 - i)) & 1);
    }
    for (int i = 6; i < 85; i++) {
        int param = bo[L9][i - 6][0];
        int bit = bo[L9][i - 6][1];
        imbe_d[i] = (char)((bvals[param] >> bit) & 1);
    }
}

static int
imbe4400_encode(mbe_imbe4400_encoder* enc, const double input[IMBE4400_ENC_SAMPLES], char imbe_d[88]) {
    struct mbe_fa_voicing_state scratch;
    struct mbe_fa_voicing_state* voicing = &enc->voicing;
    const struct mbe_fa_result* a = &enc->frame;
    unsigned char bands[MBE_FA_MAX_BANDS];
    double amplitudes[57];
    int bvals[58] = {0};

    if (!mbe_fa_push(&enc->tables, &enc->analysis, input, &enc->frame)) {
        /* While the look-ahead buffer is filling, a copy of the voicing state
         * is used, so that only analysed frames update the original. */
        imbe4400_enc_quiet_frame(&enc->frame);
        scratch = enc->voicing;
        voicing = &scratch;
    }

    double omega0 = a->omega0_hat;
    int L = mbe_fa_harmonics_count(omega0);
    if (L < 9 || L > 56) {
        return MBE_STATUS_INVALID_ARGUMENT; /* not reached: the refined pitch range limits L to 9..56 */
    }
    int k_hat = mbe_fa_determine_voicing(&enc->tables, &a->sw, omega0, a->initial_pitch_error, voicing, bands);
    mbe_fa_spectral_amplitudes(&enc->tables, &a->sw, L, k_hat, omega0, bands, amplitudes);
    /* The logarithm of a zero amplitude (digital silence) is undefined, so
     * amplitudes are floored at one 16-bit PCM step, as in OP25's
     * imbe_vocoder. */
    for (int l = 1; l <= L; l++) {
        amplitudes[l] = fmax(amplitudes[l], 1.0);
    }

    /* Eq 45: b0, the fundamental; eq 46 maps it back to the same L. */
    bvals[0] = (int)floor((4.0 * M_PI / omega0) - 39.0 + 1e-9);
    /* b1: band 1 in the most significant bit. */
    for (int k = 1; k <= k_hat; k++) {
        if (bands[k - 1]) {
            bvals[1] |= 1 << (k_hat - k);
        }
    }
    imbe4400_enc_quantize_amplitudes(L, amplitudes, &enc->prev, bvals);
    imbe4400_enc_pack(L, bvals, imbe_d);

    /* Update the prediction history by decoding the emitted frame, exactly as
     * the receiving decoder will. */
    mbe_parms cur = enc->prev;
    if (mbe_decodeImbe4400Parms(imbe_d, &cur, &enc->prev) == 0) {
        mbe_moveMbeParms(&cur, &enc->prev);
    }
    return 0;
}

/* Accept only finite samples of magnitude at most 2^20 (bit pattern
 * 0x49800000). The test inspects the IEEE 754 bit pattern rather than using
 * isfinite() or comparisons, because the library may be compiled with
 * MBELIB_ENABLE_FAST_MATH (-ffast-math, /fp:fast), under which the compiler
 * is permitted to assume that no value is NaN or infinite and may remove such
 * checks. */
static int
imbe4400_enc_samples_valid(const float* samples) {
    for (int i = 0; i < IMBE4400_ENC_SAMPLES; i++) {
        uint32_t bits;
        memcpy(&bits, &samples[i], sizeof(bits));
        if ((bits & 0x7FFFFFFFu) > 0x49800000u) {
            return 0;
        }
    }
    return 1;
}

int
mbe_encodeImbe4400Parms(mbe_imbe4400_encoder* enc, const float* samples, char imbe_d[88]) {
    double input[IMBE4400_ENC_SAMPLES];

    if (enc == NULL || samples == NULL || imbe_d == NULL || !imbe4400_enc_samples_valid(samples)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < IMBE4400_ENC_SAMPLES; i++) {
        input[i] = (double)samples[i] * IMBE4400_ENC_PCM_SCALE;
    }
    return imbe4400_encode(enc, input, imbe_d);
}

int
mbe_encodeImbe4400ParmsShort(mbe_imbe4400_encoder* enc, const short* samples, char imbe_d[88]) {
    double input[IMBE4400_ENC_SAMPLES];

    if (enc == NULL || samples == NULL || imbe_d == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < IMBE4400_ENC_SAMPLES; i++) {
        input[i] = (double)samples[i];
    }
    return imbe4400_encode(enc, input, imbe_d);
}

int
mbe_encodeImbe7200x4400Frame(const char imbe_d[88], char imbe_fr[8][23]) {
    char data[12];
    char cw[23];
    unsigned short pr[115];
    unsigned short foo = 0;

    if (imbe_d == NULL || imbe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 88; i++) {
        if (imbe_d[i] != 0 && imbe_d[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }
    memset(imbe_fr, 0, 8 * sizeof(imbe_fr[0]));

    /* Eq 81-83: u0..u3 as (23,12) Golay, u4..u6 as (15,11) Hamming, u7 raw. */
    for (int i = 0; i < 4; i++) {
        memcpy(data, imbe_d + ((size_t)12 * (size_t)i), 12);
        mbe_golay2312_encode(data, cw);
        memcpy(imbe_fr[i], cw, 23);
    }
    for (int i = 0; i < 3; i++) {
        mbe_hamming1511_encode(imbe_d + 48 + ((size_t)11 * (size_t)i), cw);
        memcpy(imbe_fr[4 + i], cw, 15);
    }
    for (int j = 0; j < 7; j++) {
        imbe_fr[7][6 - j] = imbe_d[81 + j];
    }

    /* Eq 84-94: modulate vectors 1-6 with the sequence seeded from u0. */
    for (int i = 22; i >= 11; i--) {
        foo = (unsigned short)((foo << 1) | (unsigned short)(imbe_fr[0][i] & 1));
    }
    pr[0] = (unsigned short)(16 * foo);
    for (int i = 1; i < 115; i++) {
        pr[i] = (unsigned short)(((173u * pr[i - 1]) + 13849u) & 0xFFFFu);
    }
    int k = 1;
    for (int i = 1; i < 4; i++) {
        for (int j = 22; j >= 0; j--) {
            imbe_fr[i][j] = (char)(imbe_fr[i][j] ^ (pr[k++] >> 15));
        }
    }
    for (int i = 4; i < 7; i++) {
        for (int j = 14; j >= 0; j--) {
            imbe_fr[i][j] = (char)(imbe_fr[i][j] ^ (pr[k++] >> 15));
        }
    }
    return 0;
}
