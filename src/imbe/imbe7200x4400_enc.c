// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief IMBE 7200x4400 (P25 Phase 1 full rate) speech encoder.
 *
 * Quantizes the speech model of mbe_speech_analysis.c as TIA-102.BABA chapter
 * 6 describes (equation numbers below refer to it), against the tables
 * mbe_decodeImbe4400Parms() dequantizes with:
 *
 *  - b0 from the refined fundamental (eq 45); L and K as the decoder derives
 *    them from b0 (eqs 46-48), which equal the analysis' (eqs 31, 34).
 *  - b1 from the band V/UV decisions of 5.2 (eq 49).
 *  - Spectral amplitudes by the estimator of each band's decision (5.3, eqs 43
 *    and 44), then the prediction residual (eqs 52-55), six block DCTs and
 *    the gain vector (eqs 58-61), and the quantizers of eqs 62-64.
 *  - The synchronization bit of 6.5 (eq 80).
 *
 * A frame whose mean log2 amplitude lies below the gain table's lowest level
 * is flattened toward that level (mbe_encoder_fit_floor()), as D-STAR's are;
 * quantizing it as is would play a quiet tone tens of dB too loud.
 *
 * The prediction is the decoder's, which reads log2 M_0(-1) as log2 M_1(-1)
 * where eq 56 sets M_0(-1) = 1. DVSI's decoder predicts the same way: decoding
 * DVSI's P25 test vectors with eq 56 plays 0-250 Hz about 0.8 dB further from
 * DVSI's output on every development and validation vector. The encoder
 * predicts as the decoders it is heard through do.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "imbe4400_internal.h"
#include "imbe7200x4400_const.h"
#include "mbe_ecc.h"
#include "mbe_encoder.h"
#include "mbe_speech_analysis.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"

#define IMBE_ENC_MAX_B0    207 /* 6.1: b0 208..255 are not voice frames */

/*
 * Level: IMBE_ENC_MAG_SCALE scales the analysis amplitudes so the decoded
 * level matches DVSI's encoding of the same input, both decoded by this
 * library: the median, over DVSI's P25 development vectors, of the level that
 * puts our encoding at DVSI's (as AMBE2400_ENC_MAG_SCALE for D-STAR).
 */
#define IMBE_ENC_MAG_SCALE 1.06f

struct mbe_imbe4400_encoder {
    struct mbe_encoder_frontend fe;
    unsigned char bands_prev[MBE_ANALYSIS_BANDS]; /* previous frame's band decisions, eq 37 */
    int sync;                                     /* this frame's synchronization bit, eq 80 */
};

void
mbe_imbe4400EncoderReset(mbe_imbe4400_encoder* enc) {
    if (enc == NULL) {
        return;
    }
    mbe_encoder_frontend_reset(&enc->fe);
    memset(enc->bands_prev, 0, sizeof(enc->bands_prev));
    enc->sync = 0;
}

mbe_imbe4400_encoder*
mbe_imbe4400EncoderAlloc(void) {
    mbe_imbe4400_encoder* enc = calloc(1, sizeof(*enc));
    if (enc == NULL) {
        return NULL;
    }
    if (mbe_encoder_frontend_open(&enc->fe) != 0) {
        free(enc);
        return NULL;
    }
    mbe_imbe4400EncoderReset(enc);
    return enc;
}

void
mbe_imbe4400EncoderFree(mbe_imbe4400_encoder* enc) {
    if (enc != NULL) {
        mbe_encoder_frontend_close(&enc->fe);
        free(enc);
    }
}

/* One frame's quantizer values and workspace. */
struct imbe_enc_frame {
    int L;
    int K;
    int b[58];      /* b0..b(L+1) */
    float logm[57]; /* log2 M_l */
    float Tl[57];   /* prediction residual, eq 54 */
    float C[7][11]; /* block DCT coefficients C[i][k], eq 60 */
    float G[7];     /* transformed gain vector, eq 61 */
};

/* eq 45 for the refined fundamental f0 (cycles per sample), within 6.1's range. */
static int
imbe_enc_pitch(float f0) {
    long b0 = lround(floor((2.0 / (double)f0) - 39.0)); /* 4 pi / w0 with w0 = 2 pi f0 */
    return (int)((b0 < 0) ? 0 : ((b0 > IMBE_ENC_MAX_B0) ? IMBE_ENC_MAX_B0 : b0));
}

/* b0, and L and K as the decoder derives them (eqs 46-48); every voice code
 * 0..207 gives L 9..56. */
static void
imbe_enc_set_pitch(struct imbe_enc_frame* q, int b0) {
    q->b[0] = (b0 < 0) ? 0 : ((b0 > IMBE_ENC_MAX_B0) ? IMBE_ENC_MAX_B0 : b0);
    float w0 = (float)(4 * M_PI) / ((float)q->b[0] + 39.5f);
    int L = (int)(0.9254 * (int)((M_PI / w0) + 0.25));
    q->L = (L < 9) ? 9 : ((L > 56) ? 56 : L);
    q->K = (q->L < 37) ? (q->L + 2) / 3 : 12;
}

/* eq 55, as the decoder evaluates it. */
static float
imbe_enc_rho(int L) {
    if (L <= 15) {
        return 0.4f;
    }
    if (L <= 24) {
        return (0.03f * (float)L) - 0.05f;
    }
    return 0.7f;
}

/* eqs 52-54 against the decoder's previous log2 amplitudes, with its edge
 * values (log2 M_0(-1) = log2 M_1(-1), eq 57 above the previous L). */
static void
imbe_enc_residual(struct imbe_enc_frame* q, const mbe_parms* prev_mp) {
    const int prev_L = mbe_clamp_harmonic_count(prev_mp->L);
    const float rho = imbe_enc_rho(q->L);
    float prev[57];
    float interp[57];
    memcpy(prev, prev_mp->log2Ml, sizeof(prev));
    for (int l = prev_L + 1; l <= q->L; l++) {
        prev[l] = prev[prev_L];
    }
    prev[0] = prev[1];

    float sum = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        float flokl = ((float)prev_L / (float)q->L) * (float)l;
        int intkl = (int)flokl;
        float deltal = flokl - (float)intkl;
        int upper = (intkl + 1 > MBE_MAX_HARMONIC_BANDS) ? MBE_MAX_HARMONIC_BANDS : intkl + 1;
        interp[l] = ((1.0f - deltal) * prev[intkl]) + (deltal * prev[upper]);
        sum += interp[l];
    }
    sum *= rho / (float)q->L;
    for (int l = 1; l <= q->L; l++) {
        q->Tl[l] = q->logm[l] - (rho * interp[l]) + sum;
    }
}

/* eqs 58-61: the six block DCTs, Annex J block lengths, and the 6-point DCT
 * of their first coefficients. */
static void
imbe_enc_transform(struct imbe_enc_frame* q) {
    const int L9 = q->L - 9;
    int l = 1;
    for (int i = 1; i <= 6; i++) {
        const int J = ImbeJi[L9][i - 1];
        for (int k = 1; k <= J; k++) {
            double sum = 0.0;
            for (int j = 1; j <= J; j++) {
                sum += (double)q->Tl[l + j - 1] * cos(M_PI * (double)(k - 1) * ((double)j - 0.5) / (double)J);
            }
            q->C[i][k] = (float)(sum / (double)J);
        }
        l += J;
    }
    for (int m = 1; m <= 6; m++) {
        double sum = 0.0;
        for (int i = 1; i <= 6; i++) {
            sum += (double)q->C[i][1] * cos(M_PI * (double)(m - 1) * ((double)i - 0.5) / 6.0);
        }
        q->G[m] = (float)(sum / 6.0);
    }
}

/* eqs 62-63: a uniform quantizer of `bits` bits and step `step`, whose cells
 * the decoder reconstructs at step * (b - 2^(bits-1) + 0.5). */
static int
imbe_enc_uniform(float value, int bits, float step) {
    const long half = 1L << (bits - 1);
    double cell = floor((double)value / (double)step);
    if (cell < (double)-half) {
        return 0;
    }
    if (cell >= (double)half) {
        return (int)((2 * half) - 1);
    }
    return (int)((long)cell + half);
}

/* 6.3.1 and 6.3.2: b2 (Annex E), b3..b7 (Annex F) and b8..b(L+1) (Annex G),
 * in the order [C(1,2), ..., C(1,J1), ..., C(6,2), ..., C(6,J6)] of eq 64. */
static void
imbe_enc_quantize(struct imbe_enc_frame* q) {
    const int L9 = q->L - 9;
    int best = 0;
    for (int i = 1; i < 64; i++) {
        if (fabsf(B2[i] - q->G[1]) < fabsf(B2[best] - q->G[1])) {
            best = i;
        }
    }
    q->b[2] = best;
    for (int m = 3; m <= 7; m++) {
        q->b[m] = imbe_enc_uniform(q->G[m - 1], (int)ba[L9][m - 3][0], ba[L9][m - 3][1]);
    }
    int m = 8;
    for (int i = 1; i <= 6; i++) {
        for (int k = 2; k <= ImbeJi[L9][i - 1]; k++, m++) {
            int bits = hoba[L9][m - 8];
            q->b[m] = (bits > 0) ? imbe_enc_uniform(q->C[i][k], bits, quantstep[bits - 1] * standdev[k - 2]) : 0;
        }
    }
}

/* The 88 bits in the order mbe_decodeImbe4400Parms() reads them: b0 at 0..5
 * and 85..86, the bit order table for 6..84, and the synchronization bit. */
static void
imbe_enc_pack(const struct imbe_enc_frame* q, int sync, char imbe_d[88]) {
    static const int b0_positions[8] = {0, 1, 2, 3, 4, 5, 85, 86};
    const int L9 = q->L - 9;
    for (int i = 0; i < 8; i++) {
        imbe_d[b0_positions[i]] = (char)((q->b[0] >> (7 - i)) & 1);
    }
    for (int i = 6; i < 85; i++) {
        imbe_d[i] = (char)((q->b[bo[L9][i - 6][0]] >> bo[L9][i - 6][1]) & 1);
    }
    imbe_d[87] = (char)sync;
}

void
mbe_imbe4400_quantize_amplitudes(int b0, const unsigned char bands[12], const float logm[57], const mbe_parms* prev_mp,
                                 int sync, char imbe_d[88]) {
    struct imbe_enc_frame q;
    memset(&q, 0, sizeof(q));
    imbe_enc_set_pitch(&q, b0);
    for (int k = 1; k <= q.K; k++) {
        q.b[1] |= bands[k - 1] << (q.K - k); /* eq 49 */
    }
    memcpy(q.logm, logm, sizeof(q.logm));
    imbe_enc_residual(&q, prev_mp);
    imbe_enc_transform(&q);
    imbe_enc_quantize(&q);
    imbe_enc_pack(&q, sync, imbe_d);
}

static int
imbe_encode(mbe_imbe4400_encoder* enc, const float* samples, char imbe_d[88], mbe_parms* cur_mp,
            const mbe_parms* prev_mp) {
    struct mbe_analysis_result res;
    int status = mbe_encoder_frontend_analyze(&enc->fe, samples, &res);
    if (status < 0) {
        return status;
    }
    struct imbe_enc_frame q;
    unsigned char bands[MBE_ANALYSIS_BANDS];
    float logm[57] = {0};
    memset(&q, 0, sizeof(q));
    imbe_enc_set_pitch(&q, imbe_enc_pitch(res.f0));
    mbe_analysis_band_voicing(&res, q.L, enc->bands_prev, bands);
    float mean = 0.0f;
    for (int l = 1; l <= q.L; l++) {
        int voiced = bands[mbe_analysis_band_of(l) - 1];
        float m = IMBE_ENC_MAG_SCALE * (voiced ? res.voiced_magnitude[l] : res.noise_magnitude[l]);
        logm[l] = log2f(m + 1e-12f);
        mean += logm[l];
    }
    /* G1 (eq 61) carries the mean log2 amplitude, and the gain table's lowest
     * level is the floor under it. */
    mean /= (float)q.L;
    if (mean < B2[0]) {
        (void)mbe_encoder_fit_floor(logm, q.L, mean, B2[0]);
    }
    mbe_imbe4400_quantize_amplitudes(q.b[0], bands, logm, prev_mp, enc->sync, imbe_d);

    /* The decoder's update, on a private copy of the caller's history. */
    mbe_parms history = *prev_mp;
    status = mbe_decodeImbe4400Parms(imbe_d, cur_mp, &history);
    if (status != 0) {
        return (status < 0) ? status : MBE_STATUS_INVALID_ARGUMENT;
    }
    if (!mbe_encoder_model_valid(cur_mp)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    memcpy(enc->bands_prev, bands, sizeof(enc->bands_prev));
    enc->sync ^= 1;
    return 0;
}

int
mbe_encodeImbe4400Parms(mbe_imbe4400_encoder* enc, const float* samples, char imbe_d[88], mbe_parms* cur_mp,
                        const mbe_parms* prev_mp) {
    if (enc == NULL || samples == NULL || imbe_d == NULL || cur_mp == NULL || prev_mp == NULL || cur_mp == prev_mp
        || !mbe_encoder_samples_valid(samples) || !mbe_encoder_history_valid(prev_mp)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    /* A frame that fails (a history whose model would overflow) leaves no trace. */
    const struct mbe_analysis_state saved = enc->fe.analysis;
    mbe_parms out = *cur_mp;
    int status = imbe_encode(enc, samples, imbe_d, &out, prev_mp);
    if (status < 0) {
        enc->fe.analysis = saved;
        return status;
    }
    *cur_mp = out;
    return 0;
}

int
mbe_encodeImbe4400ParmsShort(mbe_imbe4400_encoder* enc, const short* samples, char imbe_d[88], mbe_parms* cur_mp,
                             const mbe_parms* prev_mp) {
    float float_buf[MBE_ENCODER_SAMPLES];

    if (samples == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    mbe_encoder_short_to_float(samples, float_buf);
    return mbe_encodeImbe4400Parms(enc, float_buf, imbe_d, cur_mp, prev_mp);
}

/* Chapter 7: u0..u3 (12 bits each) as (23,12) Golay codes, u4..u6 (11 bits
 * each) as (15,11) Hamming codes and u7 (7 bits) unprotected, then vectors
 * 1..6 modulated by the sequence seeded from u0 (eqs 84-94), as
 * mbe_demodulateImbe7200x4400Data() removes it. */
int
mbe_encodeImbe7200x4400Frame(const char imbe_d[88], char imbe_fr[8][23]) {
    unsigned short pr[115];
    unsigned short seed = 0;

    if (imbe_d == NULL || imbe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 88; i++) {
        if (imbe_d[i] != 0 && imbe_d[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }
    char d[88];
    memcpy(d, imbe_d, sizeof(d)); /* the caller's bits may share the frame's memory */
    memset(imbe_fr, 0, 8 * sizeof(imbe_fr[0]));
    for (int i = 0; i < 4; i++) {
        mbe_golay2312_encode(d + ((size_t)12 * (size_t)i), imbe_fr[i]);
    }
    for (int i = 0; i < 3; i++) {
        mbe_hamming1511_encode(d + 48 + ((size_t)11 * (size_t)i), imbe_fr[4 + i]);
    }
    for (int j = 0; j < 7; j++) {
        imbe_fr[7][6 - j] = d[81 + j];
    }

    for (int j = 22; j >= 11; j--) {
        seed = (unsigned short)((seed << 1) | (unsigned short)imbe_fr[0][j]);
    }
    pr[0] = (unsigned short)(16 * seed);
    for (int i = 1; i < 115; i++) {
        pr[i] = (unsigned short)((173u * pr[i - 1] + 13849u) & 0xFFFFu);
    }
    int k = 1;
    for (int i = 1; i < 7; i++) {
        for (int j = (i < 4) ? 22 : 14; j >= 0; j--) {
            imbe_fr[i][j] = (char)(imbe_fr[i][j] ^ (pr[k++] >> 15));
        }
    }
    return 0;
}
