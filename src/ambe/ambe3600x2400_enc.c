// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 *
 * Copyright (C) 2010 mbelib Author
 * GPG Key ID: 0xEA5EFE2C (9E7A 5527 9CDC EBF7 BF1B  D772 4F98 E863 EA5E FE2C)
 *
 * Portions were originally under the ISC license; this mbelib-neo
 * distribution is provided under GPL-2.0-or-later. See LICENSE for details.
 */

/**
 * @file
 * @brief AMBE 3600x2400 (D-STAR DV) speech encoder.
 *
 * Performs the inverse of mbe_decodeAmbe2400Parms():
 *
 *  - Speech analysis (mbe_speech_analysis.c) follows the method of
 *    ANSI/TIA-102.BABA chapter 5: pitch from the E(P) error function with
 *    tracking and quarter-sample refinement, per-harmonic fit errors against
 *    the window spectrum, and voicing-independent spectral magnitudes
 *    (US 5,701,390).
 *  - Voicing is decided per 500 Hz column as soft decisions against the
 *    standard's frequency- and energy-dependent thresholds (US 8,595,002),
 *    then reduced to the four 1 kHz bits D-STAR transmits by an
 *    energy-weighted choice.
 *  - Every parameter is quantized against the exact pitch law, harmonic
 *    count and tables the decoder dequantizes with (AmbePlusVuv/AmbePlusDg/
 *    AmbePlusPRBA24/AmbePlusPRBA58/AmbePlusHOCb5..b8), so the
 *    reconstructed frame is bit-compatible with this library's
 *    mbe_decodeAmbe2400Parms()/mbe_processAmbe3600x2400*() path and follows
 *    the D-STAR AMBE bit layout (interleave, scrambler and Golay parity
 *    cross-checked against the MMDVM tables). Reconstruction matches DVSI's
 *    AMBE-3000 D-STAR test vectors; on-air interoperability is not verified.
 *  - The prediction state (log2Ml, gamma) is advanced with the
 *    QUANTIZED values, mirroring the decoder state so the encoder and
 *    decoder never drift apart.
 *
 * Every frame is a voice frame, as from DVSI's encoder: quiet and silent
 * input is coded as low-level voice rather than as a silence frame, and the
 * decoded level follows the input level (there is no AGC). Tone (DTMF)
 * encoding is not supported.
 */

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "ambe3600x2400_const.h"
#include "ambe3600x2400_internal.h"
#include "mbe_ecc.h"
#include "mbe_speech_analysis.h"
#include "mbe_unvoiced_fft.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"

#define AMBE2400_ENC_SAMPLES   160
#define AMBE2400_ENC_PCM_SCALE 32768.0f /* analysis runs on the 16-bit scale */

/*
 * Level: AMBE2400_ENC_MAG_SCALE maps the analysis magnitudes (A/2 for a
 * sinusoid of amplitude A on the 16-bit scale) into the codec's log2 amplitude
 * domain. DVSI's encoder keeps the input level: across its D-STAR test vectors
 * (-2 to -63 dBFS active level) its bits decode at the input level within
 * about 1 dB. The scale is the median, over the development vectors, of the
 * level that puts our encoding at DVSI's encoding of the same input, both
 * decoded by this library; the validation vectors then sit within about
 * 1 dB of DVSI's.
 */
#define AMBE2400_ENC_MAG_SCALE 1.02f

/*
 * Caller-owned encoder state: one context per stream. The analysis tables
 * are computed once at allocation and kept across resets.
 */
struct mbe_ambe2400_encoder {
    struct mbe_analysis_state analysis;
    mbe_fft_plan* fft;
    mbe_acf_plan* acf;
    struct mbe_analysis_tables tables;
};

void
mbe_ambe2400EncoderReset(mbe_ambe2400_encoder* enc) {
    if (enc == NULL) {
        return;
    }
    mbe_analysis_reset(&enc->analysis);
}

mbe_ambe2400_encoder*
mbe_ambe2400EncoderAlloc(void) {
    mbe_ambe2400_encoder* enc = calloc(1, sizeof(*enc));
    if (enc == NULL) {
        return NULL;
    }
    enc->fft = mbe_fft_plan_alloc();
    enc->acf = mbe_acf_plan_alloc();
    if (enc->fft == NULL || enc->acf == NULL) {
        mbe_fft_plan_free(enc->fft);
        mbe_acf_plan_free(enc->acf);
        free(enc);
        return NULL;
    }
    mbe_analysis_init_tables(&enc->tables);
    mbe_ambe2400EncoderReset(enc);
    return enc;
}

void
mbe_ambe2400EncoderFree(mbe_ambe2400_encoder* enc) {
    if (enc != NULL) {
        mbe_fft_plan_free(enc->fft);
        mbe_acf_plan_free(enc->acf);
        free(enc);
    }
}

/* Temporary analysis and quantization workspace; no inter-frame state. */
struct ambe2400_enc_frame {
    const struct ambe_dct_cache* cache;
    float p[57];
    float a[57], Tl[57], Tl_q[57];
    float Cik[5][18], Cik_q[5][18], Gm[9];
    int Ji[5], L, b[9];
    float f0q, gamma_q, mean_a, mean_p;
};

static void
ambe2400_enc_quantize_pitch(struct ambe2400_enc_frame* q, float f0) {
    q->b[0] = (int)lroundf(((MBE_AMBE2400_LOG2_F0_OFFSET - log2f(f0)) / MBE_AMBE2400_LOG2_F0_SLOPE) - 0.5f);
    if (q->b[0] < 0) {
        q->b[0] = 0;
    }
    if (q->b[0] > 125) {
        q->b[0] = 125;
    }

    q->f0q = mbe_ambe2400_f0(q->b[0]);
    q->L = mbe_ambe2400_harmonic_count(q->f0q);
}

/* The decoder's 500 Hz column for harmonic l at the transmitted pitch. */
static int
ambe2400_enc_column(const struct ambe2400_enc_frame* q, int l) {
    int jl = (int)((float)l * 16.0f * q->f0q);
    return (jl > 7) ? 7 : jl;
}

/* US 8,595,002 eq [3]: soft voicing from the energy-weighted fit error. */
static float
ambe2400_enc_soft_voicing(float energy, float error, float threshold) {
    if (energy <= 0.0f || threshold <= 0.0f) {
        return 0.0f;
    }
    if (error <= 0.0f) {
        return 1.0f;
    }
    float lv = 0.5f * (1.0f - log2f(error / (threshold * energy)));
    if (lv < 0.0f) {
        return 0.0f;
    }
    return (lv > 1.0f) ? 1.0f : lv;
}

/*
 * V/UV: soft voicing per 500 Hz column, then the four 1 kHz bits D-STAR sends
 * (each covers two columns) by the energy-weighted squared error of
 * US 8,595,002. Columns follow the decoder's floor(16 l f0) assignment, so a
 * decision lands on the harmonics it was measured from; ties and empty bands
 * are unvoiced.
 */
static void
ambe2400_enc_voicing(struct ambe2400_enc_frame* q, const struct mbe_analysis_state* state,
                     const struct mbe_analysis_result* res, unsigned char columns[MBE_ANALYSIS_COLUMNS]) {
    float energy[MBE_ANALYSIS_COLUMNS] = {0};
    float error[MBE_ANALYSIS_COLUMNS] = {0};
    float lv[MBE_ANALYSIS_COLUMNS];
    for (int l = 1; l <= q->L && l <= res->harmonics; l++) {
        int column = ambe2400_enc_column(q, l);
        float e = res->magnitude[l] * res->magnitude[l];
        energy[column] += e;
        error[column] += res->fit_error[l] * e;
    }
    for (int k = 0; k < MBE_ANALYSIS_COLUMNS; k++) {
        lv[k] = ambe2400_enc_soft_voicing(energy[k], error[k], mbe_analysis_column_threshold(state, res, k + 1));
    }
    q->b[1] = 0;
    for (int bit = 0; bit < 4; bit++) {
        int c0 = 2 * bit;
        int c1 = c0 + 1;
        float voiced_cost =
            (energy[c0] * (1.0f - lv[c0]) * (1.0f - lv[c0])) + (energy[c1] * (1.0f - lv[c1]) * (1.0f - lv[c1]));
        float unvoiced_cost = (energy[c0] * lv[c0] * lv[c0]) + (energy[c1] * lv[c1] * lv[c1]);
        int voiced = (energy[c0] + energy[c1] > 0.0f) && (voiced_cost < unvoiced_cost);
        columns[c0] = (unsigned char)voiced;
        columns[c1] = (unsigned char)voiced;
        q->b[1] |= voiced << (3 - bit);
    }
}

static void
ambe2400_enc_magnitudes(struct ambe2400_enc_frame* q, const struct mbe_analysis_result* res) {
    q->mean_a = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        q->a[l] = log2f((res->magnitude[l] * AMBE2400_ENC_MAG_SCALE) + 1e-12f);
        q->mean_a += q->a[l];
    }
    q->mean_a /= (float)q->L;
}

/* Energy of the envelope mean + s * (a - mean_a), in the log2 amplitude domain. */
static float
ambe2400_enc_envelope_energy(const struct ambe2400_enc_frame* q, float mean, float s) {
    float energy = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        energy += exp2f(2.0f * (mean + (s * (q->a[l] - q->mean_a))));
    }
    return energy;
}

/*
 * The decoder sets the mean log2 magnitude to the transmitted gain, which
 * cannot go below AmbePlusDg[0]. A frame whose own mean lies below that (a
 * quiet tone: one strong harmonic over empty bands) would decode with the
 * excess pushed into its peak. Scale the deviations from the mean toward a
 * flat envelope until the energy at the floor gain matches the frame's; a
 * flat envelope is the quietest the floor gain can play.
 */
static void
ambe2400_enc_fit_floor_gain(struct ambe2400_enc_frame* q) {
    const float floor_mean = q->gamma_q - (0.5f * log2f((float)q->L));
    const float target = ambe2400_enc_envelope_energy(q, q->mean_a, 1.0f);
    float lo = 0.0f;
    float hi = 1.0f;
    if (ambe2400_enc_envelope_energy(q, floor_mean, 0.0f) >= target) {
        hi = 0.0f;
    }
    for (int i = 0; i < 30 && hi > 0.0f; i++) {
        float mid = 0.5f * (lo + hi);
        if (ambe2400_enc_envelope_energy(q, floor_mean, mid) > target) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    for (int l = 1; l <= q->L; l++) {
        q->a[l] = floor_mean + (hi * (q->a[l] - q->mean_a));
    }
    q->mean_a = floor_mean;
}

static void
ambe2400_enc_quantize_gain(struct ambe2400_enc_frame* q, const mbe_parms* prev_mp) {
    /* Gain quantization */
    {
        float gamma_raw = q->mean_a + (0.5f * log2f((float)q->L));
        float target = gamma_raw - (0.5f * prev_mp->gamma);
        int best = 0;
        float best_err = 1e30f;

        for (int c = 0; c < 64; c++) {
            float err = fabsf(AmbePlusDg[c] - target);
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        q->b[2] = best;
        q->gamma_q = AmbePlusDg[q->b[2]] + (0.5f * prev_mp->gamma);
        if (target < AmbePlusDg[0]) {
            ambe2400_enc_fit_floor_gain(q);
        }
    }
}

static void
ambe2400_enc_prediction(struct ambe2400_enc_frame* q, const mbe_parms* prev_mp) {
    int prev_L;
    /* Spectral prediction and residual */
    prev_L = mbe_clamp_harmonic_count(prev_mp->L);
    q->L = mbe_clamp_harmonic_count(q->L);
    float prev_log2Ml[57];
    memcpy(prev_log2Ml, prev_mp->log2Ml, sizeof(prev_log2Ml));
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

    for (int l = 1; l <= q->L; l++) {
        q->Tl[l] = q->a[l] - q->mean_a - (MBE_AMBE2400_PREDICTION_RHO * (q->p[l] - q->mean_p));
    }
}

static void
ambe2400_enc_block_dct(struct ambe2400_enc_frame* q) {
    /* Block DCTs */
    q->Ji[1] = AmbePlusLmprbl[q->L][0];
    q->Ji[2] = AmbePlusLmprbl[q->L][1];
    q->Ji[3] = AmbePlusLmprbl[q->L][2];
    q->Ji[4] = AmbePlusLmprbl[q->L][3];

    {
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
}

static void
ambe2400_enc_prba_dct(struct ambe2400_enc_frame* q) {
    float Ri[9];
    /* Ri from block DC-ish terms */
    const float sqrt2 = 1.41421356237f;
    Ri[1] = q->Cik[1][1] + (sqrt2 * q->Cik[1][2]);
    Ri[2] = q->Cik[1][1] - (sqrt2 * q->Cik[1][2]);
    Ri[3] = q->Cik[2][1] + (sqrt2 * q->Cik[2][2]);
    Ri[4] = q->Cik[2][1] - (sqrt2 * q->Cik[2][2]);
    Ri[5] = q->Cik[3][1] + (sqrt2 * q->Cik[3][2]);
    Ri[6] = q->Cik[3][1] - (sqrt2 * q->Cik[3][2]);
    Ri[7] = q->Cik[4][1] + (sqrt2 * q->Cik[4][2]);
    Ri[8] = q->Cik[4][1] - (sqrt2 * q->Cik[4][2]);

    /* Gm: 8-point DCT of Ri (Gm[1] is the discarded DC term) */
    for (int m = 2; m <= 8; m++) {
        float sum = 0.0f;
        for (int i = 1; i <= 8; i++) {
            sum += Ri[i] * q->cache->ri_cos[m][i];
        }
        q->Gm[m] = sum / 8.0f;
    }
}

static void
ambe2400_enc_quantize_prba(struct ambe2400_enc_frame* q) {
    /* PRBA24 (b3): Gm[2..4] */
    {
        int best = 0;
        float best_err = 1e30f;
        for (int c = 0; c < 512; c++) {
            float err = 0.0f;
            for (int m = 0; m < 3; m++) {
                float d = q->Gm[2 + m] - AmbePlusPRBA24[c][m];
                err += d * d;
            }
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        q->b[3] = best;
    }

    /* PRBA58 (b4): Gm[5..8] */
    {
        int best = 0;
        float best_err = 1e30f;
        for (int c = 0; c < 128; c++) {
            float err = 0.0f;
            for (int m = 0; m < 4; m++) {
                float d = q->Gm[5 + m] - AmbePlusPRBA58[c][m];
                err += d * d;
            }
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        q->b[4] = best;
    }
}

static const float (*const ambe2400_hoc_tables[4])[4] = {AmbePlusHOCb5, AmbePlusHOCb6, AmbePlusHOCb7, AmbePlusHOCb8};

static void
ambe2400_enc_quantize_hoc(struct ambe2400_enc_frame* q) {
    /* HOC blocks (b5..b8): Cik[blk][3..min(Ji,6)] */
    {
        int codes[4];

        for (int blk = 0; blk < 4; blk++) {
            int ji = q->Ji[blk + 1];
            int kmax = (ji < 6) ? ji : 6;
            int best = 0;
            float best_err = 1e30f;

            /* Block 4 carries only bits 3..1, so only even rows are reachable. */
            for (int c = 0; c < 16; c += (blk == 3) ? 2 : 1) {
                float err = 0.0f;
                for (int k = 3; k <= kmax; k++) {
                    float d = q->Cik[blk + 1][k] - ambe2400_hoc_tables[blk][c][k - 3];
                    err += d * d;
                }
                if (err < best_err) {
                    best_err = err;
                    best = c;
                }
            }
            codes[blk] = best;
        }
        q->b[5] = codes[0];
        q->b[6] = codes[1];
        q->b[7] = codes[2];
        q->b[8] = codes[3];
    }
}

static void
ambe2400_enc_reconstruct_coefficients(struct ambe2400_enc_frame* q) {
    float Ri_q[9];
    mbe_ambe2400_reconstruct_prba(q->b[3], q->b[4], Ri_q);
    mbe_ambe2400_reconstruct_cik(Ri_q, q->b + 5, q->Ji, q->Cik_q);
}

static void
ambe2400_enc_pack(const struct ambe2400_enc_frame* q, char ambe_d[49]) {
    /* Pack the 49-bit AMBE word */
    memset(ambe_d, 0, 49);
    for (int i = 0; i < 6; i++) {
        ambe_d[i] = (char)((q->b[0] >> (6 - i)) & 1);
    }
    ambe_d[48] = (char)(q->b[0] & 1);

    ambe_d[38] = (char)((q->b[1] >> 3) & 1);
    ambe_d[39] = (char)((q->b[1] >> 2) & 1);
    ambe_d[40] = (char)((q->b[1] >> 1) & 1);
    ambe_d[41] = (char)(q->b[1] & 1);

    ambe_d[6] = (char)((q->b[2] >> 5) & 1);
    ambe_d[7] = (char)((q->b[2] >> 4) & 1);
    ambe_d[8] = (char)((q->b[2] >> 3) & 1);
    ambe_d[9] = (char)((q->b[2] >> 2) & 1);
    ambe_d[42] = (char)((q->b[2] >> 1) & 1);
    ambe_d[43] = (char)(q->b[2] & 1);

    ambe_d[10] = (char)((q->b[3] >> 8) & 1);
    ambe_d[11] = (char)((q->b[3] >> 7) & 1);
    ambe_d[12] = (char)((q->b[3] >> 6) & 1);
    ambe_d[13] = (char)((q->b[3] >> 5) & 1);
    ambe_d[14] = (char)((q->b[3] >> 4) & 1);
    ambe_d[15] = (char)((q->b[3] >> 3) & 1);
    ambe_d[16] = (char)((q->b[3] >> 2) & 1);
    ambe_d[44] = (char)((q->b[3] >> 1) & 1);
    ambe_d[45] = (char)(q->b[3] & 1);

    ambe_d[17] = (char)((q->b[4] >> 6) & 1);
    ambe_d[18] = (char)((q->b[4] >> 5) & 1);
    ambe_d[19] = (char)((q->b[4] >> 4) & 1);
    ambe_d[20] = (char)((q->b[4] >> 3) & 1);
    ambe_d[21] = (char)((q->b[4] >> 2) & 1);
    ambe_d[46] = (char)((q->b[4] >> 1) & 1);
    ambe_d[47] = (char)(q->b[4] & 1);

    ambe_d[22] = (char)((q->b[5] >> 3) & 1);
    ambe_d[23] = (char)((q->b[5] >> 2) & 1);
    ambe_d[25] = (char)((q->b[5] >> 1) & 1);
    ambe_d[26] = (char)(q->b[5] & 1);

    ambe_d[27] = (char)((q->b[6] >> 3) & 1);
    ambe_d[28] = (char)((q->b[6] >> 2) & 1);
    ambe_d[29] = (char)((q->b[6] >> 1) & 1);
    ambe_d[30] = (char)(q->b[6] & 1);

    ambe_d[31] = (char)((q->b[7] >> 3) & 1);
    ambe_d[32] = (char)((q->b[7] >> 2) & 1);
    ambe_d[33] = (char)((q->b[7] >> 1) & 1);
    ambe_d[34] = (char)(q->b[7] & 1);

    ambe_d[35] = (char)((q->b[8] >> 3) & 1);
    ambe_d[36] = (char)((q->b[8] >> 2) & 1);
    ambe_d[37] = (char)((q->b[8] >> 1) & 1);

    /* ambe_d[24] is the spare bit; filled by the FEC layer below. */
}

static void
ambe2400_enc_fill_parms(const struct ambe2400_enc_frame* q, mbe_parms* cur_mp) {
    /* Write quantized parameters to cur_mp (decoder-equivalent state) */
    cur_mp->w0 = q->f0q * (float)2 * M_PI;
    cur_mp->L = q->L;
    cur_mp->K = 0;
    cur_mp->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
    for (int l = 1; l <= q->L; l++) {
        int jl = (int)((float)l * 16.0f * q->f0q);
        if (jl > 7) {
            jl = 7;
        }
        cur_mp->Vl[l] = AmbePlusVuv[q->b[1]][jl];
        if (cur_mp->Vl[l] == 1) {
            cur_mp->K++;
        }
    }
    cur_mp->gamma = q->gamma_q;
}

/* Analyze and quantize a voice frame, reconstruct its predictor state, and
 * pack 49 bits. */
static int
ambe2400_encode_voice(mbe_ambe2400_encoder* enc, char ambe_d[49], mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    struct mbe_analysis_result res;
    int status = mbe_analysis_frame(&enc->tables, &enc->analysis, enc->fft, enc->acf, &res);
    if (status < 0) {
        return status;
    }
    struct ambe2400_enc_frame q = {0};
    unsigned char columns[MBE_ANALYSIS_COLUMNS];
    q.cache = mbe_ambe2400_get_dct_cache();
    ambe2400_enc_quantize_pitch(&q, res.f0);
    ambe2400_enc_voicing(&q, &enc->analysis, &res, columns);
    ambe2400_enc_magnitudes(&q, &res);
    ambe2400_enc_quantize_gain(&q, prev_mp);
    ambe2400_enc_prediction(&q, prev_mp);
    ambe2400_enc_block_dct(&q);
    ambe2400_enc_prba_dct(&q);
    ambe2400_enc_quantize_prba(&q);
    ambe2400_enc_quantize_hoc(&q);
    ambe2400_enc_reconstruct_coefficients(&q);
    mbe_ambe2400_inverse_dct_tl(q.Cik_q, q.Ji, q.Tl_q);
    ambe2400_enc_pack(&q, ambe_d);
    ambe2400_enc_fill_parms(&q, cur_mp);
    /* Share the decoder update without changing the caller's predictor. */
    mbe_parms prediction_prev = *prev_mp;
    mbe_ambe2400_update_spectral_amplitudes(cur_mp, &prediction_prev, q.Tl_q, 0.2046f / sqrtf(cur_mp->w0));
    mbe_analysis_commit(&enc->analysis, columns);
    return 0;
}

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of PCM into AMBE 2400
 *        parameter bits.
 *
 * @param samples Input PCM floats (160), nominal range [-1, 1]; a non-finite
 *                sample or one beyond +-2^20 rejects the frame.
 * @param ambe_d  Output parameter bits (49).
 * @param cur_mp  Output: quantized (decoder-equivalent) parameters.
 * @param prev_mp Input: previous frame state (see mbe_initMbeParms()).
 * @return 0, or negative on error (the context is unchanged).
 */
/* Finite and within +-2^20 (0x49800000), checked on the bit pattern so the
 * test survives fast-math. One bad sample would otherwise poison the DC
 * filter and analysis history for the rest of the stream. */
static bool
ambe2400_enc_samples_valid(const float* samples) {
    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        uint32_t bits;
        memcpy(&bits, &samples[i], sizeof(bits));
        if ((bits & 0x7FFFFFFFu) > 0x49800000u) {
            return false;
        }
    }
    return true;
}

int
mbe_encodeAmbe2400Parms(mbe_ambe2400_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                        const mbe_parms* prev_mp) {
    float scaled[AMBE2400_ENC_SAMPLES];

    if (enc == NULL || samples == NULL || ambe_d == NULL || cur_mp == NULL || prev_mp == NULL
        || !ambe2400_enc_samples_valid(samples)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    /* DC-filter into the analysis history, then analyze and quantize. */
    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        scaled[i] = samples[i] * AMBE2400_ENC_PCM_SCALE;
    }
    mbe_analysis_push(&enc->analysis, scaled);
    return ambe2400_encode_voice(enc, ambe_d, cur_mp, prev_mp);
}

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of 16-bit PCM into AMBE 2400
 *        parameter bits.
 *
 * @see mbe_encodeAmbe2400Parms for details.
 */
int
mbe_encodeAmbe2400ParmsShort(mbe_ambe2400_encoder* enc, const short* samples, char ambe_d[49], mbe_parms* cur_mp,
                             const mbe_parms* prev_mp) {
    float float_buf[AMBE2400_ENC_SAMPLES];

    if (samples == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        float_buf[i] = (float)samples[i] / 32768.0f;
    }
    return mbe_encodeAmbe2400Parms(enc, float_buf, ambe_d, cur_mp, prev_mp);
}

static void
ambe2400_enc_c0(const char* ambe_d, char ambe_fr[4][24]) {
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
ambe2400_enc_pr(const char ambe_fr[4][24], unsigned short pr[25]) {
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

static void
ambe2400_enc_c1(const char* ambe_d, char ambe_fr[4][24], const unsigned short pr[25]) {
    char cw[23], data12[12];
    /* C1: 12 data bits (ambe_d[12..23]) -> Golay, XOR p bits 23..1 */
    for (int i = 0; i < 12; i++) {
        data12[i] = ambe_d[12 + i];
    }
    mbe_golay2312_encode(data12, cw);
    for (int j = 0; j < 23; j++) {
        /* fr[1][j] = b word bit (j+1), scrambled with p bit (j+1) = pr[24-(j+1)] >> 15 */
        int prbit = pr[23 - j] / 32768;
        ambe_fr[1][j] = (char)((cw[j] & 1) ^ prbit);
    }

    /* Spare bit (ambe_d[24] = fr[2][10] on air): b word bit 0 = parity ^ p bit 0 */
    int ones = 0;
    for (int j = 0; j < 23; j++) {
        ones += (cw[j] & 1);
    }
    ambe_fr[2][10] = (char)(((ones & 1) ^ (pr[24] / 32768)) & 1);
}

static void
ambe2400_enc_raw(const char* ambe_d, char ambe_fr[4][24]) {
    /* C2/C3: raw bits (decoder reads them MSB-first) */
    for (int j = 0; j < 11; j++) {
        ambe_fr[2][j] = ambe_d[24 + (10 - j)];
    }
    for (int j = 0; j < 14; j++) {
        ambe_fr[3][j] = ambe_d[35 + (13 - j)];
    }
}

/**
 * @brief Encode 49 AMBE 2400 parameter bits into a 72-bit D-STAR DV
 *        data frame (FEC + interleave), in the same plane layout the
 *        decoder (mbe_decodeAmbe3600x2400Frame) consumes.
 *
 * ambe_d[24] is the spare bit; on output it carries the scrambled even
 * parity of the second Golay codeword. All other input bits round-trip
 * exactly.
 *
 * @param ambe_d  Input parameter bits (49).
 * @param ambe_fr Output frame as 4x24 bitplanes.
 * @return 0 on success.
 */
int
mbe_encodeAmbe3600x2400Frame(const char ambe_d[49], char ambe_fr[4][24]) {
    unsigned short pr[25];

    if (ambe_d == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 49; i++) {
        if (ambe_d[i] != 0 && ambe_d[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }

    memset(ambe_fr, 0, 4 * sizeof(ambe_fr[0]));

    ambe2400_enc_c0(ambe_d, ambe_fr);
    ambe2400_enc_pr((const char (*)[24])ambe_fr, pr);
    ambe2400_enc_raw(ambe_d, ambe_fr);
    ambe2400_enc_c1(ambe_d, ambe_fr, pr);

    return 0;
}
