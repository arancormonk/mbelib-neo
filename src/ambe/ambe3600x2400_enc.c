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
#include <stdlib.h>
#include <string.h>

#include "ambe3600x2400_const.h"
#include "ambe3600x2400_internal.h"
#include "ambe_encoder.h"
#include "mbe_encoder.h"
#include "mbe_speech_analysis.h"
#include "mbelib-neo/mbelib.h"

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

/* Caller-owned encoder state: one context per stream. The analysis tables
 * are computed once at allocation and kept across resets. */
struct mbe_ambe2400_encoder {
    struct mbe_encoder_frontend fe;
};

static const struct ambe_enc_tables ambe2400_enc_tables = {
    .dg = AmbePlusDg,
    .dg_count = 64,
    .prba24 = AmbePlusPRBA24,
    .prba58 = AmbePlusPRBA58,
    .hoc = {AmbePlusHOCb5, AmbePlusHOCb6, AmbePlusHOCb7, AmbePlusHOCb8},
    .hoc_rows = {16, 16, 16, 16},
    .hoc_step = {1, 1, 1, 2}, /* block 4 carries only bits 3..1 */
    .lmprbl = AmbePlusLmprbl,
    .rho = MBE_AMBE2400_PREDICTION_RHO,
    .fit_floor_gain = 1,
};

void
mbe_ambe2400EncoderReset(mbe_ambe2400_encoder* enc) {
    if (enc == NULL) {
        return;
    }
    mbe_encoder_frontend_reset(&enc->fe);
}

mbe_ambe2400_encoder*
mbe_ambe2400EncoderAlloc(void) {
    mbe_ambe2400_encoder* enc = calloc(1, sizeof(*enc));
    if (enc == NULL) {
        return NULL;
    }
    if (mbe_encoder_frontend_open(&enc->fe) != 0) {
        free(enc);
        return NULL;
    }
    return enc;
}

void
mbe_ambe2400EncoderFree(mbe_ambe2400_encoder* enc) {
    if (enc != NULL) {
        mbe_encoder_frontend_close(&enc->fe);
        free(enc);
    }
}

static void
ambe2400_enc_quantize_pitch(struct ambe_enc_frame* q, float f0) {
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
ambe2400_enc_column(const struct ambe_enc_frame* q, int l) {
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
ambe2400_enc_voicing(struct ambe_enc_frame* q, const struct mbe_analysis_state* state,
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
ambe2400_enc_magnitudes(struct ambe_enc_frame* q, const struct mbe_analysis_result* res) {
    for (int l = 1; l <= q->L; l++) {
        q->a[l] = log2f((res->magnitude[l] * AMBE2400_ENC_MAG_SCALE) + 1e-12f);
    }
    ambe_enc_mean_amplitude(q);
}

static void
ambe2400_enc_reconstruct_coefficients(struct ambe_enc_frame* q) {
    float Ri_q[9];
    mbe_ambe2400_reconstruct_prba(q->b[3], q->b[4], Ri_q);
    mbe_ambe2400_reconstruct_cik(Ri_q, q->b + 5, q->Ji, q->Cik_q);
}

static void
ambe2400_enc_pack(const struct ambe_enc_frame* q, char ambe_d[49]) {
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
ambe2400_enc_fill_parms(const struct ambe_enc_frame* q, mbe_parms* cur_mp) {
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
ambe2400_encode_voice(mbe_ambe2400_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                      const mbe_parms* prev_mp) {
    struct mbe_analysis_result res;
    int status = mbe_encoder_frontend_analyze(&enc->fe, samples, &res);
    if (status < 0) {
        return status;
    }
    struct ambe_enc_frame q = {0};
    unsigned char columns[MBE_ANALYSIS_COLUMNS];
    q.cache = mbe_ambe2400_get_dct_cache();
    ambe2400_enc_quantize_pitch(&q, res.f0);
    ambe2400_enc_voicing(&q, &enc->fe.analysis, &res, columns);
    ambe2400_enc_magnitudes(&q, &res);
    ambe_enc_quantize_gain(&q, &ambe2400_enc_tables, prev_mp);
    ambe_enc_prediction(&q, &ambe2400_enc_tables, prev_mp);
    ambe_enc_quantize_residual(&q, &ambe2400_enc_tables);
    ambe2400_enc_reconstruct_coefficients(&q);
    mbe_ambe2400_inverse_dct_tl(q.Cik_q, q.Ji, q.Tl_q);
    ambe2400_enc_pack(&q, ambe_d);
    ambe2400_enc_fill_parms(&q, cur_mp);
    /* Share the decoder update without changing the caller's predictor. */
    mbe_parms prediction_prev = *prev_mp;
    mbe_ambe2400_update_spectral_amplitudes(cur_mp, &prediction_prev, q.Tl_q, 0.2046f / sqrtf(cur_mp->w0));
    if (!mbe_encoder_model_valid(cur_mp)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    mbe_analysis_commit(&enc->fe.analysis, columns);
    return 0;
}

int
mbe_encodeAmbe2400Parms(mbe_ambe2400_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                        const mbe_parms* prev_mp) {
    if (enc == NULL || samples == NULL || ambe_d == NULL || cur_mp == NULL || prev_mp == NULL
        || !mbe_encoder_samples_valid(samples) || !mbe_encoder_history_valid(prev_mp)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    /* A frame that fails (a history whose model would overflow) leaves no trace. */
    const struct mbe_analysis_state saved = enc->fe.analysis;
    mbe_parms out = *cur_mp;
    int status = ambe2400_encode_voice(enc, samples, ambe_d, &out, prev_mp);
    if (status < 0) {
        enc->fe.analysis = saved;
        return status;
    }
    *cur_mp = out;
    return 0;
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
    float float_buf[MBE_ENCODER_SAMPLES];

    if (samples == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    mbe_encoder_short_to_float(samples, float_buf);
    return mbe_encodeAmbe2400Parms(enc, float_buf, ambe_d, cur_mp, prev_mp);
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
    if (ambe_d == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 49; i++) {
        if (ambe_d[i] != 0 && ambe_d[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }
    /* The spare bit (ambe_d[24], fr[2][10] on air) carries C1's scrambled
     * parity instead. */
    ambe_fr[2][10] = (char)ambe_enc_fec(ambe_d, ambe_fr);
    return 0;
}
