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
 * @brief AMBE+2 3600x2450 (DMR, NXDN, YSF, P25 Phase 2) speech encoder.
 *
 * C port of ham_digital_modes' float AMBE+2 half-rate encoder
 * (ambe_plus_2/encoder.rs and mbe_encode.rs):
 *
 *  - Pitch from the TIA-102.BABA analysis (mbe_frame_analysis.c), quantized
 *    to the nearest AmbeW0table entry; voicing and amplitudes are then
 *    analysed at the quantized pitch.
 *  - b1..b8 invert this library's mbe_decodeAmbe2450Parms(): the V/UV row
 *    agreeing best with the analysis (amplitude weighted), the gain step that
 *    reproduces the mean log amplitude, and nearest PRBA and higher-order
 *    coefficient vectors from the forward transforms of the decoder's.
 *  - The prediction history is updated by decoding each emitted voice frame
 *    with mbe_decodeAmbe2450Parms(), so that the encoder's prediction remains
 *    identical to the decoder's.
 *  - DTMF digits and single tones are sent as TIA-102.BABA-1 tone frames
 *    (section 7, Table 10), which the decoder plays as tones.
 *
 * Deviation from the reference: the amplitude target is taken as log2,
 * matching this library's decoder (exp2f), where the reference inverts
 * mbelib's exp(0.693 * log2Ml).
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "ambe3600x2450_const.h"
#include "mbe_ecc.h"
#include "mbe_frame_analysis.h"
#include "mbe_tone.h"
#include "mbe_tone_detect.h"
#include "mbelib-neo/mbelib.h"

#define AMBE2450_ENC_SAMPLES   160
#define AMBE2450_ENC_PCM_SCALE 32768.0 /* analysis runs on the 16-bit scale */
#define AMBE2450_ENC_RHO       0.65    /* the decoder's amplitude predictor weight */
#define AMBE2450_ENC_GAMMA_MEM 0.5     /* gamma = DG + 0.5 * gamma_prev */

struct mbe_ambe2450_encoder {
    struct mbe_fa_tables tables;
    struct mbe_fa_state analysis;
    struct mbe_fa_voicing_state voicing;
    struct mbe_fa_result frame;
    mbe_parms prev; /* the decoder's prediction history */
};

void
mbe_ambe2450EncoderReset(mbe_ambe2450_encoder* enc) {
    mbe_parms cur;
    mbe_parms prev_enhanced;

    if (enc == NULL) {
        return;
    }
    mbe_fa_reset(&enc->analysis);
    mbe_fa_voicing_reset(&enc->voicing);
    mbe_initMbeParms(&cur, &enc->prev, &prev_enhanced);
}

mbe_ambe2450_encoder*
mbe_ambe2450EncoderAlloc(void) {
    mbe_ambe2450_encoder* enc = calloc(1, sizeof(*enc));
    if (enc == NULL) {
        return NULL;
    }
    mbe_fa_init_tables(&enc->tables);
    mbe_ambe2450EncoderReset(enc);
    return enc;
}

void
mbe_ambe2450EncoderFree(mbe_ambe2450_encoder* enc) {
    free(enc);
}

/* Frames emitted before the first analysed frame are encoded as silent
 * input: a zero spectrum, the standard's initial pitch period of 100 samples,
 * and a pitch error of 1 (no periodicity). */
static void
ambe2450_enc_quiet_frame(struct mbe_fa_result* frame) {
    memset(frame, 0, sizeof(*frame));
    frame->initial_pitch = 100.0;
    frame->omega0_hat = 2.0 * M_PI / 100.0;
    frame->initial_pitch_error = 1.0;
}

static int
ambe2450_enc_quantize_pitch(double omega0) {
    double f0 = omega0 / (2.0 * M_PI);
    int best = 0;
    for (int i = 1; i < 120; i++) {
        if (fabs((double)AmbeW0table[i] - f0) < fabs((double)AmbeW0table[best] - f0)) {
            best = i;
        }
    }
    return best;
}

/* The decoder's previous log2Ml: index 0 is read as index 1, and harmonics
 * past the previous frame's last hold its value. */
static double
ambe2450_enc_prev_log2(const mbe_parms* prev, int prev_l, int j) {
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

/* Index of the table row nearest the target in its first `used`
 * coefficients (squared error); ties resolve to the lowest index. */
static int
ambe2450_enc_nearest(const float* table, int rows, int cols, const double* target, int used) {
    int best = 0;
    double best_d = 0.0;
    for (int i = 0; i < rows; i++) {
        double d = 0.0;
        for (int k = 0; k < used; k++) {
            double e = (double)table[((size_t)i * (size_t)cols) + (size_t)k] - target[k];
            d += e * e;
        }
        if (i == 0 || d < best_d) {
            best_d = d;
            best = i;
        }
    }
    return best;
}

/* One frame's quantization target and workspace; no inter-frame state. */
struct ambe2450_enc_frame {
    int L;
    double f0;                   /* quantized fundamental, cycles per sample */
    const unsigned char* voiced; /* voiced[1..L] */
    const double* ml;            /* amplitudes ml[1..L] */
    double pred[57];             /* the decoder's predicted log2Ml */
    double sum43;                /* Sum43: rho times the mean interpolated log2Ml */
    double x[57];                /* target log2Ml less the prediction */
    double mean_x;
    double cik[5][18]; /* block DCT coefficients */
};

/* b1: the V/UV row agreeing best with the target, weighted by amplitude. */
static int
ambe2450_enc_quantize_vuv(const struct ambe2450_enc_frame* q) {
    double best_score = 0.0;
    int best = 0;
    for (int row = 0; row < 32; row++) {
        double score = 0.0;
        for (int h = 1; h <= q->L; h++) {
            int jl = (int)((double)h * 16.0 * q->f0);
            if (jl > 7) {
                jl = 7;
            }
            double w = fmax(q->ml[h], 1e-6);
            score += (AmbeVuv[row][jl] == (int)q->voiced[h]) ? w : -w;
        }
        if (row == 0 || score > best_score) {
            best_score = score;
            best = row;
        }
    }
    return best;
}

/* The decoder's prediction terms, with its float harmonic mapping. */
static void
ambe2450_enc_prediction(struct ambe2450_enc_frame* q, const mbe_parms* prev) {
    int prev_l = prev->L;
    if (prev_l < 1) {
        prev_l = 1;
    }
    if (prev_l > 56) {
        prev_l = 56;
    }
    q->sum43 = 0.0;
    for (int h = 1; h <= q->L; h++) {
        float flokl = ((float)prev_l / (float)q->L) * (float)h;
        int intkl = (int)flokl;
        double delta = (double)(flokl - (float)intkl);
        double interp = ((1.0 - delta) * ambe2450_enc_prev_log2(prev, prev_l, intkl))
                        + (delta * ambe2450_enc_prev_log2(prev, prev_l, intkl + 1));
        q->sum43 += interp;
        q->pred[h] = AMBE2450_ENC_RHO * interp;
    }
    q->sum43 *= AMBE2450_ENC_RHO / (double)q->L;
}

/* Target log2Ml less the prediction, undoing the decoder's unvoiced scaling;
 * returns b2, the gain step reproducing the mean log amplitude. */
static int
ambe2450_enc_quantize_gain(struct ambe2450_enc_frame* q, const mbe_parms* prev) {
    double w0 = q->f0 * 2.0 * M_PI;
    double unvc = (double)(0.2046f / sqrtf((float)w0));
    q->mean_x = 0.0;
    for (int h = 1; h <= q->L; h++) {
        double m = fmax(q->ml[h], 1e-30);
        double eff = q->voiced[h] ? m : m / unvc;
        q->x[h] = log2(eff) - q->pred[h];
        q->mean_x += q->x[h];
    }
    q->mean_x /= (double)q->L;

    double gamma_target = q->mean_x + q->sum43 + (0.5 * log2((double)q->L));
    double delta_target = gamma_target - (AMBE2450_ENC_GAMMA_MEM * (double)prev->gamma);
    int best = 0;
    for (int i = 1; i < 32; i++) {
        if (fabs((double)AmbeDg[i] - delta_target) < fabs((double)AmbeDg[best] - delta_target)) {
            best = i;
        }
    }
    return best;
}

/* Split the zero-mean residual Tl into the decoder's four blocks and apply
 * the DCT to each. */
static void
ambe2450_enc_block_dct(struct ambe2450_enc_frame* q) {
    int start = 1;
    memset(q->cik, 0, sizeof(q->cik));
    for (int block = 0; block < 4; block++) {
        int j_len = AmbeLmprbl[q->L][block];
        int k_max = (j_len < 17) ? j_len : 17;
        for (int k = 1; k <= k_max; k++) {
            double sum = 0.0;
            for (int j = 1; j <= j_len && start + j - 1 <= q->L; j++) {
                double tl = q->x[start + j - 1] - q->mean_x;
                sum += tl * cos(M_PI * ((double)k - 1.0) * ((double)j - 0.5) / (double)j_len);
            }
            q->cik[block + 1][k] = sum / (double)j_len;
        }
        start += j_len;
    }
}

/* b3, b4: Ri from each block's first two coefficients, Gm by the 8-point
 * forward transform, then the nearest PRBA vectors. */
static void
ambe2450_enc_quantize_prba(const struct ambe2450_enc_frame* q, int b[9]) {
    double ri[9];
    double gm[9];
    for (int block = 1; block <= 4; block++) {
        double c1 = q->cik[block][1];
        double c2 = q->cik[block][2];
        int odd = (2 * block) - 1;
        ri[odd] = c1 + (M_SQRT2 * c2);
        ri[odd + 1] = c1 - (M_SQRT2 * c2);
    }
    for (int m = 1; m <= 8; m++) {
        double sum = 0.0;
        for (int i = 1; i <= 8; i++) {
            sum += ri[i] * cos(M_PI * ((double)m - 1.0) * ((double)i - 0.5) / 8.0);
        }
        gm[m] = sum / 8.0;
    }
    b[3] = ambe2450_enc_nearest(&AmbePRBA24[0][0], 512, 3, &gm[2], 3);
    b[4] = ambe2450_enc_nearest(&AmbePRBA58[0][0], 128, 4, &gm[5], 4);
}

/* b5..b8: higher-order coefficients; the decoder reads C[block][3..min(J, 6)]. */
static void
ambe2450_enc_quantize_hoc(const struct ambe2450_enc_frame* q, int b[9]) {
    static const float* const hoc_tables[4] = {&AmbeHOCb5[0][0], &AmbeHOCb6[0][0], &AmbeHOCb7[0][0], &AmbeHOCb8[0][0]};
    static const int hoc_rows[4] = {32, 16, 16, 8};
    for (int block = 0; block < 4; block++) {
        int used = AmbeLmprbl[q->L][block] - 2;
        if (used > 4) {
            used = 4;
        }
        b[5 + block] =
            (used <= 0) ? 0 : ambe2450_enc_nearest(hoc_tables[block], hoc_rows[block], 4, &q->cik[block + 1][3], used);
    }
}

/* Choose b1..b8 for a target of L harmonics at fundamental f0 (cycles per
 * sample), voicing voiced[1..L] and amplitudes ml[1..L]. */
static void
ambe2450_enc_quantize(int L, double f0, const unsigned char* voiced, const double* ml, const mbe_parms* prev,
                      int b[9]) {
    struct ambe2450_enc_frame q;
    q.L = L;
    q.f0 = f0;
    q.voiced = voiced;
    q.ml = ml;
    b[1] = ambe2450_enc_quantize_vuv(&q);
    ambe2450_enc_prediction(&q, prev);
    b[2] = ambe2450_enc_quantize_gain(&q, prev);
    ambe2450_enc_block_dct(&q);
    ambe2450_enc_quantize_prba(&q, b);
    ambe2450_enc_quantize_hoc(&q, b);
}

static void
ambe2450_enc_set_bits(char ambe_d[49], int start, int width, int value) {
    for (int i = 0; i < width; i++) {
        ambe_d[start + i] = (char)((value >> (width - 1 - i)) & 1);
    }
}

/* Scatter b0..b8 into the 49-bit layout mbe_decodeAmbe2450Parms() reads. */
static void
ambe2450_enc_pack(const int b[9], char ambe_d[49]) {
    ambe2450_enc_set_bits(ambe_d, 0, 4, b[0] >> 3);
    ambe2450_enc_set_bits(ambe_d, 37, 3, b[0] & 7);
    ambe2450_enc_set_bits(ambe_d, 4, 4, b[1] >> 1);
    ambe2450_enc_set_bits(ambe_d, 35, 1, b[1] & 1);
    ambe2450_enc_set_bits(ambe_d, 8, 4, b[2] >> 1);
    ambe2450_enc_set_bits(ambe_d, 36, 1, b[2] & 1);
    ambe2450_enc_set_bits(ambe_d, 12, 8, b[3] >> 1);
    ambe2450_enc_set_bits(ambe_d, 40, 1, b[3] & 1);
    ambe2450_enc_set_bits(ambe_d, 20, 4, b[4] >> 3);
    ambe2450_enc_set_bits(ambe_d, 41, 3, b[4] & 7);
    ambe2450_enc_set_bits(ambe_d, 24, 4, b[5] >> 1);
    ambe2450_enc_set_bits(ambe_d, 44, 1, b[5] & 1);
    ambe2450_enc_set_bits(ambe_d, 28, 3, b[6] >> 1);
    ambe2450_enc_set_bits(ambe_d, 45, 1, b[6] & 1);
    ambe2450_enc_set_bits(ambe_d, 31, 3, b[7] >> 1);
    ambe2450_enc_set_bits(ambe_d, 46, 1, b[7] & 1);
    ambe2450_enc_set_bits(ambe_d, 34, 1, b[8] >> 2);
    ambe2450_enc_set_bits(ambe_d, 47, 2, b[8] & 3);
}

/* TIA-102.BABA-1 Table 9 tone index of a detected tone; DTMF follows the
 * decoder's dual-tone table (128 + digit, A-D = 0xA-0xD, * = 0xE, # = 0xF). */
static int
ambe2450_enc_tone_id(const struct mbe_tone_detection* det) {
    static const int nibble[4][4] = {
        {0x1, 0x2, 0x3, 0xA}, {0x4, 0x5, 0x6, 0xB}, {0x7, 0x8, 0x9, 0xC}, {0xE, 0x0, 0xF, 0xD}};
    if (det->kind == MBE_TONE_DETECT_DTMF) {
        return 0x80 | nibble[det->row & 3][det->col & 3];
    }
    return det->index;
}

/* The 7-bit tone level AD for a per-tone peak amplitude on the 16-bit scale:
 * the inverse of the level law the decoder plays tone frames with
 * (mbe_tone_ambe2450_level_db, peak sqrt(2) * 32768 * 10^(dB/20) after the
 * factor of 7 applied by the float-to-16-bit output conversion). An amplitude
 * of 1000 yields AD 84, the value the DVSI AMBE-3000 encoder transmits for
 * that tone. */
static int
ambe2450_enc_tone_level(double amplitude) {
    double level_db = 20.0 * log10(fmax(amplitude, 1.0) / (M_SQRT2 * 32768.0));
    double ad = 127.0 + ((level_db - MBE_TONE_AMBE2450_DB_AT_AD127) / MBE_TONE_AMBE2450_DB_PER_STEP);
    int level = (int)lround(ad);
    return (level < 0) ? 0 : ((level > 127) ? 127 : level);
}

/*
 * Tone frame (TIA-102.BABA-1 7.2, Table 10): u0 = six ones and AD(6..1),
 * the 8-bit index repeated from bit 12 to bit 43, AD(0) at bit 44, zeros
 * after. The repetitions of the index's high nibble occupy the low bits of
 * b0, so DTMF frames read as b0 = 120 and call-progress tones as b0 = 122,
 * as transmitted by the DVSI AMBE-3000 encoder.
 */
static void
ambe2450_enc_pack_tone(int tone_id, int level, char ambe_d[49]) {
    memset(ambe_d, 0, 49);
    ambe2450_enc_set_bits(ambe_d, 0, 6, 0x3F);
    ambe2450_enc_set_bits(ambe_d, 6, 6, level >> 1);
    for (int i = 12; i < 44; i += 8) {
        ambe2450_enc_set_bits(ambe_d, i, 8, tone_id);
    }
    ambe_d[44] = (char)(level & 1);
}

static int
ambe2450_encode(mbe_ambe2450_encoder* enc, const double input[AMBE2450_ENC_SAMPLES], char ambe_d[49]) {
    struct mbe_fa_voicing_state scratch;
    struct mbe_fa_voicing_state* voicing = &enc->voicing;
    const struct mbe_fa_result* a = &enc->frame;
    unsigned char bands[MBE_FA_MAX_BANDS];
    unsigned char voiced[57] = {0};
    double amplitudes[57] = {0};
    double ml[57] = {0};
    int b[9];

    if (mbe_fa_push(&enc->tables, &enc->analysis, input, &enc->frame)) {
        struct mbe_tone_detection det;
        if (mbe_tone_detect(a->slot, &det) != MBE_TONE_DETECT_NONE) {
            ambe2450_enc_pack_tone(ambe2450_enc_tone_id(&det), ambe2450_enc_tone_level(det.amplitude), ambe_d);
            return MBE_AMBE2450_FRAME_TONE;
        }
    } else {
        /* While the look-ahead buffer is filling, a copy of the voicing state
         * is used, so that only analysed frames update the original. */
        ambe2450_enc_quiet_frame(&enc->frame);
        scratch = enc->voicing;
        voicing = &scratch;
    }

    b[0] = ambe2450_enc_quantize_pitch(a->omega0_hat);
    int L = (int)AmbeLtable[b[0]];
    double f0 = (double)AmbeW0table[b[0]];
    double w0 = f0 * 2.0 * M_PI;

    /* Voicing and amplitudes at the quantized pitch, the band decisions
     * extended to the decoder's harmonic count. */
    int decided = mbe_fa_determine_voicing(&enc->tables, &a->sw, w0, a->initial_pitch_error, voicing, bands);
    int k_hat = mbe_fa_bands_count(L);
    unsigned char fill = (decided > 0) ? bands[decided - 1] : 0;
    for (int k = decided; k < k_hat; k++) {
        bands[k] = fill;
    }
    mbe_fa_spectral_amplitudes(&enc->tables, &a->sw, L, k_hat, w0, bands, amplitudes);
    for (int h = 1; h <= L; h++) {
        int band = (h + 2) / 3;
        if (band > k_hat) {
            band = k_hat;
        }
        if (band < 1) {
            band = 1;
        }
        voiced[h] = bands[band - 1];
        /* Amplitudes are floored at 1e-3 (-60 dB) so that the logarithm of a
         * near-zero amplitude does not dominate the DCT coefficients of its
         * block. */
        ml[h] = fmax(amplitudes[h], 1e-3);
    }

    ambe2450_enc_quantize(L, f0, voiced, ml, &enc->prev, b);
    ambe2450_enc_pack(b, ambe_d);

    /* Update the prediction history by decoding the emitted frame, exactly as
     * the receiving decoder will. */
    mbe_parms cur = enc->prev;
    if (mbe_decodeAmbe2450Parms(ambe_d, &cur, &enc->prev) == MBE_AMBE2450_FRAME_VOICE) {
        mbe_moveMbeParms(&cur, &enc->prev);
    }
    return MBE_AMBE2450_FRAME_VOICE;
}

/* Accept only finite samples of magnitude at most 2^20 (bit pattern
 * 0x49800000). The test inspects the IEEE 754 bit pattern rather than using
 * isfinite() or comparisons, because the library may be compiled with
 * MBELIB_ENABLE_FAST_MATH (-ffast-math, /fp:fast), under which the compiler
 * is permitted to assume that no value is NaN or infinite and may remove such
 * checks. */
static int
ambe2450_enc_samples_valid(const float* samples) {
    for (int i = 0; i < AMBE2450_ENC_SAMPLES; i++) {
        uint32_t bits;
        memcpy(&bits, &samples[i], sizeof(bits));
        if ((bits & 0x7FFFFFFFu) > 0x49800000u) {
            return 0;
        }
    }
    return 1;
}

int
mbe_encodeAmbe2450Parms(mbe_ambe2450_encoder* enc, const float* samples, char ambe_d[49]) {
    double input[AMBE2450_ENC_SAMPLES];

    if (enc == NULL || samples == NULL || ambe_d == NULL || !ambe2450_enc_samples_valid(samples)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < AMBE2450_ENC_SAMPLES; i++) {
        input[i] = (double)samples[i] * AMBE2450_ENC_PCM_SCALE;
    }
    return ambe2450_encode(enc, input, ambe_d);
}

int
mbe_encodeAmbe2450ParmsShort(mbe_ambe2450_encoder* enc, const short* samples, char ambe_d[49]) {
    double input[AMBE2450_ENC_SAMPLES];

    if (enc == NULL || samples == NULL || ambe_d == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < AMBE2450_ENC_SAMPLES; i++) {
        input[i] = (double)samples[i];
    }
    return ambe2450_encode(enc, input, ambe_d);
}

int
mbe_encodeAmbe3600x2450Frame(const char ambe_d[49], char ambe_fr[4][24]) {
    char data12[12];
    char cw[23];
    unsigned short pr[25];
    unsigned short foo = 0;

    if (ambe_d == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 49; i++) {
        if (ambe_d[i] != 0 && ambe_d[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }
    memset(ambe_fr, 0, 4 * sizeof(ambe_fr[0]));

    /* C0: (23,12) Golay plus the even parity bit fr[0][0]. */
    memcpy(data12, ambe_d, 12);
    mbe_golay2312_encode(data12, cw);
    int ones = 0;
    for (int j = 0; j < 23; j++) {
        ambe_fr[0][j + 1] = cw[j];
        ones += cw[j];
    }
    ambe_fr[0][0] = (char)(ones & 1);

    /* C1: (23,12) Golay, scrambled by the sequence seeded from C0's data. */
    for (int i = 23; i >= 12; i--) {
        foo = (unsigned short)((foo << 1) | (unsigned short)(ambe_fr[0][i] & 1));
    }
    pr[0] = (unsigned short)(16 * foo);
    for (int i = 1; i < 25; i++) {
        pr[i] = (unsigned short)(((173u * pr[i - 1]) + 13849u) & 0xFFFFu);
    }
    memcpy(data12, ambe_d + 12, 12);
    mbe_golay2312_encode(data12, cw);
    for (int j = 0; j < 23; j++) {
        ambe_fr[1][j] = (char)((cw[j] & 1) ^ (pr[23 - j] >> 15));
    }

    /* C2 and C3: unprotected, read MSB-first by the decoder. */
    for (int j = 0; j < 11; j++) {
        ambe_fr[2][j] = ambe_d[24 + (10 - j)];
    }
    for (int j = 0; j < 14; j++) {
        ambe_fr[3][j] = ambe_d[35 + (13 - j)];
    }
    return 0;
}
