// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 *
 * Copyright (C) 2010 mbelib Author
 * GPG Key ID: 0xEA5EFE2C (9E7A 5527 9CDC EBF7 BF1B  D772 4F98 E863 EA5E FE2C)
 *
 * Portions were originally under the ISC license; this mbelib-neo
 * distribution is provided under GPL-2.0-or-later. See LICENSE for details.
 */

/**
 * @file
 * @brief AMBE 3600x2400 parameter decode, ECC, and synthesis hooks.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ambe3600x2400_const.h"
#include "ambe3600x2400_internal.h"
#include "ambe_common.h"
#include "mbe_adaptive.h"
#include "mbe_compiler.h"
#include "mbe_repeat.h"
#include "mbe_result.h"
#include "mbe_tone.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"

/**
 * @brief Thread-local cache for AMBE DCT cosine coefficients.
 *
 * Pre-computes the cosine terms used in the Ri inverse DCT and per-block
 * inverse DCT loops to eliminate repeated cosf() calls per frame.
 *
 * Ri DCT: cosf(M_PI * (m-1) * (i-0.5) / 8) for m=1..8, i=1..8
 * Per-block IDCT: cosf(M_PI * (k-1) * (j-0.5) / ji) for ji=1..17, j=1..ji, k=1..ji
 */
static MBE_THREAD_LOCAL struct ambe_dct_cache ambe_cache = {0};

/**
 * @brief Initialize or return the thread-local AMBE DCT cache.
 *
 * Fills the cosine tables on first use. Because the cache is thread-local,
 * no locking is required.
 *
 * @return Pointer to the initialized cache.
 */
const struct ambe_dct_cache*
mbe_ambe2400_get_dct_cache(void) {
    if (ambe_cache.inited) {
        return &ambe_cache;
    }

    /* Fill Ri DCT cosine table: cosf(M_PI * (m-1) * (i-0.5) / 8) */
    for (int m = 1; m <= 8; m++) {
        for (int i = 1; i <= 8; i++) {
            ambe_cache.ri_cos[m][i] = cosf((M_PI * (float)(m - 1) * ((float)i - 0.5f)) / 8.0f);
        }
    }

    /* Fill per-block IDCT cosine table: cosf(M_PI * (k-1) * (j-0.5) / ji) */
    for (int ji = 1; ji <= 17; ji++) {
        for (int j = 1; j <= ji; j++) {
            for (int k = 1; k <= ji; k++) {
                ambe_cache.idct_cos[ji][j][k] = cosf((M_PI * (float)(k - 1) * ((float)j - 0.5f)) / (float)ji);
            }
        }
    }

    ambe_cache.inited = 1;
    return &ambe_cache;
}

/**
 * @brief Print AMBE 2400 parameter bits to stderr (debug aid).
 * @param ambe_d AMBE parameter bits (49).
 */
void
mbe_dumpAmbe2400Data(const char* ambe_d) {

    int i;
    const char* ambe;

    ambe = ambe_d;
    for (i = 0; i < 49; i++) {
        fprintf(stderr, "%i", *ambe);
        ambe++;
    }
    fprintf(stderr, " ");
}

/**
 * @brief Print raw AMBE 3600x2400 frame bitplanes to stderr.
 * @param ambe_fr Frame as 4x24 bitplanes.
 */
void
mbe_dumpAmbe3600x2400Frame(const char ambe_fr[4][24]) {

    int j;

    // c0
    fprintf(stderr, "ambe_fr c0: ");
    for (j = 23; j >= 0; j--) {
        fprintf(stderr, "%i", ambe_fr[0][j]);
    }
    fprintf(stderr, " ");
    // c1
    fprintf(stderr, "ambe_fr c1: ");
    for (j = 22; j >= 0; j--) {
        fprintf(stderr, "%i", ambe_fr[1][j]);
    }
    fprintf(stderr, " ");
    // c2
    fprintf(stderr, "ambe_fr c2: ");
    for (j = 10; j >= 0; j--) {
        fprintf(stderr, "%i", ambe_fr[2][j]);
    }
    fprintf(stderr, " ");
    // c3
    fprintf(stderr, "ambe_fr c3: ");
    for (j = 13; j >= 0; j--) {
        fprintf(stderr, "%i", ambe_fr[3][j]);
    }
    fprintf(stderr, " ");
}

/**
 * @brief Apply ECC to AMBE 3600x2400 C0 and update in-place.
 * @param ambe_fr Frame as 4x24 bitplanes.
 * @return Number of corrected errors in C0.
 */
int
mbe_eccAmbe3600x2400C0(char ambe_fr[4][24]) {
    int ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }
    return mbe_eccAmbe3600C0_common(ambe_fr);
}

/**
 * @brief Apply ECC to AMBE 3600x2400 data and pack parameter bits.
 * @param ambe_fr Frame as 4x24 bitplanes.
 * @param ambe_d  Output parameter bits (49).
 * @return Number of corrected errors in protected fields.
 */
int
mbe_eccAmbe3600x2400Data(char ambe_fr[4][24], char* ambe_d) {
    if (!ambe_d) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    int ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }
    return mbe_eccAmbe3600Data_common(ambe_fr, ambe_d);
}

static int
ambe2400_decode_b0(const char* ambe_d) {
    int b0 = 0;
    b0 |= ambe_d[0] << 6;
    b0 |= ambe_d[1] << 5;
    b0 |= ambe_d[2] << 4;
    b0 |= ambe_d[3] << 3;
    b0 |= ambe_d[4] << 2;
    b0 |= ambe_d[5] << 1;
    b0 |= ambe_d[48];
    return b0;
}

/*
 * b0 125..127: a tone index or MBE_AMBE2400_FRAME_LOST (see
 * ambe3600x2400_internal.h); 0 for voice. A b0 126 tone frame plays at any
 * volume, 0 included; b0 120..124, the null frame among them, are voice.
 */
static int
ambe2400_special_frame(const char* ambe_d, int b0) {
    if (b0 < 125) {
        return 0;
    }
    if (b0 == 126) {
        int tone_index = mbe_ambe2400_tone_index(ambe_d);
        float freq1, freq2;
        return mbe_tone_lookup_dstar_freqs(tone_index, &freq1, &freq2) ? tone_index : MBE_AMBE2400_FRAME_LOST;
    }
    return MBE_AMBE2400_FRAME_LOST;
}

static void
ambe2400_setup_voice_model(mbe_parms* cur_mp, int b0, int* L, float* f0) {
    *f0 = mbe_ambe2400_f0(b0);
    cur_mp->w0 = *f0 * (float)2 * M_PI;
    *L = mbe_ambe2400_harmonic_count(*f0);
    cur_mp->L = *L;
}

static void
ambe2400_decode_vuv(const char* ambe_d, mbe_parms* cur_mp, int b0, int L, float f0) {
    (void)b0;
    int b1 = 0;
    b1 |= ambe_d[38] << 3;
    b1 |= ambe_d[39] << 2;
    b1 |= ambe_d[40] << 1;
    b1 |= ambe_d[41];

    for (int l = 1; l <= L; l++) {
        int jl = (int)((float)l * (float)16.0 * f0);
        cur_mp->Vl[l] = AmbePlusVuv[b1][jl];
#ifdef AMBE_DEBUG
        fprintf(stderr, "jl[%i]:%i Vl[%i]:%i\n", l, jl, l, cur_mp->Vl[l]);
#endif
    }
#ifdef AMBE_DEBUG
    fprintf(stderr, "\nb0:%i w0:%f L:%i b1:%i\n", b0, cur_mp->w0, L, b1);
#endif
}

static void
ambe2400_decode_gain(const char* ambe_d, mbe_parms* cur_mp) {
    int b2 = 0;
    b2 |= ambe_d[6] << 5;
    b2 |= ambe_d[7] << 4;
    b2 |= ambe_d[8] << 3;
    b2 |= ambe_d[9] << 2;
    b2 |= ambe_d[42] << 1;
    b2 |= ambe_d[43];

    cur_mp->gamma = (MBE_AMBE2400_GAIN_SCALE * AmbePlusDg[b2]) + MBE_AMBE2400_GAIN_OFFSET;
#ifdef AMBE_DEBUG
    fprintf(stderr, "b2: %i, gamma: %f\n", b2, cur_mp->gamma);
#endif
}

static void
ambe2400_decode_ri(const char* ambe_d, float Ri[9]) {
    int b3 = 0;
    b3 |= ambe_d[10] << 8;
    b3 |= ambe_d[11] << 7;
    b3 |= ambe_d[12] << 6;
    b3 |= ambe_d[13] << 5;
    b3 |= ambe_d[14] << 4;
    b3 |= ambe_d[15] << 3;
    b3 |= ambe_d[16] << 2;
    b3 |= ambe_d[44] << 1;
    b3 |= ambe_d[45];

    int b4 = 0;
    b4 |= ambe_d[17] << 6;
    b4 |= ambe_d[18] << 5;
    b4 |= ambe_d[19] << 4;
    b4 |= ambe_d[20] << 3;
    b4 |= ambe_d[21] << 2;
    b4 |= ambe_d[46] << 1;
    b4 |= ambe_d[47];

    mbe_ambe2400_reconstruct_prba(b3, b4, Ri);
}

void
mbe_ambe2400_reconstruct_prba(int b3, int b4, float Ri[9]) {
    float Gm[9];
    Gm[1] = 0;
    Gm[2] = AmbePlusPRBA24[b3][0];
    Gm[3] = AmbePlusPRBA24[b3][1];
    Gm[4] = AmbePlusPRBA24[b3][2];
    Gm[5] = AmbePlusPRBA58[b4][0];
    Gm[6] = AmbePlusPRBA58[b4][1];
    Gm[7] = AmbePlusPRBA58[b4][2];
    Gm[8] = AmbePlusPRBA58[b4][3];
#ifdef AMBE_DEBUG
    fprintf(stderr, "b3: %i Gm[2]: %f Gm[3]: %f Gm[4]: %f b4: %i Gm[5]: %f Gm[6]: %f Gm[7]: %f Gm[8]: %f\n", b3, Gm[2],
            Gm[3], Gm[4], b4, Gm[5], Gm[6], Gm[7], Gm[8]);
#endif

    mbe_ambe2400_reconstruct_ri(Gm, Ri);
}

void
mbe_ambe2400_reconstruct_ri(const float Gm[9], float Ri[9]) {
    const struct ambe_dct_cache* cache = mbe_ambe2400_get_dct_cache();
    for (int i = 1; i <= 8; i++) {
        float sum = 0;
        for (int m = 1; m <= 8; m++) {
            int am = (m == 1) ? 1 : 2;
            sum = sum + ((float)am * Gm[m] * cache->ri_cos[m][i]);
        }
        Ri[i] = sum;
#ifdef AMBE_DEBUG
        fprintf(stderr, "R%i: %f ", i, Ri[i]);
#endif
    }
#ifdef AMBE_DEBUG
    fprintf(stderr, "\n");
#endif
}

static void
ambe2400_load_hoc_block(float Cik[5][18], int block, int ji, int code, const float hoc[16][4]) {
    for (int k = 3; k <= ji; k++) {
        if (k > 6) {
            Cik[block][k] = 0;
        } else {
            Cik[block][k] = hoc[code][k - 3];
#ifdef AMBE_DEBUG
            fprintf(stderr, "C%i,%i: %f ", block, k, Cik[block][k]);
#endif
        }
    }
}

static void
ambe2400_decode_cik(const char* ambe_d, int L, const float Ri[9], float Cik[5][18], int Ji[5]) {
    int b5 = 0;
    b5 |= ambe_d[22] << 3;
    b5 |= ambe_d[23] << 2;
    b5 |= ambe_d[25] << 1;
    b5 |= ambe_d[26];

    int b6 = 0;
    b6 |= ambe_d[27] << 3;
    b6 |= ambe_d[28] << 2;
    b6 |= ambe_d[29] << 1;
    b6 |= ambe_d[30];

    int b7 = 0;
    b7 |= ambe_d[31] << 3;
    b7 |= ambe_d[32] << 2;
    b7 |= ambe_d[33] << 1;
    b7 |= ambe_d[34];

    int b8 = 0;
    b8 |= ambe_d[35] << 3;
    b8 |= ambe_d[36] << 2;
    b8 |= ambe_d[37] << 1;

    Ji[1] = AmbePlusLmprbl[L][0];
    Ji[2] = AmbePlusLmprbl[L][1];
    Ji[3] = AmbePlusLmprbl[L][2];
    Ji[4] = AmbePlusLmprbl[L][3];
#ifdef AMBE_DEBUG
    fprintf(stderr, "Ji[1]: %i Ji[2]: %i Ji[3]: %i Ji[4]: %i\n", Ji[1], Ji[2], Ji[3], Ji[4]);
    fprintf(stderr, "b5: %i b6: %i b7: %i b8: %i\n", b5, b6, b7, b8);
#endif

    const int codes[4] = {b5, b6, b7, b8};
    mbe_ambe2400_reconstruct_cik(Ri, codes, Ji, Cik);
}

void
mbe_ambe2400_reconstruct_cik(const float Ri[9], const int codes[4], const int Ji[5], float Cik[5][18]) {
    const float rconst = ((float)1 / ((float)2 * M_SQRT2));
    Cik[1][1] = (float)0.5 * (Ri[1] + Ri[2]);
    Cik[1][2] = rconst * (Ri[1] - Ri[2]);
    Cik[2][1] = (float)0.5 * (Ri[3] + Ri[4]);
    Cik[2][2] = rconst * (Ri[3] - Ri[4]);
    Cik[3][1] = (float)0.5 * (Ri[5] + Ri[6]);
    Cik[3][2] = rconst * (Ri[5] - Ri[6]);
    Cik[4][1] = (float)0.5 * (Ri[7] + Ri[8]);
    Cik[4][2] = rconst * (Ri[7] - Ri[8]);

    ambe2400_load_hoc_block(Cik, 1, Ji[1], codes[0], AmbePlusHOCb5);
    ambe2400_load_hoc_block(Cik, 2, Ji[2], codes[1], AmbePlusHOCb6);
    ambe2400_load_hoc_block(Cik, 3, Ji[3], codes[2], AmbePlusHOCb7);
    ambe2400_load_hoc_block(Cik, 4, Ji[4], codes[3], AmbePlusHOCb8);
#ifdef AMBE_DEBUG
    fprintf(stderr, "\n");
#endif
}

void
mbe_ambe2400_inverse_dct_tl(float Cik[5][18], const int Ji[5], float Tl[57]) {
    const struct ambe_dct_cache* cache = mbe_ambe2400_get_dct_cache();
    int l = 1;
    for (int i = 1; i <= 4; i++) {
        int ji = Ji[i];
        for (int j = 1; j <= ji; j++) {
            float sum = 0;
            for (int k = 1; k <= ji; k++) {
                int ak = (k == 1) ? 1 : 2;
#ifdef AMBE_DEBUG
                fprintf(stderr, "j: %i Cik[%i][%i]: %f ", j, i, k, Cik[i][k]);
#endif
                sum = sum + ((float)ak * Cik[i][k] * cache->idct_cos[ji][j][k]);
            }
            Tl[l] = sum;
#ifdef AMBE_DEBUG
            fprintf(stderr, "Tl[%i]: %f\n", l, Tl[l]);
#endif
            l++;
        }
    }
}

void
mbe_ambe2400_update_spectral_amplitudes(mbe_parms* cur_mp, mbe_parms* prev_mp, const float Tl[57], float unvc) {
    int intkl[57] = {0};
    float flokl[57] = {0.0f};
    float deltal[57] = {0.0f};
    int prev_L = mbe_clamp_harmonic_count(prev_mp->L);
    cur_mp->L = mbe_clamp_harmonic_count(cur_mp->L);

    if (cur_mp->L > prev_L) {
        for (int l = prev_L + 1; l <= cur_mp->L; l++) {
            prev_mp->Ml[l] = prev_mp->Ml[prev_L];
            prev_mp->log2Ml[l] = prev_mp->log2Ml[prev_L];
        }
    }
    prev_mp->log2Ml[0] = prev_mp->log2Ml[1];
    prev_mp->Ml[0] = prev_mp->Ml[1];

    float Sum43 = 0;
    for (int l = 1; l <= cur_mp->L; l++) {
        flokl[l] = ambe2400_prediction_position(prev_L, cur_mp->L, l);
        intkl[l] = (int)flokl[l];
#ifdef AMBE_DEBUG
        fprintf(stderr, "flok%i: %f, intk%i: %i ", l, flokl[l], l, intkl[l]);
#endif
        deltal[l] = flokl[l] - (float)intkl[l];
#ifdef AMBE_DEBUG
        fprintf(stderr, "delta%i: %f ", l, deltal[l]);
#endif
        /* intkl[l] can reach 56; keep the upper interpolation tap inside log2Ml[57]. */
        int upper = intkl[l] + 1;
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }
        Sum43 = Sum43 + ambe2400_interpolate_prediction(deltal[l], prev_mp->log2Ml[intkl[l]], prev_mp->log2Ml[upper]);
    }
    Sum43 = ((MBE_AMBE2400_PREDICTION_RHO / (float)cur_mp->L) * Sum43);
#ifdef AMBE_DEBUG
    fprintf(stderr, "\n");
    fprintf(stderr, "Sum43: %f\n", Sum43);
#endif

    float Sum42 = 0;
    for (int l = 1; l <= cur_mp->L; l++) {
        Sum42 += Tl[l];
    }
    Sum42 = Sum42 / (float)cur_mp->L;
    float BigGamma = cur_mp->gamma - (0.5f * log2f((float)cur_mp->L)) - Sum42;

    for (int l = 1; l <= cur_mp->L; l++) {
        int upper = intkl[l] + 1;
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }
        float c1 = (MBE_AMBE2400_PREDICTION_RHO * ((float)1 - deltal[l]) * prev_mp->log2Ml[intkl[l]]);
        float c2 = (MBE_AMBE2400_PREDICTION_RHO * deltal[l] * prev_mp->log2Ml[upper]);
        cur_mp->log2Ml[l] = ambe2400_reconstruct_log2Ml(Tl[l], c1, c2, Sum43, BigGamma);
        if (cur_mp->Vl[l] == 1) {
            cur_mp->Ml[l] = exp2f(cur_mp->log2Ml[l]);
        } else {
            cur_mp->Ml[l] = unvc * exp2f(cur_mp->log2Ml[l]);
        }
#ifdef AMBE_DEBUG
        fprintf(stderr, "flokl[%i]: %f, intkl[%i]: %i ", l, flokl[l], l, intkl[l]);
        fprintf(stderr, "deltal[%i]: %f ", l, deltal[l]);
        fprintf(stderr, "prev_mp->log2Ml[%i]: %f\n", l, prev_mp->log2Ml[intkl[l]]);
        fprintf(stderr, "BigGamma: %f c1: %f c2: %f Sum43: %f Tl[%i]: %f log2Ml[%i]: %f Ml[%i]: %f\n", BigGamma, c1, c2,
                Sum43, l, Tl[l], l, cur_mp->log2Ml[l], l, cur_mp->Ml[l]);
#endif
    }
}

/**
 * @brief Decode AMBE 2400 parameters from demodulated bitstream.
 * @param ambe_d  Demodulated AMBE parameter bits (49).
 * @param cur_mp  Output: current frame parameters.
 * @param prev_mp Input: previous frame parameters (for prediction).
 * @return 0 for voice, a tone index (5..122, 128..147) for a tone frame, 2 for
 *         a frame that plays as a lost frame (b0 125..127: repeat the previous
 *         frame); only voice decodes cur_mp.
 */
int
mbe_decodeAmbe2400Parms(const char* ambe_d, mbe_parms* cur_mp, mbe_parms* prev_mp) {

    int L = 0;
    int Ji[5];
    float f0;
    float Cik[5][18];
    float Tl[57] = {0};
    float Ri[9];
    int ret;

    if (!cur_mp || !prev_mp) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_validate_bits(ambe_d, 49u);
    if (ret < 0) {
        return ret;
    }

#ifdef AMBE_DEBUG
    fprintf(stderr, "\n");
#endif

    int b0 = ambe2400_decode_b0(ambe_d);
    int special = ambe2400_special_frame(ambe_d, b0);
    if (special != 0) {
        return special;
    }

    ambe2400_setup_voice_model(cur_mp, b0, &L, &f0);
    float unvc = (float)0.2046 / sqrtf(cur_mp->w0);

    ambe2400_decode_vuv(ambe_d, cur_mp, b0, L, f0);
    ambe2400_decode_gain(ambe_d, cur_mp);
    ambe2400_decode_ri(ambe_d, Ri);
    ambe2400_decode_cik(ambe_d, L, Ri, Cik, Ji);
    mbe_ambe2400_inverse_dct_tl(Cik, Ji, Tl);
    mbe_ambe2400_update_spectral_amplitudes(cur_mp, prev_mp, Tl, unvc);

    return 0;
}

/**
 * @brief Demodulate interleaved AMBE 3600x2400 data in-place.
 * @param ambe_fr Frame as 4x24 bitplanes (modified).
 */
int
mbe_demodulateAmbe3600x2400Data(char ambe_fr[4][24]) {
    int ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }
    mbe_demodulateAmbe3600Data_common(ambe_fr);
    return 0;
}

int
mbe_decodeAmbe3600x2400Frame(const char ambe_fr[4][24], char ambe_d[49], mbe_process_result* result) {
    char fr[4][24];
    int c0_errors;
    int protected_errors;
    int ret;

    if (result) {
        mbe_initProcessResult(result);
    }
    if (!ambe_d) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }

    memcpy(fr, ambe_fr, sizeof(fr));
    c0_errors = mbe_eccAmbe3600x2400C0(fr);
    ret = mbe_demodulateAmbe3600x2400Data(fr);
    if (ret < 0) {
        return ret;
    }
    protected_errors = mbe_eccAmbe3600x2400Data(fr, ambe_d);

    if (result) {
        result->c0_errors = c0_errors;
        result->protected_errors = protected_errors;
        result->total_errors = c0_errors + protected_errors;
        result->flags = MBE_PROCESS_FLAG_C0_VALID;
    }
    return c0_errors + protected_errors;
}

int
mbe_decodeAmbe3600x2400SoftFrame(const mbe_soft_bit ambe_fr[4][24], char ambe_d[49], mbe_process_result* result) {
    mbe_soft_bit fr[4][24];
    int c0_errors;
    int protected_errors;
    int ret;

    if (result) {
        mbe_initProcessResult(result);
    }
    if (!ambe_d) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_validate_soft_bits((const mbe_soft_bit*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }

    memcpy(fr, ambe_fr, sizeof(fr));
    c0_errors = mbe_eccAmbe3600C0Soft_common(fr);
    mbe_demodulateAmbe3600DataSoft_common(fr);
    protected_errors = mbe_eccAmbe3600DataSoft_common(fr, ambe_d);

    if (result) {
        result->c0_errors = c0_errors;
        result->protected_errors = protected_errors;
        result->total_errors = c0_errors + protected_errors;
        result->flags = MBE_PROCESS_FLAG_SOFT_INPUT | MBE_PROCESS_FLAG_C0_VALID;
    }
    return c0_errors + protected_errors;
}

static int
ambe2400_prepare_process(mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                         mbe_parms* prev_mp_enhanced, int* total_errors, int* c0_errors) {
    int ret;

    if (!cur_mp || !prev_mp || !prev_mp_enhanced || !total_errors || !c0_errors) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_result_resolve_total_errors(result, total_errors);
    if (ret < 0) {
        return ret;
    }
    ret = mbe_validate_bits(ambe_d, 49u);
    if (ret < 0) {
        return ret;
    }

    *c0_errors = ((result->flags & MBE_PROCESS_FLAG_C0_VALID) != 0u) ? result->c0_errors : 0;
    mbe_result_prepare_synthesis(result, *total_errors);

    /* Normalize generic (IMBE) init state to the initial AMBE model. */
    mbe_ensureAmbeDefaults_common(cur_mp, prev_mp, prev_mp_enhanced);

    /* Set AMBE-specific muting threshold (9.6% vs IMBE's 8.75%). */
    cur_mp->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;

    cur_mp->errorCountTotal = *total_errors;
    cur_mp->errorCount4 = 0; /* AMBE has no Hamming cosets */
    cur_mp->errorRate = (0.95f * prev_mp->errorRate) + (0.001064f * (float)cur_mp->errorCountTotal);
    return 0;
}

static int
ambe2400_is_tone(int bad) {
    return (bad >= 5 && bad <= 122) || (bad >= 128 && bad <= 147);
}

/* Hard-decision output mutes while the error rate is above the threshold;
 * soft decisions do not mute on it (docs/dvsi-parity.md). */
static int
ambe2400_error_rate_mutes(const mbe_process_result* result, const mbe_parms* cur_mp) {
    return ((result->flags & MBE_PROCESS_FLAG_SOFT_INPUT) == 0u) && (cur_mp->errorRate > MBE_AMBE2400_MUTE_ERROR_RATE);
}

static void
ambe2400_update_decode_state(const char ambe_d[49], int bad, int total_errors, mbe_process_result* result,
                             mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    if (bad == MBE_AMBE2400_FRAME_LOST) {
        /* Repeated like a frame with too many errors, keeping this frame's own
         * error accounting and noise state. A lost tone frame is flagged as both. */
        const unsigned tone = (ambe2400_decode_b0(ambe_d) >= 126) ? MBE_PROCESS_FLAG_TONE : 0u;
        mbe_repeat_load_model(cur_mp, prev_mp);
        cur_mp->repeatCount = mbe_repeat_next_count(prev_mp->repeatCount);
        mbe_result_set_flag(result, MBE_PROCESS_FLAG_ERASURE | MBE_PROCESS_FLAG_REPEAT | tone);
        return;
    }
    if (ambe2400_is_tone(bad)) {
        /* A played tone ends a run of repeats. */
        mbe_result_set_flag(result, MBE_PROCESS_FLAG_TONE);
        cur_mp->repeatCount = 0;
        return;
    }
    if (total_errors > 3) {
        /* Repeat the previous model; this frame's error accounting and noise state continue. */
        mbe_repeat_load_model(cur_mp, prev_mp);
        cur_mp->repeatCount = mbe_repeat_next_count(prev_mp->repeatCount);
        mbe_result_set_flag(result, MBE_PROCESS_FLAG_REPEAT);
        return;
    }

    cur_mp->repeatCount = 0;
}

/* The noise D-STAR mutes to, at the level of DVSI's muted output. */
static void
ambe2400_mute_noise(float* aout_buf) {
    mbe_synthesizeUniformNoisef(aout_buf, MBE_AMBE2400_MUTE_NOISE_AMPLITUDE * MBE_SPEECH_OUTPUT_GAIN);
}

/*
 * The repeat limit's mute (the 4th consecutive repeat). DVSI's decoder fades
 * the last frame out rather than cutting it, then plays its mute noise, and
 * keeps the prediction history: the voice after five lost frames follows
 * DVSI's best this way (docs/dvsi-parity.md). The last synthesized frame
 * becomes a silent copy of itself, so the voice after the mute fades in.
 */
static void
ambe2400_mute_repeat_limit(float* aout_buf, const mbe_process_result* result, mbe_parms* cur_mp, mbe_parms* prev_mp,
                           mbe_parms* prev_mp_enhanced) {
    mbe_parms silent = *prev_mp_enhanced;
    float noise[160];

    mbe_moveMbeParms(cur_mp, prev_mp);
    memset(silent.Ml, 0, sizeof(silent.Ml));
    silent.repeatCount = 0;
    mbe_synthesizeRepeatedSpeechf(aout_buf, &silent, prev_mp_enhanced);
    /* The next frame's synthesis draws its unvoiced noise from cur_mp, as
     * after any synthesized frame. */
    cur_mp->noiseSeed = silent.noiseSeed;
    memcpy(cur_mp->noiseOverlap, silent.noiseOverlap, sizeof(cur_mp->noiseOverlap));
    mbe_moveMbeParms(&silent, prev_mp_enhanced);
    mbe_copy_error_state(prev_mp_enhanced, cur_mp);
    ambe2400_mute_noise(noise);
    /* While the error rate mutes, the last frames' speech was never played:
     * the fade still advances the state but stays silent. */
    const int audible = !ambe2400_error_rate_mutes(result, cur_mp);
    for (int i = 0; i < 160; ++i) {
        aout_buf[i] = (audible ? aout_buf[i] : 0.0f) + noise[i];
    }
}

/*
 * Tone transitions, as on DVSI's decoder: entering a tone, the last frame's
 * speech, synthesized as a repeat of it, fades out linearly over the first
 * MBE_AMBE2400_TONE_FADE samples while the tone fades in; leaving it, the tone
 * goes on and fades out over the same span while the voice fades in. Vl[0]
 * (Vl is indexed from 1, so element 0 is otherwise unused) of the state a frame
 * commits to prev_mp records the tone it leaves sounding, whose phases that
 * state carries, and whether the error rate muted the frame: the next frame
 * may be decoded from the other kind of decisions, which mute differently.
 */
#define MBE_AMBE2400_TONE_FADE   80
#define MBE_AMBE2400_TONE_HELD   (1 << 16)
#define MBE_AMBE2400_FRAME_MUTED (1 << 17)

static void
ambe2400_hold_tone(mbe_parms* mp, int index, int volume, int muted) {
    mp->Vl[0] = MBE_AMBE2400_TONE_HELD | (muted ? MBE_AMBE2400_FRAME_MUTED : 0) | (volume << 8) | index;
}

/* Whether the error rate muted the last frame. */
static int
ambe2400_last_frame_muted(const mbe_parms* prev_mp) {
    return (prev_mp->Vl[0] & MBE_AMBE2400_FRAME_MUTED) != 0;
}

/* The tone the last frame left sounding, or 0. */
static int
ambe2400_held_tone(const mbe_parms* prev_mp, int* index, int* volume) {
    const int held = prev_mp->Vl[0];
    if ((held & MBE_AMBE2400_TONE_HELD) == 0) {
        return 0;
    }
    *index = held & 0xFF;
    *volume = (held >> 8) & 0xFF;
    return 1;
}

/* Gain of the incoming signal at sample n of a transition frame. */
static float
ambe2400_fade_in_gain(int n) {
    return (n >= MBE_AMBE2400_TONE_FADE) ? 1.0f : (float)n / (float)MBE_AMBE2400_TONE_FADE;
}

/* Whether the last frame played speech: not muted by the error rate or at the
 * repeat limit (which leaves a silent last frame). */
static int
ambe2400_last_frame_audible(const mbe_parms* prev_mp) {
    return !ambe2400_last_frame_muted(prev_mp) && prev_mp->repeatCount < MBE_MAX_FRAME_REPEATS;
}

/* Under the first tone frame, the last frame repeated, as a lost frame plays
 * it (the prediction history is unchanged), to be faded out. */
static void
ambe2400_speech_under_tone(float* speech, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_parms repeated = *cur_mp;
    mbe_repeat_load_model(&repeated, prev_mp);
    repeated.repeatCount = 0;
    float pre_enh_rm0 = mbe_spectralAmpEnhanceWithRm0(&repeated, MBE_HIGH_BAND_GAIN_AMBE2400);
    mbe_synthesizeSpeechWithPreEnhRm0f(speech, &repeated, prev_mp_enhanced, pre_enh_rm0, MBE_MUTE_NOISE_COMFORT);
    /* The next frame's synthesis draws its unvoiced noise from cur_mp. */
    cur_mp->noiseSeed = repeated.noiseSeed;
    memcpy(cur_mp->noiseOverlap, repeated.noiseOverlap, sizeof(cur_mp->noiseOverlap));
    mbe_moveMbeParms(&repeated, prev_mp_enhanced);
    mbe_copy_error_state(prev_mp_enhanced, cur_mp);
}

static void
ambe2400_synthesize_tone(float* aout_buf, mbe_process_result* result, const char ambe_d[49], int index,
                         mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    int held_index, held_volume;
    const int entering = !ambe2400_held_tone(prev_mp, &held_index, &held_volume);
    const int volume = mbe_tone_dstar_volume(ambe_d);
    float speech[160] = {0};

    if (entering && ambe2400_last_frame_audible(prev_mp)) {
        ambe2400_speech_under_tone(speech, cur_mp, prev_mp, prev_mp_enhanced);
    }
    mbe_renderToneDstarf(aout_buf, cur_mp, index, volume);
    if (entering) {
        for (int n = 0; n < 160; ++n) {
            const float gain = ambe2400_fade_in_gain(n);
            aout_buf[n] = (gain * aout_buf[n]) + ((1.0f - gain) * speech[n]);
        }
    }
    /* The error rate mutes a tone as it mutes voice, after synthesis, so the
     * tone state advances: on DVSI's dam_e10_hd vector a channel error that
     * reads as a tone stays silent there. */
    const int muted = ambe2400_error_rate_mutes(result, cur_mp);
    if (muted) {
        mbe_result_set_flag(result, MBE_PROCESS_FLAG_MUTE);
        ambe2400_mute_noise(aout_buf);
    }
    ambe2400_hold_tone(cur_mp, index, volume, muted);
    mbe_moveMbeParms(cur_mp, prev_mp);
}

/* Leaving a tone: the voice fades in under the tone's fading tail (silent if
 * the error rate muted the tone), which goes on from the tone's phases. */
struct ambe2400_tone_tail {
    int index;
    int volume;
    int audible;
    int swn;
    uint32_t tone_phase;
};

static void
ambe2400_fade_in_from_tone(float* aout_buf, const struct ambe2400_tone_tail* held) {
    float tail[160] = {0};
    if (held->audible) {
        mbe_parms tone;
        memset(&tone, 0, sizeof(tone));
        tone.swn = held->swn;
        tone.tonePhase = held->tone_phase;
        mbe_renderToneDstarf(tail, &tone, held->index, held->volume);
    }
    for (int n = 0; n < 160; ++n) {
        const float gain = ambe2400_fade_in_gain(n);
        aout_buf[n] = (gain * aout_buf[n]) + ((1.0f - gain) * tail[n]);
    }
}

static void
ambe2400_synthesize_voice(float* aout_buf, mbe_process_result* result, mbe_parms* cur_mp, mbe_parms* prev_mp,
                          mbe_parms* prev_mp_enhanced) {
    if (cur_mp->repeatCount < MBE_MAX_FRAME_REPEATS) {
        mbe_moveMbeParms(cur_mp, prev_mp);
        float pre_enh_rm0 = mbe_spectralAmpEnhanceWithRm0(cur_mp, MBE_HIGH_BAND_GAIN_AMBE2400);
        mbe_synthesizeSpeechWithPreEnhRm0f(aout_buf, cur_mp, prev_mp_enhanced, pre_enh_rm0, MBE_MUTE_NOISE_COMFORT);
        mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
        if (ambe2400_error_rate_mutes(result, cur_mp)) {
            mbe_result_set_flag(result, MBE_PROCESS_FLAG_MUTE);
            ambe2400_mute_noise(aout_buf);
        }
        return;
    }

    mbe_result_set_flag(result, MBE_PROCESS_FLAG_MUTE);
    ambe2400_mute_repeat_limit(aout_buf, result, cur_mp, prev_mp, prev_mp_enhanced);
}

static void
ambe2400_synthesize_frame(float* aout_buf, mbe_process_result* result, const char ambe_d[49], int bad,
                          mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    if (ambe2400_is_tone(bad)) {
        ambe2400_synthesize_tone(aout_buf, result, ambe_d, bad, cur_mp, prev_mp, prev_mp_enhanced);
        return;
    }
    struct ambe2400_tone_tail held = {0, 0, 0, prev_mp->swn, prev_mp->tonePhase};
    const int leaving = ambe2400_held_tone(prev_mp, &held.index, &held.volume);
    held.audible = !ambe2400_last_frame_muted(prev_mp);
    cur_mp->Vl[0] = ambe2400_error_rate_mutes(result, cur_mp) ? MBE_AMBE2400_FRAME_MUTED : 0;
    ambe2400_synthesize_voice(aout_buf, result, cur_mp, prev_mp, prev_mp_enhanced);
    if (leaving && (result->flags & MBE_PROCESS_FLAG_MUTE) == 0u) {
        ambe2400_fade_in_from_tone(aout_buf, &held);
    }
}

int
mbe_processAmbe2400Dataf(float* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                         mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;
    int bad;
    int c0_errors;
    int total_errors;
    int ret;

    if (!result) {
        mbe_initProcessResult(&local_result);
        result = &local_result;
    }
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = ambe2400_prepare_process(result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, &total_errors, &c0_errors);
    if (ret < 0) {
        return ret;
    }

    bad = mbe_decodeAmbe2400Parms(ambe_d, cur_mp, prev_mp);
    if (bad < 0) {
        return bad;
    }

    /* A tone frame that fails the error check plays as a lost frame, as one
     * whose index has no tone does: it used to reset the decoder, clearing the
     * error rate, so output unmuted mid-burst where DVSI's stays muted. */
    if (ambe2400_is_tone(bad) && (c0_errors >= 2 || total_errors >= 3)) {
        bad = MBE_AMBE2400_FRAME_LOST;
    }
    ambe2400_update_decode_state(ambe_d, bad, total_errors, result, cur_mp, prev_mp);
    ambe2400_synthesize_frame(aout_buf, result, ambe_d, bad, cur_mp, prev_mp, prev_mp_enhanced);
    return result->total_errors;
}

int
mbe_processAmbe2400Data(short* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                        mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processAmbe2400Dataf(float_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete AMBE 3600x2400 frame into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional output status populated by decode and synthesis.
 * @param ambe_fr  Input frame as 4x24 bitplanes.
 * @param ambe_d   Scratch/output parameter bits (49).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 */
int
mbe_processAmbe3600x2400Framef(float* aout_buf, mbe_process_result* result, const char ambe_fr[4][24], char ambe_d[49],
                               mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;
    int ret;
    if (!result) {
        result = &local_result;
    }
    ret = mbe_decodeAmbe3600x2400Frame(ambe_fr, ambe_d, result);
    if (ret < 0) {
        return ret;
    }
    return mbe_processAmbe2400Dataf(aout_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processAmbe3600x2400SoftFramef(float* aout_buf, mbe_process_result* result, const mbe_soft_bit ambe_fr[4][24],
                                   char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                   mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;
    int ret;
    if (!result) {
        result = &local_result;
    }
    ret = mbe_decodeAmbe3600x2400SoftFrame(ambe_fr, ambe_d, result);
    if (ret < 0) {
        return ret;
    }
    return mbe_processAmbe2400Dataf(aout_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processAmbe3600x2400SoftFrame(short* aout_buf, mbe_process_result* result, const mbe_soft_bit ambe_fr[4][24],
                                  char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processAmbe3600x2400SoftFramef(float_buf, result, ambe_fr, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete AMBE 3600x2400 frame into 16-bit PCM.
 * @see mbe_processAmbe3600x2400Framef for details.
 */
int
mbe_processAmbe3600x2400Frame(short* aout_buf, mbe_process_result* result, const char ambe_fr[4][24], char ambe_d[49],
                              mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;

    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processAmbe3600x2400Framef(float_buf, result, ambe_fr, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}
