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
 * @brief AMBE 3600x2450 parameter decode, ECC, and synthesis hooks.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ambe3600x2450_const.h"
#include "ambe_common.h"
#include "mbe_adaptive.h"
#include "mbe_compiler.h"
#include "mbe_result.h"
#include "mbe_tone.h"
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
struct ambe2450_dct_cache {
    int inited;
    float ri_cos[9][9];         /* [m][i] for m=1..8, i=1..8 (index 0 unused) */
    float idct_cos[18][18][18]; /* [ji][j][k] for ji=1..17, j=1..ji, k=1..ji */
};

static MBE_THREAD_LOCAL struct ambe2450_dct_cache ambe2450_cache = {0};

/**
 * @brief Initialize or return the thread-local AMBE 2450 DCT cache.
 *
 * Fills the cosine tables on first use. Because the cache is thread-local,
 * no locking is required.
 *
 * @return Pointer to the initialized cache.
 */
static struct ambe2450_dct_cache*
ambe2450_get_dct_cache(void) {
    if (ambe2450_cache.inited) {
        return &ambe2450_cache;
    }

    /* Fill Ri DCT cosine table: cosf(M_PI * (m-1) * (i-0.5) / 8) */
    for (int m = 1; m <= 8; m++) {
        for (int i = 1; i <= 8; i++) {
            ambe2450_cache.ri_cos[m][i] = cosf((M_PI * (float)(m - 1) * ((float)i - 0.5f)) / 8.0f);
        }
    }

    /* Fill per-block IDCT cosine table: cosf(M_PI * (k-1) * (j-0.5) / ji) */
    for (int ji = 1; ji <= 17; ji++) {
        for (int j = 1; j <= ji; j++) {
            for (int k = 1; k <= ji; k++) {
                ambe2450_cache.idct_cos[ji][j][k] = cosf((M_PI * (float)(k - 1) * ((float)j - 0.5f)) / (float)ji);
            }
        }
    }

    ambe2450_cache.inited = 1;
    return &ambe2450_cache;
}

static int
ambe2450_is_valid_tone_id(const char ambe_d[49]) {
    int id1 = 0;
    /* AMBE tone ID1 is U1[0..7] (bits 12..19 in ambe_d), not the full 12-bit U1 field. */
    for (int i = 12; i < 20; i++) {
        id1 = (id1 << 1) | (int)ambe_d[i];
    }

    return mbe_tone_id_is_valid(id1);
}

/**
 * @brief Print AMBE 2450 parameter bits to stderr (debug aid).
 * @param ambe_d AMBE parameter bits (49).
 */
void
mbe_dumpAmbe2450Data(const char* ambe_d) {

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
 * @brief Print raw AMBE 3600x2450 frame bitplanes to stderr.
 * @param ambe_fr Frame as 4x24 bitplanes.
 */
void
mbe_dumpAmbe3600x2450Frame(const char ambe_fr[4][24]) {

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
 * @brief Apply ECC to AMBE 3600x2450 C0 and update in-place.
 * @param ambe_fr Frame as 4x24 bitplanes.
 * @return Number of corrected errors in C0.
 */
int
mbe_eccAmbe3600x2450C0(char ambe_fr[4][24]) {
    return mbe_eccAmbe3600C0_common(ambe_fr);
}

/**
 * @brief Apply ECC to AMBE 3600x2450 data and pack parameter bits.
 * @param ambe_fr Frame as 4x24 bitplanes.
 * @param ambe_d  Output parameter bits (49).
 * @return Number of corrected errors in protected fields.
 */
int
mbe_eccAmbe3600x2450Data(char ambe_fr[4][24], char* ambe_d) {
    return mbe_eccAmbe3600Data_common(ambe_fr, ambe_d);
}

static void
ambe2450_read_u_fields(const char* ambe_d, int* u0, int* u1, int* u2, int* u3) {
    *u0 = 0;
    *u1 = 0;
    *u2 = 0;
    *u3 = 0;

    for (int i = 0; i < 12; i++) {
        *u0 = (*u0 << 1) | (int)ambe_d[i];
    }
    for (int i = 12; i < 24; i++) {
        *u1 = (*u1 << 1) | (int)ambe_d[i];
    }
    for (int i = 24; i < 35; i++) {
        *u2 = (*u2 << 1) | (int)ambe_d[i];
    }
    for (int i = 35; i < 49; i++) {
        *u3 = (*u3 << 1) | (int)ambe_d[i];
    }
}

static void
ambe2450_decode_vuv(const char* ambe_d, mbe_parms* cur_mp, int L, float f0, int silence, int b0) {
    (void)b0;
    int b1 = 0;
    b1 |= ambe_d[4] << 4;
    b1 |= ambe_d[5] << 3;
    b1 |= ambe_d[6] << 2;
    b1 |= ambe_d[7] << 1;
    b1 |= ambe_d[35];

    for (int l = 1; l <= L; l++) {
        int jl = (int)((float)l * (float)16.0 * f0);
        if (silence == 0) {
            cur_mp->Vl[l] = AmbeVuv[b1][jl];
        }
#ifdef AMBE_DEBUG
        fprintf(stderr, "jl[%i]:%i Vl[%i]:%i\n", l, jl, l, cur_mp->Vl[l]);
#endif
    }
#ifdef AMBE_DEBUG
    fprintf(stderr, "\nb0:%i w0:%f L:%i b1:%i\n", b0, cur_mp->w0, L, b1);
#endif
}

static void
ambe2450_decode_ri(const char* ambe_d, float Ri[9]) {
    float Gm[9];
    Gm[1] = 0;

    int b3 = 0;
    b3 |= ambe_d[12] << 8;
    b3 |= ambe_d[13] << 7;
    b3 |= ambe_d[14] << 6;
    b3 |= ambe_d[15] << 5;
    b3 |= ambe_d[16] << 4;
    b3 |= ambe_d[17] << 3;
    b3 |= ambe_d[18] << 2;
    b3 |= ambe_d[19] << 1;
    b3 |= ambe_d[40];
    Gm[2] = AmbePRBA24[b3][0];
    Gm[3] = AmbePRBA24[b3][1];
    Gm[4] = AmbePRBA24[b3][2];

    int b4 = 0;
    b4 |= ambe_d[20] << 6;
    b4 |= ambe_d[21] << 5;
    b4 |= ambe_d[22] << 4;
    b4 |= ambe_d[23] << 3;
    b4 |= ambe_d[41] << 2;
    b4 |= ambe_d[42] << 1;
    b4 |= ambe_d[43];
    Gm[5] = AmbePRBA58[b4][0];
    Gm[6] = AmbePRBA58[b4][1];
    Gm[7] = AmbePRBA58[b4][2];
    Gm[8] = AmbePRBA58[b4][3];

#ifdef AMBE_DEBUG
    fprintf(stderr, "b3: %i Gm[2]: %f Gm[3]: %f Gm[4]: %f b4: %i Gm[5]: %f Gm[6]: %f Gm[7]: %f Gm[8]: %f\n", b3, Gm[2],
            Gm[3], Gm[4], b4, Gm[5], Gm[6], Gm[7], Gm[8]);
#endif

    struct ambe2450_dct_cache* cache = ambe2450_get_dct_cache();
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
ambe2450_decode_cik(const char* ambe_d, int L, const float Ri[9], float Cik[5][18], int Ji[5]) {
    const float rconst = ((float)1 / ((float)2 * M_SQRT2));
    Cik[1][1] = (float)0.5 * (Ri[1] + Ri[2]);
    Cik[1][2] = rconst * (Ri[1] - Ri[2]);
    Cik[2][1] = (float)0.5 * (Ri[3] + Ri[4]);
    Cik[2][2] = rconst * (Ri[3] - Ri[4]);
    Cik[3][1] = (float)0.5 * (Ri[5] + Ri[6]);
    Cik[3][2] = rconst * (Ri[5] - Ri[6]);
    Cik[4][1] = (float)0.5 * (Ri[7] + Ri[8]);
    Cik[4][2] = rconst * (Ri[7] - Ri[8]);

    int b5 = 0;
    b5 |= ambe_d[24] << 4;
    b5 |= ambe_d[25] << 3;
    b5 |= ambe_d[26] << 2;
    b5 |= ambe_d[27] << 1;
    b5 |= ambe_d[44];

    int b6 = 0;
    b6 |= ambe_d[28] << 3;
    b6 |= ambe_d[29] << 2;
    b6 |= ambe_d[30] << 1;
    b6 |= ambe_d[45];

    int b7 = 0;
    b7 |= ambe_d[31] << 3;
    b7 |= ambe_d[32] << 2;
    b7 |= ambe_d[33] << 1;
    b7 |= ambe_d[46];

    int b8 = 0;
    b8 |= ambe_d[34] << 2;
    b8 |= ambe_d[47] << 1;
    b8 |= ambe_d[48];

    Ji[1] = AmbeLmprbl[L][0];
    Ji[2] = AmbeLmprbl[L][1];
    Ji[3] = AmbeLmprbl[L][2];
    Ji[4] = AmbeLmprbl[L][3];
#ifdef AMBE_DEBUG
    fprintf(stderr, "Ji[1]: %i Ji[2]: %i Ji[3]: %i Ji[4]: %i\n", Ji[1], Ji[2], Ji[3], Ji[4]);
    fprintf(stderr, "b5: %i b6: %i b7: %i b8: %i\n", b5, b6, b7, b8);
#endif

    for (int k = 3; k <= Ji[1]; k++) {
        if (k > 6) {
            Cik[1][k] = 0;
        } else {
            Cik[1][k] = AmbeHOCb5[b5][k - 3];
#ifdef AMBE_DEBUG
            fprintf(stderr, "C1,%i: %f ", k, Cik[1][k]);
#endif
        }
    }
    for (int k = 3; k <= Ji[2]; k++) {
        if (k > 6) {
            Cik[2][k] = 0;
        } else {
            Cik[2][k] = AmbeHOCb6[b6][k - 3];
#ifdef AMBE_DEBUG
            fprintf(stderr, "C2,%i: %f ", k, Cik[2][k]);
#endif
        }
    }
    for (int k = 3; k <= Ji[3]; k++) {
        if (k > 6) {
            Cik[3][k] = 0;
        } else {
            Cik[3][k] = AmbeHOCb7[b7][k - 3];
#ifdef AMBE_DEBUG
            fprintf(stderr, "C3,%i: %f ", k, Cik[3][k]);
#endif
        }
    }
    for (int k = 3; k <= Ji[4]; k++) {
        if (k > 6) {
            Cik[4][k] = 0;
        } else {
            Cik[4][k] = AmbeHOCb8[b8][k - 3];
#ifdef AMBE_DEBUG
            fprintf(stderr, "C4,%i: %f ", k, Cik[4][k]);
#endif
        }
    }
#ifdef AMBE_DEBUG
    fprintf(stderr, "\n");
#endif
}

static void
ambe2450_inverse_dct_tl(float Cik[5][18], int Ji[5], float Tl[57]) {
    struct ambe2450_dct_cache* cache = ambe2450_get_dct_cache();
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

static void
ambe2450_update_spectral_amplitudes(mbe_parms* cur_mp, mbe_parms* prev_mp, const float Tl[57], float unvc) {
    int intkl[57];
    float flokl[57], deltal[57];
    int prev_L = prev_mp->L;
    if (cur_mp->L < 1) {
        cur_mp->L = 1;
    } else if (cur_mp->L > 56) {
        cur_mp->L = 56;
    }
    if (prev_L < 1) {
        prev_L = 1;
    } else if (prev_L > 56) {
        prev_L = 56;
    }

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
        flokl[l] = ((float)prev_L / (float)cur_mp->L) * (float)l;
        intkl[l] = (int)(flokl[l]);
#ifdef AMBE_DEBUG
        fprintf(stderr, "flok%i: %f, intk%i: %i ", l, flokl[l], l, intkl[l]);
#endif
        deltal[l] = flokl[l] - (float)intkl[l];
#ifdef AMBE_DEBUG
        fprintf(stderr, "delta%i: %f ", l, deltal[l]);
#endif
        Sum43 = Sum43
                + ((((float)1 - deltal[l]) * prev_mp->log2Ml[intkl[l]]) + (deltal[l] * prev_mp->log2Ml[intkl[l] + 1]));
    }
    Sum43 = (((float)0.65 / (float)cur_mp->L) * Sum43);
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
        float c1 = ((float)0.65 * ((float)1 - deltal[l]) * prev_mp->log2Ml[intkl[l]]);
        float c2 = ((float)0.65 * deltal[l] * prev_mp->log2Ml[intkl[l] + 1]);
        cur_mp->log2Ml[l] = Tl[l] + c1 + c2 - Sum43 + BigGamma;
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

static int
ambe2450_decode_b0(const char* ambe_d) {
    int b0 = 0;
    b0 |= ambe_d[0] << 6;
    b0 |= ambe_d[1] << 5;
    b0 |= ambe_d[2] << 4;
    b0 |= ambe_d[3] << 3;
    b0 |= ambe_d[37] << 2;
    b0 |= ambe_d[38] << 1;
    b0 |= ambe_d[39];
    return b0;
}

static int
ambe2450_tone_verified(const char* ambe_d) {
    int u0, u1, u2, u3;
    ambe2450_read_u_fields(ambe_d, &u0, &u1, &u2, &u3);
    (void)u2;

    int u0_tone_check = (u0 >> 6) & 0x3f;
    int u3_tone_check = (u3 & 0xf);
    int u1_high_tone_verify = (u1 >> 8) & 0xf;
    int u1_low_tone_verify = u1 & 0xf;

#ifdef AMBE_DEBUG
    fprintf(stderr, "TONECHK u0=%d u3=%d u1h=%d u1l=%d ", u0_tone_check, u3_tone_check, u1_high_tone_verify,
            u1_low_tone_verify);
#endif

    return (u0_tone_check == 63) && ((u3_tone_check == 0) || (u1_high_tone_verify == u1_low_tone_verify));
}

static void
ambe2450_set_silence_model(mbe_parms* cur_mp, int b0, int* L, float* f0, int* silence) {
    *silence = 1;
    /* JMBE AMBEFundamentalFrequency.W124/W125: constructor uses (frequency * 2*PI), where frequency is PI/32. */
    *f0 = (float)M_PI / 32.0f;
    cur_mp->w0 = *f0 * (float)(2.0 * M_PI);
    *L = (b0 == 124) ? 15 : 14; /* JMBE maps W124->L15, W125->L14 */
    cur_mp->L = *L;
    for (int l = 1; l <= *L; l++) {
        cur_mp->Vl[l] = 0;
    }
}

static int
ambe2450_setup_frame_model(const char* ambe_d, mbe_parms* cur_mp, int total_errors, int* b0, int* L, float* f0,
                           int* silence) {
    *silence = 0;

    /* JMBE-compatible tone classification:
     * - Tone if U0 tone check passes and either U3 check is zero or U1 high/low nibbles match.
     * - If total BER is known and >= 6, do not classify as tone (JMBE AMBEFrame). */
    if (ambe2450_tone_verified(ambe_d) && total_errors < 6) {
#ifdef AMBE_DEBUG
        fprintf(stderr, "Tone Frame\n");
#endif
        return 7;
    }

    *b0 = ambe2450_decode_b0(ambe_d);

    if ((*b0 >= 120) && (*b0 <= 123)) {
#ifdef AMBE_DEBUG
        fprintf(stderr, "Erasure Frame b0 = %d\n", *b0);
#endif
        return 2;
    }
    if ((*b0 == 124) || (*b0 == 125)) {
#ifdef AMBE_DEBUG
        fprintf(stderr, "Silence Frame\n");
#endif
        ambe2450_set_silence_model(cur_mp, *b0, L, f0, silence);
        return 0;
    }
    /* If the fundamental decodes as a tone but tone verification failed above,
     * treat this as erasure to match JMBE's fallback behavior. */
    if ((*b0 == 126) || (*b0 == 127)) {
#ifdef AMBE_DEBUG
        fprintf(stderr, "Unverified tone fundamental -> erasure\n");
#endif
        return 2;
    }
    if ((*b0 < 0) || (*b0 >= 120)) {
        return 2;
    }

    *f0 = AmbeW0table[*b0];
    cur_mp->w0 = *f0 * (float)2 * M_PI;
    *L = AmbeLtable[*b0];
    cur_mp->L = *L;
    return 0;
}

/**
 * @brief Internal AMBE 2450 parameter decode with optional tone BER gate.
 *
 * @param ambe_d  Demodulated AMBE parameter bits (49).
 * @param cur_mp  Output: current frame parameters.
 * @param prev_mp Input: previous frame parameters (for prediction).
 * @param total_errors Total frame error count for JMBE tone gating, or <0 to disable.
 * @return Tone index or 0 for voice; implementation-specific non-zero for special frames.
 */
static int
mbe_decodeAmbe2450ParmsInternal(const char* ambe_d, mbe_parms* cur_mp, mbe_parms* prev_mp, int total_errors) {

    int b0, b2, L = 0;
    float f0, Cik[5][18], Tl[57] = {0}, Ri[9];
    int silence = 0;
    int Ji[5];
    float deltaGamma;
    float unvc;

#ifdef AMBE_DEBUG
    fprintf(stderr, "\n");
#endif

    // copy repeat from prev_mp
    cur_mp->repeat = prev_mp->repeat;

    int special_frame = ambe2450_setup_frame_model(ambe_d, cur_mp, total_errors, &b0, &L, &f0, &silence);
    if (special_frame != 0) {
        return special_frame;
    }

    unvc = (float)0.2046 / sqrtf(cur_mp->w0);

    // decode V/UV parameters
    ambe2450_decode_vuv(ambe_d, cur_mp, L, f0, silence, b0);

    // decode gain vector
    // load b2 from ambe_d
    b2 = 0;
    b2 |= ambe_d[8] << 4;
    b2 |= ambe_d[9] << 3;
    b2 |= ambe_d[10] << 2;
    b2 |= ambe_d[11] << 1;
    b2 |= ambe_d[36];

    deltaGamma = AmbeDg[b2];
    cur_mp->gamma = deltaGamma + ((float)0.5 * prev_mp->gamma);
#ifdef AMBE_DEBUG
    fprintf(stderr, "b2: %i, deltaGamma: %f gamma: %f gamma-1: %f\n", b2, deltaGamma, cur_mp->gamma, prev_mp->gamma);
#endif

    // decode PRBA/HOC vectors and inverse DCT
    ambe2450_decode_ri(ambe_d, Ri);
    ambe2450_decode_cik(ambe_d, L, Ri, Cik, Ji);
    ambe2450_inverse_dct_tl(Cik, Ji, Tl);

    // determine log2Ml by applying ci,j to previous log2Ml
    ambe2450_update_spectral_amplitudes(cur_mp, prev_mp, Tl, unvc);

    return (0);
}

/**
 * @brief Decode AMBE 2450 parameters from demodulated bitstream.
 * @param ambe_d  Demodulated AMBE parameter bits (49).
 * @param cur_mp  Output: current frame parameters.
 * @param prev_mp Input: previous frame parameters (for prediction).
 * @return Tone index or 0 for voice; implementation-specific non-zero for tone frames.
 */
int
mbe_decodeAmbe2450Parms(const char* ambe_d, mbe_parms* cur_mp, mbe_parms* prev_mp) {
    return mbe_decodeAmbe2450ParmsInternal(ambe_d, cur_mp, prev_mp, -1);
}

/**
 * @brief Demodulate interleaved AMBE 3600x2450 data in-place.
 * @param ambe_fr Frame as 4x24 bitplanes (modified).
 */
void
mbe_demodulateAmbe3600x2450Data(char ambe_fr[4][24]) {
    mbe_demodulateAmbe3600Data_common(ambe_fr);
}

int
mbe_decodeAmbe3600x2450Frame(const char ambe_fr[4][24], char ambe_d[49], mbe_process_result* result) {
    char fr[4][24];
    int c0_errors;
    int protected_errors;

    memcpy(fr, ambe_fr, sizeof(fr));
    c0_errors = mbe_eccAmbe3600x2450C0(fr);
    mbe_demodulateAmbe3600x2450Data(fr);
    protected_errors = mbe_eccAmbe3600x2450Data(fr, ambe_d);

    if (result) {
        mbe_initProcessResult(result);
        result->c0_errors = c0_errors;
        result->protected_errors = protected_errors;
        result->total_errors = c0_errors + protected_errors;
        result->flags = MBE_PROCESS_FLAG_C0_VALID;
    }
    return c0_errors + protected_errors;
}

int
mbe_decodeAmbe3600x2450SoftFrame(const mbe_soft_bit ambe_fr[4][24], char ambe_d[49], mbe_process_result* result) {
    mbe_soft_bit fr[4][24];
    int c0_errors;
    int protected_errors;

    memcpy(fr, ambe_fr, sizeof(fr));
    c0_errors = mbe_eccAmbe3600C0Soft_common(fr);
    mbe_demodulateAmbe3600DataSoft_common(fr);
    protected_errors = mbe_eccAmbe3600DataSoft_common(fr, ambe_d);

    if (result) {
        mbe_initProcessResult(result);
        result->c0_errors = c0_errors;
        result->protected_errors = protected_errors;
        result->total_errors = c0_errors + protected_errors;
        result->flags = MBE_PROCESS_FLAG_SOFT_INPUT | MBE_PROCESS_FLAG_C0_VALID;
    }
    return c0_errors + protected_errors;
}

/**
 * @brief Internal AMBE 2450 synthesis path with optional C0 error context.
 *
 * Frame decode paths provide C0 errors for JMBE repeat criteria parity.
 * Public Dataf calls do not have C0 context and use the historical total-error fallback.
 */
static void
mbe_processAmbe2450Dataf_internal(float* aout_buf, const int* errs2, char* err_str, const char ambe_d[49],
                                  mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced, int uvquality,
                                  int c0_errors, int c0_errors_valid) {
    int i, bad;

    /* JMBE AMBE starts from W124 defaults; normalize if caller used generic init. */
    mbe_ensureAmbeDefaults_common(cur_mp, prev_mp, prev_mp_enhanced);

    /* Set AMBE-specific muting threshold (9.6% vs IMBE's 8.75%). */
    cur_mp->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;

    /* Set error metrics for adaptive smoothing (JMBE Algorithms #55-56, #111-116). */
    cur_mp->errorCountTotal = *errs2;
    cur_mp->errorCount4 = 0; /* AMBE has no Hamming cosets */
    cur_mp->errorRate = (0.95f * prev_mp->errorRate) + (0.001064f * (float)cur_mp->errorCountTotal);

    for (i = 0; i < *errs2; i++) {
        *err_str = '=';
        err_str++;
    }

    //it should be noted that in this context, 'bad' isn't referring to bad decode, but is a return
    //value for which type of frame we should synthesize (voice, repeat, silence, erasure, or tone, etc)
    bad = mbe_decodeAmbe2450ParmsInternal(ambe_d, cur_mp, prev_mp, *errs2);
    if (bad == 2) {
        // Erasure frame
        *err_str = 'E';
        err_str++;
        cur_mp->repeat = 0;
        cur_mp->repeatCount = 0;
        /* JMBE erasure model: W120 defaults carried as previous frame state. */
        mbe_setAmbeErasureParms_common(cur_mp, prev_mp);
    } else if (bad == 3 || bad == 7) {
        // Tone Frame
        *err_str = 'T';
        err_str++;
        cur_mp->repeat = 0;
        cur_mp->repeatCount = 0;
    } else {
        int repeat_required;

        if (c0_errors_valid) {
            /* JMBE AMBE repeat criteria when C0 errors are known:
             * errors[0] >= 4 || (errors[0] >= 2 && total >= 6). */
            repeat_required = ((c0_errors >= 4) || ((c0_errors >= 2) && (*errs2 >= 6)));
        } else {
            /* Dataf callers pass parameter bits only (no C0 context); keep historical total-error fallback. */
            repeat_required = (*errs2 > 3);
        }

        if (repeat_required) {
            mbe_useLastMbeParms(cur_mp, prev_mp);
            cur_mp->repeat++;
            cur_mp->repeatCount++;
            *err_str = 'R';
            err_str++;
        } else {
            cur_mp->repeat = 0;
            cur_mp->repeatCount = 0;
        }
    }

    if (bad == 0) {
        if (cur_mp->repeat <= 3) {
            mbe_moveMbeParms(cur_mp, prev_mp);
            float pre_enh_rm0 = mbe_spectralAmpEnhanceWithRm0(cur_mp);
            mbe_synthesizeSpeechWithPreEnhRm0f(aout_buf, cur_mp, prev_mp_enhanced, uvquality, pre_enh_rm0);
            mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
        } else {
            *err_str = 'M';
            err_str++;
            mbe_synthesizeComfortNoisef(aout_buf);
            mbe_initAmbeParms_common(cur_mp, prev_mp, prev_mp_enhanced);
        }
    } else if (bad == 7) {
        /* JMBE tone behavior:
         * - valid tone IDs synthesize tone audio without replacing previous voice model
         * - invalid tone IDs synthesize prior voice frame unless repeat maxed */
        if (ambe2450_is_valid_tone_id(ambe_d)) {
            mbe_synthesizeTonef(aout_buf, ambe_d, cur_mp);
        } else if (!mbe_isMaxFrameRepeat(prev_mp)) {
            mbe_parms synth_mp;
            /* JMBE invalid-tone behavior reuses the prior VOICE model (enhanced amplitudes) while still advancing synth state. */
            mbe_moveMbeParms(prev_mp_enhanced, &synth_mp);
            mbe_synthesizeSpeechf(aout_buf, &synth_mp, prev_mp_enhanced, uvquality);
            mbe_moveMbeParms(&synth_mp, prev_mp_enhanced);
        } else {
            mbe_synthesizeComfortNoisef(aout_buf);
            mbe_initAmbeParms_common(cur_mp, prev_mp, prev_mp_enhanced);
        }
    } else if (bad == 2) {
        /* JMBE erasure behavior: synthesize white noise and keep ERASURE
         * parameters as the previous-frame context for recovery. */
        mbe_synthesizeComfortNoisef(aout_buf);
        mbe_moveMbeParms(cur_mp, prev_mp);
        mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
    } else {
        mbe_synthesizeComfortNoisef(aout_buf);
        mbe_initAmbeParms_common(cur_mp, prev_mp, prev_mp_enhanced);
    }
    *err_str = 0;
}

/**
 * @brief Process AMBE 2450 parameters into 160 float samples at 8 kHz.
 * @param aout_buf Output buffer of 160 float samples.
 * @param errs     Reserved input for API compatibility (currently ignored).
 * @param errs2    Input total/protected-field error count.
 * @param err_str  Output status trace string.
 * @param ambe_d   Demodulated parameter bits (49).
 * @param cur_mp   In/out: current frame parameters (may be enhanced).
 * @param prev_mp  In/out: previous frame parameters.
 * @param prev_mp_enhanced In/out: enhanced previous parameters for continuity.
 * @param uvquality Legacy quality knob (currently ignored; kept for API compatibility).
 */
void
mbe_processAmbe2450Dataf(float* aout_buf, const int* errs, const int* errs2, char* err_str, const char ambe_d[49],
                         mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced, int uvquality) {
    (void)errs;
    mbe_processAmbe2450Dataf_internal(aout_buf, errs2, err_str, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, uvquality, 0,
                                      0);
}

/**
 * @brief Process AMBE 2450 parameters into 160 16-bit samples at 8 kHz.
 * @see mbe_processAmbe2450Dataf for parameter details.
 */
void
mbe_processAmbe2450Data(short* aout_buf, const int* errs, const int* errs2, char* err_str, const char ambe_d[49],
                        mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced, int uvquality) {
    float float_buf[160];

    mbe_processAmbe2450Dataf(float_buf, errs, errs2, err_str, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, uvquality);
    mbe_floattoshort(float_buf, aout_buf);
}

int
mbe_processAmbe2450DatafV2(float* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                           mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced, int uvquality) {
    mbe_process_result local_result;
    char status[MBE_RESULT_STATUS_SIZE];

    if (!result) {
        mbe_initProcessResult(&local_result);
        result = &local_result;
    }

    int errs2 = result->total_errors;
    if (errs2 == 0 && (result->c0_errors != 0 || result->protected_errors != 0)) {
        errs2 = result->c0_errors + result->protected_errors;
    }
    errs2 = mbe_result_clamp_status_errors(errs2, sizeof(status));

    status[0] = '\0';
    mbe_processAmbe2450Dataf_internal(aout_buf, &errs2, status, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, uvquality,
                                      result->c0_errors, (result->flags & MBE_PROCESS_FLAG_C0_VALID) != 0u);
    result->total_errors = errs2;
    result->protected_errors = errs2 - result->c0_errors;
    mbe_result_apply_status(result, status);
    return result->total_errors;
}

int
mbe_processAmbe2450DataV2(short* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                          mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced, int uvquality) {
    float float_buf[160];
    int ret = mbe_processAmbe2450DatafV2(float_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, uvquality);
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete AMBE 3600x2450 frame into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param errs     Output corrected C0 error count.
 * @param errs2    Output total/protected-field error count.
 * @param err_str  Output status trace string.
 * @param ambe_fr  Input frame as 4x24 bitplanes.
 * @param ambe_d   Scratch/output parameter bits (49).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 * @param uvquality Legacy quality knob (currently ignored; kept for API compatibility).
 */
void
mbe_processAmbe3600x2450Framef(float* aout_buf, int* errs, int* errs2, char* err_str, char ambe_fr[4][24],
                               char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced,
                               int uvquality) {

    *errs = 0;
    *errs2 = 0;
    *errs = mbe_eccAmbe3600x2450C0(ambe_fr);
    mbe_demodulateAmbe3600x2450Data(ambe_fr);
    *errs2 = *errs;
    *errs2 += mbe_eccAmbe3600x2450Data(ambe_fr, ambe_d);

    mbe_processAmbe2450Dataf_internal(aout_buf, errs2, err_str, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, uvquality,
                                      *errs, 1);
}

int
mbe_processAmbe3600x2450SoftFramef(float* aout_buf, mbe_process_result* result, const mbe_soft_bit ambe_fr[4][24],
                                   char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced,
                                   int uvquality) {
    mbe_process_result local_result;
    if (!result) {
        result = &local_result;
    }
    (void)mbe_decodeAmbe3600x2450SoftFrame(ambe_fr, ambe_d, result);
    return mbe_processAmbe2450DatafV2(aout_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, uvquality);
}

int
mbe_processAmbe3600x2450SoftFrame(short* aout_buf, mbe_process_result* result, const mbe_soft_bit ambe_fr[4][24],
                                  char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced,
                                  int uvquality) {
    float float_buf[160];
    int ret = mbe_processAmbe3600x2450SoftFramef(float_buf, result, ambe_fr, ambe_d, cur_mp, prev_mp, prev_mp_enhanced,
                                                 uvquality);
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete AMBE 3600x2450 frame into 16-bit PCM.
 * @see mbe_processAmbe3600x2450Framef for details.
 */
void
mbe_processAmbe3600x2450Frame(short* aout_buf, int* errs, int* errs2, char* err_str, char ambe_fr[4][24],
                              char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced,
                              int uvquality) {
    float float_buf[160];

    mbe_processAmbe3600x2450Framef(float_buf, errs, errs2, err_str, ambe_fr, ambe_d, cur_mp, prev_mp, prev_mp_enhanced,
                                   uvquality);
    mbe_floattoshort(float_buf, aout_buf);
}
