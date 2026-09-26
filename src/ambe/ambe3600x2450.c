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
ambe2450_read_tone_id(const char ambe_d[49]) {
    int id1 = 0;
    /* AMBE tone ID1 is U1[0..7] (bits 12..19 in ambe_d), not the full 12-bit U1 field. */
    for (int i = 12; i < 20; i++) {
        id1 = (id1 << 1) | (int)ambe_d[i];
    }
    return id1;
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
    int ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }
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
    if (!ambe_d) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    int ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }
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

    const struct ambe2450_dct_cache* cache = ambe2450_get_dct_cache();
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
ambe2450_inverse_dct_tl(float Cik[5][18], const int Ji[5], float Tl[57]) {
    const struct ambe2450_dct_cache* cache = ambe2450_get_dct_cache();
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
        flokl[l] = ((float)prev_L / (float)cur_mp->L) * (float)l;
        intkl[l] = (int)(flokl[l]);
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
        Sum43 = Sum43 + ((((float)1 - deltal[l]) * prev_mp->log2Ml[intkl[l]]) + (deltal[l] * prev_mp->log2Ml[upper]));
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
        int upper = intkl[l] + 1;
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }
        float c1 = ((float)0.65 * ((float)1 - deltal[l]) * prev_mp->log2Ml[intkl[l]]);
        float c2 = ((float)0.65 * deltal[l] * prev_mp->log2Ml[upper]);
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

/** TIA-102.BABA-1 7: the frame is a tone frame when the first six bits of u0 equal 63. */
static int
ambe2450_is_tone_frame(const char* ambe_d) {
    for (int i = 0; i < 6; i++) {
        if (ambe_d[i] != 1) {
            return 0;
        }
    }
    return 1;
}

/**
 * Consistency check on a tone frame's redundant fields (JMBE's; the spec
 * repeats the tone index but leaves checking it to the decoder, 7.2 Table 10):
 * the last four bits of u3 are 0, or u1's trailing copy of ID(7..4) matches.
 */
static int
ambe2450_tone_fields_consistent(const char* ambe_d) {
    int u0, u1, u2, u3;
    ambe2450_read_u_fields(ambe_d, &u0, &u1, &u2, &u3);
    (void)u0;
    (void)u2;

    int u3_tone_check = (u3 & 0xf);
    int u1_high_tone_verify = (u1 >> 8) & 0xf;
    int u1_low_tone_verify = u1 & 0xf;

#ifdef AMBE_DEBUG
    fprintf(stderr, "TONECHK u3=%d u1h=%d u1l=%d ", u3_tone_check, u1_high_tone_verify, u1_low_tone_verify);
#endif

    return (u3_tone_check == 0) || (u1_high_tone_verify == u1_low_tone_verify);
}

/** A verified tone frame: the tone identifier (7) plus consistent redundant fields. */
static int
ambe2450_tone_verified(const char* ambe_d) {
    return ambe2450_is_tone_frame(ambe_d) && ambe2450_tone_fields_consistent(ambe_d);
}

#define AMBE2450_TONE_ID_ZERO_AMPLITUDE 255 /* Table 9 / 7.3 */

/** 7.3: a tone index the decoder can output, including the zero-amplitude ID 255. */
static int
ambe2450_tone_id_usable(const char* ambe_d) {
    const int tone_id = ambe2450_read_tone_id(ambe_d);
    return (tone_id == AMBE2450_TONE_ID_ZERO_AMPLITUDE) || mbe_tone_id_is_valid(tone_id);
}

/**
 * Frame type from the parameter bits alone. A tone frame (7) is recognised
 * before b0: its other bits follow Table 10, not Tables 5-8, so its b0 reads
 * 120-127 (a single tone can read as silence). A tone frame whose index is
 * invalid (Table 9), or whose redundant fields disagree, is an erasure (7.3).
 * Otherwise b0 gives the type (4.1 Table 4); b0 126/127 without the tone
 * identifier cannot be decoded as a tone and is an erasure too.
 */
static int
ambe2450_classify_frame(const char* ambe_d) {
    if (ambe2450_is_tone_frame(ambe_d)) {
        return (ambe2450_tone_fields_consistent(ambe_d) && ambe2450_tone_id_usable(ambe_d))
                   ? MBE_AMBE2450_FRAME_TONE
                   : MBE_AMBE2450_FRAME_ERASURE;
    }
    const int b0 = ambe2450_decode_b0(ambe_d);
    if (b0 < 120) {
        return MBE_AMBE2450_FRAME_VOICE;
    }
    if ((b0 == 124) || (b0 == 125)) {
        return MBE_AMBE2450_FRAME_SILENCE;
    }
    return MBE_AMBE2450_FRAME_ERASURE;
}

/* TIA-102.BABA-1 4.1 eqs 1-3: w0 = 2*pi/32, L = 14, all bands unvoiced, for
 * both silence b0 values (the encoder sends 124). The spectral amplitudes are
 * decoded like a voice frame's (4.4). JMBE's W124/W125 stored pi/32 where its
 * other entries hold f0 in cycles/sample, then scaled by 2*pi again, putting
 * harmonics 6..15 above Nyquist. */
static void
ambe2450_set_silence_model(mbe_parms* cur_mp, int* L, float* f0, int* silence) {
    *silence = 1;
    *f0 = 1.0f / 32.0f;
    cur_mp->w0 = MBE_AMBE_SILENCE_W0;
    *L = MBE_AMBE_SILENCE_L;
    cur_mp->L = *L;
    for (int l = 1; l <= *L; l++) {
        cur_mp->Vl[l] = 0;
    }
}

static int
ambe2450_setup_frame_model(const char* ambe_d, mbe_parms* cur_mp, int* b0, int* L, float* f0, int* silence) {
    const int kind = ambe2450_classify_frame(ambe_d);
    *silence = 0;
    *b0 = ambe2450_decode_b0(ambe_d);
#ifdef AMBE_DEBUG
    fprintf(stderr, "Frame type %d b0 = %d\n", kind, *b0);
#endif

    if (kind == MBE_AMBE2450_FRAME_SILENCE) {
        ambe2450_set_silence_model(cur_mp, L, f0, silence);
    } else if (kind == MBE_AMBE2450_FRAME_VOICE) {
        *f0 = AmbeW0table[*b0];
        cur_mp->w0 = *f0 * (float)2 * M_PI;
        *L = AmbeLtable[*b0];
        cur_mp->L = *L;
    }
    return kind;
}

/**
 * @brief Internal AMBE 2450 frame classification and parameter decode.
 *
 * Unlike the public wrapper, tells silence frames apart from voice frames.
 *
 * @param ambe_d  Demodulated AMBE parameter bits (49).
 * @param cur_mp  Output: current frame parameters (voice and silence frames).
 * @param prev_mp In/out: last valid voice frame (prediction history).
 * @return An MBE_AMBE2450_FRAME_* value, or a negative MBE_STATUS_* code.
 */
static int
mbe_decodeAmbe2450ParmsInternal(const char* ambe_d, mbe_parms* cur_mp, mbe_parms* prev_mp) {

    int b0, b2, L = 0;
    float f0, Cik[5][18], Tl[57] = {0}, Ri[9];
    int silence = 0;
    int Ji[5];
    float deltaGamma;
    float unvc;
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

    const int kind = ambe2450_setup_frame_model(ambe_d, cur_mp, &b0, &L, &f0, &silence);
    if ((kind != MBE_AMBE2450_FRAME_VOICE) && (kind != MBE_AMBE2450_FRAME_SILENCE)) {
        return kind;
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

    return kind;
}

/**
 * @brief Decode AMBE 2450 parameters from demodulated bitstream.
 * @param ambe_d  Demodulated AMBE parameter bits (49).
 * @param cur_mp  Output: current frame parameters.
 * @param prev_mp In/out: last valid voice frame (prediction history).
 * @return MBE_AMBE2450_FRAME_VOICE for a decoded model (voice or silence, as
 *         in 2.1), MBE_AMBE2450_FRAME_ERASURE, MBE_AMBE2450_FRAME_TONE, or a
 *         negative MBE_STATUS_* code.
 */
int
mbe_decodeAmbe2450Parms(const char* ambe_d, mbe_parms* cur_mp, mbe_parms* prev_mp) {
    const int kind = mbe_decodeAmbe2450ParmsInternal(ambe_d, cur_mp, prev_mp);
    return (kind == MBE_AMBE2450_FRAME_SILENCE) ? MBE_AMBE2450_FRAME_VOICE : kind;
}

/**
 * @brief Classify AMBE 2450 parameter bits without decoding them.
 * @param ambe_d Demodulated AMBE parameter bits (49).
 * @return An MBE_AMBE2450_FRAME_* value, or a negative MBE_STATUS_* code.
 */
int
mbe_classifyAmbe2450Frame(const char ambe_d[49]) {
    const int ret = mbe_validate_bits(ambe_d, 49u);
    if (ret < 0) {
        return ret;
    }
    return ambe2450_classify_frame(ambe_d);
}

/**
 * @brief Demodulate interleaved AMBE 3600x2450 data in-place.
 * @param ambe_fr Frame as 4x24 bitplanes (modified).
 */
int
mbe_demodulateAmbe3600x2450Data(char ambe_fr[4][24]) {
    int ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }
    mbe_demodulateAmbe3600Data_common(ambe_fr);
    return 0;
}

int
mbe_decodeAmbe3600x2450Frame(const char ambe_fr[4][24], char ambe_d[49], mbe_process_result* result) {
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
    c0_errors = mbe_eccAmbe3600x2450C0(fr);
    ret = mbe_demodulateAmbe3600x2450Data(fr);
    if (ret < 0) {
        return ret;
    }
    protected_errors = mbe_eccAmbe3600x2450Data(fr, ambe_d);

    if (result) {
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
ambe2450_prepare_process(mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                         mbe_parms* prev_mp_enhanced, int* total_errors, int* c0_errors, int* c0_errors_valid) {
    int ret;

    if (!cur_mp || !prev_mp || !prev_mp_enhanced || !total_errors || !c0_errors || !c0_errors_valid) {
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

    *c0_errors_valid = (result && ((result->flags & MBE_PROCESS_FLAG_C0_VALID) != 0u));
    *c0_errors = *c0_errors_valid ? result->c0_errors : 0;
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
ambe2450_repeat_required(const char ambe_d[49], int total_errors, int c0_errors, int c0_errors_valid) {
    if (c0_errors_valid) {
        /* TIA-102.BABA-1 5.6 repeat criteria when C0 errors are known. */
        return ((c0_errors >= 4) || ((c0_errors >= 2) && (total_errors >= 6)));
    }

    /* Dataf callers pass parameter bits only (no C0 context); keep the 2.1
     * heuristic, which also outputs a verified tone with fewer than 6 errors. */
    if (ambe2450_tone_verified(ambe_d) && (total_errors < 6)) {
        return 0;
    }
    return (total_errors > 3);
}

/*
 * AMBE 3600x2450 frame handling follows TIA-102.BABA-1 (Half-Rate Vocoder
 * Addendum). The caller-owned triplet carries the spec's two notions of
 * "previous frame":
 *  - prev_mp is the last valid voice frame: the prediction history (L, gamma,
 *    log2Ml; 4.3; 4.4.1 eq 26; 4.4.3 eqs 43-45) plus the per-frame decoder
 *    values below, which advance on every frame.
 *  - prev_mp_enhanced is the last synthesized frame and the frame-repeat
 *    source (5.6 eqs 59-64: w0, L, K, v, enhanced M).
 * Invalid frames never commit history. At frame start, cur_mp's noise state
 * equals prev_mp_enhanced's; every path below preserves that.
 */

/**
 * Not in TIA-102.BABA-1: silence the last synthesized frame (zero its
 * amplitudes and unvoiced overlap; keep its phase and noise state) after a
 * frame whose output was not synthesized speech: a mute or a tone. Whatever
 * is synthesized next then fades in from silence instead of overlapping
 * speech from before that frame, and a repeat right after it replays silence.
 * The spec does not define the synthesis state across tone frames; JMBE
 * resets to init defaults after a mute.
 */
static void
ambe2450_silence_last_synthesized(mbe_parms* prev_mp_enhanced) {
    memset(prev_mp_enhanced->Ml, 0, sizeof(prev_mp_enhanced->Ml));
    memset(prev_mp_enhanced->previousUw, 0, sizeof(prev_mp_enhanced->previousUw));
}

/**
 * 5.7: compute the 5.6 step (2) update (the repeat model, which keeps this
 * frame's error, tone and noise state), bypass synthesis and output uniform
 * noise in [-5, 5]; then silence the last synthesized frame.
 */
static void
ambe2450_mute(float* aout_buf, mbe_process_result* result, mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced) {
    mbe_repeat_load_model(cur_mp, prev_mp_enhanced);
    mbe_result_set_flag(result, MBE_PROCESS_FLAG_MUTE);
    mbe_synthesizeUniformNoisef(aout_buf, MBE_SPEC_MUTE_NOISE_AMPLITUDE);
    ambe2450_silence_last_synthesized(prev_mp_enhanced);
}

/** Enhance and synthesize a decoded frame, then make it the last synthesized frame. */
static void
ambe2450_synthesize_decoded(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced) {
    float pre_enh_rm0 = mbe_spectralAmpEnhanceWithRm0(cur_mp);
    mbe_synthesizeSpeechWithPreEnhRm0f(aout_buf, cur_mp, prev_mp_enhanced, pre_enh_rm0, MBE_MUTE_NOISE_SPEC);
    mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
}

/**
 * 5.6 steps 2-3: synthesize the repeated model as-is (no re-enhancement and
 * no adaptive smoothing, which could change its voicing and amplitudes), or
 * mute (5.7) instead of the 4th consecutive repeat.
 */
static void
ambe2450_repeat(float* aout_buf, mbe_process_result* result, mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced) {
    if (mbe_isMaxFrameRepeat(cur_mp) || mbe_requiresMuting(cur_mp)) {
        ambe2450_mute(aout_buf, result, cur_mp, prev_mp_enhanced);
        return;
    }
    /* Step 2 repeats the previous frame's model, enhanced amplitudes included
     * (eqs 59-64), so synthesis continues its WOLA and phase state. */
    mbe_repeat_load_model(cur_mp, prev_mp_enhanced);
    mbe_synthesizeRepeatedSpeechf(aout_buf, cur_mp, prev_mp_enhanced);
    mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
}

/**
 * Invalid frame: repeat criteria or erasure (5.6), or a tone frame with an
 * unusable tone index (7.3). The frame is ignored in later processing (5.6
 * step 1); the 4th consecutive invalid frame mutes instead of repeating (5.7).
 */
static void
ambe2450_process_invalid(float* aout_buf, mbe_process_result* result, unsigned extra_flags, mbe_parms* cur_mp,
                         mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_result_set_flag(result, extra_flags | MBE_PROCESS_FLAG_REPEAT);
    cur_mp->repeatCount = mbe_repeat_next_count(prev_mp->repeatCount);
    mbe_copy_error_state(prev_mp, cur_mp);
    ambe2450_repeat(aout_buf, result, cur_mp, prev_mp_enhanced);
}

/**
 * Valid voice or silence frame. A voice frame becomes the prediction history;
 * the encoder updated its predictor on it, so it commits even when the error
 * rate then mutes the output (5.7 does not say either way). A silence frame is
 * synthesized but never used for prediction (4.3; 4.4.1 eq 26; 4.4.3 eq 43):
 * only the per-frame decoder values advance.
 */
static void
ambe2450_process_valid(float* aout_buf, mbe_process_result* result, int silence, mbe_parms* cur_mp, mbe_parms* prev_mp,
                       mbe_parms* prev_mp_enhanced) {
    cur_mp->repeatCount = 0;
    if (silence) {
        mbe_result_set_flag(result, MBE_PROCESS_FLAG_SILENCE);
        mbe_copy_error_state(prev_mp, cur_mp);
    } else {
        mbe_moveMbeParms(cur_mp, prev_mp);
    }
    if (mbe_requiresMuting(cur_mp)) {
        ambe2450_mute(aout_buf, result, cur_mp, prev_mp_enhanced);
        return;
    }
    ambe2450_synthesize_decoded(aout_buf, cur_mp, prev_mp_enhanced);
}

/**
 * Tone frame with a usable index (7.3): output the tone, or silence for the
 * zero-amplitude ID 255, instead of synthesized speech. The prediction history
 * is not touched (4.3); the last synthesized frame is silenced (not in the
 * spec, see ambe2450_silence_last_synthesized).
 */
static void
ambe2450_process_tone(float* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                      mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_result_set_flag(result, MBE_PROCESS_FLAG_TONE);
    cur_mp->repeatCount = 0;
    mbe_copy_error_state(prev_mp, cur_mp);
    if (mbe_requiresMuting(cur_mp)) {
        ambe2450_mute(aout_buf, result, cur_mp, prev_mp_enhanced);
        return;
    }
    if (ambe2450_read_tone_id(ambe_d) == AMBE2450_TONE_ID_ZERO_AMPLITUDE) {
        mbe_synthesizeSilencef(aout_buf);
    } else {
        mbe_synthesizeTonef(aout_buf, ambe_d, cur_mp);
    }
    ambe2450_silence_last_synthesized(prev_mp_enhanced);
}

static void
ambe2450_process_frame(float* aout_buf, mbe_process_result* result, const char ambe_d[49], int kind, mbe_parms* cur_mp,
                       mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    if (kind == MBE_AMBE2450_FRAME_TONE) {
        ambe2450_process_tone(aout_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
    } else if (kind == MBE_AMBE2450_FRAME_ERASURE) {
        /* A tone frame with an unusable index is flagged as both. */
        const unsigned tone = ambe2450_is_tone_frame(ambe_d) ? MBE_PROCESS_FLAG_TONE : 0u;
        ambe2450_process_invalid(aout_buf, result, MBE_PROCESS_FLAG_ERASURE | tone, cur_mp, prev_mp, prev_mp_enhanced);
    } else {
        ambe2450_process_valid(aout_buf, result, kind == MBE_AMBE2450_FRAME_SILENCE, cur_mp, prev_mp, prev_mp_enhanced);
    }
}

static int
mbe_processAmbe2450Dataf_internal(float* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                                  mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    int kind;
    int c0_errors;
    int c0_errors_valid;
    int total_errors;
    int ret;

    /* Aliased sets would let a mute's silencing of prev_mp_enhanced erase the history. */
    if (!aout_buf || !mbe_parms_triplet_is_valid(cur_mp, prev_mp, prev_mp_enhanced)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = ambe2450_prepare_process(result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced, &total_errors, &c0_errors,
                                   &c0_errors_valid);
    if (ret < 0) {
        return ret;
    }

    /* 5.6: the repeat criteria apply to every frame type, before classification. */
    if (ambe2450_repeat_required(ambe_d, total_errors, c0_errors, c0_errors_valid)) {
        ambe2450_process_invalid(aout_buf, result, 0u, cur_mp, prev_mp, prev_mp_enhanced);
        return result ? result->total_errors : total_errors;
    }

    kind = mbe_decodeAmbe2450ParmsInternal(ambe_d, cur_mp, prev_mp);
    if (kind < 0) {
        return kind;
    }
    ambe2450_process_frame(aout_buf, result, ambe_d, kind, cur_mp, prev_mp, prev_mp_enhanced);
    return result ? result->total_errors : total_errors;
}

/**
 * @brief Process AMBE 2450 parameters into 160 float samples at 8 kHz.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional in/out status context.
 * @param ambe_d   Demodulated parameter bits (49).
 * @param cur_mp   In/out: current frame parameters (may be enhanced).
 * @param prev_mp  In/out: previous frame parameters.
 * @param prev_mp_enhanced In/out: enhanced previous parameters for continuity.
 */
int
mbe_processAmbe2450Dataf(float* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                         mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;

    if (!result) {
        mbe_initProcessResult(&local_result);
        result = &local_result;
    }
    return mbe_processAmbe2450Dataf_internal(aout_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processAmbe2450Data(short* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                        mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processAmbe2450Dataf(float_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete AMBE 3600x2450 frame into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional output status populated by decode and synthesis.
 * @param ambe_fr  Input frame as 4x24 bitplanes.
 * @param ambe_d   Scratch/output parameter bits (49).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 */
int
mbe_processAmbe3600x2450Framef(float* aout_buf, mbe_process_result* result, const char ambe_fr[4][24], char ambe_d[49],
                               mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;
    int ret;
    if (!result) {
        result = &local_result;
    }
    ret = mbe_decodeAmbe3600x2450Frame(ambe_fr, ambe_d, result);
    if (ret < 0) {
        return ret;
    }
    return mbe_processAmbe2450Dataf(aout_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processAmbe3600x2450SoftFramef(float* aout_buf, mbe_process_result* result, const mbe_soft_bit ambe_fr[4][24],
                                   char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                   mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;
    int ret;
    if (!result) {
        result = &local_result;
    }
    ret = mbe_decodeAmbe3600x2450SoftFrame(ambe_fr, ambe_d, result);
    if (ret < 0) {
        return ret;
    }
    return mbe_processAmbe2450Dataf(aout_buf, result, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processAmbe3600x2450SoftFrame(short* aout_buf, mbe_process_result* result, const mbe_soft_bit ambe_fr[4][24],
                                  char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processAmbe3600x2450SoftFramef(float_buf, result, ambe_fr, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete AMBE 3600x2450 frame into 16-bit PCM.
 * @see mbe_processAmbe3600x2450Framef for details.
 */
int
mbe_processAmbe3600x2450Frame(short* aout_buf, mbe_process_result* result, const char ambe_fr[4][24], char ambe_d[49],
                              mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;

    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processAmbe3600x2450Framef(float_buf, result, ambe_fr, ambe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}
