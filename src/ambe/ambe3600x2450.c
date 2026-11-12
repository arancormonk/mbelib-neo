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
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ambe3600x2450_const.h"
#include "ambe3600x2450_internal.h"
#include "ambe_common.h"
#include "mbe_adaptive.h"
#include "mbe_compiler.h"
#include "mbe_math.h"
#include "mbe_repeat.h"
#include "mbe_result.h"
#include "mbe_tone.h"
#include "mbe_unvoiced_fft.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"

/*
 * Bit positions of b0..b8 in ambe_d, most significant first, -1 ending each
 * field: TIA-102.BABA-1 5.1 Tables 5-8 as DVSI's AMBE-3000 applies them.
 * u3's bits 8..5 (ambe_d 40..43) carry b4(2), b4(1), b4(0) and then b3(0);
 * the v1.0.5 text of Table 8 prints b3(0) first. DVSI's rate-34 frames carry
 * the same 49-bit vector as this decoder reads, and with this assignment
 * DVSI's own re-quantizations of a decoded model pick the PRBA vectors this
 * library's quantizer picks (84-91% exact against 8-34% with the printed one).
 */
static const signed char ambe2450_layout[9][10] = {
    {0, 1, 2, 3, 37, 38, 39, -1},
    {4, 5, 6, 7, 35, -1},
    {8, 9, 10, 11, 36, -1},
    {12, 13, 14, 15, 16, 17, 18, 19, 43, -1},
    {20, 21, 22, 23, 40, 41, 42, -1},
    {24, 25, 26, 27, 44, -1},
    {28, 29, 30, 45, -1},
    {31, 32, 33, 46, -1},
    {34, 47, 48, -1},
};

/* Quantizer value b<field> from ambe_d. */
static int
ambe2450_field(const char* ambe_d, int field) {
    int value = 0;
    for (int j = 0; ambe2450_layout[field][j] >= 0; j++) {
        value = (value << 1) | (ambe_d[(int)ambe2450_layout[field][j]] & 1);
    }
    return value;
}

void
mbe_ambe2450_unpack(const char ambe_d[49], int b[9]) {
    for (int i = 0; i < 9; i++) {
        b[i] = ambe2450_field(ambe_d, i);
    }
}

void
mbe_ambe2450_pack(const int b[9], char ambe_d[49]) {
    for (int i = 0; i < 9; i++) {
        int width = 0;
        while (ambe2450_layout[i][width] >= 0) {
            width++;
        }
        for (int j = 0; j < width; j++) {
            ambe_d[(int)ambe2450_layout[i][j]] = (char)((b[i] >> (width - 1 - j)) & 1);
        }
    }
}

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

/* The tone amplitude field AD: the low six bits of U0, then bit 4 of U3. */
static int
ambe2450_read_tone_ad(const char ambe_d[49]) {
    int ad = 0;
    for (int i = 6; i < 12; i++) {
        ad = (ad << 1) | (int)ambe_d[i];
    }
    return (ad << 1) | (int)ambe_d[44];
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
 * @return Bit errors corrected or detected in C0, parity bits included
 *         (TIA-102.BABA-1 5.5): 0..3, or 4 when the extended parity reveals
 *         an uncorrectable pattern.
 */
int
mbe_eccAmbe3600x2450C0(char ambe_fr[4][24]) {
    int ret = mbe_validate_bits((const char*)ambe_fr, (size_t)4u * 24u);
    if (ret < 0) {
        return ret;
    }
    return mbe_eccAmbe3600C0Counted_common(ambe_fr);
}

/**
 * @brief Apply ECC to AMBE 3600x2450 data and pack parameter bits.
 * @param ambe_fr Frame as 4x24 bitplanes.
 * @param ambe_d  Output parameter bits (49).
 * @return Bit errors corrected in C1, parity bits included (TIA-102.BABA-1 5.5).
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
    return mbe_eccAmbe3600DataCounted_common(ambe_fr, ambe_d);
}

int
mbe_ambe2450_harmonic_voiced(const int row[8], int l, float f0) {
    const float column = mbe_round_to_float((float)l * 16.0f * f0);
    int own = (int)column;
    int above = (int)mbe_round_to_float(column + 0.5f);
    own = (own > 7) ? 7 : own;
    above = (above > 7) ? 7 : above;
    return row[own] | row[above];
}

/* The frame's voicing, and its mark: a silence frame is marked
 * MBE_FRAME_MARK_PREVIOUS_HELD (ambe2450_set_silence_model()), a voice frame
 * of row MBE_AMBE2450_VUV_LATE_ONSET a late onset, any other none. */
static void
ambe2450_decode_vuv(const char* ambe_d, mbe_parms* cur_mp, int L, float f0, int silence, int b0) {
    (void)b0;
    int b1 = ambe2450_field(ambe_d, 1);

    if (silence) {
        mbe_set_frame_mark(cur_mp, MBE_FRAME_MARK_PREVIOUS_HELD);
    } else {
        mbe_set_voicing_onset_delay(cur_mp, b1 == MBE_AMBE2450_VUV_LATE_ONSET);
    }
    for (int l = 1; l <= L; l++) {
        if (silence == 0) {
            cur_mp->Vl[l] = mbe_ambe2450_harmonic_voiced(AmbeVuv[b1], l, f0);
        }
#ifdef AMBE_DEBUG
        fprintf(stderr, "Vl[%i]:%i\n", l, cur_mp->Vl[l]);
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

    int b3 = ambe2450_field(ambe_d, 3);
    Gm[2] = AmbePRBA24[b3][0];
    Gm[3] = AmbePRBA24[b3][1];
    Gm[4] = AmbePRBA24[b3][2];

    int b4 = ambe2450_field(ambe_d, 4);
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

    int b5 = ambe2450_field(ambe_d, 5);
    int b6 = ambe2450_field(ambe_d, 6);
    int b7 = ambe2450_field(ambe_d, 7);
    int b8 = ambe2450_field(ambe_d, 8);

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
        flokl[l] = mbe_prediction_position(prev_L, cur_mp->L, l);
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
    return ambe2450_field(ambe_d, 0);
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
 * Independent local consistency check over the redundant tone fields in
 * TIA-102.BABA-1 7.2 Table 10 (the spec leaves checking to the decoder):
 * the last four bits of u3 are 0, or u1's trailing copy of ID(7..4) matches.
 */
static int
ambe2450_tone_fields_consistent(const char* bits) {
    if ((bits[45] | bits[46] | bits[47] | bits[48]) == 0) {
        return 1;
    }
    for (int i = 0; i < 4; ++i) {
        if (bits[12 + i] != bits[20 + i]) {
            return 0;
        }
    }
    return 1;
}

/** A verified tone frame: the tone identifier (7) plus consistent redundant fields. */
static int
ambe2450_tone_verified(const char* ambe_d) {
    return ambe2450_is_tone_frame(ambe_d) && ambe2450_tone_fields_consistent(ambe_d);
}

/**
 * 7.3: a tone index the decoder can output. Table 9 also lists ID 255 as a
 * tone of zero amplitude; DVSI's decoder does not play it as a tone but as a
 * lost frame, as it does any unusable index (measurements of an AMBE-3000R
 * published with mbevocoder: ID 255 and ID 200 play as lost frames). This
 * decoder follows DVSI (docs/dvsi-parity.md, "Rate 33 lost frames and
 * tones"): ID 255 is an erasure.
 */
static int
ambe2450_tone_id_usable(const char* ambe_d) {
    return mbe_tone_id_is_valid(ambe2450_read_tone_id(ambe_d));
}

/**
 * Frame type from the parameter bits alone. A tone frame (7) is recognised
 * before b0: its other bits follow Table 10, not Tables 5-8, so its b0 reads
 * 120-127 (a single tone can read as silence). A tone frame whose index is
 * invalid (Table 9) or the zero-amplitude ID 255 (ambe2450_tone_id_usable()),
 * or whose redundant fields disagree, is an erasure (7.3).
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

/*
 * Silence frames (b0 124 and 125; DVSI's encoder sends 124, and 125 is taken
 * to play as it does). TIA-102.BABA-1 4.1 eqs 1-3 give w0 = 2*pi/32, L = 14
 * and every band unvoiced, the spectral amplitudes decoded like a voice
 * frame's (4.4); the spec says nothing of how the frame joins the one before.
 * DVSI's decoder plays the frame on 15 bands at that pitch, its amplitudes a
 * further factor 0.2046 (13.8 dB) down, and plays the previous frame's model
 * (pitch, harmonic count, voicing, amplitudes) through the first 10 ms of the
 * silence frame, crossfading to the silence model over the second; the pulse
 * train runs on through one silence frame into a voiced frame after it. This
 * decoder follows DVSI (docs/dvsi-parity.md, "Rate 33 silence frames"): the
 * residual is reconstructed on 15 bands (block lengths AmbeLmprbl[15]) and
 * predicted from the history at 15 bands' taps, by the voice frame's law; the
 * extra factor scales Ml alone (log2Ml carries none of it); and the frame is
 * marked MBE_FRAME_MARK_PREVIOUS_HELD, whose rules in the synthesizer hold
 * the previous model to mid-frame and run the train on. A silence frame still
 * never becomes the prediction history (4.3). D-STAR's b0 124 keeps
 * MBE_AMBE_SILENCE_L.
 */
#define MBE_AMBE2450_SILENCE_L    15
#define MBE_AMBE2450_SILENCE_GAIN 0.2046f

static void
ambe2450_set_silence_model(mbe_parms* cur_mp, int* L, float* f0, int* silence) {
    *silence = 1;
    *f0 = 1.0f / 32.0f;
    cur_mp->w0 = MBE_AMBE_SILENCE_W0;
    *L = MBE_AMBE2450_SILENCE_L;
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
        *f0 = mbe_ambe2450_f0(*b0);
        cur_mp->w0 = *f0 * (float)2 * M_PI;
        *L = mbe_ambe2450_harmonic_count(*b0);
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
    if (silence) {
        unvc *= MBE_AMBE2450_SILENCE_GAIN;
    }

    // decode V/UV parameters
    ambe2450_decode_vuv(ambe_d, cur_mp, L, f0, silence, b0);

    // decode gain vector
    // load b2 from ambe_d
    b2 = ambe2450_field(ambe_d, 2);

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
    c0_errors = mbe_eccAmbe3600C0SoftCounted_common(fr);
    mbe_demodulateAmbe3600DataSoft_common(fr);
    protected_errors = mbe_eccAmbe3600DataSoftCounted_common(fr, ambe_d);

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
        /* With the 5.5 counts a frame repeats only when C0 shows four errors
         * (detected, uncorrectable). DVSI's AMBE-3000 does not apply eq
         * 57-58's second condition (eps0 >= 2 and epsT >= 6): under channel
         * errors its output repeats where eps0 is 4 and not where only that
         * condition holds, and following it puts the frames it plays closer
         * to both DVSI's output and the input speech. */
        return c0_errors >= 4;
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
 *
 * Not in TIA-102.BABA-1, which does not say what a mute or a tone leaves the
 * frames after it: DVSI's decoder silences its output there, not its model
 * (docs/dvsi-parity.md, "Rate 33 mutes"). prev_mp_enhanced keeps the model
 * the decoder would have played (the "kept model"): through a mute each valid
 * voice or silence frame's own model, enhanced and smoothed as if it were
 * synthesized, and an invalid frame's repeated one; through a tone the last
 * speech model before it. Only its unvoiced overlap is cleared. The otherwise
 * unused Vl[0] of prev_mp_enhanced (Vl is indexed from 1) holds a word that
 * says what the last frame left:
 *  - AMBE2450_TONE_HELD with the tone's ID (bits 0-7) and AD (bits 8-14): a
 *    tone is sounding;
 *  - AMBE2450_STARTS_AFRESH: the kept model was not heard at the end of the
 *    last frame (a mute or a tone came), so the next frame synthesized starts
 *    from silence. A held tone always comes with it.
 * The word is read and cleared at a frame's start (ambe2450_take_snapshot());
 * each branch of ambe2450_process_frame() produces the frame's output and the
 * word it leaves, which is written after the frame's last move of a parameter
 * set, so that it never travels with a model; cur_mp->Vl[0] is 0 on the way
 * out.
 */
#define AMBE2450_TONE_FADE     80
#define AMBE2450_TONE_HELD     (1 << 16)
#define AMBE2450_STARTS_AFRESH (1 << 17)

/* What the last frame left, read at a frame's start. */
struct ambe2450_snapshot {
    int held;            /* a tone is sounding: id and ad are its fields */
    int afresh;          /* the next frame synthesized starts from silence */
    int id;              /* tone ID (Table 9) */
    int ad;              /* tone amplitude field */
    int swn;             /* tone oscillator phases at the frame's start */
    uint32_t tone_phase; /* (the state cur_mp carries) */
};

static void
ambe2450_take_snapshot(const mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced, struct ambe2450_snapshot* in) {
    const int stored = prev_mp_enhanced->Vl[0];
    in->held = (stored & AMBE2450_TONE_HELD) != 0;
    in->afresh = in->held || ((stored & AMBE2450_STARTS_AFRESH) != 0);
    in->id = stored & 0xFF;
    in->ad = (stored >> 8) & 0x7F;
    in->swn = cur_mp->swn;
    in->tone_phase = cur_mp->tonePhase;
    prev_mp_enhanced->Vl[0] = 0;
}

/* Whether the stream has synthesized or muted a frame: DVSI's noise
 * generator then holds a 16-bit state, where a fresh stream holds a negative
 * one (read as bits, see dvsi_noise_state()). */
static int
ambe2450_has_played(const mbe_parms* prev_mp_enhanced) {
    uint32_t bits;
    memcpy(&bits, &prev_mp_enhanced->noiseSeed, sizeof(bits));
    return bits <= UINT32_C(0x477FFF00); /* 65535.0f */
}

/* The kept model's unvoiced overlap: nothing it could overlap was played. */
static void
ambe2450_clear_overlap(mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced) {
    memset(cur_mp->previousUw, 0, sizeof(cur_mp->previousUw));
    memset(prev_mp_enhanced->previousUw, 0, sizeof(prev_mp_enhanced->previousUw));
}

/*
 * The tone the last frame left sounding, going on from its phases, added to
 * the frame's first AMBE2450_TONE_FADE samples as it fades out: (1 - n / 80)
 * tail[n]. The sum is limited.
 */
static void
ambe2450_add_tone_tail(float* aout_buf, const struct ambe2450_snapshot* in) {
    float tail[160];
    mbe_parms tone;
    memset(&tone, 0, sizeof(tone));
    tone.swn = in->swn;
    tone.tonePhase = in->tone_phase;
    mbe_renderToneAmbe2450f(tail, &tone, in->id, in->ad);
    for (int n = 0; n < AMBE2450_TONE_FADE; ++n) {
        aout_buf[n] += (1.0f - ((float)n / (float)AMBE2450_TONE_FADE)) * tail[n];
    }
    mbe_clipFloatBuffer(aout_buf, AMBE2450_TONE_FADE);
}

/*
 * The frame after a mute or a tone starts from silence. TIA-102.BABA-1 does
 * not say how. DVSI's decoder synthesizes the first frame heard again from
 * the model it kept in the ordinary way, with two differences: the kept model
 * itself is not played (the first subframe rises from silence to the
 * mid-frame model, which is the ordinary one between the kept model and the
 * frame's own; no unvoiced overlap comes from before), and the pulse train
 * starts anew at the frame's start (MBE_SYNTHESIS_FROM_SILENCE). A silence
 * frame, whose hold continues what was heard, does not hold the kept model,
 * which was not: it rises from silence to its own model. On DVSI's
 * dam_e10_hd the frames leaving a mute follow DVSI's decode at a median
 * correlation of 0.95 / 0.95 / 0.96 (first half, second half, next frame)
 * this way, against 0.00 / 0.04 / 0.09 for a fade-in from nothing, and over
 * the first frames after the tones of DVSI's 52 rate-33 tone vectors the mean
 * absolute level error falls from 0.58 to 0.37 dB (docs/dvsi-parity.md, "Rate
 * 33 frames after a mute or a tone").
 */
static enum mbe_synthesis_start
ambe2450_start(int afresh) {
    return afresh ? MBE_SYNTHESIS_FROM_SILENCE : MBE_SYNTHESIS_CONTINUES;
}

/* Enhance and synthesize a decoded frame, then make it the last synthesized frame. */
static void
ambe2450_synthesize_decoded(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced, int afresh) {
    float pre_enh_rm0 = mbe_spectralAmpEnhanceWithRm0(cur_mp, MBE_HIGH_BAND_GAIN_AMBE2450);
    mbe_synthesizeSpeechWithPreEnhRm0Startf(aout_buf, cur_mp, prev_mp_enhanced, pre_enh_rm0, MBE_MUTE_NOISE_SPEC,
                                            MBE_VOICED_SYNTHESIS_DVSI, MBE_AMPLITUDE_LIMIT_NONE,
                                            ambe2450_start(afresh));
    mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
}

/*
 * The first muted frame after a frame that was heard: DVSI's decoder plays the
 * last frame once more, fading out linearly over the first 80 samples
 * (dam_e10_hd), as a repeat of it toward a silent copy of itself fades out over
 * the first subframe (as ambe2400_mute_repeat_limit() renders D-STAR's). Only
 * scratch copies are synthesized; `faded` receives the state the fade leaves.
 * Its repeat count is 0: the synthesis core mutes a model at the repeat limit.
 * The silent copy carries no frame mark: it is no decoded frame, and the last
 * frame's mark would act on the fade's mid-frame model (a silence frame's
 * would hold the last model there at full level, moving the fade into the
 * second subframe), so the mid-frame model is the ordinary one between the
 * last model and silence, whatever kind of frame the last one was. Its PHIl[0]
 * reaches nothing kept (ambe2450_keep_last() copies phases from element 1).
 */
static void
ambe2450_fade_out_last(float* aout_buf, const mbe_parms* cur_mp, const mbe_parms* prev_mp_enhanced, mbe_parms* faded) {
    mbe_parms last = *prev_mp_enhanced;
    *faded = *cur_mp;
    mbe_repeat_load_model(faded, prev_mp_enhanced);
    memset(faded->Ml, 0, sizeof(faded->Ml));
    mbe_set_frame_mark(faded, MBE_FRAME_MARK_NONE);
    faded->repeatCount = 0;
    mbe_synthesizeRepeatedSpeechf(aout_buf, faded, &last, MBE_VOICED_SYNTHESIS_DVSI);
}

/*
 * A muted frame that is not a valid voice or silence frame keeps the last
 * synthesized model (with this frame's error accounting), its pulse state and
 * harmonic phases moved on by the fade when there was one; else the noise
 * generator runs one frame on its own.
 */
static void
ambe2450_keep_last(mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced, const mbe_parms* faded) {
    mbe_repeat_load_model(cur_mp, prev_mp_enhanced);
    if (faded != NULL) {
        cur_mp->noiseSeed = faded->noiseSeed;
        memcpy(cur_mp->noiseOverlap, faded->noiseOverlap, sizeof(cur_mp->noiseOverlap));
        memcpy(&cur_mp->PSIl[1], &faded->PSIl[1], sizeof(cur_mp->PSIl) - sizeof(cur_mp->PSIl[0]));
        memcpy(&cur_mp->PHIl[1], &faded->PHIl[1], sizeof(cur_mp->PHIl) - sizeof(cur_mp->PHIl[0]));
    } else {
        mbe_dvsi_noise_skip_frame(cur_mp);
    }
    mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
}

/*
 * 5.7: mute (the error rate above 0.096, or the 4th consecutive invalid
 * frame). TIA-102.BABA-1 bypasses synthesis for uniform noise in [-5, 5].
 * DVSI's decoder plays its noise at an rms of about 0.6 (MBE_AMBE_MUTE_NOISE_AMPLITUDE,
 * as D-STAR's), and over the first muted frame after a heard one the last
 * frame fading out (ambe2450_fade_out_last()), or the tone left sounding
 * fading out under unfaded noise; later muted frames, and a mute before
 * anything was played, are the noise alone (docs/dvsi-parity.md, "Rate 33
 * mutes"). Its model goes on: a valid voice or silence frame (`own_model`)
 * becomes the kept model as an unheard synthesis of it leaves it, enhanced and
 * smoothed; any other keeps the last one. DVSI's noise generator runs one
 * frame. The sum is limited. Returns the word the frame leaves.
 */
static int
ambe2450_mute(float* aout_buf, mbe_process_result* result, const struct ambe2450_snapshot* in, int own_model,
              mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced) {
    const int fading = !in->held && !in->afresh && ambe2450_has_played(prev_mp_enhanced);
    float noise[160];
    mbe_parms faded;

    mbe_result_set_flag(result, MBE_PROCESS_FLAG_MUTE);
    if (fading) {
        ambe2450_fade_out_last(aout_buf, cur_mp, prev_mp_enhanced, &faded);
    } else {
        mbe_synthesizeSilencef(aout_buf);
        if (in->held) {
            ambe2450_add_tone_tail(aout_buf, in);
        }
    }
    if (own_model) {
        float unheard[160];
        ambe2450_synthesize_decoded(unheard, cur_mp, prev_mp_enhanced, in->afresh);
    } else {
        ambe2450_keep_last(cur_mp, prev_mp_enhanced, fading ? &faded : NULL);
    }
    ambe2450_clear_overlap(cur_mp, prev_mp_enhanced);
    mbe_synthesizeUniformNoisef(noise, MBE_AMBE_MUTE_NOISE_AMPLITUDE * MBE_SPEECH_OUTPUT_GAIN);
    for (int n = 0; n < 160; ++n) {
        aout_buf[n] += noise[n];
    }
    mbe_clipFloatBuffer(aout_buf, 160);
    return AMBE2450_STARTS_AFRESH;
}

/**
 * 5.6 steps 2-3: synthesize the repeated model as-is (no re-enhancement and
 * no adaptive smoothing, which could change its voicing and amplitudes).
 * After a mute or a tone the kept model is repeated from silence. Returns the
 * word the frame leaves.
 */
static int
ambe2450_repeat(float* aout_buf, const struct ambe2450_snapshot* in, mbe_parms* cur_mp, mbe_parms* prev_mp_enhanced) {
    /* Step 2 repeats the previous frame's model, enhanced amplitudes included
     * (eqs 59-64), so synthesis continues its WOLA and phase state. */
    mbe_repeat_load_model(cur_mp, prev_mp_enhanced);
    mbe_synthesizeRepeatedSpeechStartf(aout_buf, cur_mp, prev_mp_enhanced, MBE_VOICED_SYNTHESIS_DVSI,
                                       ambe2450_start(in->afresh));
    mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
    return 0;
}

/*
 * A lost frame while a tone sounds. TIA-102.BABA-1 5.6 repeats the last
 * synthesized frame's model, and says nothing of a tone. DVSI's decoder
 * repeats the frame last played, a tone included (measurements of an
 * AMBE-3000R published with mbevocoder: a lost frame after a tone replays the
 * tone; docs/dvsi-parity.md, "Rate 33 lost frames and tones"). Like any
 * repeat, the tone goes on at its level for up to three lost frames, the
 * repeat count running, and the fourth, the repeat limit's mute, fades it out
 * under the mute noise (ambe2450_mute()); lost frames inside a tone are not
 * heard. The held tone plays again from the oscillator phases cur_mp carries,
 * as a tone frame advances them; the kept model and the word stay as they
 * are, and the noise generator stands still, as through a tone frame.
 */
static int
ambe2450_repeat_tone(float* aout_buf, const struct ambe2450_snapshot* in, mbe_parms* cur_mp) {
    mbe_renderToneAmbe2450f(aout_buf, cur_mp, in->id, in->ad);
    return AMBE2450_TONE_HELD | AMBE2450_STARTS_AFRESH | (in->ad << 8) | in->id;
}

/**
 * Invalid frame: repeat criteria or erasure (5.6), or a tone frame with an
 * unusable tone index (7.3). The frame is ignored in later processing (5.6
 * step 1); the 4th consecutive invalid frame mutes instead of repeating (5.7).
 * While a tone sounds it is the tone that repeats. Returns the word the frame
 * leaves.
 */
static int
ambe2450_process_invalid(float* aout_buf, mbe_process_result* result, unsigned extra_flags,
                         const struct ambe2450_snapshot* in, mbe_parms* cur_mp, mbe_parms* prev_mp,
                         mbe_parms* prev_mp_enhanced) {
    mbe_result_set_flag(result, extra_flags | MBE_PROCESS_FLAG_REPEAT);
    cur_mp->repeatCount = mbe_repeat_next_count(prev_mp->repeatCount);
    mbe_copy_error_state(prev_mp, cur_mp);
    if (mbe_isMaxFrameRepeat(cur_mp) || mbe_requiresMuting(cur_mp)) {
        return ambe2450_mute(aout_buf, result, in, 0, cur_mp, prev_mp_enhanced);
    }
    if (in->held) {
        return ambe2450_repeat_tone(aout_buf, in, cur_mp);
    }
    return ambe2450_repeat(aout_buf, in, cur_mp, prev_mp_enhanced);
}

/**
 * Valid voice or silence frame. A voice frame becomes the prediction history;
 * the encoder updated its predictor on it, so it commits even when the error
 * rate then mutes the output (5.7 does not say either way). A silence frame is
 * synthesized but never used for prediction (4.3; 4.4.1 eq 26; 4.4.3 eq 43):
 * only the per-frame decoder values advance. After a tone the tone's tail
 * fades out over the frame's first 80 samples. Returns the word the frame
 * leaves.
 */
static int
ambe2450_process_valid(float* aout_buf, mbe_process_result* result, int silence, const struct ambe2450_snapshot* in,
                       mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    cur_mp->repeatCount = 0;
    if (silence) {
        mbe_result_set_flag(result, MBE_PROCESS_FLAG_SILENCE);
        mbe_copy_error_state(prev_mp, cur_mp);
    } else {
        mbe_moveMbeParms(cur_mp, prev_mp);
    }
    if (mbe_requiresMuting(cur_mp)) {
        return ambe2450_mute(aout_buf, result, in, 1, cur_mp, prev_mp_enhanced);
    }
    ambe2450_synthesize_decoded(aout_buf, cur_mp, prev_mp_enhanced, in->afresh);
    if (in->held) {
        ambe2450_add_tone_tail(aout_buf, in);
    }
    return 0;
}

/*
 * Tone transitions, as DVSI's decoder plays them (docs/dvsi-parity.md, "Rate
 * 33 tone transitions"): a tone fades in linearly over the first
 * AMBE2450_TONE_FADE samples of its first frame while the last frame's speech,
 * synthesized as a repeat of it, fades out (nothing does after a mute); a
 * tone that changes while it sounds, or comes back after lost frames, changes
 * at once; and the last tone goes on through the first AMBE2450_TONE_FADE
 * samples of the frame after it, fading out on top of that frame's output
 * (or of a mute's noise). TIA-102.BABA-1 does not say how a tone frame joins
 * the frames around it. Its phases are the tone state cur_mp carries.
 */

/* Gain of the incoming signal at sample n of a transition frame. */
static float
ambe2450_fade_in_gain(int n) {
    return (n >= AMBE2450_TONE_FADE) ? 1.0f : (float)n / (float)AMBE2450_TONE_FADE;
}

/* Crossfade from `outgoing` to the frame in aout_buf over the fade span. */
static void
ambe2450_fade_in_over(float* aout_buf, const float* outgoing) {
    for (int n = 0; n < AMBE2450_TONE_FADE; ++n) {
        const float gain = ambe2450_fade_in_gain(n);
        aout_buf[n] = (gain * aout_buf[n]) + ((1.0f - gain) * outgoing[n]);
    }
}

/* Under a tone's first frame, the last synthesized frame repeated, to be
 * faded out. Only scratch copies are synthesized: a tone frame leaves the
 * noise generator and the synthesis state where they stood. */
static void
ambe2450_speech_under_tone(float* speech, const mbe_parms* cur_mp, const mbe_parms* prev_mp_enhanced) {
    mbe_parms repeated = *cur_mp;
    mbe_parms last = *prev_mp_enhanced;
    mbe_repeat_load_model(&repeated, &last);
    repeated.repeatCount = 0;
    mbe_synthesizeRepeatedSpeechf(speech, &repeated, &last, MBE_VOICED_SYNTHESIS_DVSI);
}

/**
 * Tone frame with a usable index (7.3): output the tone instead of
 * synthesized speech. The prediction history is not touched (4.3), nor the
 * kept model (save its unvoiced overlap), and the noise generator stands
 * still. Returns the word the frame leaves.
 */
static int
ambe2450_process_tone(float* aout_buf, mbe_process_result* result, const char ambe_d[49],
                      const struct ambe2450_snapshot* in, mbe_parms* cur_mp, mbe_parms* prev_mp,
                      mbe_parms* prev_mp_enhanced) {
    float speech[160] = {0};

    mbe_result_set_flag(result, MBE_PROCESS_FLAG_TONE);
    cur_mp->repeatCount = 0;
    mbe_copy_error_state(prev_mp, cur_mp);
    if (mbe_requiresMuting(cur_mp)) {
        return ambe2450_mute(aout_buf, result, in, 0, cur_mp, prev_mp_enhanced);
    }
    if (!in->held && !in->afresh) {
        ambe2450_speech_under_tone(speech, cur_mp, prev_mp_enhanced);
    }
    mbe_synthesizeTonef(aout_buf, ambe_d, cur_mp);
    if (!in->held) {
        ambe2450_fade_in_over(aout_buf, speech);
    }
    ambe2450_clear_overlap(cur_mp, prev_mp_enhanced);
    return AMBE2450_TONE_HELD | AMBE2450_STARTS_AFRESH | (ambe2450_read_tone_ad(ambe_d) << 8)
           | ambe2450_read_tone_id(ambe_d);
}

/* A frame the repeat criteria reject before it is classified (5.6). */
#define AMBE2450_FRAME_REJECTED (-1)

/* One dispatch point: each branch gets the snapshot, writes the frame's output
 * once and returns the word the frame leaves. */
static void
ambe2450_process_frame(float* aout_buf, mbe_process_result* result, const char ambe_d[49], int kind, mbe_parms* cur_mp,
                       mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    struct ambe2450_snapshot in;
    int word;

    ambe2450_take_snapshot(cur_mp, prev_mp_enhanced, &in);
    cur_mp->Vl[0] = 0;
    if (kind == MBE_AMBE2450_FRAME_TONE) {
        word = ambe2450_process_tone(aout_buf, result, ambe_d, &in, cur_mp, prev_mp, prev_mp_enhanced);
    } else if (kind == AMBE2450_FRAME_REJECTED) {
        word = ambe2450_process_invalid(aout_buf, result, 0u, &in, cur_mp, prev_mp, prev_mp_enhanced);
    } else if (kind == MBE_AMBE2450_FRAME_ERASURE) {
        /* A tone frame with an unusable index is flagged as both. */
        const unsigned tone = ambe2450_is_tone_frame(ambe_d) ? MBE_PROCESS_FLAG_TONE : 0u;
        word = ambe2450_process_invalid(aout_buf, result, MBE_PROCESS_FLAG_ERASURE | tone, &in, cur_mp, prev_mp,
                                        prev_mp_enhanced);
    } else {
        word = ambe2450_process_valid(aout_buf, result, kind == MBE_AMBE2450_FRAME_SILENCE, &in, cur_mp, prev_mp,
                                      prev_mp_enhanced);
    }
    cur_mp->Vl[0] = 0;
    prev_mp_enhanced->Vl[0] = word;
}

static int
mbe_processAmbe2450Dataf_internal(float* aout_buf, mbe_process_result* result, const char ambe_d[49], mbe_parms* cur_mp,
                                  mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    int kind;
    int c0_errors;
    int c0_errors_valid;
    int total_errors;
    int ret;

    /* Aliased sets would let the kept model or a frame's own model overwrite the history. */
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
        ambe2450_process_frame(aout_buf, result, ambe_d, AMBE2450_FRAME_REJECTED, cur_mp, prev_mp, prev_mp_enhanced);
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
