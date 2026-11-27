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
 * @brief IMBE 7200x4400 parameter decode, ECC, and synthesis hooks.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "imbe4400_internal.h"
#include "imbe7200x4400_const.h"
#include "mbe_adaptive.h"
#include "mbe_bitpack.h"
#include "mbe_compiler.h"
#include "mbe_math.h"
#include "mbe_pulsed.h"
#include "mbe_repeat.h"
#include "mbe_result.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"

/**
 * @brief Thread-local cache for IMBE DCT cosine coefficients.
 *
 * Pre-computes the cosine terms used in the Ri inverse DCT and per-block
 * inverse DCT loops to eliminate repeated cosf() calls per frame.
 *
 * Ri DCT: cosf(M_PI * (m-1) * (i-0.5) / 6) for m=1..6, i=1..6
 * Per-block IDCT: cosf(M_PI * (k-1) * (j-0.5) / ji) for ji=1..10, j=1..ji, k=1..ji
 */
struct imbe_dct_cache {
    int inited;
    float ri_cos[7][7];         /* [m][i] for m=1..6, i=1..6 (index 0 unused) */
    float idct_cos[11][11][11]; /* [ji][j][k] for ji=1..10, j=1..ji, k=1..ji */
};

static MBE_THREAD_LOCAL struct imbe_dct_cache imbe_cache = {0};

/**
 * @brief Initialize or return the thread-local IMBE DCT cache.
 *
 * Fills the cosine tables on first use. Because the cache is thread-local,
 * no locking is required.
 *
 * @return Pointer to the initialized cache.
 */
static struct imbe_dct_cache*
imbe_get_dct_cache(void) {
    if (imbe_cache.inited) {
        return &imbe_cache;
    }

    /* Fill Ri DCT cosine table: cosf(M_PI * (m-1) * (i-0.5) / 6) */
    for (int m = 1; m <= 6; m++) {
        for (int i = 1; i <= 6; i++) {
            imbe_cache.ri_cos[m][i] = cosf((M_PI * (float)(m - 1) * ((float)i - 0.5f)) / 6.0f);
        }
    }

    /* Fill per-block IDCT cosine table: cosf(M_PI * (k-1) * (j-0.5) / ji) */
    for (int ji = 1; ji <= 10; ji++) {
        for (int j = 1; j <= ji; j++) {
            for (int k = 1; k <= ji; k++) {
                imbe_cache.idct_cos[ji][j][k] = cosf((M_PI * (float)(k - 1) * ((float)j - 0.5f)) / (float)ji);
            }
        }
    }

    imbe_cache.inited = 1;
    return &imbe_cache;
}

/* The pitch codes DVSI's P25 encoder gives a frame with no voiced band: b0 25
 * (a period of 32.25 samples, 14 harmonics) for most, and b0 128..183 for the
 * rest. */
#define IMBE_UNVOICED_CODES_FIRST 128
#define IMBE_UNVOICED_CODES_LAST  183
/* The last pitch code of a voice frame. */
#define IMBE_LAST_PITCH_CODE      207

/* Which pitch a decoded pitch code plays (imbe_decode_fundamental()). */
enum imbe_pitch {
    IMBE_PITCH_DVSI,    /* P25: the pitch DVSI's P25 decoder plays for b0, eq 46 on a 2^-19 grid */
    IMBE_PITCH_NOMINAL, /* ProVoice: TIA-102.BABA eq 46 as it stands */
};

/* Eq 46's nominal fundamental of pitch code b0 (0..207), in float: rounded to
 * float here, where the contract names it, as under x87 excess precision the
 * quotient could otherwise reach what follows unrounded (mbe_round_to_float(),
 * mbe_math.h). */
static float
imbe_nominal_fundamental(int b0) {
    return mbe_round_to_float((float)(4 * M_PI) / (float)((float)b0 + 39.5));
}

/*
 * The fundamental DVSI's P25 decoder plays for pitch code b0 (0..207;
 * docs/dvsi-parity.md, "P25 pitch"): eq 46's f0 = 4 / (2 b0 + 79) cycles per
 * sample cut to a grid of 2^-19 cycle, floor(2^21 / (2 b0 + 79)) / 2^19, as a
 * w0 in radians per sample. The grid was found first on its decodes of tones
 * its encoder sends as voice frames (P25 tones sent as voice frames) and then
 * on speech, once the mid-frame pitch (mbe_subframe.c) had been told from a
 * frame's own: the per-code drift of the pulse train over runs at one code,
 * which the half-octave law this replaces was fitted to, is that of the
 * mid-frame pitch. The floor of an integer quotient, then 2 pi times it over
 * 2^19 formed in double and rounded to float once: over the 208 codes that
 * product (with M_PI's double value) lies at least 4.9e-11 of itself from a
 * point halfway between two floats, so double arithmetic, x87's extended one
 * and fast-math's reordering (2 pi / 2^19 is exact) all round it to the same
 * float.
 */
static float
imbe_grid_fundamental(int b0) {
    const double steps = floor(2097152.0 / (double)((2 * b0) + 79));
    return mbe_round_to_float((float)(2.0 * M_PI * steps / 524288.0));
}

/* Eq 47's harmonic count of a nominal fundamental. */
static int
imbe_harmonic_count(float nominal) {
    return (int)(0.9254 * (int)((M_PI / nominal) + 0.25));
}

/*
 * b0 and the fundamental, harmonic count and band count it gives (TIA-102.BABA
 * eqs 46-48). The nominal pitch of eq 46 is computed in float, as this decoder
 * always has (4 pi rounded to float, divided by b0 + 39.5 in float, the
 * quotient rounded to float also where the arithmetic keeps excess
 * precision), and L and K, and with them the bit layout, follow it for every
 * pitch code. The pitch stored for synthesis is the one DVSI's P25 decoder
 * plays, eq 46 on the 2^-19 grid (imbe_grid_fundamental()): from 0 to 2^-19
 * cycle below eq 46's, a relative 0 to 2.26e-4 (b0 201), and the same value
 * for speech and for a tone its encoder sent as voice frames
 * (imbe4400_play_tone(), whose model keeps the pitch of its code). ProVoice
 * (IMBE 7100x4400), which reaches this decoder converted to the 7200 layout,
 * has no DVSI reference and keeps the nominal pitch (`pitch`).
 *
 * The harmonic at or next above 505 Hz (MBE_VOICING_REACH_HZ; ceil(505 Hz /
 * f0)), at which synthesis samples the first edge of the 500 Hz columns P25's
 * rendering reads its voicing in (mbe_subframe.c), is the nominal pitch's for
 * every code. Of the 208 x 14 pairs of code and edge (505 + 500 k Hz, k 0..6,
 * and the same edges 250 Hz lower), the harmonic at or next above the edge at
 * the pitch played differs from the nominal pitch's at 2; every edge lies at least
 * 5.97e-5 of a harmonic from a whole number (b0 152 at 2005 Hz;
 * test_p25_onset_edges in test_subframe_synthesis.c), so float rounding moves
 * none. Discrete steps elsewhere in synthesis (a band edge's transform bin,
 * the mid-frame model's harmonic count, the harmonic where the high-band shelf
 * starts) follow the pitch played.
 *
 * imbe_code_pitch() computes them for pitch code b0 (0..207); it returns 1,
 * with the nominal pitch in *w0 and no L or K, for a harmonic count outside
 * 9..56. One body (not inlined) serves the decode, imbe4400_take_tone(), which
 * asks a tone word's model for the pitch of its code, and the pulsed state's
 * code tests (mbe_imbe4400_pulsed_code(), mbe_imbe4400_pitch_code()), so that
 * they round alike in every build.
 */
static MBE_NOINLINE int
imbe_code_pitch(int b0, enum imbe_pitch pitch, float* w0, int* L, int* K) {
    const float nominal = imbe_nominal_fundamental(b0);
    *w0 = nominal;

    const int harmonics = imbe_harmonic_count(nominal);
    if ((harmonics > 56) || (harmonics < 9)) {
        return 1;
    }
    *L = harmonics;
    *K = (harmonics < 37) ? (int)((float)(harmonics + 2) / (float)3) : 12;
    if (pitch == IMBE_PITCH_DVSI) {
        *w0 = imbe_grid_fundamental(b0);
    }
    return 0;
}

/*
 * The reserved code (b0 128..183, IMBE_UNVOICED_CODES_*) whose pitch w0 is, as
 * mbe_pulsed.h asks: the pitch imbe_code_pitch() gives that code for P25, the
 * function the decode stores it with, compared by bit pattern. The one code
 * that can match is the one nearest 4 pi / w0 - 39.5: the played pitch lies
 * up to a relative 2.26e-4 below eq 46's (imbe_grid_fundamental()), so that
 * quotient lies from 0 to 0.046 above the code over b0 128..183 (0.054 at most
 * over all 208 codes, b0 201). w0 is
 * positive and normal, so the quotient is finite; it is bounded before it is
 * converted.
 */
int
mbe_imbe4400_pulsed_code(float w0) {
    const double position = ((4.0 * M_PI) / (double)w0) - 39.5;
    if (!(position > (double)IMBE_UNVOICED_CODES_FIRST - 0.5 && position < (double)IMBE_UNVOICED_CODES_LAST + 0.5)) {
        return -1;
    }
    const int b0 = (int)floor(position + 0.5);
    float pitch;
    int L = 0, K = 0;
    if (imbe_code_pitch(b0, IMBE_PITCH_DVSI, &pitch, &L, &K) != 0) {
        return -1;
    }
    uint32_t want, have;
    memcpy(&want, &pitch, sizeof(want));
    memcpy(&have, &w0, sizeof(have));
    return (want == have) ? b0 - IMBE_UNVOICED_CODES_FIRST : -1;
}

/*
 * The P25 pitch code (b0 0..207) whose pitch w0 is, as mbe_pulsed.h asks: the
 * pitch imbe_code_pitch() gives that code for P25, by bit pattern, or -1. As
 * for mbe_imbe4400_pulsed_code(), the one code that can match is the one
 * nearest 4 pi / w0 - 39.5: over the 208 codes the decoded pitch's quotient
 * lies from 0 to 0.054 above its own code (b0 201). w0 is positive and normal, so the quotient
 * is finite; it is bounded before it is converted.
 */
int
mbe_imbe4400_pitch_code(float w0) {
    const double position = ((4.0 * M_PI) / (double)w0) - 39.5;
    if (!(position > -0.5 && position < (double)IMBE_LAST_PITCH_CODE + 0.5)) {
        return -1;
    }
    const int b0 = (int)floor(position + 0.5);
    float pitch;
    int L = 0, K = 0;
    if (imbe_code_pitch(b0, IMBE_PITCH_DVSI, &pitch, &L, &K) != 0) {
        return -1;
    }
    uint32_t want, have;
    memcpy(&want, &pitch, sizeof(want));
    memcpy(&have, &w0, sizeof(have));
    return (want == have) ? b0 : -1;
}

static int
imbe_decode_fundamental(const char* imbe_d, enum imbe_pitch pitch, mbe_parms* cur_mp, int* out_b0, int* out_L9) {
    static const unsigned b0_indices[] = {0u, 1u, 2u, 3u, 4u, 5u, 85u, 86u};
    int b0 = mbe_bits_by_index_to_int(imbe_d, b0_indices, sizeof(b0_indices) / sizeof(b0_indices[0]));
    *out_b0 = b0;
    if (b0 > 207) {
#ifdef IMBE_DEBUG
        if ((b0 >= 216) && (b0 <= 219)) {
            fprintf(stderr, "Silence\n");
        } else {
            fprintf(stderr, "Invalid fundamental frequency\n");
        }
#endif
        return 1;
    }

    int L = 0, K = 0;
    const int invalid = imbe_code_pitch(b0, pitch, &cur_mp->w0, &L, &K);
    /* imbe_code_pitch() returns 0 only with L in 9..56; the bound is tested
     * here again, where the bit layout's tables are indexed by L - 9 (a static
     * analyzer does not look into a function that is not inlined). */
    if (invalid || (L > 56) || (L < 9)) {
#ifdef IMBE_DEBUG
        fprintf(stderr, "invalid L at b0 %i\n", b0);
#endif
        return 1;
    }
    cur_mp->L = L;
    cur_mp->K = K;
    *out_L9 = L - 9;

#ifdef IMBE_DEBUG
    fprintf(stderr, "b0:%i L:%i K:%i\n", b0, L, cur_mp->K);
#endif
    return 0;
}

static void
imbe_read_bit_layout(const char* imbe_d, int L9, char bb[58][12]) {
    const int* bo1 = bo[L9][0];
    const int* bo2 = bo1 + 1;
    for (int i = 6; i < 85; i++) {
        bb[*bo1][*bo2] = imbe_d[i];
#ifdef IMBE_DEBUG
        fprintf(stderr, "bo1: %i,bo2: %i, ", *bo1, *bo2);
#endif
        bo1 += 2;
        bo2 += 2;
    }
}

/*
 * Reads the band voicing, and marks what it tells DVSI's decoder beyond the
 * bands themselves (docs/dvsi-parity.md, "Voiced onsets and offsets" and "P25:
 * the bands that start with the frame"):
 *  - no voiced band, at a pitch code other than the ones DVSI's encoder sends
 *    only with such frames (IMBE_UNVOICED_CODES_*): the voicing of the frame
 *    before it lasts to mid-frame in the bands whose 500 Hz column, or the
 *    column 250 Hz higher, that frame covers (MBE_FRAME_MARK_PITCH_ABSENT;
 *    the pulse train runs on through the frame only when a band is held); at
 *    one of those codes no mark: the synthesizer reads the code from the
 *    model's pitch as the pulsed state (mbe_imbe4400_pulsed_code(),
 *    mbe_pulsed.h), P25's rendering only, while its validity counter allows;
 *  - voiced bands: MBE_FRAME_MARK_ONSET_COLUMNS, which P25's rendering gives
 *    every voiced model it renders but a tone model, after adaptive smoothing
 *    (the model played for a P25 tone sent as voice frames is marked
 *    MBE_FRAME_MARK_TONE instead, imbe4400_play_tone()). It chooses no band
 *    since the mid-frame columns change: P25's reads the mid-frame voicing of
 *    every voiced model but a tone model in 500 Hz columns, whatever its mark
 *    (mbe_subframe.c; docs/dvsi-parity.md,
 *    "D-STAR and P25: mid-frame voicing by columns"), which contains the onset
 *    rule the mark gave before.
 */
static void
imbe_decode_voicing(mbe_parms* cur_mp, int b0, char bb[58][12]) {
    enum mbe_frame_mark mark = MBE_FRAME_MARK_NONE;
    int j = 1;
    int k = cur_mp->K - 1;
    int top = 0;
    for (int i = 1; i <= cur_mp->L; i++) {
        /* Cast via unsigned to avoid sign-extension when 'char' is signed */
        cur_mp->Vl[i] = (int)(unsigned char)bb[1][k];
        if (cur_mp->Vl[i]) {
            top = i;
        }
        if (j == 3) {
            j = 1;
            if (k > 0) {
                k--;
            } else {
                k = 0;
            }
        } else {
            j++;
        }
    }
    if (top == 0) {
        if (b0 < IMBE_UNVOICED_CODES_FIRST || b0 > IMBE_UNVOICED_CODES_LAST) {
            mark = MBE_FRAME_MARK_PITCH_ABSENT;
        }
    } else {
        mark = MBE_FRAME_MARK_ONSET_COLUMNS;
    }
    mbe_set_frame_mark(cur_mp, mark);
}

static void
imbe_decode_gains(char bb[58][12], int L9, float Gm[7]) {
    int b2 = mbe_bits_descending_to_int(bb[2], 5, 0);
    Gm[1] = B2[b2];
#ifdef IMBE_DEBUG
    fprintf(stderr, "b2: %i g1: %e\n", b2, Gm[1]);
#endif

    const float* ba1 = ba[L9][0];
    const float* ba2 = ba1 + 1;
    for (int i = 2; i < 7; i++) {
        int bm = mbe_bits_descending_to_int(bb[i + 1], (int)*ba1 - 1, 0);
        Gm[i] = (*ba2 * ((float)bm - exp2f((*ba1) - 1.0f) + 0.5f));
#ifdef IMBE_DEBUG
        fprintf(stderr, "G%i: %e, %i, ba1: %e, ba2: %e\n", i, Gm[i], bm, *ba1, *ba2);
#endif
        ba1 += 2;
        ba2 += 2;
    }
}

static void
imbe_compute_ri(const float Gm[7], float Ri[7]) {
    const struct imbe_dct_cache* cache = imbe_get_dct_cache();
    for (int i = 1; i <= 6; i++) {
        float sum = 0;
        for (int m = 1; m <= 6; m++) {
            int am = (m == 1) ? 1 : 2;
            sum = sum + ((float)am * Gm[m] * cache->ri_cos[m][i]);
#ifdef IMBE_DEBUG
            fprintf(stderr, "sum: %e ", sum);
#endif
        }
        Ri[i] = sum;
#ifdef IMBE_DEBUG
        fprintf(stderr, "R%i: %e\n", i, Ri[i]);
#endif
    }
#ifdef IMBE_DEBUG
    fprintf(stderr, "R1: %e\n", Ri[1]);
#endif
}

static void
imbe_decode_hoc_coefficients(char bb[58][12], int L9, const float Ri[7], float Cik[7][11]) {
    int m = 8;
    for (int i = 1; i <= 6; i++) {
        Cik[i][1] = Ri[i];
        for (int k = 2; k <= ImbeJi[L9][i - 1]; k++) {
            int Bm = hoba[L9][m - 8];
            if (Bm <= 0) {
                Cik[i][k] = 0;
            } else {
                int bm = mbe_bits_descending_to_int(bb[m], Bm - 1, 0);
                Cik[i][k] = ((quantstep[Bm - 1] * standdev[k - 2]) * (((float)bm - exp2f((float)Bm - 1.0f)) + 0.5f));
            }
            m++;
        }
    }
}

static void
imbe_inverse_dct_tl(float Cik[7][11], int L9, float Tl[57]) {
    const struct imbe_dct_cache* cache = imbe_get_dct_cache();
    int l = 1;
    for (int i = 1; i <= 6; i++) {
        int ji = ImbeJi[L9][i - 1];
        for (int j = 1; j <= ji; j++) {
            float sum = 0;
            for (int k = 1; k <= ji; k++) {
                int ak = (k == 1) ? 1 : 2;
                sum = sum + ((float)ak * Cik[i][k] * cache->idct_cos[ji][j][k]);
            }
            Tl[l] = sum;
            l++;
        }
    }
#ifdef IMBE_DEBUG
    fprintf(stderr, "T1: %e\n", Tl[1]);
#endif
}

static float
imbe_spectral_rho(int L) {
    if (L <= 15) {
        return 0.4f;
    }
    if (L <= 24) {
        return (0.03f * (float)L) - 0.05f;
    }
    return 0.7f;
}

static void
imbe_update_spectral_amplitudes(mbe_parms* cur_mp, mbe_parms* prev_mp, const float Tl[57], float rho) {
    int intkl[57];
    float flokl[57], deltal[57];
    int cur_L = mbe_clamp_harmonic_count(cur_mp->L);
    int prev_L = mbe_clamp_harmonic_count(prev_mp->L);

    cur_mp->L = cur_L;

    if (cur_L > prev_L) {
        for (int l = prev_L + 1; l <= cur_L; l++) {
            prev_mp->Ml[l] = prev_mp->Ml[prev_L];
            prev_mp->log2Ml[l] = prev_mp->log2Ml[prev_L];
        }
    }
    prev_mp->log2Ml[0] = prev_mp->log2Ml[1];
    prev_mp->Ml[0] = prev_mp->Ml[1];

    float Sum77 = 0;
    for (int l = 1; l <= cur_L; l++) {
        flokl[l] = mbe_prediction_position(prev_L, cur_L, l);
        intkl[l] = (int)flokl[l];
        if (intkl[l] < 0) {
            intkl[l] = 0;
        } else if (intkl[l] > MBE_MAX_HARMONIC_BANDS) {
            intkl[l] = MBE_MAX_HARMONIC_BANDS;
        }
#ifdef IMBE_DEBUG
        fprintf(stderr, "flokl: %e, intkl: %i ", flokl[l], intkl[l]);
#endif
        deltal[l] = flokl[l] - (float)intkl[l];
#ifdef IMBE_DEBUG
        fprintf(stderr, "deltal: %e ", deltal[l]);
#endif
        int upper = intkl[l] + 1;
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }
        Sum77 = Sum77 + ((((float)1 - deltal[l]) * prev_mp->log2Ml[intkl[l]]) + (deltal[l] * prev_mp->log2Ml[upper]));
    }
    Sum77 = ((rho / (float)cur_L) * Sum77);

#ifdef IMBE_DEBUG
    fprintf(stderr, "Sum77: %e\n", Sum77);
#endif

    for (int l = 1; l <= cur_L; l++) {
        int upper = intkl[l] + 1;
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }
        float c1 = (rho * ((float)1 - deltal[l]) * prev_mp->log2Ml[intkl[l]]);
        float c2 = (rho * deltal[l] * prev_mp->log2Ml[upper]);
        cur_mp->log2Ml[l] = Tl[l] + c1 + c2 - Sum77;
        cur_mp->Ml[l] = exp2f(cur_mp->log2Ml[l]);
#ifdef IMBE_DEBUG
        fprintf(stderr, "rho: %e c1: %e c2: %e Sum77: %e T%i: %e log2M%i: %e M%i: %e\n", rho, c1, c2, Sum77, l, Tl[l],
                l, cur_mp->log2Ml[l], l, cur_mp->Ml[l]);
#endif
    }
}

/**
 * @brief Print IMBE 4400 parameter bits to stderr (debug aid).
 * @param imbe_d IMBE parameter bits (88).
 */
void
mbe_dumpImbe4400Data(const char* imbe_d) {

    int i;
    const char* imbe;

    imbe = imbe_d;
    for (i = 0; i < 88; i++) {
        fprintf(stderr, "%i", *imbe);
        imbe++;
    }
}

/**
 * @brief Print IMBE 7200x4400 parameter bits to stderr (debug aid).
 * @param imbe_d IMBE parameter bits (88).
 */
void
mbe_dumpImbe7200x4400Data(const char* imbe_d) {

    int i;
    const char* imbe;

    imbe = imbe_d;
    for (i = 0; i < 88; i++) {
        if ((i == 12) || (i == 24) || (i == 36) || (i == 48) || (i == 59) || (i == 70) || (i == 81)) {
            fprintf(stderr, " ");
        }
        fprintf(stderr, "%i", *imbe);
        imbe++;
    }
}

/**
 * @brief Print raw IMBE 7200x4400 frame bitplanes to stderr.
 * @param imbe_fr Frame as 8x23 bitplanes (last row partial).
 */
void
mbe_dumpImbe7200x4400Frame(const char imbe_fr[8][23]) {

    int i, j;

    for (i = 0; i < 4; i++) {
        for (j = 22; j >= 0; j--) {
            fprintf(stderr, "%i", imbe_fr[i][j]);
        }
        fprintf(stderr, " ");
    }
    for (i = 4; i < 7; i++) {
        for (j = 14; j >= 0; j--) {
            fprintf(stderr, "%i", imbe_fr[i][j]);
        }
        fprintf(stderr, " ");
    }
    for (j = 6; j >= 0; j--) {
        fprintf(stderr, "%i", imbe_fr[7][j]);
    }
}

/**
 * @brief Apply ECC to IMBE 7200x4400 C0 and update in-place.
 * @param imbe_fr Frame as 8x23 bitplanes.
 * @return Number of corrected errors in C0.
 */
int
mbe_eccImbe7200x4400C0(char imbe_fr[8][23]) {

    int j, errs;
    char in[23], out[23];
    int ret = mbe_validate_bits((const char*)imbe_fr, (size_t)8u * 23u);
    if (ret < 0) {
        return ret;
    }

    for (j = 0; j < 23; j++) {
        in[j] = imbe_fr[0][j];
    }
    errs = mbe_golay2312(in, out);
    for (j = 0; j < 23; j++) {
        imbe_fr[0][j] = out[j];
    }

    return errs;
}

static int
mbe_eccImbe7200x4400C0Soft(mbe_soft_bit imbe_fr[8][23]) {
    int j, errs;
    mbe_soft_bit in[23];
    char out[23];

    for (j = 0; j < 23; j++) {
        in[j] = imbe_fr[0][j];
    }
    errs = mbe_golay2312Soft(in, out);
    for (j = 0; j < 23; j++) {
        imbe_fr[0][j].bit = (uint8_t)(out[j] & 1);
    }

    return errs;
}

/**
 * @brief Internal: Apply ECC to IMBE 7200x4400 data with separate C4 error tracking.
 * @param imbe_fr Frame as 8x23 bitplanes.
 * @param imbe_d  Output parameter bits (88).
 * @param errs_c4 Output: errors in C4 (Hamming coset 4), or NULL if not needed.
 * @return Number of corrected errors in protected fields.
 */
static int
mbe_eccImbe7200x4400DataInternal(char imbe_fr[8][23], char* imbe_d, int* errs_c4) {

    int i, j, errs;
    char *imbe, gin[23], gout[23], hin[15], hout[15];

    errs = 0;
    imbe = imbe_d;
    for (i = 0; i < 4; i++) {
        if (i > 0) {
            for (j = 0; j < 23; j++) {
                gin[j] = imbe_fr[i][j];
            }
            errs += mbe_golay2312(gin, gout);
            for (j = 22; j > 10; j--) {
                *imbe = gout[j];
                imbe++;
            }
        } else {
            for (j = 22; j > 10; j--) {
                *imbe = imbe_fr[i][j];
                imbe++;
            }
        }
    }
    for (i = 4; i < 7; i++) {
        for (j = 0; j < 15; j++) {
            hin[j] = imbe_fr[i][j];
        }
        int hamming_errs = mbe_hamming1511(hin, hout);
        errs += hamming_errs;
        /* Track C4 (first Hamming coset) errors separately for adaptive smoothing */
        if (i == 4 && errs_c4 != NULL) {
            *errs_c4 = hamming_errs;
        }
        for (j = 14; j >= 4; j--) {
            *imbe = hout[j];
            imbe++;
        }
    }
    for (j = 6; j >= 0; j--) {
        *imbe = imbe_fr[7][j];
        imbe++;
    }

    return errs;
}

static int
mbe_eccImbe7200x4400DataSoftInternal(mbe_soft_bit imbe_fr[8][23], char* imbe_d, int* errs_c4) {
    int i, j, errs;
    char *imbe, gout[23], hout[15];
    mbe_soft_bit gin[23], hin[15];

    errs = 0;
    imbe = imbe_d;
    for (i = 0; i < 4; i++) {
        if (i > 0) {
            for (j = 0; j < 23; j++) {
                gin[j] = imbe_fr[i][j];
            }
            errs += mbe_golay2312Soft(gin, gout);
            for (j = 22; j > 10; j--) {
                *imbe = gout[j];
                imbe++;
            }
        } else {
            for (j = 22; j > 10; j--) {
                *imbe = (char)(imbe_fr[i][j].bit & 1u);
                imbe++;
            }
        }
    }
    for (i = 4; i < 7; i++) {
        for (j = 0; j < 15; j++) {
            hin[j] = imbe_fr[i][j];
        }
        int hamming_errs = mbe_hamming1511Soft(hin, hout);
        errs += hamming_errs;
        if (i == 4 && errs_c4 != NULL) {
            *errs_c4 = hamming_errs;
        }
        for (j = 14; j >= 4; j--) {
            *imbe = hout[j];
            imbe++;
        }
    }
    for (j = 6; j >= 0; j--) {
        *imbe = (char)(imbe_fr[7][j].bit & 1u);
        imbe++;
    }

    return errs;
}

/**
 * @brief Apply ECC to IMBE 7200x4400 data and pack parameter bits.
 * @param imbe_fr Frame as 8x23 bitplanes.
 * @param imbe_d  Output parameter bits (88).
 * @return Number of corrected errors in protected fields.
 */
int
mbe_eccImbe7200x4400Data(char imbe_fr[8][23], char* imbe_d) {
    if (!imbe_d) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    int ret = mbe_validate_bits((const char*)imbe_fr, (size_t)8u * 23u);
    if (ret < 0) {
        return ret;
    }
    return mbe_eccImbe7200x4400DataInternal(imbe_fr, imbe_d, NULL);
}

/*
 * P25 tones sent as voice frames. P25 has no tone frames: DVSI's encoder
 * sends a tone it detects as voice frames built from TIA-102.BABA-1 Annex J
 * (the pitch code b0 = floor(4 pi / w0 - 39) of the tone's fundamental, the
 * tone's harmonics l1 and l2 large, and a voicing word b1 that is the same in
 * every steady frame of the tone: 12,555 steady frames compared, while 112
 * frames at tone boundaries carry another pair), and
 * DVSI's decoder recognises such a frame and plays a tone model in place of
 * its decoded harmonics, as US 7,970,606 describes (a frame whose fundamental
 * is a table entry's, whose tone harmonics are voiced and dominate the other
 * magnitudes, is a tone frame). The rows: the distinct (b0, b1, l1, l2) of the
 * DTMF (128..143), KNOX (144..159) and single tones (5..122) as DVSI's encoder
 * sends them, b1 the K band voicing bits as transmitted, band 1 the most
 * significant. The harmonics are Annex J's (a single tone's l1 = l2 is Annex
 * J's harmonic number: 1 for IDs 5..12, 2 for 13..25, 3 for 26..38, 4 for
 * 39..51, 5 for 52..64, 6 for 65..76, 7 for 77..89, 8 for 90..102, 9 for
 * 103..115 and 10 for 116..122; not the ID / 16 + 1 by which the rate-33 and
 * D-STAR decoders play a single tone, mbe_tone.h). The pitch code is that of
 * the fundamental Annex J gives: a dual tone's printed f0 (DVSI's own 150.80
 * and 86.47 Hz for 145 and 147, mbe_tone.h), and for a single tone the
 * quotient 31.25 ID / l. For two ranges Annex J prints a multiplier instead:
 * 10.417 ID for 26..38 (31.25 / 3 rounded) and 5.2803 ID for 65..76 (31.25 /
 * 6 is 5.2083); built from the printed values, the codes of IDs 32, 65, 66,
 * 68, 69, 71, 73, 74 and 76 would come out otherwise. So
 * built, the codes equal those in DVSI's bit streams for all 154 tones. The
 * voicing word, which is not the bands of l1 and l2 under the decoder's band
 * mapping, was read from DVSI's bit streams of its tone vectors. Single tones on one code
 * and word share a row's pair at several harmonics (b0 17, 1000: l 1, 2 and
 * 3). Left out: three single tones whose harmonic lies past the frame's L
 * (120..122 at b0 2 and 3, harmonic 10 of 9), and the call-progress tones
 * (160..163), whose frames DVSI's decoder does not play as equal-amplitude
 * tone models (none of the 50 frames selected), and which play as decoded.
 * The rows are sorted by b0, and imbe_tone_first_row indexes them by it.
 */
struct imbe_tone_row {
    unsigned char b0;
    unsigned short b1;
    unsigned char l1;
    unsigned char l2;
};

static const struct imbe_tone_row imbe_tone_rows[] = {
    {1, 0x001, 7, 7},     {1, 0x001, 8, 8},     {1, 0x001, 9, 9},     {1, 0x002, 4, 4},     {1, 0x002, 5, 5},
    {1, 0x002, 6, 6},     {1, 0x004, 2, 2},     {1, 0x004, 3, 3},     {2, 0x001, 7, 7},     {2, 0x001, 8, 8},
    {2, 0x001, 9, 9},     {2, 0x002, 4, 4},     {2, 0x002, 5, 5},     {2, 0x002, 6, 6},     {2, 0x004, 3, 3},
    {3, 0x001, 7, 7},     {3, 0x001, 8, 8},     {3, 0x001, 9, 9},     {3, 0x002, 4, 4},     {3, 0x002, 5, 5},
    {3, 0x002, 6, 6},     {3, 0x004, 1, 1},     {3, 0x004, 2, 2},     {3, 0x004, 3, 3},     {4, 0x001, 10, 10},
    {4, 0x002, 7, 7},     {4, 0x003, 8, 8},     {4, 0x003, 9, 9},     {4, 0x004, 4, 4},     {4, 0x004, 5, 5},
    {4, 0x006, 6, 6},     {4, 0x00C, 3, 3},     {5, 0x001, 10, 10},   {5, 0x002, 7, 7},     {5, 0x003, 8, 8},
    {5, 0x003, 9, 9},     {5, 0x004, 4, 4},     {5, 0x004, 5, 5},     {5, 0x006, 6, 6},     {5, 0x008, 2, 2},
    {6, 0x002, 7, 7},     {6, 0x002, 8, 8},     {6, 0x004, 4, 4},     {6, 0x004, 5, 5},     {6, 0x006, 6, 6},
    {6, 0x00C, 3, 3},     {7, 0x002, 7, 7},     {7, 0x004, 4, 4},     {7, 0x004, 5, 5},     {7, 0x006, 6, 6},
    {7, 0x008, 1, 1},     {7, 0x008, 2, 2},     {7, 0x00C, 3, 3},     {8, 0x004, 4, 4},     {8, 0x004, 5, 5},
    {8, 0x006, 6, 6},     {9, 0x004, 4, 4},     {9, 0x004, 5, 5},     {9, 0x008, 2, 2},     {9, 0x00C, 3, 3},
    {10, 0x004, 4, 4},    {10, 0x004, 5, 5},    {10, 0x00C, 3, 3},    {12, 0x004, 4, 4},    {12, 0x008, 1, 1},
    {12, 0x008, 2, 2},    {12, 0x00C, 3, 3},    {13, 0x004, 4, 4},    {13, 0x00C, 3, 3},    {14, 0x008, 2, 2},
    {15, 0x00C, 3, 3},    {17, 0x008, 1, 1},    {17, 0x008, 2, 2},    {17, 0x008, 3, 3},    {20, 0x010, 3, 3},
    {21, 0x010, 2, 2},    {25, 0x010, 1, 1},    {25, 0x010, 2, 2},    {29, 0x00C, 4, 7},    {29, 0x010, 2, 2},
    {34, 0x020, 1, 1},    {34, 0x020, 2, 2},    {36, 0x038, 4, 7},    {39, 0x020, 2, 2},    {44, 0x070, 4, 7},
    {46, 0x040, 1, 1},    {53, 0x070, 4, 7},    {56, 0x060, 4, 7},    {63, 0x080, 1, 1},    {67, 0x060, 4, 7},
    {80, 0x040, 7, 9},    {80, 0x060, 7, 11},   {92, 0x1E0, 7, 11},   {97, 0x380, 7, 9},    {97, 0x3C0, 7, 11},
    {98, 0x260, 6, 14},   {106, 0x3C0, 7, 11},  {112, 0x7C0, 7, 11},  {117, 0x7F0, 8, 14},  {118, 0x6F0, 6, 14},
    {127, 0x7C8, 8, 17},  {128, 0x7C0, 7, 11},  {133, 0x7B0, 8, 14},  {144, 0x7E4, 8, 17},  {146, 0x798, 7, 15},
    {151, 0x7E4, 8, 17},  {164, 0x0F4, 12, 17}, {175, 0x7F8, 10, 14}, {186, 0x331, 12, 23}, {186, 0x3F4, 12, 17},
    {189, 0x332, 11, 21}, {189, 0x3F2, 10, 19}, {195, 0x3F4, 12, 17}, {196, 0x312, 11, 21}, {196, 0x3E4, 9, 17},
    {197, 0x3F2, 10, 19},
};

/* The rows of pitch code b0 are imbe_tone_rows[imbe_tone_first_row[b0]] up to
 * imbe_tone_first_row[b0 + 1] (the rows above are sorted by b0): the first row
 * of each code 0..197 and the end; a code past 197 has none. Derived from the
 * rows: adding or removing a row moves every entry after its code. */
#define IMBE_TONE_LAST_CODE 197
static const unsigned char imbe_tone_first_row[IMBE_TONE_LAST_CODE + 2] = {
    0,   0,   8,   15,  24,  32,  40,  46,  53,  56,  60,  63,  63,  67,  69,  70,  71,  71,  74,  74,  74,  75,  76,
    76,  76,  76,  78,  78,  78,  78,  80,  80,  80,  80,  80,  82,  82,  83,  83,  83,  84,  84,  84,  84,  84,  85,
    85,  86,  86,  86,  86,  86,  86,  86,  87,  87,  87,  88,  88,  88,  88,  88,  88,  88,  89,  89,  89,  89,  90,
    90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  92,  92,  92,  92,  92,  92,  92,  92,  92,  92,  92,
    92,  93,  93,  93,  93,  93,  95,  96,  96,  96,  96,  96,  96,  96,  96,  97,  97,  97,  97,  97,  97,  98,  98,
    98,  98,  98,  99,  100, 100, 100, 100, 100, 100, 100, 100, 100, 101, 102, 102, 102, 102, 102, 103, 103, 103, 103,
    103, 103, 103, 103, 103, 103, 103, 104, 104, 105, 105, 105, 105, 105, 106, 106, 106, 106, 106, 106, 106, 106, 106,
    106, 106, 106, 106, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 108, 108, 108, 108, 108, 108, 108, 108,
    108, 108, 108, 110, 110, 110, 112, 112, 112, 112, 112, 112, 113, 115, 116,
};

static void
imbe_tone_set_field(int* field, int value) {
    if (field != NULL) {
        *field = value;
    }
}

int
mbe_imbe4400_tone_row(int i, int* b0, int* b1, int* l1, int* l2) {
    if (i < 0 || (size_t)i >= sizeof(imbe_tone_rows) / sizeof(imbe_tone_rows[0])) {
        return 0;
    }
    const struct imbe_tone_row* row = &imbe_tone_rows[i];
    imbe_tone_set_field(b0, row->b0);
    imbe_tone_set_field(b1, row->b1);
    imbe_tone_set_field(l1, row->l1);
    imbe_tone_set_field(l2, row->l2);
    return 1;
}

int
mbe_imbe4400_tone_first_row(int b0) {
    return (b0 >= 0 && b0 <= IMBE_TONE_LAST_CODE + 1) ? imbe_tone_first_row[b0] : -1;
}

/* Others over tone: a frame passes a row when the rms of its other decoded
 * magnitudes is under a share of the rms of the tone's. For a dual tone: DVSI's
 * tone frames reach at most 0.296 (99.9% under 0.242) and the frames of the
 * speech vectors whose (b0, b1) is a dual tone's start at 0.482. */
#define IMBE_TONE_SHARE_DUAL   0.375
/* For a single tone, whose pairs also hold speech frames of high voices with a
 * dominant fundamental. On the isolated speech frames of the 25 development and
 * validation vectors that pass a single-tone row, those with a share up to
 * 0.0234 (five) follow DVSI's output as well or better through the tone path,
 * and from 0.0236 on (about 65) the speech path follows it better in nearly
 * every case, with one exception at 0.0289: the limit lies between, not at a
 * clean split. 98.9% of the single-tone frames of DVSI's tone vectors lie under
 * it (median 0.003). */
#define IMBE_TONE_SHARE_SINGLE 0.0235
/* The word a tone frame's decode leaves in the otherwise unused Vl[0] (Vl is
 * indexed from 1, and P25 keeps nothing there): the flag, b0 in bits 16-23,
 * l2 in bits 8-15 and l1 in bits 0-7. */
#define IMBE_TONE_FLAG         MBE_P25_TONE_WORD_FLAG

/* The frame's voicing word as transmitted, band 1 the most significant bit:
 * band j holds harmonics 3j - 2 .. 3j, the last band the rest
 * (imbe_decode_voicing()). */
static unsigned
imbe_voicing_word(const mbe_parms* mp) {
    unsigned word = 0u;
    for (int band = 1; band <= mp->K; band++) {
        word = (word << 1) | (mp->Vl[(3 * band) - 2] != 0 ? 1u : 0u);
    }
    return word;
}

/* A decoded magnitude squared, in double, or -1 for one that is not finite
 * (a caller's history can make it so; read from the bit pattern). */
static double
imbe_power(const float* magnitude) {
    uint32_t bits;
    memcpy(&bits, magnitude, sizeof(bits));
    if ((bits & UINT32_C(0x7FFFFFFF)) >= UINT32_C(0x7F800000)) {
        return -1.0;
    }
    return (double)*magnitude * (double)*magnitude;
}

/* A positive, finite float that is not subnormal, read from its bit pattern
 * (fast-math folds isfinite()), as the synthesis reads a magnitude
 * (positive_normal() in mbe_subframe.c): anything else is no magnitude. */
static int
imbe_positive_normal(const float* x) {
    uint32_t bits;
    memcpy(&bits, x, sizeof(bits));
    return bits - UINT32_C(0x00800000) < UINT32_C(0x7F000000);
}

/* The share of the other harmonics against the tone's, squared: the mean of
 * M_l^2 over the other harmonics over its mean over the tone's (one harmonic
 * for a single tone); -1 if a magnitude is not finite or the tone's are 0. */
static double
imbe_tone_share2(const mbe_parms* mp, int l1, int l2) {
    double tone = 0.0, others = 0.0;
    int tones = 0;
    for (int l = 1; l <= mp->L; l++) {
        const double power = imbe_power(&mp->Ml[l]);
        if (power < 0.0) {
            return -1.0;
        }
        if (l == l1 || l == l2) {
            tone += power;
            tones++;
        } else {
            others += power;
        }
    }
    if (!(tone > 0.0)) {
        return -1.0;
    }
    return (others / (double)(mp->L - tones)) / (tone / (double)tones);
}

/*
 * The tone a decoded P25 voice frame carries, as the word to keep in Vl[0], or
 * 0: among the rows of the frame's pitch code and voicing word whose harmonics
 * are within L and whose tone magnitudes are positive, finite and normal
 * (imbe_positive_normal()), the one whose tone harmonics dominate most, if the
 * others' rms is under the row's share of the tone's (IMBE_TONE_SHARE_DUAL,
 * IMBE_TONE_SHARE_SINGLE; decoded magnitudes, before enhancement). The rows of
 * one pair are all dual or all single tones. Every voice frame passes here, so
 * a code with no row returns after one lookup and the scan covers only the
 * rows the index gives for the frame's code (each row's code still checked, so
 * an index that is off can drop a row but never add one). The synthesis runs
 * it again on the model a word sits in (imbe4400_take_tone()): one body, not
 * inlined, so that the two runs on the same model give the same word in every
 * build (fast-math or x87 could round an inlined copy otherwise).
 */
static MBE_NOINLINE int
imbe_tone_word(const mbe_parms* mp, int b0) {
    if (b0 < 0 || b0 > IMBE_TONE_LAST_CODE || imbe_tone_first_row[b0] == imbe_tone_first_row[b0 + 1]) {
        return 0;
    }
    const unsigned word = imbe_voicing_word(mp);
    double best = -1.0;
    int found = 0;
    for (unsigned i = imbe_tone_first_row[b0]; i < imbe_tone_first_row[b0 + 1]; i++) {
        const struct imbe_tone_row* row = &imbe_tone_rows[i];
        if (row->b0 != b0 || row->b1 != word || row->l2 > mp->L || !imbe_positive_normal(&mp->Ml[row->l1])
            || !imbe_positive_normal(&mp->Ml[row->l2])) {
            continue;
        }
        const double limit = (row->l1 == row->l2) ? IMBE_TONE_SHARE_SINGLE : IMBE_TONE_SHARE_DUAL;
        const double share = imbe_tone_share2(mp, row->l1, row->l2);
        if (share >= 0.0 && share < limit * limit && (found == 0 || share < best)) {
            best = share;
            found = IMBE_TONE_FLAG | (b0 << 16) | (row->l2 << 8) | row->l1;
        }
    }
    return found;
}

/* mbe_decodeImbe4400Parms(), the pitch played for b0 chosen by `pitch`. */
static int
imbe4400_decode_parms(const char* imbe_d, enum imbe_pitch pitch, mbe_parms* cur_mp, mbe_parms* prev_mp) {

    int b0;
    int L9;
    float Cik[7][11] = {{0.0f}}, Tl[57] = {0.0f}, Gm[7], Ri[7];
    char bb[58][12];
    int ret;

    if (!cur_mp || !prev_mp) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_validate_bits(imbe_d, 88u);
    if (ret < 0) {
        return ret;
    }
    /* Vl[0] holds a tone word only where this frame's decode recognises a
     * tone (below): 0 for every other frame it decodes, a special or invalid
     * frame included, whatever cur_mp held before (a status error above is no
     * frame and writes nothing). prev_mp's word is not touched (the prediction
     * writes only its taps there): the process functions repeat from it, its
     * word included. */
    cur_mp->Vl[0] = 0;

    if (imbe_decode_fundamental(imbe_d, pitch, cur_mp, &b0, &L9) != 0) {
        return 1;
    }

    // read bits from imbe_d into b0..bL+1
    imbe_read_bit_layout(imbe_d, L9, bb);

    // Vl
    imbe_decode_voicing(cur_mp, b0, bb);

    imbe_decode_gains(bb, L9, Gm);

    // inverse DCT Gi to give Ri (also known as Ci,1) - using cached cosines
    imbe_compute_ri(Gm, Ri);

    // load b8..bL+1 into Ci,k
    imbe_decode_hoc_coefficients(bb, L9, Ri, Cik);

    // inverse DCT each Ci,k to give ci,j (Tl) - using cached cosines
    imbe_inverse_dct_tl(Cik, L9, Tl);

    // determine log2Ml by applying ci,j to previous log2Ml
    imbe_update_spectral_amplitudes(cur_mp, prev_mp, Tl, imbe_spectral_rho(cur_mp->L));

    /* A P25 tone sent as voice frames (ProVoice has no DVSI reference): the
     * synthesis plays its tone model (imbe4400_play_tone()). */
    cur_mp->Vl[0] = (pitch == IMBE_PITCH_DVSI) ? imbe_tone_word(cur_mp, b0) : 0;

    return 0;
}

/**
 * @brief Decode IMBE 4400 parameters from demodulated bitstream.
 *
 * The pitch stored is a fit to the one DVSI's P25 decoder plays for b0
 * (imbe_decode_fundamental()): this decoder has no frame context, so a P25
 * one. The ProVoice frame paths keep the nominal pitch. A frame that is a P25
 * tone sent as voice frames leaves its tone in the otherwise unused Vl[0]
 * (imbe_tone_word()); every other frame it decodes, special and invalid frames
 * included, leaves 0 there. A status error (MBE_STATUS_INVALID_ARGUMENT,
 * MBE_STATUS_INVALID_BITS) is no frame and writes nothing. The word of prev_mp
 * is not written (the prediction writes its taps, as before). The decoded
 * model is unchanged.
 * @param imbe_d  Demodulated IMBE parameter bits (88).
 * @param cur_mp  Output: current frame parameters.
 * @param prev_mp Input: previous frame parameters (for prediction).
 * @return 0 on voice; non-zero for special frames (implementation-specific).
 */
int
mbe_decodeImbe4400Parms(const char* imbe_d, mbe_parms* cur_mp, mbe_parms* prev_mp) {
    return imbe4400_decode_parms(imbe_d, IMBE_PITCH_DVSI, cur_mp, prev_mp);
}

/**
 * @brief Demodulate interleaved IMBE 7200x4400 data in-place.
 * @param imbe Frame as 8x23 bitplanes (modified).
 */
int
mbe_demodulateImbe7200x4400Data(char imbe[8][23]) {

    int i, j, k;
    unsigned short pr[115];
    unsigned short foo;

    int ret = mbe_validate_bits((const char*)imbe, (size_t)8u * 23u);
    if (ret < 0) {
        return ret;
    }

    // create pseudo-random modulator
    foo = (unsigned short)mbe_bits_descending_to_int(imbe[0], 22, 11);
    pr[0] = (16 * foo);
    for (i = 1; i < 115; i++) {
        pr[i] = (173 * pr[i - 1]) + 13849 - (65536 * (((173 * pr[i - 1]) + 13849) / 65536));
    }
    for (i = 1; i < 115; i++) {
        pr[i] >>= 15; /* normalize to {0,1} cheaply */
    }

    // demodulate imbe with pr
    k = 1;
    for (i = 1; i < 4; i++) {
        for (j = 22; j >= 0; j--) {
            imbe[i][j] = imbe[i][j] ^ pr[k];
            k++;
        }
    }
    for (i = 4; i < 7; i++) {
        for (j = 14; j >= 0; j--) {
            imbe[i][j] = imbe[i][j] ^ pr[k];
            k++;
        }
    }
    return 0;
}

static void
mbe_demodulateImbe7200x4400DataSoft(mbe_soft_bit imbe[8][23]) {
    int i, j, k;
    unsigned short pr[115];
    unsigned short foo;

    foo = 0;
    for (i = 22; i >= 11; i--) {
        foo <<= 1;
        foo |= (unsigned short)(imbe[0][i].bit & 1u);
    }
    pr[0] = (unsigned short)(16 * foo);
    for (i = 1; i < 115; i++) {
        pr[i] = (unsigned short)((173 * pr[i - 1]) + 13849 - (65536 * (((173 * pr[i - 1]) + 13849) / 65536)));
    }
    for (i = 1; i < 115; i++) {
        pr[i] >>= 15;
    }

    k = 1;
    for (i = 1; i < 4; i++) {
        for (j = 22; j >= 0; j--) {
            imbe[i][j].bit = (uint8_t)((imbe[i][j].bit & 1u) ^ (uint8_t)pr[k]);
            k++;
        }
    }
    for (i = 4; i < 7; i++) {
        for (j = 14; j >= 0; j--) {
            imbe[i][j].bit = (uint8_t)((imbe[i][j].bit & 1u) ^ (uint8_t)pr[k]);
            k++;
        }
    }
}

int
mbe_decodeImbe7200x4400Frame(const char imbe_fr[8][23], char imbe_d[88], mbe_process_result* result) {
    char fr[8][23];
    int c0_errors;
    int protected_errors;
    int c4_errors = 0;
    int ret;

    if (result) {
        mbe_initProcessResult(result);
    }
    if (!imbe_d) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_validate_bits((const char*)imbe_fr, (size_t)8u * 23u);
    if (ret < 0) {
        return ret;
    }

    memcpy(fr, imbe_fr, sizeof(fr));
    c0_errors = mbe_eccImbe7200x4400C0(fr);
    ret = mbe_demodulateImbe7200x4400Data(fr);
    if (ret < 0) {
        return ret;
    }
    protected_errors = mbe_eccImbe7200x4400DataInternal(fr, imbe_d, &c4_errors);

    if (result) {
        result->c0_errors = c0_errors;
        result->protected_errors = protected_errors;
        result->c4_errors = c4_errors;
        result->total_errors = c0_errors + protected_errors;
        result->flags = MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_C4_VALID;
    }
    return c0_errors + protected_errors;
}

int
mbe_decodeImbe7200x4400SoftFrame(const mbe_soft_bit imbe_fr[8][23], char imbe_d[88], mbe_process_result* result) {
    mbe_soft_bit fr[8][23];
    int c0_errors;
    int protected_errors;
    int c4_errors = 0;
    int ret;

    if (result) {
        mbe_initProcessResult(result);
    }
    if (!imbe_d) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_validate_soft_bits((const mbe_soft_bit*)imbe_fr, (size_t)8u * 23u);
    if (ret < 0) {
        return ret;
    }

    memcpy(fr, imbe_fr, sizeof(fr));
    c0_errors = mbe_eccImbe7200x4400C0Soft(fr);
    mbe_demodulateImbe7200x4400DataSoft(fr);
    protected_errors = mbe_eccImbe7200x4400DataSoftInternal(fr, imbe_d, &c4_errors);

    if (result) {
        result->c0_errors = c0_errors;
        result->protected_errors = protected_errors;
        result->c4_errors = c4_errors;
        result->total_errors = c0_errors + protected_errors;
        result->flags = MBE_PROCESS_FLAG_SOFT_INPUT | MBE_PROCESS_FLAG_C0_VALID | MBE_PROCESS_FLAG_C4_VALID;
    }
    return c0_errors + protected_errors;
}

static int
imbe4400_prepare_process(mbe_process_result* result, const char imbe_d[88], mbe_parms* cur_mp, const mbe_parms* prev_mp,
                         const mbe_parms* prev_mp_enhanced, int* total_errors, int* c0_errors, int* c0_errors_valid) {
    int c4_errors_valid;
    int ret;

    if (!cur_mp || !prev_mp || !prev_mp_enhanced || !total_errors || !c0_errors || !c0_errors_valid) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_result_resolve_total_errors(result, total_errors);
    if (ret < 0) {
        return ret;
    }
    ret = mbe_validate_bits(imbe_d, 88u);
    if (ret < 0) {
        return ret;
    }

    *c0_errors_valid = (result && ((result->flags & MBE_PROCESS_FLAG_C0_VALID) != 0u));
    c4_errors_valid = (result && ((result->flags & MBE_PROCESS_FLAG_C4_VALID) != 0u));
    *c0_errors = *c0_errors_valid ? result->c0_errors : 0;
    cur_mp->errorCount4 = c4_errors_valid ? result->c4_errors : 0;
    mbe_result_prepare_synthesis(result, *total_errors);

    cur_mp->mutingThreshold = MBE_MUTING_THRESHOLD_IMBE;
    cur_mp->errorCountTotal = *total_errors;
    cur_mp->errorRate = (0.95f * prev_mp->errorRate) + (0.000365f * (float)*total_errors);
    return 0;
}

static int
imbe4400_repeat_required(int bad, int total_errors, int c0_errors, int c0_errors_valid, float repeat_threshold) {
    if (bad == 1) {
        return 1;
    }
    if (c0_errors_valid) {
        /* TIA-102.BABA 7.7 equations 97-98 repeat criteria when C0 errors are available from frame decode. */
        return ((c0_errors >= 2) && ((float)total_errors >= repeat_threshold));
    }

    /* Dataf callers pass parameter bits only (no C0 context); keep historical total-error fallback. */
    return (total_errors > 5);
}

static void
imbe4400_apply_repeat_decision(int repeat_required, mbe_process_result* result, mbe_parms* cur_mp,
                               const mbe_parms* prev_mp) {
    if (!repeat_required) {
        cur_mp->repeatCount = 0;
        return;
    }

    if (prev_mp->repeatCount > (MBE_MAX_FRAME_REPEATS - 1)) {
        /* Library repeat-headroom policy: reset to default model after prolonged repeats. */
        mbe_resetImbeModel(cur_mp);
    } else {
        /* Repeat the previous model; this frame's error accounting and noise state continue. */
        mbe_repeat_load_model(cur_mp, prev_mp);
        cur_mp->repeatCount = mbe_repeat_next_count(prev_mp->repeatCount);
    }
    mbe_result_set_flag(result, MBE_PROCESS_FLAG_REPEAT);
}

/* Whether the parameters come from a ProVoice (IMBE 7100x4400) frame. */
static int
imbe4400_is_provoice(const mbe_process_result* result) {
    return result && ((result->flags & MBE_PROCESS_FLAG_PROVOICE) != 0u);
}

/* P25 mutes with the TIA-102.BABA 7.8 noise; ProVoice, not a TIA-102 codec, keeps library comfort noise. */
static enum mbe_mute_noise
imbe4400_mute_noise(const mbe_process_result* result) {
    return imbe4400_is_provoice(result) ? MBE_MUTE_NOISE_COMFORT : MBE_MUTE_NOISE_SPEC;
}

/* The tone model's level: the tone harmonics' rms (decoded magnitudes) times
 * this, then times DVSI's synthesis table droop (imbe4400_play_tone()). DVSI's
 * dual tones measure a median 2.0974 times that rms on the 3,574 steady DTMF
 * and KNOX frames (2.0975 with the 143 call-progress frames the tone path
 * leaves out) and its single tones 2.114, where this library renders a voiced
 * harmonic at 2.1308 times its magnitude; at
 * 0.988 the steady tone frames of knox_1..16 and dtone_1..16 play a median
 * 0.012 dB above DVSI's (0.9866 would be exact there) and those of the DTMF,
 * alert, alltone and call-progress vectors 0.022 dB (dual) and 0.031 dB
 * (single) below it. */
#define IMBE_TONE_LEVEL 0.988

/* The tone a frame's model carries (the word imbe_tone_word() left in Vl[0])
 * and its level, or 0; read before enhancement, on the decoded model (or the
 * repeated one, the decoded history). A word counts only if the decode's own
 * test, run again on the model it sits in, gives exactly that word: the model
 * at the pitch its code decodes to (w0 by its bits, L and K:
 * imbe_code_pitch()), and imbe_tone_word() on the model's voicing and decoded
 * magnitudes equal to the word. A decode's own word passes on its frame and on
 * every repeat of it; a caller's state can hold any other word, and it plays
 * as no tone. */
struct imbe_tone {
    int b0;
    int l1;
    int l2;
    double amplitude; /* the rms of the tone harmonics' decoded magnitudes, times IMBE_TONE_LEVEL */
};

/* Whether the model is at the pitch, harmonic count and band count P25 pitch
 * code b0 decodes to (the pitch by its bits). */
static int
imbe_model_at_code(const mbe_parms* mp, int b0) {
    float w0;
    int L = 0, K = 0;
    if (imbe_code_pitch(b0, IMBE_PITCH_DVSI, &w0, &L, &K) != 0 || mp->L != L || mp->K != K) {
        return 0;
    }
    uint32_t want, have;
    memcpy(&want, &w0, sizeof(want));
    memcpy(&have, &mp->w0, sizeof(have));
    return want == have;
}

static int
imbe4400_take_tone(const mbe_parms* mp, struct imbe_tone* tone) {
    const int word = mp->Vl[0];
    if ((word & ~0x00FFFFFF) != IMBE_TONE_FLAG) {
        return 0;
    }
    const int b0 = (word >> 16) & 0xFF;
    /* A code without rows returns after one lookup, before any pitch. */
    if (b0 > IMBE_TONE_LAST_CODE || imbe_tone_first_row[b0] == imbe_tone_first_row[b0 + 1]
        || !imbe_model_at_code(mp, b0) || imbe_tone_word(mp, b0) != word) {
        return 0;
    }
    /* imbe_tone_word() gave this word: l1 and l2 are within L, their
     * magnitudes positive, finite and normal. */
    tone->b0 = b0;
    tone->l2 = (word >> 8) & 0xFF;
    tone->l1 = word & 0xFF;
    const double p1 = (double)mp->Ml[tone->l1] * (double)mp->Ml[tone->l1];
    const double p2 = (double)mp->Ml[tone->l2] * (double)mp->Ml[tone->l2];
    tone->amplitude = IMBE_TONE_LEVEL * sqrt((tone->l1 == tone->l2) ? p1 : 0.5 * (p1 + p2));
    return 1;
}

/* sinc^2(pi l / 256), sinc(x) = sin(x) / x: the droop of DVSI's synthesis
 * table (US 8,200,497's K = 256 table read with linear interpolation) at
 * harmonic l. */
static double
imbe_tone_droop(int l) {
    const double x = M_PI * (double)l / 256.0;
    const double sinc = sin(x) / x;
    return sinc * sinc;
}

/*
 * The model DVSI's decoder plays for a P25 tone sent as voice frames (measured
 * on its decodes of the steady frames of the DTMF and KNOX tones, 3,574, and of
 * the single tones): its pitch is eq 46's fundamental cut to a grid of 2^-19 cycle per
 * sample, floor(2^21 / (2 b0 + 79)) / 2^19, which is the pitch of its code
 * (imbe_grid_fundamental(), for speech as for tones; imbe4400_take_tone() has
 * found the model at that pitch by its bits, so it is kept); the tone's harmonics are
 * voiced at one amplitude, the rms of their decoded magnitudes times
 * IMBE_TONE_LEVEL, times sinc^2(pi l / 256) (on DVSI's P25 decodes the upper
 * component over the lower is a median 0.9999 of the droop's ratio; the bound
 * of 0.0004 is the rate-33 and D-STAR measurement, docs/dvsi-parity.md); and every other
 * harmonic is unvoiced at 0, the library's choice, as TIA-102.BABA-1 clause 8
 * builds a tone's voice parameters (the measurement bounds them: two components
 * leave 0.0009 of the frame, so the others are at least 27 dB down; the patent
 * says 20 dB is typical). L, K and log2Ml stay the decoded frame's: the voiced
 * phases are regenerated from them (on DVSI's tone frames the frame-to-frame
 * phase error's standard deviation has a median of 0.0009 cycle over 64
 * components, around mean offsets of up to 0.07 cycle). Only this function marks a frame MBE_FRAME_MARK_TONE
 * (imbe4400_synthesize_frame() clears the mark on any other), and the mark
 * makes the synthesizer take the frame's own model at mid-frame
 * (with its own phases: on the steady frames of DVSI's knox_1..16 and
 * dtone_1..16 they follow DVSI's waveform at a median correlation of 0.9924,
 * those from the mean of the two frames' log2Ml at 0.9903). Unlike an AMBE
 * tone frame, such a frame runs DVSI's noise generator on as any voice frame
 * does: the unvoiced frames two and three after a tone run of dtmf and
 * alltone follow DVSI's waveform at a median correlation of 0.80 and 0.99
 * (154 and 7 frames) so, and at -0.01 and -0.10 with the generator held.
 */
static void
imbe4400_play_tone(mbe_parms* mp, const struct imbe_tone* tone) {
    for (int l = 1; l <= mp->L; l++) {
        const int sounds = (l == tone->l1) || (l == tone->l2);
        mp->Vl[l] = sounds ? 1 : 0;
        mp->Ml[l] = sounds ? (float)(tone->amplitude * imbe_tone_droop(l)) : 0.0f;
    }
    mbe_set_frame_mark(mp, MBE_FRAME_MARK_TONE);
}

/* Whether a model has a voiced band; none for a harmonic count outside 1..56. */
static int
imbe_has_voiced_band(const mbe_parms* mp) {
    if (!mbe_harmonic_count_is_valid(mp->L)) {
        return 0;
    }
    for (int l = 1; l <= mp->L; l++) {
        if (mp->Vl[l] == 1) {
            return 1;
        }
    }
    return 0;
}

static void
imbe4400_synthesize_frame(float* aout_buf, mbe_process_result* result, enum mbe_mute_noise mute_noise,
                          mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    int frame_muted = mbe_isMaxFrameRepeat(cur_mp) || mbe_requiresMuting(cur_mp);
    struct imbe_tone tone;

    /* Only imbe4400_play_tone() marks a frame MBE_FRAME_MARK_TONE, on the
     * frame it plays: the decoders never do, and the history never holds the
     * mark, so one here came from a caller's state and is not the frame's. It
     * gives way to the mark a decode leaves a model with a voiced band
     * (MBE_FRAME_MARK_ONSET_COLUMNS, imbe_decode_voicing()), or to none; the
     * model rendered is marked again after smoothing, by its own voicing
     * (MBE_VOICED_SYNTHESIS_DVSI_IMBE, mbe_synthesize_subframes()). A decoded
     * mark the caller's state overwrote (MBE_FRAME_MARK_PITCH_ABSENT on a frame
     * with no voiced band) is not recovered.
     * ProVoice keeps no tone word either: Vl[0] reads 0 after its frames, so
     * a P25 repeat after a ProVoice frame (a state moved between the two
     * entry points, outside the stream contract) plays the repeated model as
     * speech. */
    if (mbe_frame_mark(cur_mp) == MBE_FRAME_MARK_TONE) {
        mbe_set_frame_mark(cur_mp, imbe_has_voiced_band(cur_mp) ? MBE_FRAME_MARK_ONSET_COLUMNS : MBE_FRAME_MARK_NONE);
    }
    if (imbe4400_is_provoice(result)) {
        cur_mp->Vl[0] = 0;
    }
    mbe_moveMbeParms(cur_mp, prev_mp);
    /* The prediction history (prev_mp) keeps the decoded frame; a tone frame's
     * synthesis model is the tone model, from its decoded magnitudes. */
    const int plays_tone = !imbe4400_is_provoice(result) && imbe4400_take_tone(cur_mp, &tone);
    float pre_enh_rm0 = mbe_spectralAmpEnhanceWithRm0(cur_mp, MBE_HIGH_BAND_GAIN_IMBE);
    if (plays_tone) {
        imbe4400_play_tone(cur_mp, &tone);
    }
    /* As DVSI's P25 decoder, without eq 116's scaling; ProVoice, which plays
     * through the same path, likewise, in P25's rendering without the pulsed
     * state (no DVSI ProVoice material shows one). */
    const enum mbe_voiced_synthesis rendering =
        imbe4400_is_provoice(result) ? MBE_VOICED_SYNTHESIS_DVSI_PROVOICE : MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    mbe_synthesizeSpeechWithPreEnhRm0f(aout_buf, cur_mp, prev_mp_enhanced, pre_enh_rm0, mute_noise, rendering,
                                       MBE_AMPLITUDE_LIMIT_NONE);

    if (frame_muted) {
        mbe_result_set_flag(result, MBE_PROCESS_FLAG_MUTE);
    }

    mbe_moveMbeParms(cur_mp, prev_mp_enhanced);
}

int
mbe_processImbe4400Dataf_internal(float* aout_buf, mbe_process_result* result, const char imbe_d[88], mbe_parms* cur_mp,
                                  mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    int bad;
    int repeat_required;
    float repeat_threshold;
    int total_errors;
    int c0_errors;
    int c0_errors_valid;
    int ret;

    if (!aout_buf || !cur_mp || !prev_mp || !prev_mp_enhanced) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = imbe4400_prepare_process(result, imbe_d, cur_mp, prev_mp, prev_mp_enhanced, &total_errors, &c0_errors,
                                   &c0_errors_valid);
    if (ret < 0) {
        return ret;
    }

    /* P25 plays DVSI's pitch for b0; ProVoice, with no DVSI reference, the nominal one. */
    bad = imbe4400_decode_parms(imbe_d, imbe4400_is_provoice(result) ? IMBE_PITCH_NOMINAL : IMBE_PITCH_DVSI, cur_mp,
                                prev_mp);
    if (bad < 0) {
        return bad;
    }
    repeat_threshold = 10.0f + (40.0f * cur_mp->errorRate);
    repeat_required = imbe4400_repeat_required(bad, total_errors, c0_errors, c0_errors_valid, repeat_threshold);

    imbe4400_apply_repeat_decision(repeat_required, result, cur_mp, prev_mp);
    imbe4400_synthesize_frame(aout_buf, result, imbe4400_mute_noise(result), cur_mp, prev_mp, prev_mp_enhanced);
    return result ? result->total_errors : total_errors;
}

/**
 * @brief Process IMBE 4400 parameters into 160 float samples at 8 kHz.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional in/out status context.
 * @param imbe_d   Demodulated parameter bits (88).
 * @param cur_mp   In/out: current frame parameters (may be enhanced).
 * @param prev_mp  In/out: previous frame parameters.
 * @param prev_mp_enhanced In/out: enhanced previous parameters for continuity.
 */
int
mbe_processImbe4400Dataf(float* aout_buf, mbe_process_result* result, const char imbe_d[88], mbe_parms* cur_mp,
                         mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;

    if (!result) {
        mbe_initProcessResult(&local_result);
        result = &local_result;
    }
    return mbe_processImbe4400Dataf_internal(aout_buf, result, imbe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processImbe4400Data(short* aout_buf, mbe_process_result* result, const char imbe_d[88], mbe_parms* cur_mp,
                        mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processImbe4400Dataf(float_buf, result, imbe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete IMBE 7200x4400 frame into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional output status populated by decode and synthesis.
 * @param imbe_fr  Input frame as 8x23 bitplanes.
 * @param imbe_d   Scratch/output parameter bits (88).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 */
int
mbe_processImbe7200x4400Framef(float* aout_buf, mbe_process_result* result, const char imbe_fr[8][23], char imbe_d[88],
                               mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;
    int ret;
    if (!result) {
        result = &local_result;
    }
    ret = mbe_decodeImbe7200x4400Frame(imbe_fr, imbe_d, result);
    if (ret < 0) {
        return ret;
    }
    return mbe_processImbe4400Dataf(aout_buf, result, imbe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processImbe7200x4400SoftFramef(float* aout_buf, mbe_process_result* result, const mbe_soft_bit imbe_fr[8][23],
                                   char imbe_d[88], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                   mbe_parms* prev_mp_enhanced) {
    mbe_process_result local_result;
    int ret;
    if (!result) {
        result = &local_result;
    }
    ret = mbe_decodeImbe7200x4400SoftFrame(imbe_fr, imbe_d, result);
    if (ret < 0) {
        return ret;
    }
    return mbe_processImbe4400Dataf(aout_buf, result, imbe_d, cur_mp, prev_mp, prev_mp_enhanced);
}

int
mbe_processImbe7200x4400SoftFrame(short* aout_buf, mbe_process_result* result, const mbe_soft_bit imbe_fr[8][23],
                                  char imbe_d[88], mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    float float_buf[160];
    int ret;
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processImbe7200x4400SoftFramef(float_buf, result, imbe_fr, imbe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}

/**
 * @brief Process a complete IMBE 7200x4400 frame into 16-bit PCM.
 * @see mbe_processImbe7200x4400Framef for details.
 */
int
mbe_processImbe7200x4400Frame(short* aout_buf, mbe_process_result* result, const char imbe_fr[8][23], char imbe_d[88],
                              mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {

    float float_buf[160];
    int ret;
    if (!aout_buf) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    ret = mbe_processImbe7200x4400Framef(float_buf, result, imbe_fr, imbe_d, cur_mp, prev_mp, prev_mp_enhanced);
    if (ret < 0) {
        return ret;
    }
    mbe_floattoshort(float_buf, aout_buf);
    return ret;
}
