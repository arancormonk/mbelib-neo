// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/*
 * The AMBE log-magnitude prediction residual block lengths (D-STAR's
 * AmbePlusLmprbl, rate 33's AmbeLmprbl): TIA-102.BABA-1 Annex C but for L =
 * 20, which is DVSI's codec's 4, 5, 5, 6 where the annex prints 4, 4, 6, 6
 * (docs/dvsi-parity.md, "AMBE block lengths at 20 harmonics").
 *
 * - The tables: each row of L = 9..56 sums to L, the two decoders' tables and
 *   the two encoders' (mbe_ambe2400_enc_tables(), mbe_ambe2450_enc_tables())
 *   are equal, row 20 is 4, 5, 5, 6 and the other 47 rows are Annex C's as
 *   printed (copied below).
 * - The decoders against a block transform formed here in double: two frames
 *   built from chosen codes, at 19, 20 and 21 harmonics (every code of 20
 *   harmonics), at several V/UV rows (D-STAR's b1; rate 33's rows 0 to 31,
 *   the band-grid rows 18 to 31 among them, which only the decoder meets: the
 *   encoder never sends them), through the public parameter decoders; the
 *   second frame is predicted from the first. The chosen codes keep every
 *   other division of the frame's harmonics into blocks that moves one
 *   harmonic between two neighbouring blocks (and 5, 5, 5, 5 at L = 20) more
 *   than MARGIN away from the decode, so each such row in a table fails.
 * - The decoders against an exact result: the same frames formed here in
 *   float by the decoder's arithmetic written out (exact_frame()), compared by
 *   bit pattern where the build allows it and within EXACT_BOUND otherwise.
 * - The encoders: decoder and encoder state parity on a 20-harmonic signal;
 *   the joint refinement in effect at L = 20 (the shared quantizer and rate
 *   33's codec path); and a 20-harmonic residual that blocks of 4, 5, 5, 6
 *   represent exactly round-trips through the quantizer with the encoders'
 *   tables and not with the annex's row.
 */
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "ambe3600x2400_const.h"
#include "ambe3600x2400_internal.h"
#include "ambe3600x2450_const.h"
#include "ambe3600x2450_internal.h"
#include "ambe_encoder.h"
#include "mbelib-neo/mbelib.h"

/* TIA-102.BABA-1 v1.0.5 Annex C, J1..J4 for L = 9..56, as printed. The one
 * row the library departs from is L = 20 (4, 4, 6, 6 here). */
static const int annex_c[48][4] = {
    {2, 2, 2, 3},     {2, 2, 3, 3},     {2, 3, 3, 3},     {2, 3, 3, 4},     {3, 3, 3, 4},     {3, 3, 4, 4},
    {3, 3, 4, 5},     {3, 4, 4, 5},     {3, 4, 5, 5},     {4, 4, 5, 5},     {4, 4, 5, 6},     {4, 4, 6, 6},
    {4, 5, 6, 6},     {4, 5, 6, 7},     {5, 5, 6, 7},     {5, 5, 7, 7},     {5, 6, 7, 7},     {5, 6, 7, 8},
    {5, 6, 8, 8},     {6, 6, 8, 8},     {6, 6, 8, 9},     {6, 7, 8, 9},     {6, 7, 9, 9},     {6, 7, 9, 10},
    {7, 7, 9, 10},    {7, 8, 9, 10},    {7, 8, 10, 10},   {7, 8, 10, 11},   {8, 8, 10, 11},   {8, 9, 10, 11},
    {8, 9, 11, 11},   {8, 9, 11, 12},   {8, 9, 11, 13},   {8, 9, 12, 13},   {8, 10, 12, 13},  {9, 10, 12, 13},
    {9, 10, 12, 14},  {9, 10, 13, 14},  {9, 11, 13, 14},  {10, 11, 13, 14}, {10, 11, 13, 15}, {10, 11, 14, 15},
    {10, 12, 14, 15}, {10, 12, 14, 16}, {11, 12, 14, 16}, {11, 12, 15, 16}, {11, 12, 15, 17}, {11, 13, 15, 17},
};
static const int dvsi_row_20[4] = {4, 5, 5, 6};

/* The block lengths the library is to use at L: the annex's, but DVSI's row at 20. */
static void
expected_row(int L, int J[4]) {
    for (int i = 0; i < 4; ++i) {
        J[i] = (L == 20) ? dvsi_row_20[i] : annex_c[L - 9][i];
    }
}

static int
rows_equal(const int a[4], const int b[4]) {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

static int
test_tables(void) {
    const int (*tables[4])[4] = {AmbePlusLmprbl, AmbeLmprbl, mbe_ambe2400_enc_tables()->lmprbl,
                                 mbe_ambe2450_enc_tables()->lmprbl};
    static const char* const names[4] = {"D-STAR decoder", "rate-33 decoder", "D-STAR encoder", "rate-33 encoder"};
    int fails = 0;
    for (int t = 0; t < 4; ++t) {
        for (int L = 0; L <= 56; ++L) {
            const int* row = tables[t][L];
            if (L < 9) {
                static const int zero[4] = {0, 0, 0, 0};
                fails += !rows_equal(row, zero);
                continue;
            }
            int J[4];
            expected_row(L, J);
            if (row[0] + row[1] + row[2] + row[3] != L || !rows_equal(row, J)) {
                printf("%s table, L %d: %d %d %d %d\n", names[t], L, row[0], row[1], row[2], row[3]);
                fails++;
            }
        }
    }
    /* The annex's own rows sum to L (the copy above is checked against itself). */
    for (int L = 9; L <= 56; ++L) {
        const int* row = annex_c[L - 9];
        fails += (row[0] + row[1] + row[2] + row[3]) != L;
    }
    printf("block-length tables: %s\n", fails ? "FAIL" : "the annex's 47 rows and DVSI's row 20, in all four");
    return fails;
}

/* ---- The decoders against a block transform formed here. ---- */

struct codes {
    int b[9]; /* b0..b8 (D-STAR's b8 an even row) */
};

/* D-STAR's 49 parameter bits from b0..b8, at the positions its decoder reads them. */
static void
dstar_pack(const int b[9], char d[49]) {
    static const int b0_pos[7] = {0, 1, 2, 3, 4, 5, 48};
    static const int b1_pos[4] = {38, 39, 40, 41};
    static const int b2_pos[6] = {6, 7, 8, 9, 42, 43};
    static const int b3_pos[9] = {10, 11, 12, 13, 14, 15, 16, 44, 45};
    static const int b4_pos[7] = {17, 18, 19, 20, 21, 46, 47};
    static const int b5_pos[4] = {22, 23, 25, 26};
    static const int b6_pos[4] = {27, 28, 29, 30};
    static const int b7_pos[4] = {31, 32, 33, 34};
    static const int b8_pos[3] = {35, 36, 37};
    memset(d, 0, 49);
    const int* pos[9] = {b0_pos, b1_pos, b2_pos, b3_pos, b4_pos, b5_pos, b6_pos, b7_pos, b8_pos};
    static const int width[9] = {7, 4, 6, 9, 7, 4, 4, 4, 3};
    for (int f = 0; f < 9; ++f) {
        const int value = (f == 8) ? (b[8] >> 1) : b[f];
        for (int k = 0; k < width[f]; ++k) {
            d[pos[f][k]] = (char)((value >> (width[f] - 1 - k)) & 1);
        }
    }
}

/* The history a decode predicts from, as the test keeps it. */
struct history {
    int L;
    double log2Ml[57];
    double gamma;
};

/*
 * log2Ml of a frame of L harmonics with block lengths J from the codes b and
 * the history h, formed in double from the codebooks (TIA-102.BABA-1 4.4 eqs
 * 27-45 for rate 33; D-STAR's decoder by the same steps with its own tables,
 * gain law and prediction coefficient):
 *   G from the PRBA24 row b3 and the PRBA58 row b4, R by the 8-point inverse
 *   DCT; each block's first two coefficients from R, the next from its HOC row
 *   (none past the sixth); T by each block's inverse DCT;
 *   the prediction at the position the decoder rounds: the ratio prev_L / L
 *   rounded to float, times l, rounded to float again (mbe_prediction_position()'s
 *   arithmetic, written out here), its whole part the lower tap and the rest
 *   the weight, the history extended past its L by its last value and
 *   log2Ml[0] taken as log2Ml[1];
 *   log2Ml = T - mean(T) + rho (P - mean(P)) + gamma - 0.5 log2 L, with gamma
 *   rate 33's Dg[b2] + 0.5 gamma(-1) and D-STAR's 2 Dg[b2] - 0.09.
 * Everything else is double; the decoder works in float, so the comparison
 * allows TOLERANCE.
 */
static void
reference_frame(int dstar, const int b[9], int L, const int J[4], const struct history* h, struct history* out) {
    const float (*prba24)[3] = dstar ? AmbePlusPRBA24 : AmbePRBA24;
    const float (*prba58)[4] = dstar ? AmbePlusPRBA58 : AmbePRBA58;
    const float (*hoc[4])[4] = {dstar ? AmbePlusHOCb5 : AmbeHOCb5, dstar ? AmbePlusHOCb6 : AmbeHOCb6,
                                dstar ? AmbePlusHOCb7 : AmbeHOCb7, dstar ? AmbePlusHOCb8 : AmbeHOCb8};
    double G[9] = {0}, R[9], T[57] = {0};
    for (int m = 0; m < 3; ++m) {
        G[2 + m] = prba24[b[3]][m];
    }
    for (int m = 0; m < 4; ++m) {
        G[5 + m] = prba58[b[4]][m];
    }
    for (int i = 1; i <= 8; ++i) {
        R[i] = 0.0;
        for (int m = 1; m <= 8; ++m) {
            R[i] += (m == 1 ? 1.0 : 2.0) * G[m] * cos(M_PI * (m - 1) * (i - 0.5) / 8.0);
        }
    }
    int l = 1;
    for (int i = 1; i <= 4; ++i) {
        double C[18] = {0};
        const int Ji = J[i - 1];
        const int odd = i + i - 1; /* R[2i - 1] and R[2i] */
        C[1] = 0.5 * (R[odd] + R[odd + 1]);
        C[2] = (R[odd] - R[odd + 1]) / (2.0 * sqrt(2.0));
        for (int k = 3; k <= Ji && k <= 6; ++k) {
            C[k] = hoc[i - 1][b[4 + i]][k - 3];
        }
        for (int j = 1; j <= Ji; ++j, ++l) {
            for (int k = 1; k <= Ji; ++k) {
                T[l] += (k == 1 ? 1.0 : 2.0) * C[k] * cos(M_PI * (k - 1) * (j - 0.5) / Ji);
            }
        }
    }
    double prev[58];
    for (int k = 1; k <= 57; ++k) {
        prev[k] = h->log2Ml[(k <= h->L) ? k : h->L];
    }
    prev[0] = prev[1];
    double P[57] = {0}, meanT = 0.0, meanP = 0.0;
    for (int k = 1; k <= L; ++k) {
        volatile float ratio = (float)h->L / (float)L;
        volatile float position = ratio * (float)k;
        const int lower = (int)position;
        const double delta = (double)position - lower;
        const int upper = (lower + 1 > 56) ? 56 : lower + 1;
        P[k] = ((1.0 - delta) * prev[lower]) + (delta * prev[upper]);
        meanT += T[k] / L;
        meanP += P[k] / L;
    }
    const double rho = dstar ? (double)MBE_AMBE2400_PREDICTION_RHO : 0.65;
    out->gamma = dstar ? ((double)MBE_AMBE2400_GAIN_SCALE * AmbePlusDg[b[2]]) + (double)MBE_AMBE2400_GAIN_OFFSET
                       : (double)AmbeDg[b[2]] + (0.5 * h->gamma);
    out->L = L;
    for (int k = 0; k <= 56; ++k) {
        out->log2Ml[k] = 0.0;
    }
    for (int k = 1; k <= L; ++k) {
        out->log2Ml[k] = T[k] - meanT + (rho * (P[k] - meanP)) + out->gamma - (0.5 * log2((double)L));
    }
}

/* A float decode against the double reference. Values reach about 10 here, where a float's ulp is 9.5e-7, so the
 * bound is about 210 ulps there and more below: the decoder's single-precision sums (the 8-point DCT, each block's
 * inverse DCT, the means over 20 harmonics, the prediction) each round, and the reference rounds nothing but the
 * prediction's position. It catches a wrong row by three orders of magnitude (MARGIN) and nothing below itself;
 * exact_frame() pins the arithmetic. */
#define TOLERANCE 2e-4
#define MARGIN    0.05 /* every other division tried lies further than this from the decode */

/* Finite: at most DBL_MAX in magnitude, judged on the bit pattern so that no build's float semantics can pass a
 * NaN or an infinity. */
static int
finite_double(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return (bits & UINT64_C(0x7FFFFFFFFFFFFFFF)) < UINT64_C(0x7FF0000000000000);
}

/* The largest difference of log2Ml over the harmonics, or -1 if any difference is not finite (every caller fails
 * on a negative result, so a NaN cannot pass as a small difference nor as a large margin). */
static double
max_difference(const struct history* a, const struct history* b) {
    double worst = 0.0;
    for (int k = 1; k <= a->L; ++k) {
        const double d = fabs(a->log2Ml[k] - b->log2Ml[k]);
        if (!finite_double(d)) {
            return -1.0;
        }
        worst = fmax(worst, d);
    }
    return worst;
}

/* The divisions of L into four blocks that move one harmonic between two neighbouring blocks of J (each block
 * keeping at least one), and 5, 5, 5, 5 at L 20: n of them in alt. */
static int
alternative_rows(int L, const int J[4], int alt[8][4]) {
    int n = 0;
    for (int i = 0; i < 3; ++i) {
        for (int dir = -1; dir <= 1; dir += 2) {
            int row[4] = {J[0], J[1], J[2], J[3]};
            row[i] += dir;
            row[i + 1] -= dir;
            if (row[i] >= 1 && row[i + 1] >= 1) {
                memcpy(alt[n++], row, sizeof(row));
            }
        }
    }
    if (L == 20) {
        static const int flat[4] = {5, 5, 5, 5};
        memcpy(alt[n++], flat, sizeof(flat));
    }
    return n;
}

static void
history_from(const mbe_parms* mp, struct history* h) {
    h->L = mp->L;
    h->gamma = mp->gamma;
    for (int k = 0; k <= 56; ++k) {
        h->log2Ml[k] = mp->log2Ml[k];
    }
}

/* Decodes the frame b at L harmonics after `state`, checks it against the reference with the expected row and
 * against every alternative, and leaves the decoded state as the next frame's history. Returns failures. */
static int
check_frame(int dstar, const int b[9], int L, mbe_parms* cur, mbe_parms* prev, struct history* ref_history,
            double* least_margin) {
    char d[49];
    if (dstar) {
        dstar_pack(b, d);
    } else {
        mbe_ambe2450_pack(b, d);
    }
    const int rc = dstar ? mbe_decodeAmbe2400Parms(d, cur, prev) : mbe_decodeAmbe2450Parms(d, cur, prev);
    if (rc != 0 || cur->L != L) {
        printf("%s: b0 %d decoded rc %d L %d, wanted a voice frame of %d\n", dstar ? "D-STAR" : "rate 33", b[0], rc,
               cur->L, L);
        return 1;
    }
    struct history decoded, expected, other;
    history_from(cur, &decoded);
    int J[4], alt[8][4];
    expected_row(L, J);
    reference_frame(dstar, b, L, J, ref_history, &expected);
    int fails = 0;
    const double error = max_difference(&decoded, &expected);
    if (!(error >= 0.0 && error <= TOLERANCE)) {
        printf("%s L %d: decode differs from rows %d %d %d %d by %.3g\n", dstar ? "D-STAR" : "rate 33", L, J[0], J[1],
               J[2], J[3], error);
        fails++;
    }
    const int n = alternative_rows(L, J, alt);
    for (int a = 0; a < n; ++a) {
        reference_frame(dstar, b, L, alt[a], ref_history, &other);
        const double margin = max_difference(&expected, &other);
        if (margin >= 0.0) {
            *least_margin = fmin(*least_margin, margin);
        }
        if (!(margin > MARGIN)) {
            printf("%s L %d: rows %d %d %d %d lie only %.3g from the expected decode\n", dstar ? "D-STAR" : "rate 33",
                   L, alt[a][0], alt[a][1], alt[a][2], alt[a][3], margin);
            fails++;
        }
    }
    *ref_history = expected;
    mbe_moveMbeParms(cur, prev);
    return fails;
}

/* Codes of 19, 20 and 21 harmonics, found by decoding every code (each check asserts the decoded L): D-STAR b0
 * 51..53 (19), 54..57 (20), 58..60 (21); rate 33 b0 48..50, 51..53, 54..56. */
static int
code_of(int dstar, int L, int k) {
    static const int first[2][3] = {{48, 51, 54}, {51, 54, 58}};
    static const int count[2][3] = {{3, 3, 3}, {3, 4, 3}};
    const int c = dstar ? 1 : 0;
    return first[c][L - 19] + (k % count[c][L - 19]);
}

/* The V/UV fields the frames are built at: D-STAR's b1 (four 1 kHz bands), rate 33's rows, the band-grid rows 18
 * to 31 (w0 = 2.81 / L at the code's L) among them. */
static const int dstar_vuv[4] = {15, 0, 9, 3};
static const int r33_vuv[12] = {0, 18, 24, 31, 5, 16, 17, 13, 21, 26, 29, 3};

/* A history with a shape, so the prediction weighs: float in the decoder's state, its value in the reference. */
static void
shaped_history(mbe_parms* prev, struct history* h) {
    h->L = 20;
    h->gamma = 0.0;
    prev->L = 20;
    prev->gamma = 0.0f;
    for (int k = 0; k <= 56; ++k) {
        prev->log2Ml[k] = (float)((k == 0) ? 0.0 : 2.0 + sin(0.7 * k));
        h->log2Ml[k] = (double)prev->log2Ml[k];
    }
}

static int
test_decoder(int dstar) {
    /* Two frames, b2..b8 chosen once per codec (receivable rows: D-STAR's b8 even); b0 and b1 set per case. */
    static const int dstar_rows[2][9] = {{0, 15, 40, 301, 72, 3, 11, 9, 6}, {0, 15, 44, 117, 21, 12, 5, 14, 2}};
    static const int r33_rows[2][9] = {{0, 0, 20, 301, 72, 19, 11, 9, 5}, {0, 0, 22, 117, 21, 6, 5, 14, 2}};
    const int* vuv = dstar ? dstar_vuv : r33_vuv;
    const int nvuv = dstar ? 4 : 12;
    int fails = 0, cases = 0;
    double least = DBL_MAX;
    for (int first = 19; first <= 21; ++first) {
        for (int second = 19; second <= 21; ++second) {
            for (int v = 0; v < nvuv; ++v) {
                mbe_parms cur, prev, enhanced;
                struct history h;
                int b[9];
                mbe_initMbeParms(&cur, &prev, &enhanced);
                shaped_history(&prev, &h);
                memcpy(b, dstar ? dstar_rows[0] : r33_rows[0], sizeof(b));
                b[0] = code_of(dstar, first, v);
                b[1] = vuv[v];
                fails += check_frame(dstar, b, first, &cur, &prev, &h, &least);
                memcpy(b, dstar ? dstar_rows[1] : r33_rows[1], sizeof(b));
                b[0] = code_of(dstar, second, v + 1);
                b[1] = vuv[(v + 5) % nvuv];
                fails += check_frame(dstar, b, second, &cur, &prev, &h, &least);
                cases++;
            }
        }
    }
    printf("%s decoder: %d pairs of frames at 19, 20 and 21 harmonics and %d V/UV fields follow the block lengths "
           "(least margin to another row %.3f)\n",
           dstar ? "D-STAR" : "rate-33", cases, nvuv, least);
    return fails;
}

/* ---- The decoders against an exact result. ---- */

/* The decoder's state as the exact result keeps it: float, as the decoder holds it. */
struct fhistory {
    int L;
    float log2Ml[57];
    float gamma;
};

/* The decoders' cosine tables' entries: cosf() of the argument formed as the tables form it (M_PI in double, the
 * product in double, rounded to float on the call). */
static float
dct8_cos(int m, int i) {
    return cosf((M_PI * (float)(m - 1) * ((float)i - 0.5f)) / 8.0f);
}

static float
block_cos(int ji, int j, int k) {
    return cosf((M_PI * (float)(k - 1) * ((float)j - 0.5f)) / (float)ji);
}

/*
 * log2Ml and gamma of a frame of L harmonics with block lengths J from the codes b and the history h, formed in
 * float by the arithmetic of TIA-102.BABA-1 4.4 (rate 33; D-STAR by the same steps with its own tables, its absolute
 * gain law and its prediction coefficient) in the order the decoders evaluate it, every operation in float:
 *   gamma = Dg[b2] + 0.5 gamma(-1) (D-STAR: 2 Dg[b2] - 0.09);
 *   R_i = sum over m = 1..8, in that order, of a_m G_m cos(pi (m - 1)(i - 0.5) / 8), a_1 = 1, a_m = 2;
 *   C_i,1 = 0.5 (R_2i-1 + R_2i), C_i,2 = r (R_2i-1 - R_2i) with r = 1 / (2 sqrt 2) rounded to float, C_i,3..6 the
 *   HOC row (none past the block's length or the sixth);
 *   T_l = sum over k = 1..J_i, in that order, of a_k C_i,k cos(pi (k - 1)(j - 0.5) / J_i);
 *   the history extended past its L by its last value, log2Ml[0] = log2Ml[1]; position = (prev_L / L) l, the ratio
 *   and the product each rounded to float; P_l = (1 - d) h[k] + d h[k + 1];
 *   Sum43 = (rho / L) sum P_l (rho 0.65; D-STAR MBE_AMBE2400_PREDICTION_RHO, 0.8), Sum42 = (sum T_l) / L,
 *   Gamma = gamma - 0.5 log2f(L) - Sum42, log2Ml_l = T_l + rho (1 - d) h[k] + rho d h[k + 1] - Sum43 + Gamma, left
 *   to right.
 * The library pieces it reads are the codebooks (the PRBA, HOC and gain tables of the codecs' headers) and D-STAR's
 * three constants; it calls no library function but cosf() and log2f().
 */
static void
exact_frame(int dstar, const int b[9], int L, const int J[4], const struct fhistory* h, struct fhistory* out) {
    const float (*prba24)[3] = dstar ? AmbePlusPRBA24 : AmbePRBA24;
    const float (*prba58)[4] = dstar ? AmbePlusPRBA58 : AmbePRBA58;
    const float (*hoc[4])[4] = {dstar ? AmbePlusHOCb5 : AmbeHOCb5, dstar ? AmbePlusHOCb6 : AmbeHOCb6,
                                dstar ? AmbePlusHOCb7 : AmbeHOCb7, dstar ? AmbePlusHOCb8 : AmbeHOCb8};
    const float rho = dstar ? MBE_AMBE2400_PREDICTION_RHO : (float)0.65;
    float G[9], R[9], C[5][18], T[57] = {0}, p[57];
    memset(C, 0, sizeof(C));
    out->gamma = dstar ? (MBE_AMBE2400_GAIN_SCALE * AmbePlusDg[b[2]]) + MBE_AMBE2400_GAIN_OFFSET
                       : AmbeDg[b[2]] + ((float)0.5 * h->gamma);
    G[1] = 0.0f;
    for (int m = 0; m < 3; ++m) {
        G[2 + m] = prba24[b[3]][m];
    }
    for (int m = 0; m < 4; ++m) {
        G[5 + m] = prba58[b[4]][m];
    }
    for (int i = 1; i <= 8; ++i) {
        float sum = 0.0f;
        for (int m = 1; m <= 8; ++m) {
            sum = sum + ((float)((m == 1) ? 1 : 2) * G[m] * dct8_cos(m, i));
        }
        R[i] = sum;
    }
    const float r = ((float)1 / ((float)2 * M_SQRT2));
    for (int i = 1; i <= 4; ++i) {
        const int odd = i + i - 1;
        C[i][1] = (float)0.5 * (R[odd] + R[odd + 1]);
        C[i][2] = r * (R[odd] - R[odd + 1]);
        for (int k = 3; k <= J[i - 1]; ++k) {
            C[i][k] = (k > 6) ? 0.0f : hoc[i - 1][b[4 + i]][k - 3];
        }
    }
    int l = 1;
    for (int i = 1; i <= 4; ++i) {
        const int ji = J[i - 1];
        for (int j = 1; j <= ji; ++j, ++l) {
            float sum = 0.0f;
            for (int k = 1; k <= ji; ++k) {
                sum = sum + ((float)((k == 1) ? 1 : 2) * C[i][k] * block_cos(ji, j, k));
            }
            T[l] = sum;
        }
    }
    memcpy(p, h->log2Ml, sizeof(p));
    for (int k = h->L + 1; k <= L; ++k) {
        p[k] = p[h->L];
    }
    p[0] = p[1];
    int lower[57];
    float delta[57], sum43 = 0.0f, sum42 = 0.0f;
    for (int k = 1; k <= L; ++k) {
        const float ratio = (float)h->L / (float)L;
        const float position = ratio * (float)k;
        lower[k] = (int)position;
        delta[k] = position - (float)lower[k];
        const int up = (lower[k] + 1 > 56) ? 56 : lower[k] + 1;
        sum43 = sum43 + ((((float)1 - delta[k]) * p[lower[k]]) + (delta[k] * p[up]));
    }
    sum43 = ((rho / (float)L) * sum43);
    for (int k = 1; k <= L; ++k) {
        sum42 += T[k];
    }
    sum42 = sum42 / (float)L;
    const float big_gamma = out->gamma - (0.5f * log2f((float)L)) - sum42;
    out->L = L;
    memset(out->log2Ml, 0, sizeof(out->log2Ml));
    for (int k = 1; k <= L; ++k) {
        const int up = (lower[k] + 1 > 56) ? 56 : lower[k] + 1;
        const float c1 = (rho * ((float)1 - delta[k]) * p[lower[k]]);
        const float c2 = (rho * delta[k] * p[up]);
        out->log2Ml[k] = T[k] + c1 + c2 - sum43 + big_gamma;
    }
}

/*
 * Exact comparison: in a build whose test and library both evaluate float arithmetic as written (no fast-math in the
 * library, FLT_EVAL_METHOD 0, contraction off: the Debug builds of GCC and Clang on targets without excess precision)
 * the decode equals exact_frame() bit for bit. Elsewhere (fast-math, x87 excess precision, contraction) the two may
 * round apart; there the bound is EXACT_BOUND, measured below.
 */
#if FLT_EVAL_METHOD == 0 && !defined(MBELIB_TEST_LIBRARY_FAST_MATH) && defined(MBELIB_TEST_FP_CONTRACT_OFF)
#define EXACT_ARITHMETIC 1
#else
#define EXACT_ARITHMETIC 0
#endif
/* The largest difference measured where the comparison is not exact (log2Ml of a few units and gamma): 5.72e-6 in
 * Clang's fast-math LTO builds, 5.25e-6 in GCC's (dev-release), 1.43e-6 in the 32-bit Release build, 0 in the 32-bit
 * Debug build; the bound is about three times the largest. A wrong row lies 0.4 or more away; an error of a few ulps
 * (1e-6, say) passes there and is caught only where the comparison is exact. */
#define EXACT_BOUND 1.6e-5

static int
float_bits_equal(float a, float b) {
    uint32_t x, y;
    memcpy(&x, &a, sizeof(x));
    memcpy(&y, &b, sizeof(y));
    return x == y;
}

static int
finite_float(float x) {
    uint32_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return (bits & 0x7FFFFFFFu) < 0x7F800000u;
}

/* One frame against exact_frame(); `worst` the largest difference seen (finite values only). Returns failures. */
static int
check_exact(int dstar, int exact, const int b[9], int L, mbe_parms* cur, mbe_parms* prev, struct fhistory* h,
            double* worst) {
    char d[49];
    struct fhistory want;
    int J[4];
    if (dstar) {
        dstar_pack(b, d);
    } else {
        mbe_ambe2450_pack(b, d);
    }
    const int rc = dstar ? mbe_decodeAmbe2400Parms(d, cur, prev) : mbe_decodeAmbe2450Parms(d, cur, prev);
    if (rc != 0 || cur->L != L) {
        printf("%s exact: b0 %d decoded rc %d L %d\n", dstar ? "D-STAR" : "rate 33", b[0], rc, cur->L);
        return 1;
    }
    expected_row(L, J);
    exact_frame(dstar, b, L, J, h, &want);
    int fails = 0;
    for (int k = 0; k <= L; ++k) {
        const float got = (k == 0) ? cur->gamma : cur->log2Ml[k];
        const float ref = (k == 0) ? want.gamma : want.log2Ml[k];
        if (!finite_float(got) || !finite_float(ref)) {
            fails++;
            continue;
        }
        const double diff = fabs((double)got - (double)ref);
        *worst = fmax(*worst, diff);
        if (exact ? !float_bits_equal(got, ref) : !(diff <= EXACT_BOUND)) {
            fails++;
        }
    }
    if (fails) {
        printf("%s exact: b0 %d b1 %d L %d: %d values off the exact result\n", dstar ? "D-STAR" : "rate 33", b[0], b[1],
               L, fails);
    }
    h->L = L;
    h->gamma = cur->gamma;
    memcpy(h->log2Ml, cur->log2Ml, sizeof(h->log2Ml));
    mbe_moveMbeParms(cur, prev);
    return fails;
}

static int
test_exact(int dstar, int exact) {
    static const int dstar_rows[2][9] = {{0, 15, 40, 301, 72, 3, 11, 9, 6}, {0, 15, 44, 117, 21, 12, 5, 14, 2}};
    static const int r33_rows[2][9] = {{0, 0, 20, 301, 72, 19, 11, 9, 5}, {0, 0, 22, 117, 21, 6, 5, 14, 2}};
    const int* vuv = dstar ? dstar_vuv : r33_vuv;
    const int nvuv = dstar ? 4 : 12;
    int fails = 0, frames = 0;
    double worst = 0.0;
    for (int first = 19; first <= 21; ++first) {
        for (int second = 19; second <= 21; ++second) {
            for (int v = 0; v < nvuv; ++v) {
                mbe_parms cur, prev, enhanced;
                struct history unused;
                struct fhistory h;
                int b[9];
                mbe_initMbeParms(&cur, &prev, &enhanced);
                shaped_history(&prev, &unused);
                h.L = prev.L;
                h.gamma = prev.gamma;
                memcpy(h.log2Ml, prev.log2Ml, sizeof(h.log2Ml));
                memcpy(b, dstar ? dstar_rows[0] : r33_rows[0], sizeof(b));
                b[0] = code_of(dstar, first, v);
                b[1] = vuv[v];
                fails += check_exact(dstar, exact, b, first, &cur, &prev, &h, &worst);
                memcpy(b, dstar ? dstar_rows[1] : r33_rows[1], sizeof(b));
                b[0] = code_of(dstar, second, v + 1);
                b[1] = vuv[(v + 5) % nvuv];
                fails += check_exact(dstar, exact, b, second, &cur, &prev, &h, &worst);
                frames += 2;
            }
        }
    }
    printf("%s decoder: %d frames against the exact result %s (largest difference %.3g)\n",
           dstar ? "D-STAR" : "rate-33", frames, exact ? "by bit pattern" : "within EXACT_BOUND", worst);
    return fails;
}

/* ---- The encoders. ---- */

/* Voiced input at about 180 Hz (20 harmonics in both codecs): harmonics of f0(t) = 181 + 12 sin(0.7 t) Hz under a
 * tilt, in runs of 25 voiced frames and 25 of noise. */
static void
signal(float* x, int f, double* phase, uint32_t* seed) {
    for (int n = 0; n < 160; ++n) {
        const double t = (double)((f * 160) + n) / 8000.0;
        const double f0 = 181.0 + (12.0 * sin(0.7 * t));
        *phase += 2.0 * M_PI * f0 / 8000.0;
        double v = 0.0;
        if (((f / 25) % 2) == 0) {
            for (int k = 1; (k * f0) < 3800.0; ++k) {
                v += 0.055 * sin((k * *phase) + (0.3 * k)) / (0.6 + (0.4 * k)) * (1.0 + (0.8 * cos(k * f0 / 700.0)));
            }
        }
        *seed = (*seed * 1664525u) + 1013904223u;
        const double noise = ((double)(*seed >> 16) / 32768.0) - 1.0;
        x[n] = (float)(v + (noise * ((((f / 25) % 2) == 0) ? 0.01 : 0.03)));
    }
}

/* The fields the encoders document as the decoder's (mbelib.h, mbe_encodeAmbe2400Parms(), mbe_encodeAmbe2450Parms()):
 * w0, L, gamma, the frame mark PHIl[0], and Vl, log2Ml and Ml over harmonics 1..L, by bit pattern. K is not among
 * them (the encoder writes its count of voiced harmonics there), nor the synthesis state the encoders never touch. */
static int
states_equal(const mbe_parms* a, const mbe_parms* b) {
    if (a->L != b->L || !float_bits_equal(a->w0, b->w0) || !float_bits_equal(a->gamma, b->gamma)) {
        return 0;
    }
    for (int l = 1; l <= a->L && l <= 56; ++l) {
        if (a->Vl[l] != b->Vl[l] || !float_bits_equal(a->log2Ml[l], b->log2Ml[l])
            || !float_bits_equal(a->Ml[l], b->Ml[l])) {
            return 0;
        }
    }
    return float_bits_equal(a->PHIl[0], b->PHIl[0]);
}

/*
 * The encoder's returned model equals the public decoder's decode of the encoder's bits, frame for frame and by
 * bit pattern (states_equal()), over a signal that sends many frames of 20 harmonics. What it shows: the encoder
 * predicts from the history the decoder will hold. What it does not: that either uses the right row. Rate 33's
 * encoder takes its model from the same public decoder, and a codec's table header serves its decoder and encoder
 * alike, so a wrong row in a header (or a decoder error) passes it; only a table copied into one translation unit
 * of D-STAR's breaks it. The rows are pinned by test_tables(), test_decoder() and test_exact().
 * An encoder error stops the test before its bits are read; a tone frame (D-STAR's decoder returns the tone index,
 * rate 33's MBE_AMBE2450_FRAME_TONE) is skipped, the encoder's model then being its history; any other decoder
 * status is an error.
 */
static int
test_encoder_parity(int dstar) {
    mbe_ambe2400_encoder* e1 = dstar ? mbe_ambe2400EncoderAlloc() : NULL;
    mbe_ambe2450_encoder* e2 = dstar ? NULL : mbe_ambe2450EncoderAlloc();
    if ((dstar && e1 == NULL) || (!dstar && e2 == NULL)) {
        return 1;
    }
    mbe_parms e_cur, e_prev, e_enh, d_cur, d_prev, d_enh;
    mbe_initMbeParms(&e_cur, &e_prev, &e_enh);
    mbe_initMbeParms(&d_cur, &d_prev, &d_enh);
    double phase = 0.0;
    uint32_t seed = 0x5EEDu;
    int frames20 = 0, mismatches = 0, errors = 0, tones = 0;
    for (int f = 0; f < 1000; ++f) {
        float x[160];
        char bits[49];
        signal(x, f, &phase, &seed);
        const int rc = dstar ? mbe_encodeAmbe2400Parms(e1, x, bits, &e_cur, &e_prev)
                             : mbe_encodeAmbe2450Parms(e2, x, bits, &e_cur, &e_prev);
        if (rc < 0) {
            printf("%s parity: the encoder returned %d at frame %d\n", dstar ? "D-STAR" : "rate-33", rc, f);
            errors++;
            break;
        }
        const int dr =
            dstar ? mbe_decodeAmbe2400Parms(bits, &d_cur, &d_prev) : mbe_decodeAmbe2450Parms(bits, &d_cur, &d_prev);
        const int tone = dstar ? (dr >= 5) : (dr == MBE_AMBE2450_FRAME_TONE && rc == MBE_AMBE2450_FRAME_TONE);
        if (tone) {
            tones++;
            mbe_moveMbeParms(&e_cur, &e_prev);
            continue;
        }
        if (dr != 0 || (!dstar && rc != MBE_AMBE2450_FRAME_VOICE)) {
            printf("%s parity: frame %d encoded with status %d, decoded with %d\n", dstar ? "D-STAR" : "rate-33", f, rc,
                   dr);
            errors++;
            break;
        }
        frames20 += (d_cur.L == 20);
        mismatches += !states_equal(&e_cur, &d_cur);
        mbe_moveMbeParms(&d_cur, &d_prev);
        mbe_moveMbeParms(&e_cur, &e_prev);
    }
    if (dstar) {
        mbe_ambe2400EncoderFree(e1);
    } else {
        mbe_ambe2450EncoderFree(e2);
    }
    printf("%s encoder/decoder parity: %d mismatching frames, %d of them of 20 harmonics in 1000 (%d tone frames)\n",
           dstar ? "D-STAR" : "rate-33", mismatches, frames20, tones);
    /* The signal sends frames of 20 harmonics (measured: 147 in D-STAR, 120 in rate 33). */
    return errors != 0 || mismatches != 0 || frames20 < 100;
}

static uint32_t rng = 0xB10Cu;

static uint32_t
rnd(void) {
    rng = (rng * 1664525u) + 1013904223u;
    return rng >> 8;
}

/* Rms over l of the decoded mean-removed residual (the search's linear model of the decoder,
 * ambe_enc_model_residual(), at block lengths J) less the target's; -1 if it is not finite. */
static double
residual_error(const struct ambe_enc_tables* t, const struct ambe_enc_frame* q, const int J[5], const double* target) {
    float x[AMBE_ENC_COEFFS], model[57];
    ambe_enc_coefficients(t, q->b, J, x);
    ambe_enc_model_residual(J, q->L, x, model);
    double sum = 0.0;
    for (int l = 1; l <= q->L; ++l) {
        const double e = (double)model[l] - target[l];
        sum += e * e;
    }
    const double rms = sqrt(sum / q->L);
    return finite_double(rms) ? rms : -1.0;
}

/* Code a residual of 20 harmonics through the encoder's quantizer (nearest rows, then the joint refinement) with
 * the tables t, and return the rms error of what a decoder with t's block lengths reconstructs. */
static double
round_trip(const struct ambe_enc_tables* t, const double* target, const mbe_parms* prev, float level) {
    struct ambe_enc_frame q;
    memset(&q, 0, sizeof(q));
    q.cache = mbe_ambe2400_get_dct_cache();
    q.L = 20;
    for (int l = 1; l <= 20; ++l) {
        q.a[l] = (float)target[l] + level;
    }
    ambe_enc_mean_amplitude(&q);
    ambe_enc_prediction(&q, t, prev);
    ambe_enc_quantize_nearest(&q, t);
    ambe_enc_refine(&q, t);
    int J[5] = {0, t->lmprbl[20][0], t->lmprbl[20][1], t->lmprbl[20][2], t->lmprbl[20][3]};
    return residual_error(t, &q, J, target);
}

/*
 * A residual that blocks of 4, 5, 5, 6 represent exactly (random receivable rows b3..b8, decoded by the search's
 * model of the decoder at those lengths, mean removed: the third block's fifth coefficient is the HOC row b7's
 * third element, which a block of 6 would place elsewhere) round-trips with the encoder's own tables within
 * 1e-4 rms on every frame, and with the tables' row 20 put back to the annex's 4, 4, 6, 6 (encoder and decoder
 * both) at least 0.05 rms off on the mean of the frames.
 */
static int
test_round_trip(int dstar) {
    const struct ambe_enc_tables* own = dstar ? mbe_ambe2400_enc_tables() : mbe_ambe2450_enc_tables();
    static int annex_rows[57][4];
    memcpy(annex_rows, own->lmprbl, sizeof(annex_rows));
    memcpy(annex_rows[20], annex_c[20 - 9], sizeof(annex_rows[20]));
    struct ambe_enc_tables annex = *own;
    annex.lmprbl = (const int (*)[4])annex_rows;
    mbe_parms prev, cur, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int l = 0; l <= 56; ++l) {
        prev.log2Ml[l] = 0.0f;
    }
    prev.L = 20;
    prev.gamma = 0.0f;
    const int J[5] = {0, 4, 5, 5, 6};
    double worst_own = 0.0, sum_annex = 0.0;
    int bad = 0;
    const int frames = 200;
    for (int f = 0; f < frames; ++f) {
        struct ambe_enc_frame q;
        memset(&q, 0, sizeof(q));
        q.L = 20;
        q.b[3] = (int)(rnd() % 512u);
        q.b[4] = (int)(rnd() % 128u);
        for (int blk = 0; blk < 4; ++blk) {
            const int row = (int)(rnd() % (uint32_t)own->hoc_rows[blk]);
            q.b[5 + blk] = row - (row % own->hoc_step[blk]);
        }
        float x[AMBE_ENC_COEFFS], model[57];
        double target[57] = {0}, mean = 0.0;
        ambe_enc_coefficients(own, q.b, J, x);
        ambe_enc_model_residual(J, 20, x, model);
        for (int l = 1; l <= 20; ++l) {
            mean += (double)model[l] / 20.0;
        }
        for (int l = 1; l <= 20; ++l) {
            target[l] = (double)model[l] - mean;
        }
        const float level = 1.0f + (float)(rnd() % 400u) / 100.0f;
        /* Each error is checked before it is aggregated: fmax() would drop a NaN. */
        const double e_own = round_trip(own, target, &prev, level);
        const double e_annex = round_trip(&annex, target, &prev, level);
        if (!(e_own >= 0.0) || !(e_annex >= 0.0)) {
            bad++;
            continue;
        }
        worst_own = fmax(worst_own, e_own);
        sum_annex += e_annex;
    }
    const double mean_annex = sum_annex / frames;
    printf("%s round trip of residuals blocks of 4, 5, 5, 6 represent: worst rms error %.2g with the encoder's tables, "
           "mean %.3f with the annex's row 20\n",
           dstar ? "D-STAR" : "rate-33", worst_own, mean_annex);
    if (bad) {
        printf("%s round trip: %d frames gave an error that is not finite\n", dstar ? "D-STAR" : "rate-33", bad);
    }
    return bad != 0 || !(worst_own <= 1e-4) || !(mean_annex >= 0.05);
}

/* Uniform in [lo, hi). */
static float
uniform(float lo, float hi) {
    return lo + ((hi - lo) * ((float)(rnd() & 0xFFFFu) / 65536.0f));
}

/* Target log2 amplitudes of 20 harmonics: a few formant-like peaks over a tilt plus roughness (as
 * test_ambe_vector_search's frames), and a history of the same kind at 18 to 22 harmonics. */
static void
random_target(float a[57], mbe_parms* prev) {
    mbe_parms cur, enhanced;
    mbe_initMbeParms(&cur, prev, &enhanced);
    prev->L = 18 + (int)(rnd() % 5u);
    prev->gamma = uniform(0.0f, 8.0f);
    const float base = uniform(-1.0f, 6.0f), tilt = uniform(-0.15f, 0.0f);
    float centre[3], depth[3];
    for (int k = 0; k < 3; ++k) {
        centre[k] = uniform(0.0f, 1.0f);
        depth[k] = uniform(0.0f, 4.0f);
    }
    memset(a, 0, 57 * sizeof(a[0]));
    for (int l = 1; l <= 56; ++l) {
        float v = base + (tilt * (float)l) + uniform(-0.7f, 0.7f);
        for (int k = 0; k < 3; ++k) {
            const float u = (((float)l / 20.0f) - centre[k]) * 8.0f;
            v += depth[k] / (1.0f + (u * u));
        }
        a[l] = (l <= 20) ? v : 0.0f;
        prev->log2Ml[l] = (rnd() & 1u) ? v + uniform(-1.0f, 1.0f) : 0.0f;
    }
}

/* The refinement's weighted error of the rows q->b against q's own target (ambe_enc_search.c's definition: weights
 * min(10, max(0.1, 2^((a_l - mean a) / 2))), residuals less their means); -1 if not finite. */
static double
weighted_error(const struct ambe_enc_tables* t, const struct ambe_enc_frame* q) {
    float x[AMBE_ENC_COEFFS], model[57];
    ambe_enc_coefficients(t, q->b, q->Ji, x);
    ambe_enc_model_residual(q->Ji, q->L, x, model);
    double mean_t = 0.0, e = 0.0;
    for (int l = 1; l <= q->L; ++l) {
        mean_t += (double)q->Tl[l] / q->L;
    }
    for (int l = 1; l <= q->L; ++l) {
        const double w = fmin(10.0, fmax(0.1, pow(2.0, 0.5 * ((double)q->a[l] - (double)q->mean_a))));
        const double d = (double)model[l] - ((double)q->Tl[l] - mean_t);
        e += w * d * d;
    }
    return finite_double(e) ? e : -1.0;
}

/*
 * The joint refinement acts at L = 20. On frames of 20 harmonics that the nearest rows do not reproduce, the shared
 * quantizer (ambe_enc_quantize_residual(), which both encoders call) changes the nearest rows on most frames and
 * lowers the weighted error in all, for D-STAR's and rate 33's tables; and rate 33's codec path
 * (mbe_ambe2450_quantize_amplitudes(), at b0 52) sends the refined rows, not the nearest. D-STAR's own encode path
 * is static: a change there that bypassed the shared quantizer is not reached here.
 */
static int
test_refinement(int dstar) {
    const struct ambe_enc_tables* t = dstar ? mbe_ambe2400_enc_tables() : mbe_ambe2450_enc_tables();
    const int frames = 300;
    int changed = 0, codec_changed = 0, codec_mismatch = 0, fails = 0;
    double open_sum = 0.0, refined_sum = 0.0;
    for (int f = 0; f < frames; ++f) {
        float a[57];
        mbe_parms prev;
        struct ambe_enc_frame q;
        random_target(a, &prev);
        memset(&q, 0, sizeof(q));
        q.cache = mbe_ambe2400_get_dct_cache();
        q.L = 20;
        q.b[0] = dstar ? 55 : 52;
        q.f0q = dstar ? mbe_ambe2400_f0(55) : mbe_ambe2450_f0(52);
        memcpy(q.a, a, sizeof(q.a));
        ambe_enc_mean_amplitude(&q);
        ambe_enc_quantize_gain(&q, t, &prev);
        ambe_enc_prediction(&q, t, &prev);
        struct ambe_enc_frame open = q;
        ambe_enc_quantize_nearest(&open, t);
        ambe_enc_quantize_residual(&q, t);
        const double e_open = weighted_error(t, &open), e_refined = weighted_error(t, &q);
        if (!(e_open >= 0.0) || !(e_refined >= 0.0) || e_refined > e_open + (1e-4 * e_open) + 1e-6) {
            fails++;
            continue;
        }
        changed += memcmp(q.b, open.b, sizeof(q.b)) != 0;
        open_sum += e_open;
        refined_sum += e_refined;
        if (!dstar) {
            char d[49];
            int b[9];
            mbe_ambe2450_quantize_amplitudes(52, 0, a, &prev, d);
            mbe_ambe2450_unpack(d, b);
            codec_changed += memcmp(b + 3, open.b + 3, 6 * sizeof(int)) != 0;
            codec_mismatch += memcmp(b + 2, q.b + 2, 7 * sizeof(int)) != 0;
        }
    }
    printf("%s refinement at L 20: %d of %d frames change rows, weighted error %.4g -> %.4g",
           dstar ? "D-STAR" : "rate-33", changed, frames, open_sum, refined_sum);
    if (!dstar) {
        printf("; the codec path sends other rows than the nearest on %d, the shared quantizer's on all but %d",
               codec_changed, codec_mismatch);
    }
    printf("\n");
    /* Measured: 69 to 75% of frames change and the error falls by about 5% in all. */
    fails += changed < frames / 2 || !(refined_sum <= 0.98 * open_sum);
    if (!dstar) {
        fails += codec_changed < frames / 2 || codec_mismatch != 0;
    }
    return fails;
}

int
main(void) {
    /* Passed as an argument, not tested as a constant around its only use (MSVC's C4189). */
    const int exact = EXACT_ARITHMETIC;
    int fails = test_tables();
    for (int dstar = 1; dstar >= 0; --dstar) {
        fails += test_decoder(dstar);
        fails += test_exact(dstar, exact);
        fails += test_encoder_parity(dstar);
        fails += test_refinement(dstar);
        fails += test_round_trip(dstar);
    }
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}
