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
 * - The decoders: two frames built from chosen codes, at 19, 20 and 21
 *   harmonics, through the public parameter decoders, against the block
 *   transform and prediction formed here from the codebooks and the lengths
 *   stated here; the second frame is predicted from the first, so the
 *   prediction's use of a 20-harmonic frame is covered. The chosen codes keep
 *   every other division of the frame's harmonics into blocks that moves one
 *   harmonic between two neighbouring blocks (and 5, 5, 5, 5 at L = 20) at
 *   least MARGIN away from the decode, so each such row in a table fails.
 * - The encoders: decoder and encoder state parity on a 20-harmonic signal
 *   (D-STAR's encoder forms its model with the decoder's steps and its own
 *   table; rate 33's runs the public decoder on its bits), and a 20-harmonic
 *   residual that blocks of 4, 5, 5, 6 represent exactly round-trips through
 *   the quantizer with the encoders' tables and not with the annex's row.
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

#define TOLERANCE 2e-4 /* float decode against double reference: a few ulps of values up to about 10 */
#define MARGIN    0.05 /* every other division tried lies further than this from the decode */

/* Finite: at most DBL_MAX in magnitude, judged on the bit pattern so that no build's float semantics can pass a
 * NaN or an infinity. */
static int
finite_double(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return (bits & UINT64_C(0x7FFFFFFFFFFFFFFF)) < UINT64_C(0x7FF0000000000000);
}

/* The largest difference of log2Ml over the harmonics; a difference that is not finite counts as DBL_MAX. */
static double
max_difference(const struct history* a, const struct history* b) {
    double worst = 0.0;
    for (int k = 1; k <= a->L; ++k) {
        const double d = fabs(a->log2Ml[k] - b->log2Ml[k]);
        if (!finite_double(d)) {
            return DBL_MAX;
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
    if (!(error <= TOLERANCE)) {
        printf("%s L %d: decode differs from rows %d %d %d %d by %.3g\n", dstar ? "D-STAR" : "rate 33", L, J[0], J[1],
               J[2], J[3], error);
        fails++;
    }
    const int n = alternative_rows(L, J, alt);
    for (int a = 0; a < n; ++a) {
        reference_frame(dstar, b, L, alt[a], ref_history, &other);
        const double margin = max_difference(&expected, &other);
        *least_margin = fmin(*least_margin, margin);
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

/* Codes of 19, 20 and 21 harmonics: D-STAR b0 52 (19; 51..53 give 19), 55 (20; 54..57) and 59 (21; 58..60), rate
 * 33 b0 49 (19; 48..50), 52 (20; 51..53) and 55 (21; 54..56), found by decoding every code (each check asserts the
 * decoded L). */
static int
code_of(int dstar, int L) {
    static const int codes[2][3] = {{49, 52, 55}, {52, 55, 59}};
    return codes[dstar ? 1 : 0][L - 19];
}

static int
test_decoder(int dstar) {
    /* Two frames, b2..b8 chosen once per codec (receivable rows: D-STAR's b8 even), b1 voicing all bands. */
    static const int dstar_rows[2][9] = {{0, 15, 40, 301, 72, 3, 11, 9, 6}, {0, 15, 44, 117, 21, 12, 5, 14, 2}};
    static const int r33_rows[2][9] = {{0, 0, 20, 301, 72, 19, 11, 9, 5}, {0, 0, 22, 117, 21, 6, 5, 14, 2}};
    int fails = 0;
    double least = DBL_MAX;
    for (int first = 19; first <= 21; ++first) {
        for (int second = 19; second <= 21; ++second) {
            mbe_parms cur, prev, enhanced;
            mbe_initMbeParms(&cur, &prev, &enhanced);
            struct history h;
            h.L = 20;
            h.gamma = 0.0;
            prev.L = 20;
            prev.gamma = 0.0f;
            for (int k = 0; k <= 56; ++k) {
                /* A history with a shape, so the prediction weighs. */
                h.log2Ml[k] = (k == 0) ? 0.0 : 2.0 + sin(0.7 * k);
                prev.log2Ml[k] = (float)h.log2Ml[k];
                h.log2Ml[k] = (double)prev.log2Ml[k];
            }
            int b[9];
            memcpy(b, dstar ? dstar_rows[0] : r33_rows[0], sizeof(b));
            b[0] = code_of(dstar, first);
            fails += check_frame(dstar, b, first, &cur, &prev, &h, &least);
            memcpy(b, dstar ? dstar_rows[1] : r33_rows[1], sizeof(b));
            b[0] = code_of(dstar, second);
            fails += check_frame(dstar, b, second, &cur, &prev, &h, &least);
        }
    }
    printf("%s decoder: two frames at 19, 20 and 21 harmonics follow the block lengths (least margin to another "
           "row %.3f)\n",
           dstar ? "D-STAR" : "rate-33", least);
    return fails;
}

/* ---- The encoders. ---- */

static int
float_bits_equal(float a, float b) {
    uint32_t x, y;
    memcpy(&x, &a, sizeof(x));
    memcpy(&y, &b, sizeof(y));
    return x == y;
}

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

/* The encoder's returned model equals the public decoder's decode of the encoder's bits, frame for frame and by
 * bit pattern, over a signal that sends many frames of 20 harmonics. */
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
    int frames20 = 0, mismatches = 0, errors = 0;
    for (int f = 0; f < 1000; ++f) {
        float x[160];
        char bits[49];
        signal(x, f, &phase, &seed);
        const int rc = dstar ? mbe_encodeAmbe2400Parms(e1, x, bits, &e_cur, &e_prev)
                             : mbe_encodeAmbe2450Parms(e2, x, bits, &e_cur, &e_prev);
        const int dr =
            dstar ? mbe_decodeAmbe2400Parms(bits, &d_cur, &d_prev) : mbe_decodeAmbe2450Parms(bits, &d_cur, &d_prev);
        errors += (rc != 0);
        if (dr != 0) {
            /* A tone frame: the encoder's model is its history, which the decoder leaves alone. */
            mbe_moveMbeParms(&e_cur, &e_prev);
            continue;
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
    printf("%s encoder/decoder parity: %d mismatching frames, %d of them of 20 harmonics in 1000\n",
           dstar ? "D-STAR" : "rate-33", mismatches, frames20);
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
 * ambe_enc_model_residual(), at block lengths J) less the target's. */
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
    return sqrt(sum / q->L);
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
        worst_own = fmax(worst_own, round_trip(own, target, &prev, level));
        sum_annex += round_trip(&annex, target, &prev, level);
    }
    const double mean_annex = sum_annex / frames;
    printf("%s round trip of residuals blocks of 4, 5, 5, 6 represent: worst rms error %.2g with the encoder's tables, "
           "mean %.3f with the annex's row 20\n",
           dstar ? "D-STAR" : "rate-33", worst_own, mean_annex);
    return !(worst_own <= 1e-4) || !(mean_annex >= 0.05);
}

int
main(void) {
    int fails = test_tables();
    for (int dstar = 1; dstar >= 0; --dstar) {
        fails += test_decoder(dstar);
        fails += test_encoder_parity(dstar);
        fails += test_round_trip(dstar);
    }
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}
