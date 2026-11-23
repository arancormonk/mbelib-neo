// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Joint, amplitude-weighted choice of the AMBE encoders' six spectral
 *        vectors (b3..b8).
 *
 * TIA-102.BABA-1 4.3.1 and 4.3.2 pick each vector as the row nearest its own
 * target (minimum mean-squared error in the transform domain); this library's
 * encoders start there (ambe_enc_quantize_nearest()). An error in a transform
 * coefficient lands on the harmonics unevenly, though, and one on a spectral
 * peak counts no more than one in a valley. So the six indices are then
 * refined together against the error the decoder will produce:
 *
 *   E = sum_l w_l (That_l - mean(That) - (T_l - mean(T)))^2,
 *   w_l = min(10, max(0.1, 2^((a_l - mean(a)) / 2))),
 *
 * where T is the prediction residual to code, That the decoder's residual for
 * the chosen rows (it uses That less its mean; the gain carries the mean) and
 * a the target log2 amplitudes. Two passes visit b3, b4, b5, b6, b7 and b8 in
 * turn; each visit takes the receivable row with the smallest E, the others
 * held, keeping the current row unless another is strictly better (ties go to
 * the lower index). E never rises. docs/dvsi-parity.md, "Joint choice of the
 * AMBE spectral vectors", has the measurements.
 *
 * The decoder's residual is linear in 23 coefficients: G2..G8 of the two PRBA
 * vectors and C(i,3..6) of each block's higher-order vector. On block i it is
 * a sum of the block's cosines f_1 = 1, f_k(j) = 2 cos(pi (k-1)(j-1/2)/J_i),
 * k <= min(J_i, 6), with C(i,1) and C(i,2) from the PRBA vector and C(i,k>2)
 * from the higher-order one; the cosines k >= 2 sum to zero over a block, so
 * the residual's mean is sum_i J_i C(i,1) / L. E is therefore a quadratic in
 * the coefficients built from each block's weighted moments of its cosines,
 * and one codebook's rows are scored by a quadratic in at most four numbers.
 */

#include <math.h>
#include <string.h>

#include "ambe3600x2400_internal.h"
#include "ambe_encoder.h"

#define AMBE_ENC_G          7 /* G2..G8 at coefficients 0..6 */
#define AMBE_ENC_BLOCK_BASE 6 /* cosines per block: f_1..f_6 */
#define AMBE_ENC_MAX_DIM    4 /* widest codebook */

/* Coefficient of C(blk, k), blk 1..4, k 3..6. */
static inline int
ambe_enc_hoc_coefficient(int blk, int k) {
    return AMBE_ENC_G + (4 * (blk - 1)) + (k - 3);
}

/* Cosines of block blk: min(J, 6). */
static inline int
ambe_enc_block_basis(int ji) {
    return (ji < AMBE_ENC_BLOCK_BASE) ? ji : AMBE_ENC_BLOCK_BASE;
}

/*
 * The PRBA part of each block: C(i,1) less the frame mean (first) and C(i,2)
 * (second) per unit of G_m. R_i = 2 sum_m G_m cos(pi (m-1)(i-1/2)/8), C(i,1) =
 * (R_{2i-1} + R_{2i}) / 2 and C(i,2) = (R_{2i-1} - R_{2i}) / (2 sqrt 2), as the
 * decoders reconstruct them (mbe_ambe2400_reconstruct_ri/_cik()).
 */
struct ambe_enc_prba_map {
    float first[5][AMBE_ENC_G];
    float second[5][AMBE_ENC_G];
};

static void
ambe_enc_prba_map(const struct ambe_dct_cache* cache, const int Ji[5], int L, struct ambe_enc_prba_map* map) {
    const float inv_sqrt2 = 0.70710678118654752f;
    for (int m = 0; m < AMBE_ENC_G; m++) {
        float mean = 0.0f;
        /* R_{2i-1} and R_{2i} of block i sit at upper_i and upper_i + 1. */
        for (int blk = 1, upper_i = 1; blk <= 4; blk++, upper_i += 2) {
            const float upper = cache->ri_cos[m + 2][upper_i];
            const float lower = cache->ri_cos[m + 2][upper_i + 1];
            map->first[blk][m] = upper + lower;
            map->second[blk][m] = (upper - lower) * inv_sqrt2;
            mean += (float)Ji[blk] * map->first[blk][m];
        }
        mean /= (float)L;
        for (int blk = 1; blk <= 4; blk++) {
            map->first[blk][m] -= mean;
        }
    }
}

/* Cosine f_k of harmonic j of a block of length ji (k >= 1). */
static inline float
ambe_enc_cosine(const struct ambe_dct_cache* cache, int ji, int j, int k) {
    return (k == 1) ? 1.0f : 2.0f * cache->idct_cos[ji][j][k];
}

void
ambe_enc_coefficients(const struct ambe_enc_tables* t, const int b[9], const int Ji[5], float x[AMBE_ENC_COEFFS]) {
    memset(x, 0, AMBE_ENC_COEFFS * sizeof(x[0]));
    for (int m = 0; m < 3; m++) {
        x[m] = t->prba24[b[3]][m];
    }
    for (int m = 0; m < 4; m++) {
        x[3 + m] = t->prba58[b[4]][m];
    }
    for (int blk = 1; blk <= 4; blk++) {
        for (int k = 3; k <= ambe_enc_block_basis(Ji[blk]); k++) {
            x[ambe_enc_hoc_coefficient(blk, k)] = t->hoc[blk - 1][b[4 + blk]][k - 3];
        }
    }
}

void
ambe_enc_model_residual(const int Ji[5], int L, const float x[AMBE_ENC_COEFFS], float residual[57]) {
    const struct ambe_dct_cache* cache = mbe_ambe2400_get_dct_cache();
    struct ambe_enc_prba_map map;
    ambe_enc_prba_map(cache, Ji, L, &map);
    int l = 1;
    for (int blk = 1; blk <= 4; blk++) {
        float y[AMBE_ENC_BLOCK_BASE + 1] = {0.0f};
        for (int m = 0; m < AMBE_ENC_G; m++) {
            y[1] += map.first[blk][m] * x[m];
            y[2] += map.second[blk][m] * x[m];
        }
        const int n = ambe_enc_block_basis(Ji[blk]);
        for (int k = 3; k <= n; k++) {
            y[k] = x[ambe_enc_hoc_coefficient(blk, k)];
        }
        for (int j = 1; j <= Ji[blk]; j++, l++) {
            float sum = 0.0f;
            for (int k = 1; k <= n; k++) {
                sum += y[k] * ambe_enc_cosine(cache, Ji[blk], j, k);
            }
            residual[l] = sum;
        }
    }
}

/* E(x) = x' H x - 2 g' x + const over the 23 coefficients. */
struct ambe_enc_quadratic {
    float H[AMBE_ENC_COEFFS][AMBE_ENC_COEFFS];
    float g[AMBE_ENC_COEFFS];
};

/* Weighted moments of block blk's cosines: M[a][b] = sum_j w f_a f_b and
 * v[a] = sum_j w f_a t (a, b from 1). */
static void
ambe_enc_block_moments(const struct ambe_dct_cache* cache, int ji, const float* w, const float* t,
                       float M[AMBE_ENC_BLOCK_BASE + 1][AMBE_ENC_BLOCK_BASE + 1], float v[AMBE_ENC_BLOCK_BASE + 1]) {
    const int n = ambe_enc_block_basis(ji);
    memset(M, 0, (AMBE_ENC_BLOCK_BASE + 1) * sizeof(M[0]));
    memset(v, 0, (AMBE_ENC_BLOCK_BASE + 1) * sizeof(v[0]));
    for (int j = 1; j <= ji; j++) {
        float f[AMBE_ENC_BLOCK_BASE + 1] = {0.0f};
        for (int k = 1; k <= n; k++) {
            f[k] = ambe_enc_cosine(cache, ji, j, k);
        }
        for (int a = 1; a <= n; a++) {
            const float wf = w[j] * f[a];
            v[a] += wf * t[j];
            for (int b = a; b <= n; b++) {
                M[a][b] += wf * f[b];
            }
        }
    }
    for (int a = 1; a <= n; a++) {
        for (int b = 1; b < a; b++) {
            M[a][b] = M[b][a];
        }
    }
}

/* Add block blk's share, P' M P and P' v, where P maps the coefficients to
 * the block's cosine weights. */
static void
ambe_enc_add_block(struct ambe_enc_quadratic* e, const struct ambe_enc_prba_map* map, int blk, int n,
                   float M[AMBE_ENC_BLOCK_BASE + 1][AMBE_ENC_BLOCK_BASE + 1], const float v[AMBE_ENC_BLOCK_BASE + 1]) {
    const float* first = map->first[blk];
    const float* second = map->second[blk];
    /* Q[a][m]: cosine a's moment against the PRBA part, per G_m. */
    float Q[AMBE_ENC_BLOCK_BASE + 1][AMBE_ENC_G] = {{0.0f}};
    for (int a = 1; a <= n; a++) {
        for (int m = 0; m < AMBE_ENC_G; m++) {
            Q[a][m] = (M[a][1] * first[m]) + ((n >= 2) ? M[a][2] * second[m] : 0.0f);
        }
    }
    for (int m = 0; m < AMBE_ENC_G; m++) {
        for (int p = 0; p < AMBE_ENC_G; p++) {
            e->H[m][p] += (first[m] * Q[1][p]) + ((n >= 2) ? second[m] * Q[2][p] : 0.0f);
        }
        e->g[m] += (v[1] * first[m]) + ((n >= 2) ? v[2] * second[m] : 0.0f);
    }
    for (int k = 3; k <= n; k++) {
        const int c = ambe_enc_hoc_coefficient(blk, k);
        for (int m = 0; m < AMBE_ENC_G; m++) {
            e->H[c][m] = Q[k][m];
            e->H[m][c] = Q[k][m];
        }
        for (int k2 = 3; k2 <= n; k2++) {
            e->H[c][ambe_enc_hoc_coefficient(blk, k2)] = M[k][k2];
        }
        e->g[c] = v[k];
    }
}

/* The amplitude weight of harmonic l. */
static inline float
ambe_enc_weight(float a, float mean_a) {
    return fminf(10.0f, fmaxf(0.1f, exp2f(0.5f * (a - mean_a))));
}

static void
ambe_enc_build_quadratic(const struct ambe_enc_frame* q, struct ambe_enc_quadratic* e) {
    float w[57] = {0.0f};
    float t[57] = {0.0f};
    float mean_t = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        mean_t += q->Tl[l];
    }
    mean_t /= (float)q->L;
    for (int l = 1; l <= q->L; l++) {
        w[l] = ambe_enc_weight(q->a[l], q->mean_a);
        t[l] = q->Tl[l] - mean_t;
    }
    struct ambe_enc_prba_map map;
    ambe_enc_prba_map(q->cache, q->Ji, q->L, &map);
    memset(e, 0, sizeof(*e));
    int l = 1;
    for (int blk = 1; blk <= 4; blk++) {
        float M[AMBE_ENC_BLOCK_BASE + 1][AMBE_ENC_BLOCK_BASE + 1];
        float v[AMBE_ENC_BLOCK_BASE + 1];
        /* w and t from harmonic l on index this block's harmonics from 1. */
        ambe_enc_block_moments(q->cache, q->Ji[blk], w + l - 1, t + l - 1, M, v);
        ambe_enc_add_block(e, &map, blk, ambe_enc_block_basis(q->Ji[blk]), M, v);
        l += q->Ji[blk];
    }
}

/*
 * One codebook's error as a function of its row y, relative to the current
 * row c: with d = y - c, E(y) - E(c) = d' A d + lin' d.
 */
struct ambe_enc_row_form {
    float A[AMBE_ENC_MAX_DIM][AMBE_ENC_MAX_DIM];
    float lin[AMBE_ENC_MAX_DIM];
};

/*
 * The receivable row (every step-th of rows) that lowers E most, from its
 * first `used` elements. The current row is kept unless another is strictly
 * better, and the first of equally good rows wins; scored by its difference
 * from the current row, a row equal to it there scores exactly 0. Inline at
 * the call sites so the fixed codebook shapes specialize the loops.
 */
static inline int
ambe_enc_best_row(const float* table, int width, int rows, int step, int used, int current,
                  const struct ambe_enc_row_form* form) {
    const float* c = table + ((size_t)current * (size_t)width);
    int best = current;
    float best_delta = 0.0f;
    for (int r = 0; r < rows; r += step) {
        const float* y = table + ((size_t)r * (size_t)width);
        float d[AMBE_ENC_MAX_DIM];
        for (int a = 0; a < used; a++) {
            d[a] = y[a] - c[a];
        }
        float delta = 0.0f;
        for (int a = 0; a < used; a++) {
            float s = form->lin[a];
            for (int b = 0; b < used; b++) {
                s += form->A[a][b] * d[b];
            }
            delta += d[a] * s;
        }
        if (delta < best_delta) {
            best_delta = delta;
            best = r;
        }
    }
    return best;
}

/*
 * One visit: with every coefficient but [first, first + used) held, E as a
 * function of that codebook's row: A is H over the visited coefficients and
 * lin = 2 (H x - g) there, the gradient at the current row.
 */
static void
ambe_enc_visit_setup(const struct ambe_enc_quadratic* e, const float x[AMBE_ENC_COEFFS], int first, int used,
                     struct ambe_enc_row_form* form) {
    for (int a = 0; a < used; a++) {
        const float* row = e->H[first + a];
        float gradient = -e->g[first + a];
        for (int c = 0; c < AMBE_ENC_COEFFS; c++) {
            gradient += row[c] * x[c];
        }
        form->lin[a] = 2.0f * gradient;
        for (int b = 0; b < used; b++) {
            form->A[a][b] = row[first + b];
        }
    }
}

/* Visit codebook `index` (3..8) and store its row's coefficients in x. */
static void
ambe_enc_visit(struct ambe_enc_frame* q, const struct ambe_enc_tables* t, const struct ambe_enc_quadratic* e,
               float x[AMBE_ENC_COEFFS], int index) {
    struct ambe_enc_row_form form;
    const float* row;
    int first;
    int used;
    if (index == 3) {
        first = 0;
        used = 3;
        ambe_enc_visit_setup(e, x, first, used, &form);
        q->b[3] = ambe_enc_best_row(&t->prba24[0][0], 3, 512, 1, 3, q->b[3], &form);
        row = t->prba24[q->b[3]];
    } else if (index == 4) {
        first = 3;
        used = 4;
        ambe_enc_visit_setup(e, x, first, used, &form);
        q->b[4] = ambe_enc_best_row(&t->prba58[0][0], 4, 128, 1, 4, q->b[4], &form);
        row = t->prba58[q->b[4]];
    } else {
        const int blk = index - 4;
        used = ambe_enc_block_basis(q->Ji[blk]) - 2;
        if (used <= 0) {
            return; /* no higher-order coefficient: the index stays 0 */
        }
        first = ambe_enc_hoc_coefficient(blk, 3);
        ambe_enc_visit_setup(e, x, first, used, &form);
        q->b[index] = ambe_enc_best_row(&t->hoc[blk - 1][0][0], 4, t->hoc_rows[blk - 1], t->hoc_step[blk - 1], used,
                                        q->b[index], &form);
        row = t->hoc[blk - 1][q->b[index]];
    }
    for (int a = 0; a < used; a++) {
        x[first + a] = row[a];
    }
}

void
ambe_enc_refine(struct ambe_enc_frame* q, const struct ambe_enc_tables* t) {
    struct ambe_enc_quadratic e;
    float x[AMBE_ENC_COEFFS];
    ambe_enc_build_quadratic(q, &e);
    ambe_enc_coefficients(t, q->b, q->Ji, x);
    for (int pass = 0; pass < 2; pass++) {
        for (int index = 3; index <= 8; index++) {
            ambe_enc_visit(q, t, &e, x, index);
        }
    }
}
