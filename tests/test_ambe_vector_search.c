// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/*
 * The AMBE encoders' joint choice of their six spectral vectors
 * (ambe_enc_search.c), for D-STAR's and AMBE+2's tables: its linear model of
 * the residual is the decoders', the refinement never raises the weighted
 * error of the decoded residual, keeps the rows of a frame it can reproduce
 * exactly, and sends only rows the decoder receives. The error is measured
 * here from the decoders' own reconstruction, in double precision, with
 * tolerances that leave room for any build's float evaluation order: no check
 * depends on which of two nearly equal rows wins.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ambe3600x2400_internal.h"
#include "ambe3600x2450_internal.h"
#include "ambe_encoder.h"
#include "mbelib-neo/mbelib.h"

struct codec {
    const char* name;
    const struct ambe_enc_tables* t;
    int dstar;
    int b0_of_L[57]; /* rate 33: a pitch code with L harmonics */
};

static uint32_t rng = 0x5EEDu;

static uint32_t
rnd(void) {
    rng = (rng * 1664525u) + 1013904223u;
    return rng >> 8;
}

/* Uniform in [lo, hi). */
static float
uniform(float lo, float hi) {
    return lo + ((hi - lo) * ((float)(rnd() & 0xFFFFu) / 65536.0f));
}

/* Finite: at most 2^1023 in magnitude, a NaN failing, judged on the bit
 * pattern read from memory so that no build's float semantics can pass it. */
static int
finite_double(const double* x) {
    uint64_t bits;
    memcpy(&bits, x, sizeof(bits));
    return (bits & UINT64_C(0x7FFFFFFFFFFFFFFF)) <= UINT64_C(0x7FE0000000000000);
}

static void
block_lengths(const struct codec* c, int L, int Ji[5]) {
    Ji[0] = 0;
    for (int blk = 1; blk <= 4; blk++) {
        Ji[blk] = c->t->lmprbl[L][blk - 1];
    }
}

/* Random receivable indices b3..b8: D-STAR's b8 on even rows. */
static void
random_indices(const struct codec* c, int b[9]) {
    memset(b, 0, 9 * sizeof(b[0]));
    b[3] = (int)(rnd() % 512u);
    b[4] = (int)(rnd() % 128u);
    for (int blk = 0; blk < 4; blk++) {
        int row = (int)(rnd() % (uint32_t)c->t->hoc_rows[blk]);
        b[5 + blk] = row - (row % c->t->hoc_step[blk]);
    }
}

/*
 * The decoder's residual less its mean over l = 1..L for b3..b8. D-STAR: its
 * decoder's reconstruction steps. Rate 33 (whose steps are private to its
 * decoder): the whole parameter decoder on a frame of L harmonics, with a
 * history of log2 amplitudes 0, where log2Ml = Tl - mean(Tl) + gamma - 0.5 log2 L.
 */
static int
decoder_residual(const struct codec* c, int L, const int b[9], double r[57]) {
    if (c->dstar) {
        int Ji[5];
        float Ri[9], Cik[5][18], Tl[57] = {0};
        block_lengths(c, L, Ji);
        mbe_ambe2400_reconstruct_prba(b[3], b[4], Ri);
        mbe_ambe2400_reconstruct_cik(Ri, b + 5, Ji, Cik);
        mbe_ambe2400_inverse_dct_tl(Cik, Ji, Tl);
        double mean = 0.0;
        for (int l = 1; l <= L; l++) {
            mean += Tl[l];
        }
        mean /= L;
        for (int l = 1; l <= L; l++) {
            r[l] = (double)Tl[l] - mean;
        }
        return 0;
    }
    int bits[9];
    char d[49];
    mbe_parms cur, prev, enhanced;
    memcpy(bits, b, sizeof(bits));
    bits[0] = c->b0_of_L[L];
    bits[1] = 0;
    bits[2] = 16;
    mbe_ambe2450_pack(bits, d);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int l = 0; l <= 56; l++) {
        prev.log2Ml[l] = 0.0f;
    }
    prev.gamma = 0.0f;
    if (mbe_decodeAmbe2450Parms(d, &cur, &prev) != MBE_AMBE2450_FRAME_VOICE || cur.L != L) {
        return 1;
    }
    const double offset = (double)cur.gamma - (0.5 * log2((double)L));
    for (int l = 1; l <= L; l++) {
        r[l] = (double)cur.log2Ml[l] - offset;
    }
    return 0;
}

/* The refinement's weight of harmonic l, as src/ambe/ambe_enc_search.c defines it
 * (its file comment and ambe_enc_weight()). */
static double
weight(const struct ambe_enc_frame* q, int l) {
    double w = pow(2.0, 0.5 * ((double)q->a[l] - (double)q->mean_a));
    return fmin(10.0, fmax(0.1, w));
}

/* E = sum_l w_l (That_l - mean(That) - (T_l - mean(T)))^2 for the rows b;
 * 1 if the rows do not decode or E is not finite. */
static int
weighted_error(const struct codec* c, const struct ambe_enc_frame* q, const int b[9], double* error) {
    double r[57];
    if (decoder_residual(c, q->L, b, r) != 0) {
        return 1;
    }
    double mean_t = 0.0;
    for (int l = 1; l <= q->L; l++) {
        mean_t += q->Tl[l];
    }
    mean_t /= q->L;
    *error = 0.0;
    for (int l = 1; l <= q->L; l++) {
        const double e = r[l] - ((double)q->Tl[l] - mean_t);
        *error += weight(q, l) * e * e;
    }
    return !finite_double(error);
}

/* The decoder receives every row but D-STAR's odd b8; a block of at most two
 * harmonics carries no higher-order vector and keeps index 0. */
static int
indices_receivable(const struct codec* c, const struct ambe_enc_frame* q) {
    if (q->b[3] < 0 || q->b[3] >= 512 || q->b[4] < 0 || q->b[4] >= 128) {
        return 0;
    }
    for (int blk = 0; blk < 4; blk++) {
        const int b = q->b[5 + blk];
        if (b < 0 || b >= c->t->hoc_rows[blk] || (b % c->t->hoc_step[blk]) != 0 || (q->Ji[blk + 1] <= 2 && b != 0)) {
            return 0;
        }
    }
    return 1;
}

/* Rate 33's pitch codes reach every L from 9 to 56. */
static int
map_pitch_codes(struct codec* c) {
    for (int L = 0; L <= 56; L++) {
        c->b0_of_L[L] = -1;
    }
    for (int b0 = 0; b0 < 120; b0++) {
        c->b0_of_L[mbe_ambe2450_harmonic_count(b0)] = b0;
    }
    for (int L = 9; L <= 56; L++) {
        if (c->b0_of_L[L] < 0) {
            printf("%s: no pitch code has %d harmonics\n", c->name, L);
            return 1;
        }
    }
    return 0;
}

/*
 * The model the refinement minimises is the decoder: for random rows and
 * every L from 9 to 56, ambe_enc_model_residual() of the rows' coefficients
 * reproduces the decoder's mean-removed residual within 2e-4 (a few units in
 * the last place of residuals up to about 10; measured worst 1.4e-6).
 */
static int
test_model_is_decoder(const struct codec* c) {
    double worst = 0.0;
    for (int L = 9; L <= 56; L++) {
        int Ji[5];
        block_lengths(c, L, Ji);
        for (int trial = 0; trial < 40; trial++) {
            int b[9];
            float x[AMBE_ENC_COEFFS], model[57];
            double r[57];
            random_indices(c, b);
            ambe_enc_coefficients(c->t, b, Ji, x);
            ambe_enc_model_residual(Ji, L, x, model);
            if (decoder_residual(c, L, b, r) != 0) {
                printf("%s model: L %d did not decode\n", c->name, L);
                return 1;
            }
            for (int l = 1; l <= L; l++) {
                /* fmax() would pass over a NaN difference. */
                const double difference = fabs((double)model[l] - r[l]);
                if (!finite_double(&difference)) {
                    printf("%s model: L %d, harmonic %d is not finite\n", c->name, L, l);
                    return 1;
                }
                worst = fmax(worst, difference);
            }
        }
    }
    printf("%s: the search's linear model is the decoder's residual for L 9..56 (worst %.2g)\n", c->name, worst);
    return worst > 2e-4;
}

/*
 * A frame to code: L harmonics, an envelope of a few formant-like peaks over
 * a tilt plus roughness, and a history of the same kind (or flat), through
 * the encoder's own prediction and nearest-row search.
 */
static void
random_frame(const struct codec* c, struct ambe_enc_frame* q) {
    mbe_parms prev, cur, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    memset(q, 0, sizeof(*q));
    q->cache = mbe_ambe2400_get_dct_cache();
    q->L = 9 + (int)(rnd() % 48u);
    prev.L = 9 + (int)(rnd() % 48u);
    const float base = uniform(-1.0f, 6.0f);
    const float tilt = uniform(-0.15f, 0.0f);
    float centre[3], depth[3];
    for (int p = 0; p < 3; p++) {
        centre[p] = uniform(0.0f, 1.0f);
        depth[p] = uniform(0.0f, 4.0f);
    }
    for (int l = 1; l <= 56; l++) {
        const float f = (float)l / (float)q->L;
        float a = base + (tilt * (float)l) + uniform(-0.7f, 0.7f);
        for (int p = 0; p < 3; p++) {
            const float u = (f - centre[p]) * 8.0f;
            a += depth[p] / (1.0f + (u * u));
        }
        q->a[l] = a;
        prev.log2Ml[l] = (rnd() & 1u) ? a + uniform(-1.0f, 1.0f) : 0.0f;
    }
    ambe_enc_mean_amplitude(q);
    ambe_enc_prediction(q, c->t, &prev);
    ambe_enc_quantize_nearest(q, c->t);
}

/*
 * The refinement never raises the weighted error of the decoder's residual
 * (allowing 1e-4 of it plus 1e-6 for float evaluation order), sends only
 * receivable rows, keeps 0 for blocks without a higher-order vector, and
 * does lower the error on most frames.
 */
static int
test_never_rises(const struct codec* c) {
    const int frames = 3000;
    int changed = 0;
    double open_sum = 0.0, refined_sum = 0.0;
    for (int f = 0; f < frames; f++) {
        struct ambe_enc_frame q;
        random_frame(c, &q);
        const struct ambe_enc_frame open = q;
        ambe_enc_refine(&q, c->t);
        double e_open, e_refined;
        if (weighted_error(c, &open, open.b, &e_open) != 0 || weighted_error(c, &q, q.b, &e_refined) != 0) {
            printf("%s: frame %d did not decode to a finite error\n", c->name, f);
            return 1;
        }
        if (e_refined > e_open + (1e-4 * e_open) + 1e-6) {
            printf("%s: frame %d (L %d): weighted error rose from %.9g to %.9g\n", c->name, f, q.L, e_open, e_refined);
            return 1;
        }
        if (!indices_receivable(c, &q) || q.b[0] != open.b[0] || q.b[1] != open.b[1] || q.b[2] != open.b[2]) {
            printf("%s: frame %d sent an index the decoder cannot receive, or changed b0..b2\n", c->name, f);
            return 1;
        }
        changed += memcmp(q.b, open.b, sizeof(q.b)) != 0;
        open_sum += e_open;
        refined_sum += e_refined;
    }
    printf("%s: refinement changed %d of %d frames, weighted error %.4g -> %.4g in all\n", c->name, changed, frames,
           open_sum, refined_sum);
    /* Most frames move (72-75% when measured) and the total error falls
     * (measured 5-7%; these random envelopes' roughness is mostly beyond any
     * row): the refinement is in effect. */
    return changed < frames / 2 || refined_sum > 0.98 * open_sum;
}

/*
 * A target the decoder reproduces exactly (random rows, the decoder's
 * residual, a zero history) keeps the nearest-row search's rows, which
 * reproduce it.
 */
static int
test_exact_target_kept(const struct codec* c) {
    mbe_parms prev, cur, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int l = 0; l <= 56; l++) {
        prev.log2Ml[l] = 0.0f;
    }
    for (int f = 0; f < 2000; f++) {
        struct ambe_enc_frame q;
        int b[9];
        double r[57];
        memset(&q, 0, sizeof(q));
        q.cache = mbe_ambe2400_get_dct_cache();
        q.L = 9 + (int)(rnd() % 48u);
        random_indices(c, b);
        if (decoder_residual(c, q.L, b, r) != 0) {
            return 1;
        }
        const float level = uniform(-2.0f, 6.0f);
        for (int l = 1; l <= q.L; l++) {
            q.a[l] = (float)r[l] + level;
        }
        ambe_enc_mean_amplitude(&q);
        ambe_enc_prediction(&q, c->t, &prev);
        ambe_enc_quantize_nearest(&q, c->t);
        const struct ambe_enc_frame open = q;
        double e_open;
        if (weighted_error(c, &open, open.b, &e_open) != 0) {
            return 1;
        }
        ambe_enc_refine(&q, c->t);
        if (e_open > 1e-6 || memcmp(q.b, open.b, sizeof(q.b)) != 0) {
            printf("%s: exactly reproducible frame %d (L %d, error %.3g) changed its rows\n", c->name, f, q.L, e_open);
            return 1;
        }
    }
    printf("%s: frames the decoder reproduces exactly keep their rows\n", c->name);
    return 0;
}

int
main(void) {
    struct codec codecs[2] = {
        {"D-STAR", NULL, 1, {0}},
        {"AMBE+2", NULL, 0, {0}},
    };
    codecs[0].t = mbe_ambe2400_enc_tables();
    codecs[1].t = mbe_ambe2450_enc_tables();
    int fails = map_pitch_codes(&codecs[1]);
    for (int i = 0; i < 2; i++) {
        fails += test_model_is_decoder(&codecs[i]);
        fails += test_never_rises(&codecs[i]);
        fails += test_exact_target_kept(&codecs[i]);
    }
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}
