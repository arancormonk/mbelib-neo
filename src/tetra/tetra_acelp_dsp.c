// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Signal-processing stages of the TETRA ACELP decoder (EN 300 395-2 clause 4.2.3).
 *
 * Everything is integer fixed-point arithmetic through mbe_fixed_point.h, so the output is
 * identical on every platform and build configuration.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "mbe_fixed_point.h"
#include "tetra_acelp.h"
#include "tetra_acelp_codebooks.h"
#include "tetra_acelp_constants.h"

/* Distance the LSPs on either side of a split-VQ boundary are widened towards (Q15 cosine). */
#define TETRA_LSP_GAP_34  916
#define TETRA_LSP_GAP_67  1244
#define TETRA_SHARP       26216 /* 0.8 pitch sharpening factor (Q13 6554, scaled to Q15) */
#define TETRA_PULSE1_AMP  23168 /* 1.4142 pulse amplitude (Q12 5792, scaled to Q14) */
/* Q8 log2 offset of both gains: the Q formats of the energies, of their mantissas and of the Q12
 * gains, resolved against the reference decoder. */
#define TETRA_GAIN_OFFSET 4434
#define TETRA_GAIN_P_MAX  4915 /* 1.2 in Q12 */

/** Widen the LSP pair (i, i + 1) towards gap, moving each by half the shortfall rounded down (so the
 * gap ends at gap, or gap - 1 for an odd shortfall). */
static void
lsp_widen(int16_t lsp[TETRA_M], int i, int gap) {
    int32_t d = (int32_t)lsp[i] - lsp[i + 1];
    if (d < gap) {
        int32_t t = (int32_t)mbe_fx_asr64((int64_t)gap - d, 1);
        lsp[i] = mbe_fx_sat16((int64_t)lsp[i] + t);
        lsp[i + 1] = mbe_fx_sat16((int64_t)lsp[i + 1] - t);
    }
}

void
tetra_lsp_decode(const int lsp_index[3], int16_t lsp[TETRA_M]) {
    int i;
    for (i = 0; i < 3; i++) {
        lsp[i] = tetra_lsp_cb1[lsp_index[0]][i];
        lsp[3 + i] = tetra_lsp_cb2[lsp_index[1]][i];
    }
    for (i = 0; i < 4; i++) {
        lsp[6 + i] = tetra_lsp_cb3[lsp_index[2]][i];
    }
    /* The three sub-vectors are quantized separately, so only the LSPs on either side of a split
     * boundary can come out too close or crossed. */
    lsp_widen(lsp, 2, TETRA_LSP_GAP_34);
    lsp_widen(lsp, 5, TETRA_LSP_GAP_67);
}

int
tetra_lsp_ordered(const int16_t lsp[TETRA_M]) {
    int i;
    for (i = 0; i + 1 < TETRA_M; i++) {
        if (lsp[i] <= lsp[i + 1]) {
            return 0;
        }
    }
    return 1;
}

void
tetra_lsp_interpolate(const int16_t old_lsp[TETRA_M], const int16_t new_lsp[TETRA_M],
                      int16_t out[TETRA_SUBFRAMES][TETRA_M]) {
    static const int16_t w_old[3] = {24576, 16384, 8192};
    int k;
    int i;
    for (k = 0; k < 3; k++) {
        for (i = 0; i < TETRA_M; i++) {
            int32_t acc = mbe_fx_L_mult(old_lsp[i], w_old[k]);
            acc = mbe_fx_L_mac(acc, new_lsp[i], (int16_t)(32768 - w_old[k]));
            out[k][i] = mbe_fx_extract_h(acc);
        }
    }
    memcpy(out[3], new_lsp, sizeof(out[3]));
}

/** Coefficients f[0..5] (Q24) of the sum or difference polynomial of five LSPs (stride 2). */
static void
lsp_polynomial(const int16_t* lsp, int32_t f[6]) {
    int i;
    int j;
    f[0] = mbe_fx_L_mult(4096, 2048);
    f[1] = mbe_fx_L_msu(0, lsp[0], 512);
    for (i = 2; i <= 5; i++) {
        int16_t q = lsp[(ptrdiff_t)2 * (i - 1)];
        f[i] = f[i - 2];
        for (j = i; j > 1; j--) {
            int16_t hi;
            int16_t lo;
            int32_t t;
            mbe_fx_L_extract(f[j - 1], &hi, &lo);
            t = mbe_fx_L_shl(mbe_fx_mpy_32_16(hi, lo, q), 1);
            f[j] = mbe_fx_L_sub(mbe_fx_L_add(f[j], f[j - 2]), t);
        }
        f[1] = mbe_fx_L_msu(f[1], q, 512);
    }
}

void
tetra_lsp_to_lpc(const int16_t lsp[TETRA_M], int16_t a[TETRA_M + 1]) {
    int32_t f1[6];
    int32_t f2[6];
    int i;
    lsp_polynomial(&lsp[0], f1);
    lsp_polynomial(&lsp[1], f2);
    for (i = 5; i > 0; i--) {
        f1[i] = mbe_fx_L_add(f1[i], f1[i - 1]);
        f2[i] = mbe_fx_L_sub(f2[i], f2[i - 1]);
    }
    a[0] = 4096;
    for (i = 1; i <= 5; i++) {
        a[i] = mbe_fx_extract_l(mbe_fx_L_shr_r(mbe_fx_L_add(f1[i], f2[i]), 13));
        a[TETRA_M + 1 - i] = mbe_fx_extract_l(mbe_fx_L_shr_r(mbe_fx_L_sub(f1[i], f2[i]), 13));
    }
}

void
tetra_syn_filt(const int16_t a[TETRA_M + 1], const int16_t* x, int16_t* y, int len, int16_t mem[TETRA_M], int update) {
    int16_t buf[TETRA_M + TETRA_FRAME];
    int n;
    int i;
    memcpy(buf, mem, TETRA_M * sizeof(int16_t));
    for (n = 0; n < len; n++) {
        int32_t acc = mbe_fx_L_shl(x[n], 12);
        for (i = 1; i <= TETRA_M; i++) {
            acc = mbe_fx_L_msu0(acc, a[i], buf[TETRA_M + n - i]);
        }
        acc = mbe_fx_L_add(mbe_fx_L_shl(acc, 4), 0x8000);
        buf[TETRA_M + n] = mbe_fx_extract_h(acc);
        y[n] = buf[TETRA_M + n];
    }
    if (update) {
        memcpy(mem, &buf[len], TETRA_M * sizeof(int16_t));
    }
}

/** A(z/gamma) with gamma^i given in Q15. */
static void
weight_lpc(const int16_t a[TETRA_M + 1], const int16_t* gamma_pow, int16_t ap[TETRA_M + 1]) {
    int i;
    ap[0] = a[0];
    for (i = 1; i <= TETRA_M; i++) {
        ap[i] = mbe_fx_round(mbe_fx_L_mult(a[i], gamma_pow[i - 1]));
    }
}

/** Impulse response of F(z) = A(z/0.75) / A(z/0.85), sharpened when t0 < 60 (Q12). */
static void
shaping_response(const int16_t a[TETRA_M + 1], int t0, int16_t f[TETRA_SUBFRAME]) {
    int16_t ap1[TETRA_M + 1];
    int16_t ap2[TETRA_M + 1];
    int16_t x[TETRA_SUBFRAME];
    int16_t mem[TETRA_M];
    weight_lpc(a, tetra_gamma1_pow, ap1);
    weight_lpc(a, tetra_gamma2_pow, ap2);
    memset(x, 0, sizeof(x));
    memcpy(x, ap1, sizeof(ap1));
    memset(mem, 0, sizeof(mem));
    tetra_syn_filt(ap2, x, f, TETRA_SUBFRAME, mem, 0);
    if (t0 < TETRA_SUBFRAME) {
        int n;
        for (n = t0; n < TETRA_SUBFRAME; n++) {
            f[n] = mbe_fx_add(f[n], mbe_fx_mult(f[n - t0], TETRA_SHARP));
        }
    }
}

void
tetra_pulse_positions(const tetra_subframe_params* sp, int pos[4]) {
    pos[0] = 2 * sp->pulse1 + sp->shift;
    pos[1] = 2 + 8 * sp->pulse2 + sp->shift;
    pos[2] = 4 + 8 * sp->pulse3 + sp->shift;
    pos[3] = 6 + 8 * sp->pulse4 + sp->shift;
}

/** f(n - m), zero outside the subframe (a pulse at 60 or beyond is absent). */
static int64_t
shifted(const int16_t f[TETRA_SUBFRAME], int n, int m) {
    return (n >= m && m < TETRA_SUBFRAME) ? f[n - m] : 0;
}

void
tetra_codevector(const int16_t a[TETRA_M + 1], int t0, const tetra_subframe_params* sp, int16_t code[TETRA_SUBFRAME]) {
    int16_t f[TETRA_SUBFRAME];
    int pos[4];
    int n;
    shaping_response(a, t0, f);
    tetra_pulse_positions(sp, pos);
    for (n = 0; n < TETRA_SUBFRAME; n++) {
        int64_t acc = TETRA_PULSE1_AMP * shifted(f, n, pos[0]);
        acc += mbe_fx_lsl64(-shifted(f, n, pos[1]) + shifted(f, n, pos[2]) - shifted(f, n, pos[3]), 14);
        if (sp->sign) {
            acc = -acc;
        }
        code[n] = mbe_fx_sat16(mbe_fx_asr64(acc, 14));
    }
}

void
tetra_decode_pitch(int subframe, int index, int t1, int* t0, int* frac) {
    if (subframe == 0) {
        if (index < 197) {
            *t0 = (index + 2) / 3 + 19;
            *frac = index - 3 * *t0 + 58;
        } else {
            *t0 = index - 112;
            *frac = 0;
        }
    } else {
        int tmin = t1 - 5;
        int tmax;
        int i = (index + 2) / 3 - 1;
        if (tmin < TETRA_PIT_MIN) {
            tmin = TETRA_PIT_MIN;
        }
        tmax = tmin + 9;
        if (tmax > TETRA_PIT_MAX) {
            tmax = TETRA_PIT_MAX;
            tmin = tmax - 9;
        }
        *t0 = tmin + i;
        *frac = index - 2 - 3 * i;
    }
}

void
tetra_adaptive_vector(int16_t* exc_cur, int t0, int frac) {
    const int16_t* taps = (frac > 0) ? tetra_interp_p13 : tetra_interp_m13;
    int first = (frac > 0) ? -16 : -15;
    int n;
    int k;
    for (n = 0; n < TETRA_SUBFRAME; n++) {
        const int16_t* x = exc_cur + n - t0;
        if (frac == 0) {
            exc_cur[n] = x[0];
            continue;
        }
        {
            int32_t acc = 0;
            for (k = 0; k < 32; k++) {
                acc = mbe_fx_L_mac0(acc, x[first + k], taps[k]);
            }
            exc_cur[n] = mbe_fx_round(mbe_fx_L_add(acc, acc));
        }
    }
}

void
tetra_log2(int32_t x, int16_t* exponent, int16_t* fraction) {
    int e;
    int i;
    int16_t a;
    int32_t y;
    if (x <= 0) {
        *exponent = 0;
        *fraction = 0;
        return;
    }
    e = mbe_fx_norm_l(x);
    x = mbe_fx_L_shl(x, e);         /* now in [2^30, 2^31) */
    i = mbe_fx_L_shr(x, 25) & 0x1F; /* bits 29..25: the table interval */
    a = (int16_t)(mbe_fx_extract_l(mbe_fx_L_shr(x, 10)) & 0x7FFF);
    y = mbe_fx_L_msu(mbe_fx_L_deposit_h(tetra_log2_tab[i]), mbe_fx_sub(tetra_log2_tab[i], tetra_log2_tab[i + 1]), a);
    *exponent = (int16_t)(30 - e);
    *fraction = mbe_fx_extract_h(y);
}

int32_t
tetra_pow2(int16_t exponent, int16_t fraction) {
    /* fraction is Q15 in [0, 1): bits 14..10 pick the table interval, bits 9..0 interpolate. */
    int i = (fraction >> 10) & 0x1F;
    int16_t a = (int16_t)((fraction & 0x3FF) << 5);
    int32_t x = mbe_fx_L_deposit_h(tetra_pow2_tab[i]);
    x = mbe_fx_L_msu(x, mbe_fx_sub(tetra_pow2_tab[i], tetra_pow2_tab[i + 1]), a);
    return mbe_fx_L_shr_r(x, 30 - exponent);
}

/** log2(x) in Q15 (exponent << 15 plus the fraction). */
static int32_t
log2_q15(int32_t x) {
    int16_t e;
    int16_t f;
    tetra_log2(x, &e, &f);
    return (int32_t)e * 32768 + f;
}

/** High 16 bits of x after normalization; *shift receives the left shift that normalized it. */
static int16_t
mantissa(int32_t x, int* shift) {
    *shift = mbe_fx_norm_l(x);
    return mbe_fx_extract_h(mbe_fx_L_shl(x, *shift));
}

/** Gain (Q12) from a quantized energy (Q8) and the log2 (Q15) of the energy product of equation 31/33. */
static int16_t
gain_from_log(int16_t energy_q8, int32_t log_q15) {
    /* The log is truncated to Q8 before the subtraction (equations 38-39); z is then log2 of the
     * squared gain in Q8, i.e. log2 of the gain in Q9. */
    int32_t z = energy_q8 + TETRA_GAIN_OFFSET - (int32_t)mbe_fx_asr64(log_q15, 7);
    int32_t e = (int32_t)mbe_fx_asr64(z, 9);
    int32_t f = z - (int32_t)mbe_fx_lsl64(e, 9);
    return mbe_fx_sat16(tetra_pow2((int16_t)e, (int16_t)mbe_fx_lsl64(f, 6)));
}

void
tetra_gains(int16_t energy_p, int16_t energy_c, int32_t v_energy, int32_t c_energy, int32_t pg_energy, int16_t* gp,
            int16_t* gc) {
    int npg;
    int nv;
    int16_t mpg = mantissa(pg_energy, &npg);
    int16_t mv = mantissa(v_energy, &nv);
    /* The code energy keeps a fixed scale: its high half, not normalized (at least 43, see the tests). */
    int16_t mc = mbe_fx_extract_h(c_energy);
    /* One log2 per energy product of 16-bit mantissas; the two high halves drop 32 binary places,
     * which TETRA_GAIN_OFFSET accounts for. */
    *gp = gain_from_log(energy_p, log2_q15(mbe_fx_L_mult0(mv, mpg)) - (int32_t)mbe_fx_lsl64(nv + npg, 15));
    if (*gp > TETRA_GAIN_P_MAX) {
        *gp = TETRA_GAIN_P_MAX;
    }
    *gc = gain_from_log(energy_c, log2_q15(mbe_fx_L_mult0(mc, mpg)) - (int32_t)mbe_fx_lsl64(npg, 15));
}
