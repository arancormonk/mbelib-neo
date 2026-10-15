// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Signal-processing stages of the TETRA ACELP decoder (EN 300 395-2 clause 4.2.3).
 *
 * Integer fixed-point arithmetic through mbe_fixed_point.h and equivalent SIMD kernels,
 * with identical output on every platform and build configuration.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "mbe_compiler.h"
#include "mbe_fixed_point.h"
#include "tetra_acelp.h"
#include "tetra_acelp_codebooks.h"
#include "tetra_acelp_constants.h"
#include "tetra_simd.h" // IWYU pragma: keep (TETRA_SIMD selects the integer vector kernels)

/* Distance the LSPs on either side of a split-VQ boundary are widened towards (Q15 cosine). */
#define TETRA_LSP_GAP_34       916
#define TETRA_LSP_GAP_67       1244
#define TETRA_SHARP            26216 /* 0.8 pitch sharpening factor (Q13 6554, scaled to Q15) */
#define TETRA_PULSE1_AMP       23168 /* 1.4142 pulse amplitude (Q12 5792, scaled to Q14) */
/* Q8 log2 offset of both gains: the Q formats of the energies, of their mantissas and of the Q12
 * gains, resolved against the reference decoder. */
#define TETRA_GAIN_OFFSET      4434
#define TETRA_GAIN_P_MAX       4915 /* 1.2 in Q12 */
/* Both interpolation filters' first 18 taps have an absolute sum <= 65535.
 * Their MAC prefixes cannot overflow int32 for any int16 samples. Later taps
 * can overflow and must saturate in order. test_tetra_simd checks this bound. */
#define TETRA_INTERP_SAFE_TAPS 18

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
        i = 0;
#if defined(TETRA_SIMD)
        for (; i + 4 <= TETRA_M; i += 4) {
            tetra_v32 acc =
                tetra_vadd(tetra_vmul(old_lsp + i, w_old[k]), tetra_vmul(new_lsp + i, (int16_t)(32768 - w_old[k])));
            tetra_vstore_sat(out[k] + i, tetra_vshr(acc, 15));
        }
#endif
        for (; i < TETRA_M; i++) {
            int32_t acc = mbe_fx_L_mult(old_lsp[i], w_old[k]);
            acc = mbe_fx_L_mac(acc, new_lsp[i], (int16_t)(32768 - w_old[k]));
            out[k][i] = mbe_fx_extract_h(acc);
        }
    }
    memcpy(out[3], new_lsp, sizeof(out[3]));
}

/**
 * (hi, lo) * q, the 32 x 16-bit product of the polynomial recursion. Unlike the usual
 * double-precision product 2 hi q + 2 floor(lo q / 2^15), the product of lo's top bit (bit 14) and q
 * is kept exactly: 2 hi q + b q + 2 floor(l q / 2^15) with b = bit 14 of lo and l its lower 14 bits.
 * For an even q the two agree; for an odd q they can differ by one in the last place (resolved
 * against the reference decoder).
 */
static int32_t
lsp_mpy(int16_t hi, int16_t lo, int16_t q) {
    int64_t top = (lo >> 14) & 1;
    int64_t low = lo & 0x3FFF;
    return mbe_fx_sat32((int64_t)hi * q * 2 + top * q + 2 * mbe_fx_asr64(low * q, 15));
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
            t = mbe_fx_L_shl(lsp_mpy(hi, lo, q), 1);
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

static void
syn_filt_scalar(const int16_t a[TETRA_M + 1], const int16_t* x, int16_t* y, int len, int16_t* buf) {
    int n;
    int i;
    for (n = 0; n < len; n++) {
        int32_t acc = mbe_fx_L_shl(x[n], 12);
        for (i = 1; i <= TETRA_M; i++) {
            acc = mbe_fx_L_msu0(acc, a[i], buf[TETRA_M + n - i]);
        }
        acc = mbe_fx_L_add(mbe_fx_L_shl(acc, 4), 0x8000);
        buf[TETRA_M + n] = mbe_fx_extract_h(acc);
        y[n] = buf[TETRA_M + n];
    }
}

#if defined(TETRA_SIMD)
/* Reordering a saturating reduction is only valid when no intermediate can overflow.
 * The bound also covers the initial x * 4096, for every possible int16 input/history. */
static int
syn_filt_bounded(const int16_t a[TETRA_M + 1], int16_t reversed[TETRA_M]) {
    int32_t sum = 4096;
    for (int i = 0; i < TETRA_M; ++i) {
        int32_t v = a[TETRA_M - i];
        reversed[i] = (int16_t)v;
        sum += v < 0 ? -v : v;
    }
    return sum <= 65535;
}

static void
syn_filt_simd(const int16_t reversed[TETRA_M], const int16_t* x, int16_t* y, int len, int16_t* buf) {
    for (int n = 0; n < len; ++n) {
        int32_t acc = (int32_t)x[n] * 4096 - tetra_dot10(reversed, buf + n);
        /* sat16((acc + 2048) >> 12) is exactly the original saturated shift/round. */
        buf[TETRA_M + n] = mbe_fx_sat16(mbe_fx_asr64((int64_t)acc + 2048, 12));
        y[n] = buf[TETRA_M + n];
    }
}
#endif

void
tetra_syn_filt(const int16_t a[TETRA_M + 1], const int16_t* x, int16_t* y, int len, int16_t mem[TETRA_M], int update) {
    int16_t buf[TETRA_M + TETRA_FRAME];
    memcpy(buf, mem, TETRA_M * sizeof(int16_t));
#if defined(TETRA_SIMD)
    int16_t reversed[TETRA_M];
    if (syn_filt_bounded(a, reversed)) {
        syn_filt_simd(reversed, x, y, len, buf);
    } else
#endif
    {
        syn_filt_scalar(a, x, y, len, buf);
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
    i = 1;
#if defined(TETRA_SIMD)
    for (; i + 3 <= TETRA_M; i += 4) {
        tetra_v32 p = tetra_vmulpair(a + i, gamma_pow + i - 1);
        tetra_vstore_sat(ap + i, tetra_vshr(tetra_vadd(p, tetra_vset(0x4000)), 15));
    }
#endif
    for (; i <= TETRA_M; i++) {
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

#if !defined(TETRA_SIMD)
/** f(n - m), zero outside the subframe (a pulse at 60 or beyond is absent). */
static int64_t
shifted(const int16_t f[TETRA_SUBFRAME], int n, int m) {
    return (n >= m && m < TETRA_SUBFRAME) ? f[n - m] : 0;
}
#endif

void
tetra_codevector(const int16_t a[TETRA_M + 1], int t0, const tetra_subframe_params* sp, int16_t code[TETRA_SUBFRAME]) {
    /* A pulse position is at most 63. Zero padding makes every four-lane load
     * valid, including absent pulses and the first partial pulse block. */
#if defined(TETRA_SIMD)
    int16_t padded[64 + TETRA_SUBFRAME] = {0};
    int16_t* f = padded + 64;
#else
    int16_t f[TETRA_SUBFRAME];
#endif
    int pos[4];
    int n;
    shaping_response(a, t0, f);
    tetra_pulse_positions(sp, pos);
#if defined(TETRA_SIMD)
    for (n = 0; n < TETRA_SUBFRAME; n += 4) {
        tetra_v32 acc = tetra_vmul(f + n - pos[0], TETRA_PULSE1_AMP);
        tetra_v32 rest = tetra_vsub(tetra_vload(f + n - pos[2]), tetra_vload(f + n - pos[1]));
        rest = tetra_vsub(rest, tetra_vload(f + n - pos[3]));
        /* Modulo-2^32 adds/shifts preserve the low 16 bits after the Q14 shift,
         * even when the full pulse sum exceeds int32. No saturating narrowing. */
        acc = tetra_vadd(acc, tetra_vshl(rest, 14));
        if (sp->sign) {
            acc = tetra_vsub(tetra_vset(0), acc);
        }
        tetra_vstore_wrap(code + n, tetra_vshr(acc, 14));
    }
#else
    for (n = 0; n < TETRA_SUBFRAME; n++) {
        int64_t acc = TETRA_PULSE1_AMP * shifted(f, n, pos[0]);
        acc += mbe_fx_lsl64(-shifted(f, n, pos[1]) + shifted(f, n, pos[2]) - shifted(f, n, pos[3]), 14);
        if (sp->sign) {
            acc = -acc;
        }
        /* The reference decoder keeps the low 16 bits; only an extreme A(z) gets past them. */
        code[n] = mbe_fx_extract_l((int32_t)mbe_fx_asr64(acc, 14));
    }
#endif
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

#if defined(TETRA_SIMD) && defined(__AVX2__) && defined(MBE_SIMD_TARGET_SSE2)
static int
adaptive_avx2(int16_t* exc_cur, int t0, int first, const int16_t* taps) {
    int n = 0;
    if (t0 - (first + 31) < 8) {
        return 0; /* Short lags need the four-lane kernel's earlier stores. */
    }
    for (; n + 8 <= TETRA_SUBFRAME; n += 8) {
        __m256i acc = _mm256_setzero_si256();
        int k;
        for (k = 0; k < TETRA_INTERP_SAFE_TAPS; ++k) {
            __m128i x = _mm_loadu_si128((const __m128i*)(exc_cur + n - t0 + first + k));
            __m256i p = _mm256_mullo_epi32(_mm256_cvtepi16_epi32(x), _mm256_set1_epi32(taps[k]));
            acc = _mm256_add_epi32(acc, p);
        }
        for (; k < 32; ++k) {
            __m128i x = _mm_loadu_si128((const __m128i*)(exc_cur + n - t0 + first + k));
            __m256i p = _mm256_mullo_epi32(_mm256_cvtepi16_epi32(x), _mm256_set1_epi32(taps[k]));
            acc = tetra_vadds8(acc, p);
        }
        acc = tetra_vadds8(acc, acc);
        acc = _mm256_srai_epi32(tetra_vadds8(acc, _mm256_set1_epi32(0x8000)), 16);
        _mm_storeu_si128((__m128i*)(exc_cur + n),
                         _mm_packs_epi32(_mm256_castsi256_si128(acc), _mm256_extracti128_si256(acc, 1)));
    }
    return n;
}
#endif

void
tetra_adaptive_vector(int16_t* exc_cur, int t0, int frac) {
    const int16_t* taps = (frac > 0) ? tetra_interp_p13 : tetra_interp_m13;
    int first = (frac > 0) ? -16 : -15;
    int n;
    int k;
    n = 0;
#if defined(TETRA_SIMD) && defined(__AVX2__) && defined(MBE_SIMD_TARGET_SSE2)
    if (frac != 0) {
        n = adaptive_avx2(exc_cur, t0, first, taps);
    }
#endif
#if defined(TETRA_SIMD)
    /* Four output lanes, taps in their original order, including saturation after
     * every MAC. The newest tap must precede the block, even at the shortest lag. */
    if (frac != 0 && t0 - (first + 31) >= 4) {
        for (; n + 4 <= TETRA_SUBFRAME; n += 4) {
            tetra_v32 acc = tetra_vset(0);
            for (k = 0; k < TETRA_INTERP_SAFE_TAPS; ++k) {
                acc = tetra_vadd(acc, tetra_vmul(exc_cur + n - t0 + first + k, taps[k]));
            }
            for (; k < 32; ++k) {
                acc = tetra_vadds(acc, tetra_vmul(exc_cur + n - t0 + first + k, taps[k]));
            }
            acc = tetra_vadds(acc, acc);
            acc = tetra_vadds(acc, tetra_vset(0x8000));
            tetra_vstore_sat(exc_cur + n, tetra_vshr(acc, 16));
        }
    }
#endif
    for (; n < TETRA_SUBFRAME; n++) {
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

/** Gain (Q12, not yet narrowed, never negative) from a quantized energy (Q8) and the log2 (Q15) of
 *  the energy product of equation 31/33. */
static int32_t
gain_from_log(int16_t energy_q8, int32_t log_q15) {
    /* The log is truncated to Q8 before the subtraction (equations 38-39); z is then log2 of the
     * squared gain in Q8, i.e. log2 of the gain in Q9. */
    int32_t z = energy_q8 + TETRA_GAIN_OFFSET - (int32_t)mbe_fx_asr64(log_q15, 7);
    int32_t e = (int32_t)mbe_fx_asr64(z, 9);
    int32_t f = z - (int32_t)mbe_fx_lsl64(e, 9);
    return tetra_pow2((int16_t)e, (int16_t)mbe_fx_lsl64(f, 6));
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
    int32_t p;
    /* One log2 per energy product of 16-bit mantissas; the two high halves drop 32 binary places,
     * which TETRA_GAIN_OFFSET accounts for. */
    p = gain_from_log(energy_p, log2_q15(mbe_fx_L_mult0(mv, mpg)) - (int32_t)mbe_fx_lsl64(nv + npg, 15));
    *gp = (int16_t)(p < TETRA_GAIN_P_MAX ? p : TETRA_GAIN_P_MAX);
    /* The code gain is not capped: as in the reference decoder it keeps its low 16 bits, which a
     * gain of 8.0 or more (only from extreme energies) gets past. */
    *gc = mbe_fx_extract_l(gain_from_log(energy_c, log2_q15(mbe_fx_L_mult0(mc, mpg)) - (int32_t)mbe_fx_lsl64(npg, 15)));
}

int16_t
tetra_excitation(int16_t gp, int16_t v, int16_t gc, int16_t c) {
    /* Kept to the low 16 bits as the reference decoder does; the 32-bit sum cannot overflow, gp
     * being at most 4915. */
    return mbe_fx_extract_l(mbe_fx_L_shr_r(mbe_fx_L_add(mbe_fx_L_mult0(v, gp), mbe_fx_L_mult0(c, gc)), 12));
}

void
tetra_mix_excitation(int16_t gp, int16_t gc, int16_t v[TETRA_SUBFRAME], const int16_t c[TETRA_SUBFRAME]) {
    int n = 0;
#if defined(TETRA_SIMD)
    for (; n + 4 <= TETRA_SUBFRAME; n += 4) {
        /* gp is capped at 4915; the sum including the rounding bias fits int32. */
        tetra_v32 acc = tetra_vadd(tetra_vmul(v + n, gp), tetra_vmul(c + n, gc));
        acc = tetra_vadd(acc, tetra_vset(2048));
        tetra_vstore_wrap(v + n, tetra_vshr(acc, 12));
    }
#endif
    for (; n < TETRA_SUBFRAME; ++n) {
        v[n] = tetra_excitation(gp, v[n], gc, c[n]);
    }
}

int32_t
tetra_energy(const int16_t* x, int len, int32_t init) {
    /* Squares are nonnegative, so saturating once at the end is identical to
     * saturating after every MAC. Widen before summing: two -32768 squares
     * produce unsigned 0x80000000, and 60 full-scale samples need 36 bits. */
    int n = 0;
    int64_t sum = init;
#if defined(TETRA_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
    __m128i acc = _mm_setzero_si128();
    for (; n + 8 <= len; n += 8) {
        __m128i v = _mm_loadu_si128((const __m128i*)(x + n));
        __m128i p = _mm_madd_epi16(v, v);
        acc = _mm_add_epi64(acc, _mm_unpacklo_epi32(p, _mm_setzero_si128()));
        acc = _mm_add_epi64(acc, _mm_unpackhi_epi32(p, _mm_setzero_si128()));
    }
    int64_t lanes[2];
    _mm_storeu_si128((__m128i*)lanes, acc);
    sum += lanes[0] + lanes[1];
#elif defined(TETRA_SIMD) && defined(MBE_SIMD_TARGET_NEON)
    uint64x2_t acc = vdupq_n_u64(0);
    for (; n + 4 <= len; n += 4) {
        int16x4_t v = vld1_s16(x + n);
        acc = vpadalq_u32(acc, vreinterpretq_u32_s32(vmull_s16(v, v)));
    }
    sum += (int64_t)(vgetq_lane_u64(acc, 0) + vgetq_lane_u64(acc, 1));
#endif
    for (; n < len; ++n) {
        sum += (int64_t)((int32_t)x[n] * x[n]);
    }
    return mbe_fx_sat32(sum);
}
