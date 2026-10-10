// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Double-precision, bounded SIMD tone projection kernels. */
#include <math.h>
#include "mbe_simd.h"
#include "mbe_tone_fit.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Solve the symmetric positive definite system a x = b of order n <= 4 by
 * Cholesky factorization; returns 0 when a is not positive definite. */
static int
tone_solve(double a[4][4], const double b[4], int n, double x[4]) {
    double l[4][4] = {{0}};
    double y[4] = {0};
    if (n < 1 || n > 4) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = a[i][j];
            for (int k = 0; k < j; k++) {
                sum -= l[i][k] * l[j][k];
            }
            if (i == j) {
                if (sum <= 0.0) {
                    return 0;
                }
                l[i][i] = sqrt(sum);
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int k = 0; k < i; k++) {
            sum -= l[i][k] * y[k];
        }
        y[i] = sum / l[i][i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double sum = y[i];
        for (int k = i + 1; k < n; k++) {
            sum -= l[k][i] * x[k];
        }
        x[i] = sum / l[i][i];
    }
    return 1;
}

/* Each lane follows one residue class of sample indices. Seed consecutive
 * samples with the scalar recurrence, then rotate by WIDTH samples at a time.
 * All oscillator and projection arithmetic stays in double precision. */
struct tone_oscillator {
    mbe_vd c, s, cw, sw;
};

static struct tone_oscillator
tone_oscillator_init(double hz) {
    double w = 2.0 * M_PI * hz / 8000.0;
    double cw = cos(w), sw = sin(w);
    double c[MBE_VD_WIDTH], s[MBE_VD_WIDTH];
    double cc = 1.0, ss = 0.0;
    for (int i = 0; i < MBE_VD_WIDTH; ++i) {
        c[i] = cc;
        s[i] = ss;
        double next = cc * cw - ss * sw;
        ss = ss * cw + cc * sw;
        cc = next;
    }
    struct tone_oscillator o = {mbe_vd_load(c), mbe_vd_load(s), mbe_vd_set(cc), mbe_vd_set(ss)};
    return o;
}

static void
tone_oscillator_step(struct tone_oscillator* o) {
    mbe_vd next = mbe_vd_sub(mbe_vd_mul(o->c, o->cw), mbe_vd_mul(o->s, o->sw));
    o->s = mbe_vd_add(mbe_vd_mul(o->s, o->cw), mbe_vd_mul(o->c, o->sw));
    o->c = next;
}

static void
tone_accumulate_one(const float* x, int n, double hz, double gram[4][4], double rhs[4]) {
    struct tone_oscillator o = tone_oscillator_init(hz);
    mbe_vd cc = mbe_vd_set(0.0), cs = cc, ss = cc, cx = cc, sx = cc;
    int i = 0;
    for (; i + MBE_VD_WIDTH <= n; i += MBE_VD_WIDTH) {
        mbe_vd v = mbe_vd_load_float(x + i);
        cc = mbe_vd_add(cc, mbe_vd_mul(o.c, o.c));
        cs = mbe_vd_add(cs, mbe_vd_mul(o.c, o.s));
        ss = mbe_vd_add(ss, mbe_vd_mul(o.s, o.s));
        cx = mbe_vd_add(cx, mbe_vd_mul(o.c, v));
        sx = mbe_vd_add(sx, mbe_vd_mul(o.s, v));
        tone_oscillator_step(&o);
    }
    gram[0][0] = mbe_vd_sum(cc);
    gram[1][0] = mbe_vd_sum(cs);
    gram[1][1] = mbe_vd_sum(ss);
    rhs[0] = mbe_vd_sum(cx);
    rhs[1] = mbe_vd_sum(sx);
    double c[MBE_VD_WIDTH], s[MBE_VD_WIDTH];
    mbe_vd_store(c, o.c);
    mbe_vd_store(s, o.s);
    for (int k = 0; i < n; ++i, ++k) {
        gram[0][0] += c[k] * c[k];
        gram[1][0] += c[k] * s[k];
        gram[1][1] += s[k] * s[k];
        rhs[0] += c[k] * (double)x[i];
        rhs[1] += s[k] * (double)x[i];
    }
}

static void
tone_accumulate_two(const float* x, int n, const double hz[2], double gram[4][4], double rhs[4]) {
    struct tone_oscillator a = tone_oscillator_init(hz[0]), b = tone_oscillator_init(hz[1]);
    mbe_vd g00 = mbe_vd_set(0.0), g10 = g00, g11 = g00, g20 = g00, g21 = g00, g22 = g00;
    mbe_vd g30 = g00, g31 = g00, g32 = g00, g33 = g00, r0 = g00, r1 = g00, r2 = g00, r3 = g00;
    int i = 0;
    for (; i + MBE_VD_WIDTH <= n; i += MBE_VD_WIDTH) {
        mbe_vd v = mbe_vd_load_float(x + i);
        g00 = mbe_vd_add(g00, mbe_vd_mul(a.c, a.c));
        g10 = mbe_vd_add(g10, mbe_vd_mul(a.s, a.c));
        g11 = mbe_vd_add(g11, mbe_vd_mul(a.s, a.s));
        g20 = mbe_vd_add(g20, mbe_vd_mul(b.c, a.c));
        g21 = mbe_vd_add(g21, mbe_vd_mul(b.c, a.s));
        g22 = mbe_vd_add(g22, mbe_vd_mul(b.c, b.c));
        g30 = mbe_vd_add(g30, mbe_vd_mul(b.s, a.c));
        g31 = mbe_vd_add(g31, mbe_vd_mul(b.s, a.s));
        g32 = mbe_vd_add(g32, mbe_vd_mul(b.s, b.c));
        g33 = mbe_vd_add(g33, mbe_vd_mul(b.s, b.s));
        r0 = mbe_vd_add(r0, mbe_vd_mul(a.c, v));
        r1 = mbe_vd_add(r1, mbe_vd_mul(a.s, v));
        r2 = mbe_vd_add(r2, mbe_vd_mul(b.c, v));
        r3 = mbe_vd_add(r3, mbe_vd_mul(b.s, v));
        tone_oscillator_step(&a);
        tone_oscillator_step(&b);
    }
    gram[0][0] = mbe_vd_sum(g00);
    gram[1][0] = mbe_vd_sum(g10);
    gram[1][1] = mbe_vd_sum(g11);
    gram[2][0] = mbe_vd_sum(g20);
    gram[2][1] = mbe_vd_sum(g21);
    gram[2][2] = mbe_vd_sum(g22);
    gram[3][0] = mbe_vd_sum(g30);
    gram[3][1] = mbe_vd_sum(g31);
    gram[3][2] = mbe_vd_sum(g32);
    gram[3][3] = mbe_vd_sum(g33);
    rhs[0] = mbe_vd_sum(r0);
    rhs[1] = mbe_vd_sum(r1);
    rhs[2] = mbe_vd_sum(r2);
    rhs[3] = mbe_vd_sum(r3);
    double v[4][MBE_VD_WIDTH];
    mbe_vd_store(v[0], a.c);
    mbe_vd_store(v[1], a.s);
    mbe_vd_store(v[2], b.c);
    mbe_vd_store(v[3], b.s);
    for (int k = 0; i < n; ++i, ++k) {
        for (int r = 0; r < 4; ++r) {
            rhs[r] += v[r][k] * (double)x[i];
            for (int q = 0; q <= r; ++q) {
                gram[r][q] += v[r][k] * v[q][k];
            }
        }
    }
}

double
mbe_tone_fit(const float* samples, int n, const double hz[2], int count, double amplitude[2]) {
    if (!samples || !hz || !amplitude || n < 1 || n > 160 || count < 1 || count > 2) {
        return 0.0;
    }
    double gram[4][4] = {{0}}, rhs[4] = {0}, coef[4] = {0};
    if (count == 1) {
        tone_accumulate_one(samples, n, hz[0], gram, rhs);
    } else {
        tone_accumulate_two(samples, n, hz, gram, rhs);
    }
    int dim = 2 * count;
    for (int r = 0; r < dim; ++r) {
        for (int q = r + 1; q < dim; ++q) {
            gram[r][q] = gram[q][r];
        }
    }
    if (!tone_solve(gram, rhs, dim, coef)) {
        return 0.0;
    }
    double explained = 0.0;
    for (int r = 0; r < dim; ++r) {
        explained += rhs[r] * coef[r];
    }
    amplitude[0] = hypot(coef[0], coef[1]);
    if (count == 2) {
        amplitude[1] = hypot(coef[2], coef[3]);
    }
    return explained;
}
