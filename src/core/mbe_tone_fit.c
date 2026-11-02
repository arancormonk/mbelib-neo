// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Single-precision, bounded SIMD tone projection kernels.
 *
 * mbe_tone_fit_avx2.c builds this file a second time for AVX2/FMA (mbe_cpu.h). */
#include <math.h>
#include "mbe_cpu.h"
#include "mbe_simd.h"
#include "mbe_tone_fit.h"

/* Keep the scalar build available without changing the shared SIMD layer. */
#ifndef MBE_VF_WIDTH
typedef float mbe_vf;
#define MBE_VF_WIDTH 1

static inline mbe_vf
mbe_vf_set(float x) {
    return x;
}

static inline mbe_vf
mbe_vf_add(mbe_vf a, mbe_vf b) {
    return a + b;
}

static inline mbe_vf
mbe_vf_sub(mbe_vf a, mbe_vf b) {
    return a - b;
}

static inline mbe_vf
mbe_vf_mul(mbe_vf a, mbe_vf b) {
    return a * b;
}

static inline mbe_vf
mbe_vf_load(const float* x) {
    return *x;
}

static inline void
mbe_vf_store(float* x, mbe_vf a) {
    *x = a;
}
#endif

static float
tone_sum(mbe_vf v) {
    float lanes[MBE_VF_WIDTH];
    mbe_vf_store(lanes, v);
    float sum = 0.0f;
    for (int i = 0; i < MBE_VF_WIDTH; ++i) {
        sum += lanes[i];
    }
    return sum;
}

/* Solve the symmetric positive definite system a x = b of order n <= 4 by
 * Cholesky factorization; returns 0 when a is not positive definite. */
static int
tone_solve(float a[4][4], const float b[4], int n, float x[4]) {
    float l[4][4] = {{0}};
    float y[4] = {0};
    if (n < 1 || n > 4) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            float sum = a[i][j];
            for (int k = 0; k < j; k++) {
                sum -= l[i][k] * l[j][k];
            }
            if (i == j) {
                if (sum <= 0.0f) {
                    return 0;
                }
                l[i][i] = sqrtf(sum);
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        float sum = b[i];
        for (int k = 0; k < i; k++) {
            sum -= l[i][k] * y[k];
        }
        y[i] = sum / l[i][i];
    }
    for (int i = n - 1; i >= 0; i--) {
        float sum = y[i];
        for (int k = i + 1; k < n; k++) {
            sum -= l[k][i] * x[k];
        }
        x[i] = sum / l[i][i];
    }
    return 1;
}

/* Each lane follows one residue class of sample indices. Seed consecutive
 * samples with the scalar recurrence, then rotate by WIDTH samples at a time.
 * All oscillator and projection arithmetic stays in single precision. */
struct tone_oscillator {
    mbe_vf c, s, cw, sw;
};

static struct tone_oscillator
tone_oscillator_init(float hz) {
    float w = 2.0f * (float)M_PI * hz / MBE_TONE_FS;
    float cw = cosf(w), sw = sinf(w);
    float c[MBE_VF_WIDTH], s[MBE_VF_WIDTH];
    float cc = 1.0f, ss = 0.0f;
    for (int i = 0; i < MBE_VF_WIDTH; ++i) {
        c[i] = cc;
        s[i] = ss;
        float next = cc * cw - ss * sw;
        ss = ss * cw + cc * sw;
        cc = next;
    }
    struct tone_oscillator o = {mbe_vf_load(c), mbe_vf_load(s), mbe_vf_set(cc), mbe_vf_set(ss)};
    return o;
}

static void
tone_oscillator_step(struct tone_oscillator* o) {
    mbe_vf next = mbe_vf_sub(mbe_vf_mul(o->c, o->cw), mbe_vf_mul(o->s, o->sw));
    o->s = mbe_vf_add(mbe_vf_mul(o->s, o->cw), mbe_vf_mul(o->c, o->sw));
    o->c = next;
}

static void
tone_accumulate_one(const float* x, int n, float hz, float gram[4][4], float rhs[4]) {
    struct tone_oscillator o = tone_oscillator_init(hz);
    mbe_vf cc = mbe_vf_set(0.0f), cs = cc, ss = cc, cx = cc, sx = cc;
    int i = 0;
    for (; i + MBE_VF_WIDTH <= n; i += MBE_VF_WIDTH) {
        mbe_vf v = mbe_vf_load(x + i);
        cc = mbe_vf_add(cc, mbe_vf_mul(o.c, o.c));
        cs = mbe_vf_add(cs, mbe_vf_mul(o.c, o.s));
        ss = mbe_vf_add(ss, mbe_vf_mul(o.s, o.s));
        cx = mbe_vf_add(cx, mbe_vf_mul(o.c, v));
        sx = mbe_vf_add(sx, mbe_vf_mul(o.s, v));
        tone_oscillator_step(&o);
    }
    gram[0][0] = tone_sum(cc);
    gram[1][0] = tone_sum(cs);
    gram[1][1] = tone_sum(ss);
    rhs[0] = tone_sum(cx);
    rhs[1] = tone_sum(sx);
    float c[MBE_VF_WIDTH], s[MBE_VF_WIDTH];
    mbe_vf_store(c, o.c);
    mbe_vf_store(s, o.s);
    for (int k = 0; i < n; ++i, ++k) {
        gram[0][0] += c[k] * c[k];
        gram[1][0] += c[k] * s[k];
        gram[1][1] += s[k] * s[k];
        rhs[0] += c[k] * x[i];
        rhs[1] += s[k] * x[i];
    }
}

static void
tone_accumulate_two(const float* x, int n, const float hz[2], float gram[4][4], float rhs[4]) {
    struct tone_oscillator a = tone_oscillator_init(hz[0]), b = tone_oscillator_init(hz[1]);
    mbe_vf g00 = mbe_vf_set(0.0f), g10 = g00, g11 = g00, g20 = g00, g21 = g00, g22 = g00;
    mbe_vf g30 = g00, g31 = g00, g32 = g00, g33 = g00, r0 = g00, r1 = g00, r2 = g00, r3 = g00;
    int i = 0;
    for (; i + MBE_VF_WIDTH <= n; i += MBE_VF_WIDTH) {
        mbe_vf v = mbe_vf_load(x + i);
        g00 = mbe_vf_add(g00, mbe_vf_mul(a.c, a.c));
        g10 = mbe_vf_add(g10, mbe_vf_mul(a.s, a.c));
        g11 = mbe_vf_add(g11, mbe_vf_mul(a.s, a.s));
        g20 = mbe_vf_add(g20, mbe_vf_mul(b.c, a.c));
        g21 = mbe_vf_add(g21, mbe_vf_mul(b.c, a.s));
        g22 = mbe_vf_add(g22, mbe_vf_mul(b.c, b.c));
        g30 = mbe_vf_add(g30, mbe_vf_mul(b.s, a.c));
        g31 = mbe_vf_add(g31, mbe_vf_mul(b.s, a.s));
        g32 = mbe_vf_add(g32, mbe_vf_mul(b.s, b.c));
        g33 = mbe_vf_add(g33, mbe_vf_mul(b.s, b.s));
        r0 = mbe_vf_add(r0, mbe_vf_mul(a.c, v));
        r1 = mbe_vf_add(r1, mbe_vf_mul(a.s, v));
        r2 = mbe_vf_add(r2, mbe_vf_mul(b.c, v));
        r3 = mbe_vf_add(r3, mbe_vf_mul(b.s, v));
        tone_oscillator_step(&a);
        tone_oscillator_step(&b);
    }
    gram[0][0] = tone_sum(g00);
    gram[1][0] = tone_sum(g10);
    gram[1][1] = tone_sum(g11);
    gram[2][0] = tone_sum(g20);
    gram[2][1] = tone_sum(g21);
    gram[2][2] = tone_sum(g22);
    gram[3][0] = tone_sum(g30);
    gram[3][1] = tone_sum(g31);
    gram[3][2] = tone_sum(g32);
    gram[3][3] = tone_sum(g33);
    rhs[0] = tone_sum(r0);
    rhs[1] = tone_sum(r1);
    rhs[2] = tone_sum(r2);
    rhs[3] = tone_sum(r3);
    float v[4][MBE_VF_WIDTH];
    mbe_vf_store(v[0], a.c);
    mbe_vf_store(v[1], a.s);
    mbe_vf_store(v[2], b.c);
    mbe_vf_store(v[3], b.s);
    for (int k = 0; i < n; ++i, ++k) {
        for (int r = 0; r < 4; ++r) {
            rhs[r] += v[r][k] * x[i];
            for (int q = 0; q <= r; ++q) {
                gram[r][q] += v[r][k] * v[q][k];
            }
        }
    }
}

float
MBE_KFN(mbe_tone_fit)(const float* samples, int n, const float hz[2], int count, float amplitude[2]) {
    if (!samples || !hz || !amplitude || n < 1 || n > 160 || count < 1 || count > 2) {
        return 0.0f;
    }
    float gram[4][4] = {{0}}, rhs[4] = {0}, coef[4] = {0};
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
        return 0.0f;
    }
    float explained = 0.0f;
    for (int r = 0; r < dim; ++r) {
        explained += rhs[r] * coef[r];
    }
    amplitude[0] = hypotf(coef[0], coef[1]);
    if (count == 2) {
        amplitude[1] = hypotf(coef[2], coef[3]);
    }
    return explained;
}

#if defined(MBE_AVX2_DISPATCH)
float
mbe_tone_fit(const float* samples, int n, const float hz[2], int count, float amplitude[2]) {
    if (mbe_cpu_use_avx2()) {
        return mbe_tone_fit_avx2(samples, n, hz, count, amplitude);
    }
    return mbe_tone_fit_sse2(samples, n, hz, count, amplitude);
}
#endif
