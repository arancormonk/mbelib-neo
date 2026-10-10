// SPDX-License-Identifier: GPL-2.0-or-later
#include "mbe_analysis_kernels.h"
#include "mbe_compiler.h"
#if defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_SSE2)
#include <emmintrin.h>
typedef __m128 band_vector;
#define BAND_VECTOR 1
#define band_load   _mm_loadu_ps
#define band_set    _mm_set1_ps
#define band_add    _mm_add_ps
#define band_sub    _mm_sub_ps
#define band_mul    _mm_mul_ps
#define band_store  _mm_storeu_ps
#elif defined(MBELIB_ENABLE_SIMD) && defined(MBE_SIMD_TARGET_NEON)
#include <arm_neon.h>
typedef float32x4_t band_vector;
#define BAND_VECTOR 1
#define band_load   vld1q_f32
#define band_set    vdupq_n_f32
#define band_add    vaddq_f32
#define band_sub    vsubq_f32
#define band_mul    vmulq_f32
#define band_store  vst1q_f32
#endif

#if defined(BAND_VECTOR)
static float
band_sum(band_vector v) {
    float lane[4];
    band_store(lane, v);
    return ((lane[0] + lane[1]) + lane[2]) + lane[3];
}
#endif

void
mbe_analysis_band_sums(const float* re, const float* im, const float* w, int n, float sums[4]) {
    float cr = 0.0f, ci = 0.0f, q = 0.0f, e = 0.0f;
    int i = 0;
#if defined(BAND_VECTOR)
    band_vector vr = band_set(0.0f), vi = vr, vq = vr, ve = vr;
    for (; i + 4 <= n; i += 4) {
        band_vector r = band_load(re + i), s = band_load(im + i), t = band_load(w + i);
        vr = band_add(vr, band_mul(r, t));
        vi = band_add(vi, band_mul(s, t));
        vq = band_add(vq, band_mul(t, t));
        ve = band_add(ve, band_add(band_mul(r, r), band_mul(s, s)));
    }
    cr = band_sum(vr);
    ci = band_sum(vi);
    q = band_sum(vq);
    e = band_sum(ve);
#endif
    for (; i < n; ++i) {
        cr += re[i] * w[i];
        ci += im[i] * w[i];
        q += w[i] * w[i];
        e += re[i] * re[i] + im[i] * im[i];
    }
    sums[0] = cr;
    sums[1] = ci;
    sums[2] = q;
    sums[3] = e;
}

float
mbe_analysis_residual(const float* re, const float* im, const float* w, int n, float ar, float ai) {
    float sum = 0.0f;
    int i = 0;
#if defined(BAND_VECTOR)
    band_vector acc = band_set(0.0f), vr = band_set(ar), vi = band_set(ai);
    for (; i + 4 <= n; i += 4) {
        band_vector t = band_load(w + i);
        band_vector dr = band_sub(band_load(re + i), band_mul(vr, t));
        band_vector di = band_sub(band_load(im + i), band_mul(vi, t));
        acc = band_add(acc, band_add(band_mul(dr, dr), band_mul(di, di)));
    }
    sum = band_sum(acc);
#endif
    for (; i < n; ++i) {
        float dr = re[i] - ar * w[i];
        float di = im[i] - ai * w[i];
        sum += dr * dr + di * di;
    }
    return sum;
}
