// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Voiced oscillators: four independent sample streams per SIMD vector. */
#include <math.h>
#include "mbe_math.h"
#include "mbe_simd.h" // IWYU pragma: keep (compile-target SIMD selection)
#include "mbe_voiced.h"

/* Keep the direct phase expression for scalar targets and caller-supplied
 * phases beyond the normal decoder range, where float phase rounding dominates. */
static void
voiced_interpolated_scalar(float out[160], float phase, float prev_amp, float cur_amp, float frequency,
                           float pitch_delta, int harmonic) {
    for (int n = 0; n < 160; ++n) {
        float theta = phase + frequency * (float)n + pitch_delta * (float)(harmonic * n * n) / 320.0f;
        float amp = prev_amp + ((float)n / 160.0f) * (cur_amp - prev_amp);
        out[n] += 2.0f * amp * cosf(theta);
    }
}

#if defined(MBE_VF_WIDTH)
struct voiced_oscillator {
    mbe_vf c, s, cd, sd;
};

static void
voiced_step(struct voiced_oscillator* o) {
    mbe_vf next = mbe_vf_sub(mbe_vf_mul(o->c, o->cd), mbe_vf_mul(o->s, o->sd));
    o->s = mbe_vf_add(mbe_vf_mul(o->s, o->cd), mbe_vf_mul(o->c, o->sd));
    o->c = next;
}

static struct voiced_oscillator
voiced_init(float phase, float omega) {
    float c[4], s[4], cc, ss, cd, sd;
    mbe_sincosf(phase, &ss, &cc);
    mbe_sincosf(omega, &sd, &cd);
    for (int k = 0; k < 4; ++k) {
        c[k] = cc;
        s[k] = ss;
        float next = cc * cd - ss * sd;
        ss = ss * cd + cc * sd;
        cc = next;
    }
    /* The four-sample rotation is evaluated once, not accumulated in time. */
    mbe_sincosf(4.0f * omega, &sd, &cd);
    struct voiced_oscillator o = {mbe_vf_load(c), mbe_vf_load(s), mbe_vf_set(cd), mbe_vf_set(sd)};
    return o;
}

static void
voiced_window_one(float out[160], const float* window, const struct mbe_voiced_component* component) {
    struct voiced_oscillator o = voiced_init(component->phase, component->omega);
    mbe_vf gain = mbe_vf_set(component->gain);
    for (int n = 0; n < 160; n += 4) {
        mbe_vf v = mbe_vf_mul(mbe_vf_mul(o.c, mbe_vf_load(window + n)), gain);
        mbe_vf_store(out + n, mbe_vf_add(mbe_vf_load(out + n), v));
        voiced_step(&o);
    }
}

void
mbe_voiced_windowed(float out[160], const float window[320], const struct mbe_voiced_component* prev,
                    const struct mbe_voiced_component* cur) {
    if (!prev->active) {
        if (cur->active) {
            voiced_window_one(out, window, cur);
        }
        return;
    }
    if (!cur->active) {
        voiced_window_one(out, window + 160, prev);
        return;
    }
    struct voiced_oscillator p = voiced_init(prev->phase, prev->omega), c = voiced_init(cur->phase, cur->omega);
    mbe_vf pg = mbe_vf_set(prev->gain), cg = mbe_vf_set(cur->gain);
    for (int n = 0; n < 160; n += 4) {
        mbe_vf pv = mbe_vf_mul(mbe_vf_mul(p.c, mbe_vf_load(window + n + 160)), pg);
        mbe_vf cv = mbe_vf_mul(mbe_vf_mul(c.c, mbe_vf_load(window + n)), cg);
        mbe_vf sum = mbe_vf_add(mbe_vf_load(out + n), pv);
        mbe_vf_store(out + n, mbe_vf_add(sum, cv));
        voiced_step(&p);
        voiced_step(&c);
    }
}

/* theta(n) = phase + frequency*n + b*n*n. Each lane advances by four
 * samples, while its rotation advances by 32*b. Use double for the seeds
 * and rotations to avoid cancellation in slowly changing pitch. */
static struct voiced_oscillator
voiced_chirp_init(float phase, float frequency, double b) {
    float c[4], s[4], cd[4], sd[4];
    for (int k = 0; k < 4; ++k) {
        double theta = (double)phase + (double)frequency * k + b * k * k;
        double delta = 4.0 * frequency + b * (8 * k + 16);
        c[k] = (float)cos(theta);
        s[k] = (float)sin(theta);
        cd[k] = (float)cos(delta);
        sd[k] = (float)sin(delta);
    }
    struct voiced_oscillator o = {mbe_vf_load(c), mbe_vf_load(s), mbe_vf_load(cd), mbe_vf_load(sd)};
    return o;
}

void
mbe_voiced_interpolated(float out[160], float phase, float prev_amp, float cur_amp, float frequency, float pitch_delta,
                        int harmonic) {
    if (fabsf(phase) > 512.0f) {
        voiced_interpolated_scalar(out, phase, prev_amp, cur_amp, frequency, pitch_delta, harmonic);
        return;
    }
    double b = (double)pitch_delta * harmonic / 320.0;
    struct voiced_oscillator o = voiced_chirp_init(phase, frequency, b);
    mbe_vf delta_c = mbe_vf_set((float)cos(32.0 * b)), delta_s = mbe_vf_set((float)sin(32.0 * b));
    const float initial[4] = {0.0f, 1.0f, 2.0f, 3.0f};
    mbe_vf n = mbe_vf_load(initial), inv_n = mbe_vf_set(1.0f / 160.0f);
    mbe_vf start = mbe_vf_set(prev_amp), difference = mbe_vf_set(cur_amp - prev_amp);
    for (int i = 0; i < 160; i += 4) {
        mbe_vf amp = mbe_vf_add(start, mbe_vf_mul(mbe_vf_mul(n, inv_n), difference));
        mbe_vf v = mbe_vf_mul(mbe_vf_mul(mbe_vf_set(2.0f), amp), o.c);
        mbe_vf_store(out + i, mbe_vf_add(mbe_vf_load(out + i), v));
        voiced_step(&o);
        mbe_vf next = mbe_vf_sub(mbe_vf_mul(o.cd, delta_c), mbe_vf_mul(o.sd, delta_s));
        o.sd = mbe_vf_add(mbe_vf_mul(o.sd, delta_c), mbe_vf_mul(o.cd, delta_s));
        o.cd = next;
        n = mbe_vf_add(n, mbe_vf_set(4.0f));
    }
}
#else
/* Scalar reference retains the historical recurrence, arithmetic and order. */
static void
voiced_window_scalar(float out[160], const float* window, const struct mbe_voiced_component* component) {
    if (!component->active) {
        return;
    }
    float cc, ss, cd, sd;
    mbe_sincosf(component->phase, &ss, &cc);
    mbe_sincosf(component->omega, &sd, &cd);
    for (int n = 0; n < 160; ++n) {
        out[n] += component->gain * window[n] * cc;
        float next = cc * cd - ss * sd;
        ss = ss * cd + cc * sd;
        cc = next;
    }
}

void
mbe_voiced_windowed(float out[160], const float window[320], const struct mbe_voiced_component* prev,
                    const struct mbe_voiced_component* cur) {
    voiced_window_scalar(out, window + 160, prev);
    voiced_window_scalar(out, window, cur);
}

void
mbe_voiced_interpolated(float out[160], float phase, float prev_amp, float cur_amp, float frequency, float pitch_delta,
                        int harmonic) {
    voiced_interpolated_scalar(out, phase, prev_amp, cur_amp, frequency, pitch_delta, harmonic);
}
#endif
