// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Voiced oscillators: MBE_VF_WIDTH independent sample streams per SIMD vector. */
#include "mbe_voiced.h"

#if defined(MBE_VOICED_KERNELS)
#include <math.h>
#include "mbe_math.h"
#include "mbe_simd.h"

#define FRAME MBE_VOICED_FRAME
#define WIDTH MBE_VF_WIDTH
/* The kernels process whole vectors of one frame. The compiler checks this,
 * not the preprocessor, so analyzers that try configurations without a
 * vector width can still parse the file. */
typedef char mbe_voiced_whole_vectors[(FRAME % WIDTH == 0) ? 1 : -1];

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
    float c[WIDTH], s[WIDTH], cc, ss, cd, sd;
    mbe_sincosf(phase, &ss, &cc);
    mbe_sincosf(omega, &sd, &cd);
    for (int k = 0; k < WIDTH; ++k) {
        c[k] = cc;
        s[k] = ss;
        float next = cc * cd - ss * sd;
        ss = ss * cd + cc * sd;
        cc = next;
    }
    /* The WIDTH-sample rotation is evaluated once, not accumulated in time. */
    mbe_sincosf((float)WIDTH * omega, &sd, &cd);
    struct voiced_oscillator o = {mbe_vf_load(c), mbe_vf_load(s), mbe_vf_set(cd), mbe_vf_set(sd)};
    return o;
}

static void
voiced_window_one(float out[FRAME], const float* window, const struct mbe_voiced_component* component) {
    struct voiced_oscillator o = voiced_init(component->phase, component->omega);
    mbe_vf gain = mbe_vf_set(component->gain);
    for (int n = 0; n < FRAME; n += WIDTH) {
        mbe_vf v = mbe_vf_mul(mbe_vf_mul(o.c, mbe_vf_load(window + n)), gain);
        mbe_vf_store(out + n, mbe_vf_add(mbe_vf_load(out + n), v));
        voiced_step(&o);
    }
}

void
mbe_voiced_windowed(float out[FRAME], const float window[2 * FRAME], const struct mbe_voiced_component* prev,
                    const struct mbe_voiced_component* cur) {
    if (!prev->active) {
        if (cur->active) {
            voiced_window_one(out, window, cur);
        }
        return;
    }
    if (!cur->active) {
        voiced_window_one(out, window + FRAME, prev);
        return;
    }
    struct voiced_oscillator p = voiced_init(prev->phase, prev->omega), c = voiced_init(cur->phase, cur->omega);
    mbe_vf pg = mbe_vf_set(prev->gain), cg = mbe_vf_set(cur->gain);
    for (int n = 0; n < FRAME; n += WIDTH) {
        mbe_vf pv = mbe_vf_mul(mbe_vf_mul(p.c, mbe_vf_load(window + n + FRAME)), pg);
        mbe_vf cv = mbe_vf_mul(mbe_vf_mul(c.c, mbe_vf_load(window + n)), cg);
        mbe_vf sum = mbe_vf_add(mbe_vf_load(out + n), pv);
        mbe_vf_store(out + n, mbe_vf_add(sum, cv));
        voiced_step(&p);
        voiced_step(&c);
    }
}

/* theta(n) = phase + frequency*n + b*n*n. Each lane advances by WIDTH
 * samples, while its rotation advances by 2*WIDTH*WIDTH*b. Use double for the
 * seeds and rotations to avoid cancellation in slowly changing pitch. Each
 * lane's seed rotates the phase's own cosine and sine by its offset, so no
 * float phase, however large, absorbs the offsets; every phase takes this one
 * path. */
static struct voiced_oscillator
voiced_chirp_init(float phase, float frequency, double b) {
    float c[WIDTH], s[WIDTH], cd[WIDTH], sd[WIDTH];
    const double phase_c = cos((double)phase), phase_s = sin((double)phase);
    for (int k = 0; k < WIDTH; ++k) {
        double offset = (double)frequency * k + b * k * k;
        double delta = (double)WIDTH * frequency + b * (2 * WIDTH * k + WIDTH * WIDTH);
        c[k] = (float)(phase_c * cos(offset) - phase_s * sin(offset));
        s[k] = (float)(phase_s * cos(offset) + phase_c * sin(offset));
        cd[k] = (float)cos(delta);
        sd[k] = (float)sin(delta);
    }
    struct voiced_oscillator o = {mbe_vf_load(c), mbe_vf_load(s), mbe_vf_load(cd), mbe_vf_load(sd)};
    return o;
}

void
mbe_voiced_interpolated(float out[FRAME], float phase, float prev_amp, float cur_amp, float frequency,
                        float pitch_delta, int harmonic) {
    double b = (double)pitch_delta * harmonic / (2.0 * FRAME);
    struct voiced_oscillator o = voiced_chirp_init(phase, frequency, b);
    double rotation = (double)(2 * WIDTH * WIDTH) * b;
    mbe_vf delta_c = mbe_vf_set((float)cos(rotation)), delta_s = mbe_vf_set((float)sin(rotation));
    float initial[WIDTH];
    for (int k = 0; k < WIDTH; ++k) {
        initial[k] = (float)k;
    }
    mbe_vf n = mbe_vf_load(initial), inv_n = mbe_vf_set(1.0f / (float)FRAME);
    mbe_vf start = mbe_vf_set(prev_amp), difference = mbe_vf_set(cur_amp - prev_amp);
    for (int i = 0; i < FRAME; i += WIDTH) {
        mbe_vf amp = mbe_vf_add(start, mbe_vf_mul(mbe_vf_mul(n, inv_n), difference));
        mbe_vf v = mbe_vf_mul(mbe_vf_mul(mbe_vf_set(2.0f), amp), o.c);
        mbe_vf_store(out + i, mbe_vf_add(mbe_vf_load(out + i), v));
        voiced_step(&o);
        mbe_vf next = mbe_vf_sub(mbe_vf_mul(o.cd, delta_c), mbe_vf_mul(o.sd, delta_s));
        o.sd = mbe_vf_add(mbe_vf_mul(o.sd, delta_c), mbe_vf_mul(o.cd, delta_s));
        o.cd = next;
        n = mbe_vf_add(n, mbe_vf_set((float)WIDTH));
    }
}
#else
/* Scalar targets keep the historical loops in mbelib.c. ISO C forbids an
 * empty translation unit. */
typedef int mbe_voiced_scalar_target;
#endif
