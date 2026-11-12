// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Speech synthesis in 10 ms subframes: the mid-frame model, the pulse
 * train and its state, and the frame marks the decoders set.
 */
#include <assert.h>
#include <fenv.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "mbe_adaptive.h"
#include "mbe_subframe.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Where noiseOverlap holds the pulse phase (mbe_subframe.h). */
#define PHASE_SLOT (1 + MBE_PULSE_HISTORY)

/* Whether the floating-point overflow flag tells what the library computed:
 * in IEEE arithmetic the compiler has not rearranged (not fast-math, not
 * optimized), on the compilers whose unoptimized code keeps the flags. */
#if (defined(__GNUC__) || defined(__clang__)) && defined(FE_OVERFLOW) && !defined(__FAST_MATH__)                       \
    && !defined(__OPTIMIZE__)
#define MID_MODEL_OVERFLOW_IS_VISIBLE 1
#else
#define MID_MODEL_OVERFLOW_IS_VISIBLE 0
#endif

/* Flush-to-zero and denormals-are-zero, the floating-point environment audio
 * hosts commonly run in, where the target has the switches (SSE). Returns the
 * state to put back. */
#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && (_M_IX86_FP >= 2))
#include <xmmintrin.h>

static unsigned
flush_denormals(void) {
    const unsigned saved = _mm_getcsr();
    _mm_setcsr(saved | 0x8040u);
    return saved;
}

static void
restore_denormals(unsigned saved) {
    _mm_setcsr(saved);
}
#else
static unsigned
flush_denormals(void) {
    return 0u;
}

static void
restore_denormals(unsigned saved) {
    (void)saved;
}
#endif

static float
w0_of(double hz) {
    return (float)(2.0 * M_PI * hz / 8000.0);
}

/* A model of L harmonics at `hz`, the first `voiced` of them voiced. */
static void
set_model(mbe_parms* mp, double hz, int L, int voiced, float amplitude) {
    mp->w0 = w0_of(hz);
    mp->L = L;
    for (int l = 1; l <= 56; ++l) {
        mp->Vl[l] = (l <= voiced) ? 1 : 0;
        mp->Ml[l] = (l <= L) ? amplitude : 0.0f;
    }
    mbe_set_frame_mark(mp, MBE_FRAME_MARK_NONE);
}

/* The next frame of a stream: `last` is the frame before it as synthesized,
 * `model` the new frame's model and mark. Renders the voiced part alone (no
 * transform plan), at unit gain, and leaves the new state in `last`. */
static void
next_frame_from(float out[160], mbe_parms* last, const mbe_parms* model, int from_silence) {
    mbe_parms cur = *last;
    cur.w0 = model->w0;
    cur.L = model->L;
    memcpy(cur.Vl, model->Vl, sizeof(cur.Vl));
    memcpy(cur.Ml, model->Ml, sizeof(cur.Ml));
    mbe_set_frame_mark(&cur, mbe_frame_mark(model));
    memset(out, 0, 160 * sizeof(float));
    mbe_synthesize_subframes(out, &cur, last, NULL, 1.0f, from_silence);
    *last = cur;
}

static void
next_frame(float out[160], mbe_parms* last, const mbe_parms* model) {
    next_frame_from(out, last, model, 0);
}

static void
fresh_stream(mbe_parms* last) {
    mbe_parms prev, enhanced;
    mbe_initMbeParms(last, &prev, &enhanced);
    set_model(last, 100.0, 30, 0, 0.0f);
}

/* Cycles since the last pulse at the end of the frame `mp` holds. */
static double
pulse_phase(const mbe_parms* mp) {
    return (double)mp->noiseOverlap[PHASE_SLOT];
}

/* a - b in cycles, wrapped to [-0.5, 0.5). */
static double
cycles_apart(double a, double b) {
    const double d = a - b;
    return d - floor(d + 0.5);
}

static double
cycles_of(double hz, double samples) {
    return hz * samples / 8000.0;
}

static double
component_amplitude(const float* x, int n, double hz) {
    double re = 0.0, im = 0.0;
    for (int i = 0; i < n; ++i) {
        const double t = 2.0 * M_PI * hz * i / 8000.0;
        re += x[i] * cos(t);
        im += x[i] * sin(t);
    }
    return 2.0 * sqrt(re * re + im * im) / n;
}

/* Phase (cosine convention) of the component at `hz`, with the time origin at
 * sample `origin`. */
static double
component_phase(const float* x, int n, double hz, int origin) {
    double re = 0.0, im = 0.0;
    for (int i = 0; i < n; ++i) {
        const double t = 2.0 * M_PI * hz * (i - origin) / 8000.0;
        re += x[i] * cos(t);
        im -= x[i] * sin(t);
    }
    return atan2(im, re);
}

static int
all_finite(const float* x, int n) {
    for (int i = 0; i < n; ++i) {
        uint32_t bits;
        memcpy(&bits, &x[i], sizeof(bits));
        if ((bits & 0x7FFFFFFFu) >= 0x7F800000u) {
            return 0;
        }
    }
    return 1;
}

static int
same_bits(const float* a, const float* b, int n) {
    for (int i = 0; i < n; ++i) {
        uint32_t x, y;
        memcpy(&x, &a[i], sizeof(x));
        memcpy(&y, &b[i], sizeof(y));
        if (x != y) {
            return 0;
        }
    }
    return 1;
}

/* The frame marks: only an exact stored value counts, and one replaces another. */
static void
test_frame_marks(void) {
    static const enum mbe_frame_mark marks[] = {MBE_FRAME_MARK_NONE, MBE_FRAME_MARK_PITCH_CARRIES,
                                                MBE_FRAME_MARK_PITCH_ABSENT, MBE_FRAME_MARK_PREVIOUS_HELD,
                                                MBE_FRAME_MARK_LATE_ONSET};

    static const uint32_t others[] = {0x7FC00000u, 0x7F800000u, 0xFF800000u, 0x80000000u, 0x3F800001u, 0x40400001u,
                                      0x40800000u, 0x429E0000u, 0xBF800000u, 0xC0000000u, 0x00000001u};
    mbe_parms mp;
    memset(&mp, 0, sizeof(mp));
    assert(mbe_frame_mark(&mp) == MBE_FRAME_MARK_NONE);
    for (size_t i = 0; i < sizeof(marks) / sizeof(marks[0]); ++i) {
        mbe_set_frame_mark(&mp, marks[i]);
        assert(mbe_frame_mark(&mp) == marks[i]);
        assert(mbe_voicing_onset_delay(&mp) == ((marks[i] == MBE_FRAME_MARK_LATE_ONSET) ? MBE_VOICING_ONSET_LATE : 0));
    }
    for (size_t i = 0; i < sizeof(others) / sizeof(others[0]); ++i) {
        memcpy(&mp.PHIl[0], &others[i], sizeof(others[i]));
        assert(mbe_frame_mark(&mp) == MBE_FRAME_MARK_NONE);
        assert(mbe_voicing_onset_delay(&mp) == 0);
    }
    mbe_set_voicing_onset_delay(&mp, 1);
    assert(mbe_frame_mark(&mp) == MBE_FRAME_MARK_LATE_ONSET);
    mbe_set_voicing_onset_delay(&mp, 0);
    assert(mbe_frame_mark(&mp) == MBE_FRAME_MARK_NONE);
}

/* A steady voiced model renders as its harmonics: each at 2 Ml times the
 * roll-off's gain and the voiced gain passed in, with the phase the frame
 * stores (PHIl) plus the roll-off's at the frame's end. 187.5 Hz has a period
 * of 128/3 samples, so the pulses fall between samples. */
static void
test_steady_model_is_its_harmonics(void) {
    static float out[10 * 160];
    mbe_parms last, model;
    fresh_stream(&last);
    set_model(&model, 187.5, 18, 18, 0.0f);
    for (int l = 1; l <= model.L; ++l) {
        model.Ml[l] = 8.0f + (float)((l * 7) % 5);
    }
    for (int f = 0; f < 10; ++f) {
        next_frame(out + ((size_t)f * 160u), &last, &model);
    }
    /* Thirty periods from the third frame on. */
    for (int l = 1; l <= model.L; ++l) {
        const double hz = 187.5 * l;
        const double want = 2.0 * (double)model.Ml[l] * (double)mbe_low_band_gain((float)hz);
        const double got = component_amplitude(out + 320, 1280, hz);
        assert(fabs(got / want - 1.0) < 2e-3);
        const double phase = component_phase(&out[1344], 256, hz, 256); /* the last 256 samples */
        const double stored = (double)last.PHIl[l] + (double)mbe_low_band_phase((float)hz);
        assert(fabs(atan2(sin(phase - stored), cos(phase - stored))) < 0.02);
    }
    /* The gain argument scales the voiced part. */
    float unit[160], scaled[160];
    mbe_parms a = last, b = last, a_prev = last, b_prev = last;
    memset(unit, 0, sizeof(unit));
    memset(scaled, 0, sizeof(scaled));
    mbe_synthesize_subframes(unit, &a, &a_prev, NULL, 1.0f, 0);
    mbe_synthesize_subframes(scaled, &b, &b_prev, NULL, 1.25f, 0);
    for (int n = 0; n < 160; ++n) {
        assert(fabs((double)scaled[n] - 1.25 * (double)unit[n]) <= 1e-4 * (1.0 + fabs((double)unit[n])));
    }
}

/* The first voiced frame after one with no voiced band starts the train at
 * the frame's start, at its own pitch, whatever phases the stream held; the
 * frame's last pulse, and so its phase, follow from that alone. */
static void
test_onset_starts_the_train(void) {
    float out[2][160];
    mbe_parms stream[2], voiced;
    set_model(&voiced, 230.0, 15, 15, 16.0f);
    for (int history = 0; history < 2; ++history) {
        mbe_parms* last = &stream[history];
        fresh_stream(last);
        set_model(last, 150.0, 20, 0, 4.0f);
        for (int l = 1; l <= 56; ++l) {
            last->PSIl[l] = history ? 1.7f * (float)l : 0.0f;
            last->PHIl[l] = history ? 3.1f * (float)l : 0.0f;
        }
        /* A phase and pulses a frame with no voiced band cannot have left. */
        last->noiseOverlap[0] = history ? 3.0f : 0.0f;
        last->noiseOverlap[1] = history ? -12.0f : 0.0f;
        last->noiseOverlap[PHASE_SLOT] = history ? 0.4f : 0.0f;
        next_frame(out[history], last, &voiced);
    }
    const mbe_parms* first = &stream[0];
    const float* a = out[0];
    assert(same_bits(out[0], out[1], 160));
    assert(fabs(cycles_apart(pulse_phase(first), cycles_of(230.0, 160.0))) < 1e-5);
    /* Pulses every 8000/230 samples from the frame's start: five of them. */
    assert(first->noiseOverlap[0] == 5.0f);
    assert(first->noiseOverlap[5] == -160.0f);
    assert(fabs((double)first->noiseOverlap[4] - (-160.0 + 8000.0 / 230.0)) < 1e-3);
    /* Harmonic l's phase at the frame's end, as the next frame reads it. */
    for (int l = 1; l <= 15; ++l) {
        const double expect =
            (double)l * (2.0 * M_PI * pulse_phase(first) - (double)MBE_PULSE_PEAK * (double)voiced.w0);
        assert(fabs((double)first->PSIl[l] - expect) < 1e-3 * (1.0 + fabs(expect)));
    }
    /* Nothing sounds before the first response rises. */
    assert(a[0] == 0.0f);
    double early = 0.0, late = 0.0;
    for (int n = 0; n < 160; ++n) {
        *((n < 80) ? &early : &late) += (double)a[n] * (double)a[n];
    }
    assert(early > 0.0 && late > early);
}

/* A frame marked as a late onset after a frame with no voiced band starts its
 * voiced bands in its second subframe: no pulse before sample 80. After a
 * voiced frame the mark changes nothing. */
static void
test_late_onset(void) {
    float out[160], plain[160];
    mbe_parms last, voiced;
    fresh_stream(&last);
    set_model(&last, 150.0, 20, 0, 4.0f);
    set_model(&voiced, 230.0, 15, 15, 16.0f);
    mbe_set_frame_mark(&voiced, MBE_FRAME_MARK_LATE_ONSET);
    next_frame(out, &last, &voiced);
    assert(fabs(cycles_apart(pulse_phase(&last), cycles_of(230.0, 80.0))) < 1e-5);
    assert(last.noiseOverlap[0] == 3.0f);
    assert(last.noiseOverlap[3] == -80.0f);
    for (int n = 0; n < 80; ++n) {
        assert(out[n] == 0.0f);
    }
    assert(out[100] != 0.0f);

    mbe_parms marked = last, unmarked = last;
    next_frame(out, &marked, &voiced);
    mbe_set_frame_mark(&voiced, MBE_FRAME_MARK_NONE);
    next_frame(plain, &unmarked, &voiced);
    assert(same_bits(out, plain, 160));
}

/* Through voiced speech the pitch moves linearly from one model's to the
 * next: over a frame the pulse phase advances by the mean pitch of each
 * subframe, the mid-frame pitch being the geometric mean of the two frames'. */
static void
test_pitch_step_keeps_the_train(void) {
    float out[160];
    mbe_parms last, model;
    fresh_stream(&last);
    set_model(&model, 200.0, 18, 18, 16.0f);
    for (int f = 0; f < 4; ++f) {
        next_frame(out, &last, &model);
    }
    static const double steps[] = {230.0, 180.0, 207.0};
    double hz = 200.0;
    for (size_t i = 0; i < sizeof(steps) / sizeof(steps[0]); ++i) {
        const double before = pulse_phase(&last);
        const double mid = sqrt(hz * steps[i]);
        set_model(&model, steps[i], 16, 16, 16.0f);
        next_frame(out, &last, &model);
        const double advance = cycles_of(0.5 * (hz + mid), 80.0) + cycles_of(0.5 * (mid + steps[i]), 80.0);
        assert(fabs(cycles_apart(pulse_phase(&last), before + advance)) < 1e-4);
        /* The pulse spacing stays between the two pitch periods. */
        const int count = (int)last.noiseOverlap[0];
        const double shorter = 8000.0 / fmax(hz, steps[i]), longer = 8000.0 / fmin(hz, steps[i]);
        assert(count >= 3);
        for (int k = 1; k < count; ++k) {
            const double spacing = (double)last.noiseOverlap[k] - (double)last.noiseOverlap[k + 1];
            assert(spacing > shorter - 0.05 && spacing < longer + 0.05);
        }
        assert(all_finite(out, 160));
        hz = steps[i];
    }
}

/*
 * An independent rendering of a frame's voiced part, in double, from the rules
 * docs/dvsi-parity.md states: where the pulses fall (the pulse phase passing a
 * whole number, the pitch moving linearly across each subframe; found here by
 * bisection), and what each plays (both models of the subframe under the
 * linear crossfade, each through a window that rises over 16 samples, stays
 * whole until one pitch period after the pulse and falls over the 16 after).
 * The models are the library's, which test_mid_model checks on their own.
 */
struct reference_train {
    double phase;    /* cycles since the last pulse */
    int count;       /* pulses kept */
    double time[96]; /* in samples from the stream's start */
};

static double
reference_window(double t, double period) {
    if (t < 0.0 || t >= period + 16.0) {
        return 0.0;
    }
    if (t < 16.0) {
        return 0.5 * (1.0 - cos(M_PI * t / 16.0));
    }
    return (t < period) ? 1.0 : 0.5 * (1.0 + cos(M_PI * (t - period) / 16.0));
}

/* What a model plays t samples after a pulse. */
static double
reference_response(const struct mbe_subframe_model* model, double t) {
    double sum = 0.0;
    if (!model->any_voiced) {
        return 0.0;
    }
    for (int l = 1; l <= model->L; ++l) {
        const double w = (double)model->w0 * l;
        const double hz = w * 8000.0 / (2.0 * M_PI);
        if (model->voiced[l] != 1 || w >= M_PI) {
            continue;
        }
        sum += 2.0 * (double)model->magnitude[l] * (double)mbe_low_band_gain((float)hz)
               * cos((w * (t - 16.0)) + (double)model->phase[l] + (double)mbe_low_band_phase((float)hz));
    }
    return reference_window(t, 2.0 * M_PI / (double)model->w0) * sum;
}

/* The pulses of one subframe from sample `start`, between models a and b. */
static void
reference_advance(struct reference_train* train, double start, const struct mbe_subframe_model* a,
                  const struct mbe_subframe_model* b) {
    double from = (double)a->w0 / (2.0 * M_PI), to = (double)b->w0 / (2.0 * M_PI);
    if (!a->any_voiced && !b->any_voiced) {
        /* A new train starts after this subframe: nothing of the old one sounds on. */
        train->phase = 0.0;
        train->count = 0;
        return;
    }
    if (!b->any_voiced) {
        to = from;
    } else if (!a->any_voiced) {
        from = to;
    }
    const double first = train->phase;
    const double end = first + (80.0 * from) + ((to - from) * 80.0 * 80.0 / 160.0);
    for (int k = (first > 0.0) ? 1 : 0; (double)k < end; ++k) {
        double lo = 0.0, hi = 80.0;
        for (int step = 0; step < 60; ++step) {
            const double t = 0.5 * (lo + hi);
            if (first + (from * t) + ((to - from) * t * t / 160.0) < (double)k) {
                lo = t;
            } else {
                hi = t;
            }
        }
        assert(train->count < 96);
        train->time[train->count++] = start + hi;
    }
    train->phase = end - floor(end);
}

/* One subframe from sample `start`: every pulse so far plays model a going out and model b coming in. */
static void
reference_render(double* out, double start, const struct reference_train* train, const struct mbe_subframe_model* a,
                 const struct mbe_subframe_model* b) {
    for (int n = 0; n < 80; ++n) {
        const double fade = (double)n / 80.0;
        double sum = 0.0;
        for (int i = 0; i < train->count; ++i) {
            const double t = start + (double)n - train->time[i];
            sum += ((1.0 - fade) * reference_response(a, t)) + (fade * reference_response(b, t));
        }
        out[n] = sum;
    }
}

/* The frames a stream plays through pitch steps, changes of harmonic count,
 * level and voicing, onsets and the frame marks, against that rendering:
 * every sample, so that pulses, models, windows and crossfades have to agree
 * across the subframe and frame boundaries, not the stored state alone. */
static void
test_waveform_through_changes(void) {
    static const struct {
        double hz;
        int L, voiced;
        float amplitude;
        enum mbe_frame_mark mark;
        int from_silence;
    } frames[] = {
        {200.0, 18, 14, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* an onset */
        {200.0, 18, 14, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* steady */
        {230.0, 16, 16, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* pitch up, fewer harmonics */
        {180.0, 20, 12, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* pitch down, more of them, fewer voiced */
        {207.0, 17, 17, 40.0f, MBE_FRAME_MARK_NONE, 0},        /* a level step */
        {100.0, 30, 0, 8.0f, MBE_FRAME_MARK_PITCH_CARRIES, 0}, /* no voiced band, its pitch carried */
        {180.0, 20, 20, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* onset at the mean of the two pitches */
        {250.0, 14, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, 0},  /* the voicing before lasts to mid-frame */
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* and the train goes on */
        {150.0, 24, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0},   /* a mark that needs an unvoiced frame before it */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0},           /* voicing ends with the frame's start */
        {120.0, 28, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0},   /* voiced bands from mid-frame */
        {60.0, 56, 40, 6.0f, MBE_FRAME_MARK_NONE, 0},          /* the longest period: responses span frames */
        {60.0, 56, 40, 6.0f, MBE_FRAME_MARK_NONE, 0},
        {110.0, 30, 30, 16.0f, MBE_FRAME_MARK_NONE, 0},
        {250.0, 14, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, 0}, /* the train runs on through this frame, */
        {65.0, 52, 6, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0},   /* then starts anew at mid-frame: its pulses, */
        {65.0, 52, 30, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* under a period longer than a subframe, stay out */
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 0},
        {250.0, 15, 0, 1.6f, MBE_FRAME_MARK_PREVIOUS_HELD, 0}, /* the model before stands at mid-frame, */
        {170.0, 22, 22, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* and the train runs on into the frame after */
        {140.0, 26, 4, 16.0f, MBE_FRAME_MARK_NONE, 0},         /* voiced low only, held the same way */
        {250.0, 15, 0, 1.6f, MBE_FRAME_MARK_PREVIOUS_HELD, 0},
        {250.0, 15, 0, 1.6f, MBE_FRAME_MARK_PREVIOUS_HELD, 0}, /* a second such frame stops the train */
        {170.0, 22, 22, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* and a new one starts */
        {200.0, 18, 18, 16.0f, MBE_FRAME_MARK_NONE, 1},        /* from silence: a new train at the start, */
        {200.0, 18, 18, 16.0f, MBE_FRAME_MARK_NONE, 0},        /* and on from it */
        {160.0, 25, 9, 30.0f, MBE_FRAME_MARK_NONE, 1},         /* from silence at another pitch and level */
        {250.0, 15, 6, 4.0f, MBE_FRAME_MARK_PREVIOUS_HELD, 1}, /* the model before is not held: its own rises */
        {140.0, 26, 26, 16.0f, MBE_FRAME_MARK_NONE, 1},        /* the train it ran is dropped */
        {100.0, 30, 0, 8.0f, MBE_FRAME_MARK_NONE, 0},
        {120.0, 28, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 1}, /* a late onset from silence */
        {250.0, 14, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, 0},
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 1}, /* the marks before do not carry pulses in */
    };

    const size_t frame_count = sizeof(frames) / sizeof(frames[0]);
    struct reference_train train = {0.0, 0, {0.0}};
    mbe_parms last;
    double energy = 0.0, worst = 0.0;

    fresh_stream(&last);
    for (size_t f = 0; f < frame_count; ++f) {
        struct mbe_subframe_model models[3];
        mbe_parms model, cur;
        float out[160];
        double want[160];
        const double start = 160.0 * (double)f;

        set_model(&model, frames[f].hz, frames[f].L, frames[f].voiced, frames[f].amplitude);
        mbe_set_frame_mark(&model, frames[f].mark);
        cur = last;
        cur.w0 = model.w0;
        cur.L = model.L;
        memcpy(cur.Vl, model.Vl, sizeof(cur.Vl));
        memcpy(cur.Ml, model.Ml, sizeof(cur.Ml));
        mbe_set_frame_mark(&cur, frames[f].mark);
        mbe_subframe_own_model(&last, &models[0]);
        const int late = frames[f].mark == MBE_FRAME_MARK_LATE_ONSET && !models[0].any_voiced;
        mbe_subframe_mid_model(&last, &cur, late, &models[1]);
        mbe_subframe_own_model(&cur, &models[2]);

        /* The train starts afresh after a frame with no voiced band, unless that
         * frame lacked a pitch or held the model before it to mid-frame. A
         * frame that starts from silence plays nothing of the last frame's
         * own model and starts a new train, whatever came before. */
        if (frames[f].from_silence) {
            memset(&models[0].voiced, 2, sizeof(models[0].voiced));
            models[0].any_voiced = 0;
            /* A hold continues what was heard: from silence there is nothing
             * to hold, and the frame's own model stands at mid-frame. */
            if (frames[f].mark == MBE_FRAME_MARK_PREVIOUS_HELD) {
                models[1] = models[2];
            }
            train.phase = 0.0;
            train.count = 0;
        } else if (!models[0].any_voiced && mbe_frame_mark(&last) != MBE_FRAME_MARK_PITCH_ABSENT
                   && mbe_frame_mark(&last) != MBE_FRAME_MARK_PREVIOUS_HELD) {
            train.phase = 0.0;
            train.count = 0;
        }
        /* Pulses more than two frames back sound no longer. */
        int kept = 0;
        for (int i = 0; i < train.count; ++i) {
            if (train.time[i] > start - 320.0) {
                train.time[kept++] = train.time[i];
            }
        }
        train.count = kept;
        reference_advance(&train, start, &models[0], &models[1]);
        reference_render(want, start, &train, &models[0], &models[1]);
        reference_advance(&train, start + 80.0, &models[1], &models[2]);
        reference_render(want + 80, start + 80.0, &train, &models[1], &models[2]);

        next_frame_from(out, &last, &model, frames[f].from_silence);
        for (int n = 0; n < 160; ++n) {
            const double difference = fabs((double)out[n] - want[n]);
            worst = (difference > worst) ? difference : worst;
            energy += want[n] * want[n];
            assert(difference < 0.05 + (1.0e-4 * fabs(want[n])));
        }
        assert(fabs(cycles_apart(pulse_phase(&last), train.phase)) < 1e-4);
    }
    /* The comparison saw speech, at a level that makes the tolerance a thousandth of it. */
    assert(sqrt(energy / (160.0 * (double)frame_count)) > 50.0);
    assert(worst < 0.05);
}

/* An unvoiced frame marked as carrying its pitch gives the voiced frame after
 * it a first subframe at the geometric mean of the two pitches; unmarked, the
 * voiced frame's own pitch runs from the start. */
static void
test_unvoiced_pitch_carries(void) {
    float out[160];
    mbe_parms last, voiced;
    for (int carries = 0; carries < 2; ++carries) {
        fresh_stream(&last);
        set_model(&last, 100.0, 30, 0, 4.0f);
        mbe_set_frame_mark(&last, carries ? MBE_FRAME_MARK_PITCH_CARRIES : MBE_FRAME_MARK_NONE);
        set_model(&voiced, 180.0, 20, 20, 16.0f);
        next_frame(out, &last, &voiced);
        const double mid = carries ? sqrt(100.0 * 180.0) : 180.0;
        const double want = cycles_of(mid, 80.0) + cycles_of(0.5 * (mid + 180.0), 80.0);
        assert(fabs(cycles_apart(pulse_phase(&last), want)) < 1e-4);
    }
}

/* Whether the mid-frame model between `before` and `frame` has a voiced band. */
static int
mid_is_voiced(const mbe_parms* before, const mbe_parms* frame) {
    struct mbe_subframe_model mid;
    mbe_subframe_mid_model(before, frame, 0, &mid);
    return mid.any_voiced;
}

/* DVSI's decoders hold the voicing of the frame before a frame that has no
 * voiced band and no pitch through that frame's first half, when that voicing
 * reaches MBE_VOICING_REACH_HZ (the harmonic at or next above it): the pulse
 * train runs on at the old pitch through the whole frame, and a voiced frame
 * after it continues the train rather than starting a new one. */
static void
test_voicing_persists_without_a_pitch(void) {
    float out[160], silent[160];
    mbe_parms last, voiced, unvoiced, after;

    fresh_stream(&last);
    set_model(&voiced, 157.0, 23, 12, 16.0f); /* voiced to 1884 Hz */
    for (int f = 0; f < 3; ++f) {
        next_frame(out, &last, &voiced);
    }
    const double before = pulse_phase(&last);
    set_model(&unvoiced, 8000.0 / 32.25, 14, 0, 9.0f);
    mbe_set_frame_mark(&unvoiced, MBE_FRAME_MARK_PITCH_ABSENT);
    assert(mid_is_voiced(&last, &unvoiced));

    struct mbe_subframe_model mid;
    mbe_subframe_mid_model(&last, &unvoiced, 0, &mid);
    assert(mid.w0 == voiced.w0);
    for (int l = 1; l <= mid.L; ++l) {
        assert(mid.voiced[l] == ((l <= 12) ? 1 : 0));
    }

    mbe_parms held = last, cut = last;
    next_frame(out, &held, &unvoiced);
    assert(fabs(cycles_apart(pulse_phase(&held), before + cycles_of(157.0, 160.0))) < 1e-4);
    /* The voiced part plays through the first half and fades out by the end. */
    double first = 0.0, tail = 0.0;
    for (int n = 0; n < 160; ++n) {
        *((n < 80) ? &first : &tail) += (double)out[n] * (double)out[n];
    }
    assert(first > 0.0 && tail > 0.0 && tail < first);
    assert(fabs((double)out[159]) < 0.05 * sqrt(first / 80.0));

    /* The same frame unmarked ends the voicing at the frame's start. */
    mbe_set_frame_mark(&unvoiced, MBE_FRAME_MARK_NONE);
    assert(!mid_is_voiced(&last, &unvoiced));
    next_frame(silent, &cut, &unvoiced);
    assert(pulse_phase(&cut) == 0.0);
    for (int n = 80; n < 160; ++n) {
        assert(silent[n] == 0.0f);
    }

    /* A voiced frame next: the held train runs on, the cut one starts anew. */
    set_model(&after, 190.0, 19, 19, 16.0f);
    const double held_phase = pulse_phase(&held);
    next_frame(out, &held, &after);
    next_frame(out, &cut, &after);
    assert(fabs(cycles_apart(pulse_phase(&held), held_phase + cycles_of(190.0, 160.0))) < 1e-4);
    assert(fabs(cycles_apart(pulse_phase(&cut), cycles_of(190.0, 160.0))) < 1e-4);
    assert(fabs(cycles_apart(pulse_phase(&held), pulse_phase(&cut))) > 0.01);

    /* The voicing must reach MBE_VOICING_REACH_HZ: three voiced harmonics do
     * at a period of 47.25 samples (507.9 Hz) and do not at 47.75 (502.6 Hz);
     * voiced bands above an unvoiced band around that frequency do not count. */
    mbe_set_frame_mark(&unvoiced, MBE_FRAME_MARK_PITCH_ABSENT);
    set_model(&last, 8000.0 / 47.25, 21, 3, 16.0f);
    assert(mid_is_voiced(&last, &unvoiced));
    set_model(&last, 8000.0 / 47.75, 22, 3, 16.0f);
    assert(!mid_is_voiced(&last, &unvoiced));
    for (int l = 7; l <= 12; ++l) {
        last.Vl[l] = 1;
    }
    assert(!mid_is_voiced(&last, &unvoiced));
    /* A frame with a voiced band is not such a frame, whatever its mark. */
    set_model(&last, 157.0, 23, 12, 16.0f);
    unvoiced.Vl[2] = 1;
    mbe_subframe_mid_model(&last, &unvoiced, 0, &mid);
    assert(mid.voiced[1] != 1);
}

/* log2 of a frame's magnitude at angular frequency w, in the domain the
 * mid-frame model interpolates in: an unvoiced band's magnitude without its
 * factor over sqrt(w0), linear between the two harmonics around w. */
static double
common_log2_at(const mbe_parms* mp, double w) {
    const double unvoiced = log2((double)MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt((double)mp->w0));
    double common[57] = {0.0};
    for (int l = 1; l <= mp->L; ++l) {
        common[l] = log2((double)mp->Ml[l]) - ((mp->Vl[l] == 0) ? unvoiced : 0.0);
    }
    const double position = w / (double)mp->w0;
    if (position < 1.0) {
        return common[1];
    }
    if (position >= (double)mp->L) {
        return common[mp->L];
    }
    const int below = (int)position;
    return common[below] + ((position - below) * (common[below + 1] - common[below]));
}

/* The mid-frame model: pitch, harmonic count, voicing and magnitudes. */
static void
test_mid_model(void) {
    const double k = (double)MBE_UNVOICED_MAGNITUDE_FACTOR;
    mbe_parms prev, cur, enhanced;
    struct mbe_subframe_model mid;
    mbe_initMbeParms(&cur, &prev, &enhanced);

    /* Both voiced: the geometric-mean pitch, this frame's voicing at each
     * band's frequency, magnitudes the geometric mean of the two frames'. */
    set_model(&prev, 150.0, 24, 24, 10.0f);
    set_model(&cur, 180.0, 20, 10, 20.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, &mid);
    const double w0 = sqrt((double)prev.w0 * (double)cur.w0);
    assert(fabs((double)mid.w0 - w0) < 1e-6);
    assert(mid.L == 22 && mid.any_voiced);
    for (int l = 1; l <= mid.L; ++l) {
        const double w = (double)l * w0;
        const int harmonic = (int)((w / (double)cur.w0) + 0.5);
        const int voiced = harmonic <= 10;
        assert(mid.voiced[l] == voiced);
        /* An unvoiced band's magnitude carries its factor over sqrt(w0): out
         * of each frame's before the mean, into the mid-frame model's after. */
        const double mean = 0.5 * (common_log2_at(&prev, w) + common_log2_at(&cur, w));
        const double want = exp2(mean + (voiced ? 0.0 : log2(k / sqrt(w0))));
        assert(fabs((double)mid.magnitude[l] / want - 1.0) < 1e-4);
    }
    /* Away from the edge of this frame's voicing: a band voiced in both frames,
     * and one voiced before and unvoiced now. */
    assert(fabs((double)mid.magnitude[5] / sqrt(200.0) - 1.0) < 1e-4);
    assert(fabs((double)mid.magnitude[20] / (sqrt(10.0 * 20.0 * sqrt((double)cur.w0) / k) * k / sqrt(w0)) - 1.0)
           < 1e-4);

    /* A magnitude is read at the band's frequency, log-linear between harmonics. */
    for (int l = 1; l <= prev.L; ++l) {
        prev.Ml[l] = 10.0f * (float)l;
    }
    mbe_subframe_mid_model(&prev, &cur, 0, &mid);
    const double position = 3.0 * w0 / (double)prev.w0;
    const int below = (int)position;
    const double at = exp2(log2(10.0 * below) + (position - below) * (log2(10.0 * (below + 1)) - log2(10.0 * below)));
    assert(fabs((double)mid.magnitude[3] / sqrt(at * 20.0) - 1.0) < 1e-4);

    /* Only one frame with a voiced band: its pitch. Neither: between them. */
    set_model(&prev, 150.0, 24, 0, 10.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, &mid);
    assert(mid.w0 == cur.w0 && mid.any_voiced);
    set_model(&prev, 150.0, 24, 24, 10.0f);
    set_model(&cur, 180.0, 20, 0, 20.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, &mid);
    assert(mid.w0 == prev.w0 && !mid.any_voiced);
    set_model(&prev, 150.0, 24, 0, 10.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, &mid);
    assert(fabs((double)mid.w0 - w0) < 1e-6 && !mid.any_voiced);
    assert(fabs((double)mid.magnitude[5] / sqrt(200.0) - 1.0) < 1e-4);

    /* A late onset has no voiced band at mid-frame; its bands are noise. */
    set_model(&cur, 180.0, 20, 20, 20.0f);
    mbe_subframe_mid_model(&prev, &cur, 1, &mid);
    assert(!mid.any_voiced && mid.voiced[1] == 0 && mid.magnitude[1] > 0.0f);

    /* A band silent in either frame is silent at mid-frame. */
    prev.Ml[4] = 0.0f;
    prev.Ml[5] = 0.0f;
    mbe_subframe_mid_model(&prev, &cur, 0, &mid);
    int silent = 0;
    for (int l = 1; l <= mid.L; ++l) {
        silent += mid.magnitude[l] == 0.0f;
    }
    assert(silent >= 1 && silent <= 2);

    /* A frame's own model. */
    struct mbe_subframe_model own;
    cur.Vl[3] = 0;
    cur.Vl[4] = 7; /* neither voiced nor unvoiced: renders nothing */
    mbe_subframe_own_model(&cur, &own);
    assert(own.w0 == cur.w0 && own.L == 20 && own.any_voiced);
    assert(own.voiced[2] == 1 && own.voiced[3] == 0 && own.voiced[4] > 1 && own.magnitude[2] == 20.0f);
}

static int
same_model(const struct mbe_subframe_model* a, const struct mbe_subframe_model* b) {
    return same_bits(&a->w0, &b->w0, 1) && a->L == b->L && a->any_voiced == b->any_voiced
           && memcmp(a->voiced, b->voiced, sizeof(a->voiced)) == 0 && same_bits(a->magnitude, b->magnitude, 57)
           && same_bits(a->phase, b->phase, 57);
}

/* A frame marked MBE_FRAME_MARK_PREVIOUS_HELD (a rate-33 silence frame, as
 * DVSI's decoder plays it) has the previous frame's own model at mid-frame,
 * whatever the voicing of either frame and whatever the late-onset flag; the
 * marks P25 and D-STAR set keep their rules. The pulse train runs on through
 * such a frame at the held pitch, and rendering and advancing it leave the
 * same state. */
static void
test_previous_model_held(void) {
    static const int voiced_before[3] = {24, 3, 0};
    float out[160];
    mbe_parms prev, cur, enhanced;
    struct mbe_subframe_model mid, own;
    mbe_initMbeParms(&cur, &prev, &enhanced);

    for (int k = 0; k < 3; ++k) {
        set_model(&prev, 150.0, 24, voiced_before[k], 10.0f);
        for (int l = 1; l <= 24; ++l) {
            prev.Ml[l] = 4.0f + (float)l;
        }
        mbe_subframe_own_model(&prev, &own);
        for (int cur_voiced = 0; cur_voiced <= 15; cur_voiced += 15) {
            set_model(&cur, 250.0, 15, cur_voiced, 1.6f);
            mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PREVIOUS_HELD);
            for (int late = 0; late < 2; ++late) {
                mbe_subframe_mid_model(&prev, &cur, late, &mid);
                assert(same_model(&mid, &own));
            }
        }
        /* The same frame marked as P25's and D-STAR's decoders mark one with
         * no voiced band and no pitch meets the two frames' magnitudes. */
        set_model(&cur, 250.0, 15, 0, 1.6f);
        mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PITCH_ABSENT);
        mbe_subframe_mid_model(&prev, &cur, 0, &mid);
        assert(!same_model(&mid, &own));
    }

    /* After three voiced frames: one held frame runs the train on at 157 Hz. */
    mbe_parms last, voiced, held;
    fresh_stream(&last);
    set_model(&voiced, 157.0, 23, 12, 16.0f);
    for (int f = 0; f < 3; ++f) {
        next_frame(out, &last, &voiced);
    }
    const double before = pulse_phase(&last);
    set_model(&held, 250.0, 15, 0, 1.6f);
    mbe_set_frame_mark(&held, MBE_FRAME_MARK_PREVIOUS_HELD);

    mbe_parms rendered = last, advanced = last;
    rendered.w0 = advanced.w0 = held.w0;
    rendered.L = advanced.L = held.L;
    memcpy(rendered.Vl, held.Vl, sizeof(held.Vl));
    memcpy(advanced.Vl, held.Vl, sizeof(held.Vl));
    memcpy(rendered.Ml, held.Ml, sizeof(held.Ml));
    memcpy(advanced.Ml, held.Ml, sizeof(held.Ml));
    mbe_set_frame_mark(&rendered, MBE_FRAME_MARK_PREVIOUS_HELD);
    mbe_set_frame_mark(&advanced, MBE_FRAME_MARK_PREVIOUS_HELD);
    memset(out, 0, sizeof(out));
    mbe_synthesize_subframes(out, &rendered, &last, NULL, 1.0f, 0);
    mbe_advance_subframes(&advanced, &last, 0);
    assert(same_bits(rendered.PSIl, advanced.PSIl, 57));
    assert(same_bits(rendered.PHIl, advanced.PHIl, 57));
    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
    assert(fabs(cycles_apart(pulse_phase(&rendered), before + cycles_of(157.0, 160.0))) < 1e-4);
    assert(rendered.noiseOverlap[0] >= 1.0f);

    /* A voiced frame next continues that train. */
    mbe_parms after;
    set_model(&after, 190.0, 19, 19, 16.0f);
    const double held_phase = pulse_phase(&rendered);
    next_frame(out, &rendered, &after);
    assert(fabs(cycles_apart(pulse_phase(&rendered), held_phase + cycles_of(190.0, 160.0))) < 1e-4);
}

/* Between two models a subframe crossfades linearly. A model with the voicing
 * and pitch of its neighbour but no amplitude keeps the pulse train and makes
 * the mid-frame model silent too (a band silent in either frame is), so the
 * steady frame fades in over its second subframe alone, and out over the first. */
static void
test_subframes_crossfade_linearly(void) {
    float flat[160], in[160], out[160];
    mbe_parms last, model, silent;
    fresh_stream(&last);
    set_model(&model, 187.5, 18, 18, 16.0f);
    for (int f = 0; f < 3; ++f) {
        next_frame(flat, &last, &model);
    }
    mbe_parms from_silence = last, to_silence = last;
    memset(from_silence.Ml, 0, sizeof(from_silence.Ml));
    silent = model;
    memset(silent.Ml, 0, sizeof(silent.Ml));
    next_frame(flat, &last, &model);
    next_frame(in, &from_silence, &model);
    next_frame(out, &to_silence, &silent);
    double peak = 0.0;
    for (int n = 0; n < 160; ++n) {
        peak = fmax(peak, fabs((double)flat[n]));
    }
    assert(peak > 1.0);
    for (int n = 0; n < 160; ++n) {
        const double rising = (n < 80) ? 0.0 : (double)(n - 80) / 80.0;
        const double falling = (n < 80) ? 1.0 - ((double)n / 80.0) : 0.0;
        assert(fabs((double)in[n] - rising * (double)flat[n]) < 1e-4 * peak);
        assert(fabs((double)out[n] - falling * (double)flat[n]) < 1e-4 * peak);
    }
    /* The train ran on through both. */
    assert(same_bits(from_silence.noiseOverlap, last.noiseOverlap, 96));
    assert(same_bits(to_silence.noiseOverlap, last.noiseOverlap, 96));
}

/* A frame the core mutes moves the stream's state on exactly as a rendered
 * frame does, and leaves nothing to overlap the next. */
static void
test_muted_frame_advances(void) {
    float out[160];
    mbe_parms last, model;
    fresh_stream(&last);
    set_model(&model, 210.0, 17, 9, 16.0f);
    next_frame(out, &last, &model);
    set_model(&model, 190.0, 19, 11, 12.0f);

    mbe_parms rendered = last, advanced = last;
    rendered.w0 = advanced.w0 = model.w0;
    rendered.L = advanced.L = model.L;
    memcpy(rendered.Vl, model.Vl, sizeof(model.Vl));
    memcpy(advanced.Vl, model.Vl, sizeof(model.Vl));
    memcpy(rendered.Ml, model.Ml, sizeof(model.Ml));
    memcpy(advanced.Ml, model.Ml, sizeof(model.Ml));
    memset(out, 0, sizeof(out));
    mbe_synthesize_subframes(out, &rendered, &last, NULL, 1.0f, 0);
    mbe_advance_subframes(&advanced, &last, 0);
    assert(same_bits(rendered.PSIl, advanced.PSIl, 57));
    assert(same_bits(rendered.PHIl, advanced.PHIl, 57));
    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
    assert(rendered.noiseSeed == advanced.noiseSeed);
    for (int n = 0; n < 256; ++n) {
        assert(advanced.previousUw[n] == 0.0f);
    }

    /* Through the core: a frame muted by the error rate leaves the state a
     * rendered frame leaves, so the frame after both plays the same voiced part. */
    mbe_parms muted = last, enh_a = last, enh_b = last;
    muted.w0 = model.w0;
    muted.L = model.L;
    memcpy(muted.Vl, model.Vl, sizeof(model.Vl));
    memcpy(muted.Ml, model.Ml, sizeof(model.Ml));
    mbe_parms heard = muted;
    muted.errorRate = 1.0f;
    muted.mutingThreshold = 0.0875f;
    heard.errorRate = 0.0f;
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &muted, &enh_a, 0.0f, MBE_MUTE_NOISE_SPEC, MBE_VOICED_SYNTHESIS_DVSI,
                                       MBE_AMPLITUDE_LIMIT_NONE);
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &heard, &enh_b, 0.0f, MBE_MUTE_NOISE_SPEC, MBE_VOICED_SYNTHESIS_DVSI,
                                       MBE_AMPLITUDE_LIMIT_NONE);
    assert(same_bits(muted.noiseOverlap, heard.noiseOverlap, 96));
    assert(same_bits(muted.PSIl, heard.PSIl, 57));
    assert(muted.noiseSeed == heard.noiseSeed);
}

/* The frame a parameter set holds, given a model and mark. */
static void
take_model(mbe_parms* mp, const mbe_parms* model, enum mbe_frame_mark mark) {
    mp->w0 = model->w0;
    mp->L = model->L;
    memcpy(mp->Vl, model->Vl, sizeof(mp->Vl));
    memcpy(mp->Ml, model->Ml, sizeof(mp->Ml));
    mbe_set_frame_mark(mp, mark);
}

static float
peak_of(const float* x) {
    float m = 0.0f;
    for (int n = 0; n < 160; ++n) {
        m = fmaxf(m, fabsf(x[n]));
    }
    return m;
}

static void
held_from_silence(const mbe_parms* loud, double loud_peak) {
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    float held_out[160], plain_out[160], steady_out[160], continued_out[160];
    mbe_parms silence, own, held, plain, steady, continued;
    assert(plan != NULL);
    /* Quiet, as a silence frame is against speech. */
    set_model(&silence, 250.0, 15, 0, 0.2f);

    /* The silence frame from silence after the loud model, marked. */
    held = *loud;
    held.noiseSeed = 4321.0f;
    take_model(&held, &silence, MBE_FRAME_MARK_PREVIOUS_HELD);
    memset(held_out, 0, sizeof(held_out));
    mbe_synthesize_subframes(held_out, &held, loud, plan, 1.0f, 1);

    /* The same frame unmarked, from silence, after its own model. */
    own = *loud;
    take_model(&own, &silence, MBE_FRAME_MARK_NONE);
    plain = *loud;
    plain.noiseSeed = 4321.0f;
    take_model(&plain, &silence, MBE_FRAME_MARK_NONE);
    memset(plain_out, 0, sizeof(plain_out));
    mbe_synthesize_subframes(plain_out, &plain, &own, plan, 1.0f, 1);

    /* A steady silence frame: after its own model, continuing. */
    steady = plain;
    steady.noiseSeed = 4321.0f;
    memset(steady_out, 0, sizeof(steady_out));
    mbe_synthesize_subframes(steady_out, &steady, &plain, plan, 1.0f, 0);

    /* Continuing (no option) the mark holds the loud model to mid-frame. */
    continued = *loud;
    continued.noiseSeed = 4321.0f;
    take_model(&continued, &silence, MBE_FRAME_MARK_PREVIOUS_HELD);
    memset(continued_out, 0, sizeof(continued_out));
    mbe_synthesize_subframes(continued_out, &continued, loud, plan, 1.0f, 0);
    mbe_fft_plan_free(plan);

    const double silence_peak = (double)peak_of(steady_out);
    assert(silence_peak > 0.1 && silence_peak * 20.0 < loud_peak);
    for (int n = 0; n < 160; ++n) {
        assert(fabs((double)held_out[n] - (double)plain_out[n]) < 1e-4 * (1.0 + silence_peak));
    }
    assert(peak_of(held_out) > 0.0f);
    assert((double)peak_of(held_out) < 3.0 * silence_peak);
    assert((double)peak_of(continued_out) > 20.0 * silence_peak);
    assert(same_bits(held.noiseOverlap, plain.noiseOverlap, 96));
    assert(same_bits(&held.noiseSeed, &plain.noiseSeed, 1));
}

/*
 * Starting from silence (the first frame heard after a rate-33 mute or tone):
 *  - a steady all-voiced pair renders, over the first subframe, the ordinary
 *    rendering of the same pair started with an empty pulse train times
 *    n / 80 (the previous model's falling half is gone), and the same after it;
 *  - pulses and a phase the previous state carries neither sound nor move the
 *    new train, and no unvoiced overlap comes from it, whatever it holds;
 *  - a frame marked MBE_FRAME_MARK_PREVIOUS_HELD (a rate-33 silence frame)
 *    does not hold a previous model that was not heard: from silence it
 *    rises to its own model and plays it, as the same frame unmarked does
 *    from silence after its own model, and stays within a small multiple of
 *    its steady level after loud voiced speech; continuing, the mark still
 *    holds the previous model;
 *  - advancing (a frame the core mutes) leaves the state rendering leaves,
 *    through the core as well.
 */
static void
test_from_silence(void) {
    float steady[160], rising[160], other[160], out_scratch[160];
    mbe_parms last, model;
    fresh_stream(&last);
    set_model(&model, 187.5, 18, 18, 16.0f);
    for (int f = 0; f < 3; ++f) {
        next_frame(steady, &last, &model);
    }
    /* The pulses and phase the stream carries now, and the same state without them. */
    assert(last.noiseOverlap[0] >= 1.0f && pulse_phase(&last) > 0.0);
    mbe_parms reset = last, carried = last, emptied = last;
    memset(reset.noiseOverlap, 0, sizeof(reset.noiseOverlap));
    memset(emptied.noiseOverlap, 0, sizeof(emptied.noiseOverlap));
    mbe_subframe_state_mark(&emptied);
    const mbe_parms no_train = emptied;
    next_frame_from(steady, &reset, &model, 0);
    next_frame_from(rising, &carried, &model, 1);
    next_frame_from(other, &emptied, &model, 1);
    double peak = 0.0;
    for (int n = 0; n < 160; ++n) {
        peak = fmax(peak, fabs((double)steady[n]));
    }
    assert(peak > 1.0);
    for (int n = 0; n < 160; ++n) {
        const double gain = (n < 80) ? (double)n / 80.0 : 1.0;
        assert(fabs((double)rising[n] - gain * (double)steady[n]) < 1e-4 * peak);
    }
    assert(fabsf(rising[0]) < 1e-6f);
    assert(same_bits(rising, other, 160));
    assert(same_bits(carried.noiseOverlap, emptied.noiseOverlap, 96));
    assert(same_bits(carried.PSIl, emptied.PSIl, 57));
    assert(fabs(cycles_apart(pulse_phase(&carried), cycles_of(187.5, 160.0))) < 1e-5);

    /* The unvoiced part: an overlap left in the previous state does not play. */
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    mbe_parms unvoiced, with_overlap, without_overlap;
    fresh_stream(&unvoiced);
    set_model(&unvoiced, 150.0, 24, 6, 12.0f);
    unvoiced.noiseSeed = 1234.0f;
    with_overlap = unvoiced;
    without_overlap = unvoiced;
    for (int n = 0; n < 256; ++n) {
        with_overlap.previousUw[n] = 500.0f * (float)((n % 7) - 3);
    }
    mbe_subframe_state_mark(&with_overlap);
    mbe_subframe_state_mark(&without_overlap);
    mbe_parms a = unvoiced, b = unvoiced;
    memset(rising, 0, sizeof(rising));
    memset(other, 0, sizeof(other));
    mbe_synthesize_subframes(rising, &a, &with_overlap, plan, 1.0f, 1);
    mbe_synthesize_subframes(other, &b, &without_overlap, plan, 1.0f, 1);
    assert(same_bits(rising, other, 160));
    assert(same_bits(&a.noiseSeed, &b.noiseSeed, 1) && !same_bits(&a.noiseSeed, &unvoiced.noiseSeed, 1));
    assert(same_bits(a.previousUw, b.previousUw, 256));
    /* Continuing, the same overlap plays. */
    a = unvoiced;
    memset(other, 0, sizeof(other));
    mbe_synthesize_subframes(other, &a, &with_overlap, plan, 1.0f, 0);
    assert(!same_bits(rising, other, 160));
    mbe_fft_plan_free(plan);

    /* A silence frame after loud voiced speech (no_train): from silence it
     * rises to its own model, exactly as the unmarked frame does from silence
     * after its own model (whose ordinary mid-frame model is its own), from the
     * same noise state; its level stays that of a steady silence frame. */
    held_from_silence(&no_train, peak);

    /* Advancing and rendering leave the same state, from silence too. */
    mbe_parms rendered = carried, advanced = carried;
    set_model(&model, 210.0, 17, 9, 16.0f);
    rendered.w0 = advanced.w0 = model.w0;
    rendered.L = advanced.L = model.L;
    memcpy(rendered.Vl, model.Vl, sizeof(model.Vl));
    memcpy(advanced.Vl, model.Vl, sizeof(model.Vl));
    memcpy(rendered.Ml, model.Ml, sizeof(model.Ml));
    memcpy(advanced.Ml, model.Ml, sizeof(model.Ml));
    mbe_parms continued = advanced;
    memset(out_scratch, 0, sizeof(out_scratch));
    mbe_synthesize_subframes(out_scratch, &rendered, &carried, NULL, 1.0f, 1);
    mbe_advance_subframes(&advanced, &carried, 1);
    mbe_advance_subframes(&continued, &carried, 0);
    assert(same_bits(rendered.PSIl, advanced.PSIl, 57));
    assert(same_bits(rendered.PHIl, advanced.PHIl, 57));
    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
    assert(same_bits(&rendered.noiseSeed, &advanced.noiseSeed, 1));
    assert(!same_bits(continued.noiseOverlap, advanced.noiseOverlap, 96));

    /* Through the core: a frame it mutes moves the state on as the rendered frame does. */
    mbe_parms muted = carried, heard = carried, enh_a = carried, enh_b = carried;
    muted.w0 = heard.w0 = model.w0;
    muted.L = heard.L = model.L;
    memcpy(muted.Vl, model.Vl, sizeof(model.Vl));
    memcpy(heard.Vl, model.Vl, sizeof(model.Vl));
    memcpy(muted.Ml, model.Ml, sizeof(model.Ml));
    memcpy(heard.Ml, model.Ml, sizeof(model.Ml));
    muted.errorRate = 1.0f;
    muted.mutingThreshold = 0.0875f;
    heard.errorRate = 0.0f;
    mbe_synthesizeSpeechWithPreEnhRm0Startf(out_scratch, &muted, &enh_a, 0.0f, MBE_MUTE_NOISE_SPEC,
                                            MBE_VOICED_SYNTHESIS_DVSI, MBE_AMPLITUDE_LIMIT_NONE,
                                            MBE_SYNTHESIS_FROM_SILENCE, MBE_OUTPUT_LIMITED);
    mbe_synthesizeSpeechWithPreEnhRm0Startf(out_scratch, &heard, &enh_b, 0.0f, MBE_MUTE_NOISE_SPEC,
                                            MBE_VOICED_SYNTHESIS_DVSI, MBE_AMPLITUDE_LIMIT_NONE,
                                            MBE_SYNTHESIS_FROM_SILENCE, MBE_OUTPUT_LIMITED);
    assert(same_bits(muted.noiseOverlap, heard.noiseOverlap, 96));
    assert(same_bits(muted.PSIl, heard.PSIl, 57));
    assert(same_bits(&muted.noiseSeed, &heard.noiseSeed, 1));
}

/*
 * The output-limit option of the synthesis core: a frame left unlimited is
 * the limited frame before its limit, sample for sample (limiting it gives the
 * limited frame's bits), for a voiced, an unvoiced and a from-silence frame,
 * loud enough to pass the limit; and the state the two leave is the same.
 */
static void
test_output_limit_option(void) {
    static const struct {
        int voiced;
        enum mbe_synthesis_start start;
    } cases[] = {
        {18, MBE_SYNTHESIS_CONTINUES},
        {0, MBE_SYNTHESIS_CONTINUES},
        {18, MBE_SYNTHESIS_FROM_SILENCE},
    };

    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); ++k) {
        float out[160], limited[160], unlimited[160];
        mbe_parms last, model;
        fresh_stream(&last);
        set_model(&model, 187.5, 18, cases[k].voiced, 3000.0f);
        next_frame(out, &last, &model);
        last.noiseSeed = 777.0f;
        mbe_parms a = last, b = last, prev_a = last, prev_b = last;
        mbe_synthesizeSpeechWithPreEnhRm0Startf(limited, &a, &prev_a, 0.0f, MBE_MUTE_NOISE_SPEC,
                                                MBE_VOICED_SYNTHESIS_DVSI, MBE_AMPLITUDE_LIMIT_NONE, cases[k].start,
                                                MBE_OUTPUT_LIMITED);
        mbe_synthesizeSpeechWithPreEnhRm0Startf(unlimited, &b, &prev_b, 0.0f, MBE_MUTE_NOISE_SPEC,
                                                MBE_VOICED_SYNTHESIS_DVSI, MBE_AMPLITUDE_LIMIT_NONE, cases[k].start,
                                                MBE_OUTPUT_UNLIMITED);
        const float limit = (32767.0f * 0.95f) / 7.0f;
        int past = 0;
        for (int n = 0; n < 160; ++n) {
            past |= fabsf(unlimited[n]) > limit;
        }
        assert(past);
        assert(!same_bits(limited, unlimited, 160));
        mbe_clipFloatBuffer(unlimited, 160);
        assert(same_bits(limited, unlimited, 160));
        assert(same_bits(a.Ml, b.Ml, 57) && same_bits(a.PSIl, b.PSIl, 57) && same_bits(a.PHIl, b.PHIl, 57));
        assert(same_bits(a.noiseOverlap, b.noiseOverlap, 96) && same_bits(a.previousUw, b.previousUw, 256));
        assert(same_bits(&a.noiseSeed, &b.noiseSeed, 1) && same_bits(&a.localEnergy, &b.localEnergy, 1));
        assert(memcmp(a.Vl, b.Vl, sizeof(a.Vl)) == 0 && a.amplitudeThreshold == b.amplitudeThreshold);
    }
}

/* Pulse state that no frame stores is dropped, and never read as an index or
 * a count: a stream handed such a state plays as one starting afresh. */
static void
test_pulse_state_is_validated(void) {
    static const uint32_t counts[] = {0x42040000u /* 33 */, 0x7FC00000u,           0x7F800000u, 0xFF800000u,
                                      0xBF800000u /* -1 */, 0x4E6E6B28u /* 1e9 */, 0x80000000u, 0x7F7FFFFFu};
    static const uint32_t values[] = {
        0x7FC00000u, 0x7F800000u, 0xFF800000u, 0x3F800000u /* 1 */, 0xC3210000u /* -161 */,
        0x7F7FFFFFu, 0xFF7FFFFFu, 0x40A00000u /* 5 */};
    float reference[160], out[160];
    mbe_parms base, model, last;
    fresh_stream(&base);
    set_model(&base, 200.0, 18, 18, 16.0f);
    set_model(&model, 200.0, 18, 18, 16.0f);
    last = base;
    next_frame(reference, &last, &model);
    const mbe_parms fresh = last;

    for (size_t i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i) {
        last = base;
        memcpy(&last.noiseOverlap[0], &counts[i], sizeof(counts[i]));
        for (int k = 1; k < 96; ++k) {
            last.noiseOverlap[k] = -20.0f * (float)k;
        }
        last.noiseOverlap[PHASE_SLOT] = 0.5f;
        mbe_subframe_state_mark(&last);
        next_frame(out, &last, &model);
        assert(all_finite(out, 160));
        assert(same_bits(out, reference, 160));
        assert(same_bits(last.noiseOverlap, fresh.noiseOverlap, 96));
    }
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        last = base;
        last.noiseOverlap[0] = 32.0f;
        for (int k = 1; k < 96; ++k) {
            memcpy(&last.noiseOverlap[k], &values[i], sizeof(values[i]));
        }
        mbe_subframe_state_mark(&last);
        next_frame(out, &last, &model);
        assert(all_finite(out, 160));
        assert(same_bits(out, reference, 160));
    }
    /* A valid state is used: a phase of half a cycle moves the pulses. */
    last = base;
    last.noiseOverlap[PHASE_SLOT] = 0.5f;
    mbe_subframe_state_mark(&last);
    next_frame(out, &last, &model);
    assert(!same_bits(out, reference, 160));
    assert(fabs(cycles_apart(pulse_phase(&last), 0.5 + cycles_of(200.0, 160.0))) < 1e-5);
    /* Without this rendering's mark the same values are the other rendering's noise: not a pulse train. */
    last = base;
    last.noiseOverlap[PHASE_SLOT] = 0.5f;
    next_frame(out, &last, &model);
    assert(same_bits(out, reference, 160));
    assert(mbe_subframe_state_is_marked(&last));
}

/* Models no decoder produces render finite samples within the frame, and
 * leave a state the next frame can use. */
static void
test_invalid_models(void) {
    static const uint32_t w0_bits[] = {0x7FC00000u,          0x7F800000u, 0xFF800000u, 0x00000000u,
                                       0x80000000u,          0xBF000000u, 0x00000001u, 0x7F7FFFFFu,
                                       0x40490FDBu /* pi */, 0x40400000u /* 3 */};
    static const uint32_t ml_bits[] = {0x7FC00000u, 0x7F800000u, 0xFF800000u, 0xC1200000u, 0x00000001u, 0x7F7FFFFFu};
    mbe_parms good, last;
    float guarded[160 + 32];

    set_model(&good, 200.0, 18, 18, 16.0f);
    for (size_t i = 0; i < sizeof(w0_bits) / sizeof(w0_bits[0]); ++i) {
        for (int which = 0; which < 3; ++which) {
            mbe_parms model = good, bad = good;
            memcpy(&bad.w0, &w0_bits[i], sizeof(bad.w0));
            fresh_stream(&last);
            next_frame(guarded + 16, &last, (which == 1) ? &bad : &good);
            if (which == 2) {
                memcpy(&last.w0, &w0_bits[i], sizeof(last.w0));
            }
            for (int n = 0; n < 16; ++n) {
                guarded[n] = 1234.5f;
                guarded[176 + n] = -1234.5f;
            }
            mbe_parms cur = last;
            const mbe_parms* next = (which == 0) ? &bad : &model;
            cur.w0 = next->w0;
            cur.L = next->L;
            memcpy(cur.Vl, next->Vl, sizeof(cur.Vl));
            memcpy(cur.Ml, next->Ml, sizeof(cur.Ml));
            memset(guarded + 16, 0, 160 * sizeof(float));
            mbe_synthesize_subframes(guarded + 16, &cur, &last, NULL, 1.0f, 0);
            for (int n = 0; n < 16; ++n) {
                assert(guarded[n] == 1234.5f && guarded[176 + n] == -1234.5f);
            }
            assert(all_finite(guarded + 16, 160));
            assert(all_finite(cur.noiseOverlap, 96));
            /* The phases the frame leaves stay finite too, for the other rendering to read. */
            assert(all_finite(&cur.PSIl[1], 56) && all_finite(&cur.PHIl[1], 56));
            last = cur;
            next_frame(guarded + 16, &last, &good);
            assert(all_finite(last.noiseOverlap, 96));
        }
    }
    for (size_t i = 0; i < sizeof(ml_bits) / sizeof(ml_bits[0]); ++i) {
        mbe_parms bad = good;
        for (int l = 1; l <= 56; l += 2) {
            memcpy(&bad.Ml[l], &ml_bits[i], sizeof(bad.Ml[l]));
        }
        fresh_stream(&last);
        next_frame(guarded + 16, &last, &good);
        next_frame(guarded + 16, &last, &bad);
        assert(all_finite(guarded + 16, 160));
        next_frame(guarded + 16, &last, &good);
        assert(all_finite(last.noiseOverlap, 96));
    }
    /* Magnitudes that meet past what a float holds: a voiced band at the
     * largest float beside an unvoiced one, at a pitch whose unvoiced factor is
     * huge. The mid-frame model's stay finite and within the model's limit. */
    {
        static const uint32_t tiny_pitch = 0x00800000u; /* FLT_MIN */
        static const uint32_t huge = 0x7F7FFFFFu;       /* FLT_MAX */
        struct mbe_subframe_model model;
        mbe_parms before = good, after = good;
        memcpy(&before.w0, &tiny_pitch, sizeof(before.w0));
        for (int l = 1; l <= 56; ++l) {
            memcpy(&before.Ml[l], &huge, sizeof(before.Ml[l]));
            after.Vl[l] = 0;
            after.Ml[l] = 1000.0f;
        }
        after.w0 = 0.2f;
        for (int order = 0; order < 2; ++order) {
#if MID_MODEL_OVERFLOW_IS_VISIBLE
            /* Limited before it narrows to float, not after: no conversion overflows on the way. */
            feclearexcept(FE_ALL_EXCEPT);
            mbe_subframe_mid_model(order ? &after : &before, order ? &before : &after, 0, &model);
            assert(fetestexcept(FE_OVERFLOW) == 0);
#endif
            mbe_subframe_mid_model(order ? &after : &before, order ? &before : &after, 0, &model);
            assert(all_finite(&model.magnitude[1], 56));
            for (int l = 1; l <= model.L; ++l) {
                assert(model.magnitude[l] >= 0.0f && model.magnitude[l] <= 1.0e30f);
            }
            fresh_stream(&last);
            next_frame(guarded + 16, &last, order ? &after : &before);
            next_frame(guarded + 16, &last, order ? &before : &after);
            assert(all_finite(guarded + 16, 160));
        }
    }
    /* A harmonic count outside 1..56, in either frame: nothing is rendered,
     * the state stays as it is, and the models have no band. */
    static const int counts[] = {0, -1, 57, 1000, INT_MAX, INT_MIN};
    set_model(&good, 200.0, 18, 18, 16.0f);
    for (size_t i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i) {
        for (int which = 0; which < 2; ++which) {
            struct mbe_subframe_model model;
            fresh_stream(&last);
            next_frame(guarded + 16, &last, &good);
            mbe_parms cur = last, prev = last;
            (which ? &prev : &cur)->L = counts[i];
            unsigned char before[sizeof(mbe_parms)], after[sizeof(mbe_parms)];
            memcpy(before, &cur, sizeof(before));
            for (int n = 0; n < 192; ++n) {
                guarded[n] = 7.0f;
            }
            mbe_synthesize_subframes(guarded + 16, &cur, &prev, NULL, 1.0f, 0);
            for (int n = 0; n < 192; ++n) {
                assert(guarded[n] == 7.0f);
            }
            memcpy(after, &cur, sizeof(after));
            assert(memcmp(after, before, sizeof(after)) == 0);
            mbe_advance_subframes(&cur, &prev, 0);
            memcpy(after, &cur, sizeof(after));
            assert(memcmp(after, before, sizeof(after)) == 0);
            mbe_subframe_mid_model(&prev, &cur, 0, &model);
            assert(model.L == 0 && !model.any_voiced);
            mbe_subframe_own_model(which ? &prev : &cur, &model);
            assert(model.L == 0 && !model.any_voiced);
        }
    }

    /* A period of a few samples: an interval's pulses are capped. */
    fresh_stream(&last);
    set_model(&good, 3900.0, 9, 9, 1.0f);
    for (int f = 0; f < 3; ++f) {
        next_frame(guarded + 16, &last, &good);
        assert(all_finite(guarded + 16, 160));
        assert(last.noiseOverlap[0] <= (float)MBE_PULSE_HISTORY);
    }
}

static void
set_imbe_b0(char imbe_d[88], int b0) {
    static const int index[8] = {0, 1, 2, 3, 4, 5, 85, 86};
    for (int i = 0; i < 8; ++i) {
        imbe_d[index[i]] = (char)((b0 >> (7 - i)) & 1);
    }
}

static void
set_ambe2400(char ambe_d[49], int b0, int b1) {
    memset(ambe_d, 0, 49);
    for (int i = 0; i < 6; ++i) {
        ambe_d[i] = (char)((b0 >> (6 - i)) & 1);
    }
    ambe_d[48] = (char)(b0 & 1);
    for (int i = 0; i < 4; ++i) {
        ambe_d[38 + i] = (char)((b1 >> (3 - i)) & 1);
    }
}

/* A pitch too small to divide by is no pitch. With denormals flushed, a
 * subnormal w0 is zero to the arithmetic, and a position among the harmonics
 * computed from it is not a number: an optimized build would carry that into
 * an index. Such a model has no band, whatever the environment. */
static void
test_subnormal_pitch_is_no_pitch(void) {
    static const uint32_t w0_bits[] = {0x00000001u, 0x00000100u, 0x007FFFFFu};
    for (size_t i = 0; i < sizeof(w0_bits) / sizeof(w0_bits[0]); ++i) {
        mbe_parms voiced, unvoiced, start;
        struct mbe_subframe_model model;
        float guarded[160 + 32];

        set_model(&voiced, 200.0, 18, 18, 16.0f);
        set_model(&unvoiced, 200.0, 18, 0, 16.0f);
        memcpy(&voiced.w0, &w0_bits[i], sizeof(voiced.w0));
        mbe_subframe_own_model(&voiced, &model);
        assert(model.L == 18 && !model.any_voiced);
        for (int l = 1; l <= model.L; ++l) {
            assert(model.voiced[l] == 2 && model.magnitude[l] == 0.0f);
        }
        /* The same pitch, a magnitude as small: silence. */
        set_model(&start, 200.0, 18, 18, 16.0f);
        memcpy(&start.Ml[3], &w0_bits[i], sizeof(start.Ml[3]));
        mbe_subframe_own_model(&start, &model);
        assert(model.any_voiced && model.magnitude[3] == 0.0f && model.magnitude[4] == 16.0f);

        /* Next to a frame with no voiced band the mid-frame pitch would be
         * this frame's alone: no band there either, in both orders, and the
         * frame renders nothing outside its 160 samples. */
        for (int order = 0; order < 2; ++order) {
            const mbe_parms* before = order ? &unvoiced : &voiced;
            const mbe_parms* after = order ? &voiced : &unvoiced;
            for (int flushed = 0; flushed < 2; ++flushed) {
                const unsigned saved = flushed ? flush_denormals() : 0u;
                mbe_subframe_mid_model(before, after, 0, &model);
                assert(!model.any_voiced && model.w0 == 0.0f);
                for (int l = 1; l <= model.L; ++l) {
                    assert(model.voiced[l] == 2 && model.magnitude[l] == 0.0f);
                }
                fresh_stream(&start);
                next_frame(guarded + 16, &start, before);
                for (int n = 0; n < 192; ++n) {
                    guarded[n] = (n < 16 || n >= 176) ? 7.0f : 0.0f;
                }
                next_frame(guarded + 16, &start, after);
                for (int n = 0; n < 16; ++n) {
                    assert(guarded[n] == 7.0f && guarded[176 + n] == 7.0f);
                }
                assert(all_finite(guarded + 16, 160));
                if (flushed) {
                    restore_denormals(saved);
                }
            }
        }

        /* Through the P25 data path, as a caller's state can hold it: the last
         * synthesized frame voiced at such a pitch, then a frame with no voiced
         * band. */
        for (int flushed = 0; flushed < 2; ++flushed) {
            mbe_parms cur, prev, enhanced;
            char imbe_d[88] = {0};
            float out[160];
            const unsigned saved = flushed ? flush_denormals() : 0u;
            mbe_initMbeParms(&cur, &prev, &enhanced);
            enhanced.L = 18;
            memcpy(&enhanced.w0, &w0_bits[i], sizeof(enhanced.w0));
            for (int l = 1; l <= 56; ++l) {
                enhanced.Vl[l] = 1;
                enhanced.Ml[l] = 16.0f;
            }
            set_imbe_b0(imbe_d, 150);
            const int rc = mbe_processImbe4400Dataf(out, NULL, imbe_d, &cur, &prev, &enhanced);
            assert(rc >= 0);
            (void)rc;
            assert(all_finite(out, 160));
            if (flushed) {
                restore_denormals(saved);
            }
        }
    }
}

/* The parameter bit that voices a P25 frame's lowest band (harmonics 1..3) and
 * nothing else, for pitch code b0; -1 if there is none. */
static int
imbe_lowest_band_bit(int b0) {
    for (int bit = 6; bit < 85; ++bit) {
        mbe_parms cur, prev, enhanced;
        char imbe_d[88] = {0};
        set_imbe_b0(imbe_d, b0);
        imbe_d[bit] = 1;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        if (mbe_decodeImbe4400Parms(imbe_d, &cur, &prev) != 0) {
            continue;
        }
        int others = 0;
        for (int l = 4; l <= cur.L; ++l) {
            others |= cur.Vl[l];
        }
        if (cur.Vl[1] && cur.Vl[2] && cur.Vl[3] && !others) {
            return bit;
        }
    }
    return -1;
}

/* A P25 frame whose voicing does not reach MBE_VOICING_REACH_HZ is marked as a
 * late onset: three voiced harmonics at pitch codes 56 (502.6 Hz) and 100
 * (344 Hz) are, at codes 55 (507.9 Hz) and 30 (691 Hz) they are not. */
static void
test_p25_marks_late_onsets(void) {
    static const int codes[] = {56, 100, 207, 55, 30, 0};
    for (size_t i = 0; i < sizeof(codes) / sizeof(codes[0]); ++i) {
        mbe_parms cur, prev, enhanced;
        char imbe_d[88] = {0};
        const int bit = imbe_lowest_band_bit(codes[i]);
        assert(bit >= 0);
        set_imbe_b0(imbe_d, codes[i]);
        imbe_d[bit] = 1;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        const int kind = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
        assert(kind == 0);
        (void)kind;
        assert(mbe_frame_mark(&cur) == ((i < 3) ? MBE_FRAME_MARK_LATE_ONSET : MBE_FRAME_MARK_NONE));

        /* Through the data path, after frames with no voiced band: the pulse
         * train starts at mid-frame when marked, at the frame's start if not. */
        char unvoiced[88] = {0};
        float out[160];
        set_imbe_b0(unvoiced, 150);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int f = 0; f < 2; ++f) {
            const int rc = mbe_processImbe4400Dataf(out, NULL, unvoiced, &cur, &prev, &enhanced);
            assert(rc >= 0);
            (void)rc;
        }
        const int rc2 = mbe_processImbe4400Dataf(out, NULL, imbe_d, &cur, &prev, &enhanced);
        assert(rc2 >= 0);
        (void)rc2;
        const double hz = (double)enhanced.w0 * 8000.0 / (2.0 * M_PI);
        assert(fabs(cycles_apart(pulse_phase(&enhanced), cycles_of(hz, (i < 3) ? 80.0 : 160.0))) < 1e-4);
    }
}

/* The P25 and D-STAR decoders mark a frame with no voiced band by its pitch
 * code, replacing whatever mark the caller's state held. */
static void
test_decoders_mark_unvoiced_frames(void) {
    mbe_parms cur, prev, enhanced;
    char imbe_d[88];

    for (int b0 = 0; b0 <= 207; ++b0) {
        for (int voiced = 0; voiced < 2; ++voiced) {
            memset(imbe_d, voiced, sizeof(imbe_d));
            set_imbe_b0(imbe_d, b0);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            mbe_set_frame_mark(&cur, MBE_FRAME_MARK_LATE_ONSET);
            const int kind = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
            assert(kind == 0);
            (void)kind;
            int any = 0;
            for (int l = 1; l <= cur.L; ++l) {
                any |= cur.Vl[l];
            }
            assert(any == voiced);
            /* b0 128..183 are the codes DVSI's encoder keeps for frames with no voiced band. */
            const int own_code = b0 >= 128 && b0 <= 183;
            assert(mbe_frame_mark(&cur)
                   == ((!voiced && !own_code) ? MBE_FRAME_MARK_PITCH_ABSENT : MBE_FRAME_MARK_NONE));
        }
    }

    for (int b0 = 0; b0 <= 124; ++b0) {
        for (int b1 = 0; b1 < 16; b1 += 15) {
            char ambe_d[49];
            set_ambe2400(ambe_d, b0, b1);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            mbe_set_frame_mark(&cur, MBE_FRAME_MARK_LATE_ONSET);
            const int kind2 = mbe_decodeAmbe2400Parms(ambe_d, &cur, &prev);
            assert(kind2 == 0);
            (void)kind2;
            int any = 0;
            for (int l = 1; l <= cur.L; ++l) {
                any |= cur.Vl[l];
            }
            /* b0 64..119 are D-STAR's codes for frames with no voiced band;
             * three ranges of them carry their pitch. */
            enum mbe_frame_mark want = MBE_FRAME_MARK_NONE;
            if (!any && (b0 < 64 || b0 > 119)) {
                want = MBE_FRAME_MARK_PITCH_ABSENT;
            } else if (!any && ((b0 >= 64 && b0 <= 71) || (b0 >= 88 && b0 <= 103) || (b0 >= 112 && b0 <= 119))) {
                want = MBE_FRAME_MARK_PITCH_CARRIES;
            }
            assert(mbe_frame_mark(&cur) == want);
        }
    }
}

/* Through the P25 data path: a voiced frame, then an unvoiced frame at an
 * ordinary pitch code (b0 25 is the one DVSI's encoder sends), keeps the pulse
 * train running at the voiced frame's pitch; an unvoiced frame at one of the
 * codes kept for such frames ends it. */
static void
test_p25_voicing_persists(void) {
    static const int codes[] = {25, 26, 100, 127, 184, 128, 150, 183};
    for (size_t i = 0; i < sizeof(codes) / sizeof(codes[0]); ++i) {
        const int persists = i < 5;
        mbe_parms cur, prev, enhanced;
        char voiced[88], unvoiced[88];
        float out[160];
        mbe_initMbeParms(&cur, &prev, &enhanced);
        memset(voiced, 1, sizeof(voiced));
        set_imbe_b0(voiced, 100); /* 114.7 Hz, voiced throughout */
        memset(unvoiced, 0, sizeof(unvoiced));
        set_imbe_b0(unvoiced, codes[i]);
        for (int f = 0; f < 3; ++f) {
            const int rc = mbe_processImbe4400Dataf(out, NULL, voiced, &cur, &prev, &enhanced);
            assert(rc >= 0);
            (void)rc;
        }
        const double before = pulse_phase(&enhanced);
        const double hz = (double)enhanced.w0 * 8000.0 / (2.0 * M_PI);
        const int rc2 = mbe_processImbe4400Dataf(out, NULL, unvoiced, &cur, &prev, &enhanced);
        assert(rc2 >= 0);
        (void)rc2;
        assert(mbe_frame_mark(&enhanced) == (persists ? MBE_FRAME_MARK_PITCH_ABSENT : MBE_FRAME_MARK_NONE));
        if (persists) {
            assert(fabs(cycles_apart(pulse_phase(&enhanced), before + cycles_of(hz, 160.0))) < 1e-4);
            assert(pulse_phase(&enhanced) != 0.0);
        } else {
            assert(pulse_phase(&enhanced) == 0.0);
        }
        assert(all_finite(out, 160));
    }
}

/* P25's reset after a run of repeats replaces the frame's model, and the mark
 * the rejected frame's bits gave it goes with the model: two streams that
 * differ only in the bits of their rejected frames play alike, through the
 * reset and after it. */
static void
test_p25_reset_drops_the_mark(void) {
    static const int rejected_codes[2] = {25, 150}; /* marked as lacking a pitch, and not marked */
    mbe_parms cur[2], prev[2], enhanced[2];
    float out[2][160];
    char voiced[88];
    int resets = 0;

    memset(voiced, 1, sizeof(voiced));
    set_imbe_b0(voiced, 100);
    for (int s = 0; s < 2; ++s) {
        mbe_initMbeParms(&cur[s], &prev[s], &enhanced[s]);
    }
    for (int f = 0; f < 16; ++f) {
        const int rejected = f >= 3 && f < 12;
        for (int s = 0; s < 2; ++s) {
            mbe_process_result result;
            char imbe_d[88] = {0};
            mbe_initProcessResult(&result);
            if (rejected) {
                set_imbe_b0(imbe_d, rejected_codes[s]);
                result.total_errors = 6; /* past the repeat criterion of the parameter-bit API */
            } else {
                memcpy(imbe_d, voiced, sizeof(imbe_d));
            }
            /* A frame at the repeat limit plays mute noise from the thread's generator: the same for both. */
            mbe_setThreadRngSeed(1000u + (uint32_t)f);
            const int rc = mbe_processImbe4400Dataf(out[s], &result, imbe_d, &cur[s], &prev[s], &enhanced[s]);
            assert(rc >= 0);
            (void)rc;
            assert(((result.flags & MBE_PROCESS_FLAG_REPEAT) != 0u) == rejected);
            /* A reset frame is a repeat that is neither muted nor a copy of the frame before: count them. */
            if (s == 0 && rejected && enhanced[s].repeatCount == 0) {
                ++resets;
            }
        }
        assert(same_bits(out[0], out[1], 160));
        assert(mbe_frame_mark(&enhanced[0]) == mbe_frame_mark(&enhanced[1]));
        assert(same_bits(enhanced[0].noiseOverlap, enhanced[1].noiseOverlap, 96));
        assert(same_bits(enhanced[0].PSIl, enhanced[1].PSIl, 57));
        assert(same_bits(&enhanced[0].noiseSeed, &enhanced[1].noiseSeed, 1));
    }
    assert(resets > 0);
}

int
main(void) {
    test_frame_marks();
    test_steady_model_is_its_harmonics();
    test_onset_starts_the_train();
    test_late_onset();
    test_pitch_step_keeps_the_train();
    test_waveform_through_changes();
    test_unvoiced_pitch_carries();
    test_voicing_persists_without_a_pitch();
    test_mid_model();
    test_previous_model_held();
    test_subframes_crossfade_linearly();
    test_muted_frame_advances();
    test_from_silence();
    test_output_limit_option();
    test_pulse_state_is_validated();
    test_invalid_models();
    test_subnormal_pitch_is_no_pitch();
    test_decoders_mark_unvoiced_frames();
    test_p25_marks_late_onsets();
    test_p25_voicing_persists();
    test_p25_reset_drops_the_mark();
    return 0;
}
