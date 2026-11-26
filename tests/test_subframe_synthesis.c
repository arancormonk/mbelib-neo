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
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "ambe3600x2450_internal.h"
#include "imbe4400_internal.h"
#include "mbe_adaptive.h"
#include "mbe_subframe.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Where noiseOverlap holds the pulse phase (mbe_subframe.h). */
#define PHASE_SLOT (1 + MBE_PULSE_HISTORY)

/* Whether the library is built with fast-math (CMake defines
 * MBELIB_TEST_LIBRARY_FAST_MATH): it may then evaluate an expression otherwise
 * than this test does. */
#if defined(MBELIB_TEST_LIBRARY_FAST_MATH) || defined(__FAST_MATH__) || defined(_M_FP_FAST)
#define MID_PITCH_FAST_MATH 1
#else
#define MID_PITCH_FAST_MATH 0
#endif

/* Whether the library's evaluation of the mid-frame pitch, an expression in
 * double, may differ from this test's evaluation of the same expression by
 * the last bit: under fast-math, under excess precision (x87, or an
 * evaluation method left open), and where the target fuses multiply-adds,
 * which the compiler may contract in one translation unit and not in the
 * other. Elsewhere IEEE arithmetic gives the two the same bits. */
#if MID_PITCH_FAST_MATH || !defined(FLT_EVAL_METHOD) || (FLT_EVAL_METHOD != 0) || defined(__FP_FAST_FMA)
#define MID_PITCH_INEXACT 1
#else
#define MID_PITCH_INEXACT 0
#endif

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

/* Harmonic l's amplitude, decoded as log2Ml and rendered as Ml alike (a model
 * no enhancement has shaped; 0 for none, whose log2Ml is the initial state's). */
static void
set_amplitude(mbe_parms* mp, int l, float amplitude) {
    mp->Ml[l] = amplitude;
    mp->log2Ml[l] = (amplitude > 0.0f) ? (float)log2((double)amplitude) : 0.0f;
}

/* A model of L harmonics at `hz`, the first `voiced` of them voiced. */
static void
set_model(mbe_parms* mp, double hz, int L, int voiced, float amplitude) {
    mp->w0 = w0_of(hz);
    mp->L = L;
    for (int l = 1; l <= 56; ++l) {
        mp->Vl[l] = (l <= voiced) ? 1 : 0;
        set_amplitude(mp, l, (l <= L) ? amplitude : 0.0f);
    }
    mbe_set_frame_mark(mp, MBE_FRAME_MARK_NONE);
}

/* A model's amplitudes, as rendered (Ml) and as decoded (log2Ml, which the
 * phases are regenerated from). */
static void
copy_amplitudes(mbe_parms* dst, const mbe_parms* src) {
    memcpy(dst->Ml, src->Ml, sizeof(dst->Ml));
    memcpy(dst->log2Ml, src->log2Ml, sizeof(dst->log2Ml));
}

/* A stream's next frame: its state, with `model`'s pitch, voicing and amplitudes. */
static mbe_parms
with_model(const mbe_parms* state, const mbe_parms* model) {
    mbe_parms cur = *state;
    cur.w0 = model->w0;
    cur.L = model->L;
    memcpy(cur.Vl, model->Vl, sizeof(cur.Vl));
    copy_amplitudes(&cur, model);
    mbe_set_frame_mark(&cur, MBE_FRAME_MARK_NONE);
    return cur;
}

/* The next frame of a stream: `last` is the frame before it as synthesized,
 * `model` the new frame's model and mark. Renders the voiced part alone (no
 * transform plan), at unit gain, the way `rendering` says, and leaves the new
 * state in `last`. */
static void
next_frame_as(float out[160], mbe_parms* last, const mbe_parms* model, int from_silence,
              enum mbe_voiced_synthesis rendering) {
    mbe_parms cur = *last;
    cur.w0 = model->w0;
    cur.L = model->L;
    memcpy(cur.Vl, model->Vl, sizeof(cur.Vl));
    copy_amplitudes(&cur, model);
    mbe_set_frame_mark(&cur, mbe_frame_mark(model));
    memset(out, 0, 160 * sizeof(float));
    mbe_synthesize_subframes(out, &cur, last, NULL, 1.0f, from_silence, rendering);
    *last = cur;
}

static void
next_frame_from(float out[160], mbe_parms* last, const mbe_parms* model, int from_silence) {
    next_frame_as(out, last, model, from_silence, MBE_VOICED_SYNTHESIS_DVSI);
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

/* Whether the pulse phase a frame stored lies in [0, 1) as stored, read from
 * its bits (a comparison modulo a cycle takes 1.0 for 0). */
static int
stored_phase_in_range(const mbe_parms* mp) {
    uint32_t bits;
    memcpy(&bits, &mp->noiseOverlap[PHASE_SLOT], sizeof(bits));
    return bits < UINT32_C(0x3F800000); /* +0 up to, not including, 1.0f */
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

/* P25 parameter bits with pitch code b0. */
static void
set_imbe_b0(char imbe_d[88], int b0) {
    static const int index[8] = {0, 1, 2, 3, 4, 5, 85, 86};
    for (int i = 0; i < 8; ++i) {
        imbe_d[index[i]] = (char)((b0 >> (7 - i)) & 1);
    }
}

/* The frame marks: only an exact stored value counts, and one replaces another. */
static void
test_frame_marks(void) {
    static const enum mbe_frame_mark marks[] = {MBE_FRAME_MARK_NONE,         MBE_FRAME_MARK_PITCH_CARRIES,
                                                MBE_FRAME_MARK_PITCH_ABSENT, MBE_FRAME_MARK_PREVIOUS_HELD,
                                                MBE_FRAME_MARK_TONE,         MBE_FRAME_MARK_ONSET_COLUMNS,
                                                MBE_FRAME_MARK_MID_COLUMN_0, MBE_FRAME_MARK_MID_COLUMNS_0_1,
                                                MBE_FRAME_MARK_LATE_ONSET};

    static const uint32_t others[] = {0x7FC00000u, 0x7F800000u, 0xFF800000u, 0x80000000u, 0x3F800001u, 0x40400001u,
                                      0x40A00001u, 0x40C00001u, 0x40E00001u, 0x40BFFFFFu, 0x40900000u, 0x41000000u,
                                      0x429E0000u, 0xBF800000u, 0xC0000000u, 0x00000001u};
    mbe_parms mp;
    memset(&mp, 0, sizeof(mp));
    assert(mbe_frame_mark(&mp) == MBE_FRAME_MARK_NONE);
    for (size_t i = 0; i < sizeof(marks) / sizeof(marks[0]); ++i) {
        mbe_set_frame_mark(&mp, marks[i]);
        assert(mbe_frame_mark(&mp) == marks[i]);
        /* Stored as the mark's value in float. */
        const float stored = (float)marks[i];
        assert(same_bits(&mp.PHIl[0], &stored, 1));
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
        set_amplitude(&model, l, 8.0f + (float)((l * 7) % 5));
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
    mbe_synthesize_subframes(unit, &a, &a_prev, NULL, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
    mbe_synthesize_subframes(scaled, &b, &b_prev, NULL, 1.25f, 0, MBE_VOICED_SYNTHESIS_DVSI);
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

/* A frame marked as a late onset (rate 33's V/UV row 15) starts its voiced
 * bands in its second subframe, in the rate-33 rendering whatever the frame
 * before: after a frame with no voiced band no pulse comes before sample 80;
 * after a voiced frame the mid-frame model voices nothing too, and the pulse
 * train runs on, not restarted, through the first subframe at the previous
 * frame's pitch and through the second at this frame's. In the other
 * renderings the mark acts after a frame with no voiced band alone, as before
 * the rate-33 rows' reading: after a voiced frame it changes nothing. */
static void
test_late_onset(void) {
    const enum mbe_voiced_synthesis rate33 = MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    float out[160], plain[160];
    mbe_parms last, voiced;
    fresh_stream(&last);
    set_model(&last, 150.0, 20, 0, 4.0f);
    set_model(&voiced, 230.0, 15, 15, 16.0f);
    mbe_set_frame_mark(&voiced, MBE_FRAME_MARK_LATE_ONSET);
    mbe_parms other = last;
    float other_out[160];
    next_frame_as(out, &last, &voiced, 0, rate33);
    next_frame(other_out, &other, &voiced);
    assert(same_bits(out, other_out, 160) && same_bits(last.noiseOverlap, other.noiseOverlap, 96));
    assert(fabs(cycles_apart(pulse_phase(&last), cycles_of(230.0, 80.0))) < 1e-5);
    assert(last.noiseOverlap[0] == 3.0f);
    assert(last.noiseOverlap[3] == -80.0f);
    for (int n = 0; n < 80; ++n) {
        assert(out[n] == 0.0f);
    }
    assert(out[100] != 0.0f);

    /* After a voiced frame at 200 Hz. */
    mbe_parms before;
    set_model(&before, 200.0, 18, 18, 16.0f);
    next_frame_as(out, &last, &before, 0, rate33);
    const double phase = pulse_phase(&last);
    mbe_parms marked = last, unmarked = last, other_marked = last, other_unmarked = last;
    next_frame_as(out, &marked, &voiced, 0, rate33);
    assert(fabs(cycles_apart(pulse_phase(&marked), phase + cycles_of(200.0, 80.0) + cycles_of(230.0, 80.0))) < 1e-5);
    struct mbe_subframe_model mid;
    mbe_subframe_mid_model(&last, &voiced, 1, rate33, &mid);
    assert(!mid.any_voiced);
    next_frame(other_out, &other_marked, &voiced);
    mbe_set_frame_mark(&voiced, MBE_FRAME_MARK_NONE);
    next_frame_as(plain, &unmarked, &voiced, 0, rate33);
    assert(!same_bits(out, plain, 160));
    next_frame(plain, &other_unmarked, &voiced);
    assert(same_bits(other_out, plain, 160) && same_bits(other_marked.noiseOverlap, other_unmarked.noiseOverlap, 96));
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
 * The models' pitch, voicing and magnitudes are the library's, which
 * test_mid_model checks on their own; their phases are regenerated here from
 * the frames' log2Ml (reference_phases()).
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

/*
 * The phases US 5,701,390's edge-detection kernel regenerates from log2
 * amplitudes envelope[1..L], in double: 0 at l = 0, each value below log2(1e-6)
 * taken as that, falling by 0.72 per harmonic past L, mirrored below l = 1;
 * phase l is 0.44 times the sum over odd m up to 19 of (B[l + m] - B[l - m]) / m.
 */
static void
reference_phases(const double envelope[57], int L, double phase[57]) {
    double B[1 + 56 + 2 * 19]; /* B[19 + l], l from -19 to 75 */
    B[19] = 0.0;
    for (int l = 1; l <= L; ++l) {
        B[19 + l] = fmax(envelope[l], log2(1e-6));
    }
    for (int l = L + 1; l <= 75; ++l) {
        B[19 + l] = B[19 + L] - (0.72 * (double)(l - L));
    }
    for (int l = 1; l <= 19; ++l) {
        B[19 - l] = B[19 + l];
    }
    phase[0] = 0.0;
    for (int l = 1; l <= 56; ++l) {
        double sum = 0.0;
        for (int m = 1; m <= 19; m += 2) {
            sum += (B[19 + l + m] - B[19 + l - m]) / (double)m;
        }
        phase[l] = 0.44 * sum;
    }
}

/* A frame's log2Ml at angular frequency w: linear between the two harmonics
 * around it, the end harmonic's beyond either end. */
static double
reference_log2_at(const mbe_parms* mp, double w) {
    const double position = w / (double)mp->w0;
    if (position < 1.0) {
        return (double)mp->log2Ml[1];
    }
    if (position >= (double)mp->L) {
        return (double)mp->log2Ml[mp->L];
    }
    const int below = (int)position;
    return (double)mp->log2Ml[below]
           + ((position - below) * ((double)mp->log2Ml[below + 1] - (double)mp->log2Ml[below]));
}

/* The phases of a frame's own model, from its log2Ml as stored. */
static void
reference_own_phases(const mbe_parms* mp, double phase[57]) {
    double envelope[57] = {0.0};
    for (int l = 1; l <= mp->L; ++l) {
        envelope[l] = (double)mp->log2Ml[l];
    }
    reference_phases(envelope, mp->L, phase);
}

/* The phases of the mid-frame model of pitch w0 and L harmonics between two
 * frames: from the mean of the two frames' log2Ml at each harmonic's frequency. */
static void
reference_mid_phases(const mbe_parms* prev, const mbe_parms* cur, double w0, int L, double phase[57]) {
    double envelope[57] = {0.0};
    for (int l = 1; l <= L; ++l) {
        envelope[l] = 0.5 * (reference_log2_at(prev, w0 * l) + reference_log2_at(cur, w0 * l));
    }
    reference_phases(envelope, L, phase);
}

/* What a model plays t samples after a pulse, its harmonics at `phase`. */
static double
reference_response(const struct mbe_subframe_model* model, const double phase[57], double t) {
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
               * cos((w * (t - 16.0)) + phase[l] + (double)mbe_low_band_phase((float)hz));
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
    } else if (from / to > 1.6 && from / to < 2.4) {
        /* The new period about twice the old: the old train an octave down. */
        train->phase /= 2.0;
        from /= 2.0;
    } else if (from / to > 0.4 && from / to < 0.6) {
        /* About half: an octave up, its phase doubled within a cycle. */
        train->phase = fmod(2.0 * train->phase, 1.0);
        from *= 2.0;
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

/* One subframe from sample `start`: every pulse so far plays model a going out
 * and model b coming in, at their phases. */
static void
reference_render(double* out, double start, const struct reference_train* train, const struct mbe_subframe_model* a,
                 const double a_phase[57], const struct mbe_subframe_model* b, const double b_phase[57]) {
    for (int n = 0; n < 80; ++n) {
        const double fade = (double)n / 80.0;
        double sum = 0.0;
        for (int i = 0; i < train->count; ++i) {
            const double t = start + (double)n - train->time[i];
            sum += ((1.0 - fade) * reference_response(a, a_phase, t)) + (fade * reference_response(b, b_phase, t));
        }
        out[n] = sum;
    }
}

/* Whether a frame has a voiced band. */
static int
reference_voiced(const mbe_parms* mp) {
    for (int l = 1; l <= mp->L; ++l) {
        if (mp->Vl[l] == 1) {
            return 1;
        }
    }
    return 0;
}

/* The half-octave wave of rate 33's mid-frame pitch at pitch w0, in double:
 * h(P) = 4.244e-5 + 9.374e-5 cos t + 8.388e-4 sin t, t = 4 pi log2 P, with P
 * = 2 pi / w0 the period in samples (docs/dvsi-parity.md, "Rate 33: mid-frame
 * pitch"). */
static double
reference_r33_wave(double w0) {
    const double t = 4.0 * M_PI * log2(2.0 * M_PI / w0);
    return 4.244e-5 + (9.374e-5 * cos(t)) + (8.388e-4 * sin(t));
}

/* The factor 1 - 3.0e-4 as the library stores it, a float (rounded to float
 * also under x87 excess precision), promoted to double: the double literal
 * 1.0 - 3.0e-4 rounds otherwise for some pairs of pitches. */
static double
reference_r33_factor(void) {
    volatile float stored = 1.0f - 3.0e-4f;
    return (double)stored;
}

/* Rate 33's mid-frame pitch between two voiced frames of different pitch, in
 * double, unrounded: the geometric mean, times the stored factor, times 1
 * plus the mean of the two frames' waves, in the library's order. */
static double
reference_r33_mid_pitch(double a, double b) {
    const double wave = 0.5 * (reference_r33_wave(a) + reference_r33_wave(b));
    return sqrt(a * b) * reference_r33_factor() * (1.0 + wave);
}

/* The mid-frame model's harmonic count at its pitch w0, between two frames
 * of L harmonics `shorter` and `longer`: every harmonic up to 0.9254 of
 * Nyquist, in float as the library computes it, held between the two. */
static int
reference_mid_count(float w0, int shorter, int longer) {
    const float harmonics = floorf((0.9254f * (float)M_PI) / w0);
    return (harmonics < (float)shorter) ? shorter : ((harmonics > (float)longer) ? longer : (int)harmonics);
}

/* How far (0.9254 pi) / w0 lies from a whole number, relative to it: where a
 * build may round the pitch an ulp otherwise than this test, the harmonic
 * count it gives is decided only away from a whole number. */
static double
count_margin(float w0) {
    const double q = (0.9254 * M_PI) / (double)w0;
    return fabs(q - floor(q + 0.5)) / q;
}

/* Whether a lies within an ulp of float of b. */
static int
within_an_ulp(float a, float b) {
    const double ulp = (double)nextafterf(fabsf(b), INFINITY) - (double)fabsf(b);
    return fabs((double)a - (double)b) <= ulp;
}

/* Whether a lies within an ulp of float of the double b, b within the float
 * range: the ulp at b's nearest float, the largest finite one's at the top
 * of the range (no infinite tolerance there). Outside the range, no float
 * is within an ulp of b. */
static int
within_an_ulp_of(float a, double b) {
    if (!(fabs(b) <= (double)FLT_MAX)) {
        return 0;
    }
    const float nearest = (float)b;
    const float magnitude = fabsf(nearest);
    const double ulp = (magnitude < FLT_MAX) ? (double)nextafterf(magnitude, INFINITY) - (double)magnitude
                                             : (double)FLT_MAX - (double)nextafterf(FLT_MAX, 0.0f);
    return fabs((double)a - b) <= ulp;
}

/* The mid-frame pitch between two frames: the previous frame's own where a
 * silence frame holds it; the voiced frame's where only one is voiced, unless
 * the unvoiced one before is marked as carrying its pitch; else the geometric
 * mean (in float, sqrtf of the product, as every rendering takes it), which
 * the rate-33 rendering instead forms in double, times the stored factor and
 * the half-octave wave, between two voiced frames of different pitch
 * (reference_r33_mid_pitch(), here rounded to float once). */
static float
reference_mid_pitch(const mbe_parms* prev, const mbe_parms* cur, enum mbe_voiced_synthesis rendering) {
    const int prev_voiced = reference_voiced(prev), cur_voiced = reference_voiced(cur);
    if (mbe_frame_mark(cur) == MBE_FRAME_MARK_PREVIOUS_HELD || (prev_voiced && !cur_voiced)) {
        return prev->w0;
    }
    if (cur_voiced && !prev_voiced && mbe_frame_mark(prev) != MBE_FRAME_MARK_PITCH_CARRIES) {
        return cur->w0;
    }
    const int differ = !same_bits(&prev->w0, &cur->w0, 1);
    return (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 && prev_voiced && cur_voiced && differ)
               ? (float)reference_r33_mid_pitch((double)prev->w0, (double)cur->w0)
               : sqrtf(prev->w0 * cur->w0);
}

/* The P25 onset rule, in double, from its statement: a voiced harmonic at `hz`
 * of a frame after one with no voiced band sounds at mid-frame when the 500 Hz
 * column holding hz, or the one holding hz + 250 Hz, is covered; columns
 * [0, 505), [505, 1005), ... [3005, 3505) Hz, none above, each covered when
 * the frame is voiced at harmonic 1 (for 0 Hz) or at the harmonic at or next
 * above the edge (unvoiced past L) at both its ends. */
static int
reference_onset_sounds(const mbe_parms* mp, double hz) {
    const double f0 = (double)mp->w0 * 8000.0 / (2.0 * M_PI);
    int covered[8] = {0}, edge_voiced[8];
    edge_voiced[0] = mp->Vl[1] == 1;
    for (int k = 1; k <= 7; ++k) {
        const int n = (int)ceil((505.0 + 500.0 * (k - 1)) / f0);
        edge_voiced[k] = (n <= mp->L) && (mp->Vl[n] == 1);
    }
    for (int k = 0; k < 7; ++k) {
        covered[k] = edge_voiced[k] && edge_voiced[k + 1];
    }
    int column = 0, reach = 0;
    for (int k = 0; k < 7; ++k) {
        column += (hz >= 505.0 + 500.0 * k);
        reach += (hz + 250.0 >= 505.0 + 500.0 * k);
    }
    return covered[column] || covered[reach];
}

/* The mid-frame model's band states (0 unvoiced, 1 voiced, 2 none) at pitch
 * w0 and L harmonics, from the rules as docs/dvsi-parity.md states them: each
 * harmonic takes the voicing of the nearest harmonic of the frame whose
 * voicing holds at mid-frame (the previous frame's when this one has no voiced
 * band and is marked as lacking a pitch, and the previous one is voiced at its
 * harmonic at or next above 505 Hz; else this frame's); a frame marked as a
 * late onset voices none (reference_late()); in P25's rendering a frame
 * marked for the onset columns after one with no voiced band voices only the
 * harmonics the P25 onset rule sounds; and in the rate-33 rendering a frame of
 * V/UV row 17 or 13 voices the harmonics under Nyquist in the 500 Hz columns
 * below 1 or 2 (by the frequency of each, or of 250 Hz above it), whatever the
 * voicing around. Each mark acts in its own codec's rendering alone. */
static void
reference_mid_voicing(const mbe_parms* prev, const mbe_parms* cur, double w0, int L,
                      enum mbe_voiced_synthesis rendering, unsigned char voiced[57]) {
    const enum mbe_frame_mark mark = mbe_frame_mark(cur);
    const int onset = !reference_voiced(prev);
    const int rate33 = rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    const int late = mark == MBE_FRAME_MARK_LATE_ONSET && (rate33 || onset);
    const mbe_parms* source = cur;
    if (mark == MBE_FRAME_MARK_PITCH_ABSENT && !reference_voiced(cur)) {
        const int n = (int)ceil(505.0 * 2.0 * M_PI / 8000.0 / (double)prev->w0);
        source = (n <= prev->L && prev->Vl[n] == 1) ? prev : cur;
    }
    for (int l = 1; l <= L; ++l) {
        const double w = w0 * l;
        const int nearest = (int)fmin(fmax(floor(w / (double)source->w0 + 0.5), 1.0), (double)source->L);
        int state = (source->Vl[nearest] == 0) ? 0 : ((source->Vl[nearest] == 1) ? 1 : 2);
        if (state == 1 && late) {
            state = 0;
        }
        if (state == 1 && onset && rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE && mark == MBE_FRAME_MARK_ONSET_COLUMNS
            && !reference_onset_sounds(cur, w * 8000.0 / (2.0 * M_PI))) {
            state = 0;
        }
        int columns = (mark == MBE_FRAME_MARK_MID_COLUMN_0) ? 1 : ((mark == MBE_FRAME_MARK_MID_COLUMNS_0_1) ? 2 : 0);
        columns = rate33 ? columns : 0;
        if (columns > 0 && state != 2) {
            const double hz = w * 8000.0 / (2.0 * M_PI);
            state = (w < M_PI)
                    && (fmin(floor(hz / 500.0), 7.0) < columns || fmin(floor((hz + 250.0) / 500.0), 7.0) < columns);
        }
        voiced[l] = (unsigned char)state;
    }
}

/* The mark a frame is rendered with: in P25's rendering the onset-columns
 * mark for any model with a voiced band that is not a tone model, none for one
 * with no voiced band that carried it, as the rendering marks the model it
 * renders; in the others the frame's own. */
static enum mbe_frame_mark
reference_rendered_mark(const mbe_parms* mp, enum mbe_frame_mark mark, enum mbe_voiced_synthesis rendering) {
    if (rendering != MBE_VOICED_SYNTHESIS_DVSI_IMBE || mark == MBE_FRAME_MARK_TONE) {
        return mark;
    }
    if (reference_voiced(mp)) {
        return MBE_FRAME_MARK_ONSET_COLUMNS;
    }
    return (mark == MBE_FRAME_MARK_ONSET_COLUMNS) ? MBE_FRAME_MARK_NONE : mark;
}

/* Give frame f's model a decoded spectral shape (log2Ml, which its phases come
 * from) and render it with a high-band shelf on top (Ml), as enhancement shapes
 * a decoded model. */
static void
shape_model(mbe_parms* mp, size_t f) {
    for (int l = 1; l <= mp->L; ++l) {
        const double decoded = (double)mp->Ml[l] * exp2(0.8 * sin((0.7 * l) + (double)f));
        set_amplitude(mp, l, (float)decoded);
        mp->Ml[l] = (float)(decoded * ((2 * l > mp->L) ? 1.6 : 1.0));
    }
}

/* The frames a stream plays through pitch steps, changes of harmonic count,
 * level and voicing, onsets and the frame marks, against that rendering:
 * every sample, so that pulses, models, windows and crossfades have to agree
 * across the subframe and frame boundaries, not the stored state alone. The
 * mid-frame pitch the pulses follow is the reference's own
 * (reference_mid_pitch()), for the codec's rendering. */
static void
test_waveform_through_changes(enum mbe_voiced_synthesis rendering, int inexact) {
    static const struct {
        double hz;
        int L, voiced;
        float amplitude;
        enum mbe_frame_mark mark;
        int from_silence;
    } frames[] = {
        {200.0, 18, 14, 16.0f, MBE_FRAME_MARK_NONE, 0},          /* an onset */
        {200.0, 18, 14, 16.0f, MBE_FRAME_MARK_NONE, 0},          /* steady */
        {230.0, 16, 16, 16.0f, MBE_FRAME_MARK_NONE, 0},          /* pitch up, fewer harmonics */
        {180.0, 20, 12, 16.0f, MBE_FRAME_MARK_NONE, 0},          /* pitch down, more of them, fewer voiced */
        {207.0, 17, 17, 40.0f, MBE_FRAME_MARK_NONE, 0},          /* a level step */
        {100.0, 30, 0, 8.0f, MBE_FRAME_MARK_PITCH_CARRIES, 0},   /* no voiced band, its pitch carried */
        {180.0, 20, 20, 16.0f, MBE_FRAME_MARK_NONE, 0},          /* onset at the mean of the two pitches */
        {250.0, 14, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, 0},    /* the voicing before lasts to mid-frame */
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 0},          /* and the train goes on */
        {150.0, 24, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0},     /* after a voiced frame too: the train runs on */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0},             /* voicing ends with the frame's start */
        {127.5, 28, 23, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0}, /* a P25 onset: harmonics 20-23 from mid-frame */
        {127.5, 28, 23, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0}, /* inside the run the mark changes nothing */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0},
        {160.0, 25, 3, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0}, /* harmonics 1-3 to 480 Hz: all from mid-frame */
        {100.0, 30, 0, 8.0f, MBE_FRAME_MARK_PITCH_CARRIES, 0},
        {205.0, 17, 17, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0}, /* at the mean pitch: those under 3505 Hz */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0},
        {170.0, 22, 8, 16.0f, MBE_FRAME_MARK_NONE, 0},
        {165.0, 23, 0, 16.0f, MBE_FRAME_MARK_MID_COLUMN_0, 0},    /* rate-33 row 17: column 0 at mid-frame, */
        {150.0, 24, 9, 16.0f, MBE_FRAME_MARK_MID_COLUMNS_0_1, 0}, /* row 13 after it: columns 0-1, the train anew */
        {150.0, 24, 9, 16.0f, MBE_FRAME_MARK_MID_COLUMNS_0_1, 0}, /* and inside a voiced run */
        {140.0, 26, 10, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0},     /* row 15 after a voiced frame */
        {130.0, 28, 0, 8.0f, MBE_FRAME_MARK_MID_COLUMN_0, 0},
        {125.0, 29, 10, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0}, /* row 15 after row 17: a new train at mid-frame */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0},
        {120.0, 28, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0}, /* voiced bands from mid-frame */
        {60.0, 56, 40, 6.0f, MBE_FRAME_MARK_NONE, 0},        /* the longest period: responses span frames */
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
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 0},
        {56.0, 56, 40, 16.0f, MBE_FRAME_MARK_NONE, 0},  /* each subframe's period 1.64 times the last: octaves down */
        {170.0, 22, 20, 16.0f, MBE_FRAME_MARK_NONE, 0}, /* 0.57 times: octaves up */
        {70.8, 52, 40, 16.0f, MBE_FRAME_MARK_NONE, 0},  /* 1.55 times: a glide */
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
        shape_model(&model, f);
        mbe_set_frame_mark(&model, frames[f].mark);
        cur = last;
        cur.w0 = model.w0;
        cur.L = model.L;
        memcpy(cur.Vl, model.Vl, sizeof(cur.Vl));
        copy_amplitudes(&cur, &model);
        const enum mbe_frame_mark mark = reference_rendered_mark(&cur, frames[f].mark, rendering);
        mbe_set_frame_mark(&cur, mark);
        mbe_subframe_own_model(&last, &models[0]);
        /* A late onset: in the rate-33 rendering whatever the frame before,
         * in the others after a frame with no voiced band. */
        const int late = mark == MBE_FRAME_MARK_LATE_ONSET
                         && (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 || !reference_voiced(&last));
        mbe_subframe_mid_model(&last, &cur, late, rendering, &models[1]);
        mbe_subframe_own_model(&cur, &models[2]);
        /* The mid-frame pitch is the reference's own: bit for bit where IEEE
         * arithmetic gives the library and this test the same bits, within an
         * ulp elsewhere. So is its harmonic count, computed here from that
         * pitch (where the pitch may be an ulp off, away from a whole number
         * of harmonics only). */
        const float mid_pitch = reference_mid_pitch(&last, &cur, rendering);
        assert(same_bits(&models[1].w0, &mid_pitch, 1) || (inexact && within_an_ulp(models[1].w0, mid_pitch)));
        assert(last.L >= 1 && last.L <= 56 && cur.L >= 1 && cur.L <= 56);
        const int mid_count =
            (mark == MBE_FRAME_MARK_PREVIOUS_HELD)
                ? last.L
                : reference_mid_count(mid_pitch, (last.L < cur.L) ? last.L : cur.L, (last.L < cur.L) ? cur.L : last.L);
        assert(models[1].L == mid_count || (inexact && count_margin(mid_pitch) < 1.0e-5));
        assert(mid_count >= 1 && mid_count <= 56);
        models[1].w0 = mid_pitch;
        models[1].L = mid_count;
        /* Its voicing is the reference's own, which the library's must equal. */
        if (mark != MBE_FRAME_MARK_PREVIOUS_HELD) {
            unsigned char voiced[57] = {0};
            reference_mid_voicing(&last, &cur, (double)mid_pitch, models[1].L, rendering, voiced);
            models[1].any_voiced = 0;
            for (int l = 1; l <= models[1].L; ++l) {
                assert(models[1].voiced[l] == voiced[l]);
                models[1].voiced[l] = voiced[l];
                models[1].any_voiced |= voiced[l] == 1;
            }
        }
        /* The phases: each frame's own from its log2Ml, the mid-frame model's
         * from the two frames' met at its frequencies, except where it is the
         * last frame's own model held (and, from silence, this frame's). */
        double phases[3][57];
        reference_own_phases(&last, phases[0]);
        reference_own_phases(&cur, phases[2]);
        if (mark == MBE_FRAME_MARK_PREVIOUS_HELD) {
            memcpy(phases[1], frames[f].from_silence ? phases[2] : phases[0], sizeof(phases[1]));
        } else {
            reference_mid_phases(&last, &cur, (double)models[1].w0, models[1].L, phases[1]);
        }

        /* The train starts afresh after a frame with no voiced band, unless that
         * frame lacked a pitch or held the model before it to mid-frame. A
         * frame that starts from silence plays nothing of the last frame's
         * own model and starts a new train, whatever came before. */
        if (frames[f].from_silence) {
            memset(&models[0].voiced, 2, sizeof(models[0].voiced));
            models[0].any_voiced = 0;
            /* A hold continues what was heard: from silence there is nothing
             * to hold, and the frame's own model stands at mid-frame. */
            if (mark == MBE_FRAME_MARK_PREVIOUS_HELD) {
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
        reference_render(want, start, &train, &models[0], phases[0], &models[1], phases[1]);
        reference_advance(&train, start + 80.0, &models[1], &models[2]);
        reference_render(want + 80, start + 80.0, &train, &models[1], phases[1], &models[2], phases[2]);

        next_frame_as(out, &last, &model, frames[f].from_silence, rendering);
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

/* How expected_phase() treats a subframe's step. */
enum octave_rule {
    GLIDE,        /* the linear crossing alone */
    OCTAVE_RULE,  /* the octave step where the period ratio asks for it */
    OCTAVE_FORCED /* the octave step whatever the ratio: halved for a longer period, doubled for a shorter */
};

/* The pulse phase a voiced frame at `to` Hz leaves after one at `from` Hz
 * that left `phase`: two subframes, the mid-frame pitch the geometric mean;
 * each the octave step first as `octaves` says (the old train halved in phase
 * and rate when the new period is 1.6 to 2.4 times the old, doubled at 0.4 to
 * 0.6 times), then the linear crossing. */
static double
expected_phase(double phase, double from, double to, enum octave_rule octaves) {
    const double hz[3] = {from, sqrt(from * to), to};
    for (int k = 0; k < 2; ++k) {
        double old_rate = hz[k] / 8000.0;
        const double new_rate = hz[k + 1] / 8000.0;
        const double ratio = old_rate / new_rate;
        const int forced = octaves == OCTAVE_FORCED;
        if ((octaves == OCTAVE_RULE && ratio > 1.6 && ratio < 2.4) || (forced && ratio > 1.0)) {
            phase *= 0.5;
            old_rate *= 0.5;
        } else if ((octaves == OCTAVE_RULE && ratio > 0.4 && ratio < 0.6) || (forced && ratio < 1.0)) {
            phase = fmod(2.0 * phase, 1.0);
            old_rate *= 2.0;
        }
        phase += 40.0 * (old_rate + new_rate);
        phase -= floor(phase);
    }
    return phase;
}

/* About an octave apart the pulse train does not glide: a step whose period
 * ratio in a subframe lies within 20% of 2 (or of 1/2) takes the old train an
 * octave down (up) first. Each case steps once between two steady voiced
 * models whose pitches give each subframe the ratio named; the stored phase is
 * the rule's where it applies and the glide's elsewhere, the two predictions
 * far enough apart to tell; and advancing the frame unrendered leaves the
 * train exactly as rendering it. */
static void
test_octave_steps(void) {
    static const struct {
        double ratio; /* new period over old, per subframe */
        double from;  /* Hz */
        int octave;   /* whether the rule applies */
    } cases[] = {
        {1.70, 245.0, 1}, {0.55, 104.0, 1}, {1.55, 245.0, 0}, {2.45, 290.0, 0}, {0.62, 104.0, 0}, {0.39, 61.0, 0},
    };

    float out[160];
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        const double from = cases[i].from;
        const double to = from / (cases[i].ratio * cases[i].ratio);
        mbe_parms last, model;
        fresh_stream(&last);
        set_model(&model, from, (int)(3800.0 / from) < 56 ? (int)(3800.0 / from) : 56, 9, 16.0f);
        for (int f = 0; f < 3; ++f) {
            next_frame(out, &last, &model);
        }
        /* A stored phase where an octave step and a glide part far enough:
         * each steady frame before moves it on. */
        for (int f = 0; f < 8
                        && fabs(cycles_apart(expected_phase(pulse_phase(&last), from, to, OCTAVE_FORCED),
                                             expected_phase(pulse_phase(&last), from, to, GLIDE)))
                               < 0.05;
             ++f) {
            next_frame(out, &last, &model);
        }
        const double before = pulse_phase(&last);
        const double rule = expected_phase(before, from, to, OCTAVE_RULE);
        const double glide = expected_phase(before, from, to, GLIDE);
        const double forced = expected_phase(before, from, to, OCTAVE_FORCED);
        assert(fabs(cycles_apart(forced, glide)) >= 0.05);
        assert(fabs(cycles_apart(rule, cases[i].octave ? forced : glide)) < 1e-12);
        mbe_parms step;
        set_model(&step, to, (int)(3800.0 / to) < 56 ? (int)(3800.0 / to) : 56, 9, 16.0f);
        mbe_parms advanced = with_model(&last, &step);
        mbe_advance_subframes(&advanced, &last, 0, MBE_VOICED_SYNTHESIS_DVSI);
        next_frame(out, &last, &step);
        const double want = cases[i].octave ? rule : glide;
        assert(fabs(cycles_apart(pulse_phase(&last), want)) < 1e-5);
        assert(stored_phase_in_range(&last));
        assert(all_finite(out, 160));
        assert(same_bits(last.noiseOverlap, advanced.noiseOverlap, 96));
        assert(same_bits(&last.PSIl[1], &advanced.PSIl[1], 56) && same_bits(&last.PHIl[1], &advanced.PHIl[1], 56));
    }
}

/*
 * The pulse phase a frame stores lies in [0, 1) as stored: one whose last
 * interval ends within float's spacing below a whole cycle is stored as the
 * start of the next (0), whatever precision the arithmetic before the store
 * kept (x87 may compare the unrounded value with 1 and then store 1.0f). The
 * case: a voiced frame at w0 0.2 that left a phase of 0.9070417284965515,
 * then one at w0 0.1 marked MBE_FRAME_MARK_TONE (its own model at mid-frame):
 * the first interval an octave down (period ratio 2), the train reaches
 * 0.9999999917 cycle at the frame's end, which rounds to 1.0f. Rendered and
 * advanced unrendered alike.
 */
static void
test_stored_phase_range(void) {
    float out[160];
    mbe_parms last, model, step;
    fresh_stream(&last);
    set_model(&model, 100.0, 15, 15, 16.0f);
    model.w0 = 0.2f;
    for (int f = 0; f < 3; ++f) {
        next_frame(out, &last, &model);
    }
    last.noiseOverlap[PHASE_SLOT] = 0.9070417284965515f;
    set_model(&step, 100.0, 31, 31, 16.0f);
    step.w0 = 0.1f;
    mbe_set_frame_mark(&step, MBE_FRAME_MARK_TONE);
    mbe_parms advanced = with_model(&last, &step);
    mbe_set_frame_mark(&advanced, MBE_FRAME_MARK_TONE);
    mbe_advance_subframes(&advanced, &last, 0, MBE_VOICED_SYNTHESIS_DVSI);
    next_frame(out, &last, &step);
    assert(fabs(cycles_apart(pulse_phase(&last), 0.0)) < 1e-6);
    assert(stored_phase_in_range(&last));
    assert(stored_phase_in_range(&advanced));
    assert(same_bits(last.noiseOverlap, advanced.noiseOverlap, 96));
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
    mbe_subframe_mid_model(before, frame, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
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
    mbe_subframe_mid_model(&last, &unvoiced, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(same_bits(&mid.w0, &voiced.w0, 1));
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
    mbe_subframe_mid_model(&last, &unvoiced, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
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
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
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
        set_amplitude(&prev, l, 10.0f * (float)l);
    }
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    const double position = 3.0 * w0 / (double)prev.w0;
    const int below = (int)position;
    const double at = exp2(log2(10.0 * below) + (position - below) * (log2(10.0 * (below + 1)) - log2(10.0 * below)));
    assert(fabs((double)mid.magnitude[3] / sqrt(at * 20.0) - 1.0) < 1e-4);

    /* Only one frame with a voiced band: its pitch. Neither: between them. */
    set_model(&prev, 150.0, 24, 0, 10.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(same_bits(&mid.w0, &cur.w0, 1) && mid.any_voiced);
    set_model(&prev, 150.0, 24, 24, 10.0f);
    set_model(&cur, 180.0, 20, 0, 20.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(same_bits(&mid.w0, &prev.w0, 1) && !mid.any_voiced);
    set_model(&prev, 150.0, 24, 0, 10.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(fabs((double)mid.w0 - w0) < 1e-6 && !mid.any_voiced);
    assert(fabs((double)mid.magnitude[5] / sqrt(200.0) - 1.0) < 1e-4);

    /* A late onset has no voiced band at mid-frame; its bands are noise. */
    set_model(&cur, 180.0, 20, 20, 20.0f);
    mbe_subframe_mid_model(&prev, &cur, 1, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(!mid.any_voiced && mid.voiced[1] == 0 && mid.magnitude[1] > 0.0f);

    /* A band silent in either frame is silent at mid-frame. */
    prev.Ml[4] = 0.0f;
    prev.Ml[5] = 0.0f;
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
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
    assert(same_bits(&own.w0, &cur.w0, 1) && own.L == 20 && own.any_voiced);
    assert(own.voiced[2] == 1 && own.voiced[3] == 0 && own.voiced[4] > 1 && own.magnitude[2] == 20.0f);
}

/* The phase of the component at `hz` in the last 256 samples of ten frames,
 * against the frame's end, as test_steady_model_is_its_harmonics() reads it. */
static double
end_phase(const float* frames, double hz) {
    return component_phase(&frames[1344], 256, hz, 256);
}

static double
phase_apart(double a, double b) {
    return fabs(atan2(sin(a - b), cos(a - b)));
}

/*
 * A voiced harmonic's phase comes from the frame's decoded log2 amplitudes
 * (log2Ml), not from the amplitudes it renders (Ml), as DVSI's decoders
 * regenerate it: two steady streams with the same log2Ml, one rendered with a
 * high-band shelf (as enhancement adds one), play every harmonic at the same
 * phase and at their own amplitudes; two with the same Ml and differently
 * shaped log2Ml play the same amplitudes at other phases.
 */
static void
test_phases_follow_decoded_amplitudes(void) {
    static float out[3][10 * 160];
    mbe_parms model[3], last[3];
    struct mbe_subframe_model own[3];

    for (int k = 0; k < 3; ++k) {
        set_model(&model[k], 187.5, 18, 18, 0.0f);
        for (int l = 1; l <= 18; ++l) {
            set_amplitude(&model[k], l, (float)(12.0 * exp2(0.9 * sin(0.8 * l))));
        }
    }
    for (int l = 10; l <= 18; ++l) {
        model[1].Ml[l] *= 2.5f;
    }
    for (int l = 1; l <= 18; ++l) {
        model[2].log2Ml[l] = (float)(log2(12.0) + (0.9 * cos(1.3 * l)));
    }
    for (int k = 0; k < 3; ++k) {
        mbe_subframe_own_model(&model[k], &own[k]);
        fresh_stream(&last[k]);
        for (int f = 0; f < 10; ++f) {
            next_frame(out[k] + ((size_t)f * 160u), &last[k], &model[k]);
        }
    }
    assert(same_bits(own[0].phase, own[1].phase, 57) && !same_bits(own[0].magnitude, own[1].magnitude, 57));
    assert(same_bits(own[0].magnitude, own[2].magnitude, 57) && !same_bits(own[0].phase, own[2].phase, 57));
    assert(same_bits(last[0].PHIl, last[1].PHIl, 57));

    double moved = 0.0;
    for (int l = 1; l <= 18; ++l) {
        const double hz = 187.5 * l;
        for (int k = 0; k < 3; ++k) {
            const double want = 2.0 * (double)model[k].Ml[l] * (double)mbe_low_band_gain((float)hz);
            assert(fabs(component_amplitude(out[k] + 320, 1280, hz) / want - 1.0) < 2e-3);
            const double stored = (double)last[k].PHIl[l] + (double)mbe_low_band_phase((float)hz);
            assert(phase_apart(end_phase(out[k], hz), stored) < 0.02);
        }
        assert(phase_apart(end_phase(out[0], hz), end_phase(out[1], hz)) < 0.02);
        moved = fmax(moved, phase_apart(end_phase(out[0], hz), end_phase(out[2], hz)));
    }
    assert(moved > 0.3);
}

/* The mid-frame model's phases between frames of different pitch and harmonic
 * count: regenerated from the mean of the two frames' log2Ml, each read at the
 * band's frequency (linear between the frame's harmonics around it, the end
 * harmonic's beyond either end), whatever amplitudes the frames render. */
static void
test_mid_model_phases(void) {
    mbe_parms prev, cur, enhanced;
    struct mbe_subframe_model mid, shelved;
    double want[57];

    mbe_initMbeParms(&cur, &prev, &enhanced);
    set_model(&prev, 150.0, 24, 24, 0.0f);
    set_model(&cur, 180.0, 20, 10, 0.0f);
    for (int l = 1; l <= 24; ++l) {
        set_amplitude(&prev, l, (float)(10.0 * exp2(sin(0.5 * l))));
    }
    for (int l = 1; l <= 20; ++l) {
        set_amplitude(&cur, l, (float)(20.0 * exp2(0.7 * cos(0.9 * l))));
    }
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(mid.L == 22 && fabs((double)mid.w0 - sqrt((double)prev.w0 * (double)cur.w0)) < 1e-6);
    reference_mid_phases(&prev, &cur, (double)mid.w0, mid.L, want);
    for (int l = 1; l <= 56; ++l) {
        assert(fabs((double)mid.phase[l] - want[l]) < 1e-4);
    }
    /* Harmonic 1 lies below the current frame's first and harmonic 22 above
     * its last: its end harmonics' values hold there. */
    assert(mid.w0 < cur.w0 && 22.0f * mid.w0 > 20.0f * cur.w0);

    /* A shelf on what both frames render moves the magnitudes, not the phases. */
    for (int l = 12; l <= 24; ++l) {
        prev.Ml[l] *= 3.0f;
        cur.Ml[l] *= 3.0f;
    }
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &shelved);
    assert(same_bits(mid.phase, shelved.phase, 57) && !same_bits(mid.magnitude, shelved.magnitude, 57));
}

/* Whether a frame's models, its output and the phases it leaves are finite. */
static int
frame_is_finite(const float* out, const mbe_parms* cur, const mbe_parms* prev) {
    struct mbe_subframe_model model;
    mbe_subframe_own_model(cur, &model);
    int finite = all_finite(&model.phase[1], 56);
    mbe_subframe_mid_model(prev, cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &model);
    finite &= all_finite(&model.phase[1], 56);
    return finite && all_finite(out, 160) && all_finite(&cur->PSIl[1], 56) && all_finite(&cur->PHIl[1], 56)
           && all_finite(cur->noiseOverlap, 96) && all_finite(cur->previousUw, 256);
}

/*
 * log2Ml is a caller's state too: values no decoder produces (not numbers,
 * infinities, finite values far past any codec's) in either frame or both
 * render finite samples, models and phases, from silence or not, and leave a
 * state the next frame plays from; a value that is not a number counts as the
 * kernel's floor. Through the P25 data path too.
 */
static void
test_untrusted_log2_amplitudes(void) {
    static const uint32_t bits[] = {
        0x7FC00000u, 0xFFC00000u,           0x7F800000u, 0xFF800000u, 0x7F7FFFFFu, 0xFF7FFFFFu, 0x7149F2CAu /* 1e30 */,
        0xF149F2CAu, 0x49742400u /* 1e6 */, 0xC9742400u, 0x00000001u};
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    mbe_parms good;
    set_model(&good, 200.0, 18, 12, 16.0f);

    for (size_t i = 0; i < sizeof(bits) / sizeof(bits[0]); ++i) {
        mbe_parms bad = good;
        for (int l = 1; l <= 56; l += (i & 1u) ? 1 : 2) {
            memcpy(&bad.log2Ml[l], &bits[i], sizeof(bad.log2Ml[l]));
        }
        for (int which = 0; which < 3; ++which) {
            for (int from_silence = 0; from_silence < 2; ++from_silence) {
                mbe_parms last, cur;
                float out[160];
                fresh_stream(&last);
                next_frame(out, &last, &good);
                if (which != 1) {
                    memcpy(last.log2Ml, bad.log2Ml, sizeof(last.log2Ml));
                }
                cur = last;
                copy_amplitudes(&cur, (which != 0) ? &bad : &good);
                memset(out, 0, sizeof(out));
                mbe_synthesize_subframes(out, &cur, &last, plan, 1.0f, from_silence, MBE_VOICED_SYNTHESIS_DVSI);
                assert(frame_is_finite(out, &cur, &last));
                last = cur;
                cur = last;
                copy_amplitudes(&cur, &good);
                memset(out, 0, sizeof(out));
                mbe_synthesize_subframes(out, &cur, &last, plan, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
                assert(frame_is_finite(out, &cur, &last));
            }
        }

        /* The last synthesized frame of a P25 stream, as a caller's state can hold it. */
        mbe_parms cur, prev, enhanced;
        char imbe_d[88];
        float out[160];
        mbe_initMbeParms(&cur, &prev, &enhanced);
        memset(imbe_d, 1, sizeof(imbe_d));
        set_imbe_b0(imbe_d, 100);
        memcpy(enhanced.log2Ml, bad.log2Ml, sizeof(enhanced.log2Ml));
        for (int f = 0; f < 2; ++f) {
            const int rc = mbe_processImbe4400Dataf(out, NULL, imbe_d, &cur, &prev, &enhanced);
            assert(rc >= 0);
            (void)rc;
            assert(all_finite(out, 160));
        }
    }
    mbe_fft_plan_free(plan);

    /* Not a number: the floor, as a value below it is. */
    static const uint32_t nan_bits = 0x7FC00000u;
    struct mbe_subframe_model a, b, c;
    mbe_parms nan_state = good, floor_state = good, below = good;
    memcpy(&nan_state.log2Ml[5], &nan_bits, sizeof(nan_bits));
    floor_state.log2Ml[5] = MBE_PHASE_LOG2_FLOOR;
    below.log2Ml[5] = -40.0f;
    mbe_subframe_own_model(&nan_state, &a);
    mbe_subframe_own_model(&floor_state, &b);
    mbe_subframe_own_model(&below, &c);
    assert(same_bits(a.phase, b.phase, 57) && same_bits(b.phase, c.phase, 57));
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
 * pulse train in the same place (not the unvoiced overlap: a rendered frame
 * leaves its own subframe's second half in previousUw, an advanced one none;
 * here, with no transform plan, the voiced part alone). */
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
            set_amplitude(&prev, l, 4.0f + (float)l);
        }
        mbe_subframe_own_model(&prev, &own);
        for (int cur_voiced = 0; cur_voiced <= 15; cur_voiced += 15) {
            set_model(&cur, 250.0, 15, cur_voiced, 1.6f);
            mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PREVIOUS_HELD);
            for (int late = 0; late < 2; ++late) {
                mbe_subframe_mid_model(&prev, &cur, late, MBE_VOICED_SYNTHESIS_DVSI, &mid);
                assert(same_model(&mid, &own));
            }
        }
        /* The same frame marked as P25's and D-STAR's decoders mark one with
         * no voiced band and no pitch meets the two frames' magnitudes. */
        set_model(&cur, 250.0, 15, 0, 1.6f);
        mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PITCH_ABSENT);
        mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
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
    copy_amplitudes(&rendered, &held);
    copy_amplitudes(&advanced, &held);
    mbe_set_frame_mark(&rendered, MBE_FRAME_MARK_PREVIOUS_HELD);
    mbe_set_frame_mark(&advanced, MBE_FRAME_MARK_PREVIOUS_HELD);
    memset(out, 0, sizeof(out));
    mbe_synthesize_subframes(out, &rendered, &last, NULL, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
    mbe_advance_subframes(&advanced, &last, 0, MBE_VOICED_SYNTHESIS_DVSI);
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

/* A frame the core mutes moves the pulse train and the noise generator on
 * exactly as a rendered frame does, in each rendering, and leaves no unvoiced
 * overlap for the next, where a rendered frame with unvoiced bands leaves its
 * own subframe's second half (so the two states are not the same): with no
 * transform plan (the voiced part alone) and with one. The pitch steps
 * between two voiced frames, so that the rate-33 mid-frame pitch differs from
 * the others'. */
static void
test_muted_frame_advances(enum mbe_voiced_synthesis rendering) {
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
    copy_amplitudes(&rendered, &model);
    copy_amplitudes(&advanced, &model);
    memset(out, 0, sizeof(out));
    mbe_synthesize_subframes(out, &rendered, &last, NULL, 1.0f, 0, rendering);
    mbe_advance_subframes(&advanced, &last, 0, rendering);
    assert(same_bits(rendered.PSIl, advanced.PSIl, 57));
    assert(same_bits(rendered.PHIl, advanced.PHIl, 57));
    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
    assert(same_bits(&rendered.noiseSeed, &advanced.noiseSeed, 1));
    for (int n = 0; n < 256; ++n) {
        assert(advanced.previousUw[n] == 0.0f);
    }
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    mbe_parms whole = with_model(&last, &model);
    memset(out, 0, sizeof(out));
    mbe_synthesize_subframes(out, &whole, &last, plan, 1.0f, 0, rendering);
    mbe_fft_plan_free(plan);
    assert(same_bits(whole.PSIl, advanced.PSIl, 57) && same_bits(whole.PHIl, advanced.PHIl, 57));
    assert(same_bits(whole.noiseOverlap, advanced.noiseOverlap, 96));
    assert(same_bits(&whole.noiseSeed, &advanced.noiseSeed, 1));
    const float no_overlap[256] = {0.0f};
    assert(same_bits(whole.previousUw, no_overlap, MBE_UNVOICED_TAIL_OFFSET));
    assert(
        !same_bits(&whole.previousUw[MBE_UNVOICED_TAIL_OFFSET], &no_overlap[MBE_UNVOICED_TAIL_OFFSET], MBE_SUBFRAME));

    /* Through the core: a frame muted by the error rate leaves the pulse train
     * and the noise generator where a rendered frame leaves them, so the frame
     * after both plays the same voiced part. */
    mbe_parms muted = last, enh_a = last, enh_b = last;
    muted.w0 = model.w0;
    muted.L = model.L;
    memcpy(muted.Vl, model.Vl, sizeof(model.Vl));
    copy_amplitudes(&muted, &model);
    mbe_parms heard = muted;
    muted.errorRate = 1.0f;
    muted.mutingThreshold = 0.0875f;
    heard.errorRate = 0.0f;
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &muted, &enh_a, 0.0f, MBE_MUTE_NOISE_SPEC, rendering,
                                       MBE_AMPLITUDE_LIMIT_NONE);
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &heard, &enh_b, 0.0f, MBE_MUTE_NOISE_SPEC, rendering,
                                       MBE_AMPLITUDE_LIMIT_NONE);
    assert(same_bits(muted.noiseOverlap, heard.noiseOverlap, 96));
    assert(same_bits(muted.PSIl, heard.PSIl, 57));
    assert(same_bits(&muted.noiseSeed, &heard.noiseSeed, 1));
    /* The pulses the pitch step leaves: the mid-frame pitch's offset in the
     * rate-33 rendering (the factor and the half-octave wave) moves them, by
     * the cycles it takes from the frame. */
    mbe_parms other = last;
    other.w0 = model.w0;
    other.L = model.L;
    memcpy(other.Vl, model.Vl, sizeof(model.Vl));
    copy_amplitudes(&other, &model);
    const enum mbe_voiced_synthesis others = (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450)
                                                 ? MBE_VOICED_SYNTHESIS_DVSI
                                                 : MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    mbe_advance_subframes(&other, &last, 0, others);
    const double wave = 0.5 * (reference_r33_wave((double)w0_of(210.0)) + reference_r33_wave((double)w0_of(190.0)));
    const double lost = (1.0 - ((1.0 - 3.0e-4) * (1.0 + wave))) * cycles_of(sqrt(210.0 * 190.0), 80.0);
    const double offset = (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) ? -lost : lost;
    assert(fabs(cycles_apart(pulse_phase(&advanced), pulse_phase(&other) + offset)) < 1e-6);
}

/* Advancing a frame unrendered leaves the pulse train, the phases and the
 * noise generator exactly where rendering it does, through the frame marks
 * that shape the mid-frame model: onsets read in columns, before and inside a
 * voiced run, late onsets after unvoiced and voiced frames, and the mid-frame
 * columns of rate 33's rows 13 and 17. */
static void
test_advance_follows_the_marks(enum mbe_voiced_synthesis rendering) {
    static const struct {
        double hz;
        int L, voiced;
        enum mbe_frame_mark mark;
    } frames[] = {
        {90.0, 38, 0, MBE_FRAME_MARK_NONE},
        {127.5, 28, 23, MBE_FRAME_MARK_ONSET_COLUMNS},
        {127.5, 28, 23, MBE_FRAME_MARK_ONSET_COLUMNS},
        {90.0, 38, 0, MBE_FRAME_MARK_NONE},
        {160.0, 25, 3, MBE_FRAME_MARK_ONSET_COLUMNS},
        {160.0, 25, 9, MBE_FRAME_MARK_ONSET_COLUMNS},
        {90.0, 38, 0, MBE_FRAME_MARK_NONE},
        {150.0, 24, 3, MBE_FRAME_MARK_LATE_ONSET},
        {150.0, 24, 20, MBE_FRAME_MARK_NONE},
        {170.0, 22, 8, MBE_FRAME_MARK_NONE},
        {165.0, 23, 0, MBE_FRAME_MARK_MID_COLUMN_0},
        {150.0, 24, 9, MBE_FRAME_MARK_MID_COLUMNS_0_1},
        {140.0, 26, 10, MBE_FRAME_MARK_LATE_ONSET},
        {130.0, 28, 0, MBE_FRAME_MARK_MID_COLUMN_0},
        {125.0, 29, 10, MBE_FRAME_MARK_LATE_ONSET},
    };

    mbe_parms last, model;
    float out[160];
    fresh_stream(&last);
    for (size_t f = 0; f < sizeof(frames) / sizeof(frames[0]); ++f) {
        set_model(&model, frames[f].hz, frames[f].L, frames[f].voiced, 16.0f);
        mbe_set_frame_mark(&model, frames[f].mark);
        mbe_parms advanced = with_model(&last, &model);
        mbe_set_frame_mark(&advanced, frames[f].mark);
        mbe_advance_subframes(&advanced, &last, 0, rendering);
        next_frame_as(out, &last, &model, 0, rendering);
        assert(same_bits(last.PSIl, advanced.PSIl, 57));
        assert(same_bits(last.PHIl, advanced.PHIl, 57));
        assert(same_bits(last.noiseOverlap, advanced.noiseOverlap, 96));
        assert(same_bits(&last.noiseSeed, &advanced.noiseSeed, 1));
        assert(all_finite(out, 160));
    }
}

/* Whether a is b: bit for bit, or within an ulp where `tolerant` (the build
 * mode is passed in, not tested around the only use of a value). */
static int
same_or_an_ulp(float a, float b, int tolerant) {
    return same_bits(&a, &b, 1) || (tolerant && within_an_ulp(a, b));
}

/* The mid-frame models of the two renderings between prev and cur, and the
 * rate-33 one's pitch over the other's. */
static double
rate33_mid_pitch_ratio(const mbe_parms* prev, const mbe_parms* cur, struct mbe_subframe_model* common,
                       struct mbe_subframe_model* rate33) {
    mbe_subframe_mid_model(prev, cur, 0, MBE_VOICED_SYNTHESIS_DVSI, common);
    mbe_subframe_mid_model(prev, cur, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, rate33);
    assert(common->any_voiced == rate33->any_voiced);
    return (double)rate33->w0 / (double)common->w0;
}

/* A voiced model at rate-33 pitch code b0, with the pitch and harmonic count
 * the rate-33 decoder stores for it. */
static void
rate33_code_model(mbe_parms* mp, int b0) {
    int b[9] = {b0, 0, 20, 100, 30, 10, 5, 5, 3};
    char ambe_d[49];
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    mbe_ambe2450_pack(b, ambe_d);
    const int rc = mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev);
    assert(rc == MBE_AMBE2450_FRAME_VOICE);
    (void)rc;
    set_model(mp, 100.0, cur.L, cur.L, 10.0f);
    mp->w0 = cur.w0;
}

/* The mid-frame pitch of each rendering between two voiced frames at rate-33
 * pitch codes a and b: in the rate-33 rendering the formula
 * (reference_r33_mid_pitch(), in double) within an ulp of float, the mean of
 * the two frames' h(P) being near *wave (to 5e-5) where wave is not NULL; in
 * D-STAR's and the IMBE rendering the geometric mean in float, unchanged by
 * the wave (bit for bit in an IEEE build); returns that mean of h. */
static double
check_code_pair(int a, int b, const double* wave, int tolerant, int inexact) {
    mbe_parms prev, cur;
    struct mbe_subframe_model model;
    rate33_code_model(&prev, a);
    rate33_code_model(&cur, b);
    const double mean_wave = 0.5 * (reference_r33_wave((double)prev.w0) + reference_r33_wave((double)cur.w0));
    assert(wave == NULL || fabs(mean_wave - *wave) < 5.0e-5);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, &model);
    const double want = reference_r33_mid_pitch((double)prev.w0, (double)cur.w0);
    const float rounded = (float)want;
    assert(same_bits(&model.w0, &rounded, 1) || (inexact && within_an_ulp_of(model.w0, want)));
    /* Away from the factor alone by the wave: more than float rounding. */
    const double factor_alone = sqrt((double)prev.w0 * (double)cur.w0) * (1.0 - 3.0e-4);
    assert(fabs(mean_wave) < 1.0e-5 || fabs((double)model.w0 / factor_alone - 1.0 - mean_wave) < 1.0e-6);
    const float common = sqrtf(prev.w0 * cur.w0);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &model);
    assert(same_or_an_ulp(model.w0, common, tolerant));
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &model);
    assert(same_or_an_ulp(model.w0, common, tolerant));
    /* Equal pitches: the frame's own pitch, the wave not applied, in every
     * rendering (sqrtf(w0 * w0), which IEEE arithmetic rounds back to w0, a
     * fast-math build within an ulp). */
    for (int r = 0; r < 3; ++r) {
        static const enum mbe_voiced_synthesis renderings[3] = {
            MBE_VOICED_SYNTHESIS_DVSI, MBE_VOICED_SYNTHESIS_DVSI_IMBE, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450};

        mbe_subframe_mid_model(&cur, &cur, 0, renderings[r], &model);
        assert(same_or_an_ulp(model.w0, cur.w0, tolerant));
    }
    return mean_wave;
}

/* The wave's range over the rate-33 pitch codes, and its pitch codes nearest
 * its maximum and its minimum. */
static void
rate33_wave_extremes(int* top, int* bottom, double* high, double* low) {
    *high = -1.0;
    *low = 1.0;
    for (int b0 = 0; b0 < MBE_AMBE2450_PITCH_CODES; ++b0) {
        mbe_parms mp;
        rate33_code_model(&mp, b0);
        const double h = reference_r33_wave((double)mp.w0);
        if (h > *high) {
            *high = h;
            *top = b0;
        }
        if (h < *low) {
            *low = h;
            *bottom = b0;
        }
    }
}

/* Rate 33's mid-frame pitch (docs/dvsi-parity.md, "Rate 33: mid-frame pitch"):
 * between two voiced frames of different pitch the geometric mean times
 * 1 - 3.0e-4 times 1 + (h(Pa) + h(Pb)) / 2; every other case as in the other
 * codecs' rendering. */
static void
test_ambe2450_mid_pitch(int tolerant, int inexact) {
    struct mbe_subframe_model common, rate33;
    mbe_parms prev, cur, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);

    /* (i) Pairs of rate-33 pitch codes, the formula within an ulp. The wave
     * reaches 4.244e-5 -+ 8.440e-4 (C0 -+ sqrt(C1^2 + S1^2)); over the codes
     * its extremes lie within 1e-5 of those. */
    int top = 0, bottom = 0;
    double high = 0.0, low = 0.0;
    rate33_wave_extremes(&top, &bottom, &high, &low);
    assert(high > 8.864e-4 - 1.0e-5 && high < 8.865e-4 && low < -8.015e-4 + 1.0e-5 && low > -8.016e-4);
    /* Near the maximum (both codes beside the top) and near the minimum. */
    const double near_top = check_code_pair(top, top + 1, &high, tolerant, inexact);
    const double near_bottom = check_code_pair(bottom - 1, bottom, &low, tolerant, inexact);
    assert(near_top > 8.5e-4 && near_bottom < -7.6e-4);
    (void)near_top;
    (void)near_bottom;
    static const int pairs[][2] = {{0, 1},   {1, 0},   {10, 20}, {60, 50}, {100, 119}, {119, 100},
                                   {5, 119}, {58, 59}, {47, 48}, {30, 31}, {76, 77},   {118, 119}};
    for (size_t i = 0; i < sizeof(pairs) / sizeof(pairs[0]); ++i) {
        (void)check_code_pair(pairs[i][0], pairs[i][1], NULL, tolerant, inexact);
    }

    /* Both voiced at pitches off the code grid, the pitch stepping either way. */
    for (int order = 0; order < 2; ++order) {
        set_model(&prev, order ? 180.0 : 150.0, 24, 20, 10.0f);
        set_model(&cur, order ? 150.0 : 180.0, 20, 10, 20.0f);
        const double ratio = rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
        const double mean = sqrt((double)prev.w0 * (double)cur.w0);
        assert(fabs((double)common.w0 / mean - 1.0) < 1e-6);
        assert(within_an_ulp_of(rate33.w0, reference_r33_mid_pitch((double)prev.w0, (double)cur.w0)));
        assert(fabs(ratio - reference_r33_mid_pitch((double)prev.w0, (double)cur.w0) / mean) < 2e-7);
        (void)ratio;
    }
    /* A step of the smallest size a pitch takes. */
    set_model(&prev, 150.0, 24, 20, 10.0f);
    cur = prev;
    cur.w0 = nextafterf(prev.w0, 1.0f);
    const double step_ratio = rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
    assert(fabs(step_ratio - (1.0 - 3.0e-4) * (1.0 + reference_r33_wave((double)prev.w0))) < 2e-7);
    (void)step_ratio;

    /* (ii) Equal pitches: that pitch, in both, the factor and the wave not
     * applied. Bit for bit in an IEEE build; within an ulp with fast-math:
     * the mean is sqrtf(w0 * w0), which IEEE arithmetic rounds back to w0 but
     * a fast-math build need not. */
    set_model(&cur, 150.0, 24, 12, 20.0f);
    (void)rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
    assert(same_or_an_ulp(common.w0, cur.w0, tolerant) && same_or_an_ulp(rate33.w0, cur.w0, tolerant));

    /* (iii) One frame without a voiced band: the voiced frame's pitch, in both. */
    set_model(&prev, 150.0, 24, 0, 10.0f);
    set_model(&cur, 180.0, 20, 10, 20.0f);
    (void)rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
    assert(same_bits(&common.w0, &cur.w0, 1) && same_bits(&rate33.w0, &cur.w0, 1));
    set_model(&prev, 150.0, 24, 20, 10.0f);
    set_model(&cur, 180.0, 20, 0, 20.0f);
    (void)rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
    assert(same_bits(&common.w0, &prev.w0, 1) && same_bits(&rate33.w0, &prev.w0, 1));

    /* (iv) An unvoiced frame carrying its pitch into a voiced one, and two
     * unvoiced frames: the geometric mean, unmoved. */
    set_model(&prev, 150.0, 24, 0, 10.0f);
    mbe_set_frame_mark(&prev, MBE_FRAME_MARK_PITCH_CARRIES);
    set_model(&cur, 180.0, 20, 10, 20.0f);
    (void)rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
    assert(fabs((double)common.w0 / sqrt((double)prev.w0 * (double)cur.w0) - 1.0) < 1e-6);
    assert(same_bits(&common.w0, &rate33.w0, 1));
    set_model(&cur, 180.0, 20, 0, 20.0f);
    mbe_set_frame_mark(&prev, MBE_FRAME_MARK_NONE);
    (void)rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
    assert(same_bits(&common.w0, &rate33.w0, 1));

    /* A silence frame holds the previous frame's own model, in both. */
    set_model(&prev, 150.0, 24, 20, 10.0f);
    set_model(&cur, 250.0, 15, 0, 1.6f);
    mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PREVIOUS_HELD);
    (void)rate33_mid_pitch_ratio(&prev, &cur, &common, &rate33);
    assert(same_bits(&common.w0, &prev.w0, 1) && same_bits(&rate33.w0, &prev.w0, 1));

    /* (v) The ends of the float range a caller's state can reach. The two
     * largest floats (0x7F7FFFFF and 0x7F7FFFFE), in either order, take the
     * pitch to about 1.0005 FLT_MAX before it is rounded: the rendering keeps
     * no pitch there and a mid-frame model with no voiced band. At the lower
     * end the wave lies above 1 (about 1.0005 within 0.23% of FLT_MIN), so no
     * pair falls under FLT_MIN. Every pitch the rendering keeps is positive
     * and normal. */
    static const uint32_t ends[][2] = {
        {UINT32_C(0x7F7FFFFF), UINT32_C(0x7F7FFFFE)}, {UINT32_C(0x7F7FFFFE), UINT32_C(0x7F7FFFFF)},
        {UINT32_C(0x7F7FFFFF), UINT32_C(0x7F61B1E6)}, {UINT32_C(0x7F61B1E6), UINT32_C(0x7F7FFFFF)},
        {UINT32_C(0x00800000), UINT32_C(0x00829CBC)}, {UINT32_C(0x00829CBC), UINT32_C(0x00800000)},
        {UINT32_C(0x00800000), UINT32_C(0x00800001)}, {UINT32_C(0x00800000), UINT32_C(0x7F7FFFFF)},
        {UINT32_C(0x1E3CE508), UINT32_C(0x60AD78EC)}};
    int past_the_top = 0;
    for (size_t i = 0; i < sizeof(ends) / sizeof(ends[0]); ++i) {
        set_model(&prev, 150.0, 24, 20, 10.0f);
        set_model(&cur, 180.0, 20, 10, 20.0f);
        memcpy(&prev.w0, &ends[i][0], sizeof(prev.w0));
        memcpy(&cur.w0, &ends[i][1], sizeof(cur.w0));
        mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, &rate33);
        const double want = reference_r33_mid_pitch((double)prev.w0, (double)cur.w0);
        uint32_t bits;
        memcpy(&bits, &rate33.w0, sizeof(bits));
        assert(want >= (double)FLT_MIN);
        if (want > (double)FLT_MAX) {
            past_the_top++;
            assert(bits == 0u && !rate33.any_voiced);
            for (int l = 1; l <= rate33.L; ++l) {
                assert(rate33.voiced[l] != 1);
            }
        } else {
            const float rounded = (float)want;
            assert(same_bits(&rate33.w0, &rounded, 1) || (inexact && within_an_ulp_of(rate33.w0, want)));
        }
        /* 0, or from FLT_MIN (0x00800000) up to FLT_MAX (0x7F7FFFFF). */
        assert(bits == 0u || (bits >= UINT32_C(0x00800000) && bits <= UINT32_C(0x7F7FFFFF)));
        /* And the frame renders finite samples. */
        float out[160] = {0};
        mbe_synthesize_subframes(out, &cur, &prev, NULL, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450);
        assert(all_finite(out, 160));
    }
    /* The fallback was reached, by the pairs of the two largest floats. */
    assert(past_the_top == 2);
}

/* Every ordered pair of different rate-33 pitch codes, as the decoder stores
 * their pitches: the rate-33 mid-frame model's pitch is the formula
 * (reference_r33_mid_pitch(), rounded to float once; within an ulp where the
 * build may round it otherwise) and its harmonic count the one that pitch
 * gives (reference_mid_count()); D-STAR's and the IMBE rendering keep
 * sqrtf of the product. The wave changes the harmonic count for some pairs
 * (0 -> 94 among them: 19 harmonics, 18 at the pitch without the wave), so a
 * count taken from the pitch without the wave fails here. */
static void
test_ambe2450_mid_pitch_codes(int tolerant, int inexact) {
    static mbe_parms code[MBE_AMBE2450_PITCH_CODES];
    for (int b0 = 0; b0 < MBE_AMBE2450_PITCH_CODES; ++b0) {
        rate33_code_model(&code[b0], b0);
    }
    int moved = 0, moved_0_94 = 0;
    for (int a = 0; a < MBE_AMBE2450_PITCH_CODES; ++a) {
        for (int b = 0; b < MBE_AMBE2450_PITCH_CODES; ++b) {
            if (a == b) {
                continue;
            }
            struct mbe_subframe_model model;
            const mbe_parms* prev = &code[a];
            const mbe_parms* cur = &code[b];
            const int shorter = (prev->L < cur->L) ? prev->L : cur->L;
            const int longer = (prev->L < cur->L) ? cur->L : prev->L;
            mbe_subframe_mid_model(prev, cur, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, &model);
            const double want = reference_r33_mid_pitch((double)prev->w0, (double)cur->w0);
            const float rounded = (float)want;
            assert(same_bits(&model.w0, &rounded, 1) || (inexact && within_an_ulp_of(model.w0, want)));
            const int count = reference_mid_count(rounded, shorter, longer);
            assert(model.L == count || (inexact && count_margin(rounded) < 1.0e-5));
            /* The count at the pitch without the wave: the geometric mean
             * times the stored factor. */
            const float plain = (float)(sqrt((double)prev->w0 * (double)cur->w0) * reference_r33_factor());
            if (reference_mid_count(plain, shorter, longer) != count && count_margin(rounded) > 1.0e-5) {
                moved++;
                moved_0_94 |= (a == 0 && b == 94);
            }
            const float common = sqrtf(prev->w0 * cur->w0);
            mbe_subframe_mid_model(prev, cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &model);
            assert(same_or_an_ulp(model.w0, common, tolerant));
            mbe_subframe_mid_model(prev, cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &model);
            assert(same_or_an_ulp(model.w0, common, tolerant));
        }
    }
    assert(moved > 0 && moved_0_94);
    (void)moved;
    (void)moved_0_94;
}

/* The frame a parameter set holds, given a model and mark. */
static void
take_model(mbe_parms* mp, const mbe_parms* model, enum mbe_frame_mark mark) {
    mp->w0 = model->w0;
    mp->L = model->L;
    memcpy(mp->Vl, model->Vl, sizeof(mp->Vl));
    copy_amplitudes(mp, model);
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
    mbe_synthesize_subframes(held_out, &held, loud, plan, 1.0f, 1, MBE_VOICED_SYNTHESIS_DVSI);

    /* The same frame unmarked, from silence, after its own model. */
    own = *loud;
    take_model(&own, &silence, MBE_FRAME_MARK_NONE);
    plain = *loud;
    plain.noiseSeed = 4321.0f;
    take_model(&plain, &silence, MBE_FRAME_MARK_NONE);
    memset(plain_out, 0, sizeof(plain_out));
    mbe_synthesize_subframes(plain_out, &plain, &own, plan, 1.0f, 1, MBE_VOICED_SYNTHESIS_DVSI);

    /* A steady silence frame: after its own model, continuing. */
    steady = plain;
    steady.noiseSeed = 4321.0f;
    memset(steady_out, 0, sizeof(steady_out));
    mbe_synthesize_subframes(steady_out, &steady, &plain, plan, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);

    /* Continuing (no option) the mark holds the loud model to mid-frame. */
    continued = *loud;
    continued.noiseSeed = 4321.0f;
    take_model(&continued, &silence, MBE_FRAME_MARK_PREVIOUS_HELD);
    memset(continued_out, 0, sizeof(continued_out));
    mbe_synthesize_subframes(continued_out, &continued, loud, plan, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
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
 *  - advancing (a frame the core mutes) leaves the pulse train and the noise
 *    generator where rendering leaves them, through the core as well (not
 *    the unvoiced overlap, which an advanced frame does not leave).
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
    mbe_synthesize_subframes(rising, &a, &with_overlap, plan, 1.0f, 1, MBE_VOICED_SYNTHESIS_DVSI);
    mbe_synthesize_subframes(other, &b, &without_overlap, plan, 1.0f, 1, MBE_VOICED_SYNTHESIS_DVSI);
    assert(same_bits(rising, other, 160));
    assert(same_bits(&a.noiseSeed, &b.noiseSeed, 1) && !same_bits(&a.noiseSeed, &unvoiced.noiseSeed, 1));
    assert(same_bits(a.previousUw, b.previousUw, 256));
    /* Continuing, the same overlap plays. */
    a = unvoiced;
    memset(other, 0, sizeof(other));
    mbe_synthesize_subframes(other, &a, &with_overlap, plan, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
    assert(!same_bits(rising, other, 160));
    mbe_fft_plan_free(plan);

    /* A silence frame after loud voiced speech (no_train): from silence it
     * rises to its own model, exactly as the unmarked frame does from silence
     * after its own model (whose ordinary mid-frame model is its own), from the
     * same noise state; its level stays that of a steady silence frame. */
    held_from_silence(&no_train, peak);

    /* Advancing and rendering leave the pulse train and the noise generator in
     * the same place, from silence too (not the unvoiced overlap: a rendered
     * frame leaves its own subframe's second half, an advanced one none;
     * here, with no transform plan, the voiced part alone). */
    mbe_parms rendered = carried, advanced = carried;
    set_model(&model, 210.0, 17, 9, 16.0f);
    rendered.w0 = advanced.w0 = model.w0;
    rendered.L = advanced.L = model.L;
    memcpy(rendered.Vl, model.Vl, sizeof(model.Vl));
    memcpy(advanced.Vl, model.Vl, sizeof(model.Vl));
    copy_amplitudes(&rendered, &model);
    copy_amplitudes(&advanced, &model);
    mbe_parms continued = advanced;
    memset(out_scratch, 0, sizeof(out_scratch));
    mbe_synthesize_subframes(out_scratch, &rendered, &carried, NULL, 1.0f, 1, MBE_VOICED_SYNTHESIS_DVSI);
    mbe_advance_subframes(&advanced, &carried, 1, MBE_VOICED_SYNTHESIS_DVSI);
    mbe_advance_subframes(&continued, &carried, 0, MBE_VOICED_SYNTHESIS_DVSI);
    assert(same_bits(rendered.PSIl, advanced.PSIl, 57));
    assert(same_bits(rendered.PHIl, advanced.PHIl, 57));
    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
    assert(same_bits(&rendered.noiseSeed, &advanced.noiseSeed, 1));
    assert(!same_bits(continued.noiseOverlap, advanced.noiseOverlap, 96));

    /* Through the core: a frame it mutes moves the pulse train and the noise
     * generator on as the rendered frame does. */
    mbe_parms muted = carried, heard = carried, enh_a = carried, enh_b = carried;
    muted.w0 = heard.w0 = model.w0;
    muted.L = heard.L = model.L;
    memcpy(muted.Vl, model.Vl, sizeof(model.Vl));
    memcpy(heard.Vl, model.Vl, sizeof(model.Vl));
    copy_amplitudes(&muted, &model);
    copy_amplitudes(&heard, &model);
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

/* Whether the library is built with fast-math (CMake defines
 * MBELIB_TEST_LIBRARY_FAST_MATH): the compiler may then compile the core once
 * per output mode and round the two copies' samples apart by a few units in
 * the last place, so the comparison below allows 1e-5 of the frame's largest
 * magnitude there; in an IEEE build it is bit for bit. */
#if defined(MBELIB_TEST_LIBRARY_FAST_MATH) || defined(__FAST_MATH__) || defined(_M_FP_FAST)
#define LIBRARY_FAST_MATH 1
#else
#define LIBRARY_FAST_MATH 0
#endif

/* Whether two floats agree: bit for bit, or, where `tolerant`, within
 * `bound`. The build mode is passed in rather than tested where the bound is
 * computed, so that the bound is computed in both modes. */
static int
agree(const float* a, const float* b, double bound, int tolerant) {
    return same_bits(a, b, 1) || (tolerant && fabs((double)*a - (double)*b) <= bound);
}

static int
same_or_close(const float* a, const float* b, int n) {
    double peak = 0.0;
    for (int i = 0; i < n; ++i) {
        peak = fmax(peak, fabs((double)b[i]));
    }
    const double bound = 1e-5 * (1.0 + peak);
    for (int i = 0; i < n; ++i) {
        if (!agree(&a[i], &b[i], bound, LIBRARY_FAST_MATH)) {
            return 0;
        }
    }
    return 1;
}

/* The pulse state in noiseOverlap: the pulse count (element 0) and the tag
 * saying which rendering wrote the state (MBE_SUBFRAME_STATE_SLOT) are
 * discrete values stored as floats, compared bit for bit in every build; the
 * pulse times and the pulse phase (elements 1 to PHASE_SLOT) as
 * same_or_close() compares them; the elements between are zero, bit for bit. */
static int
pulse_state_same_or_close(const float* a, const float* b) {
    return same_bits(&a[0], &b[0], 1) && same_bits(&a[MBE_SUBFRAME_STATE_SLOT], &b[MBE_SUBFRAME_STATE_SLOT], 1)
           && same_or_close(&a[1], &b[1], PHASE_SLOT)
           && same_bits(&a[PHASE_SLOT + 1], &b[PHASE_SLOT + 1], MBE_SUBFRAME_STATE_SLOT - PHASE_SLOT - 1);
}

/*
 * The output-limit option of the synthesis core: a frame left unlimited is
 * the limited frame before its limit, sample for sample (limiting it gives the
 * limited frame's bits), for a voiced, an unvoiced and a from-silence frame,
 * loud enough to pass the limit; and the state the two leave is the same
 * (see same_or_close()).
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
        assert(same_or_close(limited, unlimited, 160));
        assert(same_or_close(a.Ml, b.Ml, 57) && same_or_close(a.PSIl, b.PSIl, 57) && same_or_close(a.PHIl, b.PHIl, 57));
        assert(pulse_state_same_or_close(a.noiseOverlap, b.noiseOverlap)
               && same_or_close(a.previousUw, b.previousUw, 256));
        assert(same_bits(&a.noiseSeed, &b.noiseSeed, 1) && same_or_close(&a.localEnergy, &b.localEnergy, 1));
        assert(memcmp(a.Vl, b.Vl, sizeof(a.Vl)) == 0 && a.amplitudeThreshold == b.amplitudeThreshold);
    }
}

/* The next frame of a stream rendered whole (voiced and unvoiced, with `plan`),
 * at unit gain, from `last`, which receives the new state. */
static void
next_whole_frame(float out[160], mbe_parms* last, const mbe_parms* model, mbe_fft_plan* plan) {
    mbe_parms cur = with_model(last, model);
    memset(out, 0, 160 * sizeof(float));
    mbe_synthesize_subframes(out, &cur, last, plan, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
    *last = cur;
}

/* Every 7th of `n` overlap samples from `start`: `value` in one copy of a
 * state, 0 in the other. */
static void
corrupt_overlap(float* bad, float* clean, int start, int n, uint32_t value) {
    for (int k = start; k < start + n; k += 7) {
        memcpy(&bad[k], &value, sizeof(value));
        clean[k] = 0.0f;
    }
}

/* Where a huge finite sample of previousUw plays as stored: at every 7th
 * output sample from `first` up to `last`, the output of the state that holds
 * it there is `want`, and that of the state with 0 there is not, at one of
 * those samples at least (a loud frame reaches the limit at some anyway). */
static void
assert_plays_huge_overlap(const float bad[160], const float clean[160], int first, int last, float want) {
    int apart = 0;
    for (int n = first; n <= last; n += 7) {
        assert(same_bits(&bad[n], &want, 1));
        apart += !same_bits(&clean[n], &want, 1);
    }
    assert(apart > 0);
}

/* The output limit of the codec paths and mbe_synthesizeSpeechf(), with the
 * sign of `bits`. */
static float
limit_with_sign_of(uint32_t bits) {
    float limit = ((bits >> 31) != 0u) ? -FLT_MAX : FLT_MAX;
    mbe_clipFloatBuffer(&limit, 1);
    return limit;
}

/* The TIA-102.BABA rendering (mbe_synthesizeSpeechf()) of a third frame of
 * `model`, from the state two frames leave (so that the noise overlap is drawn
 * noise) with `value` in every 7th sample of its previousUw (from 128) or of
 * its noise overlap, and from the same state with 0 there. */
static void
spec_frame_after(float out[2][160], const mbe_parms* model, int noise_overlap, uint32_t value) {
    mbe_parms spec, unused, spec_prev;
    mbe_initMbeParms(&spec, &spec_prev, &unused);
    for (int f = 0; f < 2; ++f) {
        mbe_parms next = with_model(&spec, model);
        mbe_synthesizeSpeechf(out[0], &next, &spec);
        spec = next;
    }
    mbe_parms spec_bad = spec, spec_clean = spec;
    if (noise_overlap) {
        corrupt_overlap(spec_bad.noiseOverlap, spec_clean.noiseOverlap, 0, 96, value);
    } else {
        corrupt_overlap(spec_bad.previousUw, spec_clean.previousUw, 128, 128, value);
    }
    for (int k = 0; k < 2; ++k) {
        mbe_parms* state = (k == 0) ? &spec_bad : &spec_clean;
        mbe_parms next = with_model(state, model);
        mbe_synthesizeSpeechf(out[k], &next, state);
    }
}

/*
 * The unvoiced overlap a frame plays from the state before it is a caller's
 * state too. The subframe rendering's overlap (previousUw from
 * MBE_UNVOICED_TAIL_OFFSET) directly and through the P25 data path, and the
 * TIA-102.BABA rendering's (mbe_synthesizeSpeechf(): previousUw and the noise
 * overlap), each with not a number of either sign, plus and minus infinity,
 * plus and minus FLT_MAX and 2^127 in it:
 *  - the frame is finite;
 *  - a sample of previousUw that is not a finite number plays as 0, and a huge
 *    finite one plays as stored: the subframe rendering alone, before any
 *    limit, gives the stored value at that sample (it drowns everything else
 *    there), and the limited paths (the P25 data path, mbe_synthesizeSpeechf())
 *    give the output limit of its sign there, where the state with 0 in its
 *    place does not (at some of those samples: the P25 frame is loud);
 *  - a sample of the noise overlap that the TIA-102.BABA generator cannot
 *    have written (it writes 0 to 53124) plays as 0: a huge one would
 *    otherwise overflow the transform of the noise buffer into NaN.
 */
static void
test_untrusted_unvoiced_overlap(void) {
    static const struct {
        uint32_t bits;
        int finite;
    } values[] = {
        {0x7FC00000u, 0}, {0xFFC00000u, 0}, {0x7F800000u, 0}, {0xFF800000u, 0}, {0x7F7FFFFFu, 1}, /* FLT_MAX */
        {0xFF7FFFFFu, 1},                                                                         /* -FLT_MAX */
        {0x7F000000u, 1},                                                                         /* 2^127 */
    };

    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    mbe_parms mixed, unused[2];
    mbe_initMbeParms(&mixed, &unused[0], &unused[1]);
    set_model(&mixed, 200.0, 18, 6, 16.0f); /* 12 unvoiced bands */
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        const uint32_t bits = values[i].bits;
        const int zeroed = !values[i].finite;
        float out[2][160];

        mbe_parms last;
        fresh_stream(&last);
        next_whole_frame(out[0], &last, &mixed, plan);
        assert(mbe_subframe_state_is_marked(&last));
        mbe_parms bad = last, clean = last;
        corrupt_overlap(bad.previousUw, clean.previousUw, MBE_UNVOICED_TAIL_OFFSET, MBE_SUBFRAME, bits);
        next_whole_frame(out[0], &bad, &mixed, plan);
        next_whole_frame(out[1], &clean, &mixed, plan);
        assert(all_finite(out[0], 160));
        assert(!zeroed || same_or_close(out[0], out[1], 160));
        float stored;
        memcpy(&stored, &bits, sizeof(stored));
        if (!zeroed) {
            /* Rendered as stored: at a sample this large nothing else counts. */
            assert_plays_huge_overlap(out[0], out[1], 0, MBE_SUBFRAME - 1, stored);
        }

        mbe_parms cur[2], prev[2], enhanced[2];
        char imbe_d[88];
        memset(imbe_d, 1, sizeof(imbe_d));
        set_imbe_b0(imbe_d, 100);
        mbe_initMbeParms(&cur[0], &prev[0], &enhanced[0]);
        int rc = mbe_processImbe4400Dataf(out[0], NULL, imbe_d, &cur[0], &prev[0], &enhanced[0]);
        assert(rc >= 0);
        assert(mbe_subframe_state_is_marked(&enhanced[0]));
        cur[1] = cur[0];
        prev[1] = prev[0];
        enhanced[1] = enhanced[0];
        corrupt_overlap(enhanced[0].previousUw, enhanced[1].previousUw, MBE_UNVOICED_TAIL_OFFSET, MBE_SUBFRAME, bits);
        for (int k = 0; k < 2; ++k) {
            rc = mbe_processImbe4400Dataf(out[k], NULL, imbe_d, &cur[k], &prev[k], &enhanced[k]);
            assert(rc >= 0);
        }
        (void)rc;
        assert(all_finite(out[0], 160));
        assert(!zeroed || same_or_close(out[0], out[1], 160));
        if (!zeroed) {
            assert_plays_huge_overlap(out[0], out[1], 0, MBE_SUBFRAME - 1, limit_with_sign_of(bits));
        }

        /* The TIA-102.BABA rendering: previousUw, then the noise overlap,
         * where every one of these values is outside what the generator writes. */
        spec_frame_after(out, &mixed, 0, bits);
        assert(all_finite(out[0], 160));
        assert(!zeroed || same_or_close(out[0], out[1], 160));
        if (!zeroed) {
            /* Where its overlap window is whole (WOLA, mbe_synthesisWindow()). */
            assert_plays_huge_overlap(out[0], out[1], 0, 55, limit_with_sign_of(bits));
        }
        spec_frame_after(out, &mixed, 1, bits);
        assert(all_finite(out[0], 160));
        assert(same_or_close(out[0], out[1], 160));
    }
    /* The bound: 53124, the largest sample the generator writes, plays as
     * stored; 53125 plays as 0. */
    float out[2][160];
    spec_frame_after(out, &mixed, 1, 0x474F8400u);
    assert(!same_bits(out[0], out[1], 160));
    spec_frame_after(out, &mixed, 1, 0x474F8500u);
    assert(same_or_close(out[0], out[1], 160));
    mbe_fft_plan_free(plan);
}

/*
 * mbe_generate_noise_with_overlap() directly: the head of the buffer is the
 * overlap, sample for sample, where a sample lies in the range the generator
 * writes (+0 to 53124, by its bit pattern; any value there, not only the
 * generator's integers), and +0 where it does not (the next float above 53124,
 * 53125, a negative value, -0, a NaN of either sign, either infinity, plus or
 * minus FLT_MAX). The last overlap slot holds 7 (-1 there is the subframe
 * rendering's mark, which clears the overlap). The head is copied as
 * memmove() copies, so a buffer that overlaps the overlap plays what the
 * copy gave before the guard. With the overlap at s and the buffer at s + 1,
 * the new overlap, written last, replaces buffer[0..94]; what survives of
 * the buffer (buffer[95] on), the new overlap and the seed are those of the
 * same call on separate buffers.
 */
static void
test_noise_overlap_guard(void) {
    static const struct {
        uint32_t bits;
        int kept;
    } values[] = {
        {0x00000000u, 1}, /* 0 */
        {0x00000001u, 1}, /* the smallest subnormal */
        {0x3F000000u, 1}, /* 0.5 */
        {0x3F800000u, 1}, /* 1 */
        {0x474F8400u, 1}, /* 53124 */
        {0x474F8401u, 0}, /* the next float above 53124 */
        {0x474F8500u, 0}, /* 53125 */
        {0xBF800000u, 0}, /* -1 */
        {0x80000000u, 0}, /* -0 */
        {0x7FC00000u, 0}, {0xFFC00000u, 0}, {0x7F800000u, 0}, {0xFF800000u, 0}, {0x7F7FFFFFu, 0}, /* FLT_MAX */
        {0xFF7FFFFFu, 0},                                                                         /* -FLT_MAX */
    };

    const float last_slot = 7.0f;
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        float overlap[MBE_NOISE_OVERLAP], buffer[MBE_FFT_SIZE];
        for (int k = 0; k < MBE_SUBFRAME_STATE_SLOT; ++k) {
            memcpy(&overlap[k], &values[i].bits, sizeof(values[i].bits));
        }
        overlap[MBE_SUBFRAME_STATE_SLOT] = last_slot;
        float seed = 1.0f;
        mbe_generate_noise_with_overlap(buffer, &seed, overlap);
        const uint32_t want = values[i].kept ? values[i].bits : 0u;
        for (int k = 0; k < MBE_SUBFRAME_STATE_SLOT; ++k) {
            uint32_t got;
            memcpy(&got, &buffer[k], sizeof(got));
            assert(got == want);
        }
        assert(same_bits(&buffer[MBE_SUBFRAME_STATE_SLOT], &last_slot, 1));
    }

    /* Overlapping buffers: overlap at s, buffer at s + 1. */
    float s[MBE_FFT_SIZE + 1], apart_overlap[MBE_NOISE_OVERLAP], apart_buffer[MBE_FFT_SIZE];
    for (int k = 0; k < MBE_NOISE_OVERLAP; ++k) {
        s[k] = (float)(k + 1);
        apart_overlap[k] = (float)(k + 1);
    }
    float seed = 1.0f, apart_seed = 1.0f;
    mbe_generate_noise_with_overlap(s + 1, &seed, s);
    mbe_generate_noise_with_overlap(apart_buffer, &apart_seed, apart_overlap);
    /* The buffer past the overlap's new contents (the overlap is written last). */
    assert(
        same_bits(&s[MBE_NOISE_OVERLAP], &apart_buffer[MBE_NOISE_OVERLAP - 1], MBE_FFT_SIZE - MBE_NOISE_OVERLAP + 1));
    assert(same_bits(s, apart_overlap, MBE_NOISE_OVERLAP));
    assert(same_bits(&seed, &apart_seed, 1));
    assert(s[MBE_NOISE_OVERLAP] > 95.5f && s[MBE_NOISE_OVERLAP] < 96.5f); /* buffer[95]: overlap[95] */
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
            copy_amplitudes(&cur, next);
            memset(guarded + 16, 0, 160 * sizeof(float));
            mbe_synthesize_subframes(guarded + 16, &cur, &last, NULL, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
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
            mbe_subframe_mid_model(order ? &after : &before, order ? &before : &after, 0, MBE_VOICED_SYNTHESIS_DVSI,
                                   &model);
            assert(fetestexcept(FE_OVERFLOW) == 0);
#endif
            mbe_subframe_mid_model(order ? &after : &before, order ? &before : &after, 0, MBE_VOICED_SYNTHESIS_DVSI,
                                   &model);
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
     * the state stays as it is, byte for byte (in P25's rendering too, which
     * marks the model it renders only on a frame it accepts), and the models
     * have no band. */
    static const int counts[] = {0, -1, 57, 1000, INT_MAX, INT_MIN};

    static const enum mbe_voiced_synthesis renderings[] = {
        MBE_VOICED_SYNTHESIS_DVSI, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_VOICED_SYNTHESIS_DVSI_IMBE};

    set_model(&good, 200.0, 18, 18, 16.0f);
    for (size_t k = 0; k < sizeof(counts) / sizeof(counts[0]) * 3; ++k) {
        const size_t i = k / 3;
        const enum mbe_voiced_synthesis rendering = renderings[k % 3];
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
            mbe_synthesize_subframes(guarded + 16, &cur, &prev, NULL, 1.0f, 0, rendering);
            for (int n = 0; n < 192; ++n) {
                assert(guarded[n] == 7.0f);
            }
            memcpy(after, &cur, sizeof(after));
            assert(memcmp(after, before, sizeof(after)) == 0);
            mbe_advance_subframes(&cur, &prev, 0, rendering);
            memcpy(after, &cur, sizeof(after));
            assert(memcmp(after, before, sizeof(after)) == 0);
            mbe_subframe_mid_model(&prev, &cur, 0, rendering, &model);
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
                mbe_subframe_mid_model(before, after, 0, MBE_VOICED_SYNTHESIS_DVSI, &model);
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

/* The parameter bit that voices band `band` of a P25 frame of pitch code b0
 * (harmonics 3 band - 2 .. 3 band, the last band the rest) and nothing else;
 * -1 if there is none. */
static int
imbe_band_bit(int b0, int band) {
    for (int bit = 6; bit < 85; ++bit) {
        mbe_parms cur, prev, enhanced;
        char imbe_d[88] = {0};
        set_imbe_b0(imbe_d, b0);
        imbe_d[bit] = 1;
        mbe_initMbeParms(&cur, &prev, &enhanced);
        if (mbe_decodeImbe4400Parms(imbe_d, &cur, &prev) != 0) {
            continue;
        }
        int match = 1;
        for (int l = 1; l <= cur.L; ++l) {
            const int in_band = (l > 3 * (band - 1)) && ((l <= 3 * band) || (band == cur.K));
            match &= cur.Vl[l] == in_band;
        }
        if (match) {
            return bit;
        }
    }
    return -1;
}

/* P25 parameter bits of pitch code b0 whose voicing word is `word` (band 1
 * first, one character per band, '1' voiced); the other fields 0. */
static void
set_imbe_voicing(char imbe_d[88], int b0, const char* word) {
    memset(imbe_d, 0, 88);
    set_imbe_b0(imbe_d, b0);
    for (int band = 1; word[band - 1] != '\0'; ++band) {
        if (word[band - 1] == '1') {
            const int bit = imbe_band_bit(b0, band);
            assert(bit >= 0);
            imbe_d[bit] = 1;
        }
    }
}

/* The P25 onset rule for harmonic l of a frame at its own pitch (any l: one
 * past L is voiced as harmonic L is). */
static int
onset_harmonic_sounds(const mbe_parms* mp, int l) {
    const int nearest = (l > mp->L) ? mp->L : l;
    return mp->Vl[nearest] == 1 && reference_onset_sounds(mp, (double)l * (double)mp->w0 * 8000.0 / (2.0 * M_PI));
}

/* A frame decoded from `imbe_d` after a P25 frame with no voiced band (b0 150,
 * a code DVSI's encoder keeps for such frames): its mid-frame model voices
 * exactly the harmonics `sounds` lists (0-terminated), which are the rule's;
 * the frame is marked for it. Returns whether any sounds. */
static int
check_p25_onset(const char imbe_d[88], const int* sounds) {
    mbe_parms cur, prev, enhanced, unvoiced;
    char before[88] = {0};
    struct mbe_subframe_model mid;
    set_imbe_b0(before, 150);
    mbe_initMbeParms(&unvoiced, &prev, &enhanced);
    const int kind = mbe_decodeImbe4400Parms(before, &unvoiced, &prev);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    const int kind2 = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
    assert(kind == 0 && kind2 == 0);
    (void)kind;
    (void)kind2;
    assert(mbe_frame_mark(&cur) == MBE_FRAME_MARK_ONSET_COLUMNS);
    mbe_subframe_mid_model(&unvoiced, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
    assert(same_bits(&mid.w0, &cur.w0, 1));
    int listed = 0, any = 0;
    for (int l = 1; l <= mid.L; ++l) {
        const int want = sounds[listed] == l;
        listed += want;
        assert(mid.voiced[l] == (want ? 1 : 0));
        assert(onset_harmonic_sounds(&cur, l) == want);
        any |= want;
    }
    assert(sounds[listed] == 0);
    return any;
}

/* Through the P25 data path, after two frames with no voiced band: the frame's
 * pulse train starts with the frame when a harmonic sounds at mid-frame and
 * at mid-frame when none does; ProVoice, which shares the path, alike. */
static void
check_p25_onset_train(const char imbe_d[88], int early, int provoice) {
    mbe_parms cur, prev, enhanced;
    char unvoiced[88] = {0};
    float out[160];
    set_imbe_b0(unvoiced, 150);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int f = 0; f < 3; ++f) {
        mbe_process_result result;
        mbe_initProcessResult(&result);
        result.flags |= provoice ? MBE_PROCESS_FLAG_PROVOICE : 0u;
        const int rc = mbe_processImbe4400Dataf(out, &result, (f < 2) ? unvoiced : imbe_d, &cur, &prev, &enhanced);
        assert(rc >= 0);
        (void)rc;
    }
    assert(mbe_frame_mark(&enhanced) == MBE_FRAME_MARK_ONSET_COLUMNS);
    const double hz = (double)enhanced.w0 * 8000.0 / (2.0 * M_PI);
    assert(fabs(cycles_apart(pulse_phase(&enhanced), cycles_of(hz, early ? 160.0 : 80.0))) < 1e-4);
    assert(all_finite(out, 160));
}

/* The bands a P25 frame after one with no voiced band starts with, as DVSI's
 * decoder plays them (MBE_FRAME_MARK_ONSET_COLUMNS): a frame voiced in its
 * first band only starts late at pitch code 56 (three harmonics up to 502.8
 * Hz at the pitch played) and early at 55 (508.1 Hz, the first column
 * covered: harmonics 1 and 2
 * sound, the third, in column 1, does not); the word 1000001100 early at 86
 * (127.5 Hz: column 5 covered, so harmonics 19, within 250 Hz below it, to
 * 23 sound, and neither band 1 nor harmonic 24, in column 6, does) and late at
 * 92 (121.7 Hz); 1011 then unvoiced bands early at 79 (column 2: harmonics 7
 * to 11) and late at 92, 99 and 102; and nothing above 3505 Hz (a frame
 * voiced throughout at code 0, whose 9th harmonic lies at 3.65 kHz). */
static void
test_p25_onset_columns(void) {
    static const struct {
        int b0;
        const char* word;
        int sounds[10];
    } cases[] = {
        {56, "10000000", {0}},
        {55, "1000000", {1, 2, 0}},
        {86, "1000001100", {19, 20, 21, 22, 23, 0}},
        {92, "1000001100", {0}},
        {79, "101100000", {7, 8, 9, 10, 11, 0}},
        {92, "1011000000", {0}},
        {99, "10110000000", {0}},
        {102, "10110000000", {0}},
        {0, "111", {1, 2, 3, 4, 5, 6, 7, 8, 0}},
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        char imbe_d[88];
        set_imbe_voicing(imbe_d, cases[i].b0, cases[i].word);
        const int early = check_p25_onset(imbe_d, cases[i].sounds);
        check_p25_onset_train(imbe_d, early, 0);
        check_p25_onset_train(imbe_d, early, 1);
    }
}

/* Where the onset rule does not act: after a voiced frame the mark changes
 * nothing (the mid-frame model is that of an unmarked frame, every voiced
 * band sounding), and a frame played as a P25 tone model keeps its own model
 * at mid-frame. */
static void
test_p25_onset_columns_only_at_onsets(void) {
    mbe_parms cur, prev, enhanced, voiced, unvoiced;
    char imbe_d[88], before[88] = {0};
    struct mbe_subframe_model marked, plain, own;
    set_imbe_voicing(imbe_d, 86, "1000001100");
    mbe_initMbeParms(&cur, &prev, &enhanced);
    const int kind = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
    assert(kind == 0);
    (void)kind;
    memset(before, 1, sizeof(before));
    set_imbe_b0(before, 80);
    mbe_initMbeParms(&voiced, &prev, &enhanced);
    const int kind2 = mbe_decodeImbe4400Parms(before, &voiced, &prev);
    assert(kind2 == 0);
    (void)kind2;
    mbe_parms unmarked = cur;
    mbe_set_frame_mark(&unmarked, MBE_FRAME_MARK_NONE);
    mbe_subframe_mid_model(&voiced, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &marked);
    mbe_subframe_mid_model(&voiced, &unmarked, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &plain);
    assert(same_model(&marked, &plain));
    /* Band 1, which no covered column holds, sounds at mid-frame. */
    assert(marked.voiced[1] == 1);

    /* The same frame as a tone model after a frame with no voiced band. */
    memset(before, 0, sizeof(before));
    set_imbe_b0(before, 150);
    mbe_initMbeParms(&unvoiced, &prev, &enhanced);
    const int kind3 = mbe_decodeImbe4400Parms(before, &unvoiced, &prev);
    assert(kind3 == 0);
    (void)kind3;
    mbe_set_frame_mark(&cur, MBE_FRAME_MARK_TONE);
    mbe_subframe_mid_model(&unvoiced, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &marked);
    mbe_subframe_own_model(&cur, &own);
    assert(same_model(&marked, &own));
}

/* The quantities the onset rule rounds up to a harmonic number, the column
 * edges 505 + 500 k Hz and the same edges 250 Hz lower over the frame's
 * fundamental, lie far from a whole number at every pitch code: at least
 * 2.5e-4 of a harmonic at the pitch the P25 decoder plays (2.86e-4 is the
 * least, code 3 at 1505 Hz) and 1.5e-4 at ProVoice's nominal one (1.55e-4,
 * code 36 at 2755 Hz), so that no build's arithmetic moves a harmonic across
 * an edge. The pitches are the decoder's, through its own paths. */
static void
test_p25_onset_edges(void) {
    double least[2] = {1.0, 1.0};
    for (int provoice = 0; provoice < 2; ++provoice) {
        for (int b0 = 0; b0 <= 207; ++b0) {
            mbe_parms cur, prev, enhanced;
            mbe_process_result result;
            char imbe_d[88] = {0};
            float out[160];
            set_imbe_b0(imbe_d, b0);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            mbe_initProcessResult(&result);
            result.flags |= provoice ? MBE_PROCESS_FLAG_PROVOICE : 0u;
            const int rc = mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &enhanced);
            assert(rc >= 0);
            (void)rc;
            for (int k = 0; k < 14; ++k) {
                const double hz = 505.0 + (500.0 * (k % 7)) - ((k < 7) ? 0.0 : 250.0);
                const double position = hz * M_PI / (4000.0 * (double)enhanced.w0);
                const double margin = fabs(position - floor(position + 0.5));
                least[provoice] = fmin(least[provoice], margin);
            }
        }
    }
    assert(least[0] > 2.5e-4 && least[0] < 3.0e-4);
    assert(least[1] > 1.5e-4 && least[1] < 2.0e-4);
}

/* Through the P25 data path: an onset frame muted by the error rate, which the
 * core advances, leaves the pulse train and the noise generator where the same
 * frame heard leaves them. The unvoiced overlap does not end alike (an advanced
 * frame leaves none, a rendered one the second half of its own subframe), so
 * the frame after plays otherwise in its first 80 samples, over which that
 * overlap is added, and alike from sample 80 on. */
static void
test_p25_muted_onset_advances(void) {
    mbe_parms cur[2], prev[2], enhanced[2];
    char unvoiced[88] = {0}, onset[88], voiced[88];
    float out[2][160];
    set_imbe_b0(unvoiced, 150);
    set_imbe_voicing(onset, 86, "1000001100");
    memset(voiced, 1, sizeof(voiced));
    set_imbe_b0(voiced, 86);
    for (int s = 0; s < 2; ++s) {
        mbe_initMbeParms(&cur[s], &prev[s], &enhanced[s]);
    }
    const char* frames[4] = {unvoiced, onset, voiced, voiced};
    for (int f = 0; f < 4; ++f) {
        for (int s = 0; s < 2; ++s) {
            mbe_process_result result;
            mbe_initProcessResult(&result);
            /* Stream 1's onset frame comes with errors enough to mute it: its
             * error rate is back under the limit for the frames after. */
            if (s == 1 && f == 1) {
                prev[s].errorRate = 1.0f;
            } else if (f == 2) {
                prev[s].errorRate = 0.0f;
            }
            const int rc = mbe_processImbe4400Dataf(out[s], &result, frames[f], &cur[s], &prev[s], &enhanced[s]);
            assert(rc >= 0);
            (void)rc;
            assert(((result.flags & MBE_PROCESS_FLAG_MUTE) != 0u) == (s == 1 && f == 1));
        }
        if (f >= 1) {
            assert(same_bits(enhanced[0].noiseOverlap, enhanced[1].noiseOverlap, 96));
            assert(same_bits(enhanced[0].PSIl, enhanced[1].PSIl, 57));
            assert(same_bits(&enhanced[0].noiseSeed, &enhanced[1].noiseSeed, 1));
        }
        if (f == 1) {
            const float none[256] = {0.0f};
            assert(same_bits(enhanced[1].previousUw, none, 256));
            assert(!same_bits(enhanced[0].previousUw, none, 256));
        }
        if (f == 2) {
            assert(!same_bits(out[0], out[1], 80));
            assert(same_bits(&out[0][80], &out[1][80], 80));
        }
    }
}

/*
 * The marks that shape the mid-frame voicing on a caller's state: any pitch
 * bit pattern in either frame (none, subnormal, the largest float, just under
 * and past Nyquist), harmonic counts 1 and 56, and voicing values no decoder
 * gives, render finite samples and state, rendered or advanced, from silence
 * or not, and memory-safely (the sanitizer build runs this).
 */
static void
test_untrusted_mid_voicing_marks(void) {
    static const enum mbe_frame_mark marks[] = {MBE_FRAME_MARK_ONSET_COLUMNS, MBE_FRAME_MARK_MID_COLUMN_0,
                                                MBE_FRAME_MARK_MID_COLUMNS_0_1, MBE_FRAME_MARK_LATE_ONSET};

    static const uint32_t w0_bits[] = {0x00800000u, 0x00000001u, 0x7F7FFFFFu, 0x40490FDAu, 0x40490FDCu,
                                       0x3C000000u, 0x7FC00000u, 0x7F800000u, 0x3E800000u};
    static const int counts[] = {1, 56, 9};
    static const int votes[] = {0, 1, 2, -1, 1000000};
    for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
        for (size_t i = 0; i < sizeof(w0_bits) / sizeof(w0_bits[0]); ++i) {
            for (size_t c = 0; c < sizeof(counts) / sizeof(counts[0]); ++c) {
                for (int v = 0; v < 5; ++v) {
                    mbe_parms last, model;
                    float out[160];
                    fresh_stream(&last);
                    set_model(&model, 110.0, 30, 0, 8.0f);
                    next_frame(out, &last, &model);
                    set_model(&model, 200.0, counts[c], counts[c], 16.0f);
                    memcpy(&model.w0, &w0_bits[i], sizeof(model.w0));
                    for (int l = 1; l <= 56; ++l) {
                        model.Vl[l] = ((l + v) % 3 == 0) ? votes[v] : votes[(v + l) % 5];
                    }
                    mbe_set_frame_mark(&model, marks[m]);
                    for (int start = 0; start < 2; ++start) {
                        for (int rendering = 0; rendering < 3; ++rendering) {
                            static const enum mbe_voiced_synthesis kinds[3] = {MBE_VOICED_SYNTHESIS_DVSI,
                                                                               MBE_VOICED_SYNTHESIS_DVSI_AMBE2450,
                                                                               MBE_VOICED_SYNTHESIS_DVSI_IMBE};
                            const enum mbe_voiced_synthesis kind = kinds[rendering];
                            mbe_parms rendered = last, advanced = with_model(&last, &model);
                            mbe_set_frame_mark(&advanced, marks[m]);
                            next_frame_as(out, &rendered, &model, start, kind);
                            assert(all_finite(out, 160) && all_finite(rendered.noiseOverlap, 96));
                            assert(all_finite(&rendered.PSIl[1], 56) && all_finite(&rendered.PHIl[1], 56));
                            mbe_advance_subframes(&advanced, &last, start, kind);
                            assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
                            /* And a frame on from it, the marked frame now the one before. */
                            set_model(&model, 150.0, 24, 12, 16.0f);
                            mbe_set_frame_mark(&model, marks[m]);
                            next_frame_as(out, &rendered, &model, 0, kind);
                            assert(all_finite(out, 160));
                            set_model(&model, 200.0, counts[c], counts[c], 16.0f);
                            memcpy(&model.w0, &w0_bits[i], sizeof(model.w0));
                            for (int l = 1; l <= 56; ++l) {
                                model.Vl[l] = ((l + v) % 3 == 0) ? votes[v] : votes[(v + l) % 5];
                            }
                            mbe_set_frame_mark(&model, marks[m]);
                        }
                    }
                }
            }
        }
    }
}

/* The P25 and D-STAR decoders mark a frame with no voiced band by its pitch
 * code, and P25 a voiced frame for the onset columns, replacing whatever mark
 * the caller's state held. */
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
            const enum mbe_frame_mark want =
                voiced ? MBE_FRAME_MARK_ONSET_COLUMNS : (own_code ? MBE_FRAME_MARK_NONE : MBE_FRAME_MARK_PITCH_ABSENT);
            assert(mbe_frame_mark(&cur) == want);
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

/* A frame of `model` marked `mark` after `last`, rendered the way `rendering`
 * says (the voiced part alone): its samples and the state it leaves. */
static void
render_marked(float out[160], mbe_parms* state, const mbe_parms* last, const mbe_parms* model, enum mbe_frame_mark mark,
              enum mbe_voiced_synthesis rendering) {
    mbe_parms marked = *model;
    *state = *last;
    mbe_set_frame_mark(&marked, mark);
    next_frame_as(out, state, &marked, 0, rendering);
}

/* Whether two frames rendered the same: their samples, phases and pulse train. */
static int
same_rendering(const float a[160], const mbe_parms* a_state, const float b[160], const mbe_parms* b_state) {
    return same_bits(a, b, 160) && same_bits(&a_state->PSIl[1], &b_state->PSIl[1], 56)
           && same_bits(&a_state->PHIl[1], &b_state->PHIl[1], 56)
           && same_bits(a_state->noiseOverlap, b_state->noiseOverlap, 96);
}

/*
 * A frame mark acts only in the rendering of the codec whose decoder writes
 * it. Rate 33's (MBE_FRAME_MARK_MID_COLUMN_0, MBE_FRAME_MARK_MID_COLUMNS_0_1,
 * and MBE_FRAME_MARK_LATE_ONSET whatever the frame before) in the rate-33
 * rendering; P25's (MBE_FRAME_MARK_ONSET_COLUMNS) in P25's. In D-STAR's (and
 * every other) rendering the marks 5, 6 and 7 play as none, and a late onset
 * acts only after a frame with no voiced band, as before these marks. (P25's
 * rendering replaces the mark of a voiced model:
 * test_p25_rendering_marks_the_model_it_renders().) A voiced frame (127.5 Hz, 23 of 28 harmonics voiced,
 * which covers the 500 Hz column 5) after an unvoiced and after a voiced one.
 */
static void
test_marks_act_in_their_own_rendering(void) {
    static const enum mbe_frame_mark marks[] = {MBE_FRAME_MARK_ONSET_COLUMNS, MBE_FRAME_MARK_MID_COLUMN_0,
                                                MBE_FRAME_MARK_MID_COLUMNS_0_1, MBE_FRAME_MARK_LATE_ONSET};

    static const enum mbe_voiced_synthesis renderings[] = {MBE_VOICED_SYNTHESIS_DVSI,
                                                           MBE_VOICED_SYNTHESIS_DVSI_AMBE2450};

    for (int voiced_before = 0; voiced_before < 2; ++voiced_before) {
        mbe_parms last, model;
        float out[160];
        fresh_stream(&last);
        set_model(&model, 150.0, 24, voiced_before ? 20 : 0, 16.0f);
        next_frame(out, &last, &model);
        set_model(&model, 127.5, 28, 23, 16.0f);
        shape_model(&model, 3);
        for (size_t r = 0; r < sizeof(renderings) / sizeof(renderings[0]); ++r) {
            const enum mbe_voiced_synthesis rendering = renderings[r];
            float plain[160], marked[160];
            mbe_parms plain_state, marked_state;
            render_marked(plain, &plain_state, &last, &model, MBE_FRAME_MARK_NONE, rendering);
            for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
                const enum mbe_frame_mark mark = marks[m];
                render_marked(marked, &marked_state, &last, &model, mark, rendering);
                int acts;
                if (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) {
                    acts = mark != MBE_FRAME_MARK_ONSET_COLUMNS;
                } else {
                    acts = mark == MBE_FRAME_MARK_LATE_ONSET && !voiced_before;
                }
                assert(same_rendering(marked, &marked_state, plain, &plain_state) == !acts);
                assert(all_finite(marked, 160));
            }
        }
    }
}

/*
 * P25's rendering marks the model it renders: a model with a voiced band that
 * is not a tone model carries MBE_FRAME_MARK_ONSET_COLUMNS whatever mark it
 * came with, so a caller's mark on it is overwritten and has no effect, the
 * late onset (80) and the held previous model (3) included, which before these
 * marks, as at the branch point, delayed the voiced bands and held the model
 * before: every one renders as the onset mark (after an unvoiced frame, where
 * the columns act, and after a voiced one), and an onset is read in the
 * columns unmarked too, unlike D-STAR's rendering. A tone model keeps
 * MBE_FRAME_MARK_TONE (its own model at mid-frame). A model with no voiced
 * band keeps any other mark (PITCH_CARRIES, PITCH_ABSENT and PREVIOUS_HELD
 * still act there) and loses the onset mark.
 */
static void
test_p25_rendering_marks_the_model_it_renders(void) {
    static const enum mbe_frame_mark marks[] = {MBE_FRAME_MARK_NONE,         MBE_FRAME_MARK_PITCH_CARRIES,
                                                MBE_FRAME_MARK_PITCH_ABSENT, MBE_FRAME_MARK_PREVIOUS_HELD,
                                                MBE_FRAME_MARK_MID_COLUMN_0, MBE_FRAME_MARK_MID_COLUMNS_0_1,
                                                MBE_FRAME_MARK_LATE_ONSET};
    const enum mbe_voiced_synthesis p25 = MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    for (int voiced_before = 0; voiced_before < 2; ++voiced_before) {
        mbe_parms last, model, unvoiced;
        float out[160], onset[160], marked[160], other[160];
        mbe_parms onset_state, marked_state, other_state;
        fresh_stream(&last);
        set_model(&model, 150.0, 24, voiced_before ? 20 : 0, 16.0f);
        next_frame(out, &last, &model);
        set_model(&model, 127.5, 28, 23, 16.0f);
        shape_model(&model, 3);
        render_marked(onset, &onset_state, &last, &model, MBE_FRAME_MARK_ONSET_COLUMNS, p25);
        assert(mbe_frame_mark(&onset_state) == MBE_FRAME_MARK_ONSET_COLUMNS);
        render_marked(other, &other_state, &last, &model, MBE_FRAME_MARK_NONE, MBE_VOICED_SYNTHESIS_DVSI);
        assert(same_rendering(other, &other_state, onset, &onset_state) == voiced_before);
        for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
            render_marked(marked, &marked_state, &last, &model, marks[m], p25);
            assert(same_rendering(marked, &marked_state, onset, &onset_state));
            assert(mbe_frame_mark(&marked_state) == MBE_FRAME_MARK_ONSET_COLUMNS);
        }
        render_marked(marked, &marked_state, &last, &model, MBE_FRAME_MARK_TONE, p25);
        assert(mbe_frame_mark(&marked_state) == MBE_FRAME_MARK_TONE);
        assert(!same_rendering(marked, &marked_state, onset, &onset_state));

        set_model(&unvoiced, 110.0, 30, 0, 8.0f);
        for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
            render_marked(marked, &marked_state, &last, &unvoiced, marks[m], p25);
            assert(mbe_frame_mark(&marked_state) == marks[m]);
        }
        render_marked(marked, &marked_state, &last, &unvoiced, MBE_FRAME_MARK_ONSET_COLUMNS, p25);
        assert(mbe_frame_mark(&marked_state) == MBE_FRAME_MARK_NONE);
    }
}

/* Three voiced P25 frames (b0 100, voiced throughout); with `onset`, then a
 * frame with no voiced band and the voiced frames' history put back, as a
 * caller's state can hold it: the history voiced, the last frame played
 * unvoiced. */
static void
p25_foreign_mark_stream(int onset, mbe_parms* cur, mbe_parms* prev, mbe_parms* enhanced) {
    char voiced[88], unvoiced[88] = {0};
    float out[160];
    memset(voiced, 1, sizeof(voiced));
    set_imbe_b0(voiced, 100);
    set_imbe_b0(unvoiced, 150);
    mbe_initMbeParms(cur, prev, enhanced);
    for (int f = 0; f < 3; ++f) {
        const int rc = mbe_processImbe4400Dataf(out, NULL, voiced, cur, prev, enhanced);
        assert(rc >= 0);
        (void)rc;
    }
    if (onset) {
        const mbe_parms history = *prev;
        const int rc = mbe_processImbe4400Dataf(out, NULL, unvoiced, cur, prev, enhanced);
        assert(rc >= 0);
        (void)rc;
        *prev = history;
    }
}

/*
 * A caller's state can hold any frame mark, and P25 plays rate 33's
 * (MBE_FRAME_MARK_MID_COLUMN_0, MBE_FRAME_MARK_MID_COLUMNS_0_1) as none: the
 * history, the last frame played or both so marked, then a frame the decoder
 * rejects (b0 255, which repeats the history over the last frame played) and a
 * voice frame play as the same stream with no mark there; after voiced frames,
 * and with the last frame played unvoiced and the history voiced (an onset).
 */
static void
test_p25_foreign_marks_in_state(void) {
    static const enum mbe_frame_mark marks[] = {MBE_FRAME_MARK_MID_COLUMN_0, MBE_FRAME_MARK_MID_COLUMNS_0_1};

    char rejected[88] = {0}, voiced[88];
    set_imbe_b0(rejected, 255);
    memset(voiced, 1, sizeof(voiced));
    set_imbe_b0(voiced, 104);
    for (int onset = 0; onset < 2; ++onset) {
        for (int target = 0; target < 3; ++target) {
            for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
                mbe_parms cur[2], prev[2], enhanced[2];
                float out[2][320];
                p25_foreign_mark_stream(onset, &cur[0], &prev[0], &enhanced[0]);
                cur[1] = cur[0];
                prev[1] = prev[0];
                enhanced[1] = enhanced[0];
                for (int s = 0; s < 2; ++s) {
                    const enum mbe_frame_mark mark = (s == 0) ? MBE_FRAME_MARK_NONE : marks[m];
                    if (target != 1) {
                        mbe_set_frame_mark(&prev[s], mark);
                    }
                    if (target != 0) {
                        mbe_set_frame_mark(&enhanced[s], mark);
                    }
                    mbe_process_result result;
                    mbe_initProcessResult(&result);
                    const int rc = mbe_processImbe4400Dataf(out[s], &result, rejected, &cur[s], &prev[s], &enhanced[s]);
                    assert(rc >= 0 && (result.flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
                    const int rc2 =
                        mbe_processImbe4400Dataf(out[s] + 160, NULL, voiced, &cur[s], &prev[s], &enhanced[s]);
                    assert(rc2 >= 0);
                    (void)rc;
                    (void)rc2;
                }
                assert(same_bits(out[0], out[1], 320));
            }
        }
    }
}

/* P25 parameter bits of pitch code b0 with no voiced band, log2 amplitude
 * `loud` at harmonics first..last and 2 elsewhere, quantized against the
 * history `prev` (the IMBE encoder's quantizer). */
static void
p25_unvoiced_frame(char imbe_d[88], int b0, int first, int last, const mbe_parms* prev) {
    const unsigned char bands[12] = {0};
    float logm[57];
    for (int l = 0; l < 57; ++l) {
        logm[l] = (l >= first && l <= last) ? 11.0f : 2.0f;
    }
    mbe_imbe4400_quantize_amplitudes(b0, bands, logm, prev, 0, imbe_d);
}

/* A P25 frame's result with errors enough for adaptive smoothing to voice its
 * loud harmonics (more than four, one in the Hamming-coded words, so that the
 * voicing threshold is TIA-102.BABA eq 112's second, the smaller), and too few
 * for the parameter-bit API to repeat the frame for them. */
static void
smoothing_result(mbe_process_result* result) {
    mbe_initProcessResult(result);
    result->flags = MBE_PROCESS_FLAG_C4_VALID;
    result->c4_errors = 1;
    result->total_errors = 5;
}

/*
 * A P25 frame that decodes with no voiced band after frames with none, where
 * the adaptive smoothing that channel errors bring voices some of its
 * harmonics: the model the synthesis renders, so voiced, is read in the onset
 * columns, as a decoded onset is. One loud harmonic voiced alone (12, at
 * 1.2 kHz at b0 120; 14, at 1.2 kHz, at b0 150) covers no column, so the pulse train starts
 * at mid-frame; loud harmonics 5 to 12 (500 to 1200 Hz at b0 120) cover the
 * column from 505 Hz, and the train starts with the frame. The decoded history
 * keeps the decoded model and its mark (b0 120 lacks a pitch, b0 150 is one of
 * the codes kept for frames with no voiced band); the model rendered carries
 * MBE_FRAME_MARK_ONSET_COLUMNS. Through a repeat of such a history (a rejected
 * frame with the same errors) likewise, and a muted frame advances to the pulse
 * train a heard one leaves.
 */
static void
test_p25_smoothed_onset(void) {
    static const struct {
        int b0, first, last, early;
    } cases[] = {{120, 12, 12, 0}, {150, 14, 14, 0}, {120, 5, 12, 1}};

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        for (int repeat = 0; repeat < 2; ++repeat) {
            mbe_parms cur, prev, enhanced;
            char imbe_d[88];
            float out[160];
            mbe_initMbeParms(&cur, &prev, &enhanced);
            for (int f = 0; f < 8; ++f) {
                p25_unvoiced_frame(imbe_d, cases[i].b0, cases[i].first, cases[i].last, &prev);
                const int rc = mbe_processImbe4400Dataf(out, NULL, imbe_d, &cur, &prev, &enhanced);
                assert(rc >= 0);
                (void)rc;
            }
            assert(!reference_voiced(&enhanced));
            const enum mbe_frame_mark decoded =
                (cases[i].b0 >= 128) ? MBE_FRAME_MARK_NONE : MBE_FRAME_MARK_PITCH_ABSENT;
            assert(mbe_frame_mark(&prev) == decoded);
            if (repeat) {
                memset(imbe_d, 1, sizeof(imbe_d)); /* b0 255: rejected */
            } else {
                p25_unvoiced_frame(imbe_d, cases[i].b0, cases[i].first, cases[i].last, &prev);
            }
            mbe_parms muted_cur = cur, muted_prev = prev, muted_enhanced = enhanced;
            mbe_process_result result;
            smoothing_result(&result);
            const int rc = mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &enhanced);
            assert(rc >= 0 && ((result.flags & MBE_PROCESS_FLAG_REPEAT) != 0u) == repeat);
            assert((result.flags & MBE_PROCESS_FLAG_MUTE) == 0u);
            (void)rc;
            /* The decoded history: no voiced band, its own mark. */
            assert(!reference_voiced(&prev) && mbe_frame_mark(&prev) == decoded);
            /* The model rendered: the loud harmonics voiced, and marked. */
            for (int l = 1; l <= enhanced.L; ++l) {
                assert(enhanced.Vl[l] == (l >= cases[i].first && l <= cases[i].last));
            }
            assert(mbe_frame_mark(&enhanced) == MBE_FRAME_MARK_ONSET_COLUMNS);
            const double hz = (double)enhanced.w0 * 8000.0 / (2.0 * M_PI);
            assert(fabs(cycles_apart(pulse_phase(&enhanced), cycles_of(hz, cases[i].early ? 160.0 : 80.0))) < 1e-4);
            assert(all_finite(out, 160));

            /* The same frame muted by the error rate: the train the heard one leaves. */
            muted_prev.errorRate = 1.0f;
            smoothing_result(&result);
            const int rc2 = mbe_processImbe4400Dataf(out, &result, imbe_d, &muted_cur, &muted_prev, &muted_enhanced);
            assert(rc2 >= 0 && (result.flags & MBE_PROCESS_FLAG_MUTE) != 0u);
            (void)rc2;
            assert(mbe_frame_mark(&muted_enhanced) == MBE_FRAME_MARK_ONSET_COLUMNS);
            assert(same_bits(muted_enhanced.noiseOverlap, enhanced.noiseOverlap, 96));
            assert(same_bits(&muted_enhanced.PSIl[1], &enhanced.PSIl[1], 56));
        }
    }
}

int
main(void) {
    test_frame_marks();
    test_steady_model_is_its_harmonics();
    test_onset_starts_the_train();
    test_late_onset();
    test_pitch_step_keeps_the_train();
    test_waveform_through_changes(MBE_VOICED_SYNTHESIS_DVSI, MID_PITCH_INEXACT);
    test_waveform_through_changes(MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MID_PITCH_INEXACT);
    test_waveform_through_changes(MBE_VOICED_SYNTHESIS_DVSI_IMBE, MID_PITCH_INEXACT);
    test_octave_steps();
    test_stored_phase_range();
    test_unvoiced_pitch_carries();
    test_voicing_persists_without_a_pitch();
    test_mid_model();
    test_phases_follow_decoded_amplitudes();
    test_mid_model_phases();
    test_untrusted_log2_amplitudes();
    test_previous_model_held();
    test_subframes_crossfade_linearly();
    test_muted_frame_advances(MBE_VOICED_SYNTHESIS_DVSI);
    test_muted_frame_advances(MBE_VOICED_SYNTHESIS_DVSI_AMBE2450);
    test_muted_frame_advances(MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    test_advance_follows_the_marks(MBE_VOICED_SYNTHESIS_DVSI);
    test_advance_follows_the_marks(MBE_VOICED_SYNTHESIS_DVSI_AMBE2450);
    test_advance_follows_the_marks(MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    test_ambe2450_mid_pitch(MID_PITCH_FAST_MATH, MID_PITCH_INEXACT);
    test_ambe2450_mid_pitch_codes(MID_PITCH_FAST_MATH, MID_PITCH_INEXACT);
    test_from_silence();
    test_output_limit_option();
    test_untrusted_unvoiced_overlap();
    test_noise_overlap_guard();
    test_pulse_state_is_validated();
    test_invalid_models();
    test_subnormal_pitch_is_no_pitch();
    test_decoders_mark_unvoiced_frames();
    test_p25_onset_columns();
    test_p25_onset_columns_only_at_onsets();
    test_p25_onset_edges();
    test_p25_muted_onset_advances();
    test_untrusted_mid_voicing_marks();
    test_p25_voicing_persists();
    test_p25_reset_drops_the_mark();
    test_p25_foreign_marks_in_state();
    test_marks_act_in_their_own_rendering();
    test_p25_rendering_marks_the_model_it_renders();
    test_p25_smoothed_onset();
    return 0;
}
