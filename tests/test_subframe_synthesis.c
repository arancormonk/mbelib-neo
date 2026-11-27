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
#include "ambe3600x2400_internal.h"
#include "ambe3600x2450_internal.h"
#include "ambe_common.h"
#include "imbe4400_internal.h"
#include "mbe_adaptive.h"
#include "mbe_low_band.h"
#include "mbe_math.h"
#include "mbe_repeat.h"
#include "mbe_subframe.h"
#include "mbe_tone.h"
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

/* Whether this test can reproduce the library's float rendering of a pulse op
 * for op and so compare it bit for bit (test_r33_pulse_times_by_bits()): the
 * library built without fast-math, float expressions evaluated in float
 * (FLT_EVAL_METHOD 0), or, in an unoptimized build, wider but stored at each
 * statement as the library stores them, and no multiply-add contracted where
 * the target fuses them (CMake defines MBELIB_TEST_FP_CONTRACT_OFF where it
 * builds library and tests with -ffp-contract=off, the Debug builds of GCC and
 * Clang). Elsewhere the test compares within a tolerance, which cannot tell
 * the pulse times formed in double from the same times formed in float. */
#if MID_PITCH_FAST_MATH || !defined(FLT_EVAL_METHOD)
#define PULSE_REPLICA_EXACT 0
#elif (FLT_EVAL_METHOD == 0 || (FLT_EVAL_METHOD == 2 && !defined(__OPTIMIZE__)))                                       \
    && (!defined(__FP_FAST_FMA) || defined(MBELIB_TEST_FP_CONTRACT_OFF))
#define PULSE_REPLICA_EXACT 1
#else
#define PULSE_REPLICA_EXACT 0
#endif

/* Whether float expressions may be evaluated wider than float (x87 excess
 * precision, or an evaluation method left open), so that a float quotient the
 * library forms need not be rounded to float before it is used. */
#if !defined(FLT_EVAL_METHOD) || (FLT_EVAL_METHOD != 0)
#define HOLD_QUOTIENT_WIDE 1
#else
#define HOLD_QUOTIENT_WIDE 0
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

/* A model of L harmonics at `hz`, the first `voiced` of them voiced. Sets
 * the model's fields only (element 0 of Vl, Ml and log2Ml to 0, as no model
 * uses it): a model handed to the library whole starts from an initialized
 * state (blank_model()). */
static void
set_model(mbe_parms* mp, double hz, int L, int voiced, float amplitude) {
    mp->w0 = w0_of(hz);
    mp->L = L;
    mp->Vl[0] = 0;
    mp->Ml[0] = 0.0f;
    mp->log2Ml[0] = 0.0f;
    for (int l = 1; l <= 56; ++l) {
        mp->Vl[l] = (l <= voiced) ? 1 : 0;
        set_amplitude(mp, l, (l <= L) ? amplitude : 0.0f);
    }
    mbe_set_frame_mark(mp, MBE_FRAME_MARK_NONE);
}

/* A model in an initialized state (mbe_initMbeParms()), so that every field
 * the library may read has a value, then set as set_model() sets it. */
static void
blank_model(mbe_parms* mp, double hz, int L, int voiced, float amplitude) {
    mbe_parms prev, enhanced;
    mbe_initMbeParms(mp, &prev, &enhanced);
    set_model(mp, hz, L, voiced, amplitude);
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

/* next_frame_as() for a model carrying a word in PHIl[0] in place of a mark
 * (a rate-33 pulsed-row word): the word goes with the model as stored. */
static void
next_frame_with_word(float out[160], mbe_parms* last, const mbe_parms* model, int from_silence,
                     enum mbe_voiced_synthesis rendering) {
    mbe_parms cur = *last;
    cur.w0 = model->w0;
    cur.L = model->L;
    memcpy(cur.Vl, model->Vl, sizeof(cur.Vl));
    copy_amplitudes(&cur, model);
    cur.PHIl[0] = model->PHIl[0];
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

/* The low band every subframe rendering plays, from its statement
 * (mbe_low_band.h, mbe_adaptive.h; docs/dvsi-parity.md, "DVSI's low band"), in
 * double. D-STAR's, P25's and rate 33's renderings: a third-order Butterworth
 * high-pass at 75 Hz, H(s) = s^3 / ((s + 1)(s^2 + s + 1)) at s = j hz / 75, and
 * the adaptive cut, (f / corner)^2 on a band under the corner; ProVoice's
 * (and every other value's): a second-order Butterworth high-pass at 100 Hz,
 * H(s) = s^2 / (s^2 + sqrt2 s + 1) at s = j hz / 100, and no cut. */
struct reference_band {
    int dvsi;      /* DVSI's filter and the cut */
    double corner; /* the cut's corner, Hz; 0 for none */
};

static int
reference_dvsi_band(enum mbe_voiced_synthesis rendering) {
    return rendering == MBE_VOICED_SYNTHESIS_DVSI || rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE
           || rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
}

/* A rendering's low band with no cut (a stream whose corner lies under every
 * harmonic the test reads). */
static struct reference_band
reference_band_of(enum mbe_voiced_synthesis rendering) {
    const struct reference_band band = {reference_dvsi_band(rendering), 0.0};
    return band;
}

static double
reference_second_order_gain(double hz) {
    if (!(hz > 0.0)) {
        return 0.0;
    }
    const double x = hz / 100.0;
    return (x * x) / sqrt((x * x * x * x) + 1.0);
}

static double
reference_second_order_phase(double hz) {
    if (!(hz > 0.0)) {
        return 0.0;
    }
    const double x = hz / 100.0;
    return M_PI - atan2(sqrt(2.0) * x, 1.0 - (x * x));
}

/* The third-order filter from its pole-zero form, |H| and arg H as the
 * product of its three factors (not the library's closed form). */
static double
reference_third_order_gain(double hz) {
    if (!(hz > 0.0)) {
        return 0.0;
    }
    const double x = hz / 75.0;
    const double first = sqrt(1.0 + (x * x));                                  /* |j x + 1| */
    const double second = sqrt(((1.0 - (x * x)) * (1.0 - (x * x))) + (x * x)); /* |(j x)^2 + j x + 1| */
    return (x * x * x) / (first * second);
}

static double
reference_third_order_phase(double hz) {
    if (!(hz > 0.0)) {
        return 0.0;
    }
    const double x = hz / 75.0;
    return (1.5 * M_PI) - atan2(x, 1.0) - atan2(x, 1.0 - (x * x));
}

/* A harmonic's frequency as the library's low band compares it with the
 * corner (the arithmetic its rule states: w0 times 8000 / (2 pi) in float,
 * rounded to float, times l, rounded again). */
static float
reference_knee_hz(float w0, int l) {
    volatile float per_harmonic = w0 * (8000.0f / (2.0f * (float)M_PI));
    volatile float hz = per_harmonic * (float)l;
    return hz;
}

/* The cut's scale on harmonic l of pitch w0: (f / corner)^2 under the corner
 * (f compared as the library forms it), 1 at or above it or with no cut. */
static double
reference_cut(const struct reference_band* band, float w0, int l) {
    const double hz = (double)reference_knee_hz(w0, l);
    if (!band->dvsi || !(band->corner > 0.0) || !(hz < band->corner)) {
        return 1.0;
    }
    return (hz / band->corner) * (hz / band->corner);
}

static double
reference_band_gain(const struct reference_band* band, double hz) {
    return band->dvsi ? reference_third_order_gain(hz) : reference_second_order_gain(hz);
}

static double
reference_band_phase(const struct reference_band* band, double hz) {
    return band->dvsi ? reference_third_order_phase(hz) : reference_second_order_phase(hz);
}

/* The adaptive low cut's tracker, from its statement (mbe_low_band.h): 110 Hz
 * at a stream's start; one frame on, a frame whose model has every band 1..L
 * voiced, a positive pitch and no tone mark moves it toward its fundamental,
 * T + 0.043 (f0 - T), in double and unrounded here; the corner a frame applies
 * is min(T, 200) of the tracker the frame before left, its own model and its
 * mid-frame model at that corner, the last frame's own model at the corner
 * that frame applied. */
struct reference_tracker {
    double track;  /* the tracker; 0 before the stream's first frame (110 then) */
    double corner; /* the corner the last frame applied; 0 before the first frame */
};

static double
reference_track(const struct reference_tracker* tracker) {
    return (tracker->track > 0.0) ? tracker->track : 110.0;
}

static double
reference_corner_now(const struct reference_tracker* tracker) {
    return fmin(reference_track(tracker), 200.0);
}

/* Whether a frame's model, as the synthesis renders it, moves the tracker. */
static int
reference_moves_tracker(const mbe_parms* mp) {
    if (!(mp->w0 > 0.0f) || mbe_frame_mark(mp) == MBE_FRAME_MARK_TONE || mp->L < 1) {
        return 0;
    }
    for (int l = 1; l <= mp->L; ++l) {
        if (mp->Vl[l] != 1) {
            return 0;
        }
    }
    return 1;
}

/* The three models' bands of a frame in `rendering` (the last frame's own, the
 * mid-frame one, this frame's own), then the tracker one frame on. */
static void
reference_tracker_frame(struct reference_tracker* tracker, const mbe_parms* cur, enum mbe_voiced_synthesis rendering,
                        struct reference_band band[3]) {
    const double now = reference_corner_now(tracker);
    const int dvsi = reference_dvsi_band(rendering);
    band[0].dvsi = band[1].dvsi = band[2].dvsi = dvsi;
    band[0].corner = dvsi ? ((tracker->corner > 0.0) ? tracker->corner : now) : 0.0;
    band[1].corner = band[2].corner = dvsi ? now : 0.0;
    if (!dvsi) {
        return;
    }
    double track = reference_track(tracker);
    if (reference_moves_tracker(cur)) {
        const double f0 = (double)cur->w0 * 8000.0 / (2.0 * M_PI);
        track = fmin(fmax(track + (0.043 * (f0 - track)), 50.0), 500.0);
    }
    tracker->track = track;
    tracker->corner = now;
}

/* Whether a frame stored the tracker and corner the reference holds: within a
 * hundred-thousandth (the library rounds the tracker to float once a frame,
 * the reference does not: the two drift apart by under a millionth), and 0 in
 * both slots in a rendering without the low cut. */
static int
stored_tracker_as(const mbe_parms* mp, const struct reference_tracker* tracker, enum mbe_voiced_synthesis rendering) {
    const double track = (double)mp->noiseOverlap[MBE_LOW_CUT_TRACK_SLOT];
    const double corner = (double)mp->noiseOverlap[MBE_LOW_CUT_CORNER_SLOT];
    if (!reference_dvsi_band(rendering)) {
        uint32_t a, b;
        memcpy(&a, &mp->noiseOverlap[MBE_LOW_CUT_TRACK_SLOT], sizeof(a));
        memcpy(&b, &mp->noiseOverlap[MBE_LOW_CUT_CORNER_SLOT], sizeof(b));
        return a == 0u && b == 0u;
    }
    return fabs(track - tracker->track) <= 1.0e-5 * tracker->track
           && fabs(corner - tracker->corner) <= 1.0e-5 * tracker->corner;
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
    /* Thirty periods from the third frame on. The low cut's corner stays under
     * the fundamental (the tracker rises toward 187.5 Hz from 110 and never
     * reaches it), so no harmonic is cut: D-STAR's filter alone. */
    const struct reference_band band = reference_band_of(MBE_VOICED_SYNTHESIS_DVSI);
    assert(last.noiseOverlap[MBE_LOW_CUT_TRACK_SLOT] < 187.0f);
    for (int l = 1; l <= model.L; ++l) {
        const double hz = 187.5 * l;
        const double want = 2.0 * (double)model.Ml[l] * reference_band_gain(&band, hz);
        const double got = component_amplitude(out + 320, 1280, hz);
        assert(fabs(got / want - 1.0) < 2e-3);
        const double phase = component_phase(&out[1344], 256, hz, 256); /* the last 256 samples */
        const double stored = (double)last.PHIl[l] + reference_band_phase(&band, hz);
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

/* Whether harmonic l of pitch w0 is rendered: below Nyquist in the library's
 * arithmetic, the product l w0 formed and compared in float (a harmonic just
 * under pi in double can round to pi in float and is then not rendered).
 * Stored through a volatile float so that this test rounds it under x87 excess
 * precision too; the library's comparison may keep the wider value there, so
 * the test that puts a harmonic on that edge runs in IEEE builds only. */
static int
below_nyquist(float w0, int l) {
    volatile float w = w0 * (float)l;
    return w < (float)M_PI;
}

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

/* What a model plays t samples after a pulse, its harmonics at `phase`, under
 * its low band (`band`: the filter, and the cut at the model's corner). */
static double
reference_response(const struct mbe_subframe_model* model, const struct reference_band* band, const double phase[57],
                   double t) {
    double sum = 0.0;
    if (!model->any_voiced) {
        return 0.0;
    }
    for (int l = 1; l <= model->L; ++l) {
        const double w = (double)model->w0 * l;
        const double hz = w * 8000.0 / (2.0 * M_PI);
        if (model->voiced[l] != 1 || !below_nyquist(model->w0, l)) {
            continue;
        }
        sum += 2.0 * (double)model->magnitude[l] * reference_band_gain(band, hz) * reference_cut(band, model->w0, l)
               * cos((w * (t - 16.0)) + phase[l] + reference_band_phase(band, hz));
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
                 const struct reference_band* a_band, const double a_phase[57], const struct mbe_subframe_model* b,
                 const struct reference_band* b_band, const double b_phase[57]) {
    for (int n = 0; n < 80; ++n) {
        const double fade = (double)n / 80.0;
        double sum = 0.0;
        for (int i = 0; i < train->count; ++i) {
            const double t = start + (double)n - train->time[i];
            sum += ((1.0 - fade) * reference_response(a, a_band, a_phase, t))
                   + (fade * reference_response(b, b_band, b_phase, t));
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

static double away_from_whole(double x);

/* The IMBE renderings, P25's and ProVoice's: P25's 500 Hz columns and marks. */
static int
reference_imbe(enum mbe_voiced_synthesis rendering) {
    return rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE || rendering == MBE_VOICED_SYNTHESIS_DVSI_PROVOICE;
}

/*
 * The pulsed state of D-STAR's and P25's renderings, from the rule's statement
 * (docs/dvsi-parity.md, "D-STAR and P25: the pulsed state"), apart from the
 * library's own helpers: the 56 patterns of US 7,970,606's Table 2 typed here
 * as the table prints them (code b0 - 128 for P25, b0 - 64 for D-STAR; eight
 * channels of 500 Hz for the frame's own subframe, then eight for its
 * mid-frame subframe, lowest first, P pulsed), and the pitch of each code as
 * the codec's public decoder gives it.
 */
static const char* const k_pulsed_rows[56] = {
    "UUUUUUUP UUUUUUUU", "UUUUUUPU UUUUUUUU", "UUUUUUPP UUUUUUUU", "UUUUUPUU UUUUUUUU", "UUUUUPPP UUUUUUUU",
    "UUUUPPPP UUUUUUUU", "UUUPUUUU UUUUUUUU", "UUUPPUUU UUUUUUUU", "UUUUUUUU PUUUUUUU", "UUUUUUUU UUUPPUUU",
    "UUUUUUUU UUUPUUUU", "UUUUUUUU UUPPPPPU", "UUUUUUUU UUPUUUUU", "UUUUUUUU UUPPUUUU", "UUUUUUUU UUPPPUUU",
    "UUUUUUUU UUUPPPUU", "UUUUUUUU UUUUUUUP", "UUUUUUUU UUUUUPPP", "UUUUUUUU UUUUUUPU", "UUUUUUUU UUUUUPUU",
    "UUUUUUUU UUUUUUPP", "UUUUUUUU UUUUPPPP", "UUUUUUUU UUUUPUUU", "UUUUUUUU UUUUPPUU", "PUUUUUUU UUUUUUUU",
    "UUPPPUUU UUUUUUUU", "PPPUUUUU UUUUUUUU", "UUUPPPUU UUUUUUUU", "PPUUUUUU UUUUUUUU", "PPPPUUUU UUUUUUUU",
    "UUUPPPPP UUUUUUUU", "UUPPUUUU UUUUUUUU", "PPPPPPPP UUUUUUUU", "UUUUPUUU UUUUUUUU", "UUPPPPPP UUUUUUUU",
    "PPPPPUUU UUUUUUUU", "UUPUUUUU UUUUUUUU", "PPPPPPUU UUUUUUUU", "UPPPPPPP UUUUUUUU", "PPPPPPPU UUUUUUUU",
    "UUUUUUUU PPPPPPPP", "UUUUUUUU PPPPPUUU", "UUUUUUUU PPPUUUUU", "UUUUUUUU UPPPPPPP", "UUUUUUUU PPUUUUUU",
    "UUUUUUUU UUUPPPPP", "UUUUUUUU UUPPPPPP", "UUUUUUUU PPPPPPUU", "PPPPPPPP PPPPPPPP", "PPPPPPPP UUPPPPPP",
    "PPPUUUUU PPPUUUUU", "UUUPPPPP UUUPPPPP", "UUPPPPPP UUPPPPPP", "PPPPPPPP PPPUUUUU", "PPUUUUUU PPUUUUUU",
    "UUPUUUUU UUPUUUUU",
};

/* Whether channel k (0..7) of a code's own (mid 0) or mid-frame (mid 1) pattern is pulsed. */
static int
reference_pattern(int code, int mid, int k) {
    return k_pulsed_rows[code][(mid ? 9 : 0) + k] == 'P';
}

static int
reference_pattern_pulses(int code, int mid) {
    int any = 0;
    for (int k = 0; k < 8 && code >= 0; ++k) {
        any |= reference_pattern(code, mid, k);
    }
    return any;
}

/* The pitch the codec's public decoder gives pitch code b0 (D-STAR's when
 * `dstar`, else P25's), for a frame with no voiced band. */
static float
decoded_pitch(int dstar, int b0) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    if (dstar) {
        char ambe_d[49] = {0};
        for (int i = 0; i < 6; ++i) {
            ambe_d[i] = (char)((b0 >> (6 - i)) & 1);
        }
        ambe_d[48] = (char)(b0 & 1);
        const int status = mbe_decodeAmbe2400Parms(ambe_d, &cur, &prev);
        assert(status == 0);
        (void)status;
    } else {
        char imbe_d[88] = {0};
        set_imbe_b0(imbe_d, b0);
        const int status = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
        assert(status == 0);
        (void)status;
    }
    return cur.w0;
}

/* The pitch of reserved code `code` (0..55) in the codec of a rendering:
 * P25's for the IMBE renderings (ProVoice's too, which has no code: so a
 * ProVoice stream meets P25's reserved pitches), D-STAR's for the others. */
static float
reserved_pitch(enum mbe_voiced_synthesis rendering, int code) {
    static float pitches[2][56];
    static int filled = 0;
    if (!filled) {
        for (int c = 0; c < 56; ++c) {
            pitches[0][c] = decoded_pitch(0, 128 + c);
            pitches[1][c] = decoded_pitch(1, 64 + c);
        }
        filled = 1;
    }
    return pitches[reference_imbe(rendering) ? 0 : 1][code];
}

/* The reserved code of a model, or -1: in D-STAR's and P25's renderings, a
 * model with no voiced band at a reserved code's decoded pitch, bit for bit;
 * in rate 33's and ProVoice's none. */
static int
reference_code(const mbe_parms* mp, enum mbe_voiced_synthesis rendering) {
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 || rendering == MBE_VOICED_SYNTHESIS_DVSI_PROVOICE) {
        return -1;
    }
    for (int l = 1; l <= mp->L; ++l) {
        if (mp->Vl[l] == 1) {
            return -1;
        }
    }
    for (int c = 0; c < 56; ++c) {
        const float pitch = reserved_pitch(rendering, c);
        if (same_bits(&pitch, &mp->w0, 1)) {
            return c;
        }
    }
    return -1;
}

/*
 * What a frame of the pulsed state leaves the next, as this test keeps it:
 * whether its mid-frame model emitted, whether its own code was read, and
 * P25's validity counter (US 7,970,606, FIG. 8): 0 for a fresh stream; a
 * model with a voiced band leaves it, one with no voiced band at the pitch of
 * the default code b0 25 or of a reserved code raises it by one to 3 at most,
 * one with no voiced band at the pitch of any other P25 code (b0 0..207, as
 * the public decoder gives it) lowers it by one to -10 at least, and one at a
 * pitch no code decodes to leaves it; P25 reads a frame's code only while the
 * counter the frame before left is 0 or more. D-STAR reads every code and
 * keeps no counter.
 */
struct reference_stream {
    int emitted, read, count;
};

/* The pitch P25's public decoder gives its default code, b0 25. */
static float
p25_default_pitch(void) {
    static float pitch = 0.0f;
    if (!(pitch > 0.0f)) {
        pitch = decoded_pitch(0, 25);
    }
    return pitch;
}

/* The codes the rendering reads: the last frame's (codes[0]) where that frame
 * read it, this frame's (codes[1]) where the counter allows; -1 for none. */
static void
reference_read_codes(const mbe_parms* last, const mbe_parms* cur, enum mbe_voiced_synthesis rendering,
                     const struct reference_stream* stream, int codes[2]) {
    const int p25 = rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    codes[0] = (!p25 || stream->read) ? reference_code(last, rendering) : -1;
    codes[1] = (!p25 || stream->count >= 0) ? reference_code(cur, rendering) : -1;
}

/* The P25 pitch code (0..207) whose pitch, as the public decoder gives it,
 * w0 is bit for bit, or -1. */
static int
reference_p25_pitch_code(float w0) {
    static float pitches[208];
    static int filled = 0;
    if (!filled) {
        for (int b0 = 0; b0 < 208; ++b0) {
            pitches[b0] = decoded_pitch(0, b0);
        }
        filled = 1;
    }
    for (int b0 = 0; b0 < 208; ++b0) {
        if (same_bits(&pitches[b0], &w0, 1)) {
            return b0;
        }
    }
    return -1;
}

/* The counter a P25 frame leaves after `count`. */
static int
reference_next_count(const mbe_parms* cur, int count) {
    if (reference_voiced(cur)) {
        return count;
    }
    const int b0 = reference_p25_pitch_code(cur->w0);
    if (b0 < 0) {
        return count;
    }
    if (b0 == 25 || (b0 >= 128 && b0 <= 183)) {
        return (count < 3) ? count + 1 : 3;
    }
    return (count > -10) ? count - 1 : -10;
}

/* The pulsed bands (1) of a model under a code's pattern, from the rule: in a
 * model with no voiced band, an unvoiced band whose 500 Hz channel (the
 * eighth open above 3500 Hz), or the channel of 250 Hz higher, is pulsed. The
 * band's frequency is l w0 4000 / pi Hz in double, as the rule states it, the
 * float pitch widened: a harmonic just under a channel's edge (749.999964 Hz)
 * lies in the channel below, where the product in float would round it onto
 * the edge. */
static int
reference_pulse_bands(const struct mbe_subframe_model* model, int code, int mid, int pulsed[57]) {
    int any = 0;
    memset(pulsed, 0, 57 * sizeof(int));
    if (code < 0 || model->any_voiced) {
        return 0;
    }
    for (int l = 1; l <= model->L; ++l) {
        const double hz = (double)l * (double)model->w0 * (4000.0 / M_PI);
        if (model->voiced[l] != 0) {
            continue;
        }
        const int own = (int)fmin(floor(hz / 500.0), 7.0);
        const int reach = (int)fmin(floor((hz + 250.0) / 500.0), 7.0);
        if (reference_pattern(code, mid, own) || reference_pattern(code, mid, reach)) {
            pulsed[l] = 1;
            any = 1;
        }
    }
    return any;
}

/* A pulsed band's magnitude, in the library's arithmetic (pulse_bands()): the
 * unvoiced factor MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt(w0) in float, rounded to
 * float (mbe_round_to_float()), the model's magnitude divided by it in double,
 * the quotient rounded to float, and the result capped at 1e30f, the largest
 * magnitude a model takes. Under fast-math the library may evaluate this
 * otherwise; the tests that use it compare within a tolerance there. */
static double
reference_pulsed_magnitude(const struct mbe_subframe_model* model, int l) {
    const float factor = mbe_round_to_float(MBE_UNVOICED_MAGNITUDE_FACTOR / sqrtf(model->w0));
    const double magnitude = (double)model->magnitude[l] / (double)factor;
    const float capped = (magnitude < (double)1.0e30f) ? (float)magnitude : 1.0e30f;
    return (double)capped;
}

/* What a model's pulsed bands play t samples after its pulse: the voiced
 * response, each band at its magnitude with the unvoiced factor taken out
 * (reference_pulsed_magnitude()). */
static double
reference_pulsed_response(const struct mbe_subframe_model* model, const struct reference_band* band,
                          const int pulsed[57], const double phase[57], double t) {
    double sum = 0.0;
    for (int l = 1; l <= model->L; ++l) {
        const double w = (double)model->w0 * l;
        const double hz = w * 8000.0 / (2.0 * M_PI);
        if (!pulsed[l] || !below_nyquist(model->w0, l)) {
            continue;
        }
        const double magnitude = reference_pulsed_magnitude(model, l);
        sum += 2.0 * magnitude * reference_band_gain(band, hz) * reference_cut(band, model->w0, l)
               * cos((w * (t - 16.0)) + phase[l] + reference_band_phase(band, hz));
    }
    return reference_window(t, 2.0 * M_PI / (double)model->w0) * sum;
}

/* The one pulse model a emits, 24 samples into the interval it starts, played
 * through its pulsed bands going out and model b's coming in, within the
 * interval's 80 samples. */
static void
reference_pulsed_render(double* out, const struct mbe_subframe_model* a, const struct reference_band* a_band,
                        const int a_pulsed[57], const double a_phase[57], const struct mbe_subframe_model* b,
                        const struct reference_band* b_band, const int b_pulsed[57], const double b_phase[57]) {
    for (int n = 24; n < 80; ++n) {
        const double fade = (double)n / 80.0;
        const double t = (double)n - 24.0;
        out[n] += ((1.0 - fade) * reference_pulsed_response(a, a_band, a_pulsed, a_phase, t))
                  + (fade * reference_pulsed_response(b, b_band, b_pulsed, b_phase, t));
    }
}

/*
 * Rate 33's pulsed rows, from the rule's statement (docs/dvsi-parity.md, "Rate
 * 33: the pulsed rows"), apart from the library's helpers: the masks of V/UV
 * rows 18 to 31 typed here as the document prints them (the frame's own
 * subframe, then the mid-frame subframe, channel k from 500 k Hz, lowest first,
 * P pulsed), the word a frame carries in PHIl[0] (1024 + 2048 held + 128 (row
 * - 18) + b0) and in the state's slot 37 (the same, not held), formed here,
 * and the pulse times computed here. What the test takes from the library for
 * them: the band-grid model of a pitch code (its w0 and L) as the public
 * decoder, mbe_decodeAmbe2450Parms(), gives it.
 */
static const char* const k_row_masks[14] = {
    "PUUUUUUU UUUUUUUU", "PPUUUUUU UUUUUUUU", "PPPPPPPP UUUUUUUU", "UUPUUUUU UUUUUUUU", "UUUUUUUP UUUUUUUU",
    "UUUUUPPP UUUUUUUU", "PPPPPPPP PPPPPPPP", "UUPPPPPP UUPPPPPP", "UUUUUUUU PUUUUUUU", "UUUUUUUU UUPUUUUU",
    "UUUUUUUU PPPPPPPP", "UUUUUUUU PPUUUUUU", "UUUUUUUU UUUUUUUP", "UUUUUUUU UUUUUUPP",
};

#define TEST_ROW_SLOT      37
#define TEST_ROW_WORD_BASE 1024

/* The word of row (18..31), pitch code b0 (0..119) and held, as a whole number. */
static int
test_row_word(int row, int b0, int held) {
    return TEST_ROW_WORD_BASE + (2048 * held) + (128 * (row - 18)) + b0;
}

/* The fields of a word this test made (test_row_word()). */
static void
test_row_fields(int word, int* row, int* b0, int* held) {
    const int fields = word - TEST_ROW_WORD_BASE;
    *held = fields / 2048;
    *row = 18 + ((fields % 2048) / 128);
    *b0 = fields % 128;
}

/* A word's bits as a float. */
static uint32_t
test_row_word_bits(int word) {
    const float value = (float)word;
    uint32_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

/* The band-grid model of rate-33 pitch code b0, as the public decoder gives a
 * voice frame of V/UV row 18 at that code. */
static void
r33_grid(int b0, float* w0, int* L) {
    static float pitches[120];
    static int counts[120];
    static int filled = 0;
    if (!filled) {
        for (int code = 0; code < 120; ++code) {
            mbe_parms cur, prev, enhanced;
            char ambe_d[49];
            int b[9] = {code, 18, 20, 100, 30, 10, 5, 5, 3};
            mbe_ambe2450_pack(b, ambe_d);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            const int status = mbe_decodeAmbe2450Parms(ambe_d, &cur, &prev);
            assert(status == MBE_AMBE2450_FRAME_VOICE);
            (void)status;
            pitches[code] = cur.w0;
            counts[code] = cur.L;
        }
        filled = 1;
    }
    *w0 = pitches[b0];
    *L = counts[b0];
}

/* Whether channel k (0..7) of a row's own (mid 0) or mid-frame (mid 1) mask is pulsed. */
static int
reference_row_channel(int row, int mid, int k) {
    return k_row_masks[row - 18][(mid ? 9 : 0) + k] == 'P';
}

/* The pulsed bands (1) of a model under a row's mask: as for the codes
 * (reference_pulse_bands()), an unvoiced band of a model with no voiced band
 * whose channel, or the channel of 250 Hz more, is pulsed. */
static int
reference_row_bands(const struct mbe_subframe_model* model, int row, int mid, int pulsed[57]) {
    int any = 0;
    memset(pulsed, 0, 57 * sizeof(int));
    if (model->any_voiced) {
        return 0;
    }
    for (int l = 1; l <= model->L; ++l) {
        const double hz = (double)l * (double)model->w0 * (4000.0 / M_PI);
        if (model->voiced[l] != 0) {
            continue;
        }
        const int own = (int)fmin(floor(hz / 500.0), 7.0);
        const int reach = (int)fmin(floor((hz + 250.0) / 500.0), 7.0);
        if (reference_row_channel(row, mid, own) || reference_row_channel(row, mid, reach)) {
            pulsed[l] = 1;
            any = 1;
        }
    }
    return any;
}

/* The time of the pulse of model `index` (0 the last frame's own, at sample 0;
 * 1 the mid-frame model, at 80; 2 this frame's own, at 160), in samples from
 * the frame's start, as the rule states it: point + 1.016 b0 - 80.89 for rows
 * 18 to 23 and 26 to 31; for rows 24 and 25 point - 75.6 + 11.875 i, i b0 mod
 * 11 for the mid-frame model and b0 div 11 for an own one, and b0 div 11 at
 * both points for a repeat (held). Formed in thousandths of a sample, whole
 * numbers, and rounded to float once, as the rule asks the library to (the
 * exact time lies far from any float's rounding edge, so the two agree). */
static double
reference_row_time(int row, int b0, int held, int index) {
    const long point = 80000L * index;
    long thousandths;
    if (row == 24 || row == 25) {
        const int grid = (index == 1 && !held) ? (b0 % 11) : (b0 / 11);
        thousandths = point - 75600L + (11875L * grid);
    } else {
        thousandths = point + (1016L * b0) - 80890L;
    }
    volatile float rounded = (float)((double)thousandths / 1000.0);
    return (double)rounded;
}

/* Model `index`'s pulse (at `time`, frame samples) through its pulsed bands
 * alone, under its own crossfade weight: rising in the interval that ends at
 * its point (index - 1), falling in the one that starts there (index), from
 * its time (or the interval's start) until its response or the interval ends. */
static void
reference_row_render(double out[160], const struct mbe_subframe_model* model, const struct reference_band* band,
                     const int pulsed[57], const double phase[57], int index, double time) {
    for (int k = index - 1; k <= index; ++k) {
        if (k < 0 || k > 1) {
            continue;
        }
        const int rising = k == index - 1;
        for (int n = 0; n < 80; ++n) {
            const double fade = (double)n / 80.0;
            const double t = (80.0 * k) + (double)n - time;
            if (t < 0.0) {
                continue;
            }
            out[(80 * k) + n] +=
                (rising ? fade : (1.0 - fade)) * reference_pulsed_response(model, band, pulsed, phase, t);
        }
    }
}

/* What a stream of the rate-33 rendering carries for the pulsed rows: the row
 * and code of a frame whose own model had a pulsed band (row 0: none). */
struct reference_rows {
    int row, b0;
};

/*
 * The reference's pulsed rows of a frame's three models, in the rate-33
 * rendering alone (in every other no model has a pulsed row): the last frame's
 * own model under the own mask of the row `rows` carries (none from silence),
 * and, for a frame carrying a word (`word`, 0 for none) with no voiced band
 * (cur) on the band-grid model of its code, the mid-frame model under the
 * row's mid-frame mask (its own for a held word) and this frame's own model
 * under its own; each model with a pulsed band emits one pulse, at
 * reference_row_time(). Updates `rows` for the next frame.
 */
static void
reference_row_frame(const mbe_parms* cur, int word, const struct mbe_subframe_model models[3],
                    enum mbe_voiced_synthesis rendering, int from_silence, struct reference_rows* rows,
                    int pulsed[3][57], int emits[3], double time[3]) {
    const int r33 = rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    int row = 0, b0 = 0, held = 0;
    for (int i = 0; i < 3; ++i) {
        emits[i] = 0;
        time[i] = 0.0;
        memset(pulsed[i], 0, 57 * sizeof(int));
    }
    if (r33 && !from_silence && rows->row != 0) {
        emits[0] = reference_row_bands(&models[0], rows->row, 0, pulsed[0]);
        time[0] = reference_row_time(rows->row, rows->b0, 0, 0);
    }
    rows->row = 0;
    if (!r33 || word < TEST_ROW_WORD_BASE || reference_voiced(cur)) {
        return;
    }
    test_row_fields(word, &row, &b0, &held);
    float grid_w0;
    int grid_L;
    r33_grid(b0, &grid_w0, &grid_L);
    if (!same_bits(&grid_w0, &cur->w0, 1) || grid_L != cur->L) {
        return;
    }
    emits[1] = reference_row_bands(&models[1], row, !held, pulsed[1]);
    emits[2] = reference_row_bands(&models[2], row, 0, pulsed[2]);
    time[1] = reference_row_time(row, b0, held, 1);
    time[2] = reference_row_time(row, b0, held, 2);
    if (emits[2]) {
        rows->row = row;
        rows->b0 = b0;
    }
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
reference_mid_pitch(const mbe_parms* prev, const mbe_parms* cur, enum mbe_voiced_synthesis rendering,
                    const int codes[2]) {
    const int prev_voiced = reference_voiced(prev), cur_voiced = reference_voiced(cur);
    if (mbe_frame_mark(cur) == MBE_FRAME_MARK_TONE) {
        return cur->w0;
    }
    if (mbe_frame_mark(cur) == MBE_FRAME_MARK_PREVIOUS_HELD || (prev_voiced && !cur_voiced)) {
        return prev->w0;
    }
    if (cur_voiced && !prev_voiced && mbe_frame_mark(prev) != MBE_FRAME_MARK_PITCH_CARRIES) {
        return cur->w0;
    }
    /* Between two frames with no voiced band of which one has a reserved
     * code that is read: P25 this frame's pitch; D-STAR the frame before's,
     * unless this frame's own pattern pulses (then this frame's, or the
     * geometric mean if the frame before's own pattern pulses too). */
    const int cur_code = codes[1], prev_code = codes[0];
    if (!prev_voiced && !cur_voiced && (cur_code >= 0 || prev_code >= 0)) {
        if (rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE) {
            return cur->w0;
        }
        if (!reference_pattern_pulses(cur_code, 0)) {
            return prev->w0;
        }
        return reference_pattern_pulses(prev_code, 0) ? sqrtf(prev->w0 * cur->w0) : cur->w0;
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

/* x, a quantity a rule rounds to a whole number, away from one: where a build
 * may round it otherwise than this test, the rule's reading is decided only
 * there. */
static double
away_from_whole(double x) {
    assert(fabs(x - floor(x + 0.5)) > 1.0e-6 * fmax(1.0, fabs(x)));
    return x;
}

/* Whether a P25 frame reads its last band as unvoiced in the columns: more
 * than 36 harmonics, harmonic 34 voiced and harmonic 33 not. */
static int
reference_last_band_excepted(const mbe_parms* mp) {
    return mp->L > 36 && mp->Vl[34] == 1 && mp->Vl[33] != 1;
}

/* A frame's voicing sampled at the mid-frame rule's columns, s[0..7] and s[8]
 * past the last (unvoiced), from the rule's statement, in double: a frame with
 * no voiced band unvoiced throughout; D-STAR (`dstar`) both columns of 1 kHz
 * band j at the band's lowest harmonic at or above 1000 j Hz (harmonic 1 for
 * band 0), unvoiced past L or where the band holds no harmonic; P25 harmonic 1
 * for column 0 and the harmonic at or next above 505 + 500 (k - 1) Hz for
 * column k, unvoiced past L, and the last band unvoiced where
 * reference_last_band_excepted() says. */
static void
reference_samples(const mbe_parms* mp, int dstar, int s[9]) {
    const double f0 = (double)mp->w0 * 8000.0 / (2.0 * M_PI);
    memset(s, 0, 9 * sizeof(int));
    if (!reference_voiced(mp)) {
        return;
    }
    if (dstar) {
        for (int j = 0; j < 4; ++j) {
            const int n = (j == 0) ? 1 : (int)ceil(away_from_whole(1000.0 * j / f0));
            const int in_band = (j == 3) || ((double)n * f0 < 1000.0 * (j + 1));
            const int column = j + j;
            s[column] = s[column + 1] = in_band && n <= mp->L && mp->Vl[n] == 1;
        }
        return;
    }
    const int last = reference_last_band_excepted(mp) ? 33 : mp->L;
    s[0] = mp->Vl[1] == 1;
    for (int k = 1; k < 8; ++k) {
        const int n = (int)ceil(away_from_whole((505.0 + 500.0 * (k - 1)) / f0));
        s[k] = n <= last && mp->Vl[n] == 1;
    }
}

/* The column (0..7) of a frequency: D-STAR's 500 k Hz, P25's edges 505 + 500 k Hz. */
static int
reference_column(double hz, int dstar) {
    if (dstar) {
        return (int)fmin(floor(away_from_whole(hz / 500.0)), 7.0);
    }
    int k = 0;
    for (int i = 0; i < 7; ++i) {
        k += hz >= 505.0 + 500.0 * i;
    }
    assert(k == 7 || fabs(hz - (505.0 + 500.0 * k)) > 1.0e-3);
    return k;
}

/* M(k): a column voiced in both frames, or one that changes under a column that changes too. */
static int
reference_column_voiced(const int p[9], const int c[9], int k) {
    return (p[k] && c[k]) || ((p[k] != c[k]) && (p[k + 1] != c[k + 1]));
}

/* A frame's band state (0 unvoiced, 1 voiced, 2 none) at its harmonic nearest w. */
static int
reference_state_at(const mbe_parms* mp, double w) {
    const int nearest = (int)fmin(fmax(floor(w / (double)mp->w0 + 0.5), 1.0), (double)mp->L);
    return (mp->Vl[nearest] == 0) ? 0 : ((mp->Vl[nearest] == 1) ? 1 : 2);
}

/* The rate-33 rendering's hold through a frame lacking a pitch: whether the
 * frame before is voiced at its harmonic at or next above 505 Hz. That hold
 * is the rendering's behaviour from before the columns, kept byte for byte, so
 * it is stated in the library's own float arithmetic: the position
 * ceil(505 Hz in radians per sample over w0), each step rounded to float
 * (volatile, also under x87), then made an index on its bit pattern, a position
 * at or past L (or infinite) read at L and one under 1 at 1, so that no
 * conversion sees a value out of range (at w0 = FLT_MIN the position is about
 * 3.4e37). A double quotient takes another harmonic at some pitches: at the
 * pitch whose bits are 0x3e0761b7 the float quotient is exactly 3 and the double
 * one 3.000000018. */
static int
reference_holds_505(const mbe_parms* prev) {
    const volatile float reach = 505.0f * ((2.0f * (float)M_PI) / 8000.0f);
    const volatile float quotient = reach / prev->w0;
    const volatile float position = ceilf(quotient);
    const float at = position;
    const float last = (float)prev->L;
    uint32_t bits, last_bits, one_bits;
    const float one = 1.0f;
    memcpy(&bits, &at, sizeof(bits));
    memcpy(&last_bits, &last, sizeof(last_bits));
    memcpy(&one_bits, &one, sizeof(one_bits));
    int n = 1;
    if (bits >= last_bits && bits <= UINT32_C(0x7F800000)) {
        n = prev->L;
    } else if (bits >= one_bits && bits < last_bits) {
        n = (int)at;
    }
    return prev->Vl[n] == 1;
}

/* Whether D-STAR's and P25's mid-frame rule of 500 Hz columns sets the
 * mid-frame model's voicing: in their renderings, for a frame with a voiced
 * band or marked as lacking a pitch, whatever the frame before (a tone model
 * included). */
static int
reference_by_columns(const mbe_parms* prev, const mbe_parms* cur, enum mbe_voiced_synthesis rendering) {
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) {
        return 0;
    }
    (void)prev;
    return reference_voiced(cur) || mbe_frame_mark(cur) == MBE_FRAME_MARK_PITCH_ABSENT;
}

/* The state of the mid-frame model's band at angular frequency w by the
 * columns' rule: none where either frame's band has none; the same state in
 * both frames kept; else voiced if and only if the column holding the band's
 * frequency, or the one holding it plus 250 Hz, is voiced at mid-frame. */
static int
reference_column_state(const mbe_parms* prev, const mbe_parms* cur, double w, int dstar) {
    const int before = reference_state_at(prev, w), now = reference_state_at(cur, w);
    if (before == 2 || now == 2) {
        return 2;
    }
    if (before == now) {
        return now;
    }
    int p[9], c[9];
    reference_samples(prev, dstar, p);
    reference_samples(cur, dstar, c);
    const double hz = w * 8000.0 / (2.0 * M_PI);
    return reference_column_voiced(p, c, reference_column(hz, dstar))
           || reference_column_voiced(p, c, reference_column(hz + 250.0, dstar));
}

/* The mid-frame model's band states (0 unvoiced, 1 voiced, 2 none) at pitch
 * w0 and L harmonics, from the rules as docs/dvsi-parity.md states them: in
 * D-STAR's and P25's renderings by the columns' rule
 * (reference_column_state()) where it applies (reference_by_columns()), else
 * each harmonic takes the voicing of the nearest harmonic of the frame whose
 * voicing holds at mid-frame (in the rate-33 rendering the previous frame's
 * when this one has no voiced band and is marked as lacking a pitch, and the
 * previous one is voiced at its harmonic at or next above 505 Hz; else this
 * frame's); a frame marked as a late onset voices none (in the rate-33
 * rendering whatever the frame before, in the others after a frame with no
 * voiced band); and in the rate-33 rendering a frame of V/UV row 17 or 13
 * voices the harmonics under Nyquist in the 500 Hz columns below 1 or 2 (by
 * the frequency of each, or of 250 Hz above it), whatever the voicing around.
 * The marks of the mid-frame voicing (the late onset's reading whatever the
 * frame before, rows 17 and 13) act in their own codec's rendering alone; the
 * columns' rule reads no mark of the frame before, PITCH_ABSENT only on this
 * frame. */
static void
reference_mid_voicing(const mbe_parms* prev, const mbe_parms* cur, double w0, int L,
                      enum mbe_voiced_synthesis rendering, unsigned char voiced[57]) {
    const enum mbe_frame_mark mark = mbe_frame_mark(cur);
    const int rate33 = rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    const int late = mark == MBE_FRAME_MARK_LATE_ONSET && (rate33 || !reference_voiced(prev));
    const int by_columns = reference_by_columns(prev, cur, rendering);
    const mbe_parms* source = cur;
    if (rate33 && mark == MBE_FRAME_MARK_PITCH_ABSENT && !reference_voiced(cur)) {
        source = reference_holds_505(prev) ? prev : cur;
    }
    for (int l = 1; l <= L; ++l) {
        const double w = w0 * l;
        int state = by_columns ? reference_column_state(prev, cur, w, !reference_imbe(rendering))
                               : reference_state_at(source, w);
        if (state == 1 && late) {
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
    if (!reference_imbe(rendering) || mark == MBE_FRAME_MARK_TONE) {
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
 * (reference_mid_pitch()), for the codec's rendering. In D-STAR's and P25's
 * renderings, which read the voicing in 500 Hz columns, every pitch is played
 * 0.1% above the one listed, so that no position the columns' rule rounds (a
 * band's lowest harmonic, a column edge) falls on a whole number, where this
 * test's arithmetic and the library's could round apart. */
struct wave_frame {
    double hz;
    int L, voiced;
    float amplitude;
    enum mbe_frame_mark mark;
    int from_silence;
    int reserved; /* 1 + a reserved code (0..55) whose pitch the frame takes, as its codec decodes it; 0 for none;
                     negative for a pitch of edge_pitch(); a rate-33 pulsed-row word (test_row_word(), 1024 or
                     more), which the frame carries in PHIl[0] on the band-grid model of its code */
};

/* The pitch of a wave_frame whose `reserved` is negative: -1 and -2 the bit
 * patterns 0x3E490FDA (an ulp under D-STAR's b0 124 pitch, 2 pi / 32: its 3rd
 * harmonic lies at 749.999964 Hz, under the 750 Hz a float product rounds
 * it to) and 0x3E3D3C19 (its 17th harmonic is under pi in double, at pi in
 * float), -3 the codec's default code for a frame with no voiced band (P25 b0
 * 25, D-STAR b0 34) and -4 the code this library's encoder sends such a frame
 * with (P25 b0 127, D-STAR b0 123), as its public decoder gives them. */
static float
edge_pitch(enum mbe_voiced_synthesis rendering, int which) {
    const uint32_t bits = (which == -1) ? UINT32_C(0x3E490FDA) : UINT32_C(0x3E3D3C19);
    float w0;
    if (which == -3) {
        return decoded_pitch(!reference_imbe(rendering), reference_imbe(rendering) ? 25 : 34);
    }
    if (which == -4) {
        return decoded_pitch(!reference_imbe(rendering), reference_imbe(rendering) ? 127 : 123);
    }
    memcpy(&w0, &bits, sizeof(w0));
    return w0;
}

/* The rms over the stream of what the reference rendered, all of it and its
 * pulsed state's part. */
struct waveform_levels {
    double rms;
    double pulsed_rms;
};

static struct waveform_levels run_waveform(enum mbe_voiced_synthesis rendering, int inexact,
                                           const struct wave_frame* frames, size_t frame_count);

static void
test_waveform_through_changes(enum mbe_voiced_synthesis rendering, int inexact) {
    static const struct wave_frame frames[] = {
        {200.0, 18, 14, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},          /* an onset */
        {200.0, 18, 14, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},          /* steady */
        {230.0, 16, 16, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},          /* pitch up, fewer harmonics */
        {180.0, 20, 12, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},          /* pitch down, more of them, fewer voiced */
        {207.0, 17, 17, 40.0f, MBE_FRAME_MARK_NONE, 0, 0},          /* a level step */
        {100.0, 30, 0, 8.0f, MBE_FRAME_MARK_PITCH_CARRIES, 0, 0},   /* no voiced band, its pitch carried */
        {180.0, 20, 20, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},          /* onset at the mean of the two pitches */
        {250.0, 14, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, 0, 0},    /* the voicing before lasts to mid-frame */
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},          /* and the train goes on */
        {150.0, 24, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0, 0},     /* after a voiced frame too: the train runs on */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0, 0},             /* voicing ends with the frame's start */
        {127.5, 28, 23, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0, 0}, /* a P25 onset: harmonics 20-23 from mid-frame */
        {127.5, 28, 23, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0, 0}, /* inside the run the mark changes nothing */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {160.0, 25, 3, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0, 0}, /* harmonics 1-3 to 480 Hz: all from mid-frame */
        {100.0, 30, 0, 8.0f, MBE_FRAME_MARK_PITCH_CARRIES, 0, 0},
        {205.0, 17, 17, 16.0f, MBE_FRAME_MARK_ONSET_COLUMNS, 0, 0}, /* at the mean pitch: those under 3505 Hz */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {170.0, 22, 8, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {165.0, 23, 0, 16.0f, MBE_FRAME_MARK_MID_COLUMN_0, 0, 0},    /* rate-33 row 17: column 0 at mid-frame, */
        {150.0, 24, 9, 16.0f, MBE_FRAME_MARK_MID_COLUMNS_0_1, 0, 0}, /* row 13 after it: columns 0-1, the train on */
        {150.0, 24, 9, 16.0f, MBE_FRAME_MARK_MID_COLUMNS_0_1, 0, 0}, /* and inside a voiced run */
        {140.0, 26, 10, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0, 0},     /* row 15 after a voiced frame */
        {130.0, 28, 0, 8.0f, MBE_FRAME_MARK_MID_COLUMN_0, 0, 0},
        {125.0, 29, 10, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0, 0}, /* row 15 after row 17: a new train at mid-frame */
        {90.0, 38, 0, 8.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {120.0, 28, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0, 0}, /* voiced bands from mid-frame */
        {60.0, 56, 40, 6.0f, MBE_FRAME_MARK_NONE, 0, 0},        /* the longest period: responses span frames */
        {60.0, 56, 40, 6.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {110.0, 30, 30, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {250.0, 14, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, 0, 0}, /* the train runs on through this frame, */
        {65.0, 52, 6, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 0, 0},   /* then starts anew at mid-frame: its pulses, */
        {65.0, 52, 30, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},        /* under a period longer than a subframe, stay out */
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {250.0, 15, 0, 1.6f, MBE_FRAME_MARK_PREVIOUS_HELD, 0, 0}, /* the model before stands at mid-frame, */
        {170.0, 22, 22, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},        /* and the train runs on into the frame after */
        {140.0, 26, 4, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},         /* voiced low only, held the same way */
        {250.0, 15, 0, 1.6f, MBE_FRAME_MARK_PREVIOUS_HELD, 0, 0},
        {250.0, 15, 0, 1.6f, MBE_FRAME_MARK_PREVIOUS_HELD, 0, 0}, /* a second such frame stops the train */
        {170.0, 22, 22, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},        /* and a new one starts */
        {200.0, 18, 18, 16.0f, MBE_FRAME_MARK_NONE, 1, 0},        /* from silence: a new train at the start, */
        {200.0, 18, 18, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},        /* and on from it */
        {160.0, 25, 9, 30.0f, MBE_FRAME_MARK_NONE, 1, 0},         /* from silence at another pitch and level */
        {250.0, 15, 6, 4.0f, MBE_FRAME_MARK_PREVIOUS_HELD, 1, 0}, /* the model before is not held: its own rises */
        {140.0, 26, 26, 16.0f, MBE_FRAME_MARK_NONE, 1, 0},        /* the train it ran is dropped */
        {100.0, 30, 0, 8.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {120.0, 28, 3, 16.0f, MBE_FRAME_MARK_LATE_ONSET, 1, 0}, /* a late onset from silence */
        {250.0, 14, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, 0, 0},
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 1, 0}, /* the marks before do not carry pulses in */
        {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},
        {56.0, 56, 40, 16.0f, MBE_FRAME_MARK_NONE, 0, 0}, /* each subframe's period 1.64 times the last: octaves down */
        {170.0, 22, 20, 16.0f, MBE_FRAME_MARK_NONE, 0, 0}, /* 0.57 times: octaves up */
        {70.8, 52, 40, 16.0f, MBE_FRAME_MARK_NONE, 0, 0},  /* 1.55 times: a glide */
    };
    const struct waveform_levels levels = run_waveform(rendering, inexact, frames, sizeof(frames) / sizeof(frames[0]));
    /* The comparison saw speech, at a level that makes the tolerance a
     * thousandth of it; no frame here has a reserved code. */
    assert(levels.rms > 50.0);
    assert(levels.pulsed_rms <= 0.0);
}

/* The reserved codes of each kind: own pattern only, mid-frame pattern only,
 * both. */
static const int k_own_codes[] = {0,  1,  2,  3,  4,  5,  6,  7,  24, 25, 26, 27,
                                  28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39};
static const int k_mid_codes[] = {8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                  20, 21, 22, 23, 40, 41, 42, 43, 44, 45, 46, 47};
static const int k_both_codes[] = {48, 49, 50, 51, 52, 53, 54, 55};

/* A frame of kind `kind` for the pulsed state's waveform stream: 0 a frame
 * with no voiced band and no code (lacking a pitch, as its decoder marks
 * one; at a pitch no code decodes to), 1 to 3 a frame with a code whose own
 * pattern only, mid-frame pattern only or both pulse (the `n`th of its kind,
 * cycling), 4 a voiced frame, 5 a tone model, 6 a frame of kind 3 marked as a
 * tone model, 7 a frame at the default code, 8 a frame with no voiced band at
 * the code this library's encoder sends such a frame with. A
 * frame with a code is marked as D-STAR's decoder marks it in D-STAR's
 * rendering (MBE_FRAME_MARK_PITCH_CARRIES when its own pattern pulses) and
 * carries no mark elsewhere. */
/* A frame with a reserved code of kind 1 to 3 (pulsed_wave_frame()'s kinds), the
 * `n`th of its kind, cycling. */
static struct wave_frame
code_frame(enum mbe_voiced_synthesis rendering, int code, int from_silence) {
    struct wave_frame frame = {125.0, 30, 0, 8.0f, MBE_FRAME_MARK_NONE, from_silence, 0};
    const float w0 = reserved_pitch(rendering, code);
    const int L = (int)floor(0.9254 * M_PI / (double)w0);
    frame.reserved = 1 + code;
    frame.L = (L > 56) ? 56 : L;
    frame.mark = (rendering == MBE_VOICED_SYNTHESIS_DVSI && reference_pattern_pulses(code, 0))
                     ? MBE_FRAME_MARK_PITCH_CARRIES
                     : MBE_FRAME_MARK_NONE;
    return frame;
}

static struct wave_frame
code_wave_frame(enum mbe_voiced_synthesis rendering, int kind, size_t n, int from_silence) {
    const int* codes = (kind == 1) ? k_own_codes : ((kind == 2) ? k_mid_codes : k_both_codes);
    const size_t count = (kind == 3) ? 8u : 24u;
    return code_frame(rendering, codes[(n * 7u) % count], from_silence);
}

static struct wave_frame
pulsed_wave_frame(enum mbe_voiced_synthesis rendering, int kind, size_t n, int from_silence) {
    struct wave_frame frame = {125.0, 30, 0, 8.0f, MBE_FRAME_MARK_PITCH_ABSENT, from_silence, 0};
    if (kind == 5) {
        /* A tone model: two voiced harmonics, its own model at mid-frame. */
        frame.hz = 200.0;
        frame.L = 18;
        frame.voiced = 2;
        frame.amplitude = 20.0f;
        frame.mark = MBE_FRAME_MARK_TONE;
        return frame;
    }
    if (kind == 6) {
        /* A caller's TONE mark on a frame with no voiced band at a reserved
         * pitch: its own model has the code, the model it puts at mid-frame
         * takes no pattern. */
        frame = code_wave_frame(rendering, 3, n, from_silence);
        frame.mark = MBE_FRAME_MARK_TONE;
        return frame;
    }
    if (kind == 8) {
        /* An ordinary code's pitch, outside the table: P25's counter falls. */
        const int L = (int)floor(0.9254 * M_PI / (double)edge_pitch(rendering, -4));
        frame.L = (L > 56) ? 56 : L;
        frame.mark = MBE_FRAME_MARK_NONE;
        frame.reserved = -4;
        return frame;
    }
    if (kind == 7) {
        /* The default code: no pattern, but P25's validity counter counts it. */
        frame.L = 14;
        frame.mark = MBE_FRAME_MARK_NONE;
        frame.reserved = -3;
        return frame;
    }
    if (kind == 4) {
        frame.hz = 150.0;
        frame.L = 24;
        frame.voiced = 20;
        frame.amplitude = 16.0f;
        frame.mark = MBE_FRAME_MARK_NONE;
        return frame;
    }
    return (kind == 0) ? frame : code_wave_frame(rendering, kind, n, from_silence);
}

/* The independent rendering through the pulsed state (test_waveform_through_changes()'s
 * reference with the rule of docs/dvsi-parity.md, "D-STAR and P25: the pulsed
 * state"): frames with no voiced band without a code, and with codes whose own
 * pattern, mid-frame pattern or both pulse, in every order of two and three (a
 * de Bruijn sequence), voiced frames before and after each kind, tone models
 * (MBE_FRAME_MARK_TONE, their own model at mid-frame) before and after frames
 * with codes, a caller's tone mark on a frame with a code (its mid-frame model
 * takes no pattern), and frames from silence; frames at the default code
 * (kind 7) first and with P25's validity counter walked down to its floor by
 * frames with no voiced band at an ordinary code's pitch (kind 8), up through
 * the codes it does not read and to its ceiling by default-code frames, and
 * kept where it stands by frames at a pitch no code decodes to (kind 0; the
 * reference keeps the counter, reference_next_count()). In rate 33's rendering and in ProVoice's
 * the same pitches have no code, and the reference plays no pulsed band
 * there. */
static void
test_waveform_through_pulsed_codes(enum mbe_voiced_synthesis rendering, int inexact) {
    static const char kinds[] = "777000100200301101201302102202303103203311121131221231321332223233300"
                                "4140441140424044224043404433404"
                                "1515525535051615262536361"
                                "888888888883333333333388888888888777777777723"
                                "8880000000000032";
    static const int silent[][2] = {{1, 1}, {3, 1}, {4, 1}, {2, 0}, {3, 1}, {2, 1}, {0, 0}};
    struct wave_frame frames[sizeof(kinds) + sizeof(silent) / sizeof(silent[0])];
    size_t count = 0;
    for (size_t i = 0; kinds[i] != '\0'; ++i) {
        frames[count] = pulsed_wave_frame(rendering, kinds[i] - '0', i, 0);
        count++;
    }
    for (size_t i = 0; i < sizeof(silent) / sizeof(silent[0]); ++i) {
        frames[count] = pulsed_wave_frame(rendering, silent[i][0], i + 3, silent[i][1]);
        count++;
    }
    const struct waveform_levels levels = run_waveform(rendering, inexact, frames, count);
    /* The comparison saw the pulses, at a level that makes the tolerance about
     * a thousandth of them, where the rendering has the pulsed state; where
     * it has none the same pitches play no pulse. */
    assert(levels.rms > 35.0);
    const int pulsed = rendering == MBE_VOICED_SYNTHESIS_DVSI || rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    assert(pulsed ? levels.pulsed_rms > 30.0 : levels.pulsed_rms <= 0.0);
}

/* A rate-33 frame of V/UV row `row` at pitch code b0 (held: a repeat's word),
 * on the band-grid model of the code, with no voiced band (voiced: the first
 * `voiced` harmonics voiced, a word on a model it is not the word of). */
static struct wave_frame
row_wave_frame(int row, int b0, int held, int from_silence, int voiced) {
    float w0;
    int L;
    r33_grid(b0, &w0, &L);
    const struct wave_frame frame = {
        125.0, L, voiced, 8.0f, MBE_FRAME_MARK_NONE, from_silence, test_row_word(row, b0, held)};
    return frame;
}

/* The rows' waveform stream, built frame by frame. */
struct row_stream {
    struct wave_frame frames[3600];
    size_t count;
    int last_row, last_b0; /* the last pulsed-row frame appended (for its repeats) */
    size_t n;              /* how many of each kind so far, for cycling */
};

static void
row_push(struct row_stream* s, struct wave_frame frame) {
    assert(s->count < sizeof(s->frames) / sizeof(s->frames[0]));
    s->frames[s->count++] = frame;
}

static void
row_push_row(struct row_stream* s, int row, int b0, int from_silence) {
    row_push(s, row_wave_frame(row, b0, 0, from_silence, 0));
    s->last_row = row;
    s->last_b0 = b0;
}

/*
 * One frame of kind `kind`: P a frame with no voiced band and no word (V/UV
 * row 16, its coded pitch), O a frame of a row pulsing its own subframe only
 * (18-23), M one pulsing its mid-frame subframe only (26-31), B one pulsing
 * both (24, 25), V a voiced frame, S a silence frame (MBE_FRAME_MARK_PREVIOUS_HELD),
 * H a repeat of the last pulsed-row frame (its held word), W a word of rows
 * 18-25 on a model with a voiced band; lower case o, m, b, p and v the same
 * from silence. Rows and codes cycle so that the stream visits many.
 */
static void
row_push_kind(struct row_stream* s, char kind) {
    static const int codes[] = {0, 13, 37, 55, 64, 79, 80, 87, 98, 105, 119, 47, 71, 26, 110, 3};
    const size_t n = s->n++;
    const int from_silence = kind >= 'a' && kind <= 'z';
    const int b0 = codes[n % (sizeof(codes) / sizeof(codes[0]))];
    struct wave_frame plain = {125.0, 30, 0, 8.0f, MBE_FRAME_MARK_NONE, from_silence, 0};
    struct wave_frame voiced = {150.0, 24, 20, 16.0f, MBE_FRAME_MARK_NONE, from_silence, 0};
    const struct wave_frame silence = {250.0, 15, 0, 1.6f, MBE_FRAME_MARK_PREVIOUS_HELD, 0, 0};
    switch (kind) {
        case 'P':
        case 'p': row_push(s, plain); break;
        case 'V':
        case 'v': row_push(s, voiced); break;
        case 'S': row_push(s, silence); break;
        case 'O':
        case 'o': row_push_row(s, 18 + (int)(n % 6u), b0, from_silence); break;
        case 'M':
        case 'm': row_push_row(s, 26 + (int)(n % 6u), b0, from_silence); break;
        case 'B':
        case 'b': row_push_row(s, 24 + (int)(n % 2u), b0, from_silence); break;
        case 'H': row_push(s, row_wave_frame(s->last_row, s->last_b0, 1, 0, 0)); break;
        default: row_push(s, row_wave_frame(18 + (int)(n % 8u), b0, 0, 0, 5)); break;
    }
}

/*
 * The independent rendering through rate 33's pulsed rows
 * (test_waveform_through_changes()'s reference with the rule of
 * docs/dvsi-parity.md, "Rate 33: the pulsed rows"): each of the fourteen rows
 * at every code 0..119 (all 1,680 pairs, each after and before a frame with
 * no voiced band, so that each pair's playback is checked sample for sample,
 * not only its word), so that the pulses fall before their model's window (b0
 * 0: the mid-frame pulse at -0.89, before the frame; the own pulse at 79.11,
 * before its window), in the rising half, on either side of the point (b0 79
 * puts the single-row pulse 0.626 sample before it and b0 80 0.390 after it:
 * no legal code puts a pulse exactly on a model's point), in the falling half
 * and, for an own model, in the next frame (b0 80 and above), b0 0 and 119
 * among them (L 9 and 56, the grid's two ends), rows 24 and 25 at every pair of
 * grid indices that occurs; plain, own-pulsed, mid-pulsed and
 * both-pulsed frames in every order of two and three (a de Bruijn sequence);
 * voiced frames before and after each kind (the train exactly as without the
 * word: the reference's train is the one with no pulsed row); silence frames
 * after pulsed frames; repeats (held words) of each kind, twice for own-only
 * and both-pulsed rows, and the frame after; frames from silence of each kind
 * and after pulsed frames; and words on models with voiced bands, which play
 * as no word. In D-STAR's, P25's and ProVoice's renderings the same frames,
 * words included, play no pulsed row.
 */
static void
test_waveform_through_pulsed_rows(enum mbe_voiced_synthesis rendering, int inexact) {
    static const char* const sequences[] = {
        "PPPOPPMPPBPOOPOMPOBPMOPMMPMBPBOPBMPBBOOOMOOBOMMOMBOBMOBBMMMBMBBBPP",
        "VOVVMVVBVVPVVOOVVMMVVBBVOVMVBVP",
        "OSPBSPMSPOSSPBSVMSV",
        "OHHPBHHPMHPOHVBHSMHHP",
        "OoPBbPMmPOpPBvPMbPoMPbOP",
        "WPOWPWWVWBP",
    };
    static struct row_stream stream;
    memset(&stream, 0, sizeof(stream));
    row_push_kind(&stream, 'P');
    for (int row = 18; row <= 31; ++row) {
        for (int b0 = 0; b0 < 120; ++b0) {
            row_push_row(&stream, row, b0, 0);
            row_push_kind(&stream, 'P');
        }
    }
    for (size_t q = 0; q < sizeof(sequences) / sizeof(sequences[0]); ++q) {
        for (const char* k = sequences[q]; *k != '\0'; ++k) {
            row_push_kind(&stream, *k);
        }
    }
    const struct waveform_levels levels = run_waveform(rendering, inexact, stream.frames, stream.count);
    assert(levels.rms > 1.0);
    assert((rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) ? levels.pulsed_rms > 1.0 : levels.pulsed_rms <= 0.0);
}

/*
 * D-STAR's pulsed bands at a channel's edge and at Nyquist, against the
 * independent rendering: after a frame with no voiced band at 0x3E490FDA (no
 * code), a code-12 frame (b0 76: no own pattern, so the mid-frame model takes
 * the frame before's pitch; its mid-frame pattern pulses channel 2, which the
 * 250 Hz reach makes [750, 1500) Hz): the mid-frame model's 3rd harmonic lies
 * at 749.999964 Hz and is not pulsed, where its frequency rounded to float
 * (750) would pulse it. In IEEE builds also, after a frame of 17 harmonics at
 * 0x3E3D3C19, a code-16 frame (b0 80, its mid-frame pattern channel 7): the
 * mid-frame model, at that pitch with 17 harmonics, renders no 17th harmonic,
 * which is under pi in double but at it in float.
 */
static void
test_pulsed_channel_edges(int inexact) {
    struct wave_frame frames[6];
    size_t count = 0;
    const struct wave_frame edge = {125.0, 14, 0, 8.0f, MBE_FRAME_MARK_NONE, 0, -1};
    frames[count++] = edge;
    frames[count++] = code_frame(MBE_VOICED_SYNTHESIS_DVSI, 12, 0);
    if (!inexact) {
        const struct wave_frame nyquist = {125.0, 17, 0, 8.0f, MBE_FRAME_MARK_NONE, 0, -2};
        frames[count++] = nyquist;
        frames[count++] = code_frame(MBE_VOICED_SYNTHESIS_DVSI, 16, 0);
    }
    const struct waveform_levels levels = run_waveform(MBE_VOICED_SYNTHESIS_DVSI, inexact, frames, count);
    assert(levels.pulsed_rms > 1.0);
}

/* The reference's pulsed state of a frame's three models (pulsed[i][l]) and
 * which of the first two emit, from the rule: each model's code from its own
 * pitch and voicing (the last frame's own model, this frame's own model, and
 * this frame's for the mid-frame one) where the rendering reads it
 * (reference_read_codes(): `codes`), the own pattern for the frames' own
 * models and the mid-frame pattern for the mid-frame model, which takes none
 * where the frame's mark puts another model there (the model before held, a
 * tone's own); a model with a pulsed band emits unless the model before it
 * did. `stream` carries what the frame leaves the next. */
static void
reference_pulsed_frame(const mbe_parms* cur, const struct mbe_subframe_model models[3],
                       enum mbe_voiced_synthesis rendering, const int codes[2], int pulsed[3][57], int emits[2],
                       struct reference_stream* stream) {
    const enum mbe_frame_mark mark = mbe_frame_mark(cur);
    const int replaced = mark == MBE_FRAME_MARK_PREVIOUS_HELD || mark == MBE_FRAME_MARK_TONE;
    const int any0 = reference_pulse_bands(&models[0], codes[0], 0, pulsed[0]);
    const int any1 = reference_pulse_bands(&models[1], replaced ? -1 : codes[1], 1, pulsed[1]);
    (void)reference_pulse_bands(&models[2], codes[1], 0, pulsed[2]);
    emits[0] = any0 && !stream->emitted;
    emits[1] = any1 && !emits[0];
    stream->emitted = emits[1];
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE) {
        stream->read = codes[1] >= 0;
        stream->count = reference_next_count(cur, stream->count);
    }
}

/* The noiseOverlap elements the pulsed state is stored in, and the values,
 * as this test pins them: whether the mid-frame model emitted (1.0f), whether
 * the frame's code was read (1.0f) and P25's validity counter (a whole number
 * as a float), the values every other case leaves 0. */
#define TEST_EMITTED_SLOT 34
#define TEST_READ_SLOT    35
#define TEST_COUNT_SLOT   36
#define TEST_ONE_BITS     UINT32_C(0x3F800000)

/* Whether a frame's stored pulsed-row word (slot 37) is the reference's: the
 * word, not held, of the row and code `rows` carries, or 0. */
static int
stored_row_as(const mbe_parms* mp, const struct reference_rows* rows) {
    uint32_t bits;
    memcpy(&bits, &mp->noiseOverlap[TEST_ROW_SLOT], sizeof(bits));
    return bits == ((rows->row != 0) ? test_row_word_bits(test_row_word(rows->row, rows->b0, 0)) : 0u);
}

/* Whether a frame's stored pulsed state is the reference's. */
static int
stored_as(const mbe_parms* mp, const struct reference_stream* stream, enum mbe_voiced_synthesis rendering) {
    uint32_t emitted, read, count;
    const float want_count = (rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE) ? (float)stream->count : 0.0f;
    uint32_t want_count_bits;
    memcpy(&emitted, &mp->noiseOverlap[TEST_EMITTED_SLOT], sizeof(emitted));
    memcpy(&read, &mp->noiseOverlap[TEST_READ_SLOT], sizeof(read));
    memcpy(&count, &mp->noiseOverlap[TEST_COUNT_SLOT], sizeof(count));
    memcpy(&want_count_bits, &want_count, sizeof(want_count_bits));
    const int p25 = rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    return emitted == (stream->emitted ? TEST_ONE_BITS : 0u) && read == ((p25 && stream->read) ? TEST_ONE_BITS : 0u)
           && count == want_count_bits;
}

static struct waveform_levels
run_waveform(enum mbe_voiced_synthesis rendering, int inexact, const struct wave_frame* frames, size_t frame_count) {
    struct reference_train train = {0.0, 0, {0.0}};
    mbe_parms last;
    double energy = 0.0, pulsed_energy = 0.0, worst = 0.0;
    struct reference_stream stream = {0, 0, 0};
    struct reference_rows rows = {0, 0};
    struct reference_tracker tracker = {0.0, 0.0};

    fresh_stream(&last);
    for (size_t f = 0; f < frame_count; ++f) {
        struct mbe_subframe_model models[3];
        mbe_parms model, cur;
        float out[160];
        double want[160];
        const double start = 160.0 * (double)f;

        const double shift = (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) ? 1.0 : 1.001;
        set_model(&model, frames[f].hz * shift, frames[f].L, frames[f].voiced, frames[f].amplitude);
        const int word = (frames[f].reserved >= TEST_ROW_WORD_BASE) ? frames[f].reserved : 0;
        if (word != 0) {
            int row, b0, held, grid_L;
            test_row_fields(word, &row, &b0, &held);
            r33_grid(b0, &model.w0, &grid_L);
            assert(grid_L == model.L);
        } else if (frames[f].reserved > 0) {
            model.w0 = reserved_pitch(rendering, frames[f].reserved - 1);
        } else if (frames[f].reserved < 0) {
            model.w0 = edge_pitch(rendering, frames[f].reserved);
        }
        shape_model(&model, f);
        mbe_set_frame_mark(&model, frames[f].mark);
        cur = last;
        cur.w0 = model.w0;
        cur.L = model.L;
        memcpy(cur.Vl, model.Vl, sizeof(cur.Vl));
        copy_amplitudes(&cur, &model);
        const enum mbe_frame_mark mark = reference_rendered_mark(&cur, frames[f].mark, rendering);
        mbe_set_frame_mark(&cur, mark);
        if (word != 0) {
            /* The word in place of a mark: no mark acts (P25's rendering
             * marks a voiced model all the same). */
            assert(mark == MBE_FRAME_MARK_NONE || mark == MBE_FRAME_MARK_ONSET_COLUMNS);
            if (mark == MBE_FRAME_MARK_NONE) {
                cur.PHIl[0] = (float)word;
            }
            model.PHIl[0] = (float)word;
        }
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
        int codes[2];
        reference_read_codes(&last, &cur, rendering, &stream, codes);
        const float mid_pitch = reference_mid_pitch(&last, &cur, rendering, codes);
        assert(same_bits(&models[1].w0, &mid_pitch, 1) || (inexact && within_an_ulp(models[1].w0, mid_pitch)));
        assert(last.L >= 1 && last.L <= 56 && cur.L >= 1 && cur.L <= 56);
        /* A tone frame's own model stands at mid-frame (MBE_FRAME_MARK_TONE). */
        const int tone = mark == MBE_FRAME_MARK_TONE;
        const int mid_count = (mark == MBE_FRAME_MARK_PREVIOUS_HELD)
                                  ? last.L
                                  : (tone ? cur.L
                                          : reference_mid_count(mid_pitch, (last.L < cur.L) ? last.L : cur.L,
                                                                (last.L < cur.L) ? cur.L : last.L));
        assert(models[1].L == mid_count || (inexact && count_margin(mid_pitch) < 1.0e-5));
        assert(mid_count >= 1 && mid_count <= 56);
        models[1].w0 = mid_pitch;
        models[1].L = mid_count;
        /* Its voicing is the reference's own, which the library's must equal. */
        if (mark != MBE_FRAME_MARK_PREVIOUS_HELD) {
            unsigned char voiced[57] = {0};
            reference_mid_voicing(&last, &cur, (double)mid_pitch, models[1].L, rendering, voiced);
            for (int l = 1; l <= models[1].L && tone; ++l) {
                voiced[l] = (unsigned char)((cur.Vl[l] == 0) ? 0 : ((cur.Vl[l] == 1) ? 1 : 2));
            }
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
        } else if (tone) {
            memcpy(phases[1], phases[2], sizeof(phases[1]));
        } else {
            reference_mid_phases(&last, &cur, (double)models[1].w0, models[1].L, phases[1]);
        }

        /* The train starts afresh after a frame with no voiced band, unless that
         * frame lacked a pitch or held the model before it to mid-frame, or, in
         * the rate-33 rendering alone, voiced column 0 at mid-frame (row 17).
         * A frame that starts from silence plays nothing of the last frame's
         * own model and starts a new train, whatever came before. */
        const enum mbe_frame_mark last_mark = mbe_frame_mark(&last);
        const int runs_on =
            last_mark == MBE_FRAME_MARK_PITCH_ABSENT || last_mark == MBE_FRAME_MARK_PREVIOUS_HELD
            || (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 && last_mark == MBE_FRAME_MARK_MID_COLUMN_0);
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
        } else if (!models[0].any_voiced && !runs_on) {
            train.phase = 0.0;
            train.count = 0;
        }
        /* The pulsed state: its own pulses, in pulsed bands only, beside the
         * voiced train, which models with pulsed bands and no voiced band leave
         * as models with no voiced band leave it. */
        int pulsed[3][57], emits[2];
        reference_pulsed_frame(&cur, models, rendering, codes, pulsed, emits, &stream);
        int row_pulsed[3][57], row_emits[3];
        double row_time[3];
        reference_row_frame(&cur, word, models, rendering, frames[f].from_silence, &rows, row_pulsed, row_emits,
                            row_time);
        /* Pulses more than two frames back sound no longer. */
        int kept = 0;
        for (int i = 0; i < train.count; ++i) {
            if (train.time[i] > start - 320.0) {
                train.time[kept++] = train.time[i];
            }
        }
        train.count = kept;
        /* The low band: the filter of the rendering, and the cut at each
         * model's corner, from the reference's own tracker. */
        struct reference_band band[3];
        reference_tracker_frame(&tracker, &cur, rendering, band);
        reference_advance(&train, start, &models[0], &models[1]);
        reference_render(want, start, &train, &models[0], &band[0], phases[0], &models[1], &band[1], phases[1]);
        reference_advance(&train, start + 80.0, &models[1], &models[2]);
        reference_render(want + 80, start + 80.0, &train, &models[1], &band[1], phases[1], &models[2], &band[2],
                         phases[2]);
        double pulses[160] = {0.0};
        if (emits[0]) {
            reference_pulsed_render(pulses, &models[0], &band[0], pulsed[0], phases[0], &models[1], &band[1], pulsed[1],
                                    phases[1]);
        }
        if (emits[1]) {
            reference_pulsed_render(pulses + 80, &models[1], &band[1], pulsed[1], phases[1], &models[2], &band[2],
                                    pulsed[2], phases[2]);
        }
        for (int i = 0; i < 3; ++i) {
            if (row_emits[i]) {
                reference_row_render(pulses, &models[i], &band[i], row_pulsed[i], phases[i], i, row_time[i]);
            }
        }
        for (int n = 0; n < 160; ++n) {
            want[n] += pulses[n];
            pulsed_energy += pulses[n] * pulses[n];
        }

        if (word != 0) {
            next_frame_with_word(out, &last, &model, frames[f].from_silence, rendering);
        } else {
            next_frame_as(out, &last, &model, frames[f].from_silence, rendering);
        }
        for (int n = 0; n < 160; ++n) {
            const double difference = fabs((double)out[n] - want[n]);
            worst = (difference > worst) ? difference : worst;
            energy += want[n] * want[n];
            assert(difference < 0.05 + (1.0e-4 * fabs(want[n])));
        }
        assert(fabs(cycles_apart(pulse_phase(&last), train.phase)) < 1e-4);
        assert(stored_as(&last, &stream, rendering));
        assert(stored_row_as(&last, &rows));
        assert(stored_tracker_as(&last, &tracker, rendering));
    }
    assert(worst < 0.05);
    const struct waveform_levels levels = {sqrt(energy / (160.0 * (double)frame_count)),
                                           sqrt(pulsed_energy / (160.0 * (double)frame_count))};
    return levels;
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

static int
same_model(const struct mbe_subframe_model* a, const struct mbe_subframe_model* b) {
    return same_bits(&a->w0, &b->w0, 1) && a->L == b->L && a->any_voiced == b->any_voiced
           && memcmp(a->voiced, b->voiced, sizeof(a->voiced)) == 0 && same_bits(a->magnitude, b->magnitude, 57)
           && same_bits(a->phase, b->phase, 57);
}

/* Whether the mid-frame model between `before` and `frame`, rendered the way
 * `rendering` says, has a voiced band. */
static int
mid_is_voiced_as(const mbe_parms* before, const mbe_parms* frame, enum mbe_voiced_synthesis rendering) {
    struct mbe_subframe_model mid;
    mbe_subframe_mid_model(before, frame, 0, rendering, &mid);
    return mid.any_voiced;
}

static int
mid_is_voiced(const mbe_parms* before, const mbe_parms* frame) {
    return mid_is_voiced_as(before, frame, MBE_VOICED_SYNTHESIS_DVSI);
}

/* DVSI's D-STAR and P25 decoders hold the voicing of the frame before a frame
 * that has no voiced band and no pitch through that frame's first half in the
 * 500 Hz columns that frame covers (voiced at both their ends): the pulse
 * train runs on at the old pitch through the whole frame, and a voiced frame
 * after it continues the train rather than starting a new one. In D-STAR's
 * rendering (bands of 1 kHz, both columns of each taking its bit) a frame
 * voiced to 1884 Hz at 157 Hz, b1 1100 as the decoder would give it, covers
 * columns 0 to 2: harmonics 1 to 9 (to 1413 Hz, the 9th's column 2), not 10 to
 * 12 (1570 Hz and up, in column 3, the top of the run, or 250 Hz below it). */
static void
test_voicing_held_by_columns_without_a_pitch(void) {
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
        assert(mid.voiced[l] == ((l <= 9) ? 1 : 0));
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

    /* P25's columns start at 505 Hz (MBE_VOICING_REACH_HZ): three voiced
     * harmonics cover column 0 at a period of 47.25 samples (the third at
     * 507.9 Hz), which holds harmonics 1 and 2, and none at 47.75 (502.6 Hz).
     * There voiced harmonics 7 to 12 above the unvoiced ones (1173 to 2010 Hz)
     * cover column 3 (from 1505 Hz, harmonic 9, to 2005 Hz, harmonic 12) and
     * hold harmonics 8 (1340 Hz, within 250 Hz below it) to 11. D-STAR's
     * rendering holds column 0 from three harmonics at either period. */
    const enum mbe_voiced_synthesis p25 = MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    mbe_set_frame_mark(&unvoiced, MBE_FRAME_MARK_PITCH_ABSENT);
    set_model(&last, 8000.0 / 47.25, 21, 3, 16.0f);
    mbe_subframe_mid_model(&last, &unvoiced, 0, p25, &mid);
    for (int l = 1; l <= mid.L; ++l) {
        assert(mid.voiced[l] == ((l <= 2) ? 1 : 0));
    }
    assert(mid_is_voiced(&last, &unvoiced));
    set_model(&last, 8000.0 / 47.75, 22, 3, 16.0f);
    assert(!mid_is_voiced_as(&last, &unvoiced, p25));
    assert(mid_is_voiced(&last, &unvoiced));
    for (int l = 7; l <= 12; ++l) {
        last.Vl[l] = 1;
    }
    mbe_subframe_mid_model(&last, &unvoiced, 0, p25, &mid);
    for (int l = 1; l <= mid.L; ++l) {
        assert(mid.voiced[l] == ((l >= 8 && l <= 11) ? 1 : 0));
    }
    /* A frame with a voiced band is not such a frame, whatever its mark: its
     * mid-frame model is that of the same frame unmarked. */
    struct mbe_subframe_model plain;
    set_model(&last, 157.0, 23, 12, 16.0f);
    unvoiced.Vl[2] = 1;
    mbe_subframe_mid_model(&last, &unvoiced, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    mbe_set_frame_mark(&unvoiced, MBE_FRAME_MARK_NONE);
    mbe_subframe_mid_model(&last, &unvoiced, 0, MBE_VOICED_SYNTHESIS_DVSI, &plain);
    assert(mid.any_voiced && same_model(&mid, &plain));
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

    /* Both voiced: the geometric-mean pitch, the voicing of the columns' rule
     * (reference_mid_voicing()), magnitudes the geometric mean of the two
     * frames'. */
    set_model(&prev, 151.0, 24, 24, 10.0f);
    set_model(&cur, 181.0, 20, 10, 20.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    const double w0 = sqrt((double)prev.w0 * (double)cur.w0);
    assert(fabs((double)mid.w0 - w0) < 1e-6);
    assert(mid.L == 22 && mid.any_voiced);
    unsigned char states[57] = {0};
    reference_mid_voicing(&prev, &cur, (double)mid.w0, mid.L, MBE_VOICED_SYNTHESIS_DVSI, states);
    for (int l = 1; l <= mid.L; ++l) {
        const double w = (double)l * w0;
        const int voiced = states[l] == 1;
        assert(mid.voiced[l] == voiced);
        /* An unvoiced band's magnitude carries its factor over sqrt(w0): out
         * of each frame's before the mean, into the mid-frame model's after. */
        const double mean = 0.5 * (common_log2_at(&prev, w) + common_log2_at(&cur, w));
        const double want = exp2(mean + (voiced ? 0.0 : log2(k / sqrt(w0))));
        assert(fabs((double)mid.magnitude[l] / want - 1.0) < 1e-4);
    }
    /* A band voiced in both frames; one voiced before and unvoiced now that
     * its column holds (harmonic 20, 3306 Hz: column 6 changes, as column 7
     * above it does); and one that its column lets go (harmonic 22, 3637 Hz:
     * column 7, the top of the run of changing columns). */
    assert(mid.voiced[5] == 1 && mid.voiced[20] == 1 && mid.voiced[22] == 0);
    assert(fabs((double)mid.magnitude[5] / sqrt(200.0) - 1.0) < 1e-4);
    assert(fabs((double)mid.magnitude[20] / sqrt(10.0 * 20.0 * sqrt((double)cur.w0) / k) - 1.0) < 1e-4);
    assert(fabs((double)mid.magnitude[22] / (sqrt(10.0 * 20.0 * sqrt((double)cur.w0) / k) * k / sqrt(w0)) - 1.0)
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
    set_model(&prev, 151.0, 24, 0, 10.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(same_bits(&mid.w0, &cur.w0, 1) && mid.any_voiced);
    set_model(&prev, 151.0, 24, 24, 10.0f);
    set_model(&cur, 181.0, 20, 0, 20.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(same_bits(&mid.w0, &prev.w0, 1) && !mid.any_voiced);
    set_model(&prev, 151.0, 24, 0, 10.0f);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(fabs((double)mid.w0 - w0) < 1e-6 && !mid.any_voiced);
    assert(fabs((double)mid.magnitude[5] / sqrt(200.0) - 1.0) < 1e-4);

    /* A late onset has no voiced band at mid-frame; its bands are noise. */
    set_model(&cur, 181.0, 20, 20, 20.0f);
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

    /* D-STAR's filter; the corner stays under 187.5 Hz, so no harmonic is cut. */
    const struct reference_band band = reference_band_of(MBE_VOICED_SYNTHESIS_DVSI);
    double moved = 0.0;
    for (int l = 1; l <= 18; ++l) {
        const double hz = 187.5 * l;
        for (int k = 0; k < 3; ++k) {
            assert(last[k].noiseOverlap[MBE_LOW_CUT_TRACK_SLOT] < 187.0f);
            const double want = 2.0 * (double)model[k].Ml[l] * reference_band_gain(&band, hz);
            assert(fabs(component_amplitude(out[k] + 320, 1280, hz) / want - 1.0) < 2e-3);
            const double stored = (double)last[k].PHIl[l] + reference_band_phase(&band, hz);
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
    /* The train is dropped from silence; the low cut's tracker, no part of it,
     * goes on from the state's (its corner under the fundamental in both, so
     * the samples are the same). */
    mbe_parms carried_train = carried, emptied_train = emptied;
    for (int i = MBE_LOW_CUT_TRACK_SLOT; i <= MBE_LOW_CUT_CORNER_SLOT; ++i) {
        carried_train.noiseOverlap[i] = 0.0f;
        emptied_train.noiseOverlap[i] = 0.0f;
    }
    assert(same_bits(carried_train.noiseOverlap, emptied_train.noiseOverlap, 96));
    assert(last.noiseOverlap[MBE_LOW_CUT_TRACK_SLOT] > 110.5f);
    assert(carried.noiseOverlap[MBE_LOW_CUT_TRACK_SLOT] > last.noiseOverlap[MBE_LOW_CUT_TRACK_SLOT]);
    assert(emptied.noiseOverlap[MBE_LOW_CUT_TRACK_SLOT] < last.noiseOverlap[MBE_LOW_CUT_TRACK_SLOT]);
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

/* The onset-columns mark changes nothing after a voiced frame (the mid-frame
 * model is that of an unmarked frame: the columns read every voiced model of
 * P25's rendering, and here every band the frame voices is voiced in the frame
 * before too and sounds), and a frame played as a P25 tone model keeps its own
 * model at mid-frame. */
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

/* ---- D-STAR's and P25's mid-frame voicing by 500 Hz columns ---- */

/* M(k) of the mid-frame rule from two frames' column samples s[0..8]. */
static int
rule_column(const int p[9], const int c[9], int k) {
    return (p[k] && c[k]) || ((p[k] != c[k]) && (p[k + 1] != c[k + 1]));
}

/* Bit `band` (0 for 0-1 kHz) of a D-STAR voicing word b1. */
static int
dstar_bit(int b1, int band) {
    return (b1 >> (3 - band)) & 1;
}

/* A D-STAR frame's voicing at its harmonic nearest angular frequency w, from
 * its word as the decoder voices a harmonic: its 1 kHz band's bit or that of
 * the band 250 Hz above it. */
static int
dstar_word_state(const mbe_parms* mp, int b1, double w) {
    const int n = (int)fmin(fmax(floor(w / (double)mp->w0 + 0.5), 1.0), (double)mp->L);
    const double hz = (double)n * (double)mp->w0 * 8000.0 / (2.0 * M_PI);
    const int own = (int)fmin(floor(away_from_whole(hz / 1000.0)), 3.0);
    const int above = (int)fmin(floor(away_from_whole((hz + 250.0) / 1000.0)), 3.0);
    return dstar_bit(b1, own) | dstar_bit(b1, above);
}

/* The mid-frame state of harmonic l of the mid-frame model of pitch w0
 * between a D-STAR frame of word `before` and one of word `now`, from the
 * rule's statement on the words: both columns of each 1 kHz band take its
 * bit; a frame with no voiced band after a voiced one holds by the columns
 * only when it lacks a pitch (`pitch_absent`), else is unvoiced; a band in the
 * same state in both frames keeps it; one that changes is voiced if and only
 * if M(k) or M(k'), k and k' the columns of its frequency and of 250 Hz more. */
static int
dstar_rule_state(const mbe_parms* prev, int before, const mbe_parms* cur, int now, int pitch_absent, double w0, int l) {
    const double w = w0 * l;
    const int a = dstar_word_state(prev, before, w), b = dstar_word_state(cur, now, w);
    if ((now == 0 && !pitch_absent) || a == b) {
        return b;
    }
    int p[9] = {0}, c[9] = {0};
    for (int k = 0; k < 8; ++k) {
        p[k] = dstar_bit(before, k / 2);
        c[k] = dstar_bit(now, k / 2);
    }
    const double hz = w * 8000.0 / (2.0 * M_PI);
    const int k = (int)fmin(floor(away_from_whole(hz / 500.0)), 7.0);
    const int k2 = (int)fmin(floor(away_from_whole((hz + 250.0) / 500.0)), 7.0);
    return rule_column(p, c, k) || rule_column(p, c, k2);
}

/* A D-STAR frame of pitch code b0 and word b1 as the decoder gives it. */
static void
dstar_frame(mbe_parms* mp, int b0, int b1) {
    mbe_parms prev, enhanced;
    char ambe_d[49];
    set_ambe2400(ambe_d, b0, b1);
    mbe_initMbeParms(mp, &prev, &enhanced);
    const int kind = mbe_decodeAmbe2400Parms(ambe_d, mp, &prev);
    assert(kind == 0);
    (void)kind;
}

/* The pitch codes of the D-STAR pairs: the same code (the highest pitch, 418
 * Hz at b0 0, where a band holds two or three harmonics, and lower ones), and
 * steps of pitch between the two frames. The second frame of each pair at a
 * code outside 64..119, where a frame with no voiced band lacks a pitch. */
static const int k_dstar_pairs[][2] = {{0, 0}, {1, 1},   {12, 12}, {40, 40}, {60, 60},  {123, 123},
                                       {0, 9}, {20, 34}, {50, 30}, {35, 60}, {123, 34}, {90, 120}};

/* Every pair of D-STAR words at each pair of pitch codes: the mid-frame
 * model's voicing is the rule's on the words (dstar_rule_state()), in D-STAR's
 * rendering; in the rate-33 rendering it is this frame's voicing, held from the
 * frame before through a frame lacking a pitch when that frame is voiced at its
 * harmonic at or next above 505 Hz, as before the columns. The words' pairs
 * include every one of the table in docs/dvsi-parity.md, an onset of each of
 * the fifteen voiced words, and each followed by a frame with no voiced band. */
static void
test_dstar_mid_columns_by_words(void) {
    for (size_t i = 0; i < sizeof(k_dstar_pairs) / sizeof(k_dstar_pairs[0]); ++i) {
        for (int before = 0; before < 16; ++before) {
            for (int now = 0; now < 16; ++now) {
                mbe_parms prev, cur;
                struct mbe_subframe_model mid, r33;
                dstar_frame(&prev, k_dstar_pairs[i][0], before);
                dstar_frame(&cur, k_dstar_pairs[i][1], now);
                const int pitch_absent = now == 0;
                assert((mbe_frame_mark(&cur) == MBE_FRAME_MARK_PITCH_ABSENT) == pitch_absent);
                mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
                mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, &r33);
                for (int l = 1; l <= mid.L; ++l) {
                    assert(mid.voiced[l] == dstar_rule_state(&prev, before, &cur, now, 1, (double)mid.w0, l));
                }
                const int held = pitch_absent && reference_holds_505(&prev);
                for (int l = 1; l <= r33.L; ++l) {
                    const double w = (double)r33.w0 * l;
                    assert(r33.voiced[l]
                           == (held ? dstar_word_state(&prev, before, w) : dstar_word_state(&cur, now, w)));
                }
            }
        }
    }
}

/* Whether every harmonic of the mid-frame model between [lo, hi) Hz has the
 * state `voiced`, and one at least lies there. */
static int
mid_voiced_in(const struct mbe_subframe_model* mid, double lo, double hi, int voiced) {
    int seen = 0, all = 1;
    for (int l = 1; l <= mid->L; ++l) {
        const double hz = (double)mid->w0 * l * 8000.0 / (2.0 * M_PI);
        if (hz >= lo && hz < hi) {
            seen = 1;
            all &= mid->voiced[l] == voiced;
        }
    }
    return seen && all;
}

/* D-STAR's mid-frame model between two words at pitch code b0. */
static void
dstar_mid(struct mbe_subframe_model* mid, int b0, int before, int now) {
    mbe_parms prev, cur;
    dstar_frame(&prev, b0, before);
    dstar_frame(&cur, b0, now);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, mid);
}

/* The cases of the table that the rule turns on, at b0 45 (212 Hz), each
 * frame's voicing reaching 250 Hz below its voiced bands as the decoder's
 * does: a column that changes under one that changes too is voiced at
 * mid-frame, and the top column of the run of changing columns is not,
 * whatever stands above it (1101 -> 1111 and 1111 -> 1101: column 5, under
 * column 6 voiced in both frames, starts late, or is let go, from 2500 to 2750
 * Hz, while column 4 below it starts early, or is held); a column whose upper
 * neighbour changes the other way is voiced (1000 -> 0100: column 1, where the
 * frame's own voicing ends, held to 750 Hz; columns 2 and 3 of the band
 * turning on start early and late). An onset starts the top column of each
 * voiced region late (0000 -> 1110: column 5); a frame lacking a pitch holds
 * the columns before it but the top one (1100 -> 0000: column 3 let go). */
static void
test_dstar_mid_columns_table(void) {
    struct mbe_subframe_model mid;
    dstar_mid(&mid, 45, 0xD, 0xF);
    assert(mid_voiced_in(&mid, 2000.0, 2500.0, 1) && mid_voiced_in(&mid, 2500.0, 2750.0, 0));
    assert(mid_voiced_in(&mid, 2750.0, 3700.0, 1));
    dstar_mid(&mid, 45, 0xF, 0xD);
    assert(mid_voiced_in(&mid, 2000.0, 2500.0, 1) && mid_voiced_in(&mid, 2500.0, 2750.0, 0));
    dstar_mid(&mid, 45, 0x8, 0x4);
    assert(mid_voiced_in(&mid, 0.0, 1500.0, 1) && mid_voiced_in(&mid, 1500.0, 2000.0, 0));
    dstar_mid(&mid, 45, 0x0, 0xE);
    assert(mid_voiced_in(&mid, 0.0, 2500.0, 1) && mid_voiced_in(&mid, 2500.0, 3000.0, 0));
    dstar_mid(&mid, 45, 0xC, 0x0);
    assert(mid_voiced_in(&mid, 0.0, 1500.0, 1) && mid_voiced_in(&mid, 1500.0, 2000.0, 0));
    /* A frame with no voiced band at a code DVSI's encoder keeps for such
     * frames (b0 100) holds nothing. */
    mbe_parms prev, cur;
    dstar_frame(&prev, 45, 0xF);
    dstar_frame(&cur, 100, 0x0);
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    assert(!mid.any_voiced && mid.L >= 1);
}

/* The band (1..K) of a P25 frame's harmonic n: three harmonics each, the last
 * band the rest. */
static int
p25_band(const mbe_parms* mp, int n) {
    const int band = (n + 2) / 3;
    return (band > mp->K) ? mp->K : band;
}

/* A P25 frame's voicing at harmonic n by its word (band 1 first; an empty
 * word for none), unvoiced past L. */
static int
p25_word_voiced(const mbe_parms* mp, const char* word, int n) {
    if (n < 1 || n > mp->L) {
        return 0;
    }
    const int band = p25_band(mp, n);
    return band >= 1 && (size_t)band <= strlen(word) && word[band - 1] == '1';
}

/* A P25 frame's column samples from its word: harmonic 1, then the harmonic
 * at or next above 505 + 500 (k - 1) Hz, unvoiced past L; in a frame of more
 * than 36 harmonics whose last band is voiced and band 11 is not, the last
 * band read as unvoiced. */
static void
p25_word_samples(const mbe_parms* mp, const char* word, int s[9]) {
    const double f0 = (double)mp->w0 * 8000.0 / (2.0 * M_PI);
    const int excepted = mp->L > 36 && strlen(word) == 12 && word[11] == '1' && word[10] != '1';
    memset(s, 0, 9 * sizeof(int));
    s[0] = p25_word_voiced(mp, word, 1);
    for (int k = 1; k < 8; ++k) {
        const int n = (int)ceil(away_from_whole((505.0 + 500.0 * (k - 1)) / f0));
        s[k] = p25_word_voiced(mp, word, n) && !(excepted && p25_band(mp, n) == 12);
    }
}

/* The mid-frame states of P25's rule on the words of two frames (an empty
 * word for a frame with no voiced band), the mid-frame model of pitch w0 and
 * L harmonics, into want[1..L]. */
static void
p25_rule_states(const mbe_parms* prev, const char* before, const mbe_parms* cur, const char* now, double w0, int L,
                int want[57]) {
    int p[9], c[9];
    p25_word_samples(prev, before, p);
    p25_word_samples(cur, now, c);
    for (int l = 1; l <= L; ++l) {
        const double w = w0 * l;
        const int na = (int)fmin(fmax(floor(w / (double)prev->w0 + 0.5), 1.0), (double)prev->L);
        const int nb = (int)fmin(fmax(floor(w / (double)cur->w0 + 0.5), 1.0), (double)cur->L);
        const int a = p25_word_voiced(prev, before, na), b = p25_word_voiced(cur, now, nb);
        const double hz = w * 8000.0 / (2.0 * M_PI);
        int k = 0, k2 = 0;
        for (int i = 0; i < 7; ++i) {
            k += hz >= 505.0 + 500.0 * i;
            k2 += hz + 250.0 >= 505.0 + 500.0 * i;
        }
        want[l] = (a == b) ? b : (rule_column(p, c, k) || rule_column(p, c, k2));
    }
}

/* A P25 frame of pitch code b0 and voicing word `word` ("" or all '0': none). */
static void
p25_frame(mbe_parms* mp, int b0, const char* word) {
    mbe_parms prev, enhanced;
    char imbe_d[88] = {0};
    if (strchr(word, '1') != NULL) {
        set_imbe_voicing(imbe_d, b0, word);
    } else {
        set_imbe_b0(imbe_d, b0);
    }
    mbe_initMbeParms(mp, &prev, &enhanced);
    const int kind = mbe_decodeImbe4400Parms(imbe_d, mp, &prev);
    assert(kind == 0);
    (void)kind;
}

/* Two P25 frames: pitch codes and words (an empty word: no voiced band), and
 * whether the mid-frame model has a voiced band. */
static const struct {
    int b0[2];
    const char* word[2];
    int any;
} k_p25_pairs[] = {
    /* Onsets: column 5 covered, its harmonics and those 250 Hz below it
     * sound; a frame voiced to 502.8 Hz covers nothing; the isolated last band
     * of a short frame (L 9, word 001, a single tone's at b0 2) covers column 5. */
    {{150, 86}, {"", "1000001100"}, 1},
    {{150, 56}, {"", "10000000"}, 0},
    {{150, 2}, {"", "001"}, 1},
    /* P25's exception: an isolated voiced last band in a frame of 37
     * harmonics (b0 120) reads unvoiced in the columns, so that column 6 is not
     * covered; not with band 11 voiced too, nor in a frame of 36 (b0 116). */
    {{150, 120}, {"", "111111111101"}, 1},
    {{150, 120}, {"", "111111111111"}, 1},
    {{150, 116}, {"", "111111111101"}, 1},
    {{120, 120}, {"100000000000", "100000000001"}, 1},
    {{116, 116}, {"100000000000", "100000000001"}, 1},
    /* Frames lacking a pitch (b0 25): a frame voiced from column 2 (its 505 Hz
     * harmonic unvoiced) holds its covered columns; one voiced in isolated
     * bands only covers none and holds nothing. */
    {{86, 25}, {"0011110000", ""}, 1},
    {{97, 25}, {"10000000000", ""}, 0},
    {{76, 25}, {"100010000", ""}, 0},
    {{120, 25}, {"111111111101", ""}, 1},
    /* Inside voiced runs, both ways, with pitch steps. */
    {{86, 86}, {"1111110000", "1111111100"}, 1},
    {{86, 86}, {"1111111100", "1111110000"}, 1},
    {{60, 80}, {"11111100", "111011110"}, 1},
    {{100, 40}, {"10101010101", "011101"}, 1},
};

/* P25's mid-frame voicing by columns, against the rule's statement on the
 * voicing words (p25_rule_states()); in the rate-33 rendering this frame's
 * voicing, held from the frame before through a frame lacking a pitch when
 * that frame is voiced at its harmonic at or next above 505 Hz, as before
 * the columns. */
static void
test_p25_mid_columns_by_words(void) {
    for (size_t i = 0; i < sizeof(k_p25_pairs) / sizeof(k_p25_pairs[0]); ++i) {
        mbe_parms prev, cur;
        struct mbe_subframe_model mid, r33;
        int want[57];
        p25_frame(&prev, k_p25_pairs[i].b0[0], k_p25_pairs[i].word[0]);
        p25_frame(&cur, k_p25_pairs[i].b0[1], k_p25_pairs[i].word[1]);
        mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
        p25_rule_states(&prev, k_p25_pairs[i].word[0], &cur, k_p25_pairs[i].word[1], (double)mid.w0, mid.L, want);
        for (int l = 1; l <= mid.L; ++l) {
            assert(mid.voiced[l] == want[l]);
        }
        assert(mid.any_voiced == k_p25_pairs[i].any);
        mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, &r33);
        const int held = !reference_voiced(&cur) && mbe_frame_mark(&cur) == MBE_FRAME_MARK_PITCH_ABSENT
                         && reference_holds_505(&prev);
        const mbe_parms* source = held ? &prev : &cur;
        const char* word = held ? k_p25_pairs[i].word[0] : k_p25_pairs[i].word[1];
        for (int l = 1; l <= r33.L; ++l) {
            const int m = (int)fmin(fmax(floor((double)r33.w0 * l / (double)source->w0 + 0.5), 1.0), (double)source->L);
            assert(r33.voiced[l] == p25_word_voiced(source, word, m));
        }
    }
}

/* The exception's harmonics, and where the rule does not apply in P25's
 * rendering. At an onset of word 111111111101 (band 11 unvoiced) a frame of 37
 * harmonics (b0 120, 100.3 Hz) reads its last band as unvoiced: column 6 is not
 * covered, so harmonic 29 (2909 Hz, column 5) sounds at mid-frame and 30 to
 * 34 (from 3009 Hz) do not; with band 11 voiced too it reads its last band,
 * column 6 is covered, and harmonics 30 to 34 sound; a frame of 36 (b0 116,
 * 103.1 Hz) reads it, and harmonic 30 (3093 Hz) sounds. A frame after a tone
 * model takes the rule as after any model (the mark of the model before is
 * caller state, and the columns' rule reads no mark of that model), and a frame
 * with no voiced band that does not lack a pitch holds nothing. */
static void
test_p25_mid_columns_cases(void) {
    mbe_parms prev, cur, tone;
    struct mbe_subframe_model mid, own;
    p25_frame(&prev, 150, "");
    p25_frame(&cur, 120, "111111111101");
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
    assert(cur.L == 37 && mid.voiced[29] == 1 && mid.voiced[30] == 0 && mid.voiced[34] == 0);
    p25_frame(&cur, 120, "111111111111");
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
    assert(mid.voiced[30] == 1 && mid.voiced[34] == 1);
    p25_frame(&cur, 116, "111111111101");
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
    assert(cur.L == 36 && mid.voiced[30] == 1);

    /* The frame before marked as a tone model is read in the columns as any
     * model: a voiced frame after it, and a frame lacking a pitch (held by the
     * columns the model covers), take the mid-frame model of the same pair
     * with the mark not there, and the rule's voicing. */
    p25_frame(&tone, 86, "1111110000");
    mbe_parms plain = tone;
    mbe_set_frame_mark(&tone, MBE_FRAME_MARK_TONE);

    static const struct {
        int b0;
        const char* word;
    } after[3] = {{86, "1111111100"}, {25, ""}, {150, ""}};

    for (int i = 0; i < 3; ++i) {
        int want[57];
        p25_frame(&cur, after[i].b0, after[i].word);
        mbe_subframe_mid_model(&tone, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
        mbe_subframe_mid_model(&plain, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &own);
        assert(same_model(&mid, &own));
        p25_rule_states(&plain, "1111110000", &cur, after[i].word, (double)mid.w0, mid.L, want);
        const int by_columns = i < 2; /* b0 150 is a code kept for frames with no voiced band */
        for (int l = 1; l <= mid.L; ++l) {
            assert(!by_columns || mid.voiced[l] == want[l]);
        }
        assert(mid.any_voiced == by_columns);
    }
}

/* The positions D-STAR's columns are read at, l 16 f0 with f0 the stored
 * pitch over 2 pi in float (as the synthesis forms it), lie far from a whole
 * number and from a half at every D-STAR pitch code the decoder plays voiced
 * (b0 0..123; b0 124's silence model has no voiced band, and the columns are
 * not read at its pitch: next to a voiced frame the mid-frame pitch is that
 * frame's), at least 6.0e-5 of a column from a whole number (6.58e-5, b0 80,
 * harmonic 4) and 1.5e-5 from a half (b0 122, harmonic 49: 2.00e-5 in IEEE
 * builds, 1.81e-5 where the decoder's stored pitch rounds an ulp otherwise,
 * with fast-math and in 32-bit x87 builds), and give
 * the columns the decoder's own f0 (mbe_ambe2400_f0()) gives, from which the
 * frame's voicing was set; so no build's arithmetic moves a harmonic to
 * another column, and the model's voicing and the columns' reading cannot
 * round apart. The pitches are the decoder's, through its own path. */
static void
test_dstar_column_edges(void) {
    double least_whole = 1.0, least_half = 1.0;
    for (int b0 = 0; b0 <= 123; ++b0) {
        mbe_parms mp;
        dstar_frame(&mp, b0, 0xF);
        const volatile float f0 = mp.w0 / (2.0f * (float)M_PI);
        const float decoder_f0 = mbe_ambe2400_f0(b0);
        for (int l = 1; l <= mp.L; ++l) {
            const volatile float position = (float)l * 16.0f * f0;
            const volatile float decoder_position = (float)l * 16.0f * decoder_f0;
            const volatile float above = position + 0.5f;
            const volatile float decoder_above = decoder_position + 0.5f;
            const int own = (int)position, decoder_own = (int)decoder_position;
            const int reach = (int)above, decoder_reach = (int)decoder_above;
            assert(own == decoder_own && reach == decoder_reach);
            const double x = (double)position;
            least_whole = fmin(least_whole, fabs(x - floor(x + 0.5)));
            least_half = fmin(least_half, fabs(x + 0.5 - floor(x + 1.0)));
        }
    }
    assert(least_whole > 6.0e-5 && least_whole < 7.0e-5);
    assert(least_half > 1.5e-5 && least_half < 2.2e-5);
}

/* A stream of decoded models, each rendered (the voiced part, at unit gain)
 * and, from the same state, advanced unrendered: the two leave the same pulse
 * train and mark after every frame. The pulse phase after each frame goes to
 * phase[]. */
static void
render_and_advance(const mbe_parms* frames, int count, enum mbe_voiced_synthesis rendering, double phase[]) {
    mbe_parms rendered, advanced;
    float out[160];
    fresh_stream(&rendered);
    advanced = rendered;
    for (int f = 0; f < count; ++f) {
        mbe_parms next = with_model(&advanced, &frames[f]);
        mbe_set_frame_mark(&next, mbe_frame_mark(&frames[f]));
        mbe_advance_subframes(&next, &advanced, 0, rendering);
        advanced = next;
        next_frame_as(out, &rendered, &frames[f], 0, rendering);
        assert(all_finite(out, 160));
        assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
        assert(same_bits(&rendered.PHIl[0], &advanced.PHIl[0], 1));
        phase[f] = pulse_phase(&rendered);
    }
}

static double
hz_of(const mbe_parms* mp) {
    return (double)mp->w0 * 8000.0 / (2.0 * M_PI);
}

/* A frame voiced only above its first band, then a frame with no voiced band
 * lacking a pitch, then a voiced frame (D-STAR: b1 0100 at b0 123, b0 34 with
 * b1 0000, b1 0111 at b0 115, as DVSI's vectors hold them before the onsets
 * that showed the hold): the voiced frame's column 2 holds through the first
 * half of the frame lacking a pitch, so the pulse train runs on through it at
 * the voiced frame's pitch and the frame after continues it. With the frame
 * between at a code DVSI's encoder keeps for frames with no voiced band
 * (b0 80), nothing holds and the train starts afresh. P25 alike: a frame
 * voiced from 1005 Hz (word 0011110000 at b0 86) holds its covered columns
 * through a frame lacking a pitch (b0 25); one voiced in isolated bands only
 * (word 10000000000 at b0 97) covers none, and the train starts afresh, as
 * DVSI's P25 decoder starts it after such frames. Rendering and advancing
 * leave the same train throughout. */
static void
test_train_through_held_columns(void) {
    double phase[4];
    mbe_parms frames[4];
    for (int cut = 0; cut < 2; ++cut) {
        dstar_frame(&frames[0], 123, 0x4);
        frames[1] = frames[0];
        dstar_frame(&frames[2], cut ? 80 : 34, 0x0);
        dstar_frame(&frames[3], 115, 0x7);
        render_and_advance(frames, 4, MBE_VOICED_SYNTHESIS_DVSI, phase);
        const double through = cut ? 0.0 : phase[1] + cycles_of(hz_of(&frames[1]), 160.0);
        assert(fabs(cycles_apart(phase[2], through)) < 1e-4);
        assert(fabs(cycles_apart(phase[3], through + cycles_of(hz_of(&frames[3]), 160.0))) < 1e-4);
        assert(cut || fabs(cycles_apart(phase[2], 0.0)) > 0.01);
    }

    static const struct {
        int b0;
        const char* word;
        int held;
    } p25[2] = {{86, "0011110000", 1}, {97, "10000000000", 0}};

    for (int i = 0; i < 2; ++i) {
        p25_frame(&frames[0], p25[i].b0, p25[i].word);
        frames[1] = frames[0];
        p25_frame(&frames[2], 25, "");
        p25_frame(&frames[3], 86, "1111111100");
        render_and_advance(frames, 4, MBE_VOICED_SYNTHESIS_DVSI_IMBE, phase);
        const double through = p25[i].held ? phase[1] + cycles_of(hz_of(&frames[1]), 160.0) : 0.0;
        assert(fabs(cycles_apart(phase[2], through)) < 1e-4);
        assert(fabs(cycles_apart(phase[3], through + cycles_of(hz_of(&frames[3]), 160.0))) < 1e-4);
        assert(!p25[i].held || fabs(cycles_apart(phase[2], 0.0)) > 0.01);
    }
}

/* The mid-frame model's band nearest harmonic l of a 151 Hz frame, between
 * `prev` and `cur` (both at 151 Hz), the way `rendering` says. */
static int
mid_band(const mbe_parms* prev, const mbe_parms* cur, enum mbe_voiced_synthesis rendering, int l) {
    struct mbe_subframe_model mid;
    mbe_subframe_mid_model(prev, cur, 0, rendering, &mid);
    assert(mid.L == cur->L);
    return mid.voiced[l];
}

/* A band that plays nothing (BAND_NONE: a caller's voicing that is neither 0
 * nor 1) in either frame plays nothing at mid-frame where the columns decide
 * (the rule has no state to compare), so that wherever the frame whose voicing
 * stood at mid-frame before had such a band, it still plays nothing: inside a
 * voiced run, at an onset (the frame before, with no voiced band, holding one)
 * and through a frame lacking a pitch. In the rate-33 rendering this frame's
 * band stands, as before. Pitches that are not positive and normal, and
 * harmonic counts outside 1..56, give a mid-frame model with no voiced band. */
static void
test_mid_columns_bands_with_no_state(void) {
    static const enum mbe_voiced_synthesis renderings[2] = {MBE_VOICED_SYNTHESIS_DVSI, MBE_VOICED_SYNTHESIS_DVSI_IMBE};

    for (int r = 0; r < 2; ++r) {
        for (int which = 0; which < 3; ++which) {
            mbe_parms prev, cur;
            set_model(&prev, 151.0, 24, 24, 16.0f);
            set_model(&cur, 151.0, 24, 12, 16.0f);
            prev.Vl[18] = (which != 1) ? 7 : 1;
            cur.Vl[18] = (which != 0) ? -1 : 0;
            assert(mid_band(&prev, &cur, renderings[r], 18) == 2);
            assert(mid_band(&prev, &cur, renderings[r], 17) != 2);
            assert(mid_band(&prev, &cur, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, 18) == ((which != 0) ? 2 : 0));
        }
        mbe_parms prev, cur;
        set_model(&prev, 151.0, 24, 0, 16.0f);
        prev.Vl[5] = 1000000;
        set_model(&cur, 151.0, 24, 24, 16.0f);
        assert(mid_band(&prev, &cur, renderings[r], 5) == 2 && mid_band(&prev, &cur, renderings[r], 4) == 1);
        set_model(&prev, 151.0, 24, 24, 16.0f);
        set_model(&cur, 151.0, 24, 0, 16.0f);
        cur.Vl[3] = 2;
        mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PITCH_ABSENT);
        assert(mid_band(&prev, &cur, renderings[r], 3) == 2 && mid_band(&prev, &cur, renderings[r], 2) == 1);
    }
    static const uint32_t w0_bits[] = {0x7FC00000u, 0x7F800000u, 0x00000000u, 0x00000001u, 0x007FFFFFu, 0xBF000000u};
    static const int counts[] = {0, -1, 57, 1000};
    for (int r = 0; r < 3; ++r) {
        static const enum mbe_voiced_synthesis all[3] = {MBE_VOICED_SYNTHESIS_DVSI, MBE_VOICED_SYNTHESIS_DVSI_IMBE,
                                                         MBE_VOICED_SYNTHESIS_DVSI_AMBE2450};

        for (int which = 0; which < 2; ++which) {
            for (size_t i = 0; i < sizeof(w0_bits) / sizeof(w0_bits[0]); ++i) {
                mbe_parms frames[2];
                struct mbe_subframe_model mid;
                set_model(&frames[0], 151.0, 24, 24, 16.0f);
                set_model(&frames[1], 183.0, 20, 6, 16.0f);
                memcpy(&frames[which].w0, &w0_bits[i], sizeof(float));
                mbe_subframe_mid_model(&frames[0], &frames[1], 0, all[r], &mid);
                assert(!mid.any_voiced);
            }
            for (size_t i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i) {
                mbe_parms frames[2];
                struct mbe_subframe_model mid;
                set_model(&frames[0], 151.0, 24, 24, 16.0f);
                set_model(&frames[1], 183.0, 20, 6, 16.0f);
                frames[which].L = counts[i];
                mbe_subframe_mid_model(&frames[0], &frames[1], 0, all[r], &mid);
                assert(!mid.any_voiced && mid.L == 0);
            }
        }
    }
}

/* Every mark (and a stored value that is none) in the frame before, in this
 * frame or in both, in each subframe rendering, over pairs of frames the
 * columns read (a voiced run with bands turning on and off and a pitch step;
 * an onset; a frame with no voiced band after a voiced one; a P25-like frame of
 * 37 harmonics with an isolated last band), at pitches from 101.3 to 183 Hz
 * (pitches at the edges of what a float holds are test_untrusted_mid_voicing_marks'
 * and test_mid_columns_bands_with_no_state's): the mid-frame voicing is the reference's
 * (reference_mid_voicing(); a frame marked as a tone model or as holding the
 * previous model takes a frame's own model), and the frame renders finite
 * samples and state, rendered and advanced alike (the sanitizer build runs
 * this). */
static void
test_marks_with_columns(void) {
    static const float marks[] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 80.0f, 5.5f};

    static const struct {
        double hz[2];
        int L[2], voiced[2];
    } pairs[] = {{{151.0, 183.0}, {24, 20}, {20, 9}},
                 {{183.0, 151.0}, {20, 24}, {9, 20}},
                 {{151.0, 151.0}, {24, 24}, {0, 17}},
                 {{151.0, 151.0}, {24, 24}, {13, 0}},
                 {{101.3, 101.3}, {37, 37}, {6, 37}}};

    static const enum mbe_voiced_synthesis renderings[3] = {MBE_VOICED_SYNTHESIS_DVSI, MBE_VOICED_SYNTHESIS_DVSI_IMBE,
                                                            MBE_VOICED_SYNTHESIS_DVSI_AMBE2450};

    for (size_t p = 0; p < sizeof(pairs) / sizeof(pairs[0]); ++p) {
        for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
            for (int where = 0; where < 3; ++where) {
                for (int r = 0; r < 3; ++r) {
                    mbe_parms last, model;
                    float out[160];
                    fresh_stream(&last);
                    set_model(&model, pairs[p].hz[0], pairs[p].L[0], pairs[p].voiced[0], 16.0f);
                    next_frame_as(out, &last, &model, 0, renderings[r]);
                    set_model(&model, pairs[p].hz[1], pairs[p].L[1], pairs[p].voiced[1], 16.0f);
                    if (p == 4) {
                        model.Vl[33] = 0;
                    }
                    if (where != 1) {
                        memcpy(&last.PHIl[0], &marks[m], sizeof(float));
                    }
                    if (where != 0) {
                        memcpy(&model.PHIl[0], &marks[m], sizeof(float));
                    }
                    mbe_parms cur = with_model(&last, &model);
                    memcpy(&cur.PHIl[0], &model.PHIl[0], sizeof(float));
                    struct mbe_subframe_model mid, own;
                    const enum mbe_frame_mark mark = mbe_frame_mark(&cur);
                    /* A late onset as the synthesis passes it: in the rate-33
                     * rendering whatever the frame before, in the others after a
                     * frame with no voiced band. */
                    const int late =
                        mark == MBE_FRAME_MARK_LATE_ONSET
                        && (renderings[r] == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 || !reference_voiced(&last));
                    mbe_subframe_mid_model(&last, &cur, late, renderings[r], &mid);
                    if (mark == MBE_FRAME_MARK_TONE || mark == MBE_FRAME_MARK_PREVIOUS_HELD) {
                        mbe_subframe_own_model((mark == MBE_FRAME_MARK_TONE) ? &cur : &last, &own);
                        assert(same_model(&mid, &own));
                    } else {
                        unsigned char want[57] = {0};
                        reference_mid_voicing(&last, &cur, (double)mid.w0, mid.L, renderings[r], want);
                        for (int l = 1; l <= mid.L; ++l) {
                            assert(mid.voiced[l] == want[l]);
                        }
                    }
                    mbe_parms advanced = cur;
                    mbe_advance_subframes(&advanced, &last, 0, renderings[r]);
                    mbe_parms rendered = cur;
                    memset(out, 0, sizeof(out));
                    mbe_synthesize_subframes(out, &rendered, &last, NULL, 1.0f, 0, renderings[r]);
                    assert(all_finite(out, 160) && all_finite(rendered.noiseOverlap, 96));
                    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
                    assert(same_bits(&rendered.PHIl[0], &advanced.PHIl[0], 1));
                }
            }
        }
    }
}

/* One frame of a codec's data path: 0 D-STAR, 1 P25, 2 rate 33. */
static void
codec_frame(int codec, const char* bits, float out[160], mbe_parms* cur, mbe_parms* prev, mbe_parms* enhanced) {
    mbe_process_result result;
    mbe_initProcessResult(&result);
    int rc;
    if (codec == 0) {
        rc = mbe_processAmbe2400Dataf(out, &result, bits, cur, prev, enhanced);
    } else if (codec == 1) {
        rc = mbe_processImbe4400Dataf(out, &result, bits, cur, prev, enhanced);
    } else {
        rc = mbe_processAmbe2450Dataf(out, &result, bits, cur, prev, enhanced);
    }
    assert(rc >= 0);
    (void)rc;
    assert(all_finite(out, 160));
}

/* A rate-33 voice frame of pitch code b0 and V/UV row `row`. */
static void
set_ambe2450(char ambe_d[49], int b0, int row) {
    int b[9] = {b0, row, 20, 100, 30, 10, 5, 5, 3};
    mbe_ambe2450_pack(b, ambe_d);
}

/* Every mark (and a stored value that is none) in each of a stream's three
 * states (cur_mp, prev_mp, prev_mp_enhanced), through each codec's data path:
 * after voiced frames, a frame with no voiced band lacking a pitch (D-STAR and
 * P25) or of V/UV row 16 (rate 33), then a voiced frame, play finite samples,
 * and play exactly as the stream with no mark there, in every slot. Here no
 * mark can act: the decoder replaces cur_mp's mark, prev_mp's is read only by a
 * repeat, and the last frame played is voiced, so its mark of an unvoiced frame
 * (PITCH_ABSENT, PITCH_CARRIES, PREVIOUS_HELD) does not change how the train
 * loads or the mid-frame pitch (mid_pitch() reads a previous PITCH_CARRIES,
 * the train's loading a previous PITCH_ABSENT or PREVIOUS_HELD, and in rate
 * 33's rendering MID_COLUMN_0, all only after a frame with no voiced band;
 * test_row17_mark_in_codec_states() forges 6 after one), and the columns' rule reads no mark of
 * the model before (a forged MBE_FRAME_MARK_TONE there changed P25's columns
 * before the rule was applied after a tone model). */
static void
test_marks_in_codec_states(void) {
    static const float marks[] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 80.0f, 5.5f};
    char voiced[3][88], unvoiced[3][88], after[3][88];
    set_ambe2400(voiced[0], 45, 0xC);
    set_ambe2400(unvoiced[0], 34, 0x0);
    set_ambe2400(after[0], 50, 0xF);
    set_imbe_voicing(voiced[1], 86, "0011110000");
    memset(unvoiced[1], 0, 88);
    set_imbe_b0(unvoiced[1], 25);
    set_imbe_voicing(after[1], 86, "1111111100");
    set_ambe2450(voiced[2], 40, 8);
    set_ambe2450(unvoiced[2], 40, 16);
    set_ambe2450(after[2], 45, 0);
    for (int codec = 0; codec < 3; ++codec) {
        for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
            for (int state = 0; state < 3; ++state) {
                mbe_parms s[2][3];
                float out[2][320];
                mbe_initMbeParms(&s[0][0], &s[0][1], &s[0][2]);
                for (int f = 0; f < 3; ++f) {
                    codec_frame(codec, voiced[codec], out[0], &s[0][0], &s[0][1], &s[0][2]);
                }
                memcpy(s[1], s[0], sizeof(s[1]));
                memcpy(&s[1][state].PHIl[0], &marks[m], sizeof(float));
                for (int k = 0; k < 2; ++k) {
                    codec_frame(codec, unvoiced[codec], out[k], &s[k][0], &s[k][1], &s[k][2]);
                    codec_frame(codec, after[codec], out[k] + 160, &s[k][0], &s[k][1], &s[k][2]);
                }
                assert(same_bits(out[0], out[1], 320));
            }
        }
    }
}

/* A P25 pitch step whose mid-frame harmonic lies just under halfway between
 * two harmonics of the frame before: b0 134, word 000000000010, then b0 161,
 * word 111111111101. The mid-frame model's harmonic 36 lies at 33.4999986 of
 * the first frame's harmonics (nearest 33, voiced) and 38.69 of the second's
 * (39, voiced): voiced in both, it stays voiced. (A float quotient rounds the
 * first position to 33.5 and took harmonic 34, unvoiced, for it, and the
 * columns silenced the band.) */
static void
test_p25_nearest_harmonic_across_a_step(void) {
    mbe_parms prev, cur;
    struct mbe_subframe_model mid;
    int want[57];
    p25_frame(&prev, 134, "000000000010");
    p25_frame(&cur, 161, "111111111101");
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
    const double position = (double)mid.w0 * 36.0 / (double)prev.w0;
    assert(position > 33.4999 && position < 33.5);
    assert(prev.Vl[33] == 1 && prev.Vl[34] == 0 && cur.Vl[39] == 1);
    assert(mid.voiced[36] == 1);
    p25_rule_states(&prev, "000000000010", &cur, "111111111101", (double)mid.w0, mid.L, want);
    for (int l = 1; l <= mid.L; ++l) {
        assert(mid.voiced[l] == want[l]);
    }
}

/* The voicing words of the P25 sweep for a frame of K bands (at most 12
 * characters and the terminator). */
static void
p25_sweep_word(int K, int pattern, char word[13]) {
    for (int band = 0; band < K; ++band) {
        int voiced = 0;
        switch (pattern) {
            case 0: voiced = 1; break;                  /* every band */
            case 1: voiced = band < (K + 1) / 2; break; /* the lower half */
            case 2: voiced = (band % 2) == 0; break;    /* alternate bands */
            default: voiced = (band % 3) != 1 || band == K - 1; break;
        }
        word[band] = voiced ? '1' : '0';
    }
    word[K] = '\0';
}

/* P25's columns across pitch steps, against the rule's statement on the
 * voicing words (p25_rule_states(): each frame's nearest harmonic in double,
 * ties upward): every pair of a code among every third from 0, with 134 and 161
 * (72 codes), before one among every fourth from 1, with 134 (53 codes), the 19
 * pairs of a code with itself left out (3,797 pairs, steps of up to 205 codes
 * up and 206 down), each with four voicing words in each frame (16 pairs of words:
 * 60,752 frame pairs; about a second in a Debug build), the frames decoded by
 * the P25 decoder. Among all ordered pairs of different P25 codes, only 134 and
 * 161, in either order, put a mid-frame harmonic where the float quotient and
 * the double one take different nearest harmonics (harmonic 36:
 * test_p25_nearest_harmonic_across_a_step), so the sweep holds both. The
 * voicing words are set from each band's parameter bit, found once per code. */
static void
test_p25_mid_columns_across_steps(void) {
    static int band_bits[208][13];
    static mbe_parms frames[208][4];
    static char words[208][4][13];
    static int ready[208];
    long pairs = 0;
    for (int b0 = 0; b0 < 208; ++b0) {
        if (b0 % 3 != 0 && b0 % 4 != 1 && b0 != 134) {
            continue;
        }
        mbe_parms probe;
        p25_frame(&probe, b0, "");
        for (int band = 1; band <= probe.K; ++band) {
            band_bits[b0][band] = imbe_band_bit(b0, band);
            assert(band_bits[b0][band] >= 0);
        }
        for (int w = 0; w < 4; ++w) {
            char imbe_d[88] = {0};
            mbe_parms prev, enhanced;
            p25_sweep_word(probe.K, w, words[b0][w]);
            set_imbe_b0(imbe_d, b0);
            for (int band = 1; band <= probe.K; ++band) {
                imbe_d[band_bits[b0][band]] = (char)(words[b0][w][band - 1] == '1');
            }
            mbe_initMbeParms(&frames[b0][w], &prev, &enhanced);
            const int kind = mbe_decodeImbe4400Parms(imbe_d, &frames[b0][w], &prev);
            assert(kind == 0);
            (void)kind;
        }
        ready[b0] = 1;
    }
    for (int a = 0; a < 208; ++a) {
        for (int b = 0; b < 208; ++b) {
            const int first = (a % 3 == 0) || a == 134 || a == 161;
            const int second = (b % 4 == 1) || b == 134;
            if (a == b || !first || !second) {
                continue;
            }
            assert(ready[a] && ready[b]);
            ++pairs;
            for (int wa = 0; wa < 4; ++wa) {
                for (int wb = 0; wb < 4; ++wb) {
                    struct mbe_subframe_model mid;
                    int want[57];
                    const mbe_parms* prev = &frames[a][wa];
                    const mbe_parms* cur = &frames[b][wb];
                    mbe_subframe_mid_model(prev, cur, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
                    p25_rule_states(prev, words[a][wa], cur, words[b][wb], (double)mid.w0, mid.L, want);
                    for (int l = 1; l <= mid.L; ++l) {
                        assert(mid.voiced[l] == want[l]);
                    }
                }
            }
        }
    }
    assert(pairs == 3797);
}

/* The state of a stream (all three states and the frame's samples) as bits. */
static int
same_stream_state(const mbe_parms* a, const mbe_parms* b) {
    int same = same_bits(&a->w0, &b->w0, 1) && a->L == b->L && a->K == b->K && a->swn == b->swn
               && a->tonePhase == b->tonePhase && a->amplitudeThreshold == b->amplitudeThreshold
               && a->errorCountTotal == b->errorCountTotal && a->errorCount4 == b->errorCount4
               && a->repeatCount == b->repeatCount;
    for (int l = 0; l < 57; ++l) {
        same &= a->Vl[l] == b->Vl[l];
    }
    same &= same_bits(a->Ml, b->Ml, 57) && same_bits(a->log2Ml, b->log2Ml, 57) && same_bits(a->PHIl, b->PHIl, 57);
    same &= same_bits(a->PSIl, b->PSIl, 57) && same_bits(&a->gamma, &b->gamma, 1)
            && same_bits(&a->localEnergy, &b->localEnergy, 1) && same_bits(&a->errorRate, &b->errorRate, 1)
            && same_bits(&a->mutingThreshold, &b->mutingThreshold, 1);
    same &= same_bits(a->previousUw, b->previousUw, 256) && same_bits(&a->noiseSeed, &b->noiseSeed, 1)
            && same_bits(a->noiseOverlap, b->noiseOverlap, 96);
    return same;
}

/* A Golay (23,12) word whose data bits 11..22 are set in `word`: the parity
 * bits 0..10 for which the decoder corrects nothing and gives the data back. */
static void
golay_parity(char word[23]) {
    for (int parity = 0; parity < 2048; ++parity) {
        char out[23];
        for (int j = 0; j < 11; ++j) {
            word[j] = (char)((parity >> j) & 1);
        }
        if (mbe_golay2312(word, out) == 0) {
            int same = 1;
            for (int j = 11; j < 23; ++j) {
                same &= out[j] == word[j];
            }
            if (same) {
                return;
            }
        }
    }
    assert(0);
}

/* A ProVoice (15,11) Hamming word with data bits 4..14 set: its parity bits. */
static void
hamming7100_parity(char word[15]) {
    for (int parity = 0; parity < 16; ++parity) {
        char out[15];
        for (int j = 0; j < 4; ++j) {
            word[j] = (char)((parity >> j) & 1);
        }
        if (mbe_7100x4400hamming1511(word, out) == 0) {
            return;
        }
    }
    assert(0);
}

/* The ProVoice (IMBE 7100x4400) layout of parameter bits given in the 7200
 * layout the ProVoice decoder converts to (mbe_convertImbe7100to7200(), whose
 * permutation this runs backwards). */
static void
imbe7100_layout(const char t[88], char d[88]) {
    const int b0 =
        (t[0] << 7) | (t[1] << 6) | (t[2] << 5) | (t[3] << 4) | (t[4] << 3) | (t[5] << 2) | (t[85] << 1) | t[86];
    const float w0 = (float)(4 * M_PI) / (float)((float)b0 + 39.5);
    const int L = (int)(0.9254 * (int)((M_PI / w0) + 0.25));
    const int K = (L < 37) ? (int)((float)(L + 2) / (float)3) : 12;
    memset(d, 0, 88);
    d[0] = t[87];
    d[42] = t[48 + K];
    d[43] = t[49 + K];
    for (int i = 0; i < K; ++i) {
        d[44 + i] = t[48 + i];
    }
    int j = 0, k = 1;
    while (j < 87) {
        d[k] = t[j];
        if (++j == 48) {
            j += K + 2;
        }
        if (++k == 42) {
            k += K + 2;
        }
    }
}

/* A ProVoice frame (7 x 24 bits, as received) that decodes, with no error, to
 * parameter bits `t` in the 7200 layout: the 7100 layout's fields placed in the
 * code words as mbe_decodeImbe7100x4400Frame() reads them, their parity bits
 * found by the decoders themselves, and the code words C1 to C5 modulated with
 * the sequence seeded from C0's data, which the decoder takes off again. */
static void
pack_imbe7100(const char t[88], char fr[7][24]) {
    char d[88], word[23], hword[15];
    imbe7100_layout(t, d);
    memset(fr, 0, sizeof(char[7][24]));
    memset(word, 0, sizeof(word));
    for (int j = 18; j > 11; --j) {
        word[j - 1] = d[18 - j];
    }
    golay_parity(word);
    for (int j = 0; j < 18; ++j) {
        fr[0][j + 1] = word[j];
    }
    int n = 7;
    for (int i = 1; i < 4; ++i) {
        memset(word, 0, sizeof(word));
        for (int j = 22; j > 10; --j) {
            word[j] = d[n++];
        }
        golay_parity(word);
        for (int j = 0; j < 23; ++j) {
            fr[i][(i == 1) ? j + 1 : j] = word[j];
        }
    }
    for (int i = 4; i < 6; ++i) {
        memset(hword, 0, sizeof(hword));
        for (int j = 14; j >= 4; --j) {
            hword[j] = d[n++];
        }
        hamming7100_parity(hword);
        for (int j = 0; j < 15; ++j) {
            fr[i][j] = hword[j];
        }
    }
    for (int j = 22; j >= 0; --j) {
        fr[6][j] = d[n++];
    }
    assert(n == 88);
    /* The modulation is an exclusive or, so modulating is demodulating. */
    const int rc = mbe_demodulateImbe7100x4400Data(fr);
    assert(rc == 0);
    (void)rc;
    char back[88];
    mbe_process_result result;
    const int errors = mbe_decodeImbe7100x4400Frame((const char (*)[24])fr, back, &result);
    assert(errors == 0 && (result.flags & MBE_PROCESS_FLAG_PROVOICE) != 0u);
    (void)errors;
    for (int i = 0; i < 88; ++i) {
        assert(back[i] == t[i]);
    }
}

/* A tone mark forged in a caller's P25 state changes nothing: three frames of
 * b0 86, word 0011110000, then a frame with no voiced band at b0 25 (lacking a
 * pitch), through mbe_processImbe4400Dataf() (P25, and ProVoice with its flag)
 * and through ProVoice's own entry point, mbe_processImbe7100x4400Framef(), the
 * same parameter bits packed in ProVoice frames (pack_imbe7100()). With
 * MBE_FRAME_MARK_TONE in prev_mp_enhanced (the last frame played), in prev_mp
 * or in cur_mp before the last frame, its samples and all three states after it
 * are exactly those of the stream without the mark. (Before the rule was
 * applied after a tone model, the mark in prev_mp_enhanced changed 159 of the
 * 160 samples and the pulse phase.) */
static void
test_p25_forged_tone_in_state(void) {
    char voiced[88], unvoiced[88] = {0};
    char voiced_fr[7][24], unvoiced_fr[7][24];
    const float tone = 4.0f;
    set_imbe_voicing(voiced, 86, "0011110000");
    set_imbe_b0(unvoiced, 25);
    pack_imbe7100(voiced, voiced_fr);
    pack_imbe7100(unvoiced, unvoiced_fr);
    /* 0: P25 through the data path; 1: ProVoice through the data path with its
     * flag; 2: ProVoice through its own frame entry point, which sets it. */
    for (int entry = 0; entry < 3; ++entry) {
        for (int slot = 0; slot < 3; ++slot) {
            mbe_parms s[2][3];
            float out[2][160];
            mbe_initMbeParms(&s[0][0], &s[0][1], &s[0][2]);
            for (int f = 0; f < 4; ++f) {
                if (f == 3) {
                    memcpy(s[1], s[0], sizeof(s[1]));
                    memcpy(&s[1][slot].PHIl[0], &tone, sizeof(tone));
                }
                for (int k = 0; k < ((f == 3) ? 2 : 1); ++k) {
                    mbe_process_result result;
                    char imbe_d[88];
                    int rc;
                    mbe_initProcessResult(&result);
                    if (entry == 2) {
                        rc = mbe_processImbe7100x4400Framef(out[k], &result,
                                                            (const char (*)[24])((f == 3) ? unvoiced_fr : voiced_fr),
                                                            imbe_d, &s[k][0], &s[k][1], &s[k][2]);
                        assert((result.flags & MBE_PROCESS_FLAG_PROVOICE) != 0u);
                    } else {
                        result.flags |= (entry == 1) ? MBE_PROCESS_FLAG_PROVOICE : 0u;
                        rc = mbe_processImbe4400Dataf(out[k], &result, (f == 3) ? unvoiced : voiced, &s[k][0], &s[k][1],
                                                      &s[k][2]);
                    }
                    assert(rc >= 0);
                    (void)rc;
                }
            }
            assert(same_bits(out[0], out[1], 160));
            for (int i = 0; i < 3; ++i) {
                assert(same_stream_state(&s[0][i], &s[1][i]));
            }
        }
    }
}

/* The paths test_row17_mark_in_codec_states() plays. */
enum row17_path { ROW17_DSTAR, ROW17_P25, ROW17_RATE33, ROW17_PROVOICE, ROW17_PATHS };

/* A frame's bits for each path: parameter bits for the data paths, ProVoice
 * frames (7 x 24, pack_imbe7100()) for ProVoice's own entry point. */
struct row17_frame {
    char bits[88];
    char fr[7][24];
};

/* One frame of a path, returning its result's flags: D-STAR, P25 and rate 33
 * through their data paths, ProVoice through mbe_processImbe7100x4400Framef(). */
static unsigned
row17_path_frame(enum row17_path path, const struct row17_frame* frame, float out[160], mbe_parms s[3]) {
    mbe_process_result result;
    char imbe_d[88];
    int rc;
    mbe_initProcessResult(&result);
    switch (path) {
        case ROW17_DSTAR: rc = mbe_processAmbe2400Dataf(out, &result, frame->bits, &s[0], &s[1], &s[2]); break;
        case ROW17_P25: rc = mbe_processImbe4400Dataf(out, &result, frame->bits, &s[0], &s[1], &s[2]); break;
        case ROW17_RATE33: rc = mbe_processAmbe2450Dataf(out, &result, frame->bits, &s[0], &s[1], &s[2]); break;
        default:
            rc = mbe_processImbe7100x4400Framef(out, &result, (const char (*)[24])frame->fr, imbe_d, &s[0], &s[1],
                                                &s[2]);
            assert((result.flags & MBE_PROCESS_FLAG_PROVOICE) != 0u);
            break;
    }
    assert(rc >= 0);
    (void)rc;
    assert(all_finite(out, 160));
    return result.flags;
}

/* The last two frames of the stream test_row17_mark_in_codec_states() plays:
 * the state `s` with `mark` in slot `slot` (cur_mp, prev_mp,
 * prev_mp_enhanced), then, with `lost`, a lost frame, and a voiced frame. The
 * samples go to `out` and the three states the path returns to `t`. */
static void
row17_stream_tail(enum row17_path path, const mbe_parms s[3], int slot, enum mbe_frame_mark mark, int lost,
                  const struct row17_frame* lost_frame, const struct row17_frame* after, float out[320],
                  mbe_parms between[3], mbe_parms t[3]) {
    memcpy(t, s, 3 * sizeof(mbe_parms));
    mbe_set_frame_mark(&t[slot], mark);
    if (lost) {
        const unsigned flags = row17_path_frame(path, lost_frame, out, t);
        assert((flags & MBE_PROCESS_FLAG_REPEAT) != 0u && (flags & MBE_PROCESS_FLAG_MUTE) == 0u);
        (void)flags;
    } else {
        memset(out, 0, 160 * sizeof(float));
    }
    memcpy(between, t, 3 * sizeof(mbe_parms));
    row17_path_frame(path, after, out + 160, t);
}

/* Two returned states alike in everything but the mark word (PHIl[0]). */
static int
same_state_but_mark(const mbe_parms* a, const mbe_parms* b) {
    mbe_parms x = *a, y = *b;
    mbe_set_frame_mark(&x, MBE_FRAME_MARK_NONE);
    mbe_set_frame_mark(&y, MBE_FRAME_MARK_NONE);
    return same_stream_state(&x, &y);
}

/* Frame bits for a path: D-STAR b0 and b1, P25 (and ProVoice) b0 and a voicing
 * word ("" for none), rate 33 b0 and V/UV row. */
static void
row17_frame_for(enum row17_path path, int b0, int voicing, const char* word, struct row17_frame* frame) {
    memset(frame, 0, sizeof(*frame));
    if (path == ROW17_DSTAR) {
        set_ambe2400(frame->bits, b0, voicing);
    } else if (path == ROW17_RATE33) {
        set_ambe2450(frame->bits, b0, voicing);
    } else {
        if (word[0] != '\0') {
            set_imbe_voicing(frame->bits, b0, word);
        } else {
            set_imbe_b0(frame->bits, b0);
        }
        if (path == ROW17_PROVOICE) {
            pack_imbe7100(frame->bits, frame->fr);
        }
    }
}

/*
 * MBE_FRAME_MARK_MID_COLUMN_0 in a caller's state, in each of the three slots
 * (cur_mp, prev_mp, prev_mp_enhanced), through D-STAR's, P25's and rate 33's
 * data paths and ProVoice's own entry point (mbe_processImbe7100x4400Framef(),
 * frames packed for it). The stream: three voiced frames, a frame with no
 * voiced band through whose first half the train ran on (rate 33: V/UV row
 * 17, which its decoder marks MBE_FRAME_MARK_MID_COLUMN_0; D-STAR b0 34 and
 * P25 and ProVoice b0 25, which lack a pitch and are marked
 * MBE_FRAME_MARK_PITCH_ABSENT, the voicing before held to mid-frame in the
 * columns), then a voiced frame of another pitch, directly or after a lost
 * frame (D-STAR b0 125, P25 and ProVoice b0 255, rate 33 b0 121), which
 * repeats the last frame played or the history. Against the same state with
 * no mark in that slot (MBE_FRAME_MARK_NONE in its place):
 * - the samples: on the last frame played (prev_mp_enhanced) of the rate-33
 *   path, before the voiced frame, the train that frame stored is loaded under
 *   the mark, as under MBE_FRAME_MARK_PITCH_ABSENT, whose stream it plays bit
 *   for bit; everywhere else (every slot of D-STAR's, P25's and ProVoice's
 *   paths, rate 33's other slots, and rate 33's last frame played before a
 *   lost frame, whose repeat drops the mark and starts afresh) the samples are
 *   bit for bit those with no mark there;
 * - the three states returned: after the voiced frame, wherever the samples
 *   are alike, all three exactly those of the stream with no mark (and in the
 *   rate-33 case where the mark acts, those of the stream under
 *   PITCH_ABSENT); after a lost frame alike in everything but the mark word,
 *   which is state like any other: the forged 6 in prev_mp travels on where
 *   the path copies it (D-STAR's, P25's and ProVoice's repeat loads the
 *   history, so all three returned states hold 6 where the stream with no mark
 *   holds 0; rate 33's repeat keeps its history in prev_mp, which still holds
 *   6), and every other mark word is the stream's with no mark.
 * The mark in the slot replaces the mark the state held: on D-STAR's and P25's
 * last frame played it replaces MBE_FRAME_MARK_PITCH_ABSENT, whose train the
 * next frame would otherwise take on, so 6 there plays as no mark does, not as
 * the stream with its own mark (shown below).
 */
static void
test_row17_mark_in_codec_states(void) {
    static const enum mbe_frame_mark own_marks[ROW17_PATHS] = {MBE_FRAME_MARK_PITCH_ABSENT, MBE_FRAME_MARK_PITCH_ABSENT,
                                                               MBE_FRAME_MARK_MID_COLUMN_0,
                                                               MBE_FRAME_MARK_PITCH_ABSENT};

    for (int p = 0; p < ROW17_PATHS; ++p) {
        const enum row17_path path = (enum row17_path)p;
        struct row17_frame voiced, held, after, lost;
        row17_frame_for(path, (path == ROW17_RATE33) ? 40 : ((path == ROW17_DSTAR) ? 45 : 86),
                        (path == ROW17_DSTAR) ? 0xC : 0, "1111111100", &voiced);
        row17_frame_for(path, (path == ROW17_RATE33) ? 60 : ((path == ROW17_DSTAR) ? 34 : 25),
                        (path == ROW17_RATE33) ? 17 : 0, "", &held);
        row17_frame_for(path, (path == ROW17_RATE33) ? 50 : ((path == ROW17_DSTAR) ? 50 : 92),
                        (path == ROW17_DSTAR) ? 0xF : 0, "1111111100", &after);
        row17_frame_for(path, (path == ROW17_RATE33) ? 121 : ((path == ROW17_DSTAR) ? 125 : 255),
                        (path == ROW17_DSTAR) ? 0xF : 0, "", &lost);
        mbe_parms s[3];
        float out[160];
        mbe_initMbeParms(&s[0], &s[1], &s[2]);
        for (int f = 0; f < 3; ++f) {
            row17_path_frame(path, &voiced, out, s);
        }
        row17_path_frame(path, &held, out, s);
        /* The last frame played: no voiced band, its decoder's mark, a train
         * that ran on through its first half. */
        assert(!reference_voiced(&s[2]) && mbe_frame_mark(&s[2]) == own_marks[p]);
        assert(s[2].noiseOverlap[0] >= 1.0f && pulse_phase(&s[2]) > 0.0);
        for (int lost_first = 0; lost_first < 2; ++lost_first) {
            for (int slot = 0; slot < 3; ++slot) {
                float plain[320], marked[320];
                mbe_parms plain_s[3], marked_s[3], plain_b[3], marked_b[3];
                row17_stream_tail(path, s, slot, MBE_FRAME_MARK_NONE, lost_first, &lost, &after, plain, plain_b,
                                  plain_s);
                row17_stream_tail(path, s, slot, MBE_FRAME_MARK_MID_COLUMN_0, lost_first, &lost, &after, marked,
                                  marked_b, marked_s);
                const int acts = path == ROW17_RATE33 && slot == 2 && !lost_first;
                assert(same_bits(marked, plain, 320) == !acts);
                if (acts) {
                    /* As the marks already honoured there load it: the samples
                     * and the states those of the stream under PITCH_ABSENT. */
                    float own[320];
                    mbe_parms own_b[3], own_s[3];
                    row17_stream_tail(path, s, slot, MBE_FRAME_MARK_PITCH_ABSENT, 0, &lost, &after, own, own_b, own_s);
                    assert(same_bits(marked, own, 320));
                    for (int i = 0; i < 3; ++i) {
                        assert(same_stream_state(&marked_s[i], &own_s[i]));
                    }
                    continue;
                }
                for (int i = 0; i < 3; ++i) {
                    /* After the voiced frame: every state as with no mark. */
                    assert(same_stream_state(&marked_s[i], &plain_s[i]));
                    if (!lost_first) {
                        continue;
                    }
                    /* After the lost frame: alike but for the mark word, which
                     * the forged prev_mp carries on where the path copies it:
                     * D-STAR's, P25's and ProVoice's repeat loads the history
                     * into all three states; rate 33's keeps its history in
                     * prev_mp and repeats the last frame played. */
                    assert(same_state_but_mark(&marked_b[i], &plain_b[i]));
                    const int travels = slot == 1 && (path != ROW17_RATE33 || i == 1);
                    assert(mbe_frame_mark(&marked_b[i])
                           == (travels ? MBE_FRAME_MARK_MID_COLUMN_0 : mbe_frame_mark(&plain_b[i])));
                    assert(!travels || mbe_frame_mark(&plain_b[i]) == MBE_FRAME_MARK_NONE);
                }
            }
        }
        if (path != ROW17_RATE33) {
            /* The train the last frame played stored is heard when it runs on:
             * under the decoder's own PITCH_ABSENT the samples differ from
             * those with no mark (and so with 6) in its place. */
            float own[320], plain[320];
            mbe_parms own_b[3], own_s[3], plain_b[3], plain_s[3];
            row17_stream_tail(path, s, 2, MBE_FRAME_MARK_PITCH_ABSENT, 0, &lost, &after, own, own_b, own_s);
            row17_stream_tail(path, s, 2, MBE_FRAME_MARK_NONE, 0, &lost, &after, plain, plain_b, plain_s);
            assert(!same_bits(own, plain, 320));
        }
    }
}

/* A D-STAR harmonic exactly on a column edge belongs to the upper column, as
 * the decoder's whole part takes it: at 250 Hz (w0 2 pi / 32, so f0 is 1/32
 * and harmonic l lies at l / 2 columns exactly), 14 harmonics, 1 to 3 voiced,
 * after a frame with no voiced band, harmonic 2 at 500 Hz lies in column 1,
 * which the onset does not cover (its band, 0 to 1 kHz, is covered in column 0
 * only: column 1 is the top of the voiced region), and is unvoiced at
 * mid-frame; harmonic 1 sounds, and harmonic 3 (750 Hz, column 1) does not. */
static void
test_dstar_exact_column_edge(void) {
    mbe_parms prev, cur;
    struct mbe_subframe_model mid;
    set_model(&prev, 250.0, 14, 0, 16.0f);
    set_model(&cur, 250.0, 14, 3, 16.0f);
    prev.w0 = (2.0f * (float)M_PI) / 32.0f;
    cur.w0 = prev.w0;
    mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    const float f0 = cur.w0 / (2.0f * (float)M_PI);
    const float position = 2.0f * 16.0f * f0;
    assert(same_bits(&position, &(const float){1.0f}, 1));
    assert(mid.voiced[1] == 1 && mid.voiced[2] == 0 && mid.voiced[3] == 0);
}

/* The rate-33 rendering keeps the hold through a frame lacking a pitch at the
 * harmonic at or next above 505 Hz, found in float as the library always found
 * it, a harmonic past L taken at L: a frame of one harmonic at 100 Hz, voiced,
 * holds (its 6th harmonic, past L, read at harmonic 1); at the pitch whose bits
 * are 0x3e0761b7 (14 harmonics) the harmonic is the 3rd, where a double
 * quotient would read the 4th; at FLT_MIN it is L. The reference agrees with
 * the library in each, but at 0x3e0761b7 in builds that evaluate float wider
 * than float (x87), where the library's quotient (an expression unchanged since
 * before the columns) is not rounded to float before ceilf() and can read the
 * 4th harmonic: those two cases are not compared there. */
static void
test_rate33_hold_past_l(void) {
    static const struct {
        uint32_t w0_bits; /* 0: 100 Hz */
        int L, voiced_harmonic, held;
        int on_a_whole; /* the quotient lies on a whole number */
    } cases[] = {
        {0u, 1, 1, 1, 0},            /* the 6th harmonic, past L, read at harmonic 1 */
        {0x3e0761b7u, 14, 3, 1, 1},  /* the float quotient exactly 3: harmonic 3 */
        {0x3e0761b7u, 14, 4, 0, 1},  /* harmonic 4, which a double quotient would read, is not it */
        {0x00800000u, 14, 14, 1, 0}, /* FLT_MIN: a position of about 3.4e37, read at L */
        {0x00800000u, 14, 13, 0, 0},
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        mbe_parms prev, cur;
        struct mbe_subframe_model mid;
        unsigned char want[57] = {0};
        set_model(&prev, 100.0, cases[i].L, 0, 16.0f);
        prev.Vl[cases[i].voiced_harmonic] = 1;
        if (cases[i].w0_bits != 0u) {
            memcpy(&prev.w0, &cases[i].w0_bits, sizeof(prev.w0));
        }
        set_model(&cur, 100.0, 10, 0, 16.0f);
        mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PITCH_ABSENT);
        assert(reference_holds_505(&prev) == cases[i].held);
        mbe_subframe_mid_model(&prev, &cur, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, &mid);
        /* Where float is evaluated in excess precision (x87), the library's
         * quotient is not rounded to float before ceilf(): on a whole number
         * there it can read the next harmonic, as this hold always has. */
        if (cases[i].on_a_whole && HOLD_QUOTIENT_WIDE) {
            continue;
        }
        assert(mid.any_voiced == cases[i].held);
        reference_mid_voicing(&prev, &cur, (double)mid.w0, mid.L, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, want);
        for (int l = 1; l <= mid.L; ++l) {
            assert(mid.voiced[l] == want[l]);
        }
    }
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
 * rendering; P25's (MBE_FRAME_MARK_ONSET_COLUMNS) none since the mid-frame
 * columns change, which read every voiced model in P25's rendering. In
 * D-STAR's (and every other) rendering the marks 5, 6 and 7 play as none, and
 * a late onset acts only after a frame with no voiced band, as before these
 * marks. (P25's rendering replaces the mark of a voiced model:
 * test_p25_rendering_marks_the_model_it_renders().) A voiced frame (127.5 Hz, 23 of 28 harmonics voiced,
 * which covers the 500 Hz column 5) after an unvoiced and after a voiced one.
 * On the previous model, MBE_FRAME_MARK_MID_COLUMN_0 loads the train that
 * model stored in the rate-33 rendering, as MBE_FRAME_MARK_PITCH_ABSENT and
 * MBE_FRAME_MARK_PREVIOUS_HELD do in every rendering, and in D-STAR's and
 * P25's renders the frame after as no mark does: its samples, phases and
 * train (previous_marks_in_their_own_rendering()).
 */
static void previous_marks_in_their_own_rendering(void);

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
    previous_marks_in_their_own_rendering();
}

/* The frame after a model with no voiced band through whose first half the
 * train ran on (a 140 Hz frame lacking a pitch after a 150 Hz one voiced to
 * harmonic 20, the voicing before held to mid-frame), with that model's mark
 * replaced: a voiced 127.5 Hz frame after it, its samples, phases and train.
 * Under MBE_FRAME_MARK_MID_COLUMN_0 it renders as under
 * MBE_FRAME_MARK_PITCH_ABSENT and MBE_FRAME_MARK_PREVIOUS_HELD in the rate-33
 * rendering and as under no mark in D-STAR's and P25's; the two differ in
 * each (the train that model stored is heard). */
static void
previous_marks_in_their_own_rendering(void) {
    static const enum mbe_voiced_synthesis renderings[] = {
        MBE_VOICED_SYNTHESIS_DVSI, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_VOICED_SYNTHESIS_DVSI_IMBE};

    static const enum mbe_frame_mark marks[] = {MBE_FRAME_MARK_NONE, MBE_FRAME_MARK_MID_COLUMN_0,
                                                MBE_FRAME_MARK_PITCH_ABSENT, MBE_FRAME_MARK_PREVIOUS_HELD};

    for (size_t r = 0; r < sizeof(renderings) / sizeof(renderings[0]); ++r) {
        mbe_parms last, model, onset;
        float out[160], rendered[4][160];
        mbe_parms states[4];
        fresh_stream(&last);
        set_model(&model, 150.0, 24, 20, 16.0f);
        next_frame_as(out, &last, &model, 0, renderings[r]);
        set_model(&model, 140.0, 26, 0, 8.0f);
        mbe_set_frame_mark(&model, MBE_FRAME_MARK_PITCH_ABSENT);
        next_frame_as(out, &last, &model, 0, renderings[r]);
        assert(!reference_voiced(&last) && last.noiseOverlap[0] >= 1.0f && pulse_phase(&last) > 0.0);
        set_model(&onset, 127.5, 28, 23, 16.0f);
        shape_model(&onset, 5);
        for (size_t m = 0; m < sizeof(marks) / sizeof(marks[0]); ++m) {
            states[m] = last;
            mbe_set_frame_mark(&states[m], marks[m]);
            next_frame_as(rendered[m], &states[m], &onset, 0, renderings[r]);
            assert(all_finite(rendered[m], 160));
        }
        const int rate33 = renderings[r] == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
        assert(same_rendering(rendered[1], &states[1], rendered[rate33 ? 2 : 0], &states[rate33 ? 2 : 0]));
        assert(same_rendering(rendered[2], &states[2], rendered[3], &states[3]));
        assert(!same_rendering(rendered[0], &states[0], rendered[2], &states[2]));
    }
}

/*
 * P25's rendering marks the model it renders: a model with a voiced band that
 * is not a tone model carries MBE_FRAME_MARK_ONSET_COLUMNS whatever mark it
 * came with, so a caller's mark on it is overwritten and has no effect, the
 * late onset (80) and the held previous model (3) included, which before these
 * marks, as at the branch point, delayed the voiced bands and held the model
 * before: every one renders as the onset mark (after an unvoiced frame and
 * after a voiced one, at 151 Hz: off the pitches whose harmonics fall on a
 * column edge), the columns reading every voiced model. D-STAR's
 * rendering reads the same frame in its own columns: alike at this onset,
 * otherwise after the voiced frame (until the mid-frame columns change the two
 * rendered alike after a voiced frame, each taking this frame's voicing). A tone model keeps
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
        set_model(&model, 151.0, 24, voiced_before ? 20 : 0, 16.0f);
        next_frame(out, &last, &model);
        set_model(&model, 127.5, 28, 23, 16.0f);
        shape_model(&model, 3);
        mbe_parms onset_model = with_model(&last, &model);
        mbe_set_frame_mark(&onset_model, MBE_FRAME_MARK_ONSET_COLUMNS);
        render_marked(onset, &onset_state, &last, &model, MBE_FRAME_MARK_ONSET_COLUMNS, p25);
        assert(mbe_frame_mark(&onset_state) == MBE_FRAME_MARK_ONSET_COLUMNS);
        /* D-STAR's rendering reads the voicing in its own columns: at this
         * onset both cover columns 0 to 4 (to 2505 Hz in P25's, 2500 Hz in
         * D-STAR's) and render alike; after the voiced frame they differ (the
         * mid-frame model's harmonic 22, 3053 Hz, turns unvoiced: D-STAR's
         * column 6 holds it, the column above changing too, P25's does not). */
        render_marked(other, &other_state, &last, &model, MBE_FRAME_MARK_NONE, MBE_VOICED_SYNTHESIS_DVSI);
        assert(same_rendering(other, &other_state, onset, &onset_state) == !voiced_before);
        if (voiced_before) {
            struct mbe_subframe_model dstar_mid, p25_mid;
            mbe_subframe_mid_model(&last, &onset_model, 0, MBE_VOICED_SYNTHESIS_DVSI, &dstar_mid);
            mbe_subframe_mid_model(&last, &onset_model, 0, p25, &p25_mid);
            assert(dstar_mid.voiced[22] == 1 && p25_mid.voiced[22] == 0);
        }
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

/* ---- the pulsed state through the codecs' data paths ---- */

/* A purpose-specific generator for the tests' random fields (xorshift32). */
static uint32_t
pulsed_rng(uint32_t* x) {
    *x ^= *x << 13;
    *x ^= *x >> 17;
    *x ^= *x << 5;
    return *x;
}

/* P25 parameter bits of pitch code b0 with no voiced band: log2 amplitudes
 * random about 5 (within 1.5), quantized against the stream's prediction
 * history `history` by the IMBE encoder's own quantizer. */
static void
p25_unvoiced_bits(char imbe_d[88], int b0, uint32_t* rng, const mbe_parms* history) {
    const unsigned char bands[12] = {0};
    float logm[57] = {0.0f};
    for (int l = 1; l < 57; ++l) {
        logm[l] = 3.5f + (3.0f * (float)(pulsed_rng(rng) % 1024u) / 1023.0f);
    }
    mbe_imbe4400_quantize_amplitudes(b0, bands, logm, history, 0, imbe_d);
}

/* D-STAR parameter bits of pitch code b0 with no voiced band (b1 0), the gain
 * index b2 (bits 6..9, 42 and 43) random in 16..47, the other fields random. */
static void
dstar_unvoiced_bits(char ambe_d[49], int b0, uint32_t* rng) {
    static const int gain_bits[6] = {6, 7, 8, 9, 42, 43};
    for (int i = 0; i < 49; ++i) {
        ambe_d[i] = (char)(pulsed_rng(rng) & 1u);
    }
    const int b2 = 16 + (int)(pulsed_rng(rng) % 32u);
    for (int i = 0; i < 6; ++i) {
        ambe_d[gain_bits[i]] = (char)((b2 >> (5 - i)) & 1);
    }
    for (int i = 0; i < 6; ++i) {
        ambe_d[i] = (char)((b0 >> (6 - i)) & 1);
    }
    ambe_d[48] = (char)(b0 & 1);
    for (int i = 38; i < 42; ++i) {
        ambe_d[i] = 0;
    }
}

/* One frame through a codec's data path: the last frame as synthesized
 * (prev_mp_enhanced before the call), the prediction history before it
 * (prev_mp), this frame as synthesized (prev_mp_enhanced after it), the noise
 * state it drew from, and what it played. */
struct path_step {
    mbe_parms before, history, after;
    float seed;
    float out[160];
};

enum path_codec { PATH_DSTAR, PATH_P25, PATH_PROVOICE };

/* Run one frame of parameter bits (49 for D-STAR, 88 in the 7200 layout for
 * P25 and ProVoice, whose frames are packed for its own entry point) through
 * the codec's data path on the stream `s` (cur_mp, prev_mp, prev_mp_enhanced). */
static void
path_frame(enum path_codec codec, const char* bits, mbe_parms s[3], struct path_step* step) {
    int status;
    step->before = s[2];
    step->history = s[1];
    step->seed = s[0].noiseSeed;
    if (codec == PATH_DSTAR) {
        status = mbe_processAmbe2400Dataf(step->out, NULL, bits, &s[0], &s[1], &s[2]);
    } else if (codec == PATH_P25) {
        status = mbe_processImbe4400Dataf(step->out, NULL, bits, &s[0], &s[1], &s[2]);
    } else {
        char fr[7][24], imbe_d[88];
        mbe_process_result result;
        pack_imbe7100(bits, fr);
        status =
            mbe_processImbe7100x4400Framef(step->out, &result, (const char (*)[24])fr, imbe_d, &s[0], &s[1], &s[2]);
    }
    assert(status >= 0);
    (void)status;
    step->after = s[2];
}

/* The rendering a codec's data path synthesizes with. */
static enum mbe_voiced_synthesis
path_rendering(enum path_codec codec) {
    if (codec == PATH_DSTAR) {
        return MBE_VOICED_SYNTHESIS_DVSI;
    }
    return (codec == PATH_P25) ? MBE_VOICED_SYNTHESIS_DVSI_IMBE : MBE_VOICED_SYNTHESIS_DVSI_PROVOICE;
}

/*
 * What a frame a data path rendered with no voiced band in it or in the frame
 * before plays by the rule, on the output scale (want[]), from `last` (the
 * frame before as the path renders from it) to `cur` (this frame as
 * synthesized) on the noise state `seed`: the pulses this test renders
 * (reference_pulsed_frame(): which bands are pulsed, which model emits; the
 * mid-frame pitch of clause 5 from reference_mid_pitch(); the response, the
 * window and the low-band roll-off from this test's own arithmetic, at the
 * codec's voiced gain) plus the noise of the bands it leaves unpulsed. The
 * library's own pieces: the three models' pitch, harmonic count, magnitudes and
 * regenerated phases (mbe_subframe_own_model(), mbe_subframe_mid_model(); the
 * mid-frame pitch and harmonic count checked against this test's), and the
 * unvoiced synthesis of the unpulsed bands (mbe_synthesizeUnvoicedSubframes(),
 * DVSI's generator). Checks that `cur` stored what the frame leaves the next
 * (stored_as()) and returns the energy of the pulses.
 */
static double
reference_path_frame(const mbe_parms* last, const mbe_parms* cur, float seed, enum mbe_voiced_synthesis rendering,
                     struct reference_stream* stream, mbe_fft_plan* plan, double want[160]) {
    struct mbe_subframe_model models[3];
    int pulsed[3][57], emits[2], codes[2];
    double phase[3][57], pulses[160] = {0.0}, energy = 0.0;

    assert(!reference_voiced(last) && !reference_voiced(cur));
    mbe_subframe_own_model(last, &models[0]);
    mbe_subframe_mid_model(last, cur, 0, rendering, &models[1]);
    mbe_subframe_own_model(cur, &models[2]);
    reference_read_codes(last, cur, rendering, stream, codes);
    const float mid_pitch = reference_mid_pitch(last, cur, rendering, codes);
    assert(same_bits(&models[1].w0, &mid_pitch, 1) || (MID_PITCH_INEXACT && within_an_ulp(models[1].w0, mid_pitch)));
    const int shorter = (last->L < cur->L) ? last->L : cur->L, longer = (last->L < cur->L) ? cur->L : last->L;
    assert(models[1].L == reference_mid_count(mid_pitch, shorter, longer)
           || (MID_PITCH_INEXACT && count_margin(mid_pitch) < 1.0e-5));
    reference_pulsed_frame(cur, models, rendering, codes, pulsed, emits, stream);
    for (int i = 0; i < 3; ++i) {
        for (int l = 0; l < 57; ++l) {
            phase[i][l] = (double)models[i].phase[l];
        }
    }
    /* The low band: the rendering's filter and, in D-STAR's and P25's, the cut
     * at the corner of the tracker the state before stored (its rule is
     * test_low_cut_tracker()'s; a frame with no voiced band leaves it, so the
     * last frame's own model takes the same corner). */
    struct reference_band band = reference_band_of(rendering);
    if (band.dvsi) {
        const double track =
            mbe_subframe_state_is_marked(last) ? (double)last->noiseOverlap[MBE_LOW_CUT_TRACK_SLOT] : 0.0;
        band.corner = fmin((track >= 50.0 && track <= 500.0) ? track : 110.0, 200.0);
    }
    if (emits[0]) {
        reference_pulsed_render(pulses, &models[0], &band, pulsed[0], phase[0], &models[1], &band, pulsed[1], phase[1]);
    }
    if (emits[1]) {
        reference_pulsed_render(pulses + 80, &models[1], &band, pulsed[1], phase[1], &models[2], &band, pulsed[2],
                                phase[2]);
    }
    struct mbe_unvoiced_subframe sub[2];
    memset(sub, 0, sizeof(sub));
    for (int k = 0; k < 2; ++k) {
        const struct mbe_subframe_model* model = &models[k + 1];
        sub[k].w0 = model->w0;
        sub[k].L = model->L;
        sub[k].dvsi_low_band = band.dvsi;
        sub[k].cut_hz = (float)band.corner;
        for (int l = 1; l <= model->L; ++l) {
            if (model->voiced[l] == 0 && !pulsed[k + 1][l]) {
                sub[k].amplitude[l] = MBE_UNVOICED_SUBFRAME_SCALE * model->magnitude[l];
            }
        }
    }
    float noise[160] = {0.0f};
    mbe_parms drawn = *cur;
    drawn.noiseSeed = seed;
    mbe_synthesizeUnvoicedSubframes(noise, &drawn, last, plan, &sub[0], &sub[1]);
    for (int n = 0; n < 160; ++n) {
        want[n] = (double)MBE_SPEECH_OUTPUT_GAIN * ((double)MBE_RENDER_VOICED_GAIN_DVSI * pulses[n] + (double)noise[n]);
        assert(fabs(want[n]) < 4000.0); /* far under the output limit */
        energy += pulses[n] * pulses[n];
    }
    assert(stored_as(cur, stream, rendering));
    return energy;
}

/* Whether a sample is the reference's, within `noise_bound` more where the
 * path adds mute noise. */
static void
assert_path_sample(float out, double want, double noise_bound) {
    assert(fabs((double)out - want) <= noise_bound + 1.0e-3 + (1.0e-4 * fabs(want)));
}

/*
 * A frame a data path rendered with no voiced band in it or in the frame
 * before, against the rule (reference_path_frame()), from sample `first` on
 * (160: none). Returns the energy of the pulses.
 */
static double
check_path_frame_from(const struct path_step* step, enum mbe_voiced_synthesis rendering,
                      struct reference_stream* stream, mbe_fft_plan* plan, double noise_bound, int first) {
    double want[160];
    const double energy = reference_path_frame(&step->before, &step->after, step->seed, rendering, stream, plan, want);
    for (int n = first; n < 160; ++n) {
        assert_path_sample(step->out[n], want[n], noise_bound);
    }
    return energy;
}

static double
check_unvoiced_path_frame(const struct path_step* step, enum mbe_voiced_synthesis rendering,
                          struct reference_stream* stream, mbe_fft_plan* plan, double noise_bound) {
    return check_path_frame_from(step, rendering, stream, plan, noise_bound, 0);
}

/*
 * A D-STAR frame whose first MBE_AMBE2400_TONE_FADE (80) samples crossfade the
 * frame's speech with `other` (on the output scale): the speech (the rule's,
 * reference_path_frame(), from the frame before as synthesized to the model
 * the frame left as synthesized) fading out under a tone that fades in
 * (`speech_fades_out`: a tone's first frame, whose speech is the repeat of the
 * last frame, `other` the tone) or fading in under a tone's tail that fades
 * out (the voice frame after a tone, `other` the tail). Linearly, the speech's
 * gain 1 - n / 80 or n / 80 at sample n; after the crossfade the frame plays
 * `other` alone or the speech alone. Returns the energy of the pulses.
 */
static double
check_path_crossfade(const struct path_step* step, struct reference_stream* stream, mbe_fft_plan* plan,
                     const float other[160], int speech_fades_out) {
    double want[160];
    const double energy =
        reference_path_frame(&step->before, &step->after, step->seed, MBE_VOICED_SYNTHESIS_DVSI, stream, plan, want);
    for (int n = 0; n < 160; ++n) {
        const double fade = (n < 80) ? (double)n / 80.0 : 1.0;
        const double gain = speech_fades_out ? 1.0 - fade : fade;
        assert_path_sample(step->out[n], (gain * want[n]) + ((1.0 - gain) * (double)other[n]), 0.0);
    }
    return energy;
}

/* The pitch code of a reserved code (0..55) in a codec, or, for -1, of the
 * default code DVSI's encoders give a frame with no voiced band (P25 b0 25,
 * D-STAR b0 34). */
static int
code_b0(enum path_codec codec, int code) {
    if (codec == PATH_DSTAR) {
        return (code < 0) ? 34 : 64 + code;
    }
    return (code < 0) ? 25 : 128 + code;
}

/* Parameter bits of a frame with no voiced band at pitch code b0, in the
 * codec's layout. */
static void
unvoiced_frame_bits(enum path_codec codec, int b0, uint32_t* rng, const mbe_parms* history, char bits[88]) {
    if (codec == PATH_DSTAR) {
        dstar_unvoiced_bits(bits, b0, rng);
    } else {
        p25_unvoiced_bits(bits, b0, rng, history);
    }
}

/* One frame with no voiced band at pitch code b0 through a codec's data path:
 * its random fields drawn again (from the stream's state as it was) while the
 * frame plays loud enough to come near the output limit, which the checks
 * keep clear of. The bits it played go to `played` (NULL: nowhere). */
static void
path_unvoiced_b0(enum path_codec codec, int b0, uint32_t* rng, mbe_parms s[3], struct path_step* step,
                 char played[88]) {
    mbe_parms saved[3];
    memcpy(saved, s, sizeof(saved));
    for (int attempt = 0; attempt < 64; ++attempt) {
        char bits[88] = {0};
        float peak = 0.0f;
        memcpy(s, saved, sizeof(saved));
        unvoiced_frame_bits(codec, b0, rng, &s[1], bits);
        path_frame(codec, bits, s, step);
        for (int n = 0; n < 160; ++n) {
            peak = fmaxf(peak, fabsf(step->out[n]));
        }
        if (peak < 2000.0f) {
            if (played != NULL) {
                memcpy(played, bits, sizeof(bits));
            }
            return;
        }
    }
    assert(0 && "no frame quiet enough");
}

/* path_unvoiced_b0() at a reserved code (0..55) or the default code (-1). */
static void
path_unvoiced_frame(enum path_codec codec, int code, uint32_t* rng, mbe_parms s[3], struct path_step* step) {
    path_unvoiced_b0(codec, code_b0(codec, code), rng, s, step, NULL);
}

/* A fresh stream for a codec's data path. D-STAR's first call replaces the
 * generic initial state with the AMBE one before it synthesizes, so that the
 * state captured before it is not the one it renders from: a D-STAR stream
 * first plays one default-code frame (not checked), after which every frame
 * renders from the state captured before it. */
static void
path_stream(enum path_codec codec, uint32_t* rng, mbe_parms s[3]) {
    mbe_initMbeParms(&s[0], &s[1], &s[2]);
    if (codec == PATH_DSTAR) {
        struct path_step step;
        path_unvoiced_frame(PATH_DSTAR, -1, rng, s, &step);
    }
}

/*
 * Every reserved code of D-STAR and P25 through its codec's data path, in a
 * stream of frames with no voiced band: the default code, the reserved code
 * three times, the default code, the reserved code and the default code
 * again, each frame's amplitude fields random. Each of the seven frames is
 * checked against the rule (check_unvoiced_path_frame()): the pattern's bands
 * silent in the noise, the pulse where clause 3 puts it, emitted or not as the
 * model before decides (the reserved code's second and third frames emit at
 * sample 24 only when the frame before's mid-frame model did not), the
 * mid-frame pitch of clause 5, the state stored. ProVoice, through its own
 * entry point, plays the same bits as P25 (the P25 stream's, frame for frame)
 * at eq 46's pitch with no pulsed band. D-STAR and P25 emit pulses in every
 * code's stream.
 */
static void
test_pulsed_codes_through_data_paths(void) {
    static const int sequence[7] = {-1, 0, 0, 0, -1, 0, -1}; /* -1 the default code, 0 the reserved code */
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    for (int code = 0; code < 56; ++code) {
        for (int codec = PATH_DSTAR; codec <= PATH_P25; ++codec) {
            uint32_t rng = 0x9E3779B9u ^ (uint32_t)(code + 1);
            mbe_parms s[3];
            struct path_step step;
            struct reference_stream stream = {0, 0, 0};
            char played[7][88];
            double energy = 0.0;
            path_stream((enum path_codec)codec, &rng, s);
            for (int f = 0; f < 7; ++f) {
                path_unvoiced_b0((enum path_codec)codec, code_b0((enum path_codec)codec, (sequence[f] < 0) ? -1 : code),
                                 &rng, s, &step, played[f]);
                energy += check_unvoiced_path_frame(&step, path_rendering((enum path_codec)codec), &stream, plan, 0.0);
            }
            assert(energy > 1.0);
            if (codec != PATH_P25) {
                continue;
            }
            struct reference_stream provoice = {0, 0, 0};
            double provoice_energy = 0.0;
            mbe_initMbeParms(&s[0], &s[1], &s[2]);
            for (int f = 0; f < 7; ++f) {
                path_frame(PATH_PROVOICE, played[f], s, &step);
                provoice_energy +=
                    check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI_PROVOICE, &provoice, plan, 0.0);
            }
            assert(provoice_energy <= 0.0);
        }
    }
    mbe_fft_plan_free(plan);
}

/* Bits of a voiced P25 frame at b0 97, every band voiced. */
static void
p25_voiced_bits(char imbe_d[88], uint32_t* rng, const mbe_parms* history) {
    unsigned char bands[12];
    float logm[57] = {0.0f};
    memset(bands, 1, sizeof(bands));
    for (int l = 1; l < 57; ++l) {
        logm[l] = 3.5f + (3.0f * (float)(pulsed_rng(rng) % 1024u) / 1023.0f);
    }
    mbe_imbe4400_quantize_amplitudes(97, bands, logm, history, 0, imbe_d);
}

/* The validity counter a frame stored, by its bits: -10..3, or 99 for any other value. */
static int
stored_count_of(const mbe_parms* mp) {
    for (int count = -10; count <= 3; ++count) {
        const float value = (float)count;
        if (same_bits(&mp->noiseOverlap[TEST_COUNT_SLOT], &value, 1)) {
            return count;
        }
    }
    return 99;
}

/* One frame of the counter's walk through P25's data path at pitch code b0
 * (all bands unvoiced), checked against the rule, with the counter it must
 * leave and whether its pulses sound. */
static void
walk_frame(int b0, int count, int sounds, uint32_t* rng, mbe_parms s[3], struct reference_stream* stream,
           mbe_fft_plan* plan) {
    struct path_step step;
    path_unvoiced_b0(PATH_P25, b0, rng, s, &step, NULL);
    const double energy = check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI_IMBE, stream, plan, 0.0);
    assert(stored_count_of(&step.after) == count);
    assert(sounds ? energy > 1.0 : energy <= 0.0);
}

/* One P25 frame of pitch code b0 with no voiced band (b0 220: no usable
 * code) through the data path, every log2 amplitude 5 (level, so that the
 * adaptive smoothing an error rate brings voices no band), the error rate set
 * to `error_rate` before it (1 mutes the frame, 0 does not); returns the
 * frame's result flags. */
static uint32_t
p25_frame_at_rate(int b0, float error_rate, mbe_parms s[3]) {
    char bits[88] = {0};
    float out[160];
    mbe_process_result result;
    if (b0 > 207) {
        set_imbe_b0(bits, b0);
    } else {
        const unsigned char bands[12] = {0};
        float logm[57];
        for (int l = 0; l < 57; ++l) {
            logm[l] = 5.0f;
        }
        mbe_imbe4400_quantize_amplitudes(b0, bands, logm, &s[1], 0, bits);
    }
    mbe_initProcessResult(&result);
    s[1].errorRate = error_rate;
    const int status = mbe_processImbe4400Dataf(out, &result, bits, &s[0], &s[1], &s[2]);
    assert(status >= 0);
    (void)status;
    return result.flags;
}

static int
muted_flags(uint32_t flags) {
    return (flags & MBE_PROCESS_FLAG_MUTE) != 0u;
}

/* Whether the state says its frame's code was read. */
static int
stored_read(const mbe_parms* mp) {
    uint32_t read;
    memcpy(&read, &mp->noiseOverlap[TEST_READ_SLOT], sizeof(read));
    return read == TEST_ONE_BITS;
}

/*
 * Fifteen lost P25 frames (b0 220) after a frame of pitch code `last_b0`,
 * each checked for the model it plays and what it leaves: losses 1 to 4
 * repeat the last decoded model (repeat counts 1 to 4; the 4th is muted and
 * advances that model), loss 5 installs the repeat limit's reset model (eq 46's
 * nominal pitch of b0 134, which no code decodes to) with a repeat count of 0
 * and renders it with the mute flag clear, and the run goes on so, muting at
 * losses 9 and 14. The replayed model counts by its own code (raising a
 * counter at 3 to 3, lowering one at -10 to -10, and read only at 0 or more:
 * `count` stays), the reset model leaves the counter and reads nothing.
 */
static void
p25_lost_run(mbe_parms s[3], int count, int last_read) {
    static const int repeats[15] = {1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0};
    const float decoded = s[2].w0;
    mbe_parms prev, enhanced, reset;
    mbe_initMbeParms(&reset, &prev, &enhanced);
    for (int f = 0; f < 15; ++f) {
        const uint32_t flags = p25_frame_at_rate(220, 0.0f, s);
        assert((flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
        assert(muted_flags(flags) == (f == 3 || f == 8 || f == 13));
        assert(s[2].repeatCount == repeats[f]);
        assert(same_bits(&s[2].w0, (f < 4) ? &decoded : &reset.w0, 1));
        assert(stored_count_of(&s[2]) == count);
        assert(stored_read(&s[2]) == ((f < 4) ? last_read : 0));
    }
}

/*
 * P25's counter where the data path repeats and mutes frames. It counts the
 * model each frame renders or advances, and moves only for a model at the
 * pitch some P25 code decodes to (the document counts received frames). Lost
 * frames leave the counter where the replayed model and the reset model leave
 * it (p25_lost_run()), so the first code frame after them is read only if the
 * counter stood at 0 or more before the run: three frames at b0 176 (code 48)
 * take it to 3, fifteen losses leave it at 3, and the next b0 176 frame is
 * read; twelve frames at b0 127 take it to -10, fifteen losses leave it there,
 * and the next b0 176 frame is not read and raises it to -9. Frames muted by
 * the error rate are decoded and advanced as rendered, counting by the model
 * advanced, after adaptive smoothing: three muted frames at b0 127 (with level
 * amplitudes, which the smoothing leaves unvoiced) lower the counter from 1 to
 * -2, so the next two b0 176 frames are not read and the third is.
 */
static void
test_p25_counter_through_a_mute(void) {
    uint32_t rng = 0xA11CEu;
    mbe_parms s[3];
    path_stream(PATH_P25, &rng, s);
    for (int f = 0; f < 3; ++f) {
        assert(!muted_flags(p25_frame_at_rate(176, 0.0f, s)));
        assert(stored_count_of(&s[2]) == f + 1 && stored_read(&s[2]));
    }
    p25_lost_run(s, 3, 1);
    assert(!muted_flags(p25_frame_at_rate(176, 0.0f, s)));
    assert(stored_count_of(&s[2]) == 3 && stored_read(&s[2]));

    path_stream(PATH_P25, &rng, s);
    for (int f = 0; f < 12; ++f) {
        assert(!muted_flags(p25_frame_at_rate(127, 0.0f, s)));
        assert(stored_count_of(&s[2]) == ((f < 10) ? -(f + 1) : -10));
    }
    p25_lost_run(s, -10, 0);
    assert(!muted_flags(p25_frame_at_rate(176, 0.0f, s)));
    assert(stored_count_of(&s[2]) == -9 && !stored_read(&s[2]));

    path_stream(PATH_P25, &rng, s);
    assert(!muted_flags(p25_frame_at_rate(176, 0.0f, s)));
    assert(stored_count_of(&s[2]) == 1);
    for (int f = 0; f < 3; ++f) {
        assert(muted_flags(p25_frame_at_rate(127, 1.0f, s)));
        assert(!reference_voiced(&s[2]));
        assert(stored_count_of(&s[2]) == -f && !stored_read(&s[2]));
    }
    for (int f = 0; f < 3; ++f) {
        assert(!muted_flags(p25_frame_at_rate(176, 0.0f, s)));
        assert(stored_count_of(&s[2]) == -1 + f && stored_read(&s[2]) == (f == 2));
    }
}

/*
 * P25's validity counter through its data path (US 7,970,606, FIG. 8), every
 * frame against the rule and the counter it leaves stated here: a fresh
 * stream reads its first code (code 48: both patterns pulse, so every frame
 * whose code is read sounds a pulse); frames with no voiced band at b0 127,
 * the code this library's own P25 encoder sends them with, take the counter
 * down to -10 and keep it there; the next ten code-48 frames are not read
 * (their pulses do not sound) and bring it back to 0, the eleventh is read,
 * three more take it to its ceiling of 3; down to -10 again, ten frames at the
 * default code b0 25 bring it to 0 and the next code frame is read. A voiced
 * frame (b0 97) leaves it: a code frame after one at -1 with a voiced frame
 * between is still not read. D-STAR has no counter: after twelve frames with
 * no voiced band at b0 123 (outside its table) every code frame is read.
 */
static void
test_p25_validity_counter_walk(void) {
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    const int code48 = 128 + 48;
    uint32_t rng = 0x5EEDu;
    mbe_parms s[3];
    struct reference_stream stream = {0, 0, 0};
    int count = 1;
    assert(plan != NULL);
    path_stream(PATH_P25, &rng, s);
    walk_frame(code48, 1, 1, &rng, s, &stream, plan);
    for (int f = 0; f < 12; ++f) {
        count = (count > -10) ? count - 1 : -10;
        walk_frame(127, count, 0, &rng, s, &stream, plan);
    }
    assert(count == -10);
    for (int f = 0; f < 9; ++f) {
        walk_frame(code48, ++count, 0, &rng, s, &stream, plan);
    }
    /* At -1: a voiced frame between leaves it, and the code after it is still
     * not read (the frame after a voiced one is not the rule's alone to
     * check: its stored state is). */
    struct path_step step;
    char bits[88];
    const float zero = 0.0f;
    p25_voiced_bits(bits, &rng, &s[1]);
    path_frame(PATH_P25, bits, s, &step);
    assert(reference_voiced(&step.after) && stored_count_of(&step.after) == -1);
    path_unvoiced_b0(PATH_P25, code48, &rng, s, &step, NULL);
    assert(stored_count_of(&step.after) == 0);
    assert(same_bits(&step.after.noiseOverlap[TEST_READ_SLOT], &zero, 1));
    assert(same_bits(&step.after.noiseOverlap[TEST_EMITTED_SLOT], &zero, 1));
    stream.emitted = 0;
    stream.read = 0;
    stream.count = 0;
    walk_frame(code48, 1, 1, &rng, s, &stream, plan);
    walk_frame(code48, 2, 1, &rng, s, &stream, plan);
    walk_frame(code48, 3, 1, &rng, s, &stream, plan);
    walk_frame(code48, 3, 1, &rng, s, &stream, plan);
    count = 3;
    for (int f = 0; f < 13; ++f) {
        count = (count > -10) ? count - 1 : -10;
        walk_frame(127, count, 0, &rng, s, &stream, plan);
    }
    for (int f = 0; f < 10; ++f) {
        walk_frame(25, ++count, 0, &rng, s, &stream, plan);
    }
    assert(count == 0);
    walk_frame(code48, 1, 1, &rng, s, &stream, plan);

    /* D-STAR reads every code, whatever came before, and stores no counter. */
    struct reference_stream dstar = {0, 0, 0};
    path_stream(PATH_DSTAR, &rng, s);
    for (int f = 0; f < 12; ++f) {
        path_unvoiced_b0(PATH_DSTAR, 123, &rng, s, &step, NULL);
        (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI, &dstar, plan, 0.0);
    }
    for (int f = 0; f < 11; ++f) {
        path_unvoiced_frame(PATH_DSTAR, 48, &rng, s, &step);
        assert(check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI, &dstar, plan, 0.0) > 1.0);
        assert(stored_count_of(&step.after) == 0);
    }
    mbe_fft_plan_free(plan);
}

/* A model with no voiced band at the pitch of reserved code `code` in the
 * codec of `rendering` (reserved_pitch()), every band at `amplitude`, with
 * as many harmonics as fit under 0.9254 of Nyquist, in an initialized state
 * (blank_model()). */
static void
pulsed_model(mbe_parms* mp, enum mbe_voiced_synthesis rendering, int code, float amplitude) {
    const float w0 = reserved_pitch(rendering, code);
    const int L = (int)floor(0.9254 * M_PI / (double)w0);
    blank_model(mp, 100.0, (L > 56) ? 56 : L, 0, amplitude);
    mp->w0 = w0;
}

/* Whether any sample of a frame is not 0. */
static int
frame_sounds(const float out[160]) {
    for (int n = 0; n < 160; ++n) {
        if (fabsf(out[n]) > 0.0f) {
            return 1;
        }
    }
    return 0;
}

/* Float bit patterns a caller's state can hold in a slot of the pulsed state:
 * not-a-number (quiet, negative, signalling), the infinities, subnormals, the
 * zeros, the largest float, 0.5, 1.0 and its neighbours, every whole number a
 * counter can hold and its neighbours, and whole numbers past its range. */
static const uint32_t k_slot_values[] = {
    0x7FC00000u,
    0xFFC00000u,
    0x7F800001u,
    0x7F800000u,
    0xFF800000u,
    0x00000001u,
    0x007FFFFFu,
    0x80000001u,
    0x00000000u,
    0x80000000u,
    0x7F7FFFFFu,
    0xFF7FFFFFu,
    0x3F000000u,
    0x3F800000u,
    0x3F800001u,
    0x3F7FFFFFu,
    0xBF800000u,
    0xBF800001u,
    0x40000000u,
    0x40400000u,
    0x40400001u,
    0x403FFFFFu,
    0x40600000u,
    0x40800000u,
    0xC0000000u,
    0xC0400000u,
    0xC0800000u,
    0xC0A00000u,
    0xC0C00000u,
    0xC0E00000u,
    0xC1000000u,
    0xC1100000u,
    0xC1200000u,
    0xC1200001u,
    0xC11FFFFFu,
    0xC1280000u,
    0xC1300000u,
    0x41200000u,
    0x4B000000u,
    0xCB000000u,
    /* the neighbours of every whole number the counter holds (and so of 1.0) not listed above */
    0x3FFFFFFFu,
    0x40000001u,
    0xBF7FFFFFu,
    0xBFFFFFFFu,
    0xC0000001u,
    0xC03FFFFFu,
    0xC0400001u,
    0xC07FFFFFu,
    0xC0800001u,
    0xC09FFFFFu,
    0xC0A00001u,
    0xC0BFFFFFu,
    0xC0C00001u,
    0xC0DFFFFFu,
    0xC0E00001u,
    0xC0FFFFFFu,
    0xC1000001u,
    0xC10FFFFFu,
    0xC1100001u,
};

/* Whether k_slot_values holds `bits`. */
static int
slot_value_listed(uint32_t bits) {
    for (size_t i = 0; i < sizeof(k_slot_values) / sizeof(k_slot_values[0]); ++i) {
        if (k_slot_values[i] == bits) {
            return 1;
        }
    }
    return 0;
}

/* What a caller's value in each slot reads as: exactly 1.0 for "emitted" and
 * "read", one of the whole numbers -10..3 for the counter (else 0). */
static float
slot_reading(int slot, uint32_t bits) {
    float value;
    memcpy(&value, &bits, sizeof(value));
    if (slot != TEST_COUNT_SLOT) {
        return (bits == TEST_ONE_BITS) ? 1.0f : 0.0f;
    }
    for (int count = -10; count <= 3; ++count) {
        const float whole = (float)count;
        if (same_bits(&whole, &value, 1)) {
            return whole;
        }
    }
    return 0.0f;
}

/* A code-48 frame (both patterns pulse) after `last`, played and stored. */
static void
code48_after(const mbe_parms* last, enum mbe_voiced_synthesis rendering, float out[160], mbe_parms* stored) {
    mbe_parms frame;
    pulsed_model(&frame, rendering, 48, 8.0f);
    *stored = *last;
    next_frame_as(out, stored, &frame, 0, rendering);
}

/*
 * The pulsed state's slots are caller state, read by their bit patterns: a
 * code-48 frame after a state holding any of k_slot_values in any of the three
 * slots plays, and stores, exactly what it does after the same state holding
 * the value it reads as (slot_reading()); in D-STAR's rendering the read flag
 * and the counter are not read at all (the frame is the same whatever they
 * hold, and stores 0 in both). The readings themselves differ where the rule
 * says they do: emitted or not, read or not (P25), and a counter of -1 against
 * 0 (P25), which decides whether the code is read; and the counter each value
 * reads as is what the frame raises by one (to 3 at most). A state without
 * this rendering's subframe mark reads as a fresh stream's: not emitted, not
 * read, a counter of 0.
 */
static void
test_untrusted_pulsed_state(enum mbe_voiced_synthesis rendering) {
    static const int slots[3] = {TEST_EMITTED_SLOT, TEST_READ_SLOT, TEST_COUNT_SLOT};
    const int p25 = rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    mbe_parms first, fresh, stored, clean_stored;
    float out[160], clean[160];

    /* The list holds every value a frame stores and both its neighbours. */
    for (int count = -10; count <= 3; ++count) {
        const float value = (float)count;
        const float below = nextafterf(value, -INFINITY), above = nextafterf(value, INFINITY);
        uint32_t bits[3];
        memcpy(&bits[0], &value, sizeof(bits[0]));
        memcpy(&bits[1], &below, sizeof(bits[1]));
        memcpy(&bits[2], &above, sizeof(bits[2]));
        assert(slot_value_listed(bits[0]) && slot_value_listed(bits[1]) && slot_value_listed(bits[2]));
    }
    fresh_stream(&fresh);
    code48_after(&fresh, rendering, out, &first);
    for (int i = 0; i < 3; ++i) {
        for (size_t v = 0; v < sizeof(k_slot_values) / sizeof(k_slot_values[0]); ++v) {
            mbe_parms forged = first, canonical = first;
            /* The read flag and the counter matter to a frame whose last
             * model did not emit (the emitted slot is then the one forged). */
            forged.noiseOverlap[TEST_EMITTED_SLOT] = 0.0f;
            canonical.noiseOverlap[TEST_EMITTED_SLOT] = 0.0f;
            memcpy(&forged.noiseOverlap[slots[i]], &k_slot_values[v], sizeof(float));
            canonical.noiseOverlap[slots[i]] = (p25 || slots[i] == TEST_EMITTED_SLOT)
                                                   ? slot_reading(slots[i], k_slot_values[v])
                                                   : first.noiseOverlap[slots[i]];
            code48_after(&forged, rendering, out, &stored);
            code48_after(&canonical, rendering, clean, &clean_stored);
            assert(same_bits(out, clean, 160));
            assert(same_bits(stored.noiseOverlap, clean_stored.noiseOverlap, 96));
            if (slots[i] == TEST_COUNT_SLOT && p25) {
                const int reads_as = (int)slot_reading(TEST_COUNT_SLOT, k_slot_values[v]);
                assert(stored_count_of(&stored) == ((reads_as < 3) ? reads_as + 1 : 3));
            } else if (!p25) {
                assert(stored_count_of(&stored) == 0);
            }
        }
    }
    /* The readings differ where the rule says. */
    float a[160], b[160];
    mbe_parms x = first, y = first;
    x.noiseOverlap[TEST_EMITTED_SLOT] = 1.0f;
    y.noiseOverlap[TEST_EMITTED_SLOT] = 0.0f;
    code48_after(&x, rendering, a, &stored);
    code48_after(&y, rendering, b, &stored);
    assert(!same_bits(a, b, 160));
    if (p25) {
        x = first;
        y = first;
        x.noiseOverlap[TEST_EMITTED_SLOT] = 0.0f;
        y.noiseOverlap[TEST_EMITTED_SLOT] = 0.0f;
        y.noiseOverlap[TEST_READ_SLOT] = 0.0f;
        code48_after(&x, rendering, a, &stored);
        code48_after(&y, rendering, b, &stored);
        assert(!same_bits(a, b, 160));
        x = first;
        y = first;
        x.noiseOverlap[TEST_COUNT_SLOT] = 0.0f;
        y.noiseOverlap[TEST_COUNT_SLOT] = -1.0f;
        code48_after(&x, rendering, a, &stored);
        code48_after(&y, rendering, b, &stored);
        assert(frame_sounds(a) && !frame_sounds(b));
    }
    /* Without the subframe mark: a fresh stream's reading. */
    mbe_parms unmarked = first, zeroed = first;
    unmarked.noiseOverlap[MBE_SUBFRAME_STATE_SLOT] = 0.0f;
    for (int i = 0; i < 3; ++i) {
        unmarked.noiseOverlap[slots[i]] = 1.0f;
        zeroed.noiseOverlap[slots[i]] = 0.0f;
    }
    unmarked.noiseOverlap[TEST_COUNT_SLOT] = -1.0f;
    code48_after(&unmarked, rendering, a, &stored);
    code48_after(&zeroed, rendering, b, &clean_stored);
    assert(same_bits(a, b, 160));
    assert(same_bits(&stored.noiseOverlap[TEST_EMITTED_SLOT], &clean_stored.noiseOverlap[TEST_EMITTED_SLOT], 3));
}

/* A model of a kind the counter steps on, in P25's codec: 0 a voiced model,
 * 1 a tone model (MBE_FRAME_MARK_TONE, two voiced harmonics), 2 no voiced band
 * at the default code b0 25's pitch, 3 at a reserved code's (48), 4 at b0
 * 127's, 5 at an ulp above code 48's (no code's), 6 at the pitch of the
 * initial state and the repeat limit's reset model (eq 46's nominal pitch of
 * b0 134, no code's), 7 and 8 at the floats just below and just above the
 * default code's pitch (no code's). */
static void
count_kind_model(mbe_parms* mp, int kind) {
    if (kind == 0 || kind == 1) {
        blank_model(mp, 150.0, 24, (kind == 0) ? 20 : 2, 16.0f);
        mbe_set_frame_mark(mp, (kind == 1) ? MBE_FRAME_MARK_TONE : MBE_FRAME_MARK_NONE);
        return;
    }
    pulsed_model(mp, MBE_VOICED_SYNTHESIS_DVSI_IMBE, 48, 8.0f);
    if (kind == 2) {
        mp->w0 = p25_default_pitch();
        mp->L = 14;
    } else if (kind == 4) {
        mp->w0 = decoded_pitch(0, 127);
        mp->L = 37;
    } else if (kind == 5) {
        mp->w0 = nextafterf(mp->w0, 1.0f);
    } else if (kind == 6) {
        mbe_parms prev, enhanced, initial;
        mbe_initMbeParms(&initial, &prev, &enhanced);
        mp->w0 = initial.w0;
    } else if (kind == 7 || kind == 8) {
        mp->w0 = nextafterf(p25_default_pitch(), (kind == 7) ? 0.0f : 1.0f);
        mp->L = 14;
    }
}

/*
 * The counter's steps, frame by frame, in P25's rendering, from every value
 * it can hold (a state that stored it), rendered and advanced alike (a muted
 * frame, mbe_advance_subframes()): a voiced model and a tone model leave it,
 * a model with no voiced band at the default code's pitch or a reserved
 * code's raises it by one to 3 at most, one at another code's pitch (b0
 * 127's) lowers it by one to -10 at least, one at a pitch no code decodes to
 * (an ulp off a reserved code's or off the default code's either way, the
 * initial state's and the reset model's) leaves it; the reserved code is
 * read (stored as read, and its pulses sound) only from 0 up. The advanced
 * state stores exactly what the rendered one does. D-STAR's rendering, from
 * the same states, stores no counter and no read flag, and reads the code
 * whatever the counter says.
 */
static void
test_pulsed_count_steps(void) {
    static const int expected_step[9] = {0, 0, 1, 1, -1, 0, 0, 0, 0};
    for (int count = -10; count <= 3; ++count) {
        mbe_parms base, out_state;
        float out[160];
        fresh_stream(&base);
        code48_after(&base, MBE_VOICED_SYNTHESIS_DVSI_IMBE, out, &out_state);
        out_state.noiseOverlap[TEST_COUNT_SLOT] = (float)count;
        for (int kind = 0; kind < 9; ++kind) {
            mbe_parms model;
            count_kind_model(&model, kind);
            mbe_parms rendered = out_state;
            next_frame_as(out, &rendered, &model, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
            mbe_parms advanced = with_model(&out_state, &model);
            mbe_set_frame_mark(&advanced, mbe_frame_mark(&model));
            mbe_advance_subframes(&advanced, &out_state, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
            int want = count + expected_step[kind];
            want = (want > 3) ? 3 : ((want < -10) ? -10 : want);
            assert(stored_count_of(&rendered) == want && stored_count_of(&advanced) == want);
            assert(same_bits(&rendered.noiseOverlap[TEST_EMITTED_SLOT], &advanced.noiseOverlap[TEST_EMITTED_SLOT], 3));
            uint32_t read;
            memcpy(&read, &rendered.noiseOverlap[TEST_READ_SLOT], sizeof(read));
            assert(read == ((kind == 3 && count >= 0) ? TEST_ONE_BITS : 0u));
            if (kind == 3) {
                assert(frame_sounds(out) == (count >= 0));
            }
            mbe_parms dstar = out_state;
            next_frame_as(out, &dstar, &model, 0, MBE_VOICED_SYNTHESIS_DVSI);
            assert(stored_count_of(&dstar) == 0);
            memcpy(&read, &dstar.noiseOverlap[TEST_READ_SLOT], sizeof(read));
            assert(read == 0u);
        }
        /* D-STAR reads its own code 48 at any counter value in the state. */
        mbe_parms dstar_model, dstar_state = out_state;
        pulsed_model(&dstar_model, MBE_VOICED_SYNTHESIS_DVSI, 48, 8.0f);
        next_frame_as(out, &dstar_state, &dstar_model, 0, MBE_VOICED_SYNTHESIS_DVSI);
        assert(frame_sounds(out));
    }
    /* A fresh stream (no subframe mark) reads its first code. */
    mbe_parms fresh, model;
    float out[160];
    fresh_stream(&fresh);
    pulsed_model(&model, MBE_VOICED_SYNTHESIS_DVSI_IMBE, 48, 8.0f);
    next_frame_as(out, &fresh, &model, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    assert(frame_sounds(out) && stored_count_of(&fresh) == 1);
}

/*
 * What P25 reads of the frame before is what that frame read: a code-48 model
 * whose state says its code was not read plays no pulse of its own (the frame
 * after plays as after the same model an ulp off the code, within the change
 * of pitch), and takes no part in clause 5: between it and a code frame whose
 * code is not read either (a counter of -1), the mid-frame pitch is the
 * geometric mean, and with the read flag set it is this frame's pitch. In
 * D-STAR the flag is not read: the mid-frame pitch is the same either way.
 */
static void
test_pulsed_read_flag(void) {
    mbe_parms fresh, last, model, off_last;
    float out[160], read_out[160], off_out[160];
    fresh_stream(&fresh);
    code48_after(&fresh, MBE_VOICED_SYNTHESIS_DVSI_IMBE, out, &last);
    last.noiseOverlap[TEST_EMITTED_SLOT] = 0.0f;
    last.noiseOverlap[TEST_COUNT_SLOT] = 0.0f;
    mbe_parms unread = last;
    unread.noiseOverlap[TEST_READ_SLOT] = 0.0f;
    off_last = unread;
    off_last.w0 = nextafterf(last.w0, 1.0f);
    pulsed_model(&model, MBE_VOICED_SYNTHESIS_DVSI_IMBE, 48, 8.0f);
    mbe_parms a = unread, b = last, c = off_last;
    next_frame_as(out, &a, &model, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    next_frame_as(read_out, &b, &model, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    next_frame_as(off_out, &c, &model, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    double energy = 0.0, difference = 0.0, apart = 0.0;
    for (int n = 0; n < 160; ++n) {
        energy += (double)off_out[n] * (double)off_out[n];
        difference += ((double)out[n] - (double)off_out[n]) * ((double)out[n] - (double)off_out[n]);
        apart += ((double)out[n] - (double)read_out[n]) * ((double)out[n] - (double)read_out[n]);
    }
    assert(energy > 1.0 && difference < 1.0e-6 * energy && apart > 1.0e-2 * energy);
    /* Clause 5. The current frame at code 8's pitch, so that the two pitches differ. */
    for (int dstar = 0; dstar <= 1; ++dstar) {
        const enum mbe_voiced_synthesis rendering = dstar ? MBE_VOICED_SYNTHESIS_DVSI : MBE_VOICED_SYNTHESIS_DVSI_IMBE;
        mbe_parms prev, cur;
        struct mbe_subframe_model mid_read, mid_unread;
        fresh_stream(&fresh);
        code48_after(&fresh, rendering, out, &prev);
        prev.noiseOverlap[TEST_READ_SLOT] = 1.0f;
        prev.noiseOverlap[TEST_COUNT_SLOT] = -1.0f;
        pulsed_model(&model, rendering, 8, 8.0f);
        cur = with_model(&prev, &model);
        mbe_subframe_mid_model(&prev, &cur, 0, rendering, &mid_read);
        prev.noiseOverlap[TEST_READ_SLOT] = 0.0f;
        mbe_subframe_mid_model(&prev, &cur, 0, rendering, &mid_unread);
        const float mean = sqrtf(prev.w0 * cur.w0);
        if (dstar) {
            /* Code 8's own pattern pulses nothing: the frame before's pitch, flag or not. */
            assert(same_bits(&mid_read.w0, &prev.w0, 1) && same_bits(&mid_unread.w0, &prev.w0, 1));
        } else {
            assert(same_bits(&mid_read.w0, &cur.w0, 1));
            assert(same_bits(&mid_unread.w0, &mean, 1) || (MID_PITCH_INEXACT && within_an_ulp(mid_unread.w0, mean)));
        }
    }
}

/*
 * Only D-STAR's and P25's renderings have the pulsed state: in rate 33's,
 * ProVoice's and the TIA rendering's (MBE_VOICED_SYNTHESIS_SPEC, passed to the
 * subframe synthesizer, which the public TIA path never does) a model at
 * either codec's code-48 pitch with no voiced band plays nothing (no voiced
 * band, no transform plan, no pulse), rendered or advanced, from silence or
 * not, and stores 0 in the three slots; in D-STAR's and P25's it pulses.
 */
static void
test_pulsed_renderings(void) {
    static const enum mbe_voiced_synthesis none[3] = {MBE_VOICED_SYNTHESIS_DVSI_AMBE2450,
                                                      MBE_VOICED_SYNTHESIS_DVSI_PROVOICE, MBE_VOICED_SYNTHESIS_SPEC};

    const float zero[3] = {0.0f, 0.0f, 0.0f};
    for (int codec = 0; codec < 2; ++codec) {
        const enum mbe_voiced_synthesis own = codec ? MBE_VOICED_SYNTHESIS_DVSI_IMBE : MBE_VOICED_SYNTHESIS_DVSI;
        mbe_parms model, last;
        float out[160];
        pulsed_model(&model, own, 48, 8.0f);
        fresh_stream(&last);
        next_frame_as(out, &last, &model, 0, own);
        assert(frame_sounds(out));
        for (int r = 0; r < 3; ++r) {
            for (int from_silence = 0; from_silence <= 1; ++from_silence) {
                fresh_stream(&last);
                mbe_parms advanced = with_model(&last, &model), prev = last;
                mbe_advance_subframes(&advanced, &prev, from_silence, none[r]);
                assert(same_bits(&advanced.noiseOverlap[TEST_EMITTED_SLOT], zero, 3));
                for (int f = 0; f < 2; ++f) {
                    next_frame_as(out, &last, &model, from_silence, none[r]);
                    assert(!frame_sounds(out));
                    assert(same_bits(&last.noiseOverlap[TEST_EMITTED_SLOT], zero, 3));
                }
            }
        }
    }
}

/* A frame of pitch code b0 with no voiced band through the codec's public
 * parameter decoder (D-STAR's when `dstar`, P25's else), every band at
 * amplitude 8 (the decode's own amplitudes are those of all-zero bits);
 * returns 0 if the code does not decode to a voice frame with no voiced
 * band. */
static int
decoded_unvoiced_model(int dstar, int b0, mbe_parms* model) {
    mbe_parms prev, enhanced;
    mbe_initMbeParms(model, &prev, &enhanced);
    int status;
    if (dstar) {
        char ambe_d[49] = {0};
        for (int i = 0; i < 6; ++i) {
            ambe_d[i] = (char)((b0 >> (6 - i)) & 1);
        }
        ambe_d[48] = (char)(b0 & 1);
        status = mbe_decodeAmbe2400Parms(ambe_d, model, &prev);
    } else {
        char imbe_d[88] = {0};
        set_imbe_b0(imbe_d, b0);
        status = mbe_decodeImbe4400Parms(imbe_d, model, &prev);
    }
    if (status != 0 || reference_voiced(model) || !(model->w0 > 0.0f) || model->L < 1 || model->L > 56) {
        return 0;
    }
    for (int l = 1; l <= 56; ++l) {
        set_amplitude(model, l, (l <= model->L) ? 8.0f : 0.0f);
    }
    mbe_set_frame_mark(model, MBE_FRAME_MARK_NONE);
    return 1;
}

/*
 * Only the reserved codes carry a pattern, and every code steps P25's counter
 * as the rule says: every P25 code (0..207) and D-STAR code (0..123) the public
 * decoders give as a voice frame with no voiced band, and the test requires
 * every one of them to decode so (208 and 124, none skipped). An ordinary code
 * plays nothing in its codec's rendering (no voiced band, no transform plan: a
 * pattern would sound a pulse) in a stream's first frame or second, while each
 * reserved code sounds in its second. From a fresh stream, a P25 code of the
 * table (b0 25 and the 56) leaves the counter at 1 and then 2, every other P25
 * code at -1 and then -2; D-STAR leaves it 0.
 */
static void
test_ordinary_codes_carry_no_pattern(void) {
    for (int dstar = 0; dstar <= 1; ++dstar) {
        const enum mbe_voiced_synthesis rendering = dstar ? MBE_VOICED_SYNTHESIS_DVSI : MBE_VOICED_SYNTHESIS_DVSI_IMBE;
        const int first = dstar ? 64 : 128;
        const int codes = dstar ? 124 : 208;
        int ordinary = 0, reserved = 0;
        for (int b0 = 0; b0 < codes; ++b0) {
            mbe_parms model, last;
            float out[2][160];
            const int decoded = decoded_unvoiced_model(dstar, b0, &model);
            assert(decoded);
            (void)decoded;
            fresh_stream(&last);
            const int is_reserved = b0 >= first && b0 < first + 56;
            const int step = dstar ? 0 : ((is_reserved || b0 == 25) ? 1 : -1);
            for (int f = 0; f < 2; ++f) {
                next_frame_as(out[f], &last, &model, 0, rendering);
                assert(stored_count_of(&last) == step * (f + 1));
            }
            /* A code whose own pattern alone pulses emits in its second frame. */
            assert(is_reserved ? frame_sounds(out[1]) : (!frame_sounds(out[0]) && !frame_sounds(out[1])));
            ordinary += !is_reserved;
            reserved += is_reserved;
        }
        assert(reserved == 56 && ordinary == codes - 56);
    }
}

/*
 * What has a code is the model itself, and nothing else: a pitch one ulp from
 * a reserved code's has none (the frame plays nothing at all with no voiced
 * band and no transform plan), the reserved pitch itself pulses; a band with
 * no state (a voicing neither 0 nor 1) is not pulsed, so it plays as an
 * unvoiced band of no amplitude would (the same log2Ml, so the same phases);
 * and a voiced band takes the code away. No frame mark of either frame
 * changes which bands are pulsed, but those that put another model at
 * mid-frame (MBE_FRAME_MARK_PREVIOUS_HELD, MBE_FRAME_MARK_TONE), whose
 * mid-frame model takes no pattern: a code-40 frame (its mid-frame pattern
 * alone pulses) after another plays only its mid-frame pulse, and nothing
 * with those two marks. A caller's model at a reserved pitch is a frame with
 * that code: there is nothing to forge.
 */
static void
test_pulsed_state_is_the_models(enum mbe_voiced_synthesis rendering) {
    static const float marks[] = {0.0f, 1.0f, 2.0f, 5.0f, 6.0f, 7.0f, 80.0f, 1.5f, -1.0f, 1000.0f};
    mbe_parms last, frame;
    float out[160], reference[160];

    for (int code = 0; code < 56; ++code) {
        pulsed_model(&frame, rendering, code, 8.0f);
        fresh_stream(&last);
        next_frame_as(out, &last, &frame, 0, rendering);
        next_frame_as(out, &last, &frame, 0, rendering);
        assert(frame_sounds(out));
        for (int step = -1; step <= 1; step += 2) {
            mbe_parms off = frame;
            off.w0 = nextafterf(frame.w0, (step < 0) ? 0.0f : 1.0f);
            fresh_stream(&last);
            next_frame_as(out, &last, &off, 0, rendering);
            next_frame_as(out, &last, &off, 0, rendering);
            assert(!frame_sounds(out));
        }
        /* One band voiced (as adaptive smoothing could leave it): no code. */
        mbe_parms voiced = frame;
        voiced.Vl[3] = 1;
        voiced.Ml[3] = 0.0f;
        fresh_stream(&last);
        next_frame_as(out, &last, &voiced, 0, rendering);
        next_frame_as(out, &last, &voiced, 0, rendering);
        assert(!frame_sounds(out));
    }
    /* A band with no state is not pulsed: it plays as an unvoiced band of no
     * amplitude would, within 1e-6 of the frame's energy. Not bit for bit: the
     * float quotient of a mid-frame harmonic's frequency over a frame's pitch
     * need not be the whole number it stands for, and the log-domain
     * interpolation then gives the neighbouring band a weight of the order of
     * 1e-7, under which that band's log2 magnitude (3 for an amplitude of 8,
     * or the floor of a band of no amplitude, -1000) moves the mid-frame
     * magnitudes by about 1e-4. Code 40: its mid-frame pattern pulses every
     * channel. */
    pulsed_model(&frame, rendering, 40, 8.0f);
    mbe_parms none = frame, silent = frame;
    for (int l = 2; l <= frame.L; l += 3) {
        none.Vl[l] = 2;
        silent.Ml[l] = 0.0f;
    }
    fresh_stream(&last);
    next_frame_as(out, &last, &none, 0, rendering);
    next_frame_as(out, &last, &none, 0, rendering);
    fresh_stream(&last);
    next_frame_as(reference, &last, &silent, 0, rendering);
    next_frame_as(reference, &last, &silent, 0, rendering);
    double energy = 0.0, difference = 0.0;
    for (int n = 0; n < 160; ++n) {
        energy += (double)reference[n] * (double)reference[n];
        difference += ((double)out[n] - (double)reference[n]) * ((double)out[n] - (double)reference[n]);
    }
    assert(energy > 1.0 && difference < 1.0e-6 * energy);
    /* Marks. */
    pulsed_model(&frame, rendering, 40, 8.0f);
    fresh_stream(&last);
    next_frame_as(reference, &last, &frame, 0, rendering);
    const mbe_parms first = last;
    next_frame_as(reference, &last, &frame, 0, rendering);
    assert(frame_sounds(reference));
    for (size_t i = 0; i < sizeof(marks) / sizeof(marks[0]); ++i) {
        for (size_t j = 0; j < sizeof(marks) / sizeof(marks[0]); ++j) {
            mbe_parms marked = frame;
            last = first;
            last.PHIl[0] = marks[i];
            marked.PHIl[0] = marks[j];
            /* next_frame_as() takes the mark from the model's mbe_frame_mark(). */
            mbe_parms cur = last;
            cur.w0 = marked.w0;
            cur.L = marked.L;
            memcpy(cur.Vl, marked.Vl, sizeof(cur.Vl));
            copy_amplitudes(&cur, &marked);
            cur.PHIl[0] = marks[j];
            memset(out, 0, sizeof(out));
            mbe_synthesize_subframes(out, &cur, &last, NULL, 1.0f, 0, rendering);
            assert(same_bits(out, reference, 160));
        }
    }
    for (int k = 0; k < 2; ++k) {
        mbe_parms cur = first;
        cur.w0 = frame.w0;
        cur.L = frame.L;
        memcpy(cur.Vl, frame.Vl, sizeof(cur.Vl));
        copy_amplitudes(&cur, &frame);
        mbe_set_frame_mark(&cur, (k == 0) ? MBE_FRAME_MARK_PREVIOUS_HELD : MBE_FRAME_MARK_TONE);
        last = first;
        memset(out, 0, sizeof(out));
        mbe_synthesize_subframes(out, &cur, &last, NULL, 1.0f, 0, rendering);
        assert(!frame_sounds(out));
    }
}

/*
 * A caller's MBE_FRAME_MARK_PITCH_ABSENT on a frame with a code after a voiced
 * frame holds that frame's covered columns to mid-frame, as on any frame with
 * no voiced band (the rules of the mid-frame voicing are a code's as a frame's
 * with no voiced band), and a mid-frame model with a voiced band takes no
 * pulsed band: the frame (code 40, whose mid-frame pattern alone pulses)
 * plays as the same frame one ulp off the code does, within the change of
 * pitch, and without the mark its mid-frame pulse sounds.
 */
static void
test_pulsed_held_by_a_forged_mark(enum mbe_voiced_synthesis rendering) {
    mbe_parms start, voiced, frame;
    float out[160], off_out[160], plain_out[160];
    fresh_stream(&start);
    blank_model(&voiced, 150.0 * 1.001, 24, 24, 16.0f);
    shape_model(&voiced, 3);
    next_frame_as(out, &start, &voiced, 0, rendering);
    next_frame_as(out, &start, &voiced, 0, rendering);
    pulsed_model(&frame, rendering, 40, 8.0f);
    mbe_set_frame_mark(&frame, MBE_FRAME_MARK_PITCH_ABSENT);
    mbe_parms last = start;
    next_frame_as(out, &last, &frame, 0, rendering);
    mbe_parms off = frame;
    off.w0 = nextafterf(frame.w0, 1.0f);
    last = start;
    next_frame_as(off_out, &last, &off, 0, rendering);
    double energy = 0.0, difference = 0.0;
    for (int n = 0; n < 160; ++n) {
        energy += (double)out[n] * (double)out[n];
        difference += ((double)out[n] - (double)off_out[n]) * ((double)out[n] - (double)off_out[n]);
    }
    assert(energy > 1.0 && difference < 1.0e-6 * energy);
    mbe_set_frame_mark(&frame, MBE_FRAME_MARK_NONE);
    last = start;
    next_frame_as(plain_out, &last, &frame, 0, rendering);
    double apart = 0.0;
    for (int n = 80; n < 160; ++n) {
        apart += ((double)plain_out[n] - (double)out[n]) * ((double)plain_out[n] - (double)out[n]);
    }
    assert(apart > 1.0e-2 * energy);
}

/*
 * A muted frame (mbe_advance_subframes()) stores whether its mid-frame model
 * emitted, as rendering it does: a code-48 frame after a default-code frame
 * (no pattern) emits from its mid-frame model and so stores it, and the next
 * code-48 frame's first model then does not emit; rendered after the advanced
 * frame it plays exactly as after the rendered one (no transform plan, so no
 * unvoiced overlap is lost). The advanced state starts from a slot that says
 * otherwise. (The counter and the read flag: test_pulsed_count_steps().)
 */
static void
test_pulsed_advance_stores_the_emission(enum mbe_voiced_synthesis rendering) {
    mbe_parms rendered, advanced, frame, plain;
    float out[160], after_rendered[160], after_advanced[160];
    fresh_stream(&rendered);
    blank_model(&plain, 125.0 * 1.001, 14, 0, 8.0f);
    plain.w0 = edge_pitch(rendering, -3);
    next_frame_as(out, &rendered, &plain, 0, rendering);
    advanced = rendered;
    pulsed_model(&frame, rendering, 48, 8.0f);
    next_frame_as(out, &rendered, &frame, 0, rendering);
    mbe_parms cur = with_model(&advanced, &frame);
    cur.noiseOverlap[TEST_EMITTED_SLOT] = 0.25f;
    mbe_advance_subframes(&cur, &advanced, 0, rendering);
    advanced = cur;
    const uint32_t one_bits = TEST_ONE_BITS;
    float one;
    memcpy(&one, &one_bits, sizeof(one));
    assert(same_bits(&rendered.noiseOverlap[TEST_EMITTED_SLOT], &one, 1));
    assert(same_bits(&advanced.noiseOverlap[TEST_EMITTED_SLOT], &one, 1));
    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
    next_frame_as(after_rendered, &rendered, &frame, 0, rendering);
    next_frame_as(after_advanced, &advanced, &frame, 0, rendering);
    assert(frame_sounds(after_rendered) && same_bits(after_rendered, after_advanced, 160));
}

/*
 * Adaptive smoothing can voice a band of a frame decoded with none; the frame
 * then has no code. A P25 frame at a reserved pitch with one loud band and the
 * error count that makes the smoothing voice it, through the synthesis core:
 * the frame it renders has a voiced band and plays as the same frame one ulp
 * off the code does, within the change of pitch, and stores no emission, no
 * read code and the counter it found (0: a voiced model leaves it). The same
 * frame with an error rate over the muting threshold, which the core mutes and
 * advances, is voiced by the smoothing as well and stores the same.
 */
static void
test_pulsed_smoothing_takes_the_code(void) {
    mbe_parms prev, frame, off, first, second;
    float out[160], off_out[160];
    mbe_initMbeParms(&prev, &first, &second);
    pulsed_model(&frame, MBE_VOICED_SYNTHESIS_DVSI_IMBE, 11, 4.0f);
    frame.Ml[5] = 4000.0f;
    frame.errorCountTotal = 5;
    frame.errorRate = 0.002f;
    frame.localEnergy = 1.0f;
    frame.mutingThreshold = 0.0875f;
    off = frame;
    off.w0 = nextafterf(frame.w0, 1.0f);
    first = frame;
    second = off;
    mbe_parms prev_a = prev, prev_b = prev;
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &first, &prev_a, 1.0f, MBE_MUTE_NOISE_SPEC, MBE_VOICED_SYNTHESIS_DVSI_IMBE,
                                       MBE_AMPLITUDE_LIMIT_NONE);
    mbe_synthesizeSpeechWithPreEnhRm0f(off_out, &second, &prev_b, 1.0f, MBE_MUTE_NOISE_SPEC,
                                       MBE_VOICED_SYNTHESIS_DVSI_IMBE, MBE_AMPLITUDE_LIMIT_NONE);
    assert(reference_voiced(&first));
    double energy = 0.0, difference = 0.0;
    for (int n = 0; n < 160; ++n) {
        energy += (double)out[n] * (double)out[n];
        difference += ((double)out[n] - (double)off_out[n]) * ((double)out[n] - (double)off_out[n]);
    }
    assert(energy > 1.0 && difference < 1.0e-6 * energy);
    const float zeros[3] = {0.0f, 0.0f, 0.0f};
    assert(same_bits(&first.noiseOverlap[TEST_EMITTED_SLOT], zeros, 3));
    mbe_parms muted = frame, prev_c = prev;
    float muted_out[160];
    muted.errorRate = 0.2f;
    mbe_synthesizeSpeechWithPreEnhRm0f(muted_out, &muted, &prev_c, 1.0f, MBE_MUTE_NOISE_SPEC,
                                       MBE_VOICED_SYNTHESIS_DVSI_IMBE, MBE_AMPLITUDE_LIMIT_NONE);
    /* Muted: the specification's mute noise, far under the rendered frame. */
    assert(peak_of(muted_out) < 0.01f * peak_of(out));
    assert(reference_voiced(&muted) && mbe_subframe_state_is_marked(&muted));
    assert(same_bits(&muted.noiseOverlap[TEST_EMITTED_SLOT], zeros, 3));
}

/* D-STAR bits of a lost frame (b0 125), or of a tone frame (b0 126) of tone
 * index `index` at volume 0x40. */
static void
dstar_special_bits(char ambe_d[49], int b0, int index) {
    memset(ambe_d, 0, 49);
    for (int i = 0; i < 6; ++i) {
        ambe_d[i] = (char)((b0 >> (6 - i)) & 1);
    }
    ambe_d[48] = (char)(b0 & 1);
    if (b0 == 126) {
        mbe_ambe2400_set_tone_index(ambe_d, index);
        mbe_tone_dstar_set_volume(ambe_d, 0x40);
    }
}

/*
 * The state D-STAR's first voice frame after the repeat limit's mute renders
 * from, as its decoder builds it (ambe2400_resume_after_mute() in
 * ambe3600x2400.c, restated here): the last frame as the mute left it, given
 * the prediction history's model (pitch, harmonic count, voicing, log2Ml) at
 * 0.3 times the amplitudes the decoder builds from log2Ml (an unvoiced band's
 * with the factor 0.2046 / sqrt(w0)), then the library's spectral enhancement
 * (mbe_spectralAmpEnhanceWithRm0()), and the unvoiced overlap a frame of that
 * model leaves, from the noise block after the one the next frame draws: two
 * scratch repeats of it by the library's repeat synthesis
 * (mbe_synthesizeRepeatedSpeechf()), the second's overlap kept, as the
 * decoder forms it (ambe2400_resume_overlap()).
 */
static void
dstar_resumed_state(const struct path_step* step, mbe_parms* resumed) {
    const mbe_parms* history = &step->history;
    mbe_parms enhanced = step->before;
    const float unvoiced_scale = 0.2046f / sqrtf(history->w0);
    enhanced.L = history->L;
    enhanced.w0 = history->w0;
    memcpy(enhanced.log2Ml, history->log2Ml, sizeof(enhanced.log2Ml));
    memset(enhanced.Ml, 0, sizeof(enhanced.Ml));
    for (int l = 1; l <= history->L; ++l) {
        enhanced.Vl[l] = history->Vl[l];
        enhanced.Ml[l] = 0.3f * ((history->Vl[l] == 1) ? 1.0f : unvoiced_scale) * exp2f(history->log2Ml[l]);
    }
    (void)mbe_spectralAmpEnhanceWithRm0(&enhanced, MBE_HIGH_BAND_GAIN_AMBE2400);
    mbe_parms first = enhanced, second = enhanced, prev = enhanced;
    float scratch[160];
    first.repeatCount = 0;
    second.repeatCount = 0;
    prev.repeatCount = 0;
    mbe_synthesizeRepeatedSpeechf(scratch, &first, &prev, MBE_VOICED_SYNTHESIS_DVSI);
    second.noiseSeed = first.noiseSeed;
    memcpy(second.noiseOverlap, first.noiseOverlap, sizeof(second.noiseOverlap));
    prev = enhanced;
    prev.repeatCount = 0;
    mbe_synthesizeRepeatedSpeechf(scratch, &second, &prev, MBE_VOICED_SYNTHESIS_DVSI);
    *resumed = step->before;
    memcpy(resumed->previousUw, second.previousUw, sizeof(resumed->previousUw));
    mbe_subframe_state_adopt(resumed);
    memcpy(resumed->Ml, enhanced.Ml, sizeof(resumed->Ml));
    memcpy(resumed->log2Ml, enhanced.log2Ml, sizeof(resumed->log2Ml));
    memcpy(&resumed->Vl[1], &enhanced.Vl[1], sizeof(enhanced.Vl) - sizeof(enhanced.Vl[0]));
    resumed->L = enhanced.L;
    resumed->w0 = enhanced.w0;
}

/* The tone D-STAR index 130 at volume 0x40 plays from fresh phases (a tone's
 * first frame), or from the phases a state holds (the tail after it), by the
 * library's tone renderer. */
static void
dstar_tone_samples(float out[160], const mbe_parms* phases) {
    mbe_parms tone;
    memset(&tone, 0, sizeof(tone));
    memset(out, 0, 160 * sizeof(float));
    if (phases == NULL) {
        mbe_startToneDstar(&tone, 130);
    } else {
        tone.swn = phases->swn;
        tone.tonePhase = phases->tonePhase;
    }
    mbe_renderToneDstarf(out, &tone, 130, 0x40);
}

/*
 * The pulsed state through D-STAR's lost frames, its repeat limit's mute, the
 * voice after it and its tone crossfades, and P25's repeats, each frame
 * against the rule (reference_path_frame()) within one stream each: a lost
 * frame repeats a model with its code, the alternation going on; the repeat
 * limit's 4th frame fades the last frame out (whose own model emits at sample
 * 24 unless its frame's mid-frame model did) under DVSI's mute noise, within
 * that noise's amplitude; the next voice frame resumes from the prediction
 * history at 0.3 of its amplitudes (dstar_resumed_state(), its unvoiced
 * overlap rebuilt with the library's repeat synthesis), compared wholly, and
 * the two after it wholly.
 * Into a tone (DTMF, D-STAR index 130): the tone's first frame crossfades the
 * repeat of the last frame (the rule's, to the model it left) out under the
 * tone over its first 80 samples and plays the tone alone after them; the
 * voice frame after the tone crossfades in under the tone's tail over its
 * first 80, the one after it plays wholly (check_path_crossfade(); the tone
 * and its tail from the library's tone renderer). A P25 frame with no usable
 * pitch code (b0 220) repeats the last model, code and all.
 */
static void
test_pulsed_repeats_mutes_and_tones(void) {
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    const double mute = (double)MBE_AMBE_MUTE_NOISE_AMPLITUDE * (double)MBE_SPEECH_OUTPUT_GAIN;
    const int codes[] = {5, 11, 48, 50};
    assert(plan != NULL);
    for (size_t c = 0; c < sizeof(codes) / sizeof(codes[0]); ++c) {
        uint32_t rng = 0xC0FFEEu + (uint32_t)c;
        mbe_parms s[3], resumed;
        struct path_step step;
        struct reference_stream stream = {0, 0, 0};
        char bits[88];
        double want[160];
        /* Lost frames up to the repeat limit's mute, and the voice after it. */
        path_stream(PATH_DSTAR, &rng, s);
        path_unvoiced_frame(PATH_DSTAR, codes[c], &rng, s, &step);
        (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI, &stream, plan, 0.0);
        for (int lost = 1; lost <= 4; ++lost) {
            dstar_special_bits(bits, 125, 0);
            path_frame(PATH_DSTAR, bits, s, &step);
            (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI, &stream, plan,
                                            (lost == 4) ? mute + 1.0e-6 : 0.0);
        }
        path_unvoiced_frame(PATH_DSTAR, codes[c], &rng, s, &step);
        dstar_resumed_state(&step, &resumed);
        (void)reference_path_frame(&resumed, &step.after, step.seed, MBE_VOICED_SYNTHESIS_DVSI, &stream, plan, want);
        for (int n = 0; n < 160; ++n) {
            assert_path_sample(step.out[n], want[n], 0.0);
        }
        for (int f = 0; f < 2; ++f) {
            path_unvoiced_frame(PATH_DSTAR, codes[c], &rng, s, &step);
            (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI, &stream, plan, 0.0);
        }
        /* Into and out of a tone. */
        float tone[160];
        memset(&stream, 0, sizeof(stream));
        path_stream(PATH_DSTAR, &rng, s);
        for (int f = 0; f < 2; ++f) {
            path_unvoiced_frame(PATH_DSTAR, codes[c], &rng, s, &step);
            (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI, &stream, plan, 0.0);
        }
        dstar_special_bits(bits, 126, 130);
        path_frame(PATH_DSTAR, bits, s, &step);
        dstar_tone_samples(tone, NULL);
        (void)check_path_crossfade(&step, &stream, plan, tone, 1);
        path_unvoiced_frame(PATH_DSTAR, codes[c], &rng, s, &step);
        dstar_tone_samples(tone, &step.history);
        (void)check_path_crossfade(&step, &stream, plan, tone, 0);
        path_unvoiced_frame(PATH_DSTAR, codes[c], &rng, s, &step);
        (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI, &stream, plan, 0.0);
    }
    for (size_t c = 0; c < sizeof(codes) / sizeof(codes[0]); ++c) {
        uint32_t rng = 0xFACADEu + (uint32_t)c;
        mbe_parms s[3];
        struct path_step step;
        char bits[88];
        struct reference_stream stream = {0, 0, 0};
        path_stream(PATH_P25, &rng, s);
        path_unvoiced_frame(PATH_P25, -1, &rng, s, &step);
        (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &stream, plan, 0.0);
        path_unvoiced_frame(PATH_P25, codes[c], &rng, s, &step);
        (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &stream, plan, 0.0);
        for (int lost = 1; lost <= 2; ++lost) {
            memset(bits, 0, sizeof(bits));
            set_imbe_b0(bits, 220);
            path_frame(PATH_P25, bits, s, &step);
            (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &stream, plan, 0.0);
        }
        path_unvoiced_frame(PATH_P25, codes[c], &rng, s, &step);
        (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &stream, plan, 0.0);
    }
    mbe_fft_plan_free(plan);
}

/* Frame marks a caller's state can hold: every mark (0..7), 80 (a P25 frame
 * voiced only below 505 Hz, as P25's decoder once marked it), and values
 * that are no mark (1.5, -1, a large one, not-a-number). */
static const uint32_t k_mark_values[] = {0x00000000u, 0x3F800000u, 0x40000000u, 0x40400000u, 0x40800000u,
                                         0x40A00000u, 0x40C00000u, 0x40E00000u, 0x42A00000u, 0x3FC00000u,
                                         0xBF800000u, 0x447A0000u, 0x7FC00000u};

/*
 * ProVoice has no pulsed state, whatever its states hold: a ProVoice frame
 * after a caller's model at a P25 reserved pitch with no voiced band (code 32,
 * whose own pattern pulses every channel), with each value of k_slot_values
 * in each of the three slots of the pulsed state (the others as a P25 frame
 * of that code leaves them: emitted 0, read 1.0, a counter of 1) and with
 * each mark of k_mark_values, plays as the rule without the pulsed state says
 * (noise in every unvoiced band, no pulse, the geometric mean at mid-frame),
 * and stores 0 in the three slots.
 */
static void
test_provoice_has_no_pulsed_state(void) {
    static const int slots[3] = {TEST_EMITTED_SLOT, TEST_READ_SLOT, TEST_COUNT_SLOT};
    const size_t values = sizeof(k_slot_values) / sizeof(k_slot_values[0]);
    const size_t marks = sizeof(k_mark_values) / sizeof(k_mark_values[0]);
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    for (size_t i = 0; i < 3 * values + marks; ++i) {
        uint32_t rng = 0xBEEFu + (uint32_t)i;
        mbe_parms s[3];
        struct path_step step;
        struct reference_stream stream = {0, 0, 0};
        mbe_initMbeParms(&s[0], &s[1], &s[2]);
        path_unvoiced_frame(PATH_PROVOICE, 32, &rng, s, &step);
        s[2].w0 = reserved_pitch(MBE_VOICED_SYNTHESIS_DVSI_IMBE, 32);
        s[2].noiseOverlap[TEST_EMITTED_SLOT] = 0.0f;
        s[2].noiseOverlap[TEST_READ_SLOT] = 1.0f;
        s[2].noiseOverlap[TEST_COUNT_SLOT] = 1.0f;
        if (i < 3 * values) {
            memcpy(&s[2].noiseOverlap[slots[i / values]], &k_slot_values[i % values], sizeof(float));
        } else {
            memcpy(&s[2].PHIl[0], &k_mark_values[i - 3 * values], sizeof(float));
        }
        path_unvoiced_frame(PATH_PROVOICE, 32, &rng, s, &step);
        (void)check_unvoiced_path_frame(&step, MBE_VOICED_SYNTHESIS_DVSI_PROVOICE, &stream, plan, 0.0);
        assert(!stream.emitted);
    }
    mbe_fft_plan_free(plan);
}

/* ---- rate 33's pulsed rows through the decoder and the data path ---- */

/* Rate-33 tone frame bits (TIA-102.BABA-1 7.2 Table 10) of tone `id` at
 * amplitude field `ad`. */
static void
r33_tone_bits(char d[49], int id, int ad) {
    const int fields[][3] = {{0, 63, 6},       {6, ad >> 1, 6}, {12, id, 8}, {20, id >> 4, 4}, {24, id, 4},
                             {28, id >> 1, 7}, {35, id, 1},     {36, id, 8}, {44, ad, 1},      {45, 0, 4}};
    for (size_t f = 0; f < sizeof(fields) / sizeof(fields[0]); ++f) {
        for (int i = 0; i < fields[f][2]; ++i) {
            d[fields[f][0] + i] = (char)((fields[f][1] >> (fields[f][2] - 1 - i)) & 1);
        }
    }
}

/* A PHIl[0] (or slot) holds the word of row, b0 and held, by its bits. */
static int
holds_word(const float* stored, int row, int b0, int held) {
    uint32_t bits;
    memcpy(&bits, stored, sizeof(bits));
    return bits == test_row_word_bits(test_row_word(row, b0, held));
}

/* Whether a stored value is at least the smallest word, 1024.0f, by its bits:
 * no mark (0..80) is. */
static int
holds_a_word(const float* stored) {
    uint32_t bits;
    memcpy(&bits, stored, sizeof(bits));
    return bits >= test_row_word_bits(TEST_ROW_WORD_BASE) && bits < UINT32_C(0x80000000);
}

/* One rate-33 frame through its data path, with C0's error count when
 * `c0_errors` is 0 or more (4 makes the frame a rejected repeat). */
static unsigned
r33_frame(const char bits[49], float out[160], mbe_parms s[3], int c0_errors) {
    mbe_process_result result;
    mbe_initProcessResult(&result);
    if (c0_errors >= 0) {
        result.flags |= MBE_PROCESS_FLAG_C0_VALID;
        result.c0_errors = c0_errors;
    }
    const int rc = mbe_processAmbe2450Dataf(out, &result, bits, &s[0], &s[1], &s[2]);
    assert(rc >= 0);
    (void)rc;
    assert(all_finite(out, 160));
    return result.flags;
}

/*
 * The decoder's word. Every V/UV row 0..31 at every voice code 0..119, decoded
 * from a fresh state by the public decoder and through the data path: a row
 * of 18..31 gives exactly the word of that row and code in cur_mp (the public
 * decoder) and in all three returned states (the data path), which
 * mbe_frame_mark() reads as no mark, on the band-grid model of the code; no
 * other row gives a word. No silence frame (b0 124, 125) at any row, no
 * erasure (b0 120..123, 126, 127), no tone frame and no rejected repeat gives
 * one from a fresh state; after a frame of row 20 (code 50) a rejected repeat,
 * whatever row it was coded with, holds the held word of the repeated frame,
 * not a word of its own bits, and the prediction history keeps the word as
 * decoded.
 */
static void
test_r33_rows_decoder_word(void) {
    for (int row = 0; row < 32; ++row) {
        for (int b0 = 0; b0 < 120; ++b0) {
            char bits[49];
            float out[160];
            mbe_parms cur, prev, enhanced, s[3];
            set_ambe2450(bits, b0, row);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            const int kind = mbe_decodeAmbe2450Parms(bits, &cur, &prev);
            assert(kind == MBE_AMBE2450_FRAME_VOICE);
            (void)kind;
            assert(mbe_frame_mark(&cur) == MBE_FRAME_MARK_NONE || row < 18);
            assert(holds_a_word(&cur.PHIl[0]) == (row >= 18));
            if (row >= 18) {
                float w0;
                int L;
                r33_grid(b0, &w0, &L);
                assert(holds_word(&cur.PHIl[0], row, b0, 0) && same_bits(&cur.w0, &w0, 1) && cur.L == L);
            }
            mbe_initMbeParms(&s[0], &s[1], &s[2]);
            (void)r33_frame(bits, out, s, -1);
            for (int i = 0; i < 3; ++i) {
                assert(holds_a_word(&s[i].PHIl[0]) == (row >= 18));
                assert(row < 18 || holds_word(&s[i].PHIl[0], row, b0, 0));
            }
        }
    }
    static const int specials[] = {120, 121, 122, 123, 124, 125, 126, 127};
    for (size_t k = 0; k < sizeof(specials) / sizeof(specials[0]); ++k) {
        for (int row = 0; row < 32; ++row) {
            char bits[49];
            float out[160];
            mbe_parms s[3];
            set_ambe2450(bits, specials[k], row);
            mbe_initMbeParms(&s[0], &s[1], &s[2]);
            (void)r33_frame(bits, out, s, -1);
            for (int i = 0; i < 3; ++i) {
                assert(!holds_a_word(&s[i].PHIl[0]));
            }
        }
    }
    char tone[49], rejected[49], word[49];
    float out[160];
    mbe_parms s[3];
    r33_tone_bits(tone, 130, 64);
    mbe_initMbeParms(&s[0], &s[1], &s[2]);
    const unsigned tone_flags = r33_frame(tone, out, s, -1);
    assert((tone_flags & MBE_PROCESS_FLAG_TONE) != 0u);
    (void)tone_flags;
    for (int i = 0; i < 3; ++i) {
        assert(!holds_a_word(&s[i].PHIl[0]));
    }
    set_ambe2450(rejected, 70, 25);
    mbe_initMbeParms(&s[0], &s[1], &s[2]);
    const unsigned rejected_flags = r33_frame(rejected, out, s, 4);
    assert((rejected_flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
    (void)rejected_flags;
    for (int i = 0; i < 3; ++i) {
        assert(!holds_a_word(&s[i].PHIl[0]));
    }
    set_ambe2450(word, 50, 20);
    mbe_initMbeParms(&s[0], &s[1], &s[2]);
    (void)r33_frame(word, out, s, 0);
    const unsigned repeat_flags = r33_frame(rejected, out, s, 4);
    assert((repeat_flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
    (void)repeat_flags;
    assert(holds_word(&s[2].PHIl[0], 20, 50, 1) && holds_word(&s[0].PHIl[0], 20, 50, 1));
    assert(holds_word(&s[1].PHIl[0], 20, 50, 0));
}

/* The stored pulsed-row word (slot 37) of a state, as test_row_word() bits,
 * or 0. */
static uint32_t
stored_row_bits(const mbe_parms* mp) {
    uint32_t bits;
    memcpy(&bits, &mp->noiseOverlap[TEST_ROW_SLOT], sizeof(bits));
    return bits;
}

/* A stream that reaches a frame of `row` at code b0 after voiced frames. */
static void
r33_stream_to_row(mbe_parms s[3], int row, int b0) {
    char voiced[49], frame[49];
    float out[160];
    set_ambe2450(voiced, 40, 0);
    set_ambe2450(frame, b0, row);
    mbe_initMbeParms(&s[0], &s[1], &s[2]);
    for (int f = 0; f < 3; ++f) {
        (void)r33_frame(voiced, out, s, 0);
    }
    (void)r33_frame(frame, out, s, 0);
}

/*
 * A lost frame's repeat through the data path, clause 7: after a frame of row
 * 19 (own subframe pulsed), 24 (both) or 27 (mid-frame only), two lost frames
 * (b0 121) and a voiced frame. A repeat of rows 19 and 24 holds the held word
 * (both repeats) and stores the word, not held, for the frame after (its own
 * model pulsed); one of row 27 clears the word and stores none. Against the
 * same stream with the word cleared from the last frame played before the lost
 * frames (no word to repeat): the repeats of rows 19 and 24 sound otherwise,
 * and so does the frame after them (the second repeat's own model's pulse
 * carried into it); the repeat of row 27 and the frame after it are the same
 * bit for bit, as are the states but for the word (its own mask pulses
 * nothing). The repeat's own pulsed bands are its own mask: tested on the
 * rendering itself in test_waveform_through_pulsed_rows().
 */
static void
test_r33_rows_repeat(void) {
    static const int rows[3] = {19, 24, 27};
    char lost[49], after[49];
    set_ambe2450(lost, 121, 0);
    set_ambe2450(after, 45, 0);
    for (int r = 0; r < 3; ++r) {
        const int row = rows[r], b0 = 82;
        const int own = row <= 25;
        mbe_parms s[3], plain[3];
        float out[3][160], plain_out[3][160];
        r33_stream_to_row(s, row, b0);
        assert(holds_word(&s[2].PHIl[0], row, b0, 0));
        assert(stored_row_bits(&s[2]) == (own ? test_row_word_bits(test_row_word(row, b0, 0)) : 0u));
        memcpy(plain, s, sizeof(plain));
        mbe_set_frame_mark(&plain[2], MBE_FRAME_MARK_NONE);
        for (int f = 0; f < 2; ++f) {
            const unsigned flags = r33_frame(lost, out[f], s, 0);
            const unsigned plain_flags = r33_frame(lost, plain_out[f], plain, 0);
            assert((flags & MBE_PROCESS_FLAG_REPEAT) != 0u && (plain_flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
            (void)flags;
            (void)plain_flags;
            assert(own ? holds_word(&s[2].PHIl[0], row, b0, 1) : !holds_a_word(&s[2].PHIl[0]));
            assert(stored_row_bits(&s[2]) == (own ? test_row_word_bits(test_row_word(row, b0, 0)) : 0u));
            assert(same_bits(out[f], plain_out[f], 160) == !own);
        }
        (void)r33_frame(after, out[2], s, 0);
        (void)r33_frame(after, plain_out[2], plain, 0);
        assert(same_bits(out[2], plain_out[2], 160) == !own);
        assert(stored_row_bits(&s[2]) == 0u);
        if (!own) {
            for (int i = 0; i < 3; ++i) {
                assert(same_state_but_mark(&s[i], &plain[i]));
            }
        }
    }
}

/*
 * The repeat under a tone's first frame loads the last frame as a lost frame's
 * repeat does: after a frame of row 20, the tone's first frame plays the
 * repeat with the held word under its first 80 samples. Against the same state
 * with the word cleared from the last frame played (no pulse), the tone frame
 * differs in its first 80 samples only, and there by 1 - n / 80 of the
 * difference between the repeat with the held word and the repeat without it
 * (each rendered here as the tone frame renders its outgoing speech, with the
 * library's mbe_repeat_load_model() and mbe_synthesizeRepeatedSpeechStartf(),
 * unlimited: what the test takes from the library); with the held word written
 * there by the caller it is the same bit for bit (the load gives the held form
 * either way). After a frame of row 27 (no own pulse) the three are the same.
 */
static void
test_r33_rows_under_tone(void) {
    static const int rows[2] = {20, 27};
    char tone[49];
    r33_tone_bits(tone, 130, 64);
    for (int r = 0; r < 2; ++r) {
        mbe_parms s[3], none[3], held[3];
        float out[160], none_out[160], held_out[160];
        r33_stream_to_row(s, rows[r], 30);
        memcpy(none, s, sizeof(none));
        memcpy(held, s, sizeof(held));
        mbe_set_frame_mark(&none[2], MBE_FRAME_MARK_NONE);
        held[2].PHIl[0] = (float)test_row_word(rows[r], 30, 1);
        /* The outgoing speech, with the held word and without any. */
        float speech[2][160];
        for (int k = 0; k < 2; ++k) {
            mbe_parms repeated = (k == 0) ? s[0] : none[0];
            mbe_parms last = (k == 0) ? s[2] : none[2];
            mbe_repeat_load_model(&repeated, &last);
            repeated.PHIl[0] = (k == 0 && rows[r] == 20) ? (float)test_row_word(20, 30, 1) : 0.0f;
            repeated.repeatCount = 0;
            mbe_synthesizeRepeatedSpeechStartf(speech[k], &repeated, &last, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450,
                                               MBE_SYNTHESIS_CONTINUES, MBE_OUTPUT_UNLIMITED);
        }
        const unsigned flags = r33_frame(tone, out, s, 0);
        assert((flags & MBE_PROCESS_FLAG_TONE) != 0u);
        (void)flags;
        (void)r33_frame(tone, none_out, none, 0);
        (void)r33_frame(tone, held_out, held, 0);
        assert(same_bits(out, held_out, 160));
        assert(same_bits(out + 80, none_out + 80, 80));
        assert(same_bits(out, none_out, 80) == (rows[r] != 20));
        double scale = 0.0;
        for (int n = 0; n < 80; ++n) {
            scale = fmax(scale, fmax(fabs((double)speech[0][n]), fabs((double)out[n])));
        }
        for (int n = 0; n < 80; ++n) {
            const double want = (1.0 - ((double)n / 80.0)) * ((double)speech[0][n] - (double)speech[1][n]);
            assert(fabs(((double)out[n] - (double)none_out[n]) - want) <= 1.0e-5 * scale + 1.0e-6);
        }
        assert(rows[r] != 20 || scale > 1.0e-3);
    }
}

/*
 * A frame advanced unrendered (mbe_advance_subframes(), as the synthesis core
 * advances a frame it mutes) stores the pulsed-row state rendering leaves:
 * after a plain frame, a frame of row 23 (own pulsed, code 60) rendered and the
 * same frame advanced, the second starting from a slot that says otherwise,
 * leave the same noiseOverlap, the word of row 23 among it, and the next frame
 * (row 23 again) plays the same samples after either, its last frame's own
 * pulse in them.
 */
static void
test_r33_rows_advance(void) {
    const enum mbe_voiced_synthesis r33 = MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    mbe_parms rendered, advanced, frame, plain;
    float out[160], after_rendered[160], after_advanced[160];
    float w0;
    int L;
    fresh_stream(&rendered);
    blank_model(&plain, 125.0, 30, 0, 8.0f);
    next_frame_as(out, &rendered, &plain, 0, r33);
    advanced = rendered;
    r33_grid(60, &w0, &L);
    blank_model(&frame, 125.0, L, 0, 8.0f);
    frame.w0 = w0;
    frame.PHIl[0] = (float)test_row_word(23, 60, 0);
    next_frame_with_word(out, &rendered, &frame, 0, r33);
    mbe_parms cur = with_model(&advanced, &frame);
    cur.PHIl[0] = frame.PHIl[0];
    cur.noiseOverlap[TEST_ROW_SLOT] = (float)test_row_word(19, 7, 0);
    mbe_advance_subframes(&cur, &advanced, 0, r33);
    advanced = cur;
    assert(stored_row_bits(&rendered) == test_row_word_bits(test_row_word(23, 60, 0)));
    assert(same_bits(rendered.noiseOverlap, advanced.noiseOverlap, 96));
    next_frame_with_word(after_rendered, &rendered, &frame, 0, r33);
    next_frame_with_word(after_advanced, &advanced, &frame, 0, r33);
    assert(frame_sounds(after_rendered) && same_bits(after_rendered, after_advanced, 160));
}

/*
 * Adaptive smoothing can voice a band of a frame decoded with none; such a
 * frame is no pulsed-row frame. A grid model of row 20 (code 40) carrying its
 * word, with one loud band and the error count that makes the smoothing voice
 * it, through the synthesis core: the frame it renders has a voiced band and
 * plays, and stores, exactly what the same frame without the word does, the
 * slot 0; muted (an error rate over the threshold) the core advances it voiced
 * all the same and stores 0 there.
 */
static void
test_r33_rows_smoothing(void) {
    mbe_parms prev, frame, plain, first, second;
    float out[160], plain_out[160];
    float w0;
    int L;
    mbe_initMbeParms(&prev, &first, &second);
    r33_grid(40, &w0, &L);
    blank_model(&frame, 125.0, L, 0, 4.0f);
    frame.w0 = w0;
    frame.Ml[5] = 4000.0f;
    frame.errorCountTotal = 5;
    frame.errorRate = 0.002f;
    frame.localEnergy = 1.0f;
    frame.mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
    plain = frame;
    frame.PHIl[0] = (float)test_row_word(20, 40, 0);
    first = frame;
    second = plain;
    mbe_parms prev_a = prev, prev_b = prev;
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &first, &prev_a, 1.0f, MBE_MUTE_NOISE_SPEC,
                                       MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_AMPLITUDE_LIMIT_NONE);
    mbe_synthesizeSpeechWithPreEnhRm0f(plain_out, &second, &prev_b, 1.0f, MBE_MUTE_NOISE_SPEC,
                                       MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_AMPLITUDE_LIMIT_NONE);
    assert(reference_voiced(&first) && frame_sounds(out));
    assert(same_bits(out, plain_out, 160) && same_state_but_mark(&first, &second));
    assert(stored_row_bits(&first) == 0u);
    /* Without the error count the same frame is pulsed and sounds otherwise. */
    mbe_parms clean = frame, clean_plain = plain, prev_c = prev, prev_d = prev;
    clean.errorCountTotal = 0;
    clean_plain.errorCountTotal = 0;
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &clean, &prev_c, 1.0f, MBE_MUTE_NOISE_SPEC,
                                       MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_AMPLITUDE_LIMIT_NONE);
    mbe_synthesizeSpeechWithPreEnhRm0f(plain_out, &clean_plain, &prev_d, 1.0f, MBE_MUTE_NOISE_SPEC,
                                       MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_AMPLITUDE_LIMIT_NONE);
    assert(!reference_voiced(&clean) && !same_bits(out, plain_out, 160));
    assert(stored_row_bits(&clean) == test_row_word_bits(test_row_word(20, 40, 0)));
    mbe_parms muted = frame, prev_e = prev;
    muted.errorRate = 0.2f;
    mbe_synthesizeSpeechWithPreEnhRm0f(out, &muted, &prev_e, 1.0f, MBE_MUTE_NOISE_SPEC,
                                       MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MBE_AMPLITUDE_LIMIT_NONE);
    assert(reference_voiced(&muted) && mbe_subframe_state_is_marked(&muted));
    assert(stored_row_bits(&muted) == 0u);
}

/*
 * Pulsed bands play no noise: two frames of row 24 (every channel pulsed in
 * both subframes, codes 30 and 75) after a plain frame, rendered with the
 * unvoiced part (a transform plan) and without it: the second plays the same
 * samples either way (its two models have no unvoiced band and the first left
 * no unvoiced overlap), while a plain frame does not.
 */
static void
test_r33_pulsed_bands_play_no_noise(void) {
    const enum mbe_voiced_synthesis r33 = MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    mbe_parms with_plan, without, plain;
    float a[160], b[160];
    fresh_stream(&with_plan);
    blank_model(&plain, 125.0, 30, 0, 8.0f);
    for (int f = 0; f < 3; ++f) {
        mbe_parms model = plain;
        float w0;
        int L;
        if (f > 0) {
            r33_grid((f == 1) ? 30 : 75, &w0, &L);
            blank_model(&model, 125.0, L, 0, 8.0f);
            model.w0 = w0;
            model.PHIl[0] = (float)test_row_word(24, (f == 1) ? 30 : 75, 0);
        }
        without = with_plan;
        mbe_parms cur_a = with_model(&with_plan, &model), cur_b = with_model(&without, &model);
        cur_a.PHIl[0] = model.PHIl[0];
        cur_b.PHIl[0] = model.PHIl[0];
        memset(a, 0, sizeof(a));
        memset(b, 0, sizeof(b));
        mbe_synthesize_subframes(a, &cur_a, &with_plan, plan, 1.0f, 0, r33);
        mbe_synthesize_subframes(b, &cur_b, &without, NULL, 1.0f, 0, r33);
        assert(same_bits(a, b, 160) == (f == 2));
        with_plan = cur_a;
    }
    mbe_fft_plan_free(plan);
}

/* The samples and the state a frame of `model` (its PHIl[0] as given) plays
 * after `last` in the rate-33 rendering. */
static void
r33_after(const mbe_parms* last, const mbe_parms* model, float out[160], mbe_parms* state) {
    *state = *last;
    next_frame_with_word(out, state, model, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450);
}

/* Float bit patterns of each class a caller's state can hold, beside the
 * words: k_slot_values (not-a-number, the infinities, subnormals, the zeros,
 * the largest float, small whole numbers and their neighbours) and these. */
static const uint32_t k_row_values[] = {
    0x44800000u, /* 1024.0f, the first word */
    0x447FFFFFu, /* just under it */
    0x44800001u, /* just over it */
    0x44801000u, /* 1024.5 */
    0xC4800000u, /* -1024.0 */
    0x447FC000u, /* 1023.0 */
    0x457F7000u, /* 4087.0, the last word (held row 25, code 119) */
    0x457F7001u, /* just over it */
    0x457F8000u, /* 4088.0: held row 25 with code 120, which is no code */
    0x45800000u, /* 4096.0: held row 26, code 0, which nothing writes */
    0x4597B800u, /* 4855.0: held row 31, code 119, which nothing writes */
    0x4597C000u, /* 4856.0 */
    0x459C4000u, /* 5000.0 */
    0x44F00000u, /* 1920.0: row 25's code 0 (1024 + 7 * 128) */
    0x44AF0000u, /* 1400.0: row 20's code 120 (1024 + 2 * 128 + 120) */
    0x44AFE000u, /* 1407.0: row 20's code 127 */
    0x45300000u, /* 2816.0: row 32's code 0 */
    0x46000000u, /* 8192.0 */
    0x4B000000u, /* 2^23 */
};

/*
 * Words in a caller's state, by their bit patterns. The reader
 * (mbe_pulsed_row_read()) takes exactly the 2,640 values a writer gives (rows
 * 18..31 at codes 0..119 not held, the decoder's 1,680; rows 18..25 held, the
 * repeat's 960) with their fields, and neither neighbour of any of them, nor a
 * held form of rows 26..31, nor any of k_slot_values and k_row_values that is
 * not one of them. On a grid model of row 21's code 40 after a plain frame,
 * each of those values in PHIl[0] plays and stores exactly what the frame with
 * no word does, or, a word of a code with the model's L (whose grid model is
 * the same: the model cannot tell the codes apart), plays otherwise than the
 * frame with no word and stores the word of its row (that each pair plays as
 * the rule says is the waveform test's); a word whose code has another L, a
 * word on the model with a voiced band, on a pitch an ulp off the grid's or on
 * the grid pitch with another harmonic count, a held word of rows 26..31 play
 * as no word. In slot 37 of the last frame's state, after a frame of row 21 at
 * code 40 (which stores its word), each value plays the next frame (row 16) as
 * with 0 there, but for a word as a frame stores it of a code with that
 * model's L, whose pulse sounds where its response reaches the frame; the
 * state each leaves is the one with 0. A held word, a word of rows 26..31, a
 * word of another L, the stored word with the last model an ulp off its grid
 * pitch and a state without the subframe mark are all as 0.
 */
static void
test_untrusted_row_words(void) {
    /* The reader. */
    for (int row = 18; row <= 31; ++row) {
        for (int b0 = 0; b0 < 120; ++b0) {
            for (int held = 0; held < 2; ++held) {
                const float value = (float)test_row_word(row, b0, held);
                const float below = nextafterf(value, -INFINITY), above = nextafterf(value, INFINITY);
                struct mbe_pulsed_row word = {0, 0, 0};
                /* A held word is a repeat's of a row with an own mask, 18..25. */
                const int is_word = !held || row <= 25;
                const int read = mbe_pulsed_row_read(&value, &word);
                assert(read == is_word);
                assert(!is_word || (word.row == row && word.b0 == b0 && word.held == held));
                const int read_below = mbe_pulsed_row_read(&below, &word);
                const int read_above = mbe_pulsed_row_read(&above, &word);
                assert(!read_below && !read_above);
                (void)read;
                (void)read_below;
                (void)read_above;
            }
        }
    }
    uint32_t values[sizeof(k_slot_values) / sizeof(k_slot_values[0]) + sizeof(k_row_values) / sizeof(k_row_values[0])];
    size_t count = 0;
    for (size_t i = 0; i < sizeof(k_slot_values) / sizeof(k_slot_values[0]); ++i) {
        values[count++] = k_slot_values[i];
    }
    for (size_t i = 0; i < sizeof(k_row_values) / sizeof(k_row_values[0]); ++i) {
        values[count++] = k_row_values[i];
    }
    for (size_t i = 0; i < count; ++i) {
        float value;
        struct mbe_pulsed_row word;
        memcpy(&value, &values[i], sizeof(value));
        const int is_word = values[i] == 0x44800000u || values[i] == 0x457F7000u || values[i] == 0x44F00000u;
        const int read = mbe_pulsed_row_read(&value, &word);
        assert(read == is_word);
        (void)read;
    }
    /* In PHIl[0]. */
    float w0;
    int L;
    r33_grid(40, &w0, &L);
    mbe_parms last, model, plain_state, state, canonical_state;
    float plain_out[160], out[160], canonical[160];
    fresh_stream(&last);
    blank_model(&model, 125.0, 30, 0, 8.0f);
    next_frame_as(out, &last, &model, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450);
    blank_model(&model, 125.0, L, 0, 8.0f);
    model.w0 = w0;
    r33_after(&last, &model, plain_out, &plain_state);
    const float own_word = (float)test_row_word(21, 40, 0);
    model.PHIl[0] = own_word;
    r33_after(&last, &model, out, &state);
    assert(!same_bits(out, plain_out, 160) && holds_word(&state.noiseOverlap[TEST_ROW_SLOT], 21, 40, 0));
    const float own_neighbours[2] = {nextafterf(own_word, -INFINITY), nextafterf(own_word, INFINITY)};
    for (size_t i = 0; i < count + 2 + 3360; ++i) {
        float value;
        if (i < count) {
            memcpy(&value, &values[i], sizeof(value));
        } else if (i < count + 2) {
            value = own_neighbours[i - count];
        } else {
            const int k = (int)(i - count - 2);
            value = (float)test_row_word(18 + (k / 240), (k / 2) % 120, k % 2);
        }
        struct mbe_pulsed_row word;
        int same_L = 0;
        if (mbe_pulsed_row_read(&value, &word)) {
            float other_w0;
            int other_L;
            r33_grid(word.b0, &other_w0, &other_L);
            same_L = other_L == L;
            assert(same_bits(&other_w0, &w0, 1) == same_L);
        }
        model.PHIl[0] = value;
        r33_after(&last, &model, out, &state);
        if (same_L) {
            /* The word of a code of the model's L is taken (the grid model is the
             * same) and plays its pulses, so the frame differs from the one with no
             * word, and the frame stores the word of a row with an own mask; that
             * each pair plays as the rule says, its times included, is checked
             * against the independent rendering in test_waveform_through_pulsed_rows(). */
            const uint32_t stored_want =
                (word.row <= 25) ? test_row_word_bits(test_row_word(word.row, word.b0, 0)) : 0u;
            assert(!same_bits(out, plain_out, 160) && stored_row_bits(&state) == stored_want);
            (void)stored_want;
        } else {
            assert(same_bits(out, plain_out, 160) && same_state_but_mark(&state, &plain_state));
        }
    }
    /* A word on a model it is not the word of (with a voiced band; a pitch an
     * ulp off the grid's; a code whose L differs), and a held word of a row
     * whose own mask pulses nothing. */
    float other_w0;
    int other_L;
    r33_grid(0, &other_w0, &other_L);
    assert(other_L != L);
    mbe_parms checks[5][2];
    for (int c = 0; c < 4; ++c) {
        checks[c][0] = model;
        checks[c][0].PHIl[0] = own_word;
    }
    checks[0][0].Vl[3] = 1;
    checks[1][0].w0 = nextafterf(w0, 1.0f);
    checks[2][0].PHIl[0] = (float)test_row_word(21, 0, 0);
    checks[3][0].PHIl[0] = (float)test_row_word(28, 40, 1);
    /* The grid pitch of code 40 with a harmonic count one above its own. */
    blank_model(&checks[4][0], 125.0, L + 1, 0, 8.0f);
    checks[4][0].w0 = w0;
    checks[4][0].PHIl[0] = own_word;
    for (int c = 0; c < 5; ++c) {
        float with_word[160], without_word[160];
        mbe_parms with_state, without_state;
        checks[c][1] = checks[c][0];
        mbe_set_frame_mark(&checks[c][1], MBE_FRAME_MARK_NONE);
        r33_after(&last, &checks[c][0], with_word, &with_state);
        r33_after(&last, &checks[c][1], without_word, &without_state);
        assert(same_bits(with_word, without_word, 160) && same_state_but_mark(&with_state, &without_state));
    }
    /* A word on a model whose one voiced band is harmonic 1, after a voiced frame at about four times
     * its pitch: the mid-frame model, at about twice it, has no voiced band (its harmonics fall on
     * the frame's even ones), so the word's mid-frame mask could pulse it; the frame is no pulsed-row
     * frame all the same, and plays as with no word. */
    int b40 = 0;
    float w40;
    int L40 = 0;
    while (L40 < 40) {
        r33_grid(++b40, &w40, &L40);
    }
    mbe_parms high, top, top_plain, high_state, top_state, top_plain_state;
    float high_out[160], top_out[160], top_plain_out[160];
    fresh_stream(&high_state);
    blank_model(&high, 125.0, 10, 10, 8.0f);
    high.w0 = 4.0f * w40;
    next_frame_as(high_out, &high_state, &high, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450);
    blank_model(&top, 125.0, L40, 1, 8.0f);
    top.w0 = w40;
    top_plain = top;
    top.PHIl[0] = (float)test_row_word(28, b40, 0);
    struct mbe_subframe_model mid_of_top;
    mbe_subframe_mid_model(&high_state, &top, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, &mid_of_top);
    assert(!mid_of_top.any_voiced && mid_of_top.L >= 10);
    r33_after(&high_state, &top, top_out, &top_state);
    r33_after(&high_state, &top_plain, top_plain_out, &top_plain_state);
    assert(same_bits(top_out, top_plain_out, 160) && same_state_but_mark(&top_state, &top_plain_state));
    /* In slot 37 of the last frame's state. */
    mbe_parms stored = state;
    model.PHIl[0] = own_word;
    r33_after(&last, &model, out, &stored);
    blank_model(&model, 125.0, 30, 0, 8.0f);
    mbe_parms zero = stored;
    zero.noiseOverlap[TEST_ROW_SLOT] = 0.0f;
    r33_after(&zero, &model, plain_out, &plain_state);
    r33_after(&stored, &model, out, &state);
    assert(!same_bits(out, plain_out, 160));
    for (size_t i = 0; i < count + 2 + 3360; ++i) {
        float value;
        if (i < count) {
            memcpy(&value, &values[i], sizeof(value));
        } else if (i < count + 2) {
            value = own_neighbours[i - count];
        } else {
            const int k = (int)(i - count - 2);
            value = (float)test_row_word(18 + (k / 240), (k / 2) % 120, k % 2);
        }
        struct mbe_pulsed_row word;
        int takes = 0;
        if (mbe_pulsed_row_read(&value, &word) && !word.held && word.row <= 25) {
            float other_pitch;
            int other_count;
            r33_grid(word.b0, &other_pitch, &other_count);
            takes = other_count == L;
        }
        mbe_parms forged = stored;
        forged.noiseOverlap[TEST_ROW_SLOT] = value;
        r33_after(&forged, &model, out, &state);
        /* A word taken plays the last frame's own pulse, at its time less 160,
         * through that model's bands: audible where its response reaches the
         * frame (its time plus the model's period and 16 samples past 0), and
         * leaving the stored state as it is. */
        const int sounds =
            takes && reference_row_time(word.row, word.b0, 0, 0) + (2.0 * M_PI / (double)w0) + 16.0 > 0.0;
        assert(same_bits(out, plain_out, 160) == !sounds && same_stream_state(&state, &plain_state));
        (void)sounds;
    }
    /* A word as a frame stores it, with the last frame's model a unit in the
     * last place off its code's grid pitch: no slot is read. */
    mbe_parms off = stored, off_zero = zero;
    off.w0 = nextafterf(stored.w0, 1.0f);
    off_zero.w0 = off.w0;
    r33_after(&off, &model, out, &state);
    r33_after(&off_zero, &model, canonical, &canonical_state);
    assert(same_bits(out, canonical, 160) && same_stream_state(&state, &canonical_state));
    /* Without the subframe mark: no slot is read. */
    mbe_parms unmarked = stored, unmarked_zero = zero;
    unmarked.noiseOverlap[MBE_SUBFRAME_STATE_SLOT] = 0.0f;
    unmarked_zero.noiseOverlap[MBE_SUBFRAME_STATE_SLOT] = 0.0f;
    r33_after(&unmarked, &model, out, &state);
    r33_after(&unmarked_zero, &model, plain_out, &plain_state);
    assert(same_bits(out, plain_out, 160) && same_stream_state(&state, &plain_state));
}

/*
 * The state another rendering wrote, read by rate 33's and the other way
 * round. A P25 frame at reserved code 48 (both patterns) leaves its emission,
 * read flag and counter (slots 34 to 36: 1.0, 1.0, 1.0) and 0 in slot 37; a
 * rate-33 frame after that state plays as after the same state with those
 * slots 0 (they are no row). A rate-33 frame of row 21 leaves its word in slot
 * 37; a D-STAR, P25 or ProVoice frame after that state (its rendering) plays
 * and stores as after the same state with slot 37 at 0, and stores 0 there.
 */
static void
test_row_slot_between_renderings(void) {
    static const enum mbe_voiced_synthesis others[3] = {MBE_VOICED_SYNTHESIS_DVSI, MBE_VOICED_SYNTHESIS_DVSI_IMBE,
                                                        MBE_VOICED_SYNTHESIS_DVSI_PROVOICE};

    mbe_parms fresh, p25, zeroed, a_state, b_state, model;
    float out[160], a[160], b[160];
    fresh_stream(&fresh);
    code48_after(&fresh, MBE_VOICED_SYNTHESIS_DVSI_IMBE, out, &p25);
    assert(p25.noiseOverlap[TEST_EMITTED_SLOT] > 0.0f || p25.noiseOverlap[TEST_READ_SLOT] > 0.0f);
    zeroed = p25;
    for (int slot = TEST_EMITTED_SLOT; slot <= TEST_ROW_SLOT; ++slot) {
        zeroed.noiseOverlap[slot] = 0.0f;
    }
    float w0;
    int L;
    r33_grid(40, &w0, &L);
    blank_model(&model, 125.0, L, 0, 8.0f);
    model.w0 = w0;
    model.PHIl[0] = (float)test_row_word(21, 40, 0);
    r33_after(&p25, &model, a, &a_state);
    r33_after(&zeroed, &model, b, &b_state);
    assert(same_bits(a, b, 160) && same_stream_state(&a_state, &b_state));
    mbe_parms r33_state;
    r33_after(&fresh, &model, out, &r33_state);
    assert(holds_word(&r33_state.noiseOverlap[TEST_ROW_SLOT], 21, 40, 0));
    mbe_parms r33_zero = r33_state;
    r33_zero.noiseOverlap[TEST_ROW_SLOT] = 0.0f;
    for (int i = 0; i < 3; ++i) {
        mbe_parms next;
        blank_model(&next, 125.0, 30, 0, 8.0f);
        a_state = r33_state;
        b_state = r33_zero;
        next_frame_as(a, &a_state, &next, 0, others[i]);
        next_frame_as(b, &b_state, &next, 0, others[i]);
        assert(same_bits(a, b, 160) && same_stream_state(&a_state, &b_state));
        assert(stored_row_bits(&a_state) == 0u);
    }
}

/* One frame of D-STAR's, P25's or ProVoice's path (row17_path_frame()), its
 * PHIl[0] in slot `slot` and noiseOverlap[37] of prev_mp_enhanced forged
 * first. */
static void
row_word_stream_tail(enum row17_path path, const mbe_parms s[3], int slot, float word, float row_slot,
                     const struct row17_frame* lost, const struct row17_frame* after, float out[320], mbe_parms t[3]) {
    memcpy(t, s, 3 * sizeof(mbe_parms));
    t[slot].PHIl[0] = word;
    t[2].noiseOverlap[TEST_ROW_SLOT] = row_slot;
    row17_path_frame(path, lost, out, t);
    row17_path_frame(path, after, out + 160, t);
}

/*
 * The word in the other codecs' paths, set in caller state: D-STAR's, P25's
 * and ProVoice's streams (their own entry points: the data paths and
 * mbe_processImbe7100x4400Framef()), with a pulsed-row word (held or not, of a
 * row with an own mask and of one without) in each of the three states'
 * PHIl[0], and the word in noiseOverlap[37] of the last frame played, then a
 * lost frame and a voiced frame: the samples are bit for bit those of the same
 * state with no mark (0) in that slot and 0 in noiseOverlap[37], and the
 * returned states alike in everything but the mark word, which travels where
 * the path copies the state (a lost frame's repeat loads the history from
 * prev_mp in these codecs): the word stands where the stream with no word
 * holds none, or where a decoder has since written its own mark.
 */
static void
test_row_words_in_codec_states(void) {
    static const enum row17_path paths[3] = {ROW17_DSTAR, ROW17_P25, ROW17_PROVOICE};

    static const int words[4][3] = {{20, 50, 0}, {20, 50, 1}, {28, 3, 0}, {31, 119, 1}};
    for (int p = 0; p < 3; ++p) {
        const enum row17_path path = paths[p];
        struct row17_frame voiced, unvoiced, after, lost;
        row17_frame_for(path, (path == ROW17_DSTAR) ? 45 : 86, (path == ROW17_DSTAR) ? 0xC : 0, "1111111100", &voiced);
        row17_frame_for(path, (path == ROW17_DSTAR) ? 76 : 150, 0, "", &unvoiced);
        row17_frame_for(path, (path == ROW17_DSTAR) ? 50 : 92, (path == ROW17_DSTAR) ? 0xF : 0, "1111111100", &after);
        row17_frame_for(path, (path == ROW17_DSTAR) ? 125 : 255, (path == ROW17_DSTAR) ? 0xF : 0, "", &lost);
        mbe_parms s[3];
        float out[160];
        mbe_initMbeParms(&s[0], &s[1], &s[2]);
        for (int f = 0; f < 3; ++f) {
            row17_path_frame(path, &voiced, out, s);
        }
        row17_path_frame(path, &unvoiced, out, s);
        for (int w = 0; w < 4; ++w) {
            const float word = (float)test_row_word(words[w][0], words[w][1], words[w][2]);
            for (int slot = 0; slot < 3; ++slot) {
                float plain[320], forged[320];
                mbe_parms plain_s[3], forged_s[3];
                mbe_parms base[3];
                memcpy(base, s, sizeof(base));
                row_word_stream_tail(path, base, slot, 0.0f, 0.0f, &lost, &after, plain, plain_s);
                row_word_stream_tail(path, base, slot, word, word, &lost, &after, forged, forged_s);
                assert(same_bits(plain, forged, 320));
                for (int i = 0; i < 3; ++i) {
                    assert(same_state_but_mark(&forged_s[i], &plain_s[i]));
                    assert(stored_row_bits(&forged_s[i]) == stored_row_bits(&plain_s[i]));
                }
            }
        }
    }
}

/*
 * Finding 1 of the review: the repeat loader reads the word only on the
 * band-grid model of its code. Through mbe_processAmbe2450Dataf(): a frame of
 * row 16 at b0 40 (a coded pitch, no grid model), a caller's value in the
 * enhanced model's PHIl[0], then an erasure (b0 121). For each of a word of
 * row 20 (1320.0f), one of row 28, a held word of row 20, and a word of code
 * 40 on the grid model of another harmonic count (code 0's, L 9), the three
 * returned PHIl[0] are the branch point's values, the value itself in cur_mp
 * and prev_mp_enhanced and 0 in prev_mp ([1320, 0, 1320] for the first), and
 * the samples and the states but for that word are those of the same stream
 * without the forged value.
 */
static void
test_r33_repeat_leaves_foreign_words(void) {
    const float values[4] = {(float)test_row_word(20, 40, 0), (float)test_row_word(28, 40, 0),
                             (float)test_row_word(20, 40, 1), (float)test_row_word(20, 40, 0)};
    char frame[49], lost[49];
    set_ambe2450(frame, 40, 16);
    set_ambe2450(lost, 121, 0);
    float grid0_w0;
    int grid0_L;
    r33_grid(0, &grid0_w0, &grid0_L);
    for (int c = 0; c < 4; ++c) {
        mbe_parms s[3], plain[3];
        float out[2][160];
        mbe_initMbeParms(&s[0], &s[1], &s[2]);
        (void)r33_frame(frame, out[0], s, 0);
        if (c == 3) {
            s[2].w0 = grid0_w0;
            s[2].L = grid0_L;
        }
        memcpy(plain, s, sizeof(plain));
        s[2].PHIl[0] = values[c];
        const unsigned flags = r33_frame(lost, out[0], s, -1);
        const unsigned plain_flags = r33_frame(lost, out[1], plain, -1);
        assert(flags == plain_flags && (flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
        (void)flags;
        (void)plain_flags;
        assert(same_bits(&s[0].PHIl[0], &values[c], 1) && same_bits(&s[2].PHIl[0], &values[c], 1));
        assert(test_row_word_bits(0) == 0u);
        uint32_t history;
        memcpy(&history, &s[1].PHIl[0], sizeof(history));
        assert(history == 0u);
        assert(same_bits(out[0], out[1], 160));
        for (int i = 0; i < 3; ++i) {
            assert(same_state_but_mark(&s[i], &plain[i]));
        }
    }
}

/*
 * The repeat loader tests the grid model (pitch by bit pattern, harmonic
 * count), not voicing: on the grid model of the word's code the word is read
 * whatever the voicing, as it is for a decoded frame whose adaptive smoothing
 * voiced a band. Through mbe_processAmbe2450Dataf(): a frame of row 16 at b0
 * 40, the enhanced model given code 40's grid pitch and harmonic count with
 * harmonic 1 voiced, a word of row 20 (1320.0f) or of row 28 in its PHIl[0],
 * then four erasures (three repeats and the repeat limit's mute). The returned
 * words are the library's (row 20: its held form, 3368.0f, in cur_mp and
 * prev_mp_enhanced, 0 in prev_mp; row 28: cleared, where the branch point
 * returned the value itself), and the samples, slot 37 and every other field
 * of the three states are those of the same stream without the word (the mute
 * noise's generator seeded alike for both streams).
 */
static void
test_r33_repeat_on_a_voiced_grid_model(void) {
    const float values[2] = {(float)test_row_word(20, 40, 0), (float)test_row_word(28, 40, 0)};
    const float returned[2] = {(float)test_row_word(20, 40, 1), 0.0f};
    char frame[49], lost[49];
    set_ambe2450(frame, 40, 16);
    set_ambe2450(lost, 121, 0);
    float w0;
    int L;
    r33_grid(40, &w0, &L);
    for (int c = 0; c < 2; ++c) {
        mbe_parms s[3], plain[3];
        float out[2][160];
        mbe_initMbeParms(&s[0], &s[1], &s[2]);
        (void)r33_frame(frame, out[0], s, 0);
        s[2].w0 = w0;
        s[2].L = L;
        s[2].Vl[1] = 1;
        memcpy(plain, s, sizeof(plain));
        s[2].PHIl[0] = values[c];
        for (int f = 0; f < 4; ++f) {
            mbe_setThreadRngSeed(12345u);
            const unsigned flags = r33_frame(lost, out[0], s, -1);
            mbe_setThreadRngSeed(12345u);
            const unsigned plain_flags = r33_frame(lost, out[1], plain, -1);
            assert(flags == plain_flags && (flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
            assert(f < 3 || (flags & MBE_PROCESS_FLAG_MUTE) != 0u);
            (void)flags;
            (void)plain_flags;
            assert(same_bits(&s[0].PHIl[0], &returned[c], 1) && same_bits(&s[2].PHIl[0], &returned[c], 1));
            uint32_t history;
            memcpy(&history, &s[1].PHIl[0], sizeof(history));
            assert(history == 0u);
            (void)history;
            assert(same_bits(out[0], out[1], 160));
            assert(same_bits(&s[2].noiseOverlap[TEST_ROW_SLOT], &plain[2].noiseOverlap[TEST_ROW_SLOT], 1));
            for (int i = 0; i < 3; ++i) {
                assert(same_state_but_mark(&s[i], &plain[i]));
            }
        }
    }
}

/*
 * What is not promised: a caller's state that holds a word on the band-grid
 * model of its code, or the word in slot 37 with such a last model, is to the
 * library the history of a frame of that row, which its bits never held. The
 * representation cannot tell the two apart (and two codes of one harmonic count
 * share a grid model). After a frame of row 16 at b0 40, the enhanced model set
 * to the grid model of code 40 with the word of row 20 at code 40: the erasure
 * after it repeats a pulsed-row frame (the held word in the returned states,
 * its samples other than the same stream's without the word: 159 of 160 in
 * dev-debug); the same grid model with no word but that word in slot 37: the
 * next frame of row 16 plays the last own model's pulse (14 samples in
 * dev-debug). Pinned as the limit of the representation.
 */
static void
test_r33_forged_history_plays(void) {
    char frame[49], lost[49];
    set_ambe2450(frame, 40, 16);
    set_ambe2450(lost, 121, 0);
    float w0;
    int L;
    r33_grid(40, &w0, &L);
    for (int slot = 0; slot < 2; ++slot) {
        mbe_parms s[3], plain[3];
        float out[2][160];
        mbe_initMbeParms(&s[0], &s[1], &s[2]);
        (void)r33_frame(frame, out[0], s, 0);
        s[2].w0 = w0;
        s[2].L = L;
        memcpy(plain, s, sizeof(plain));
        if (slot) {
            s[2].noiseOverlap[TEST_ROW_SLOT] = (float)test_row_word(20, 40, 0);
        } else {
            s[2].PHIl[0] = (float)test_row_word(20, 40, 0);
        }
        (void)r33_frame(slot ? frame : lost, out[0], s, -1);
        (void)r33_frame(slot ? frame : lost, out[1], plain, -1);
        int differ = 0;
        for (int n = 0; n < 160; ++n) {
            differ += !same_bits(&out[0][n], &out[1][n], 1);
        }
        assert(differ > 0);
        assert(slot || (holds_word(&s[0].PHIl[0], 20, 40, 1) && holds_word(&s[2].PHIl[0], 20, 40, 1)));
        (void)differ;
    }
}

/*
 * The cap on a pulsed band's magnitude. A grid model of row 24 at code 30
 * (every channel pulsed in both halves) with every Ml at 1e30, the largest
 * magnitude a model takes, after a frame of no amplitude: its own model's
 * pulsed bands take 1e30 divided by the unvoiced factor, over the cap, and so
 * the cap itself; the mid-frame model, between a silent frame and this one,
 * has no amplitude. The frame (no unvoiced part) is the independent rendering
 * of the own model's pulse at the capped magnitude
 * (reference_pulsed_magnitude(), which states the library's arithmetic), to a
 * thousandth of its peak; the models are the library's
 * (mbe_subframe_own_model()), the phases the test's own.
 */
static void
test_r33_pulsed_magnitude_cap(void) {
    float w0;
    int L;
    r33_grid(30, &w0, &L);
    mbe_parms last, model;
    float out[160];
    fresh_stream(&last);
    blank_model(&model, 125.0, L, 0, 1.0f);
    model.w0 = w0;
    for (int l = 1; l <= L; ++l) {
        model.Ml[l] = 1.0e30f;
        model.log2Ml[l] = (float)(0.3 * sin(0.9 * l));
    }
    model.PHIl[0] = (float)test_row_word(24, 30, 0);
    struct mbe_subframe_model own;
    mbe_parms cur = with_model(&last, &model);
    cur.PHIl[0] = model.PHIl[0];
    mbe_subframe_own_model(&cur, &own);
    int pulsed[57];
    const int any_pulsed = reference_row_bands(&own, 24, 0, pulsed);
    assert(any_pulsed);
    (void)any_pulsed;
    double phase[57], want[160] = {0.0};
    reference_own_phases(&cur, phase);
    /* A fresh stream's corner, 110 Hz, and rate 33's filter. */
    const struct reference_band band = {1, 110.0};
    reference_row_render(want, &own, &band, pulsed, phase, 2, reference_row_time(24, 30, 0, 2));
    next_frame_with_word(out, &last, &model, 0, MBE_VOICED_SYNTHESIS_DVSI_AMBE2450);
    double peak = 0.0;
    for (int n = 0; n < 160; ++n) {
        peak = fmax(peak, fabs(want[n]));
    }
    assert(peak > 1.0e30);
    for (int n = 0; n < 160; ++n) {
        assert(fabs((double)out[n] - want[n]) <= 1.0e-3 * peak);
    }
}

/*
 * The pulse times bit for bit, where the test can reproduce the library's
 * arithmetic (PULSE_REPLICA_EXACT; elsewhere the library and this test may
 * round apart, and a tolerance cannot tell the times in float from the times
 * in double). A replica of the library's float rendering of one
 * pulse (prepare_render_dvsi(), start_oscillators(), render_pulse() and
 * response_window() in mbe_subframe.c, op for op), fed the time the rule gives,
 * formed in double and rounded to float once, reproduces the library's frame
 * sample for sample. Taken from the library: the three models
 * (mbe_subframe_own_model(), mbe_subframe_mid_model()), their pulsed bands
 * chosen here as the rule says, and the helpers the rendering itself calls
 * (mbe_sincosf(), mbe_low_band_dvsi_gain(), mbe_low_band_dvsi_phase(),
 * mbe_low_cut_gain(), at a fresh stream's corner, 110 Hz: no frame of these
 * streams has a voiced band, so the tracker stays where it starts). Every row at every
 * code after a frame with no voiced band and no pulse, and the plain frame
 * after it, which plays the rest of an own pulse that falls past the frame;
 * every pair sounds in one of the two. Computing the grid rows' times in float,
 * which gives other floats at 132 of rows 24 and 25's 480 time evaluations (and
 * the single law's at 212 of its 360), moves samples here.
 */
struct exact_render {
    int count;
    float w0, period;
    float real[56], imag[56], step_c[56], step_s[56];
    short harmonic[56];
};

/* prepare_render_dvsi() op for op, at the corner `corner`. */
static void
exact_prepare(const struct mbe_subframe_model* model, const int pulsed[57], float corner, struct exact_render* r) {
    const float hz_per_harmonic = mbe_round_to_float(model->w0 * (8000.0f / (2.0f * (float)M_PI)));
    const float factor = mbe_round_to_float(MBE_UNVOICED_MAGNITUDE_FACTOR / sqrtf(model->w0));
    r->count = 0;
    r->w0 = model->w0;
    r->period = (2.0f * (float)M_PI) / model->w0;
    for (int l = 1; l <= model->L; l++) {
        const float w = model->w0 * (float)l;
        if (!pulsed[l] || !(w < (float)M_PI)) {
            continue;
        }
        const double magnitude = (double)model->magnitude[l] / (double)factor;
        const float m = (magnitude < (double)1.0e30f) ? (float)magnitude : 1.0e30f;
        const float hz = mbe_round_to_float(hz_per_harmonic * (float)l);
        const float cut = mbe_low_cut_gain(hz, corner);
        const float amplitude = 2.0f * m * 1.0f * mbe_low_band_dvsi_gain(hz) * cut;
        const int k = r->count++;
        float sn, cs;
        mbe_sincosf(model->phase[l] + mbe_low_band_dvsi_phase(hz), &sn, &cs);
        r->real[k] = amplitude * cs;
        r->imag[k] = amplitude * sn;
        mbe_sincosf(w, &r->step_s[k], &r->step_c[k]);
        r->harmonic[k] = (short)l;
    }
}

static float
exact_window(float t, float period) {
    float sn, cs;
    if (t < 16.0f) {
        mbe_sincosf((float)M_PI * (t / 16.0f), &sn, &cs);
        return 0.5f * (1.0f - cs);
    }
    if (t < period) {
        return 1.0f;
    }
    mbe_sincosf((float)M_PI * ((t - period) / 16.0f), &sn, &cs);
    return 0.5f * (1.0f + cs);
}

static void
exact_pulse(float* out, const struct exact_render* r, float pulse, int rising) {
    const float last = pulse + r->period + 16.0f;
    float real[56], imag[56];
    int first = (pulse > 0.0f) ? (int)ceilf(pulse) : 0;
    if (r->count == 0 || !(last > 0.0f) || first >= 80) {
        return;
    }
    float turn_s, turn_c, sn = 0.0f, cs = 1.0f;
    int harmonic = 0;
    const float t = (float)first - pulse;
    mbe_sincosf(r->w0 * (t - 16.0f), &turn_s, &turn_c);
    for (int k = 0; k < r->count; k++) {
        for (; harmonic < r->harmonic[k]; harmonic++) {
            const float next = (cs * turn_c) - (sn * turn_s);
            sn = (sn * turn_c) + (cs * turn_s);
            cs = next;
        }
        real[k] = (r->real[k] * cs) - (r->imag[k] * sn);
        imag[k] = (r->imag[k] * cs) + (r->real[k] * sn);
    }
    for (int n = first; n < 80 && (float)n < last; n++) {
        const float fade = (float)n * (1.0f / 80.0f);
        float sum = 0.0f;
        for (int k = 0; k < r->count; k++) {
            sum += real[k];
            const float next = (real[k] * r->step_c[k]) - (imag[k] * r->step_s[k]);
            imag[k] = (imag[k] * r->step_c[k]) + (real[k] * r->step_s[k]);
            real[k] = next;
        }
        out[n] += exact_window((float)n - pulse, r->period) * (rising ? fade : (1.0f - fade)) * sum;
    }
}

static void
test_r33_pulse_times_by_bits(void) {
    const enum mbe_voiced_synthesis r33 = MBE_VOICED_SYNTHESIS_DVSI_AMBE2450;
    mbe_parms last, plain;
    float out[160];
    fresh_stream(&last);
    blank_model(&plain, 125.0, 30, 0, 8.0f);
    next_frame_as(out, &last, &plain, 0, r33);
    int checked = 0, heard = 0;
    for (int row = 18; row <= 31; ++row) {
        for (int b0 = 0; b0 < 120; ++b0) {
            float w0;
            int L;
            mbe_parms model, state = last;
            r33_grid(b0, &w0, &L);
            blank_model(&model, 125.0, L, 0, 8.0f);
            model.w0 = w0;
            shape_model(&model, (size_t)b0);
            model.PHIl[0] = (float)test_row_word(row, b0, 0);
            mbe_parms cur = with_model(&last, &model);
            cur.PHIl[0] = model.PHIl[0];
            struct mbe_subframe_model models[3];
            mbe_subframe_mid_model(&last, &cur, 0, r33, &models[1]);
            mbe_subframe_own_model(&cur, &models[2]);
            int pulsed[3][57];
            const int any_mid = reference_row_bands(&models[1], row, 1, pulsed[1]);
            const int any_own = reference_row_bands(&models[2], row, 0, pulsed[2]);
            struct exact_render render[3];
            exact_prepare(&models[1], pulsed[1], MBE_LOW_CUT_START_HZ, &render[1]);
            exact_prepare(&models[2], pulsed[2], MBE_LOW_CUT_START_HZ, &render[2]);
            const float t1 = (float)reference_row_time(row, b0, 0, 1);
            const float t2 = (float)reference_row_time(row, b0, 0, 2);
            float want[160] = {0.0f};
            if (any_mid) {
                exact_pulse(want, &render[1], t1, 1);
                exact_pulse(want + 80, &render[1], t1 - 80.0f, 0);
            }
            if (any_own) {
                exact_pulse(want + 80, &render[2], t2 - 80.0f, 1);
            }
            /* The frame, then a plain frame after it, in which the frame's own
             * model plays the rest of its pulse falling, at its time less 160
             * (formed as every time is, in double and rounded once). */
            float next_want[160] = {0.0f}, next_out[160];
            if (any_own) {
                exact_pulse(next_want, &render[2], (float)reference_row_time(row, b0, 0, 0), 0);
            }
            next_frame_with_word(out, &state, &model, 0, r33);
            next_frame_as(next_out, &state, &plain, 0, r33);
            /* Elsewhere (fast-math, optimized x87, contracted multiply-adds)
             * the library may order or round the oscillators' sums otherwise:
             * within a thousandth of the pair's peak, a tolerance relative to
             * the parts, which the float-times mutation passes. */
            float peak = 0.0f;
            for (int n = 0; n < 160; ++n) {
                peak = fmaxf(peak, fmaxf(fabsf(want[n]), fabsf(next_want[n])));
            }
            for (int n = 0; n < 160; ++n) {
                if (PULSE_REPLICA_EXACT) {
                    assert(same_bits(&out[n], &want[n], 1) && same_bits(&next_out[n], &next_want[n], 1));
                } else {
                    assert(fabsf(out[n] - want[n]) <= 1.0e-3f * peak
                           && fabsf(next_out[n] - next_want[n]) <= 1.0e-3f * peak);
                }
            }
            checked += any_mid || any_own;
            heard += peak > 1.0f;
        }
    }
    assert(checked == 14 * 120 && heard == 14 * 120);
    (void)checked;
    (void)heard;
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
    test_waveform_through_changes(MBE_VOICED_SYNTHESIS_DVSI_PROVOICE, MID_PITCH_INEXACT);
    test_waveform_through_pulsed_codes(MBE_VOICED_SYNTHESIS_DVSI, MID_PITCH_INEXACT);
    test_waveform_through_pulsed_codes(MBE_VOICED_SYNTHESIS_DVSI_IMBE, MID_PITCH_INEXACT);
    test_waveform_through_pulsed_codes(MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MID_PITCH_INEXACT);
    test_waveform_through_pulsed_codes(MBE_VOICED_SYNTHESIS_DVSI_PROVOICE, MID_PITCH_INEXACT);
    test_pulsed_channel_edges(MID_PITCH_INEXACT);
    test_waveform_through_pulsed_rows(MBE_VOICED_SYNTHESIS_DVSI_AMBE2450, MID_PITCH_INEXACT);
    test_waveform_through_pulsed_rows(MBE_VOICED_SYNTHESIS_DVSI, MID_PITCH_INEXACT);
    test_waveform_through_pulsed_rows(MBE_VOICED_SYNTHESIS_DVSI_IMBE, MID_PITCH_INEXACT);
    test_waveform_through_pulsed_rows(MBE_VOICED_SYNTHESIS_DVSI_PROVOICE, MID_PITCH_INEXACT);
    test_octave_steps();
    test_stored_phase_range();
    test_unvoiced_pitch_carries();
    test_voicing_held_by_columns_without_a_pitch();
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
    test_dstar_mid_columns_by_words();
    test_dstar_mid_columns_table();
    test_p25_mid_columns_by_words();
    test_p25_mid_columns_cases();
    test_dstar_column_edges();
    test_train_through_held_columns();
    test_mid_columns_bands_with_no_state();
    test_marks_with_columns();
    test_marks_in_codec_states();
    test_p25_nearest_harmonic_across_a_step();
    test_p25_mid_columns_across_steps();
    test_p25_forged_tone_in_state();
    test_row17_mark_in_codec_states();
    test_dstar_exact_column_edge();
    test_rate33_hold_past_l();
    test_p25_muted_onset_advances();
    test_untrusted_mid_voicing_marks();
    test_p25_voicing_persists();
    test_p25_reset_drops_the_mark();
    test_p25_foreign_marks_in_state();
    test_marks_act_in_their_own_rendering();
    test_p25_rendering_marks_the_model_it_renders();
    test_p25_smoothed_onset();
    test_pulsed_codes_through_data_paths();
    test_p25_validity_counter_walk();
    test_p25_counter_through_a_mute();
    test_untrusted_pulsed_state(MBE_VOICED_SYNTHESIS_DVSI);
    test_untrusted_pulsed_state(MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    test_pulsed_count_steps();
    test_pulsed_read_flag();
    test_pulsed_renderings();
    test_ordinary_codes_carry_no_pattern();
    test_pulsed_state_is_the_models(MBE_VOICED_SYNTHESIS_DVSI);
    test_pulsed_state_is_the_models(MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    test_pulsed_held_by_a_forged_mark(MBE_VOICED_SYNTHESIS_DVSI);
    test_pulsed_held_by_a_forged_mark(MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    test_pulsed_advance_stores_the_emission(MBE_VOICED_SYNTHESIS_DVSI);
    test_pulsed_advance_stores_the_emission(MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    test_pulsed_smoothing_takes_the_code();
    test_pulsed_repeats_mutes_and_tones();
    test_provoice_has_no_pulsed_state();
    test_r33_rows_decoder_word();
    test_r33_rows_repeat();
    test_r33_rows_under_tone();
    test_r33_rows_advance();
    test_r33_rows_smoothing();
    test_r33_pulsed_bands_play_no_noise();
    test_untrusted_row_words();
    test_row_slot_between_renderings();
    test_row_words_in_codec_states();
    test_r33_repeat_leaves_foreign_words();
    test_r33_repeat_on_a_voiced_grid_model();
    test_r33_forged_history_plays();
    test_r33_pulsed_magnitude_cap();
    test_r33_pulse_times_by_bits();
    return 0;
}
