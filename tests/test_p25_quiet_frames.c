// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief P25's quiet frames: the quiet-frame gain and where P25's data path
 * applies it (docs/dvsi-parity.md, "P25: quiet frames").
 *
 * On DVSI's decodes of its vectors DVSI's P25 decoder plays a frame whose
 * decoded magnitudes are low lower still; the library multiplies every
 * magnitude a P25 frame renders by g = 2^clamp(K (x - X0), FLOOR, 0), x the
 * mean of the frame's decoded log2Ml[1..L], after adaptive smoothing, and
 * rounds P25's 16-bit output to nearest (test_frame_paths pins the 16-bit entry
 * points). The gain here is this test's own arithmetic (QUIET_K, QUIET_X0,
 * QUIET_FLOOR below, the documented constants), each step a stored binary64
 * operation, as the library's; the arithmetic tests (test_exact_arithmetic())
 * hold expected values written out as hexadecimal literals, computed outside
 * the test.
 * The frames are built with the IMBE encoder's quantizer
 * (mbe_imbe4400_quantize_amplitudes()) and played through the public P25 entry
 * points; the reference rendering ("replica") takes the model the data path
 * committed to the history (prev_mp after the frame) and renders it with the
 * library's enhancement (mbe_spectralAmpEnhanceWithRm0()), the library's
 * smoothing where a case has channel errors
 * (mbe_applyAdaptiveSmoothingWithRm0Limit()) and the library's synthesizer
 * without a gain (mbe_synthesizeSpeechWithPreEnhRm0f(), or
 * mbe_synthesizeRepeatedSpeechStartf() after the test's own smoothing), the
 * test's gain applied to the magnitudes in between. The decode is the
 * library's (mbe_decodeImbe4400Parms()).
 *
 * Where it compares exactly: the gain, by its bits, in every build (the
 * library's and this test's are stored binary64 steps), and the arithmetic
 * tests' gains and products in every build; the rendered magnitudes and the
 * samples of a scaled frame against the replica by their bits in builds
 * without fast-math whose float expressions are evaluated in float
 * (QUIET_EXACT: dev-debug, asan-ubsan-debug, notones-debug, the SIMD and the
 * x86-64-v3 Debug builds); elsewhere (dev-release, both Clang LTO builds, both
 * 32-bit builds) each magnitude within two float ulps (two steps of its bit
 * pattern) and each sample within 1e-5 of the frame's peak, a non-finite
 * magnitude only as the same bit pattern: the replica's enhancement and
 * synthesis are other compiled copies than the library's there. Those
 * tolerances cannot catch a magnitude or a sample an ulp off; the gain's own
 * mutations are caught by the exact gain comparisons in every build. A loud
 * frame, a frame of no gain (a tone model, ProVoice, a non-finite mean) and a
 * repeat's magnitudes against the frame it repeats are compared by their bits
 * in every build: they run the same code as their reference.
 */
#include <assert.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "imbe4400_internal.h"
#include "mbe_adaptive.h"
#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* The documented constants (docs/dvsi-parity.md, "P25: quiet frames"). */
#define QUIET_K     1.84
#define QUIET_X0    (-1.73)
#define QUIET_FLOOR (-2.00)

#if defined(MBELIB_TEST_LIBRARY_FAST_MATH) || defined(__FAST_MATH__) || defined(_M_FP_FAST)
#define QUIET_FAST_MATH 1
#else
#define QUIET_FAST_MATH 0
#endif
#if !QUIET_FAST_MATH && defined(FLT_EVAL_METHOD) && (FLT_EVAL_METHOD == 0)
#define QUIET_EXACT 1
#else
#define QUIET_EXACT 0
#endif

/* MBE_FRAME_MARK_TONE in PHIl[0], by its bits (4.0f). */
#define TONE_MARK_BITS 0x40800000u

static uint32_t
bits_of(float x) {
    uint32_t b;
    memcpy(&b, &x, sizeof(b));
    return b;
}

static uint64_t
bits_of_double(double x) {
    uint64_t b;
    memcpy(&b, &x, sizeof(b));
    return b;
}

/* Whether a gain is exactly 1, by its bits. */
static int
is_one(double g) {
    return bits_of_double(g) == bits_of_double(1.0);
}

static int
same_float(float a, float b) {
    return bits_of(a) == bits_of(b);
}

static int
same_floats(const float* a, const float* b, int n) {
    for (int i = 0; i < n; ++i) {
        if (!same_float(a[i], b[i])) {
            return 0;
        }
    }
    return 1;
}

static int
same_ints(const int* a, const int* b, int n) {
    for (int i = 0; i < n; ++i) {
        if (a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

/* Two states alike, every field by its bits. With uw_and_mute 0 two fields are left out, for the muted case alone:
 * previousUw, because a frame the core mutes is advanced, not rendered, and advancing leaves no unvoiced overlap
 * (mbe_advance_subframes(), by design: the muted frame's own subframe is never heard), while the replica renders
 * the frame and fills it; and mutingThreshold, because the replica sets it to 1e9 so that it renders the frame
 * the stream muted. */
static int
same_parms(const mbe_parms* a, const mbe_parms* b, int uw_and_mute) {
    return same_float(a->w0, b->w0) && a->L == b->L && a->K == b->K && same_ints(a->Vl, b->Vl, 57)
           && same_floats(a->Ml, b->Ml, 57) && same_floats(a->log2Ml, b->log2Ml, 57)
           && same_floats(a->PHIl, b->PHIl, 57) && same_floats(a->PSIl, b->PSIl, 57) && same_float(a->gamma, b->gamma)
           && a->tonePhase == b->tonePhase && a->swn == b->swn && same_float(a->localEnergy, b->localEnergy)
           && a->amplitudeThreshold == b->amplitudeThreshold && same_float(a->errorRate, b->errorRate)
           && a->errorCountTotal == b->errorCountTotal && a->errorCount4 == b->errorCount4
           && a->repeatCount == b->repeatCount && (!uw_and_mute || same_float(a->mutingThreshold, b->mutingThreshold))
           && (!uw_and_mute || same_floats(a->previousUw, b->previousUw, 256)) && same_float(a->noiseSeed, b->noiseSeed)
           && same_floats(a->noiseOverlap, b->noiseOverlap, 96);
}

static int
finite_float(float x) {
    return (bits_of(x) & 0x7FFFFFFFu) < 0x7F800000u;
}

/* The test's own gain: 1 for L outside 1..56 or a log2Ml in 1..L beyond 2^20 in
 * magnitude (a NaN or an infinity among them), or a mean at or above X0; else
 * 2^clamp(K (x - X0), FLOOR, 0), each step a stored binary64 operation (the
 * sum in order, harmonic 1 to L). */
/* 2^20, the largest |log2Ml| the gain reads, as a float bit pattern without its sign. */
#define QUIET_LOG_BOUND 0x49800000u

static double
own_gain(const mbe_parms* mp) {
    if (mp->L < 1 || mp->L > 56) {
        return 1.0;
    }
    volatile double sum = 0.0;
    for (int l = 1; l <= mp->L; ++l) {
        if ((bits_of(mp->log2Ml[l]) & 0x7FFFFFFFu) > QUIET_LOG_BOUND) {
            return 1.0;
        }
        sum = sum + (double)mp->log2Ml[l];
    }
    volatile double mean = sum / (double)mp->L;
    if (!(mean < QUIET_X0)) {
        return 1.0;
    }
    volatile double below = mean - QUIET_X0;
    volatile double exponent = QUIET_K * below;
    volatile double g = exp2((exponent < QUIET_FLOOR) ? QUIET_FLOOR : exponent);
    return g;
}

/* The test's own product: one stored binary64 multiplication, then float. */
static float
own_product(float magnitude, double gain) {
    volatile double product = (double)magnitude * gain;
    volatile float rounded = (float)product;
    return rounded;
}

static double
mean_log2(const mbe_parms* mp) {
    double sum = 0.0;
    for (int l = 1; l <= mp->L; ++l) {
        sum += (double)mp->log2Ml[l];
    }
    return sum / (double)mp->L;
}

/* The library's gain and this test's agree by their bits, in every build; and
 * the return says whether g < 1. */
static void
check_gain(const mbe_parms* mp, const char* what) {
    double gain = -1.0;
    const int quiet = mbe_imbe4400_quiet_gain(mp, &gain);
    const double want = own_gain(mp);
    if (bits_of_double(gain) != bits_of_double(want)) {
        fprintf(stderr, "%s: gain %.17g, want %.17g\n", what, gain, want);
    }
    assert(bits_of_double(gain) == bits_of_double(want));
    assert(quiet == (want < 1.0));
    if (!quiet) {
        assert(is_one(gain));
    }
}

/* A model of L harmonics, log2Ml[1..L] = x, the slots past L NaN (which must not be read). */
static mbe_parms
flat_model(int L, float x) {
    mbe_parms cur, prev, enh;
    mbe_initMbeParms(&cur, &prev, &enh);
    cur.L = L;
    for (int l = 1; l <= 56; ++l) {
        cur.log2Ml[l] = (l <= L) ? x : NAN;
    }
    return cur;
}

static double
gain_at(int L, float x) {
    const mbe_parms mp = flat_model(L, x);
    double gain = -1.0;
    (void)mbe_imbe4400_quiet_gain(&mp, &gain);
    check_gain(&mp, "flat model");
    return gain;
}

static float
float_step(float x, int steps) {
    uint32_t b = bits_of(x);
    /* x < 0 here, where a larger pattern is more negative: steps > 0 move toward 0, steps < 0 away from it. */
    b = (uint32_t)((int64_t)b - steps);
    float y;
    memcpy(&y, &b, sizeof(y));
    return y;
}

/*
 * The gain function: exactly 1 at and above the knee X0 (the magnitudes bit
 * for bit untouched), falling below it as 2^(K (x - X0)), the floor 2^FLOOR
 * exactly from X0 + FLOOR / K down (the lowest gain level the bits express,
 * B2[0] = -2.842205 and -2^20 among it), no gain for a log2 magnitude beyond
 * 2^20 in magnitude (a NaN or an infinity among them; -1e30 and -FLT_MAX), or
 * a harmonic count outside 1..56, the bound's two sides at 2^20 and the floats
 * beside it, and monotone and continuous across both knees (sweeps of 4,000
 * float neighbours each).
 */
static void
test_gain_function(void) {
    const double floor_gain = exp2(QUIET_FLOOR);
    const float loud[] = {0.0f, 5.0f, -1.3f, FLT_MAX};
    for (size_t i = 0; i < sizeof(loud) / sizeof(loud[0]); ++i) {
        const double g = gain_at(20, loud[i]);
        assert(is_one(g));
    }
    /* The two floats beside X0 (-1.73 is no float): the upper one no gain, the lower one a gain just under 1. */
    float upper = (float)QUIET_X0, lower = float_step(upper, -1);
    if ((double)upper < QUIET_X0) {
        lower = upper;
        upper = float_step(upper, 1);
    }
    const double at_upper = gain_at(20, upper);
    assert(is_one(at_upper));
    const double below = gain_at(20, lower);
    assert(below < 1.0 && below > 1.0 - 1e-6);
    /* Between the knees. */
    const double g2 = gain_at(20, -2.0f), g22 = gain_at(20, -2.4f);
    assert(fabs(g2 - exp2(QUIET_K * (-2.0 - QUIET_X0))) < 1e-6);
    assert(fabs(g22 - exp2(QUIET_K * (-2.4 - QUIET_X0))) < 1e-6);
    assert(g22 < g2 && g2 < 1.0 && g22 > floor_gain);
    /* At and below the lower knee: the floor exactly. */
    const float knee = (float)(QUIET_X0 + (QUIET_FLOOR / QUIET_K));
    const float under[] = {float_step(knee, -2), -2.9f, -2.842205f, -10.0f, -0x1p20f};
    for (size_t i = 0; i < sizeof(under) / sizeof(under[0]); ++i) {
        const double g = gain_at(20, under[i]);
        assert(bits_of_double(g) == UINT64_C(0x3FD0000000000000)); /* 1/4 */
    }
    /* Beyond the bound in every band: no gain, 1 exactly. */
    const float beyond[] = {float_step(-0x1p20f, -1), -1e30f, -FLT_MAX};
    for (size_t i = 0; i < sizeof(beyond) / sizeof(beyond[0]); ++i) {
        const double g = gain_at(20, beyond[i]);
        assert(is_one(g));
    }
    /* The bound's two sides in one band of a quiet model of 56 harmonics (the others -2.9), the pair +a, -b keeping
     * the mean under the knee: |a| and |b| at 2^20 and the float under it take the gain, the float over it none. */
    {
        const float sides[] = {float_step(-0x1p20f, 1), -0x1p20f, float_step(-0x1p20f, -1)};
        for (int k = 0; k < 3; ++k) {
            for (int sign = 0; sign < 2; ++sign) {
                mbe_parms mp = flat_model(56, -2.9f);
                mp.log2Ml[5] = sign ? -sides[k] : -0x1p20f;
                mp.log2Ml[9] = sign ? -0x1p20f : sides[k];
                check_gain(&mp, "the bound's sides");
                double g = -1.0;
                const int quiet = mbe_imbe4400_quiet_gain(&mp, &g);
                assert(quiet == (k < 2));
            }
        }
    }
    const double over_knee = gain_at(20, float_step(knee, 2));
    assert(over_knee > floor_gain);
    /* Harmonic counts: 1 and 56 take the gain, 0, 57 and -1 do not. */
    const double g1 = gain_at(1, -2.4f), g56 = gain_at(56, -2.4f);
    const double g0 = gain_at(0, -2.4f), g57 = gain_at(57, -2.4f), gm1 = gain_at(-1, -2.4f);
    assert(g1 < 1.0 && g56 < 1.0 && is_one(g0) && is_one(g57) && is_one(gm1));
    /* A mean that is not finite: NaN or either infinity in one band of 1..L; the slots past L are NaN in every
     * flat model above and are not read. */
    const float bad[] = {NAN, -NAN, INFINITY, -INFINITY};
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        for (int at = 1; at <= 20; at += 19) {
            mbe_parms mp = flat_model(20, -2.9f);
            mp.log2Ml[at] = bad[i];
            double g = -1.0;
            const int quiet = mbe_imbe4400_quiet_gain(&mp, &g);
            assert(quiet == 0 && is_one(g));
        }
    }
    /* Huge finite values, beyond the bound: a mix of +FLT_MAX and -FLT_MAX gives no gain; subnormals average to
     * about 0. */
    {
        mbe_parms mp = flat_model(20, FLT_MAX);
        for (int l = 1; l <= 10; ++l) {
            mp.log2Ml[l] = -FLT_MAX;
        }
        check_gain(&mp, "a +-FLT_MAX mix");
        assert(is_one(own_gain(&mp)));
        mp = flat_model(20, 1e-45f);
        check_gain(&mp, "subnormals");
        double g = -1.0;
        const int quiet = mbe_imbe4400_quiet_gain(&mp, &g);
        assert(quiet == 0);
        const int none = mbe_imbe4400_quiet_gain(NULL, &g);
        assert(none == 0 && is_one(g));
    }
    /* Every class of float in every harmonic 1..56 of a quiet model of 56 harmonics (the others -2.9, a gain at
     * the floor): NaNs of either sign and both infinities give no gain, exactly 1, with any other content;
     * zeros of either sign, subnormals, FLT_MIN, +-FLT_MAX and +-1e30 the test's gain (the last four beyond the
     * bound: none, 1 exactly; the others the gain of their mean). */
    {
        static const uint32_t none_classes[] = {0x7FC00000u, 0xFFC00000u, 0x7F800001u,
                                                0xFFFFFFFFu, 0x7F800000u, 0xFF800000u};
        static const uint32_t finite_classes[] = {0x00000000u, 0x80000000u, 0x00000001u, 0x807FFFFFu, 0x00800000u,
                                                  0x7F7FFFFFu, 0xFF7FFFFFu, 0x7149F2CAu, 0xF149F2CAu};
        int none = 0, finite = 0;
        for (int at = 1; at <= 56; ++at) {
            for (size_t i = 0; i < sizeof(none_classes) / sizeof(none_classes[0]); ++i) {
                mbe_parms mp = flat_model(56, -2.9f);
                memcpy(&mp.log2Ml[at], &none_classes[i], sizeof(float));
                double g = -1.0;
                const int quiet = mbe_imbe4400_quiet_gain(&mp, &g);
                assert(quiet == 0 && is_one(g));
                ++none;
            }
            for (size_t i = 0; i < sizeof(finite_classes) / sizeof(finite_classes[0]); ++i) {
                mbe_parms mp = flat_model(56, -2.9f);
                memcpy(&mp.log2Ml[at], &finite_classes[i], sizeof(float));
                check_gain(&mp, "a class of float in one harmonic");
                ++finite;
            }
        }
        printf("quiet-frame gain: %d states with a NaN or infinity in one harmonic give 1 exactly, %d with another "
               "class of float the test's gain (none beyond 2^20)\n",
               none, finite);
    }
    {
        const double g = gain_at(20, -0.0f);
        assert(is_one(g));
    }
    /* Monotone and continuous across both knees, float neighbour by float neighbour. */
    const float centres[] = {(float)QUIET_X0, knee};
    for (size_t c = 0; c < 2; ++c) {
        double last = gain_at(20, float_step(centres[c], -2000));
        for (int s = -1999; s <= 2000; ++s) {
            const double g = gain_at(20, float_step(centres[c], s));
            assert(g >= last);
            assert(g - last < 1e-6);
            last = g;
        }
    }
    printf("quiet-frame gain: 1 at and above X0, 2^(K (x - X0)) between (x -2.0: %.6f, -2.4: %.6f), the floor %.6f "
           "from x %.6f down; none beyond 2^20 in a band or L outside 1..56\n",
           g2, g22, floor_gain, (double)knee);
}

/* ---- P25 frames ---- */

static int
harmonics_of(int b0, int* K) {
    const float w0 = (float)(4 * M_PI) / (float)((float)b0 + 39.5);
    const int L = (int)(0.9254 * (int)((M_PI / w0) + 0.25));
    *K = (L < 37) ? (L + 2) / 3 : 12;
    return L;
}

/* A P25 voice frame at pitch code b0, voicing word `word` over the K bands (band 1 the most significant bit), log2
 * amplitude base + swing cos(0.4 l), quantized against the decoder's history. */
static void
build_frame(char imbe_d[88], int b0, unsigned word, float base, float swing, const mbe_parms* history) {
    int K;
    const int L = harmonics_of(b0, &K);
    unsigned char bands[12] = {0};
    float logm[57] = {0.0f};
    for (int band = 1; band <= K; ++band) {
        bands[band - 1] = (unsigned char)((word >> (K - band)) & 1u);
    }
    for (int l = 1; l <= L; ++l) {
        logm[l] = base + (swing * (float)cos(0.4 * (double)l));
    }
    mbe_imbe4400_quantize_amplitudes(b0, bands, logm, history, 0, imbe_d);
}

/* A lost frame: b0 255 (no pitch). */
static void
lost_frame(char imbe_d[88]) {
    static const int idx[8] = {0, 1, 2, 3, 4, 5, 85, 86};
    memset(imbe_d, 0, 88);
    for (int i = 0; i < 8; ++i) {
        imbe_d[idx[i]] = 1;
    }
}

/* play() with its status checked outside assert(). */
#define PLAY(s, d, pv, o, fl)                                                                                          \
    do {                                                                                                               \
        const int play_rc = play((s), (d), (pv), (o), (fl));                                                           \
        assert(play_rc >= 0);                                                                                          \
        (void)play_rc;                                                                                                 \
    } while (0)

struct stream {
    mbe_parms cur, prev, enh;
};

static void
stream_init(struct stream* s) {
    mbe_initMbeParms(&s->cur, &s->prev, &s->enh);
}

/* One frame through the public float data path; status and flags returned. */
static int
play(struct stream* s, const char imbe_d[88], int provoice, float out[160], unsigned* flags) {
    mbe_process_result r;
    mbe_initProcessResult(&r);
    if (provoice) {
        r.flags |= MBE_PROCESS_FLAG_PROVOICE;
    }
    const int rc = mbe_processImbe4400Dataf(out, &r, imbe_d, &s->cur, &s->prev, &s->enh);
    *flags = r.flags;
    return rc;
}

/*
 * The replica of the frame a stream just played: the model the data path
 * committed (`history`, the stream's prev_mp after the frame), enhanced by the
 * library, smoothed by the library where `smooth` (errors), times `gain` (1:
 * none), rendered by the library's synthesizer from `enh_before` (the stream's
 * prev_mp_enhanced before the frame). Leaves the rendered model in *rendered
 * and the new prev_mp_enhanced in *enh_after.
 */
static void
replica(float out[160], const mbe_parms* history, const mbe_parms* enh_before, double gain, int provoice, int smooth,
        mbe_parms* rendered, mbe_parms* enh_after) {
    mbe_parms cur = *history;
    *enh_after = *enh_before;
    const float rm0 = mbe_spectralAmpEnhanceWithRm0(&cur, MBE_HIGH_BAND_GAIN_IMBE);
    const enum mbe_voiced_synthesis rendering =
        provoice ? MBE_VOICED_SYNTHESIS_DVSI_PROVOICE : MBE_VOICED_SYNTHESIS_DVSI_IMBE;
    if (smooth) {
        mbe_applyAdaptiveSmoothingWithRm0Limit(&cur, enh_after, rm0, MBE_AMPLITUDE_LIMIT_NONE);
    }
    if (gain < 1.0) {
        for (int l = 1; l <= cur.L; ++l) {
            cur.Ml[l] = own_product(cur.Ml[l], gain);
        }
    }
    if (smooth) {
        /* Rendered, not muted, whatever the error rate: the muted case's state is compared with this one. */
        cur.mutingThreshold = 1e9f;
        mbe_synthesizeRepeatedSpeechStartf(out, &cur, enh_after, rendering, MBE_SYNTHESIS_CONTINUES,
                                           MBE_OUTPUT_LIMITED);
    } else {
        mbe_synthesizeSpeechWithPreEnhRm0f(out, &cur, enh_after, rm0,
                                           provoice ? MBE_MUTE_NOISE_COMFORT : MBE_MUTE_NOISE_SPEC, rendering,
                                           MBE_AMPLITUDE_LIMIT_NONE);
    }
    *rendered = cur;
    *enh_after = cur;
}

static double
peak(const float* x, int n) {
    double m = 0.0;
    for (int i = 0; i < n; ++i) {
        m = fmax(m, fabs((double)x[i]));
    }
    return m;
}

/* Samples equal by their bits (exact), or within 1e-5 of the frame's peak. */
static void
check_samples(const float* a, const float* b, int exact, const char* what) {
    if (exact) {
        for (int i = 0; i < 160; ++i) {
            if (!same_float(a[i], b[i])) {
                fprintf(stderr, "%s: sample %d %.9g against %.9g\n", what, i, (double)a[i], (double)b[i]);
            }
            assert(same_float(a[i], b[i]));
        }
        return;
    }
    /* Within 1e-5 of the larger of the two frames' peaks (the largest magnitude of the 160 samples); frames that
     * are silent must be equal. */
    const double bound = 1e-5 * fmax(peak(a, 160), peak(b, 160));
    for (int i = 0; i < 160; ++i) {
        /* A non-finite sample fails here whatever the other (inf <= inf would pass the bound). */
        assert(finite_float(a[i]) && finite_float(b[i]));
        assert(fabs((double)a[i] - (double)b[i]) <= bound);
    }
}

/* Float ulps between two finite floats of the same sign (the distance of their bit patterns); -1 otherwise. */
static int64_t
ulp_distance(float a, float b) {
    const uint32_t x = bits_of(a), y = bits_of(b);
    if (!finite_float(a) || !finite_float(b) || ((x ^ y) & 0x80000000u) != 0u) {
        return (x == y) ? 0 : -1;
    }
    return (x > y) ? (int64_t)(x - y) : (int64_t)(y - x);
}

/* Magnitudes 1..L equal by their bits (exact), or at most two ulps apart, two steps of the bit pattern of the
 * same sign (a non-finite one only as the same pattern, which a caller's NaN carried through both renderings is). */
static void
check_magnitudes_from(const mbe_parms* a, const mbe_parms* b, int exact, const mbe_parms* source, const char* what) {
    assert(a->L == b->L);
    /* A source model (the history a forged state repeats) holding a non-finite magnitude in 1..L: the enhancement
     * spreads it, and the two renderings must then hold the same bit patterns where either is not finite. */
    int source_nonfinite = 0;
    for (int l = 1; source != NULL && l <= source->L && l <= 56; ++l) {
        source_nonfinite |= !finite_float(source->Ml[l]);
    }
    for (int l = 1; l <= a->L; ++l) {
        if (exact) {
            if (!same_float(a->Ml[l], b->Ml[l])) {
                fprintf(stderr, "%s: Ml[%d] %.9g against %.9g\n", what, l, (double)a->Ml[l], (double)b->Ml[l]);
            }
            assert(same_float(a->Ml[l], b->Ml[l]));
        } else {
            /* Both finite first (ulp_distance() takes two identical non-finite patterns as 0), then at most two
             * float steps apart with the same sign; non-finite only from a non-finite source, then identical. */
            if (source_nonfinite && (!finite_float(a->Ml[l]) || !finite_float(b->Ml[l]))) {
                assert(same_float(a->Ml[l], b->Ml[l]));
                continue;
            }
            assert(finite_float(a->Ml[l]) && finite_float(b->Ml[l]));
            const int64_t d = ulp_distance(a->Ml[l], b->Ml[l]);
            assert(d >= 0 && d <= 2);
        }
    }
}

static void
check_magnitudes(const mbe_parms* a, const mbe_parms* b, int exact, const char* what) {
    check_magnitudes_from(a, b, exact, NULL, what);
}

/* The decoded model in the history: what the public parameter decoder gives for the frame's bits from the
 * history before it (w0, L, voicing, Ml and log2Ml by their bits). */
static void
check_history(const mbe_parms* history, const char imbe_d[88], const mbe_parms* prev_before) {
    mbe_parms cur = *prev_before, prev = *prev_before;
    const int rc = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
    assert(rc == 0);
    assert(history->L == cur.L && same_float(history->w0, cur.w0));
    for (int l = 1; l <= cur.L; ++l) {
        assert(history->Vl[l] == cur.Vl[l]);
        assert(same_float(history->Ml[l], cur.Ml[l]) && same_float(history->log2Ml[l], cur.log2Ml[l]));
    }
}

/* Every magnitude the frame rendered over the enhanced one it would render without the gain: g in every band,
 * voiced and unvoiced alike (to 1e-6), and how many of each. */
static void
check_ratio(const mbe_parms* rendered, const mbe_parms* history, double gain, int* voiced, int* unvoiced) {
    mbe_parms plain = *history;
    (void)mbe_spectralAmpEnhanceWithRm0(&plain, MBE_HIGH_BAND_GAIN_IMBE);
    for (int l = 1; l <= plain.L; ++l) {
        if (plain.Ml[l] > 0.0f) {
            assert(fabs(((double)rendered->Ml[l] / (double)plain.Ml[l]) - gain) <= 1e-6);
        }
        *voiced += rendered->Vl[l] == 1;
        *unvoiced += rendered->Vl[l] != 1;
    }
}

struct quiet_target {
    int b0;
    unsigned word;
    float base, swing;
};

/*
 * A P25 stream of quiet frames (at the lowest gain index and across the
 * curve, with and without voiced bands) and loud ones: each frame's samples
 * and rendered magnitudes against the replica with the test's gain (bits where
 * QUIET_EXACT, else the tolerances above; loud frames by their bits in every
 * build), every band's magnitude g times the unscaled one, the history the
 * decoder's model (Ml and log2Ml not scaled, so the next frame decodes from
 * the decoded values), and log2Ml of the rendered model the decoded one.
 */
static void
test_quiet_stream(void) {
    static const struct quiet_target targets[] = {
        {60, 0x000u, -4.0f, 0.0f},  {60, 0x000u, -3.2f, 0.6f},  {60, 0x7E0u, -3.0f, 0.8f},  {60, 0x000u, -2.7f, 0.4f},
        {100, 0x000u, -2.5f, 0.8f}, {100, 0xF00u, -2.4f, 0.0f}, {20, 0x00u, -2.2f, 0.8f},   {20, 0x3Fu, -2.1f, 0.5f},
        {60, 0x000u, -2.0f, 0.8f},  {60, 0xFFFu, -1.9f, 0.3f},  {100, 0x000u, -1.6f, 0.8f}, {60, 0x000u, 1.0f, 0.8f},
        {60, 0xFF0u, 2.0f, 0.5f},   {20, 0x3Fu, 3.0f, 0.0f},    {100, 0x000u, -3.5f, 0.8f},
    };
    struct stream s;
    stream_init(&s);
    int floor_frames = 0, slope_frames = 0, loud_frames = 0, voiced = 0, unvoiced = 0;
    for (int round = 0; round < 3; ++round) {
        for (size_t t = 0; t < sizeof(targets) / sizeof(targets[0]); ++t) {
            char imbe_d[88];
            float out[160], want[160];
            unsigned flags;
            const struct stream before = s;
            build_frame(imbe_d, targets[t].b0, targets[t].word, targets[t].base, targets[t].swing, &s.prev);
            PLAY(&s, imbe_d, 0, out, &flags);
            assert((flags & (MBE_PROCESS_FLAG_REPEAT | MBE_PROCESS_FLAG_MUTE)) == 0);
            check_history(&s.prev, imbe_d, &before.prev);
            const double g = own_gain(&s.prev);
            mbe_parms rendered, enh;
            replica(want, &s.prev, &before.enh, g, 0, 0, &rendered, &enh);
            const int exact = QUIET_EXACT || is_one(g);
            check_samples(out, want, exact, "quiet stream");
            check_magnitudes(&s.enh, &rendered, exact, "quiet stream");
            check_ratio(&s.enh, &s.prev, g, &voiced, &unvoiced);
            for (int l = 1; l <= s.prev.L; ++l) {
                assert(same_float(s.enh.log2Ml[l], s.prev.log2Ml[l]));
            }
            if (is_one(g)) {
                ++loud_frames;
                assert(same_parms(&s.enh, &enh, 1));
            } else if (g <= exp2(QUIET_FLOOR) * (1.0 + 1e-12)) {
                ++floor_frames;
            } else {
                ++slope_frames;
            }
        }
    }
    assert(floor_frames >= 6 && slope_frames >= 6 && loud_frames >= 6);
    printf("quiet P25 stream: %d frames at the floor, %d on the slope, %d loud ones untouched; %d voiced and %d "
           "unvoiced bands rendered (%s)\n",
           floor_frames, slope_frames, loud_frames, voiced, unvoiced, QUIET_EXACT ? "by bits" : "within tolerances");
}

/* A quiet frame then n lost frames; returns the stream. Checks each repeat against a second rendering of the
 * stored model (its magnitudes the first rendering's by their bits, its samples the replica's), the 4th muted, and
 * the reset model after it played without a gain. */
static void
test_repeats(void) {
    struct stream s;
    stream_init(&s);
    char imbe_d[88];
    float out[160];
    unsigned flags;
    for (int f = 0; f < 4; ++f) {
        build_frame(imbe_d, 60, 0x7E0u, -2.4f, 0.6f, &s.prev);
        PLAY(&s, imbe_d, 0, out, &flags);
    }
    const double g = own_gain(&s.prev);
    assert(g < 1.0);
    const mbe_parms stored = s.prev, first = s.enh;
    for (int lost = 1; lost <= 7; ++lost) {
        const struct stream before = s;
        lost_frame(imbe_d);
        PLAY(&s, imbe_d, 0, out, &flags);
        assert((flags & MBE_PROCESS_FLAG_REPEAT) != 0);
        if (lost <= MBE_MAX_FRAME_REPEATS) {
            /* The repeat loads the history: the stored model again, the gain read from it once. */
            assert(s.prev.L == stored.L && same_float(s.prev.w0, stored.w0));
            for (int l = 1; l <= stored.L; ++l) {
                assert(same_float(s.prev.Ml[l], stored.Ml[l]) && same_float(s.prev.log2Ml[l], stored.log2Ml[l]));
                assert(same_float(s.enh.Ml[l], first.Ml[l]));
            }
        } else {
            /* The repeat limit's reset model (and its repeats): log2Ml 0, no gain. */
            for (int l = 1; l <= s.prev.L; ++l) {
                assert(bits_of(s.prev.log2Ml[l]) == 0u);
            }
            assert(is_one(own_gain(&s.prev)));
        }
        const int muted = (flags & MBE_PROCESS_FLAG_MUTE) != 0;
        assert(muted == (lost == MBE_MAX_FRAME_REPEATS));
        mbe_parms rendered, enh;
        float want[160];
        const double gr = own_gain(&s.prev);
        replica(want, &s.prev, &before.enh, gr, 0, 0, &rendered, &enh);
        check_magnitudes(&s.enh, &rendered, QUIET_EXACT || is_one(gr), "repeat");
        if (!muted) {
            check_samples(out, want, QUIET_EXACT || is_one(gr), "repeat");
        }
    }
    printf("lost frames after a quiet frame (gain %.4f): repeats 1 to 3 play the stored model at its gain once, the "
           "4th is muted and advances it, the reset model after takes none\n",
           g);
}

/*
 * A quiet frame the error rate mutes (the history's error rate set high): the
 * state it leaves is the one rendering the scaled model leaves, its magnitudes
 * scaled, every field by its bits but two (same_parms() says which and why:
 * previousUw, which advancing a muted frame leaves empty by design, and the
 * muting threshold the replica raises to render the frame). A quiet frame whose smoothing voices
 * harmonics it decoded unvoiced (a moderate error rate, the local energy at
 * its lowest): the rendered magnitudes are what smoothing left, times the gain
 * read from the decoded log2Ml.
 */
/* A quiet frame at b0 60 with no voiced band, log2 amplitude -3.6 but at harmonic 5. */
static void
spiked_frame(char imbe_d[88], float spike, const mbe_parms* history) {
    int K;
    const int L = harmonics_of(60, &K);
    unsigned char bands[12] = {0};
    float logm[57] = {0.0f};
    for (int l = 1; l <= L; ++l) {
        logm[l] = (l == 5) ? spike : -3.6f;
    }
    mbe_imbe4400_quantize_amplitudes(60, bands, logm, history, 0, imbe_d);
}

#define SPIKE 9.0f

static void
test_errors(void) {
    for (int c = 0; c < 2; ++c) {
        struct stream s;
        stream_init(&s);
        char imbe_d[88];
        float out[160], want[160];
        unsigned flags;
        /* Ten frames of the target, so that the prediction carries the spike. */
        for (int f = 0; f < 10; ++f) {
            spiked_frame(imbe_d, (c == 0) ? -1.0f : SPIKE, &s.prev);
            PLAY(&s, imbe_d, 0, out, &flags);
        }
        s.prev.errorRate = (c == 0) ? 0.5f : 0.03f;
        s.enh.localEnergy = 10000.0f;
        const struct stream before = s;
        spiked_frame(imbe_d, (c == 0) ? -1.0f : SPIKE, &s.prev);
        PLAY(&s, imbe_d, 0, out, &flags);
        const double g = own_gain(&s.prev);
        assert(g < 1.0);
        mbe_parms rendered, enh;
        replica(want, &s.prev, &before.enh, g, 0, 1, &rendered, &enh);
        check_magnitudes(&s.enh, &rendered, QUIET_EXACT, "errors");
        int newly_voiced = 0;
        for (int l = 1; l <= s.prev.L; ++l) {
            assert(s.prev.Vl[l] == 0);
            newly_voiced += s.enh.Vl[l] == 1;
            assert(s.enh.Vl[l] == rendered.Vl[l]);
        }
        if (c == 0) {
            assert((flags & MBE_PROCESS_FLAG_MUTE) != 0);
            if (QUIET_EXACT) {
                assert(same_parms(&s.enh, &enh, 0));
            }
            printf("a quiet frame muted by the error rate advances the scaled model: its state is the rendered "
                   "case's%s\n",
                   QUIET_EXACT ? " byte for byte (previousUw aside)" : "'s magnitudes within tolerance");
        } else {
            assert((flags & MBE_PROCESS_FLAG_MUTE) == 0);
            assert(newly_voiced > 0);
            check_samples(out, want, QUIET_EXACT, "smoothed");
            printf("a quiet frame smoothing voices (%d of %d harmonics): what smoothing left, times the gain of the "
                   "decoded log2Ml\n",
                   newly_voiced, s.prev.L);
        }
    }
}

static double
droop(int l) {
    const double x = M_PI * (double)l / 256.0;
    return (sin(x) / x) * (sin(x) / x);
}

/* A DTMF tone (DTMF 1: b0 53, word 1110000, harmonics 4 and 7) sent as voice frames at low magnitudes, its mean
 * log2 magnitude under X0: it plays its tone model at its own level, no gain (the tone harmonics at 0.988 times
 * their decoded rms times the droop). */
static void
test_tone(void) {
    struct stream s;
    stream_init(&s);
    char imbe_d[88];
    float out[160];
    unsigned flags;
    int K;
    const int L = harmonics_of(53, &K);
    for (int f = 0; f < 10; ++f) {
        unsigned char bands[12] = {0};
        float logm[57] = {0.0f};
        for (int band = 1; band <= K; ++band) {
            bands[band - 1] = (unsigned char)((0x70u >> (K - band)) & 1u);
        }
        for (int l = 1; l <= L; ++l) {
            logm[l] = (l == 4 || l == 7) ? 1.0f : -3.4f;
        }
        mbe_imbe4400_quantize_amplitudes(53, bands, logm, &s.prev, 0, imbe_d);
        PLAY(&s, imbe_d, 0, out, &flags);
    }
    assert(bits_of(s.enh.PHIl[0]) == TONE_MARK_BITS);
    assert(own_gain(&s.prev) < 1.0);
    const double rms =
        sqrt(0.5 * ((double)s.prev.Ml[4] * (double)s.prev.Ml[4] + (double)s.prev.Ml[7] * (double)s.prev.Ml[7]));
    for (int l = 4; l <= 7; l += 3) {
        const double want = 0.988 * rms * droop(l);
        assert(fabs((double)s.enh.Ml[l] - want) <= 1e-6 * want);
    }
    printf("a P25 tone at a mean log2 magnitude of %.3f (gain %.4f were it speech) plays its tone model unscaled\n",
           mean_log2(&s.prev), own_gain(&s.prev));
}

/* The DTMF 1 tone stream of test_tone(), ten frames: the stream after it. */
static void
tone_stream(struct stream* s) {
    char imbe_d[88];
    float out[160];
    unsigned flags;
    int K;
    const int L = harmonics_of(53, &K);
    stream_init(s);
    for (int f = 0; f < 10; ++f) {
        unsigned char bands[12] = {0};
        float logm[57] = {0.0f};
        for (int band = 1; band <= K; ++band) {
            bands[band - 1] = (unsigned char)((0x70u >> (K - band)) & 1u);
        }
        for (int l = 1; l <= L; ++l) {
            logm[l] = (l == 4 || l == 7) ? 1.0f : -3.4f;
        }
        mbe_imbe4400_quantize_amplitudes(53, bands, logm, &s->prev, 0, imbe_d);
        PLAY(s, imbe_d, 0, out, &flags);
    }
}

/*
 * What decides "a tone model, no gain" against forged states: the quiet DTMF
 * stream's history (the repeat source) altered, then a lost frame, through the
 * public data path. Whatever the state, the frame takes no gain exactly when it
 * plays a tone model (the TONE mark in the model it rendered), and otherwise
 * the test's gain of the history it repeats, against the replica: the real
 * word (a tone model, at the tone's own level); the word removed, another
 * tone code's word, the word with its flag alone, the TONE mark alone (each
 * speech, at the gain of its mean); the tone harmonics' magnitudes at
 * impossible levels (1 and 100) and NaN, infinity, subnormal and negative-zero
 * magnitudes and log2 magnitudes, in one harmonic or in all.
 */
static void
test_tone_forgeries(void) {
    enum { CASES = 14 };

    static const char* const names[CASES] = {
        "the real word",    "no word",       "another code's word", "the flag alone",     "the TONE mark alone",
        "levels 1 and 100", "NaN log2Ml[1]", "NaN log2Ml in all",   "-inf log2Ml[2]",     "subnormal log2Ml",
        "-0 log2Ml",        "NaN Ml[4]",     "subnormal Ml[7]",     "-0 Ml in the others"};
    int tones = 0, speech = 0;
    for (int c = 0; c < CASES; ++c) {
        struct stream s;
        tone_stream(&s);
        assert(bits_of(s.enh.PHIl[0]) == TONE_MARK_BITS);
        mbe_parms* h = &s.prev;
        const int word = h->Vl[0];
        switch (c) {
            case 1: h->Vl[0] = 0; break;
            case 2: h->Vl[0] = (word & ~0x00FF0000) | (20 << 16); break;
            case 3: h->Vl[0] = word & ~0x00FFFFFF; break;
            case 4:
                h->Vl[0] = 0;
                h->PHIl[0] = 4.0f;
                break;
            case 5:
                h->Ml[4] = 1.0f;
                h->Ml[7] = 100.0f;
                break;
            case 6: h->log2Ml[1] = NAN; break;
            case 7:
                for (int l = 1; l <= h->L; ++l) {
                    h->log2Ml[l] = NAN;
                }
                break;
            case 8: h->log2Ml[2] = -INFINITY; break;
            case 9:
                for (int l = 1; l <= h->L; ++l) {
                    h->log2Ml[l] = 1e-45f;
                }
                break;
            case 10:
                for (int l = 1; l <= h->L; ++l) {
                    h->log2Ml[l] = -0.0f;
                }
                break;
            case 11: h->Ml[4] = NAN; break;
            case 12: h->Ml[7] = 1e-45f; break;
            case 13:
                for (int l = 1; l <= h->L; ++l) {
                    if (l != 4 && l != 7) {
                        h->Ml[l] = -0.0f;
                    }
                }
                break;
            default: break;
        }
        const struct stream before = s;
        char imbe_d[88];
        float out[160], want[160];
        unsigned flags;
        lost_frame(imbe_d);
        PLAY(&s, imbe_d, 0, out, &flags);
        for (int i = 0; i < 160; ++i) {
            assert(finite_float(out[i]));
        }
        const int tone = bits_of(s.enh.PHIl[0]) == TONE_MARK_BITS;
        if (tone) {
            /* The tone model at the tone's own level, from the magnitudes it repeats: no gain. */
            const double p4 = (double)s.prev.Ml[4] * (double)s.prev.Ml[4],
                         p7 = (double)s.prev.Ml[7] * (double)s.prev.Ml[7];
            for (int l = 4; l <= 7; l += 3) {
                const double want_l = 0.988 * sqrt(0.5 * (p4 + p7)) * droop(l);
                assert(fabs((double)s.enh.Ml[l] - want_l) <= 1e-6 * want_l);
            }
            ++tones;
        } else {
            const double g = own_gain(&s.prev);
            mbe_parms rendered, enh;
            replica(want, &s.prev, &before.enh, g, 0, 0, &rendered, &enh);
            check_samples(out, want, QUIET_EXACT || is_one(g), names[c]);
            check_magnitudes_from(&s.enh, &rendered, QUIET_EXACT || is_one(g), &s.prev, names[c]);
            ++speech;
        }
        /* The real word plays the tone; a word or mark that is not the decode's plays speech at the gain. */
        if (c == 0) {
            assert(tone);
        }
        if (c >= 1 && c <= 4) {
            assert(!tone && own_gain(&s.prev) < 1.0);
        }
        printf("  forged tone state, %s: %s\n", names[c],
               tone ? "a tone model, no gain" : "speech, the gain of its mean");
    }
    printf("forged tone states: %d play a tone model with no gain, %d speech at the gain of their mean\n", tones,
           speech);
}

/* A model of L harmonics at pitch w0, Ml[1..L] 1, no voiced band, log2Ml[1..L] from `logs` (the rest 0). */
static mbe_parms
arithmetic_model(int L, float w0, const float* logs, int n) {
    mbe_parms cur, prev, enh;
    mbe_initMbeParms(&cur, &prev, &enh);
    cur.L = L;
    cur.K = (L < 37) ? (L + 2) / 3 : 12;
    cur.w0 = w0;
    for (int l = 1; l <= 56; ++l) {
        cur.Vl[l] = 0;
        cur.Ml[l] = 1.0f;
        cur.log2Ml[l] = (l <= n) ? logs[l - 1] : 0.0f;
    }
    return cur;
}

static void
expect_gain(const mbe_parms* mp, int quiet_want, uint64_t gain_bits, const char* what) {
    double gain = -1.0;
    const int quiet = mbe_imbe4400_quiet_gain(mp, &gain);
    if (quiet != quiet_want || bits_of_double(gain) != gain_bits) {
        fprintf(stderr, "%s: flag %d gain %a, want %d and bits %llx\n", what, quiet, gain, quiet_want,
                (unsigned long long)gain_bits);
    }
    assert(quiet == quiet_want && bits_of_double(gain) == gain_bits);
}

/*
 * The gain's arithmetic, exact in every build (no tolerance): each expected
 * value below was computed outside the test, in Python's IEEE binary64 (the
 * gains through glibc's exp2), and is written as its bit pattern.
 * - The mean exactly on the knee: L 4, log2Ml -0x1.bae148p+2, 0x1.47ae14p-24,
 *   0x1p-49, 0 (summed in order, exactly 4 X0, the mean X0): no gain; L 51,
 *   -0x1.60eb86p+6, 0x1.c28f5cp-19, 0x1p-46 and 48 zeros (the mean X0 by its
 *   one rounded division, where K S / L - K X0, the order Clang's fast-math
 *   LTO builds once took, is -2^-51 and would take the gain): no gain.
 * - The doubles beside the knee: the next above (0x1.8p-49 in place of the
 *   last part), no gain; the next below (0x1p-50), the gain 0x1.ffffffffffffdp-1.
 * - L 55, forty-four -2 and eleven -1 (the mean -99/55): the gain
 *   0x1.d44562b64079ep-1.
 * - Cancelling histories beyond the bound, L 10: {2^58, -23, -2^58, 0, ...}
 *   (-32 in binary64, -23 in an extended register) and the review's
 *   {2^58, 32 + 2^-12, -2^58, -23, 0, ...} (41 in binary64, -23 under x87
 *   after double rounding, a gain of 1 against 0.483): no gain, 1 exactly, at
 *   the helper and through the public repeat path, where the rendered
 *   magnitudes equal those of the same model with every log2Ml 0 by their
 *   bits, and differ from those of the model with every log2Ml -3.2 (the
 *   floor).
 * - A history within the bound, {2^20, -2^-40, -2^20, -3.2 x 7}, whose
 *   binary64 sum loses the -2^-40 that an extended register keeps: the gain
 *   0x1.0b2ad7ecc9836p-1 (by its bits where FLT_EVAL_METHOD is 0, within one
 *   place under x87; a sum kept in an extended register gives ...9614p-1).
 * - A decoded frame on which x87 rounds twice: DVSI's `dtmf`, frame 3185 (L 14,
 *   its log2Ml as the decode stores them): the gain 0x1.32d283f8e55f5p-2 by
 *   its bits where FLT_EVAL_METHOD is 0 (the eight 64-bit configurations), and
 *   within one place of it where it is not (the two 32-bit x87 builds, which
 *   give 0x1.32d283f8e55f6p-2).
 * - The product: Ml 0x1.4affa2p+0 times the gain 0x1.2e561292a499p-1 through
 *   the scaled synthesis, 0x1.86e8dcp-1 (the binary64 product 0x1.86e8dbp-1,
 *   exactly halfway, rounded to even; one rounding of the exact product from
 *   an extended register gives 0x1.86e8dap-1).
 * - A mean of -2^20, the bound: the floor, and no overflow, underflow or
 *   invalid exception raised (the exponent is clamped before exp2(); exp2() of
 *   the unclamped exponent underflows); FE_INEXACT is not checked (the mean's
 *   and the exponent's roundings may raise it). A mean of -FLT_MAX: beyond the
 *   bound, no gain.
 */
static void
test_exact_arithmetic(void) {
    const float knee4[] = {-0x1.bae148p+2f, 0x1.47ae14p-24f, 0x1p-49f, 0.0f};
    const float above4[] = {-0x1.bae148p+2f, 0x1.47ae14p-24f, 0x1.8p-49f, 0.0f};
    const float below4[] = {-0x1.bae148p+2f, 0x1.47ae14p-24f, 0x1p-50f, 0.0f};
    const float knee51[] = {-0x1.60eb86p+6f, 0x1.c28f5cp-19f, 0x1p-46f};
    mbe_parms mp = arithmetic_model(4, 0.3f, knee4, 4);
    expect_gain(&mp, 0, UINT64_C(0x3FF0000000000000), "the mean on the knee, L 4");
    mp = arithmetic_model(4, 0.3f, above4, 4);
    expect_gain(&mp, 0, UINT64_C(0x3FF0000000000000), "the double above the knee");
    mp = arithmetic_model(4, 0.3f, below4, 4);
    expect_gain(&mp, 1, UINT64_C(0x3FEFFFFFFFFFFFFD), "the double below the knee");
    mp = arithmetic_model(51, 0.06f, knee51, 3);
    expect_gain(&mp, 0, UINT64_C(0x3FF0000000000000), "the mean on the knee, L 51");
    float logs55[55];
    for (int l = 0; l < 55; ++l) {
        logs55[l] = (l < 44) ? -2.0f : -1.0f;
    }
    mp = arithmetic_model(55, 0.05f, logs55, 55);
    expect_gain(&mp, 1, UINT64_C(0x3FED44562B64079E), "L 55, a mean of -99/55");
    const float cancel[] = {0x1p58f, -23.0f, -0x1p58f};
    const float review[] = {0x1p58f, 0x1.00008p+5f, -0x1p58f, -23.0f};
    mp = arithmetic_model(10, 0.28f, cancel, 3);
    expect_gain(&mp, 0, UINT64_C(0x3FF0000000000000), "the cancelling history");
    mp = arithmetic_model(10, 0.28f, review, 4);
    expect_gain(&mp, 0, UINT64_C(0x3FF0000000000000), "the review's cancelling history");
    {
        float bound[20], huge[20];
        for (int l = 0; l < 20; ++l) {
            bound[l] = -0x1p20f;
            huge[l] = -FLT_MAX;
        }
        mp = arithmetic_model(20, 0.15f, bound, 20);
#if defined(FE_UNDERFLOW) && defined(FE_OVERFLOW) && defined(FE_INVALID)
        feclearexcept(FE_ALL_EXCEPT);
        expect_gain(&mp, 1, UINT64_C(0x3FD0000000000000), "a mean of -2^20");
        assert(!fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID));
#else
        expect_gain(&mp, 1, UINT64_C(0x3FD0000000000000), "a mean of -2^20");
#endif
        mp = arithmetic_model(20, 0.15f, huge, 20);
        expect_gain(&mp, 0, UINT64_C(0x3FF0000000000000), "a mean of -FLT_MAX");
    }
    /* A history within the bound whose binary64 sum loses a term an extended register keeps: 2^20, -2^-40,
     * -2^20 and seven -3.2f (L 10). In binary64 the -2^-40 is lost in the first sum (the mean
     * -0x1.1eb8523333333p+1, the gain 0x1.0b2ad7ecc9836p-1); a sum kept in an 80-bit register keeps it (the
     * gain 0x1.0b2ad7ecc9614p-1). By its bits where FLT_EVAL_METHOD is 0, within one place under x87 (each
     * stored step is exact there, but exp2() is the platform's). */
    {
        float logs[10];
        logs[0] = 0x1p20f;
        logs[1] = -0x1p-40f;
        logs[2] = -0x1p20f;
        for (int l = 3; l < 10; ++l) {
            logs[l] = -3.2f;
        }
        mp = arithmetic_model(10, 0.28f, logs, 10);
        double g = -1.0;
        const int quiet = mbe_imbe4400_quiet_gain(&mp, &g);
        const uint64_t want = UINT64_C(0x3FE0B2AD7ECC9836);
        assert(quiet == 1);
#if defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD == 0
        assert(bits_of_double(g) == want);
#else
        assert(bits_of_double(g) == want || bits_of_double(g) == want + 1u || bits_of_double(g) == want - 1u);
#endif
    }
    /* DVSI's dtmf, frame 3185, as decoded (the stored floats by their bits). */
    {
        static const uint32_t dtmf3185[14] = {0xbf9b19ddu, 0xc054b2c2u, 0xc0281f94u, 0xc05e9ddcu, 0xc0250bf8u,
                                              0xc0549754u, 0xc0361568u, 0xc025a23au, 0xc048461eu, 0xc0068d26u,
                                              0xc0289808u, 0xc0197ab0u, 0xc051a41eu, 0xbff7c669u};
        float logs[14];
        memcpy(logs, dtmf3185, sizeof(logs));
        mp = arithmetic_model(14, 0.3f, logs, 14);
        double g = -1.0;
        const int quiet = mbe_imbe4400_quiet_gain(&mp, &g);
        const uint64_t want = UINT64_C(0x3FD32D283F8E55F5);
        assert(quiet == 1);
#if defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD == 0
        assert(bits_of_double(g) == want);
#else
        assert(bits_of_double(g) == want || bits_of_double(g) == want + 1u || bits_of_double(g) == want - 1u);
#endif
    }
    /* The cancelling histories through a public repeat. */
    {
        float flat0[10], flat32[10];
        for (int l = 0; l < 10; ++l) {
            flat0[l] = 0.0f;
            flat32[l] = -3.2f;
        }
        const float* histories[4] = {cancel, review, flat0, flat32};
        const int counts[4] = {3, 4, 10, 10};
        mbe_parms rendered[4];
        for (int h = 0; h < 4; ++h) {
            struct stream s;
            stream_init(&s);
            s.prev = arithmetic_model(10, 0.28f, histories[h], counts[h]);
            char imbe_d[88];
            float out[160];
            unsigned flags;
            lost_frame(imbe_d);
            PLAY(&s, imbe_d, 0, out, &flags);
            assert((flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
            for (int i = 0; i < 160; ++i) {
                assert(finite_float(out[i]));
            }
            rendered[h] = s.enh;
        }
        assert(same_floats(rendered[0].Ml, rendered[2].Ml, 57));
        assert(same_floats(rendered[1].Ml, rendered[2].Ml, 57));
        assert(!same_floats(rendered[0].Ml, rendered[3].Ml, 57));
    }
    /* The product, through the scaled synthesis (no smoothing change: the model has no error state). */
    {
        mbe_parms cur = arithmetic_model(10, 0.28f, cancel, 0);
        mbe_parms prev = cur;
        cur.Ml[1] = 0x1.4affa2p+0f;
        float out[160];
        uint64_t gb = UINT64_C(0x3FE2E561292A4990);
        double gain;
        memcpy(&gain, &gb, sizeof(gain));
        mbe_synthesizeSpeechScaledf(out, &cur, &prev, 10.0f, gain, MBE_MUTE_NOISE_SPEC, MBE_VOICED_SYNTHESIS_DVSI_IMBE,
                                    MBE_AMPLITUDE_LIMIT_NONE);
        assert(same_float(cur.Ml[1], 0x1.86e8dcp-1f));
        assert(same_float(own_product(0x1.4affa2p+0f, gain), 0x1.86e8dcp-1f));
    }
    printf("the gain's arithmetic: the knee, the doubles beside it, two cancelling histories beyond the bound (at the "
           "helper and through a public repeat), the product, a mean of -2^20 and of -FLT_MAX, by their bits; dtmf "
           "frame 3185%s\n",
           (FLT_EVAL_METHOD == 0) ? " by its bits" : " within one place");
}

/*
 * The model the data path returns and a scaled frame the error rate mutes:
 * after a quiet frame the returned cur_mp is the rendered (scaled) model,
 * every field the model left in prev_mp_enhanced; a quiet frame muted by the
 * error rate plays the P25 mute noise (TIA-102.BABA 7.8's +-5 on s(n), the
 * thread's generator seeded alike), not silence.
 */
static void
test_returned_and_muted(void) {
    struct stream s;
    stream_init(&s);
    char imbe_d[88];
    float out[160];
    unsigned flags;
    for (int f = 0; f < 4; ++f) {
        build_frame(imbe_d, 60, 0x000u, -2.6f, 0.5f, &s.prev);
        PLAY(&s, imbe_d, 0, out, &flags);
        assert(own_gain(&s.prev) < 1.0);
        assert(same_parms(&s.cur, &s.enh, 1));
    }
    s.prev.errorRate = 0.5f;
    build_frame(imbe_d, 60, 0x000u, -2.6f, 0.5f, &s.prev);
    mbe_setThreadRngSeed(0x0517u);
    PLAY(&s, imbe_d, 0, out, &flags);
    assert((flags & MBE_PROCESS_FLAG_MUTE) != 0u && own_gain(&s.prev) < 1.0);
    float want[160];
    mbe_setThreadRngSeed(0x0517u);
    mbe_synthesizeUniformNoisef(want, MBE_SPEC_MUTE_NOISE_AMPLITUDE * MBE_SPEECH_OUTPUT_GAIN);
    double energy = 0.0;
    for (int i = 0; i < 160; ++i) {
        energy += (double)out[i] * (double)out[i];
    }
    assert(energy > 0.0);
    check_samples(out, want, !QUIET_FAST_MATH, "the scaled frame's mute noise");
    printf("the returned model is the rendered one; a scaled frame the error rate mutes plays the mute noise\n");
}

/* ---- ProVoice frames (IMBE 7100x4400), packed from parameter bits ---- */

/* A Golay (23,12) word whose data bits 11..22 are set in `word`: the parity bits 0..10 for which the decoder
 * corrects nothing and gives the data back (as in test_subframe_synthesis.c). */
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

/* The ProVoice layout of parameter bits given in the 7200 layout (mbe_convertImbe7100to7200() run backwards). */
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

/* A ProVoice frame (7 x 24 bits, as received) that decodes with no error to parameter bits `t` (7200 layout). */
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
    const int rc = mbe_demodulateImbe7100x4400Data(fr);
    assert(rc == 0);
    (void)rc;
}

/*
 * ProVoice through its own frame entry points (mbe_processImbe7100x4400Framef()
 * and its 16-bit twin) on quiet frames, lost frames up to the reset model and
 * a frame the error rate mutes: every frame's float samples and states those
 * of the replica without a gain, by their bits in every build (the same code
 * path), and the 16-bit samples the float ones truncated; at least half of
 * the frames are ones P25 would play lower.
 */
static void
test_provoice_frames(void) {
    struct stream s;
    stream_init(&s);
    int lower = 0, frames = 0, muted = 0;
    for (int f = 0; f < 18; ++f) {
        char imbe_d[88], fr[7][24], scratch[88];
        float out[160], want[160];
        short pcm[160];
        if (f >= 8 && f < 14) {
            lost_frame(imbe_d);
        } else {
            build_frame(imbe_d, 60, (f % 2) ? 0x7E0u : 0x000u, -2.6f + 0.1f * (float)(f % 5), 0.6f, &s.prev);
        }
        if (f == 16) {
            s.prev.errorRate = 0.5f;
        }
        pack_imbe7100(imbe_d, fr);
        const struct stream before = s;
        struct stream twin = s;
        mbe_process_result r;
        mbe_initProcessResult(&r);
        /* The mute noise is the thread's generator's: seeded alike before each of the three renderings. */
        mbe_setThreadRngSeed(0x7100u + (unsigned)f);
        const int rc =
            mbe_processImbe7100x4400Framef(out, &r, (const char (*)[24])fr, scratch, &s.cur, &s.prev, &s.enh);
        assert(rc >= 0 && (r.flags & MBE_PROCESS_FLAG_PROVOICE) != 0u);
        mbe_process_result r2;
        mbe_initProcessResult(&r2);
        mbe_setThreadRngSeed(0x7100u + (unsigned)f);
        const int rc2 =
            mbe_processImbe7100x4400Frame(pcm, &r2, (const char (*)[24])fr, scratch, &twin.cur, &twin.prev, &twin.enh);
        assert(rc2 >= 0);
        mbe_parms rendered, enh;
        mbe_setThreadRngSeed(0x7100u + (unsigned)f);
        replica(want, &s.prev, &before.enh, 1.0, 1, 0, &rendered, &enh);
        assert(same_parms(&twin.enh, &s.enh, 1));
        /* Heard or muted (comfort noise), the samples, the 16-bit samples truncated and every field of the state. */
        check_samples(out, want, 1, "provoice frames");
        assert(same_parms(&s.enh, &enh, 1));
        short trunc[160];
        mbe_floattoshort(out, trunc);
        for (int i = 0; i < 160; ++i) {
            assert(trunc[i] == pcm[i]);
        }
        muted += (r.flags & MBE_PROCESS_FLAG_MUTE) != 0u;
        lower += own_gain(&s.prev) < 1.0;
        ++frames;
    }
    assert(2 * lower >= frames && muted >= 2);
    printf("ProVoice through its own frame entry points, %d frames (%d P25 would play lower, %d muted; repeats and the "
           "reset model among them): as without the gain\n",
           frames, lower, muted);
}

/* ProVoice through the P25 data path with its flag, on the same quiet frames, lost frames and an error-rate mute:
 * every frame's samples and states those of the replica without a gain, by their bits in every build. */
static void
test_provoice(void) {
    struct stream s;
    stream_init(&s);
    int frames = 0;
    for (int f = 0; f < 16; ++f) {
        char imbe_d[88];
        float out[160], want[160];
        unsigned flags;
        const struct stream before = s;
        if (f >= 8 && f < 13) {
            lost_frame(imbe_d);
        } else {
            build_frame(imbe_d, 60, (f % 2) ? 0x7E0u : 0x000u, -2.6f + 0.1f * (float)(f % 5), 0.6f, &s.prev);
        }
        if (f == 14) {
            s.prev.errorRate = 0.5f;
        }
        mbe_setThreadRngSeed(0x4400u + (unsigned)f);
        PLAY(&s, imbe_d, 1, out, &flags);
        mbe_parms rendered, enh;
        mbe_setThreadRngSeed(0x4400u + (unsigned)f);
        replica(want, &s.prev, &before.enh, 1.0, 1, 0, &rendered, &enh);
        /* Heard or muted (comfort noise, the generator seeded alike), the samples and every field of the state. */
        check_samples(out, want, 1, "provoice");
        assert(same_parms(&s.enh, &enh, 1));
        frames += own_gain(&s.prev) < 1.0;
    }
    assert(frames >= 8);
    printf("ProVoice on %d frames P25 would play lower (quiet, repeated, muted): as without the gain\n", frames);
}

/*
 * Caller states: a lost frame repeats the history, whose log2Ml a caller can
 * set to anything. NaN, either infinity in one band or in all, subnormals,
 * +FLT_MAX or -FLT_MAX in all bands or in one, -2^20 in all bands; values past
 * L (NaN) that must not be read; L out of range. The samples are finite; where
 * a log2Ml lies beyond 2^20 (and L out of range) they are the replica's
 * without a gain by their bits; -2^20 plays at the floor; NaN past L changes
 * nothing.
 */
static void
test_caller_states(void) {
    static const struct {
        const char* name;
        float value;
        int all;
    } cases[] = {{"NaN in band 3", NAN, 0},     {"+inf in band 3", INFINITY, 0},     {"-inf in band 3", -INFINITY, 0},
                 {"NaN in every band", NAN, 1}, {"subnormals", 1e-45f, 1},           {"+FLT_MAX", FLT_MAX, 1},
                 {"-FLT_MAX", -FLT_MAX, 1},     {"-FLT_MAX in band 3", -FLT_MAX, 0}, {"-2^20", -0x1p20f, 1}};

    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]) + 2; ++c) {
        struct stream s;
        stream_init(&s);
        char imbe_d[88];
        float out[160], want[160];
        unsigned flags;
        for (int f = 0; f < 3; ++f) {
            build_frame(imbe_d, 60, 0x000u, -2.4f, 0.6f, &s.prev);
            PLAY(&s, imbe_d, 0, out, &flags);
        }
        if (c < sizeof(cases) / sizeof(cases[0])) {
            for (int l = 1; l <= s.prev.L; ++l) {
                if (cases[c].all || l == 3) {
                    s.prev.log2Ml[l] = cases[c].value;
                }
            }
        } else if (c == sizeof(cases) / sizeof(cases[0])) {
            s.prev.L = 57;
        } else {
            s.prev.L = 0;
        }
        struct stream past = s;
        for (int l = s.prev.L + 1; l <= 56 && s.prev.L >= 0; ++l) {
            past.prev.log2Ml[l] = NAN;
        }
        const struct stream before = s;
        lost_frame(imbe_d);
        PLAY(&s, imbe_d, 0, out, &flags);
        float out_past[160];
        PLAY(&past, imbe_d, 0, out_past, &flags);
        check_samples(out, out_past, 1, "NaN past L");
        for (int i = 0; i < 160; ++i) {
            assert(finite_float(out[i]));
        }
        const double g = own_gain(&s.prev);
        if (s.prev.L >= 1 && s.prev.L <= 56) {
            mbe_parms rendered, enh;
            replica(want, &s.prev, &before.enh, g, 0, 0, &rendered, &enh);
            check_samples(out, want, QUIET_EXACT || is_one(g), "caller state");
            if (c == 6 || c == 7) {
                assert(is_one(g));
            }
            if (c == 8) {
                assert(bits_of_double(g) == UINT64_C(0x3FD0000000000000));
            }
        }
    }
    printf("caller states (non-finite, subnormal and huge log2Ml, NaN past L, L 0 and 57): finite samples, no gain "
           "beyond 2^20, the floor for -2^20\n");
}

/* The public mbe_synthesizeSpeechf() takes no gain: a model of quiet log2Ml plays as the same model with log2Ml 0. */
static void
test_public_synthesizer(void) {
    mbe_parms a, b, pa, pb, e;
    mbe_initMbeParms(&a, &pa, &e);
    a.w0 = 0.105f;
    a.L = 27;
    for (int l = 1; l <= a.L; ++l) {
        a.Vl[l] = (l % 3) == 0;
        a.Ml[l] = 0.12f + 0.01f * (float)(l % 5);
        a.log2Ml[l] = -2.8f;
    }
    pa = a;
    b = a;
    for (int l = 1; l <= b.L; ++l) {
        b.log2Ml[l] = 0.0f;
    }
    pb = pa;
    for (int l = 1; l <= pb.L; ++l) {
        pb.log2Ml[l] = 0.0f;
    }
    float oa[160], ob[160];
    mbe_setThreadRngSeed(0x5157u);
    mbe_synthesizeSpeechf(oa, &a, &pa);
    mbe_setThreadRngSeed(0x5157u);
    mbe_synthesizeSpeechf(ob, &b, &pb);
    check_samples(oa, ob, 1, "public synthesizer");
    for (int l = 1; l <= a.L; ++l) {
        assert(same_float(a.Ml[l], b.Ml[l]));
    }
    printf("mbe_synthesizeSpeechf(): a quiet model plays as the same model with log2Ml 0\n");
}

int
main(void) {
    test_gain_function();
    test_exact_arithmetic();
    test_quiet_stream();
    test_repeats();
    test_errors();
    test_tone();
    test_provoice();
    test_provoice_frames();
    test_returned_and_muted();
    test_tone_forgeries();
    test_caller_states();
    test_public_synthesizer();
    printf("P25 quiet frames: ok\n");
    return 0;
}
