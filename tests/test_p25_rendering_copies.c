// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief P25's copies of the shared mid-frame and own-model code against the
 * shared code, on the behaviour they share.
 *
 * P25's rendering has its own copies of the mid-frame filler and of the column
 * helpers (p25_fill_mid_bands(), p25_column_band_state(), p25_columns_voice()
 * and p25_own_model() in src/core/mbe_subframe.c), so that the other
 * renderings' code keeps its inlining: one body serving both moved D-STAR's
 * and rate 33's samples in the 32-bit Release build. The copies differ from
 * the shared code only where P25's rules do (the phase step on unvoiced bands,
 * each frame's state as its own model plays it by the columns). This test
 * compiles that file's static functions here (the library's own object is
 * not linked: every symbol it defines is defined by this translation unit)
 * and checks, over every ordered pair of the 208 P25 pitch codes, each pair
 * with one of 64 pairs of voicing classes (every class pair on many code
 * pairs), that the copies give what the shared code gives where P25's rules
 * do not apply (band states exactly, floats by their bits; in an optimized
 * build with excess precision, the 32-bit Release one, or with fast-math, the
 * two bodies may round apart, and floats there are compared within 1e-5;
 * 1.34e-6 is the largest difference seen, in the 32-bit Release build): the
 * column tests on P25's samples, each band's
 * mid-frame state with neither frame's own model read by the columns, the
 * mid-frame bands with no unvoiced band carrying a finite log2Ml (no step)
 * and no column rule (or one with neither frame's own model read by it), and
 * the own model with no step, the columns applied to the shared one. It also
 * runs the floored train's grid test on every point of the grid, 1 to 2^18
 * steps, and on both floats beside each (check_grid_points()), and the
 * phase a floored interval stores at a whole cycle, then an octave step
 * (check_whole_cycle_then_octave()). What it
 * does not cover: it compiles mbe_subframe.c into its own translation unit,
 * with the test's flags, so in a library built with fast-math it compares the
 * two bodies as the test's flags compile them, not the library's object (the
 * differential harness and the waveform tests see the library's).
 */
#include "../src/core/mbe_subframe.c" // NOLINT(bugprone-suspicious-include): the test compiles the library's static functions

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbe_adaptive.h"
#include "mbe_math.h"
#include "mbe_subframe.h"
#include "mbelib-neo/mbelib.h"

/* Whether two floats have the same bits. */
static int
same_float(const float* a, const float* b) {
    uint32_t x, y;
    memcpy(&x, a, sizeof(x));
    memcpy(&y, b, sizeof(y));
    return x == y;
}

/* Whether the two copies, compiled here as separate functions, can round
 * their float arithmetic apart: under excess precision in an optimized build
 * (x87: which intermediates stay in 80-bit registers depends on each body's
 * register allocation), and where the build may reorder (fast-math). */
#if !defined(FLT_EVAL_METHOD) || (FLT_EVAL_METHOD != 0 && defined(__OPTIMIZE__)) || defined(__FAST_MATH__)
#define COPIES_INEXACT 1
#else
#define COPIES_INEXACT 0
#endif

/* Whether a float is finite, by its bit pattern. */
static int
finite_float(const float* x) {
    uint32_t bits;
    memcpy(&bits, x, sizeof(bits));
    return (bits & 0x7FFFFFFFu) < 0x7F800000u;
}

/* How far two floats lie apart, relative to the larger magnitude (or to 1):
 * +infinity where either is not finite, identical bits included (so that
 * neither mode accepts an infinity or a NaN); 0 for the same bits; otherwise
 * at least DBL_MIN, so that two different bit patterns at no distance (+0
 * and -0) fail the exact mode and pass only a tolerance. */
static double
float_apart(const float* a, const float* b) {
    if (!finite_float(a) || !finite_float(b)) {
        return INFINITY;
    }
    if (same_float(a, b)) {
        return 0.0;
    }
    const double d = fabs((double)*a - (double)*b);
    return fmax(DBL_MIN, d / fmax(1.0, fmax(fabs((double)*a), fabs((double)*b))));
}

/* The comparator's own cases: an infinity or a NaN against its own bits fails
 * both modes, +0 against -0 fails the exact mode, the same finite bits pass. */
static void
check_comparator(void) {
    static const uint32_t cases[][2] = {
        {0x7F800000u, 0x7F800000u}, {0xFF800000u, 0xFF800000u}, {0x7FC01234u, 0x7FC01234u}, {0x00000000u, 0x80000000u}};
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        float a, b;
        memcpy(&a, &cases[i][0], sizeof(a));
        memcpy(&b, &cases[i][1], sizeof(b));
        const double apart = float_apart(&a, &b);
        assert(apart > 0.0);
        assert((i < 3) ? !(apart <= 1.0e-5) : (apart <= 1.0e-5));
    }
    const float one = 1.0f;
    const double same = float_apart(&one, &one);
    assert(!(same > 0.0));
}

/* How far two subframe models lie apart: -1 where the pitch, the counts or a
 * band state differ; else the largest float_apart() of their magnitudes and
 * phases. */
static double
subframe_apart(const struct mbe_subframe_model* a, const struct mbe_subframe_model* b) {
    if (!same_float(&a->w0, &b->w0) || a->L != b->L || a->any_voiced != b->any_voiced
        || a->any_pulsed != b->any_pulsed) {
        return -1.0;
    }
    double apart = 0.0;
    for (int l = 0; l < 57; ++l) {
        if (a->voiced[l] != b->voiced[l]) {
            return -1.0;
        }
        apart =
            fmax(apart, fmax(float_apart(&a->magnitude[l], &b->magnitude[l]), float_apart(&a->phase[l], &b->phase[l])));
    }
    return apart;
}

static double worst_apart;

/* Whether two subframe models are the same: pitch, counts and band states
 * exactly; magnitudes and phases by their bits, or, where the copies can
 * round apart (COPIES_INEXACT), within 1e-5 of the larger magnitude (or of 1),
 * both finite. Records the largest difference seen. */
static int
same_subframe(const struct mbe_subframe_model* a, const struct mbe_subframe_model* b) {
    const double apart = subframe_apart(a, b);
    if (apart > worst_apart) {
        worst_apart = apart;
    }
    return apart == 0.0 || (COPIES_INEXACT && apart > 0.0 && apart <= 1.0e-5);
}

/* The model P25's decoder gives pitch code b0 (w0, L and K), its bands in
 * voicing class `kind`, magnitudes and log2 amplitudes that vary by band. */
static void
code_model(mbe_parms* mp, int b0, int kind) {
    static const unsigned b0_index[8] = {0, 1, 2, 3, 4, 5, 85, 86};
    char imbe_d[88] = {0};
    mbe_parms prev, enhanced;
    for (int i = 0; i < 8; ++i) {
        imbe_d[b0_index[i]] = (char)((b0 >> (7 - i)) & 1);
    }
    mbe_initMbeParms(mp, &prev, &enhanced);
    const int rc = mbe_decodeImbe4400Parms(imbe_d, mp, &prev);
    assert(rc == 0);
    (void)rc;
    const int L = mp->L;
    for (int l = 1; l <= 56; ++l) {
        int v;
        switch (kind) {
            case 0: v = 0; break;                                 /* none */
            case 1: v = 1; break;                                 /* all */
            case 2: v = l <= L / 2; break;                        /* a low region */
            case 3: v = ((l - 1) / 3) % 2 == 0; break;            /* alternate bands of three */
            case 4: v = l == 7; break;                            /* a lone harmonic */
            case 5: v = 3 * l > 2 * L; break;                     /* the top only */
            case 6: v = (l == 5) ? 2 : (l <= (2 * L) / 3); break; /* a band with no state */
            default: v = l == L; break;                           /* the last band alone */
        }
        mp->Vl[l] = v;
        mp->Ml[l] = (float)(1.0 + 0.5 * sin(0.37 * l + b0));
        mp->log2Ml[l] = (float)(2.0 + 1.5 * cos(0.21 * l + 0.1 * b0));
    }
}

/* A copy of a model with no unvoiced band carrying a finite log2Ml: the
 * unvoiced bands' log2Ml not a number, so that P25's phase source adds no
 * step and enters the floor, as the shared one does. */
static void
without_step(mbe_parms* mp) {
    const uint32_t nan_bits = 0x7FC00000u;
    for (int l = 1; l <= mp->L; ++l) {
        if (mp->Vl[l] == 0) {
            memcpy(&mp->log2Ml[l], &nan_bits, sizeof(nan_bits));
        }
    }
}

/* The mid-frame model's pitch and count as mbe_subframe_mid_model() sets
 * them, before its bands are filled. */
static int
mid_frame(const mbe_parms* prev, const mbe_parms* cur, struct mbe_subframe_model* model) {
    memset(model, 0, sizeof(*model));
    const int shorter = (prev->L < cur->L) ? prev->L : cur->L;
    const int longer = (prev->L < cur->L) ? cur->L : prev->L;
    const float w0 = mid_pitch(prev, cur, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    if (!positive_normal(&w0)) {
        return 0;
    }
    const float harmonics = floorf((0.9254f * (float)M_PI) / w0);
    model->w0 = w0;
    model->L = (harmonics < (float)shorter) ? shorter : ((harmonics > (float)longer) ? longer : (int)harmonics);
    return 1;
}

/* The column tests and the band states of one pair. */
static long
check_columns(const mbe_parms* prev, const mbe_parms* cur, const struct mbe_subframe_model* model) {
    struct column_rule rule;
    memset(&rule, 0, sizeof(rule));
    if (!column_rule(prev, cur, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &rule)) {
        return 0;
    }
    assert(!rule.dstar);
    const int readers = p25_reads_by_columns(prev) || p25_reads_by_columns(cur);
    for (int l = 1; l <= model->L; ++l) {
        const int copy = p25_columns_voice(&rule, model->w0, l), shared = columns_voice(&rule, model->w0, l);
        assert(copy == shared);
        if (!readers) {
            const int copy_state = p25_column_band_state(prev, cur, &rule, model, l);
            const int shared_state = column_band_state(prev, cur, &rule, model, l);
            assert(copy_state == shared_state);
        }
    }
    return model->L;
}

/* The mid-frame bands of one pair, with no step: no column rule, and the rule
 * where neither frame's own model is read by it. */
static int
check_fill(const mbe_parms* prev, const mbe_parms* cur, const struct mbe_subframe_model* model, int late) {
    mbe_parms stepless_prev = *prev, stepless_cur = *cur;
    without_step(&stepless_prev);
    without_step(&stepless_cur);
    /* const objects: what the copies are handed cannot change under them */
    const mbe_parms a = stepless_prev, b = stepless_cur;
    if (a.L < 1 || a.L > 56 || b.L < 1 || b.L > 56 || model->L < 1 || model->L > 56) {
        return 0; /* every code's L is 9..56; the bound tells the analyzer */
    }
    struct mbe_subframe_model shared = *model, copy = *model;
    fill_mid_bands(&a, &b, &b, late, NULL, NULL, &shared);
    p25_fill_mid_bands(&a, &b, &b, late, NULL, NULL, &copy);
    const int alike = same_subframe(&shared, &copy);
    assert(alike);
    struct column_rule rule;
    memset(&rule, 0, sizeof(rule));
    if (column_rule(&a, &b, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &rule) && !p25_reads_by_columns(&a)
        && !p25_reads_by_columns(&b)) {
        shared = *model;
        copy = *model;
        fill_mid_bands(&a, &b, &b, late, &rule, NULL, &shared);
        p25_fill_mid_bands(&a, &b, &b, late, &rule, NULL, &copy);
        const int ruled_alike = same_subframe(&shared, &copy);
        assert(ruled_alike);
        return 1;
    }
    return 0;
}

/* The own model, with no step: P25's against the shared one with the columns
 * applied to it. */
static void
check_own(const mbe_parms* mp) {
    mbe_parms stepless = *mp;
    without_step(&stepless);
    const mbe_parms a = stepless;
    if (a.L < 1 || a.L > 56) {
        return; /* every code's L is 9..56; the bound tells the analyzer */
    }
    struct mbe_subframe_model shared, copy;
    mbe_subframe_own_model(&a, &shared);
    p25_column_voicing(&a, &shared);
    p25_own_model(&a, &copy);
    const int alike = same_subframe(&shared, &copy);
    assert(alike);
}

/* The grid test of P25's floored pulse train (p25_grid_steps()): every whole
 * number of steps n from 1 to 2^18 (the grid's Nyquist point, pi, included)
 * passes with its own canonical float, 2 pi n / 2^19 rounded to float once,
 * and gives n; both floats beside it fail; so do 0, a negative pitch and the
 * float past the last point. Returns the count of points checked. */
static long
check_grid_points(void) {
    long checked = 0;
    for (long n = 1; n <= P25_GRID_NYQUIST; ++n) {
        const float w0 = mbe_round_to_float((float)(2.0 * M_PI * (double)n / 524288.0));
        long steps = 0;
        const int on = p25_grid_steps(w0, &steps);
        assert(on && steps == n);
        uint32_t bits;
        memcpy(&bits, &w0, sizeof(bits));
        for (int side = 0; side < 2; ++side) {
            const uint32_t beside_bits = side ? bits + 1u : bits - 1u;
            float beside;
            memcpy(&beside, &beside_bits, sizeof(beside));
            long ignored = 0;
            const int beside_on = p25_grid_steps(beside, &ignored);
            assert(!beside_on);
        }
        ++checked;
    }
    const float last = mbe_round_to_float((float)M_PI);
    uint32_t past_bits;
    memcpy(&past_bits, &last, sizeof(past_bits));
    ++past_bits;
    float past, zero = 0.0f, negative = -mbe_round_to_float((float)(2.0 * M_PI / 524288.0));
    memcpy(&past, &past_bits, sizeof(past));
    long ignored = 0;
    const int past_on = p25_grid_steps(past, &ignored), zero_on = p25_grid_steps(zero, &ignored),
              negative_on = p25_grid_steps(negative, &ignored);
    assert(!past_on && !zero_on && !negative_on);
    return checked;
}

/* A voiced subframe model on the grid at n steps (harmonic 1 voiced). */
static void
grid_end(struct mbe_subframe_model* model, long n) {
    memset(model, 0, sizeof(*model));
    model->w0 = mbe_round_to_float((float)(2.0 * M_PI * (double)n / 524288.0));
    model->L = 1;
    model->voiced[1] = 1;
    model->any_voiced = 1;
}

/* The phase a floored interval leaves is the double it stores, tested as
 * stored: a constant interval at 6553 steps entered at 48 2^-19 - 2^-56
 * reaches 1 - 2^-56 cycle, which rounds to 1.0 as a double, so the phase left
 * is exactly 0 in every build (an x87 build that tested the 80-bit value kept
 * 1.0); the octave step after it (6553 to 3276, 5 n_a between 8 n_b and 12 n_b)
 * then halves 0 and launches its pulse at the interval's start, sample 80. */
static void
check_whole_cycle_then_octave(void) {
    struct mbe_subframe_model a, b;
    grid_end(&a, 6553);
    grid_end(&b, 3276);
    assert(!p25_octave_steps(6553, 6553) && p25_octave_steps(6553, 3276));
    struct pulse_train train;
    memset(&train, 0, sizeof(train));
    train.phase = (48.0 / 524288.0) - ldexp(1.0, -56);
    p25_advance_interval(&train, 0.0f, &a, &a);
    const double zero = 0.0;
    uint64_t have, want;
    memcpy(&have, &train.phase, sizeof(have));
    memcpy(&want, &zero, sizeof(want));
    assert(have == want);
    const int before = train.count;
    p25_advance_interval(&train, 80.0f, &a, &b);
    const float sample80 = 80.0f;
    assert(train.count == before + 1 && same_float(&train.time[before], &sample80));
}

int
main(void) {
    check_comparator();
    check_whole_cycle_then_octave();
    const long grid_points = check_grid_points();
    assert(grid_points == P25_GRID_NYQUIST);
    static mbe_parms models[8][208];
    for (int kind = 0; kind < 8; ++kind) {
        for (int b0 = 0; b0 < 208; ++b0) {
            code_model(&models[kind][b0], b0, kind);
            check_own(&models[kind][b0]);
        }
    }
    long bands = 0, pairs = 0, ruled = 0;
    int classes[64] = {0};
    for (int a = 0; a < 208; ++a) {
        for (int b = 0; b < 208; ++b) {
            const int c = ((a * 7) + b) % 64;
            const mbe_parms* prev = &models[c / 8][a];
            mbe_parms cur = models[c % 8][b];
            if ((c % 8) == 0 && ((a + b) % 2) == 0) {
                /* A frame with no voiced band lacking a pitch: the column rule
                 * applies, and with no voiced band before it neither own
                 * model is read by it. */
                mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PITCH_ABSENT);
            }
            struct mbe_subframe_model model;
            if (!mid_frame(prev, &cur, &model)) {
                continue;
            }
            bands += check_columns(prev, &cur, &model);
            ruled += check_fill(prev, &cur, &model, (a + b) % 3 == 0);
            ++pairs;
            classes[c] = 1;
        }
    }
    for (int c = 0; c < 64; ++c) {
        assert(classes[c]);
    }
    assert(pairs == 208L * 208L && bands > 0 && ruled > 0);
    printf("p25 rendering copies: %ld code pairs, %ld mid-frame bands by the columns, %ld fills under the rule, "
           "alike (%s; the largest float difference %.3g of its value); %ld grid points and their neighbours\n",
           pairs, bands, ruled, COPIES_INEXACT ? "within 1e-5" : "bit for bit", worst_apart, grid_points);
    return 0;
}
