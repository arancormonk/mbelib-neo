// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Numerics of the MBE speech analysis used by the encoders.
 */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbe_speech_analysis.h"
#include "mbelib-neo/mbelib.h"

static struct mbe_analysis_tables tables;
static uint32_t rng = 0x5EEDu;

static double
gauss(void) {
    double sum = 0.0;
    for (int i = 0; i < 12; i++) {
        rng = (rng * 1664525u) + 1013904223u;
        sum += (double)(rng >> 8) / 16777216.0;
    }
    return sum - 6.0;
}

static double
bessel_i0(double x) {
    double sum = 1.0;
    double term = 1.0;
    for (int k = 1; k < 60; k++) {
        term *= (x / (2.0 * k)) * (x / (2.0 * k));
        sum += term;
    }
    return sum;
}

static int
near(double a, double b, double tolerance) {
    return fabs(a - b) <= tolerance;
}

/* The windows are the Kaiser windows the standard's Annex B and C tabulate. */
static void
test_windows(void) {
    double sum_sq = 0.0;
    for (int n = 0; n <= 2 * MBE_ANALYSIS_PITCH_HALF; n++) {
        sum_sq += (double)tables.w_i[n] * tables.w_i[n];
    }
    assert(near(sum_sq, 1.0, 1e-5)); /* eq (6) */
    double peak = tables.w_i[MBE_ANALYSIS_PITCH_HALF];
    assert(near(tables.w_i[0] / peak, 1.0 / bessel_i0(5.25), 1e-6));
    assert(near(tables.w_r[MBE_ANALYSIS_WR_HALF], 1.0, 1e-7));
    assert(near(tables.w_r[0], 1.0 / bessel_i0(6.0), 1e-7));
    assert(near(tables.w_r[(size_t)2 * MBE_ANALYSIS_WR_HALF], tables.w_r[0], 1e-9));

    double lpf_sum = 0.0;
    for (int n = 0; n <= 2 * MBE_ANALYSIS_LPF_HALF; n++) {
        lpf_sum += tables.lpf[n];
        assert(near(tables.lpf[n], tables.lpf[(2 * MBE_ANALYSIS_LPF_HALF) - n], 1e-9));
    }
    assert(near(lpf_sum, 1.0, 1e-6));
}

/* W_R(nu) matches the signed-index DTFT of w_R. */
static void
test_window_transform(void) {
    const double offsets[] = {0.0, 0.5, 1.37, 2.0, 3.2, -4.75, 6.5};
    const int centre = MBE_ANALYSIS_WTAB_BINS * MBE_ANALYSIS_WTAB_STEPS;
    for (size_t i = 0; i < sizeof(offsets) / sizeof(offsets[0]); i++) {
        int index = centre + (int)lround(offsets[i] * MBE_ANALYSIS_WTAB_STEPS);
        double nu = (double)(index - centre) / MBE_ANALYSIS_WTAB_STEPS; /* the table's grid point */
        double direct = 0.0;
        for (int n = -MBE_ANALYSIS_WR_HALF; n <= MBE_ANALYSIS_WR_HALF; n++) {
            direct += tables.w_r[n + MBE_ANALYSIS_WR_HALF] * cos(2.0 * M_PI * nu * n / 256.0);
        }
        assert(near(tables.w_r_dtft[index], direct, 1e-4 * tables.w_r_sum));
    }
    assert(near(tables.w_r_dtft[centre], tables.w_r_sum, 1e-3));
}

struct harmonic_signal {
    double f0;  /* Hz */
    double amp; /* 16-bit scale */
    int count;
    double phase[64];
};

static void
synthesize(struct harmonic_signal* s, float out[MBE_ANALYSIS_FRAME]) {
    for (int i = 0; i < MBE_ANALYSIS_FRAME; i++) {
        double v = 0.0;
        for (int h = 1; h <= s->count; h++) {
            s->phase[h] += 2.0 * M_PI * h * s->f0 / 8000.0;
            v += (s->amp / h) * sin(s->phase[h]);
        }
        out[i] = (float)v;
    }
}

/* Off-grid harmonics with arbitrary phases: pitch within a quarter sample,
 * magnitudes A/2 within 1 dB for resolved harmonics, small fit errors. */
static void
test_harmonic_fit(mbe_fft_plan* fft) {
    const double f0s[] = {137.3, 181.9, 233.7};
    for (size_t c = 0; c < sizeof(f0s) / sizeof(f0s[0]); c++) {
        struct mbe_analysis_state state;
        struct mbe_analysis_result r;
        struct harmonic_signal s = {.f0 = f0s[c], .amp = 4000.0, .count = (int)(3500.0 / f0s[c])};
        for (int h = 1; h <= s.count; h++) {
            s.phase[h] = 0.7 * h * h;
        }
        mbe_analysis_reset(&state);
        for (int f = 0; f < 8; f++) {
            float in[MBE_ANALYSIS_FRAME];
            float filtered[MBE_ANALYSIS_FRAME];
            synthesize(&s, in);
            mbe_analysis_push(&state, in, filtered);
            assert(mbe_analysis_frame(&tables, &state, fft, 1.0f, &r) == 0);
        }
        double period = 8000.0 / s.f0;
        assert(near(1.0 / r.f0, period, 0.3));
        for (int h = 1; h <= s.count && h <= 8; h++) {
            double db = 20.0 * log10(r.magnitude[h] / ((s.amp / h) / 2.0));
            assert(fabs(db) < 1.0);
            assert(r.fit_error[h] < 0.05f);
        }
        assert(r.pitch_error < 0.25f);
    }
}

/* White noise: E(P) near 1 on average, fit errors large, and no
 * alternating-width ripple across harmonics on average. */
static void
test_noise(mbe_fft_plan* fft) {
    struct mbe_analysis_state state;
    double error_sum = 0.0;
    double fit_sum = 0.0;
    double ripple[MBE_ANALYSIS_HARMONICS + 1] = {0};
    int fit_count = 0;
    int frames = 0;
    mbe_analysis_reset(&state);
    for (int f = 0; f < 300; f++) {
        float in[MBE_ANALYSIS_FRAME];
        float filtered[MBE_ANALYSIS_FRAME];
        struct mbe_analysis_result r;
        for (int i = 0; i < MBE_ANALYSIS_FRAME; i++) {
            in[i] = (float)(2000.0 * gauss());
        }
        mbe_analysis_push(&state, in, filtered);
        assert(mbe_analysis_frame(&tables, &state, fft, 1.0f, &r) == 0);
        if (f < 4) {
            continue;
        }
        frames++;
        error_sum += r.pitch_error;
        const double f0_bins = r.f0 * 256.0;
        for (int l = 2; l < r.harmonics; l++) {
            fit_sum += r.fit_error[l];
            fit_count++;
            /* M_l^2 = sigma^2 f0_bins / 256 for flat noise (DC-filter gain ~1) */
            ripple[l] += (r.magnitude[l] * r.magnitude[l]) / (2000.0 * 2000.0 * f0_bins / 256.0);
        }
    }
    assert(error_sum / frames > 0.6);
    assert(fit_sum / fit_count > 0.2); /* periodic signals stay below 0.05 */
    double lo = 1e9;
    double hi = 0.0;
    for (int l = 4; l <= 12; l++) {
        double v = ripple[l] / frames;
        lo = (v < lo) ? v : lo;
        hi = (v > hi) ? v : hi;
    }
    assert(10.0 * log10(hi / lo) < 1.0);
}

/* eq (37) at omega0 0.1309 with the hysteresis, the frequency taper and the
 * periodicity gate. */
static void
test_thresholds(void) {
    struct mbe_analysis_state state;
    struct mbe_analysis_result r;
    mbe_analysis_reset(&state);
    memset(&r, 0, sizeof(r));
    r.pitch_error = 0.2f;
    r.energy_factor = 1.0f;
    assert(near(mbe_analysis_column_threshold(&state, &r, 1), 0.45, 1e-6));
    state.columns_prev[0] = 1;
    assert(near(mbe_analysis_column_threshold(&state, &r, 1), 0.5625, 1e-6));
    double taper = 1.0 - (0.3096 * 7.0 * 0.1309);
    assert(near(mbe_analysis_column_threshold(&state, &r, 8), 0.45 * taper, 1e-6));
    r.energy_factor = 0.5f;
    assert(near(mbe_analysis_column_threshold(&state, &r, 1), 0.28125, 1e-6));
    r.pitch_error = 0.45f;
    for (int k = 1; k <= MBE_ANALYSIS_COLUMNS; k++) {
        assert(mbe_analysis_column_threshold(&state, &r, k) == 0.0f);
    }
}

static void
test_reset_and_commit(void) {
    struct mbe_analysis_state state;
    const unsigned char voiced[MBE_ANALYSIS_COLUMNS] = {1, 1, 0, 0, 1, 1, 0, 0};
    mbe_analysis_reset(&state);
    assert(near(state.pitch_prev[0], 100.0, 0.0) && near(state.pitch_prev[1], 100.0, 0.0));
    assert(state.trusted == 0 && near(state.xi_max, 20000.0, 0.0));
    mbe_analysis_commit(&state, voiced, 1);
    assert(memcmp(state.columns_prev, voiced, sizeof(voiced)) == 0);
    state.trusted = 1;
    state.pitch_prev[0] = 40.0f;
    mbe_analysis_commit(&state, NULL, 0);
    assert(state.trusted == 0 && near(state.pitch_prev[0], 100.0, 0.0));
    for (int k = 0; k < MBE_ANALYSIS_COLUMNS; k++) {
        assert(state.columns_prev[k] == 0);
    }
    assert(mbe_analysis_frame(NULL, &state, NULL, 1.0f, NULL) == MBE_STATUS_INVALID_ARGUMENT);
}

int
main(void) {
    mbe_fft_plan* fft = mbe_fft_plan_alloc();
    assert(fft != NULL);
    mbe_analysis_init_tables(&tables);
    test_windows();
    test_window_transform();
    test_harmonic_fit(fft);
    test_noise(fft);
    test_thresholds();
    test_reset_and_commit();
    mbe_fft_plan_free(fft);
    puts("speech analysis: ok");
    return 0;
}
