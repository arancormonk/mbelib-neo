// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Numerics of the MBE speech analysis used by the encoders.
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "mbe_speech_analysis.h"
#include "mbe_unvoiced_fft.h"
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
test_harmonic_fit(mbe_fft_plan* fft, mbe_acf_plan* acf) {
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
            synthesize(&s, in);
            mbe_analysis_push(&state, in);
            assert(mbe_analysis_frame(&tables, &state, fft, acf, &r) == 0);
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

/* White noise: E(P) near 1 on average and no alternating-width magnitude
 * ripple. The global pitch gate leaves unused fit errors at 1; this is an
 * unvoiced default, not a measured noise-disagreement value. */
static void
test_noise(mbe_fft_plan* fft, mbe_acf_plan* acf) {
    struct mbe_analysis_state state;
    double error_sum = 0.0;
    double fit_sum = 0.0;
    double ripple[MBE_ANALYSIS_HARMONICS + 1] = {0};
    int fit_count = 0;
    int frames = 0;
    int gated_frames = 0;
    mbe_analysis_reset(&state);
    for (int f = 0; f < 300; f++) {
        float in[MBE_ANALYSIS_FRAME];
        struct mbe_analysis_result r;
        for (int i = 0; i < MBE_ANALYSIS_FRAME; i++) {
            in[i] = (float)(2000.0 * gauss());
        }
        mbe_analysis_push(&state, in);
        assert(mbe_analysis_frame(&tables, &state, fft, acf, &r) == 0);
        if (f < 4) {
            continue;
        }
        frames++;
        error_sum += r.pitch_error;
        if (r.pitch_error > 0.4f) { /* the existing global voicing gate */
            gated_frames++;
            for (int l = 1; l <= r.harmonics; l++) {
                assert(r.fit_error[l] == 1.0f);
            }
            /* Even a previously voiced band cannot pass the closed gate. */
            struct mbe_analysis_state voiced_state = state;
            memset(voiced_state.columns_prev, 1, sizeof(voiced_state.columns_prev));
            for (int k = 1; k <= MBE_ANALYSIS_COLUMNS; k++) {
                assert(mbe_analysis_column_threshold(&voiced_state, &r, k) == 0.0f);
            }
            unsigned char prev[MBE_ANALYSIS_BANDS];
            unsigned char bands[MBE_ANALYSIS_BANDS];
            memset(prev, 1, sizeof(prev));
            mbe_analysis_band_voicing(&r, r.harmonics, prev, bands);
            for (int k = 0; k < MBE_ANALYSIS_BANDS; k++) {
                assert(bands[k] == 0);
            }
        }
        const double f0_bins = r.f0 * 256.0;
        for (int l = 2; l < r.harmonics; l++) {
            fit_sum += r.fit_error[l];
            fit_count++;
            /* M_l^2 = sigma^2 f0_bins / 256 for flat noise (DC-filter gain ~1) */
            ripple[l] += (r.magnitude[l] * r.magnitude[l]) / (2000.0 * 2000.0 * f0_bins / 256.0);
        }
    }
    assert(gated_frames > 0);
    assert(error_sum / frames > 0.6);
    /* Includes the unvoiced default on gated frames. Measured periodic fits
     * remain covered by test_harmonic_fit() and the hum/voice regressions. */
    assert(fit_sum / fit_count > 0.2);
    double lo = 1e9;
    double hi = 0.0;
    for (int l = 4; l <= 12; l++) {
        double v = ripple[l] / frames;
        lo = (v < lo) ? v : lo;
        hi = (v > hi) ? v : hi;
    }
    assert(10.0 * log10(hi / lo) < 1.0);
}

/* Feed one frame and analyse it. */
static void
analyse(mbe_fft_plan* fft, mbe_acf_plan* acf, struct mbe_analysis_state* state, const float in[MBE_ANALYSIS_FRAME],
        struct mbe_analysis_result* r) {
    mbe_analysis_push(state, in);
    assert(mbe_analysis_frame(&tables, state, fft, acf, r) == 0);
}

/* Count bands wholly above min_hz, and their voiced decisions. */
static int
voiced_bands(const struct mbe_analysis_result* r, double min_hz, int* bands_out) {
    unsigned char prev[MBE_ANALYSIS_BANDS] = {0};
    unsigned char bands[MBE_ANALYSIS_BANDS];
    const int L = (int)(0.9254 * 0.5 / r->f0);
    const int K = mbe_analysis_band_voicing(r, L, prev, bands);
    int voiced = 0;
    *bands_out = 0;
    for (int k = 0; k < K; k++) {
        if ((3 * k + 1) * r->f0 * 8000.0 >= min_hz) {
            voiced += bands[k];
            (*bands_out)++;
        }
    }
    return voiced;
}

/* Mains hum at the level of #149's quiet stretches (about 75 LSB, with weak
 * 2nd and 3rd harmonics) over its noise floor (13 LSB RMS), and over a nearly
 * silent one: at most one upper band in ten is voiced, as in the encoder
 * noise regressions. The hum itself is periodic; it must not turn unrelated
 * noise above 1 kHz into harmonics. */
static void
test_hum_unvoiced(mbe_fft_plan* fft, mbe_acf_plan* acf) {
    const double mains[] = {60.0, 50.0, 60.0};
    const double floor_rms[] = {13.0, 13.0, 1.0};
    for (size_t c = 0; c < sizeof(mains) / sizeof(mains[0]); c++) {
        struct mbe_analysis_state state;
        double phase = 0.0;
        int voiced = 0;
        int bands = 0;
        mbe_analysis_reset(&state);
        for (int f = 0; f < 150; f++) {
            float in[MBE_ANALYSIS_FRAME];
            struct mbe_analysis_result r;
            for (int i = 0; i < MBE_ANALYSIS_FRAME; i++) {
                phase += 2.0 * M_PI * mains[c] / 8000.0;
                in[i] = (float)((75.0 * sin(phase)) + (7.0 * sin((2.0 * phase) + 1.0))
                                + (7.0 * sin((3.0 * phase) + 2.0)) + (floor_rms[c] * gauss()));
            }
            analyse(fft, acf, &state, in, &r);
            if (f >= 50) {
                int count;
                voiced += voiced_bands(&r, 1000.0, &count);
                bands += count;
            }
        }
        printf("hum %.0f Hz, noise %.0f LSB: %d/%d upper bands voiced\n", mains[c], floor_rms[c], voiced, bands);
        assert(voiced * 10 <= bands);
    }
}

/* Steady voices at the bottom of the pitch range, the lowest the codecs send
 * (62.5 Hz) included, keep their pitch, pass the periodicity gate and stay voiced. */
static void
test_low_voices(mbe_fft_plan* fft, mbe_acf_plan* acf) {
    const double f0s[] = {62.5, 66.0, 75.0};
    for (size_t c = 0; c < sizeof(f0s) / sizeof(f0s[0]); c++) {
        struct mbe_analysis_state state;
        struct mbe_analysis_result r;
        struct harmonic_signal s = {.f0 = f0s[c], .amp = 4000.0, .count = (int)(3500.0 / f0s[c])};
        for (int h = 1; h <= s.count; h++) {
            s.phase[h] = 0.7 * h * h;
        }
        mbe_analysis_reset(&state);
        for (int f = 0; f < 100; f++) {
            float in[MBE_ANALYSIS_FRAME];
            synthesize(&s, in);
            analyse(fft, acf, &state, in, &r);
        }
        assert(near(1.0 / r.f0, 8000.0 / s.f0, 0.5));
        assert(r.pitch_error < 0.25f);
        int bands;
        const int voiced = voiced_bands(&r, 0.0, &bands);
        assert(voiced * 10 >= bands * 9);
    }
}

/* White noise at the given RMS for `frames` frames. */
static void
feed_noise(mbe_fft_plan* fft, mbe_acf_plan* acf, struct mbe_analysis_state* state, double rms, int frames,
           struct mbe_analysis_result* r) {
    for (int f = 0; f < frames; f++) {
        float in[MBE_ANALYSIS_FRAME];
        for (int i = 0; i < MBE_ANALYSIS_FRAME; i++) {
            in[i] = (float)(rms * gauss());
        }
        analyse(fft, acf, state, in, r);
    }
}

/* A sustained low voice with weak upper harmonics (1/h^2), in noise 10 dB
 * down, retains voicing both on its own and after a vowel 13 dB louder. */
static void
test_loud_low_voice(mbe_fft_plan* fft, mbe_acf_plan* acf) {
    for (int history = 0; history <= 1; history++) {
        struct mbe_analysis_state state;
        struct mbe_analysis_result r;
        double phase[64] = {0};
        int voiced = 0;
        int frames = 0;
        mbe_analysis_reset(&state);
        if (history) {
            struct harmonic_signal loud = {.f0 = 150.0, .amp = 16000.0, .count = 23};
            for (int f = 0; f < 100; f++) {
                float in[MBE_ANALYSIS_FRAME];
                synthesize(&loud, in);
                analyse(fft, acf, &state, in, &r);
            }
        }
        for (int f = 0; f < 200; f++) {
            float in[MBE_ANALYSIS_FRAME];
            for (int i = 0; i < MBE_ANALYSIS_FRAME; i++) {
                double v = 0.0;
                for (int h = 1; h <= 53; h++) {
                    phase[h] += 2.0 * M_PI * h * 66.0 / 8000.0;
                    v += (4000.0 / (h * h)) * sin(phase[h] + (0.7 * h * h));
                }
                in[i] = (float)(v + (900.0 * gauss()));
            }
            analyse(fft, acf, &state, in, &r);
            if (f >= 50) {
                int bands;
                voiced += voiced_bands(&r, 0.0, &bands) > 0;
                frames++;
            }
        }
        assert(voiced * 20 >= frames * 19);
    }
}

/* Noise keeps its spectral magnitude; a voice starting abruptly out of it
 * is voiced from the frame that fills the analysis window. */
static void
test_voice_onset(mbe_fft_plan* fft, mbe_acf_plan* acf) {
    struct mbe_analysis_state state;
    struct mbe_analysis_result r;
    mbe_analysis_reset(&state);
    feed_noise(fft, acf, &state, 300.0, 100, &r);
    double magnitude = 0.0;
    int count = 0;
    for (int f = 0; f < 50; f++) {
        feed_noise(fft, acf, &state, 300.0, 1, &r);
        for (int l = 2; l < r.harmonics; l++) {
            /* M_l^2 = sigma^2 f0_bins / 256 for flat noise, as in test_noise() */
            magnitude += r.magnitude[l] * r.magnitude[l] / (300.0 * 300.0 * r.f0);
            count++;
        }
    }
    assert(fabs(10.0 * log10(magnitude / count)) < 1.5); /* magnitudes untouched */

    struct harmonic_signal s = {.f0 = 150.0, .amp = 4000.0, .count = 23};
    for (int f = 0; f < 10; f++) {
        float in[MBE_ANALYSIS_FRAME];
        synthesize(&s, in);
        for (int i = 0; i < MBE_ANALYSIS_FRAME; i++) {
            in[i] += (float)(300.0 * gauss());
        }
        analyse(fft, acf, &state, in, &r);
        if (f >= 1) { /* once the onset fills the analysis window */
            int bands;
            const int voiced = voiced_bands(&r, 0.0, &bands);
            assert(voiced * 2 >= bands);
            for (int h = 1; h <= 4; h++) {
                assert(fabs(20.0 * log10(r.magnitude[h] / ((s.amp / h) / 2.0))) < 1.5);
            }
        }
    }
}

/* A quieter sustained vowel following loud noise, directly or after silence,
 * promptly recovers its voicing and magnitude. */
static void
test_noise_recovery(mbe_fft_plan* fft, mbe_acf_plan* acf) {
    for (int silence = 0; silence <= 1; silence++) {
        struct mbe_analysis_state state;
        struct mbe_analysis_result r;
        mbe_analysis_reset(&state);
        feed_noise(fft, acf, &state, 700.0, 100, &r); /* -33 dBFS */
        if (silence) {
            feed_noise(fft, acf, &state, 0.0, 10, &r);
        }
        struct harmonic_signal s = {.f0 = 150.0, .amp = 400.0, .count = 23}; /* under the noise */
        for (int f = 0; f < 60; f++) {
            float in[MBE_ANALYSIS_FRAME];
            synthesize(&s, in);
            analyse(fft, acf, &state, in, &r);
            if (f >= 3) {
                int bands;
                const int voiced = voiced_bands(&r, 0.0, &bands);
                assert(voiced * 2 >= bands);
            }
        }
        int bands;
        const int voiced = voiced_bands(&r, 0.0, &bands);
        assert(voiced * 2 >= bands); /* most bands voiced */
        for (int h = 1; h <= 8; h++) {
            assert(fabs(20.0 * log10(r.magnitude[h] / ((s.amp / h) / 2.0))) < 1.5);
        }
    }
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
    mbe_analysis_commit(&state, voiced);
    assert(memcmp(state.columns_prev, voiced, sizeof(voiced)) == 0);
    mbe_analysis_commit(&state, NULL);
    assert(memcmp(state.columns_prev, voiced, sizeof(voiced)) == 0);
    assert(mbe_analysis_frame(NULL, &state, NULL, NULL, NULL) == MBE_STATUS_INVALID_ARGUMENT);
}

/* Check the 256-point ordered spectrum against a direct DFT, including DC
 * and Nyquist packing. This exercises the production FFT integration only. */
static void
test_fft_spectrum(mbe_fft_plan* fft) {
    float input[MBE_FFT_SIZE];
    float spectrum[MBE_FFT_SIZE];
    float again[MBE_FFT_SIZE];
    for (int n = 0; n < MBE_FFT_SIZE; ++n) {
        input[n] = (float)(0.3 + 0.7 * cos(2.0 * M_PI * 13.0 * n / MBE_FFT_SIZE)
                           + 0.2 * sin(2.0 * M_PI * 31.0 * n / MBE_FFT_SIZE));
        input[n] += n % 2 ? -0.1f : 0.1f;
    }
    assert(mbe_fft_forward_real(fft, input, spectrum) == 0);
    assert(mbe_fft_forward_real(fft, input, again) == 0);
    for (int n = 0; n < MBE_FFT_SIZE; ++n) {
        uint32_t first_bits;
        uint32_t second_bits;
        memcpy(&first_bits, &spectrum[n], sizeof(first_bits));
        memcpy(&second_bits, &again[n], sizeof(second_bits));
        assert(first_bits == second_bits);
    }
    for (int k = 0; k <= MBE_FFT_SIZE / 2; ++k) {
        double real = 0.0;
        double imag = 0.0;
        for (int n = 0; n < MBE_FFT_SIZE; ++n) {
            double phase = 2.0 * M_PI * k * n / MBE_FFT_SIZE;
            real += input[n] * cos(phase);
            imag -= input[n] * sin(phase);
        }
        int index = k == MBE_FFT_SIZE / 2 ? 1 : 2 * k;
        assert(fabs(spectrum[index] - real) < 1e-4);
        if (k > 0 && k < MBE_FFT_SIZE / 2) {
            assert(fabs(spectrum[index + 1] - imag) < 1e-4);
        }
    }
}

/* FFT autocorrelation matches the direct sum for every lag. */
static void
test_autocorrelation(mbe_acf_plan* acf) {
    float x[301];
    double scale = 0.0;
    for (int n = 0; n < 301; ++n) {
        x[n] = (float)(gauss() * 1000.0);
    }
    float r[153];
    assert(mbe_acf_compute(acf, x, 301, r, 152) == 0);
    for (int n = 0; n < 301; ++n) {
        scale += (double)x[n] * x[n];
    }
    for (int lag = 0; lag <= 152; ++lag) {
        double want = 0.0;
        for (int n = 0; n + lag < 301; ++n) {
            want += (double)x[n] * x[n + lag];
        }
        assert(fabs((double)r[lag] - want) <= 1e-5 * scale);
    }
    assert(mbe_acf_compute(acf, x, 400, r, 152) == MBE_STATUS_INVALID_ARGUMENT); /* 400 + 152 > 512 */
    assert(mbe_acf_compute(NULL, x, 301, r, 152) == MBE_STATUS_INVALID_ARGUMENT);
}

/* Counts at and past the transform size: lags up to 511, len + lags up to 512,
 * and sums that would overflow int are rejected before any access. */
static void
test_autocorrelation_bounds(mbe_acf_plan* acf) {
    static float x[MBE_ACF_SIZE];
    static float r[MBE_ACF_SIZE];
    for (int n = 0; n < MBE_ACF_SIZE; ++n) {
        x[n] = 1.0f;
        r[n] = -1.0f;
    }
    assert(mbe_acf_compute(acf, x, 0, r, MBE_ACF_SIZE - 1) == 0);
    for (int lag = 0; lag < MBE_ACF_SIZE; ++lag) {
        assert(fabsf(r[lag]) < 1e-6f); /* empty input */
    }
    assert(mbe_acf_compute(acf, x, MBE_ACF_SIZE, r, 0) == 0 && fabsf(r[0] - (float)MBE_ACF_SIZE) < 1e-2f);
    assert(mbe_acf_compute(acf, x, 301, r, MBE_ACF_SIZE - 301) == 0);
    assert(fabsf(r[MBE_ACF_SIZE - 301] - 90.0f) < 1e-3f); /* 301 - 211 overlapping ones */
    assert(mbe_acf_compute(acf, x, 0, r, MBE_ACF_SIZE) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(acf, x, 301, r, MBE_ACF_SIZE - 300) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(acf, x, MBE_ACF_SIZE + 1, r, 0) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(acf, x, INT_MAX, r, INT_MAX) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(acf, x, -1, r, 10) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(acf, x, 10, r, -1) == MBE_STATUS_INVALID_ARGUMENT);
}

int
main(void) {
    mbe_fft_plan* fft = mbe_fft_plan_alloc();
    assert(fft != NULL);
    test_fft_spectrum(fft);
    mbe_acf_plan* acf = mbe_acf_plan_alloc();
    assert(acf != NULL);
    test_autocorrelation(acf);
    test_autocorrelation_bounds(acf);
    mbe_analysis_init_tables(&tables);
    test_windows();
    test_window_transform();
    test_harmonic_fit(fft, acf);
    test_noise(fft, acf);
    test_hum_unvoiced(fft, acf);
    test_low_voices(fft, acf);
    test_voice_onset(fft, acf);
    test_loud_low_voice(fft, acf);
    test_noise_recovery(fft, acf);
    test_thresholds();
    test_reset_and_commit();
    mbe_acf_plan_free(acf);
    mbe_fft_plan_free(fft);
    puts("speech analysis: ok");
    return 0;
}
