// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief MBE speech analysis: pitch, per-harmonic voicing measures and
 *        spectral magnitudes for the encoders.
 *
 * Implements the analysis method of ANSI/TIA-102.BABA chapter 5 from its
 * equations (numbers below refer to that document):
 *
 *  - DC removal, eq (3).
 *  - Initial pitch from the error function E(P), eqs (5)-(9), with look-back
 *    tracking, eqs (10)-(12), sub-multiple checks, eqs (18)-(20), and the
 *    decision rules, eqs (21)-(23).
 *  - Pitch refinement to a quarter sample, eqs (24)-(28).
 *  - Per-harmonic fit error and the energy factor M(xi), eqs (35)-(42), and the
 *    three-harmonic band V/UV decisions of 5.2.
 *  - The voiced and unvoiced spectral amplitudes of eqs (43) and (44) next to
 *    the voicing-independent one.
 *
 * The analysis windows are generated rather than copied: the standard's
 * Annex B and C tables are Kaiser windows (beta 5.25 over +-150 samples and
 * beta 6 over +-110 samples; they match to the printed precision). The pitch
 * low-pass is this file's own 21-tap Kaiser-windowed sinc.
 *
 * Deviations, all documented where they occur:
 *  - Look-ahead tracking needs two future frames (40 ms more delay). It is
 *    replaced by CE_F(P) = 3 E(P), which keeps the 3-frame scale of CE_B so the
 *    standard's thresholds apply. Look-back is only used when the previous
 *    frame was periodic; otherwise acquisition always runs the sub-multiple
 *    checks.
 *  - The pitch range is extended from 21-122 to 20-128 to cover D-STAR's.
 *  - Magnitudes use the voicing-independent estimate of US 5,701,390
 *    (eqs 3-4, expired), which the AMBE decoders' unvoiced scaling expects.
 *  - The periodicity gate of eq (37) is 0.4 and covers every column and band;
 *    see mbe_analysis_column_threshold() and mbe_analysis_band_voicing().
 *  - Background noise is tracked, and at low pitches a 500 Hz column that
 *    holds little more than it is not voiced; see "Background noise" below.
 *    The standard's analysis has no such step.
 */
#include "mbe_speech_analysis.h"

#include <math.h>
#include <string.h>

#include "mbe_simd.h" // IWYU pragma: keep (MBE_VF_WIDTH selects the SIMD pitch filter)
#include "mbelib-neo/mbelib.h"

#define ANALYSIS_FFT       256
#define ANALYSIS_BINS      128 /* ANALYSIS_FFT / 2 */
#define PITCH_GRID         2   /* candidates per sample of pitch */
#define PITCH_MIN_GRID     40  /* PITCH_MIN in grid steps */
#define PITCH_MIN          ((float)PITCH_MIN_GRID / (float)PITCH_GRID)
#define PITCH_STEP         (1.0f / (float)PITCH_GRID)
#define LPF_CUTOFF         0.1757 /* cycles per sample (about 1.4 kHz) */
#define LPF_BETA           4.18
#define W_I_BETA           5.25
#define W_R_BETA           6.0
#define REFINE_FIRST_BIN   50 /* eq (24) */
#define REFINE_CANDIDATES  10
#define XI_MAX_FLOOR       20000.0f
#define PERIODIC_ERROR_MAX 0.5f /* periodic enough to trust the pitch track */
#define VOICING_ERROR_MAX  0.4f /* see mbe_analysis_column_threshold() */
#define SAMPLE_CENTRE      (MBE_ANALYSIS_HISTORY - 1)

/* Modified Bessel function of the first kind, order 0 (power series). */
static double
bessel_i0(double x) {
    double sum = 1.0;
    double term = 1.0;
    for (int k = 1; k < 50; ++k) {
        term *= (x / (2.0 * (double)k)) * (x / (2.0 * (double)k));
        sum += term;
        if (term < sum * 1e-17) {
            break;
        }
    }
    return sum;
}

static double
kaiser(int n, int half, double beta) {
    double r = (double)n / (double)half;
    return bessel_i0(beta * sqrt(fmax(0.0, 1.0 - (r * r)))) / bessel_i0(beta);
}

static void
init_windows(struct mbe_analysis_tables* t) {
    double sum_sq = 0.0;
    for (int n = -MBE_ANALYSIS_PITCH_HALF; n <= MBE_ANALYSIS_PITCH_HALF; ++n) {
        double w = kaiser(n, MBE_ANALYSIS_PITCH_HALF, W_I_BETA);
        sum_sq += w * w;
    }
    double w_i4 = 0.0;
    for (int n = -MBE_ANALYSIS_PITCH_HALF; n <= MBE_ANALYSIS_PITCH_HALF; ++n) {
        double w = kaiser(n, MBE_ANALYSIS_PITCH_HALF, W_I_BETA) / sqrt(sum_sq); /* eq (6): sum w_I^2 = 1 */
        t->w_i[n + MBE_ANALYSIS_PITCH_HALF] = (float)w;
        w_i4 += w * w * w * w;
    }
    t->w_i4_sum = (float)w_i4;

    double w_r_sum = 0.0;
    double w_r_sq = 0.0;
    for (int n = -MBE_ANALYSIS_WR_HALF; n <= MBE_ANALYSIS_WR_HALF; ++n) {
        double w = kaiser(n, MBE_ANALYSIS_WR_HALF, W_R_BETA);
        t->w_r[n + MBE_ANALYSIS_WR_HALF] = (float)w;
        w_r_sum += w;
        w_r_sq += w * w;
    }
    t->w_r_sum = (float)w_r_sum;
    t->w_r_sq_sum = (float)w_r_sq;
}

/* The covariance of white noise's windowed spectrum between bins d apart. */
static void
init_noise_cov(struct mbe_analysis_tables* t) {
    for (int d = 0; d < 16; ++d) {
        double c = 0.0;
        for (int n = -MBE_ANALYSIS_WR_HALF; n <= MBE_ANALYSIS_WR_HALF; ++n) {
            double w = t->w_r[n + MBE_ANALYSIS_WR_HALF];
            c += w * w * cos(2.0 * M_PI * (double)d * (double)n / ANALYSIS_FFT);
        }
        t->noise_cov[d] = (float)c;
    }
}

static void
init_lowpass(struct mbe_analysis_tables* t) {
    double taps[(2 * MBE_ANALYSIS_LPF_HALF) + 1];
    double sum = 0.0;
    for (int n = -MBE_ANALYSIS_LPF_HALF; n <= MBE_ANALYSIS_LPF_HALF; ++n) {
        double x = 2.0 * M_PI * LPF_CUTOFF * (double)n;
        double sinc = (n == 0) ? 1.0 : sin(x) / x;
        taps[n + MBE_ANALYSIS_LPF_HALF] = 2.0 * LPF_CUTOFF * sinc * kaiser(n, MBE_ANALYSIS_LPF_HALF, LPF_BETA);
        sum += taps[n + MBE_ANALYSIS_LPF_HALF];
    }
    for (int i = 0; i < (2 * MBE_ANALYSIS_LPF_HALF) + 1; ++i) {
        t->lpf[i] = (float)(taps[i] / sum); /* unit gain at DC */
    }
}

/* W_R(nu) of eq (30) for nu in bins of the 256-point DFT: real, since w_R is
 * symmetric about n = 0. Zero beyond the table, where bands never reach. */
static void
init_window_transform(struct mbe_analysis_tables* t) {
    const int centre = MBE_ANALYSIS_WTAB_BINS * MBE_ANALYSIS_WTAB_STEPS;
    /* Both w_R and its transform are even: fold the sum and mirror the table. */
    for (int i = 0; i <= centre; ++i) {
        double nu = (double)i / (double)MBE_ANALYSIS_WTAB_STEPS;
        double sum = (double)t->w_r[MBE_ANALYSIS_WR_HALF];
        for (int n = 1; n <= MBE_ANALYSIS_WR_HALF; ++n) {
            sum += 2.0 * (double)t->w_r[n + MBE_ANALYSIS_WR_HALF] * cos(2.0 * M_PI * nu * (double)n / ANALYSIS_FFT);
        }
        t->w_r_dtft[centre + i] = (float)sum;
        t->w_r_dtft[centre - i] = (float)sum;
    }
}

void
mbe_analysis_init_tables(struct mbe_analysis_tables* tables) {
    if (tables == NULL) {
        return;
    }
    init_windows(tables);
    init_lowpass(tables);
    init_window_transform(tables);
    init_noise_cov(tables);
}

void
mbe_analysis_reset(struct mbe_analysis_state* state) {
    if (state == NULL) {
        return;
    }
    memset(state, 0, sizeof(*state));
    for (int m = 0; m < MBE_ANALYSIS_BINS; ++m) {
        state->noise_cur[m] = -1.0f; /* no noise-like frame taken */
    }
    state->pitch_prev[0] = 100.0f; /* 5.1.3 initial values */
    state->pitch_prev[1] = 100.0f;
    state->xi_max = XI_MAX_FLOOR;
}

/*
 * The high-pass behind the noise estimator's periodicity test: a 6th-order
 * Butterworth at 100 Hz (bilinear transform, prewarped), as three sections
 * with Q 0.5176, 0.7071 and 1.9319, each b0 (1 - 2 z^-1 + z^-2) /
 * (1 + a1 z^-1 + a2 z^-2). It takes 50 and 60 Hz down 36 and 27 dB and 120 Hz
 * 0.5 dB.
 */
#define GATE_SECTIONS 3
static const float GATE_HP[GATE_SECTIONS][3] = {
    {0.92812044f, -1.85337538f, 0.859106383f}, /* b0, a1, a2 */
    {0.945976856f, -1.88903308f, 0.894874345f},
    {0.97858679f, -1.95415227f, 0.960194894f},
};

/* One section in transposed direct form II; z holds its two states. */
static float
gate_section(const float c[3], float* z, float x) {
    float y = (c[0] * x) + z[0];
    z[0] = (-2.0f * c[0] * x) - (c[1] * y) + z[1];
    z[1] = (c[0] * x) - (c[2] * y);
    return y;
}

void
mbe_analysis_push(struct mbe_analysis_state* state, const float input[MBE_ANALYSIS_FRAME]) {
    memmove(state->buf, state->buf + MBE_ANALYSIS_FRAME, MBE_ANALYSIS_HISTORY * sizeof(float));
    memmove(state->gate_buf, state->gate_buf + MBE_ANALYSIS_FRAME, MBE_ANALYSIS_HISTORY * sizeof(float));
    float* frame = state->buf + MBE_ANALYSIS_HISTORY;
    float* gate = state->gate_buf + MBE_ANALYSIS_HISTORY;
    for (int i = 0; i < MBE_ANALYSIS_FRAME; ++i) {
        /* eq (3): H(z) = (1 - z^-1) / (1 - 0.99 z^-1) */
        float y = input[i] - state->hp_in + (0.99f * state->hp_out);
        state->hp_in = input[i];
        state->hp_out = y;
        frame[i] = y;
        for (int k = 0; k < GATE_SECTIONS; ++k) {
            y = gate_section(GATE_HP[k], state->gate_filter + (2 * (size_t)k), y);
        }
        gate[i] = y;
    }
}

/* ---- Initial pitch, eqs (5)-(9) ---- */

struct pitch_work {
    float s_lpf[(2 * MBE_ANALYSIS_PITCH_HALF) + 1];
    float u[(2 * MBE_ANALYSIS_PITCH_HALF) + 1]; /* s_LPF(j) w_I^2(j) */
    float r[MBE_ANALYSIS_PITCH_HALF + 2];
    float error[MBE_ANALYSIS_PITCH_COUNT];
};

static float
pitch_of(int index) {
    return PITCH_MIN + (PITCH_STEP * (float)index);
}

#if defined(MBE_VF_WIDTH) && !defined(__FAST_MATH__) && !defined(_M_FP_FAST)
#define MBE_VECTOR_PITCH_LOWPASS 1

/* Independent output samples occupy lanes; every lane still adds the 21
 * filter taps in order. This also vectorizes IEEE builds, where an automatic
 * vectorizer cannot reassociate the surrounding energy reduction. Fast-math
 * builds let the compiler vectorize the filter and energy in a single pass. */
static void
pitch_lowpass(const struct mbe_analysis_tables* t, const float* buf, float* out) {
    const int span = (2 * MBE_ANALYSIS_PITCH_HALF) + 1;
    const float* input = buf + SAMPLE_CENTRE - MBE_ANALYSIS_PITCH_HALF;
    int j = 0;
    for (; j + MBE_VF_WIDTH <= span; j += MBE_VF_WIDTH) {
        mbe_vf acc = mbe_vf_set(0.0f);
        for (int k = -MBE_ANALYSIS_LPF_HALF; k <= MBE_ANALYSIS_LPF_HALF; ++k) {
            mbe_vf tap = mbe_vf_set(t->lpf[k + MBE_ANALYSIS_LPF_HALF]);
            acc = mbe_vf_add(acc, mbe_vf_mul(mbe_vf_load(input + j - k), tap));
        }
        mbe_vf_store(out + j, acc);
    }
    for (; j < span; ++j) {
        float acc = 0.0f;
        for (int k = -MBE_ANALYSIS_LPF_HALF; k <= MBE_ANALYSIS_LPF_HALF; ++k) {
            acc += input[j - k] * t->lpf[k + MBE_ANALYSIS_LPF_HALF];
        }
        out[j] = acc;
    }
}
#endif

/* eqs (7) and (9) over one history: buf or gate_buf of the state. */
static int
pitch_autocorrelation(const struct mbe_analysis_tables* t, const float* buf, mbe_acf_plan* acf, struct pitch_work* w,
                      float* energy) {
    const int span = (2 * MBE_ANALYSIS_PITCH_HALF) + 1;
    float e = 0.0f;
#if defined(MBE_VECTOR_PITCH_LOWPASS)
    pitch_lowpass(t, buf, w->s_lpf);
#endif
    for (int j = 0; j < span; ++j) {
#if defined(MBE_VECTOR_PITCH_LOWPASS)
        float acc = w->s_lpf[j];
#else
        const float* x = buf + SAMPLE_CENTRE - MBE_ANALYSIS_PITCH_HALF + j; /* sample j - 150 */
        float acc = 0.0f;
        for (int k = -MBE_ANALYSIS_LPF_HALF; k <= MBE_ANALYSIS_LPF_HALF; ++k) {
            acc += x[-k] * t->lpf[k + MBE_ANALYSIS_LPF_HALF]; /* eq (9) */
        }
        w->s_lpf[j] = acc;
#endif
        float wi2 = t->w_i[j] * t->w_i[j];
        w->u[j] = acc * wi2;
        e += acc * acc * wi2;
    }
    *energy = e;
    return mbe_acf_compute(acf, w->u, span, w->r, MBE_ANALYSIS_PITCH_HALF + 1); /* eq (7) */
}

/* Candidate i's pitch is (PITCH_MIN_GRID + i) / PITCH_GRID samples, so
 * quotient and remainder below select eq (8)'s bins and fraction exactly. */
#define PITCH_SPAN_GRID (MBE_ANALYSIS_PITCH_HALF * PITCH_GRID) /* r's last lag, in grid steps */

/* Pitch candidates in one group share their autocorrelation span, the n with
 * n * pitch <= MBE_ANALYSIS_PITCH_HALF. Every tap is then at or below r[HALF]
 * and the +1 interpolation tap at or below r[HALF + 1]. */
static void
pitch_error_group(const struct mbe_analysis_tables* t, struct pitch_work* w, float energy, int lo, int hi, int span) {
    float sums[MBE_ANALYSIS_PITCH_COUNT];
    for (int i = lo; i < hi; ++i) {
        sums[i] = w->r[0];
    }
    for (int n = 1; n <= span; ++n) {
        for (int i = lo; i < hi; ++i) {
            int lag = n * (PITCH_MIN_GRID + i); /* n * pitch_of(i), in grid steps */
            int base = lag / PITCH_GRID;
            float frac = (float)(lag % PITCH_GRID) / (float)PITCH_GRID;
            sums[i] += 2.0f * (((1.0f - frac) * w->r[base]) + (frac * w->r[base + 1])); /* r(-t) = r(t) */
        }
    }
    for (int i = lo; i < hi; ++i) {
        float p = pitch_of(i);
        float denominator = energy * (1.0f - (p * t->w_i4_sum));
        /* eq (5); window effects can make it slightly negative, as in the standard. */
        w->error[i] = (denominator > 1e-6f) ? (energy - (p * sums[i])) / denominator : 1.0f;
    }
}

static void
pitch_errors(const struct mbe_analysis_tables* t, struct pitch_work* w, float energy) {
    int lo = 0;
    while (lo < MBE_ANALYSIS_PITCH_COUNT) {
        /* span = floor(HALF / pitch); the group ends at the last candidate
         * with the same span, so no candidate's taps pass r[HALF + 1]. */
        int span = PITCH_SPAN_GRID / (PITCH_MIN_GRID + lo);
        int hi = (span > 0) ? (PITCH_SPAN_GRID / span) - PITCH_MIN_GRID + 1 : MBE_ANALYSIS_PITCH_COUNT;
        if (hi > MBE_ANALYSIS_PITCH_COUNT) {
            hi = MBE_ANALYSIS_PITCH_COUNT;
        }
        pitch_error_group(t, w, energy, lo, hi, span);
        lo = hi;
    }
}

static int
pitch_index(float p) {
    int i = (int)lrintf((p - PITCH_MIN) / PITCH_STEP);
    if (i < 0) {
        return 0;
    }
    return (i >= MBE_ANALYSIS_PITCH_COUNT) ? MBE_ANALYSIS_PITCH_COUNT - 1 : i;
}

static int
argmin_error(const float* error, int lo, int hi) {
    int best = lo;
    for (int i = lo + 1; i <= hi; ++i) {
        if (error[i] < error[best]) {
            best = i;
        }
    }
    return best;
}

/* The better of the two grid points bracketing a sub-multiple. The standard
 * takes the closest one, but on a half-sample grid that can sit 0.25 samples
 * off a near-perfect period, enough to fail the ratio tests and leave an
 * octave error (a steady 150 Hz voice tracked at 75 Hz). */
static int
submultiple_index(const float* error, float pitch) {
    int lo = (int)floorf((pitch - PITCH_MIN) / PITCH_STEP);
    if (lo < 0) {
        return 0;
    }
    if (lo >= MBE_ANALYSIS_PITCH_COUNT - 1) {
        return MBE_ANALYSIS_PITCH_COUNT - 1;
    }
    return (error[lo + 1] < error[lo]) ? lo + 1 : lo;
}

/* Sub-multiple checks, eqs (18)-(20), with CE_F(P) = 3 E(P): every
 * sub-multiple at or above the pitch floor, smallest first (5.1.4). In the
 * ratio tests CE_F(P0) is floored at the eq (20) threshold, so a near-perfect
 * fit at a multiple of the period cannot block every sub-multiple. */
static int
forward_pitch(const float* error) {
    int p0 = argmin_error(error, 0, MBE_ANALYSIS_PITCH_COUNT - 1);
    float ce0 = 3.0f * error[p0];
    if (ce0 < 0.05f) {
        ce0 = 0.05f;
    }
    for (int n = (int)(pitch_of(p0) / PITCH_MIN); n >= 2; --n) {
        int i = submultiple_index(error, pitch_of(p0) / (float)n);
        float ce = 3.0f * error[i];
        if ((ce <= 0.85f && ce <= 1.7f * ce0) || (ce <= 0.4f && ce <= 3.5f * ce0) || ce <= 0.05f) {
            return i;
        }
    }
    return p0;
}

static int
track_pitch(const struct mbe_analysis_state* s, const float* error) {
    int forward = forward_pitch(error);
    if (!s->trusted) {
        return forward;
    }
    /* Look-back, eqs (10)-(12), then the decision rules (21)-(23). */
    int lo = pitch_index(0.8f * s->pitch_prev[0]);
    int hi = pitch_index(1.2f * s->pitch_prev[0]);
    int back = argmin_error(error, lo, hi);
    float ce_back = error[back] + s->error_prev[0] + s->error_prev[1];
    if (ce_back <= 0.48f || ce_back <= 3.0f * error[forward]) {
        return back;
    }
    return forward;
}

/* ---- Spectrum and harmonic fitting, eqs (24)-(29) ---- */

struct spectrum {
    float re[ANALYSIS_BINS + 1];
    float im[ANALYSIS_BINS + 1];
};

static int
analysis_spectrum(const struct mbe_analysis_tables* t, const struct mbe_analysis_state* s, mbe_fft_plan* fft,
                  struct spectrum* out) {
    float in[ANALYSIS_FFT] = {0};
    float packed[ANALYSIS_FFT];
    /* eq (29): s(n) w_R(n) for n = -110..110 placed with n = 0 at index 0. */
    for (int n = -MBE_ANALYSIS_WR_HALF; n <= MBE_ANALYSIS_WR_HALF; ++n) {
        in[(n + ANALYSIS_FFT) % ANALYSIS_FFT] = s->buf[SAMPLE_CENTRE + n] * t->w_r[n + MBE_ANALYSIS_WR_HALF];
    }
    int status = mbe_fft_forward_real(fft, in, packed);
    if (status < 0) {
        return status;
    }
    out->re[0] = packed[0];
    out->im[0] = 0.0f;
    out->re[ANALYSIS_BINS] = packed[1];
    out->im[ANALYSIS_BINS] = 0.0f;
    for (int m = 1; m < ANALYSIS_BINS; ++m) {
        out->re[m] = packed[2 * (size_t)m];
        out->im[m] = packed[(2 * (size_t)m) + 1];
    }
    return 0;
}

static float
window_at(const struct mbe_analysis_tables* t, float offset_bins) {
    int i =
        (int)lrintf(offset_bins * (float)MBE_ANALYSIS_WTAB_STEPS) + (MBE_ANALYSIS_WTAB_BINS * MBE_ANALYSIS_WTAB_STEPS);
    return (i < 0 || i >= MBE_ANALYSIS_WTAB_LEN) ? 0.0f : t->w_r_dtft[i];
}

/* Bins one harmonic band can span: refinement goes down to a pitch of
 * PITCH_MIN - 9/8, where f0 is 13.6 bins, so a band covers at most 14 bins. */
#define BAND_BINS_MAX 16

struct band_fit {
    float a_re; /* A_l of eq (28) */
    float a_im;
    float energy;           /* sum |S|^2 over the band */
    float window_energy;    /* sum W_R^2 over the band */
    float error;            /* sum |S - A W|^2 over the band */
    int bins;               /* bins fitted, from lo */
    float w[BAND_BINS_MAX]; /* W_R at each fitted bin */
};

/* Least-squares fit of one harmonic over bins [lo, hi). */
static void
fit_band(const struct mbe_analysis_tables* t, const struct spectrum* sp, float centre, int lo, int hi,
         struct band_fit* fit) {
    float c_re = 0.0f;
    float c_im = 0.0f;
    float q = 0.0f;
    float energy = 0.0f;
    fit->bins = 0;
    for (int m = lo; m < hi && m <= ANALYSIS_BINS && fit->bins < BAND_BINS_MAX; ++m) {
        float w = window_at(t, (float)m - centre);
        fit->w[fit->bins++] = w;
        c_re += sp->re[m] * w;
        c_im += sp->im[m] * w;
        q += w * w;
        energy += (sp->re[m] * sp->re[m]) + (sp->im[m] * sp->im[m]);
    }
    fit->energy = energy;
    fit->window_energy = q;
    fit->a_re = (q > 0.0f) ? c_re / q : 0.0f;
    fit->a_im = (q > 0.0f) ? c_im / q : 0.0f;
    float explained = (q > 0.0f) ? ((c_re * c_re) + (c_im * c_im)) / q : 0.0f;
    float error = energy - explained;
    fit->error = (error > 0.0f) ? error : 0.0f; /* guard rounding */
}

static void
band_edges(float f0_bins, int l, int* lo, int* hi) {
    *lo = (int)ceilf(((float)l - 0.5f) * f0_bins);
    *hi = (int)ceilf(((float)l + 0.5f) * f0_bins);
    if (*lo < 0) {
        *lo = 0;
    }
    if (*hi > ANALYSIS_BINS + 1) {
        *hi = ANALYSIS_BINS + 1;
    }
}

/* eq (24): fit error over bins first .. upper for one candidate pitch. */
static float
refinement_error(const struct mbe_analysis_tables* t, const struct spectrum* sp, float pitch, int first) {
    const float f0_bins = (float)ANALYSIS_FFT / pitch;
    const int top_harmonic = (int)((0.9254f * pitch * 0.5f) - 0.5f);
    const int upper = (int)((float)top_harmonic * f0_bins);
    float total = 0.0f;
    for (int l = 1; l <= top_harmonic; ++l) {
        int lo;
        int hi;
        band_edges(f0_bins, l, &lo, &hi);
        if (hi <= first || lo > upper) {
            continue;
        }
        struct band_fit fit;
        fit_band(t, sp, (float)l * f0_bins, lo, hi, &fit);
        /* The fit's window samples cover bins lo .. lo + bins - 1. */
        for (int k = 0; k < fit.bins && k < BAND_BINS_MAX && lo + k <= upper; ++k) {
            const int m = lo + k;
            if (m < first) {
                continue;
            }
            float w = fit.w[k];
            float dr = sp->re[m] - (fit.a_re * w);
            float di = sp->im[m] - (fit.a_im * w);
            total += (dr * dr) + (di * di);
        }
    }
    return total;
}

/* Octave-down check. When the pitch is twice the true period, its odd
 * harmonics fall midway between the real ones, about one window main-lobe
 * width from each, and hold almost no energy; at the right pitch they are real
 * harmonics. E(P) itself can prefer the doubled period on very steady voices
 * (a 150 Hz vowel: E 0.063 at the period, 0.017 at twice it), and the ratio
 * tests of eqs (18)-(19) then miss by a hair. */
#define OCTAVE_ODD_RATIO 0.1f

static float
line_energy(const struct spectrum* sp, float centre) {
    int m = (int)lrintf(centre);
    if (m < 1 || m > ANALYSIS_BINS) {
        return 0.0f;
    }
    return (sp->re[m] * sp->re[m]) + (sp->im[m] * sp->im[m]);
}

static int
is_octave_low(const struct spectrum* sp, float pitch) {
    const float f0_bins = (float)ANALYSIS_FFT / pitch;
    float odd = 0.0f;
    float even = 0.0f;
    for (int h = 1; (float)h * f0_bins < (float)ANALYSIS_BINS; ++h) {
        float e = line_energy(sp, (float)h * f0_bins);
        if (h & 1) {
            odd += e;
        } else {
            even += e;
        }
    }
    return even > 0.0f && odd < OCTAVE_ODD_RATIO * even;
}

/* The standard scores bins from 1.56 kHz up, where pitch errors show most.
 * When that region holds almost no energy (a low tone, a muffled input) the
 * score would only rank rounding noise, so start from the first bin. */
static int
refinement_first_bin(const struct spectrum* sp) {
    float total = 0.0f;
    float high = 0.0f;
    for (int m = 1; m <= ANALYSIS_BINS; ++m) {
        float e = (sp->re[m] * sp->re[m]) + (sp->im[m] * sp->im[m]);
        total += e;
        high += (m >= REFINE_FIRST_BIN) ? e : 0.0f;
    }
    return (high >= 0.01f * total) ? REFINE_FIRST_BIN : 1;
}

static float
refine_pitch(const struct mbe_analysis_tables* t, const struct spectrum* sp, float initial) {
    const int first = refinement_first_bin(sp);
    float best_pitch = initial;
    float best_error = 0.0f;
    for (int c = 0; c < REFINE_CANDIDATES; ++c) {
        float pitch = initial + ((float)((2 * c) - 9) / 8.0f); /* initial - 9/8 .. + 9/8 */
        float error = refinement_error(t, sp, pitch, first);
        if (c == 0 || error < best_error) {
            best_error = error;
            best_pitch = pitch;
        }
    }
    return best_pitch;
}

/* ---- Magnitudes (US 5,701,390 eqs 3-4) and fit errors ---- */

/* Tapered band weight G: 1 inside the band, falling linearly to 0 across one
 * bin around each band edge. */
static float
band_weight(float distance, float half_width) {
    float d = fabsf(distance);
    if (d <= half_width - 0.5f) {
        return 1.0f;
    }
    if (d >= half_width + 0.5f) {
        return 0.0f;
    }
    return 0.5f - (d - half_width);
}

static float
harmonic_magnitude(const struct mbe_analysis_tables* t, const struct spectrum* sp, float centre, float f0_bins) {
    const float half = 0.5f * f0_bins;
    int lo = (int)floorf(centre - half - 0.5f);
    int hi = (int)ceilf(centre + half + 0.5f);
    if (lo < 0) {
        lo = 0;
    }
    if (hi > ANALYSIS_BINS) {
        hi = ANALYSIS_BINS;
    }
    float sum = 0.0f;
    float weight_sum = 0.0f;
    for (int m = lo; m <= hi; ++m) {
        float g = band_weight((float)m - centre, half);
        sum += g * ((sp->re[m] * sp->re[m]) + (sp->im[m] * sp->im[m]));
        weight_sum += g;
    }
    if (weight_sum <= 0.0f) {
        return -1.0f; /* band entirely above Nyquist */
    }
    return sqrtf(sum / ((float)ANALYSIS_FFT * t->w_r_sq_sum));
}

/* The normalized fit error white noise leaves in this band, expected: 1 minus
 * the share of its energy the window-shaped line explains, w' C w / (w'w
 * bins C(0)) with C the noise covariance between bins. */
static float
noise_fit_error(const struct mbe_analysis_tables* t, const struct band_fit* fit) {
    float explained = 0.0f;
    for (int a = 0; a < fit->bins; ++a) {
        for (int b = 0; b < fit->bins; ++b) {
            explained += fit->w[a] * fit->w[b] * t->noise_cov[(a > b) ? a - b : b - a];
        }
    }
    const float total = fit->window_energy * (float)fit->bins * t->noise_cov[0];
    return (total > 0.0f) ? 1.0f - (explained / total) : 1.0f;
}

static void
harmonic_measures(const struct mbe_analysis_tables* t, const struct spectrum* sp, struct mbe_analysis_result* r) {
    const float f0_bins = r->f0 * (float)ANALYSIS_FFT;
    float last = 0.0f;
    r->harmonics = 0;
    for (int l = 1; l <= MBE_ANALYSIS_HARMONICS; ++l) {
        float centre = (float)l * f0_bins;
        float m = harmonic_magnitude(t, sp, centre, f0_bins);
        /* Past Nyquist the decoder still expects a value; repeat the last. */
        r->magnitude[l] = (m < 0.0f) ? last : m;
        last = r->magnitude[l];
        int lo;
        int hi;
        band_edges(f0_bins, l, &lo, &hi);
        if (centre < (float)ANALYSIS_BINS && lo < hi) {
            struct band_fit fit;
            fit_band(t, sp, centre, lo, hi, &fit);
            r->fit_error[l] = (fit.energy > 0.0f) ? fit.error / fit.energy : 1.0f;
            r->fit_energy[l] = fit.energy;
            r->noise_fit_error[l] = noise_fit_error(t, &fit);
            /* eq (43): band energy relative to the window's over the same bins. */
            r->voiced_magnitude[l] = (fit.window_energy > 0.0f) ? sqrtf(fit.energy / fit.window_energy) : 0.0f;
            /* eq (44): mean band energy per bin, scaled by 1 / sum w_R. */
            r->noise_magnitude[l] = (fit.bins > 0) ? sqrtf(fit.energy / (float)fit.bins) / t->w_r_sum : 0.0f;
            r->harmonics = l;
        } else {
            r->fit_error[l] = 1.0f;
            r->noise_fit_error[l] = 1.0f;
            r->voiced_magnitude[l] = r->voiced_magnitude[l - 1];
            r->noise_magnitude[l] = r->noise_magnitude[l - 1];
        }
    }
}

/* eqs (38)-(42), with energies on the 16-bit input scale. */
static float
energy_factor(const struct mbe_analysis_tables* t, struct mbe_analysis_state* s, const struct spectrum* sp) {
    const float norm = 1.0f / (t->w_r_sum * t->w_r_sum);
    float lf = 0.0f;
    float hf = 0.0f;
    for (int m = 0; m <= ANALYSIS_BINS; ++m) {
        float e = ((sp->re[m] * sp->re[m]) + (sp->im[m] * sp->im[m])) * norm;
        if (m < 64) {
            lf += e;
        } else {
            hf += e;
        }
    }
    float xi0 = lf + hf;
    if (xi0 > s->xi_max) {
        s->xi_max = (0.5f * s->xi_max) + (0.5f * xi0);
    } else {
        float decayed = (0.99f * s->xi_max) + (0.01f * xi0);
        s->xi_max = (decayed > XI_MAX_FLOOR) ? decayed : XI_MAX_FLOOR;
    }
    float factor = ((0.0025f * s->xi_max) + xi0) / ((0.01f * s->xi_max) + xi0);
    if (lf < 5.0f * hf) {
        factor *= sqrtf(lf / (5.0f * hf));
    }
    return factor;
}

/* ---- Background noise ---- */

/*
 * At the lowest pitches a harmonic band spans about two bins, and the voicing
 * test passes noise as voiced (#149: 60 Hz mains hum read as a voice at
 * 65 Hz, every band voiced). The analysis therefore estimates the background
 * noise per bin, and the voicing decisions treat a 500 Hz column at such a
 * pitch as unvoiced when it holds little more than that estimate
 * (mbe_analysis_noise_dominated()). The magnitudes are not changed. DVSI's
 * AMBE-3000 encoders suppress background noise by default
 * (docs/dvsi-parity.md, "Encoder voicing over background noise"); this only
 * borrows an estimate of it to keep it from being voiced.
 *
 * The estimate per bin is the lowest smoothed power of the last 48 or so
 * noise-like frames, doubled (after minimum statistics, R. Martin, IEEE Trans.
 * Speech Audio Process. 9(5), 2001), and at most their mean. A frame is
 * noise-like when it does not repeat in the voice band (E at the tracked pitch,
 * on the input high-passed at 100 Hz, above the periodicity gate; mains hum
 * without its fundamental repeats at no pitch the search sends), its total
 * power is within 6 dB of the previous frame's, which keeps most of a speech
 * onset out, and it is quiet (under -30 dBFS) or 15 dB under the loudest
 * sound tracked so far. A quiet voice whose harmonics are weak beside its
 * fundamental is periodic in the voice band only faintly and can be taken for
 * noise as hum is; a loud one only right after a sound 15 dB louder.
 */
#define NOISE_SMOOTH        0.5f     /* first-order smoothing of the spectrum power */
#define NOISE_WINDOW_FRAMES 6        /* noise-like frames per sub-window; 8 sub-windows */
#define NOISE_IDLE_FRAMES   50       /* frames without a noise-like one before the estimate is forgotten */
#define NOISE_BIAS          2.0f     /* noise estimate over the minimum */
#define NOISE_STEADY        3.98107f /* 6 dB, as a power ratio */
#define NOISE_BELOW_PEAK    0.0316f  /* 15 dB under the tracked peak energy, eq (41) */
#define NOISE_CEILING       9.1e5f   /* xi_0 of white noise at -30 dBFS */

/* The estimate from the sub-windows' minima: NOISE_BIAS times the lowest,
 * but no more than the mean of the noise-like frames they took. With few such
 * frames their minimum is close to their mean, and doubling it would put
 * speech under the estimate. */
static void
noise_estimate(struct mbe_analysis_state* s) {
    int taken = 0;
    for (int u = 0; u <= s->noise_windows && u <= MBE_ANALYSIS_NOISE_WINDOWS; ++u) {
        taken += s->noise_taken[u];
    }
    for (int m = 0; m < MBE_ANALYSIS_BINS; ++m) {
        float low = s->noise_cur[m];
        float sum = 0.0f;
        for (int u = 0; u <= s->noise_windows && u <= MBE_ANALYSIS_NOISE_WINDOWS; ++u) {
            sum += s->noise_sum[u][m];
        }
        for (int u = 0; u < s->noise_windows; ++u) {
            float v = s->noise_min[u][m];
            if (v >= 0.0f && (low < 0.0f || v < low)) {
                low = v;
            }
        }
        s->noise[m] = (low > 0.0f && taken > 0) ? fminf(NOISE_BIAS * low, sum / (float)taken) : 0.0f;
    }
}

/* Forget the noise estimate. */
static void
noise_clear(struct mbe_analysis_state* s) {
    memset(s->noise_sum, 0, sizeof(s->noise_sum));
    memset(s->noise_taken, 0, sizeof(s->noise_taken));
    memset(s->noise, 0, sizeof(s->noise));
    for (int m = 0; m < MBE_ANALYSIS_BINS; ++m) {
        s->noise_cur[m] = -1.0f;
    }
    s->noise_count = 0;
    s->noise_windows = 0;
}

/* Take the frame's smoothed power into the open sub-window, and close it
 * after NOISE_WINDOW_FRAMES noise-like frames. */
static void
noise_take(struct mbe_analysis_state* s) {
    for (int m = 0; m < MBE_ANALYSIS_BINS; ++m) {
        if (s->noise_cur[m] < 0.0f || s->noise_power[m] < s->noise_cur[m]) {
            s->noise_cur[m] = s->noise_power[m];
        }
        s->noise_sum[0][m] += s->noise_power[m];
    }
    s->noise_taken[0]++;
    if (++s->noise_count < NOISE_WINDOW_FRAMES) {
        return;
    }
    memmove(s->noise_min[1], s->noise_min[0], (MBE_ANALYSIS_NOISE_WINDOWS - 1) * sizeof(s->noise_min[0]));
    memcpy(s->noise_min[0], s->noise_cur, sizeof(s->noise_cur));
    s->noise_windows += (s->noise_windows < MBE_ANALYSIS_NOISE_WINDOWS);
    s->noise_count = 0;
    memmove(s->noise_sum[1], s->noise_sum[0], MBE_ANALYSIS_NOISE_WINDOWS * sizeof(s->noise_sum[0]));
    memmove(s->noise_taken + 1, s->noise_taken, MBE_ANALYSIS_NOISE_WINDOWS * sizeof(s->noise_taken[0]));
    memset(s->noise_sum[0], 0, sizeof(s->noise_sum[0]));
    s->noise_taken[0] = 0;
    for (int m = 0; m < MBE_ANALYSIS_BINS; ++m) {
        s->noise_cur[m] = -1.0f; /* nothing taken yet */
    }
}

/* Smooth the frame's spectrum power into the state and, when the frame is
 * noise-like, take it into the open sub-window and refresh the estimate.
 * Sub-windows close every NOISE_WINDOW_FRAMES noise-like frames, so the
 * minimum spans the last 50 or so of them however much speech lies between.
 * Once NOISE_IDLE_FRAMES pass without one (sustained voiced sound, or loud
 * noise that stopped under a vowel) the estimate is forgotten, so a stale one
 * suppresses at most that long. */
static void
noise_track(const struct mbe_analysis_tables* t, struct mbe_analysis_state* s, const struct spectrum* sp,
            int aperiodic) {
    float total = 0.0f;
    for (int m = 0; m < MBE_ANALYSIS_BINS; ++m) {
        float p = (sp->re[m] * sp->re[m]) + (sp->im[m] * sp->im[m]);
        s->noise_power[m] =
            (s->noise_frames == 0) ? p : (NOISE_SMOOTH * s->noise_power[m]) + ((1.0f - NOISE_SMOOTH) * p);
        total += p;
    }
    /* Within 6 dB of the previous frame, digital silence after silence included. */
    const int steady = (total <= NOISE_STEADY * s->noise_prev_total && s->noise_prev_total <= NOISE_STEADY * total);
    s->noise_prev_total = total;
    /* Background is quiet, or lies under the signal around it: the frame's
     * energy, as xi_0 of eq (40), under that of white noise at -30 dBFS (DVSI's
     * encoder takes noise at -34 dBFS for noise from the start of a stream),
     * or 15 dB under xi_max of eq (41) (its value before this frame), as
     * background lies under speech (18 dB and more in DVSI's noisy vectors).
     * A loud sustained sound does not become its own noise, unless a sound
     * 15 dB louder came shortly before (xi_max falls about 0.04 dB a frame).
     * The first frame sees mostly the zeroed history. */
    const float xi0 = total / (t->w_r_sum * t->w_r_sum);
    const int quiet = xi0 < fmaxf(NOISE_CEILING, NOISE_BELOW_PEAK * s->xi_max);
    const int take = aperiodic && steady && quiet && s->noise_frames > 0;
    if (s->noise_frames < NOISE_IDLE_FRAMES) {
        s->noise_frames++;
    }
    if (!take) {
        if (s->noise_idle < NOISE_IDLE_FRAMES) {
            s->noise_idle++;
            if (s->noise_idle == NOISE_IDLE_FRAMES) {
                noise_clear(s);
            }
        }
        return;
    }
    s->noise_idle = 0;
    noise_take(s);
    noise_estimate(s);
}

/* The noise estimate summed over each analysed harmonic's band. */
static void
noise_bands(const struct mbe_analysis_state* s, struct mbe_analysis_result* r) {
    const float f0_bins = r->f0 * (float)ANALYSIS_FFT;
    for (int l = 1; l <= r->harmonics; ++l) {
        int lo;
        int hi;
        band_edges(f0_bins, l, &lo, &hi);
        float noise = 0.0f;
        for (int m = lo; m < hi && m < lo + BAND_BINS_MAX && m < MBE_ANALYSIS_BINS; ++m) {
            noise += s->noise[m];
        }
        r->noise_energy[l] = noise;
    }
}

/*
 * E at the tracked pitch candidate on the history high-passed at 100 Hz
 * (gate_buf), for the noise estimator; w is reused as scratch.
 */
static int
gate_error(const struct mbe_analysis_tables* t, const struct mbe_analysis_state* s, mbe_acf_plan* acf,
           struct pitch_work* w, int index, float* error) {
    float energy;
    int status = pitch_autocorrelation(t, s->gate_buf, acf, w, &energy);
    if (status < 0) {
        return status;
    }
    pitch_error_group(t, w, energy, index, index + 1, PITCH_SPAN_GRID / (PITCH_MIN_GRID + index));
    *error = w->error[index];
    return 0;
}

int
mbe_analysis_frame(const struct mbe_analysis_tables* tables, struct mbe_analysis_state* state, mbe_fft_plan* fft,
                   mbe_acf_plan* acf, struct mbe_analysis_result* result) {
    if (tables == NULL || state == NULL || fft == NULL || acf == NULL || result == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    struct pitch_work work;
    float energy;
    int acf_status = pitch_autocorrelation(tables, state->buf, acf, &work, &energy);
    if (acf_status < 0) {
        return acf_status;
    }
    pitch_errors(tables, &work, energy);
    int initial = track_pitch(state, work.error);
    float initial_pitch = pitch_of(initial);

    struct spectrum sp;
    int status = analysis_spectrum(tables, state, fft, &sp);
    if (status < 0) {
        return status;
    }
    if (initial_pitch * 0.5f >= PITCH_MIN && is_octave_low(&sp, initial_pitch)) {
        initial = submultiple_index(work.error, initial_pitch * 0.5f);
        initial_pitch = pitch_of(initial);
    }
    const float tracked_error = work.error[initial];
    float band_error;
    status = gate_error(tables, state, acf, &work, initial, &band_error);
    if (status < 0) {
        return status;
    }
    float pitch = refine_pitch(tables, &sp, initial_pitch);

    memset(result, 0, sizeof(*result));
    result->f0 = 1.0f / pitch;
    result->pitch_error = tracked_error;
    harmonic_measures(tables, &sp, result);
    noise_track(tables, state, &sp, band_error > VOICING_ERROR_MAX);
    noise_bands(state, result);
    result->energy_factor = energy_factor(tables, state, &sp);

    state->pitch_prev[1] = state->pitch_prev[0];
    state->pitch_prev[0] = initial_pitch;
    state->error_prev[1] = state->error_prev[0];
    state->error_prev[0] = tracked_error;
    state->trusted = tracked_error <= PERIODIC_ERROR_MAX;
    return 0;
}

/* Theta of eq (37) for 1-based band k at fundamental w0 (radians per sample)
 * and energy factor M(xi), before the periodicity gate. */
static float
voicing_threshold(int prev_voiced, int k, float w0, float factor) {
    float base = prev_voiced ? 0.5625f : 0.45f;
    return base * (1.0f - (0.3096f * (float)(k - 1) * w0)) * factor;
}

float
mbe_analysis_column_threshold(const struct mbe_analysis_state* state, const struct mbe_analysis_result* result,
                              int column) {
    /* eq (37) at omega0 = 0.1309, where three harmonics span one 500 Hz
     * column (US 8,595,002). The standard unvoices bands above the first when
     * E(P_I) > 0.5. Here a frame with E(P_I) > 0.4 is unvoiced in every
     * column: on DVSI's D-STAR test vectors this raises agreement with DVSI's
     * voicing on clean speech (0.83 to 0.85 in the lowest 1 kHz band) and
     * brings voicing in noisy speech down to DVSI's, where the standard's
     * gate voiced 0.1-0.15 more of each band. */
    if (result->pitch_error > VOICING_ERROR_MAX) {
        return 0.0f;
    }
    return voicing_threshold(state->columns_prev[column - 1], column, 0.1309f, result->energy_factor);
}

/* 500 Hz column of harmonic l at the analysed fundamental. */
static int
analysis_column(const struct mbe_analysis_result* result, int l) {
    int column = (int)((float)l * 16.0f * result->f0);
    return (column < MBE_ANALYSIS_COLUMNS) ? column : MBE_ANALYSIS_COLUMNS - 1;
}

int
mbe_analysis_band_voicing(const struct mbe_analysis_result* result, int L, const unsigned char prev[MBE_ANALYSIS_BANDS],
                          unsigned char bands[MBE_ANALYSIS_BANDS]) {
    float error[MBE_ANALYSIS_BANDS + 1] = {0};
    float energy[MBE_ANALYSIS_BANDS + 1] = {0};
    float column_energy[MBE_ANALYSIS_COLUMNS] = {0};
    float column_noise[MBE_ANALYSIS_COLUMNS] = {0};
    float column_noise_error[MBE_ANALYSIS_COLUMNS] = {0};
    if (L < 1) {
        L = 1;
    }
    if (L > MBE_ANALYSIS_HARMONICS) {
        L = MBE_ANALYSIS_HARMONICS;
    }
    const int K = mbe_analysis_band_of(L); /* eq (34) */
    for (int l = 1; l <= L && l <= result->harmonics; l++) {
        int k = mbe_analysis_band_of(l);
        energy[k] += result->fit_energy[l];
        error[k] += result->fit_error[l] * result->fit_energy[l];
        int column = analysis_column(result, l);
        column_energy[column] += result->fit_energy[l];
        column_noise[column] += result->noise_energy[l];
        column_noise_error[column] += result->noise_fit_error[l] * result->fit_energy[l];
    }
    const float w0 = 2.0f * (float)M_PI * result->f0;
    for (int k = 1; k <= MBE_ANALYSIS_BANDS; k++) {
        float theta = 0.0f;
        /* The gate of mbe_analysis_column_threshold(); see the header. */
        if (k <= K && result->pitch_error <= VOICING_ERROR_MAX) {
            theta = voicing_threshold(prev[k - 1], k, w0, result->energy_factor);
        }
        /* A band in a column of mostly noise is unvoiced: the column of its
         * middle harmonic (band K runs from harmonic 3K - 2 to L). */
        const int first = (3 * k) - 2;
        const int last = (k < K) ? 3 * k : L;
        int column = analysis_column(result, (first + last) / 2);
        if (mbe_analysis_noise_dominated(column_energy[column], column_noise[column], column_noise_error[column])) {
            theta = 0.0f;
        }
        /* eqs (35)-(36); an empty band is unvoiced. */
        bands[k - 1] = (unsigned char)(energy[k] > 0.0f && error[k] < theta * energy[k]);
    }
    return K;
}

void
mbe_analysis_commit(struct mbe_analysis_state* state, const unsigned char columns[MBE_ANALYSIS_COLUMNS]) {
    if (columns != NULL) {
        memcpy(state->columns_prev, columns, sizeof(state->columns_prev));
    }
}
