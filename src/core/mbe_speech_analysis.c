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
 *  - Local fractional-period refinement of the initial pitch.
 *  - Per-harmonic periodicity and the energy factor M(xi), eqs (38)-(42), and the
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
 *  - Voicing confidence measures agreement between successive periods,
 *    rather than fitting an amplitude to a few correlated spectral bins.
 *    Fine pitch locates a bounded period minimum before interpolation, with
 *    low-frequency emphasis while retaining upper-band timing.
 *  - Magnitudes use the voicing-independent estimate of US 5,701,390
 *    (eqs 3-4, expired), which the AMBE decoders' unvoiced scaling expects.
 *  - The periodicity gate of eq (37) is 0.4 and covers every column and band;
 *    see mbe_analysis_column_threshold() and mbe_analysis_band_voicing().
 */
#include "mbe_speech_analysis.h"

#include <fenv.h>
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
#define XI_MAX_FLOOR       20000.0f
#define PERIODIC_ERROR_MAX 0.5f /* periodic enough to trust the pitch track */
#define VOICING_ERROR_MAX  0.4f /* see mbe_analysis_column_threshold() */
#define SAMPLE_CENTRE      (MBE_ANALYSIS_HISTORY - 1)

/* Fixed windows and filter coefficients are generated offline. Copying their
 * float representations preserves the accepted coefficients without startup
 * double-precision math or per-context transcendental evaluation. */
#include "mbe_analysis_tables.h"
#include "mbe_analysis_window_index.h"

void
mbe_analysis_init_tables(struct mbe_analysis_tables* tables) {
    if (tables != NULL) {
        *tables = analysis_tables;
    }
}

void
mbe_analysis_reset(struct mbe_analysis_state* state) {
    if (state == NULL) {
        return;
    }
    memset(state, 0, sizeof(*state));
    state->pitch_prev[0] = 100.0f; /* 5.1.3 initial values */
    state->pitch_prev[1] = 100.0f;
    state->xi_max = XI_MAX_FLOOR;
}

void
mbe_analysis_push(struct mbe_analysis_state* state, const float input[MBE_ANALYSIS_FRAME]) {
    memmove(state->buf, state->buf + MBE_ANALYSIS_FRAME, MBE_ANALYSIS_HISTORY * sizeof(float));
    float* frame = state->buf + MBE_ANALYSIS_HISTORY;
    for (int i = 0; i < MBE_ANALYSIS_FRAME; ++i) {
        /* eq (3): H(z) = (1 - z^-1) / (1 - 0.99 z^-1) */
        float y = input[i] - state->hp_in + (0.99f * state->hp_out);
        state->hp_in = input[i];
        state->hp_out = y;
        frame[i] = y;
    }
}

/* ---- Initial pitch, eqs (5)-(9) ---- */

struct pitch_work {
    float refinement_samples[(2 * MBE_ANALYSIS_PITCH_HALF) + 1]; /* low-pass plus DC-filtered input */
    float u[(2 * MBE_ANALYSIS_PITCH_HALF) + 1];                  /* s_LPF(j) w_I^2(j) */
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
pitch_lowpass(const struct mbe_analysis_tables* t, const struct mbe_analysis_state* s, float* out) {
    const int span = (2 * MBE_ANALYSIS_PITCH_HALF) + 1;
    const float* input = s->buf + SAMPLE_CENTRE - MBE_ANALYSIS_PITCH_HALF;
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

static int
pitch_autocorrelation(const struct mbe_analysis_tables* t, const struct mbe_analysis_state* s, mbe_acf_plan* acf,
                      struct pitch_work* w, float* energy) {
    const int span = (2 * MBE_ANALYSIS_PITCH_HALF) + 1;
    float e = 0.0f;
#if defined(MBE_VECTOR_PITCH_LOWPASS)
    pitch_lowpass(t, s, w->refinement_samples);
#endif
    for (int j = 0; j < span; ++j) {
#if defined(MBE_VECTOR_PITCH_LOWPASS)
        float acc = w->refinement_samples[j];
#else
        const float* x = s->buf + SAMPLE_CENTRE - MBE_ANALYSIS_PITCH_HALF + j; /* sample j - 150 */
        float acc = 0.0f;
        for (int k = -MBE_ANALYSIS_LPF_HALF; k <= MBE_ANALYSIS_LPF_HALF; ++k) {
            acc += x[-k] * t->lpf[k + MBE_ANALYSIS_LPF_HALF]; /* eq (9) */
        }
#endif
        /* Coarse pitch keeps acc alone. Refinement retains upper-band
         * timing with low-frequency emphasis, without another filter. */
        w->refinement_samples[j] = acc + s->buf[SAMPLE_CENTRE - MBE_ANALYSIS_PITCH_HALF + j];
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

/* ---- Spectrum and harmonic energies ---- */

struct spectrum {
    /* All analysis consumers need power, not complex phase. Form it once. */
    float power[ANALYSIS_BINS + 1];
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
    out->power[0] = packed[0] * packed[0] + 0.0f;
    out->power[ANALYSIS_BINS] = packed[1] * packed[1] + 0.0f;
    for (int m = 1; m < ANALYSIS_BINS; ++m) {
        float re = packed[2 * (size_t)m];
        float im = packed[(2 * (size_t)m) + 1];
        out->power[m] = re * re + im * im;
    }
    return 0;
}

static float
window_at(const struct mbe_analysis_tables* t, float offset_bins) {
    int i =
        (int)lrintf(offset_bins * (float)MBE_ANALYSIS_WTAB_STEPS) + (MBE_ANALYSIS_WTAB_BINS * MBE_ANALYSIS_WTAB_STEPS);
    return (i < 0 || i >= MBE_ANALYSIS_WTAB_LEN) ? 0.0f : t->w_r_dtft[i];
}

struct band_energy {
    float energy;
    float window_energy;
    int bins;
};

/* Signal and window energies over the same bins, for eqs (43)-(44). */
static void
measure_band_energy(const struct mbe_analysis_tables* t, const struct spectrum* sp, float centre, int lo, int hi,
                    struct band_energy* band) {
    float window_energy = 0.0f;
    float energy = 0.0f;
    band->bins = 0;
    for (int m = lo; m < hi && m <= ANALYSIS_BINS; ++m) {
        float w = window_at(t, (float)m - centre);
        ++band->bins;
        window_energy += w * w;
        energy += sp->power[m];
    }
    band->energy = energy;
    band->window_energy = window_energy;
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
    return sp->power[m];
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

/* ---- Magnitudes (US 5,701,390 eqs 3-4) ---- */

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
        sum += g * sp->power[m];
        weight_sum += g;
    }
    if (weight_sum <= 0.0f) {
        return -1.0f; /* band entirely above Nyquist */
    }
    return sqrtf(sum / ((float)ANALYSIS_FFT * t->w_r_sq_sum));
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
            struct band_energy band;
            measure_band_energy(t, sp, centre, lo, hi, &band);
            r->fit_energy[l] = band.energy;
            /* eq (43): band energy relative to the window's over the same bins. */
            r->voiced_magnitude[l] = (band.window_energy > 0.0f) ? sqrtf(band.energy / band.window_energy) : 0.0f;
            /* eq (44): mean band energy per bin, scaled by 1 / sum w_R. */
            r->noise_magnitude[l] = (band.bins > 0) ? sqrtf(band.energy / (float)band.bins) / t->w_r_sum : 0.0f;
            r->harmonics = l;
        } else {
            r->fit_error[l] = 1.0f;
            r->voiced_magnitude[l] = r->voiced_magnitude[l - 1];
            r->noise_magnitude[l] = r->noise_magnitude[l - 1];
        }
    }
}

/* ---- Period-based pitch refinement and voicing ---- */

/* The supported pitch range gives half=84..110 in both window builders.
 * Keep the original calculation when the caller changes the rounding mode. */
static const unsigned char*
window_indices(int half) {
    return fegetround() == FE_TONEAREST
               ? mbe_window_indices + (half * half - MBE_WINDOW_INDEX_HALF_MIN * MBE_WINDOW_INDEX_HALF_MIN)
               : NULL;
}

/* Refine only the period used by harmonic analysis; the coarse pitch and its
 * error still drive the tracker. Integer-lag comparisons avoid introducing
 * a different interpolation filter at each candidate period. The low-pass
 * plus input signal emphasizes low frequencies while retaining upper-band
 * timing. Squared Kaiser weights match the energy of tapered observations. */
static float
refine_pitch(const struct mbe_analysis_tables* t, const struct pitch_work* work, float initial) {
    enum { REFINE_RADIUS = 3, REFINE_PROBES = 2 * REFINE_RADIUS + 1 };

    int coarse_lag = (int)lrintf(initial);
    /* Size the common support for the largest probe. Changing the support
     * between probes would compare different portions of evolving speech. */
    int half = (2 * MBE_ANALYSIS_PITCH_HALF - (coarse_lag + REFINE_RADIUS)) / 2;
    if (half > MBE_ANALYSIS_WR_HALF) {
        half = MBE_ANALYSIS_WR_HALF;
    }
    float weights[2 * MBE_ANALYSIS_WR_HALF + 1];
    const unsigned char* indices = window_indices(half);
    for (int n = -half; n <= half; ++n) {
        int wi = indices ? indices[n + half]
                         : (int)lrintf((float)n * (float)MBE_ANALYSIS_WR_HALF / (float)half) + MBE_ANALYSIS_WR_HALF;
        weights[n + half] = t->w_r[wi] * t->w_r[wi];
    }
    float score[REFINE_PROBES];
    for (int probe = 0; probe < REFINE_PROBES; ++probe) {
        int lag = coarse_lag + probe - REFINE_RADIUS;
        int earlier_center = MBE_ANALYSIS_PITCH_HALF - lag / 2;
        float error = 0.0f;
        float energy = 0.0f;
        for (int n = -half; n <= half; ++n) {
            float a = work->refinement_samples[earlier_center + lag + n];
            float b = work->refinement_samples[earlier_center + n];
            float w = weights[n + half];
            float d = a - b;
            error += w * d * d;
            energy += w * (a * a + b * b);
        }
        score[probe] = energy > 0.0f ? error / energy : 1.0f;
    }
    /* Search the five interior lags; the outer two supply their neighbors.
     * Start at the coarse lag so a tie with it needs no adjustment. A winner
     * at the search boundary need not bracket a true minimum. */
    int best = REFINE_RADIUS;
    for (int probe = 1; probe < REFINE_PROBES - 1; ++probe) {
        if (score[probe] < score[best]) {
            best = probe;
        }
    }
    float left = score[best - 1];
    float middle = score[best];
    float right = score[best + 1];
    float curvature = left + right - 2.0f * middle;
    float offset = curvature > 0.0f ? 0.5f * (left - right) / curvature : 0.0f;
    /* Keep the interpolated result within half a sample of the winning lag,
     * then quantize to eighths for the fractional-period phase correction. */
    offset = fmaxf(-0.5f, fminf(0.5f, offset));
    return (float)(coarse_lag + best - REFINE_RADIUS) + 0.125f * roundf(8.0f * offset);
}

/* Compare one harmonic's bins after correcting the fractional-period phase.
 * The earlier spectrum needs rotation by 2*pi*harmonic*f0*lag. Subtracting
 * the harmonic's whole cycles gives -2*pi*harmonic*f0*(period - lag).
 * Use the harmonic frequency, not each FFT bin's frequency: leakage from a
 * harmonic shares its phase. No input resampling is needed.
 *
 * Overlapping windows share noise samples. Their normalized correlation is
 * overlap, with a phase that varies by bin. Divide by the mean expected
 * noise disagreement, 1 - overlap*cos(angle), to keep white noise near unit
 * error even at short periods. Periodic input approaches zero; errors above
 * one are valid and are left for the existing voicing thresholds to reject.
 * The caller supplies an analysed harmonic with a nonempty bin interval. */
static float
harmonic_period_error(const float* later, const float* earlier, float f0, int harmonic, float period, int lag,
                      float overlap) {
    float phase = -2.0f * (float)M_PI * (float)harmonic * f0 * (period - (float)lag);
    float phase_cos = cosf(phase);
    float phase_sin = sinf(phase);
    int lo;
    int hi;
    band_edges(f0 * (float)ANALYSIS_FFT, harmonic, &lo, &hi);
    float energy = 0.0f;
    float error = 0.0f;
    float noise_fraction = 0.0f;
    for (int m = lo; m < hi; ++m) {
        float later_re, later_im, earlier_re, earlier_im;
        /* Real FFT packing: DC and Nyquist occupy slots 0 and 1; the
         * remaining bins are interleaved real/imaginary pairs. */
        if (m == 0 || m == ANALYSIS_BINS) {
            int k = m == 0 ? 0 : 1;
            later_re = later[k];
            later_im = 0.0f;
            earlier_re = earlier[k];
            earlier_im = 0.0f;
        } else {
            const size_t k = 2 * (size_t)m;
            later_re = later[k];
            later_im = later[k + 1];
            earlier_re = earlier[k];
            earlier_im = earlier[k + 1];
        }
        float dr = later_re - (phase_cos * earlier_re - phase_sin * earlier_im);
        float di = later_im - (phase_sin * earlier_re + phase_cos * earlier_im);
        error += dr * dr + di * di;
        energy += later_re * later_re + later_im * later_im + earlier_re * earlier_re + earlier_im * earlier_im;
        float angle = 2.0f * (float)M_PI * ((float)m / (float)ANALYSIS_FFT - (float)harmonic * f0) * (float)lag;
        noise_fraction += 1.0f - overlap * cosf(angle);
    }
    noise_fraction /= (float)(hi - lo);
    return energy > 0.0f && noise_fraction > 0.0f ? error / (energy * noise_fraction) : 1.0f;
}

/* Two tapered observations one rounded period apart distinguish repeating
 * harmonics from noise that fits a few correlated bins in a single spectrum.
 * Keep the longest window that fits the existing history for both observations;
 * this limits leakage without adding latency or per-stream state. The windows
 * use the same weights, and their normalized overlap supplies the expected
 * noise correlation used by harmonic_period_error(). */
static int
harmonic_periodicity(const struct mbe_analysis_tables* t, const struct mbe_analysis_state* s, mbe_fft_plan* fft,
                     struct mbe_analysis_result* r) {
    /* The existing coarse-pitch gate already rules out voicing in every band.
     * Leave the same unmeasured error used outside the analysed harmonics;
     * neither classifier needs a periodicity measurement once this gate closes. */
    if (r->pitch_error > VOICING_ERROR_MAX) {
        for (int l = 1; l <= r->harmonics; ++l) {
            r->fit_error[l] = 1.0f;
        }
        return 0;
    }
    float period = 1.0f / r->f0;
    int lag = (int)lrintf(period);
    int half = (MBE_ANALYSIS_BUFFER - 1 - lag) / 2;
    if (half > MBE_ANALYSIS_WR_HALF) {
        half = MBE_ANALYSIS_WR_HALF;
    }
    float window[2 * MBE_ANALYSIS_WR_HALF + 1];
    const unsigned char* indices = window_indices(half);
    float window_energy = 0.0f;
    float overlap = 0.0f;
    for (int i = 0; i <= 2 * half; ++i) {
        int wi =
            indices ? indices[i]
                    : (int)lrintf((float)(i - half) * (float)MBE_ANALYSIS_WR_HALF / (float)half) + MBE_ANALYSIS_WR_HALF;
        window[i] = t->w_r[wi];
        window_energy += window[i] * window[i];
        if (i >= lag) {
            overlap += window[i] * window[i - lag];
        }
    }
    overlap /= window_energy;
    /* Straddle the analysis centre. Integer indexing puts an odd-lag pair's
     * midpoint half a sample later; both observations keep the same taper. */
    int earlier_center = SAMPLE_CENTRE - lag / 2;
    float later[ANALYSIS_FFT] = {0};
    float earlier[ANALYSIS_FFT] = {0};
    float later_fft[ANALYSIS_FFT];
    float earlier_fft[ANALYSIS_FFT];
    for (int n = -half; n <= half; ++n) {
        float w = window[n + half];
        int m = (n + ANALYSIS_FFT) % ANALYSIS_FFT;
        later[m] = w * s->buf[earlier_center + lag + n];
        earlier[m] = w * s->buf[earlier_center + n];
    }
    int status = mbe_fft_forward_real(fft, later, later_fft);
    if (status < 0) {
        return status;
    }
    status = mbe_fft_forward_real(fft, earlier, earlier_fft);
    if (status < 0) {
        return status;
    }
    for (int l = 1; l <= r->harmonics; ++l) {
        r->fit_error[l] = harmonic_period_error(later_fft, earlier_fft, r->f0, l, period, lag, overlap);
    }
    return 0;
}

/* eqs (38)-(42), with energies on the 16-bit input scale. */
static float
energy_factor(const struct mbe_analysis_tables* t, struct mbe_analysis_state* s, const struct spectrum* sp) {
    const float norm = 1.0f / (t->w_r_sum * t->w_r_sum);
    float lf = 0.0f;
    float hf = 0.0f;
    for (int m = 0; m <= ANALYSIS_BINS; ++m) {
        float e = sp->power[m] * norm;
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

int
mbe_analysis_frame(const struct mbe_analysis_tables* tables, struct mbe_analysis_state* state, mbe_fft_plan* fft,
                   mbe_acf_plan* acf, struct mbe_analysis_result* result) {
    if (tables == NULL || state == NULL || fft == NULL || acf == NULL || result == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    struct pitch_work work;
    float energy;
    int acf_status = pitch_autocorrelation(tables, state, acf, &work, &energy);
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
    float pitch = refine_pitch(tables, &work, initial_pitch);

    memset(result, 0, sizeof(*result));
    result->f0 = 1.0f / pitch;
    result->pitch_error = work.error[initial];
    harmonic_measures(tables, &sp, result);
    status = harmonic_periodicity(tables, state, fft, result);
    if (status < 0) {
        return status;
    }
    result->energy_factor = energy_factor(tables, state, &sp);

    state->pitch_prev[1] = state->pitch_prev[0];
    state->pitch_prev[0] = initial_pitch;
    state->error_prev[1] = state->error_prev[0];
    state->error_prev[0] = work.error[initial];
    state->trusted = work.error[initial] <= PERIODIC_ERROR_MAX;
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
     * column (US 8,595,002). Retain the existing global gate: E(P_I) > 0.4
     * unvoices every column, rather than only bands above the first at 0.5
     * as in the standard. This is the coarse pitch error, distinct from the
     * per-harmonic period disagreement used below the gate. */
    if (result->pitch_error > VOICING_ERROR_MAX) {
        return 0.0f;
    }
    return voicing_threshold(state->columns_prev[column - 1], column, 0.1309f, result->energy_factor);
}

int
mbe_analysis_band_voicing(const struct mbe_analysis_result* result, int L, const unsigned char prev[MBE_ANALYSIS_BANDS],
                          unsigned char bands[MBE_ANALYSIS_BANDS]) {
    float error[MBE_ANALYSIS_BANDS + 1] = {0};
    float energy[MBE_ANALYSIS_BANDS + 1] = {0};
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
    }
    const float w0 = 2.0f * (float)M_PI * result->f0;
    for (int k = 1; k <= MBE_ANALYSIS_BANDS; k++) {
        float theta = 0.0f;
        /* The gate of mbe_analysis_column_threshold(); see the header. */
        if (k <= K && result->pitch_error <= VOICING_ERROR_MAX) {
            theta = voicing_threshold(prev[k - 1], k, w0, result->energy_factor);
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
