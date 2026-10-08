// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 *
 * Based on Bruce Perens' ham_digital_modes
 * (https://github.com/BrucePerens/hams_open).
 */

/**
 * @file
 * @brief Encoder-side DTMF and single-tone detection.
 *
 * Ported from ham_digital_modes' float tone_detect.rs, whose behaviour was
 * calibrated against the DVSI AMBE-3000 chip: DTMF digits are detected in
 * their first frame, and single tones from 400 Hz to 3800 Hz and at 200 Hz,
 * with index round(f / 31.25 Hz).
 */

#include "mbe_tone_detect.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double mbe_dtmf_row_hz[4] = {697.0, 770.0, 852.0, 941.0};
static const double mbe_dtmf_col_hz[4] = {1209.0, 1336.0, 1477.0, 1633.0};

/* Fraction of a frame's energy one sinusoid must explain for the frame to be
 * classified as a single tone. A threshold of 0.9 classified voiced speech
 * dominated by its fundamental as a tone. */
#define MBE_TONE_SINGLE_MIN_EXPLAINED 0.99
/* Explained-energy fraction required of a 200 Hz tone, the only tone below
 * 400 Hz that is reported. */
#define MBE_TONE_LOW_MIN_EXPLAINED    0.999

/* Amplitude of the least-squares sinusoid of frequency hz in x. */
static double
mbe_tone_fit_sinusoid(const double* x, double hz) {
    double w = 2.0 * M_PI * hz / 8000.0;
    double sc = 0.0;
    double cc = 0.0;
    double sn = 0.0;
    double xc = 0.0;
    double xs = 0.0;
    for (int i = 0; i < MBE_TONE_DETECT_FRAME; i++) {
        double s = sin(w * (double)i);
        double c = cos(w * (double)i);
        cc += c * c;
        sn += s * s;
        sc += s * c;
        xc += x[i] * c;
        xs += x[i] * s;
    }
    double det = (cc * sn) - (sc * sc);
    if (fabs(det) < 1e-9) {
        return 0.0;
    }
    double a = ((xc * sn) - (xs * sc)) / det;
    double b = ((xs * cc) - (xc * sc)) / det;
    return sqrt((a * a) + (b * b));
}

/* Fraction of the frame's energy explained by independently fitted sinusoids. */
static double
mbe_tone_explained_fraction(const double* x, const double* hzs, int count) {
    double total = 0.0;
    double explained = 0.0;
    for (int i = 0; i < MBE_TONE_DETECT_FRAME; i++) {
        total += x[i] * x[i];
    }
    for (int i = 0; i < count; i++) {
        double a = mbe_tone_fit_sinusoid(x, hzs[i]);
        explained += a * a / 2.0 * (double)MBE_TONE_DETECT_FRAME;
    }
    if (total > 0.0) {
        return fmin(explained / total, 1.5);
    }
    return 0.0;
}

/* Index of the largest value; ties resolve to the highest index. */
static int
mbe_tone_arg_max(const double v[4]) {
    int best = 0;
    for (int i = 1; i < 4; i++) {
        if (v[i] >= v[best]) {
            best = i;
        }
    }
    return best;
}

/* DTMF: the strongest row and column, together explaining the frame. */
static int
mbe_tone_detect_dtmf(const double* frame, struct mbe_tone_detection* out) {
    double row_amps[4];
    double col_amps[4];
    for (int i = 0; i < 4; i++) {
        row_amps[i] = mbe_tone_fit_sinusoid(frame, mbe_dtmf_row_hz[i]);
        col_amps[i] = mbe_tone_fit_sinusoid(frame, mbe_dtmf_col_hz[i]);
    }
    int r = mbe_tone_arg_max(row_amps);
    int c = mbe_tone_arg_max(col_amps);
    double ra = row_amps[r];
    double ca = col_amps[c];
    if (!(ra > 100.0 && ca > 100.0 && fmax(ra / ca, ca / ra) < 3.0)) {
        return 0;
    }
    const double hzs[2] = {mbe_dtmf_row_hz[r], mbe_dtmf_col_hz[c]};
    if (!(mbe_tone_explained_fraction(frame, hzs, 2) > 0.85)) {
        return 0;
    }
    out->kind = MBE_TONE_DETECT_DTMF;
    out->row = r;
    out->col = c;
    out->amplitude = (ra + ca) / 2.0;
    return 1;
}

/* Single tone: coarse scan for the strongest sinusoid, then refine. */
static int
mbe_tone_detect_single(const double* frame, struct mbe_tone_detection* out) {
    double best_amp = 0.0;
    double best_hz = 0.0;
    /* 150 to 3900 Hz in 12.5 Hz steps, then +-12.5 Hz in 1 Hz steps. Both
     * steps are exactly representable in binary floating point, so these
     * values equal those of the reference's accumulated loop variables. */
    for (int step = 0; step <= 300; step++) {
        double hz = 150.0 + (12.5 * (double)step);
        double a = mbe_tone_fit_sinusoid(frame, hz);
        if (a > best_amp) {
            best_amp = a;
            best_hz = hz;
        }
    }
    double amp = best_amp;
    double hz = best_hz;
    for (int step = 0; step <= 25; step++) {
        double f = best_hz - 12.5 + (double)step;
        double a = mbe_tone_fit_sinusoid(frame, f);
        if (a > amp) {
            amp = a;
            hz = f;
        }
    }
    double frac = mbe_tone_explained_fraction(frame, &hz, 1);
    if (!(amp > 100.0 && frac > MBE_TONE_SINGLE_MIN_EXPLAINED)) {
        return 0;
    }
    int index = (int)round(hz / 31.25);
    /* Below 400 Hz, voiced speech is itself close to a pure sinusoid, so only a
     * component at 200 Hz that explains nearly all of the frame's energy is
     * accepted there. */
    int in_range =
        (index >= 13 && index <= 122) || (index == 6 && fabs(hz - 200.0) < 3.0 && frac > MBE_TONE_LOW_MIN_EXPLAINED);
    if (!in_range) {
        return 0;
    }
    out->kind = MBE_TONE_DETECT_SINGLE;
    out->index = index;
    out->hz = hz;
    out->amplitude = amp;
    return 1;
}

enum mbe_tone_detect_kind
mbe_tone_detect(const double frame[MBE_TONE_DETECT_FRAME], struct mbe_tone_detection* out) {
    double energy = 0.0;
    memset(out, 0, sizeof(*out));

    for (int i = 0; i < MBE_TONE_DETECT_FRAME; i++) {
        energy += frame[i] * frame[i];
    }
    energy /= (double)MBE_TONE_DETECT_FRAME;
    if (energy < 100.0 * 100.0 / 2.0) {
        return MBE_TONE_DETECT_NONE; /* below the mean power of a sinusoid of amplitude 100 */
    }
    if (mbe_tone_detect_dtmf(frame, out) || mbe_tone_detect_single(frame, out)) {
        return out->kind;
    }
    return MBE_TONE_DETECT_NONE;
}
