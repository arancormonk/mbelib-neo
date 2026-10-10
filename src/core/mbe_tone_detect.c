// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Tone detection for the AMBE+2 and D-STAR encoders' tone frames.
 *
 * TIA-102.BABA-1 7.1 leaves the method open. This one works on one 20 ms span:
 *
 *  - The tone must cover at least 45% of the span: 20-sample blocks within
 *    -10 dB of the strongest one mark the active part, which carries the
 *    rest of the test. DVSI's AMBE-3000 sends tone frames from about the
 *    same coverage at a tone's start and end. The part is trimmed to the
 *    samples the tone plays in, since a tone that starts or stops inside a
 *    block leaves silent samples there that a steady fit cannot explain.
 *  - A Hann-windowed spectrum must put most of its energy in two peaks, which
 *    rejects most speech cheaply and seeds the frequencies.
 *  - One or two sinusoids are fitted jointly by least squares (the projection
 *    onto their sine and cosine terms, so the result does not depend on phase
 *    and close pairs such as 440/480 Hz are fitted together), and their
 *    frequencies refined to maximize the explained energy.
 *  - A single tone is accepted when it explains nearly all of the active
 *    energy, at index round(f / 31.25 Hz); below 400 Hz, where a voiced
 *    frame's fundamental can dominate, nearly all means 99.7%. A pair is
 *    matched to the nearest DTMF or KNOX pair within 2% and 10 dB of twist
 *    (DVSI takes 1.5% and 8 dB, and about half the frames at 2% and 9 dB), or
 *    call-progress pair within 15% (DVSI sends 340 + 550 Hz as 350 + 490 Hz).
 *  - The level of a pair is the mean of its components' levels in dB, which
 *    is how DVSI's AD follows twisted DTMF.
 *  - A call-progress pair must also explain the whole span, not just its
 *    active part, so a span counts toward a call-progress tone only when the
 *    tone fills it (see mbe_tone_track()).
 *
 * The thresholds come from DVSI's rate-33 tone vectors, where the encoders
 * send DVSI's tone index on 98% of DVSI's tone frames, and its speech vectors,
 * where DVSI sends no tone frames and neither does this detector.
 */

#include "mbe_tone_detect.h"

#include <math.h>
#include <string.h>

#include "mbe_tone.h"

#define TONE_FS                   8000.0
#define TONE_BLOCK                20
#define TONE_BLOCKS               (MBE_TONE_SPAN / TONE_BLOCK)
#define TONE_MIN_ACTIVE           72    /* 45% of the span */
#define TONE_FLOOR_AMPLITUDE      40.0  /* peak, 16-bit scale */
#define TONE_ACTIVE_RATIO         0.1   /* block energy relative to the strongest block */
#define TONE_EDGE_RATIO           0.1   /* edge sample magnitude relative to the active part's largest */
#define TONE_PEAK_SHARE           0.6   /* energy within two bins of the two largest peaks */
#define TONE_PURITY_SINGLE_LOW    0.997 /* index below 13 (400 Hz) */
#define TONE_PURITY_SINGLE        0.98
#define TONE_PURITY_DUAL          0.93
#define TONE_PURITY_CALL_PROGRESS 0.98
#define TONE_DUAL_TOLERANCE       0.02   /* DTMF and KNOX, relative */
#define TONE_CALL_TOLERANCE       0.15   /* call progress, relative */
#define TONE_MAX_TWIST            3.1623 /* 10 dB between the components' amplitudes */
#define TONE_MIN_SEPARATION_HZ    15.0
#define TONE_MIN_HZ               60.0
#define TONE_MAX_HZ               3950.0
#define TONE_SINGLE_STEP_HZ       31.25
#define TONE_FIRST_SINGLE         5
#define TONE_LAST_SINGLE          122
#define TONE_FIRST_DUAL           128
#define TONE_LAST_DUAL            163
#define TONE_FIRST_CALL_PROGRESS  160

/* The part of the span the tone occupies. */
struct tone_span {
    const float* x;
    int n;
    double energy;
};

/* One or two sinusoids and how much of the span's energy they explain. */
struct tone_fit {
    int count;
    double hz[2];
    double amplitude[2];
    double explained;
};

/* Solve the symmetric positive definite system a x = b of order n <= 4 by
 * Cholesky factorization; returns 0 when a is not positive definite. */
static int
tone_solve(double a[4][4], const double b[4], int n, double x[4]) {
    double l[4][4] = {{0}};
    double y[4] = {0};
    if (n < 1 || n > 4) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = a[i][j];
            for (int k = 0; k < j; k++) {
                sum -= l[i][k] * l[j][k];
            }
            if (i == j) {
                if (sum <= 0.0) {
                    return 0;
                }
                l[i][i] = sqrt(sum);
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int k = 0; k < i; k++) {
            sum -= l[i][k] * y[k];
        }
        y[i] = sum / l[i][i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double sum = y[i];
        for (int k = i + 1; k < n; k++) {
            sum -= l[k][i] * x[k];
        }
        x[i] = sum / l[i][i];
    }
    return 1;
}

/* Least-squares fit of fit->count sinusoids at fit->hz: the energy of the
 * projection onto their cosine and sine terms, and each one's amplitude. */
static void
tone_project(const struct tone_span* s, struct tone_fit* fit) {
    double gram[4][4] = {{0}};
    double rhs[4] = {0};
    double coef[4] = {0};
    double c[2] = {1.0, 1.0}, sn[2] = {0.0, 0.0}, cw[2], sw[2];
    const int dim = 2 * fit->count;
    for (int k = 0; k < fit->count; k++) {
        cw[k] = cos(2.0 * M_PI * fit->hz[k] / TONE_FS);
        sw[k] = sin(2.0 * M_PI * fit->hz[k] / TONE_FS);
    }
    for (int i = 0; i < s->n; i++) {
        double v[4];
        for (int k = 0; k < fit->count; k++) {
            v[(size_t)2 * (size_t)k] = c[k];
            v[((size_t)2 * (size_t)k) + 1] = sn[k];
            double next = (c[k] * cw[k]) - (sn[k] * sw[k]);
            sn[k] = (sn[k] * cw[k]) + (c[k] * sw[k]);
            c[k] = next;
        }
        for (int r = 0; r < dim; r++) {
            rhs[r] += v[r] * (double)s->x[i];
            for (int q = 0; q <= r; q++) {
                gram[r][q] += v[r] * v[q];
            }
        }
    }
    for (int r = 0; r < dim; r++) {
        for (int q = r + 1; q < dim; q++) {
            gram[r][q] = gram[q][r];
        }
    }
    fit->explained = 0.0;
    if (!tone_solve(gram, rhs, dim, coef)) {
        return;
    }
    for (int r = 0; r < dim; r++) {
        fit->explained += rhs[r] * coef[r];
    }
    for (int k = 0; k < fit->count; k++) {
        fit->amplitude[k] = hypot(coef[(size_t)2 * (size_t)k], coef[((size_t)2 * (size_t)k) + 1]);
    }
}

/* Coordinate search on each frequency, halving the step, keeping any move
 * that explains more energy. */
static void
tone_refine(const struct tone_span* s, struct tone_fit* fit) {
    double step = 8.0;
    tone_project(s, fit);
    for (int iteration = 0; iteration < 4; iteration++, step *= 0.5) {
        for (int k = 0; k < 2 * fit->count; k++) {
            struct tone_fit trial = *fit;
            trial.hz[k / 2] += (k & 1) ? step : -step;
            if (trial.hz[k / 2] <= TONE_MIN_HZ || trial.hz[k / 2] >= TONE_MAX_HZ
                || (fit->count == 2 && fabs(trial.hz[0] - trial.hz[1]) < TONE_MIN_SEPARATION_HZ)) {
                continue;
            }
            tone_project(s, &trial);
            if (trial.explained > fit->explained) {
                *fit = trial;
            }
        }
    }
}

/* Samples [*lo, *hi) of active blocks first..last, less the samples of the
 * first and last block before the tone starts or after it ends: those below
 * TONE_EDGE_RATIO of the largest. A steady fit cannot explain them, so a tone
 * that starts or stops inside a block would otherwise fail its purity test. */
static void
tone_trim_edges(const float span[MBE_TONE_SPAN], int first, int last, int* lo, int* hi) {
    double top = 0.0;
    *lo = first * TONE_BLOCK;
    *hi = (last + 1) * TONE_BLOCK;
    for (int i = *lo; i < *hi; i++) {
        top = fmax(top, fabs((double)span[i]));
    }
    const double edge = TONE_EDGE_RATIO * top;
    while (*lo < (first + 1) * TONE_BLOCK && fabs((double)span[*lo]) < edge) {
        (*lo)++;
    }
    while (*hi > last * TONE_BLOCK && fabs((double)span[*hi - 1]) < edge) {
        (*hi)--;
    }
}

/* The active part of the span: from the first to the last 20-sample block
 * within TONE_ACTIVE_RATIO of the strongest, trimmed to the samples the tone
 * plays in. Returns 0 if it is too short. */
static int
tone_active_span(const float span[MBE_TONE_SPAN], struct tone_span* s) {
    double block[TONE_BLOCKS] = {0};
    double peak = 0.0;
    for (int b = 0; b < TONE_BLOCKS; b++) {
        for (int i = 0; i < TONE_BLOCK; i++) {
            double v = span[(b * TONE_BLOCK) + i];
            block[b] += v * v;
        }
        peak = fmax(peak, block[b]);
    }
    /* The strongest block always qualifies, so first <= last. */
    int first = 0;
    while (first < TONE_BLOCKS - 1 && block[first] < TONE_ACTIVE_RATIO * peak) {
        first++;
    }
    int last = TONE_BLOCKS - 1;
    while (last > first && block[last] < TONE_ACTIVE_RATIO * peak) {
        last--;
    }
    int lo = 0;
    int hi = 0;
    tone_trim_edges(span, first, last, &lo, &hi);
    s->x = span + lo;
    s->n = hi - lo;
    s->energy = 0.0;
    for (int i = lo; i < hi; i++) {
        s->energy += (double)span[i] * (double)span[i];
    }
    return s->n >= TONE_MIN_ACTIVE && s->energy > 0.0;
}

#define TONE_BINS ((MBE_FFT_SIZE / 2) + 1)

/* Power spectrum of the Hann-windowed span; returns its total, or a negative
 * MBE_STATUS_* value from the FFT. */
static double
tone_power_spectrum(mbe_fft_plan* fft, const float span[MBE_TONE_SPAN], double power[TONE_BINS]) {
    float in[MBE_FFT_SIZE] = {0};
    float packed[MBE_FFT_SIZE];
    for (int i = 0; i < MBE_TONE_SPAN; i++) {
        double w = 0.5 - (0.5 * cos(2.0 * M_PI * (double)(i + 1) / (double)(MBE_TONE_SPAN + 1)));
        in[i] = (float)(w * (double)span[i]);
    }
    int status = mbe_fft_forward_real(fft, in, packed);
    if (status < 0) {
        return (double)status;
    }
    /* Packed real FFT: DC and Nyquist first, then (re, im) of bins 1..127. */
    double total = 0.0;
    power[0] = (double)packed[0] * (double)packed[0];
    power[TONE_BINS - 1] = (double)packed[1] * (double)packed[1];
    for (int m = 1; m < TONE_BINS - 1; m++) {
        const size_t re = (size_t)2 * (size_t)m;
        power[m] = ((double)packed[re] * (double)packed[re]) + ((double)packed[re + 1] * (double)packed[re + 1]);
    }
    for (int m = 0; m < TONE_BINS; m++) {
        total += power[m];
    }
    return total;
}

/* The two largest local maxima of the spectrum (0 where there is none). */
static void
tone_top_bins(const double power[TONE_BINS], int top[2]) {
    top[0] = 0;
    top[1] = 0;
    for (int m = 1; m < TONE_BINS - 1; m++) {
        if (power[m] < power[m - 1] || power[m] < power[m + 1]) {
            continue;
        }
        if (top[0] == 0 || power[m] > power[top[0]]) {
            top[1] = top[0];
            top[0] = m;
        } else if (top[1] == 0 || power[m] > power[top[1]]) {
            top[1] = m;
        }
    }
}

/* The frequency of the peak at bin m, by parabolic interpolation. */
static double
tone_peak_hz(const double power[TONE_BINS], int m) {
    double den = power[m - 1] - (2.0 * power[m]) + power[m + 1];
    double offset = (den != 0.0) ? 0.5 * (power[m - 1] - power[m + 1]) / den : 0.0;
    return ((double)m + offset) * TONE_FS / (double)MBE_FFT_SIZE;
}

/* The two largest peaks of the Hann-windowed power spectrum, interpolated,
 * and the share of the energy within two bins of them. Returns how many
 * peaks there are (0..2) or a negative MBE_STATUS_* value. */
static int
tone_peaks(mbe_fft_plan* fft, const float span[MBE_TONE_SPAN], double hz[2], double* share) {
    double power[TONE_BINS];
    double total = tone_power_spectrum(fft, span, power);
    if (total < 0.0) {
        return (int)total;
    }
    int top[2];
    unsigned char near_peak[TONE_BINS] = {0};
    double peak_energy = 0.0;
    int found = 0;
    tone_top_bins(power, top);
    for (int k = 0; k < 2 && top[k] > 0; k++, found++) {
        hz[k] = tone_peak_hz(power, top[k]);
        for (int j = top[k] - 2; j <= top[k] + 2; j++) {
            if (j >= 0 && j < TONE_BINS && !near_peak[j]) {
                near_peak[j] = 1;
                peak_energy += power[j];
            }
        }
    }
    *share = (total > 0.0) ? peak_energy / total : 0.0;
    return found;
}

/* Frequencies of the two sinusoids that a 4th-order symmetric linear
 * predictor fits: x(n) + x(n-4) = s (x(n-1) + x(n-3)) - p x(n-2), whose
 * characteristic roots are 2 cos(w1) and 2 cos(w2). Resolves pairs closer
 * than the spectrum can. Returns 0 when the fit has no such roots. */
static int
tone_pair_from_predictor(const struct tone_span* s, double hz[2]) {
    double vv = 0.0, vw = 0.0, ww = 0.0, vu = 0.0, wu = 0.0;
    for (int i = 4; i < s->n; i++) {
        double u = (double)s->x[i] + (double)s->x[i - 4];
        double v = (double)s->x[i - 1] + (double)s->x[i - 3];
        double w = s->x[i - 2];
        vv += v * v;
        vw += v * w;
        ww += w * w;
        vu += v * u;
        wu += w * u;
    }
    double det = (vv * ww) - (vw * vw);
    if (det <= 0.0) {
        return 0;
    }
    double sum = ((vu * ww) - (vw * wu)) / det;  /* s */
    double prod = ((vu * vw) - (vv * wu)) / det; /* p */
    double disc = (sum * sum) - (4.0 * (prod - 2.0));
    if (disc < 0.0) {
        return 0;
    }
    const double r[2] = {0.5 * (sum + sqrt(disc)), 0.5 * (sum - sqrt(disc))};
    for (int k = 0; k < 2; k++) {
        if (fabs(r[k]) > 2.0) {
            return 0;
        }
        hz[k] = acos(0.5 * r[k]) * TONE_FS / (2.0 * M_PI);
    }
    return 1;
}

/* The single-tone index of a fitted sinusoid, or -1. */
static int
tone_single_id(const struct tone_span* s, const struct tone_fit* fit) {
    long id = lround(fit->hz[0] / TONE_SINGLE_STEP_HZ);
    if (id < TONE_FIRST_SINGLE || id > TONE_LAST_SINGLE) {
        return -1;
    }
    double purity = (id < 13) ? TONE_PURITY_SINGLE_LOW : TONE_PURITY_SINGLE;
    return (fit->explained >= purity * s->energy) ? (int)id : -1;
}

/* The dual-tone index nearest a fitted pair, or -1. */
static int
tone_dual_id(const struct tone_span* s, const struct tone_fit* fit) {
    double lo = fmin(fit->hz[0], fit->hz[1]);
    double hi = fmax(fit->hz[0], fit->hz[1]);
    double a_lo = fmin(fit->amplitude[0], fit->amplitude[1]);
    double a_hi = fmax(fit->amplitude[0], fit->amplitude[1]);
    if (a_hi > TONE_MAX_TWIST * a_lo) {
        return -1;
    }
    int best = -1;
    double best_error = 0.0;
    for (int id = TONE_FIRST_DUAL; id <= TONE_LAST_DUAL; id++) {
        float f_hi, f_lo;
        (void)mbe_tone_lookup_freqs(id, &f_hi, &f_lo);
        double tolerance = (id >= TONE_FIRST_CALL_PROGRESS) ? TONE_CALL_TOLERANCE : TONE_DUAL_TOLERANCE;
        double e_lo = fabs(lo - f_lo) / f_lo;
        double e_hi = fabs(hi - f_hi) / f_hi;
        if (e_lo <= tolerance && e_hi <= tolerance && (best < 0 || e_lo + e_hi < best_error)) {
            best = id;
            best_error = e_lo + e_hi;
        }
    }
    if (best < 0) {
        return -1;
    }
    double purity = (best >= TONE_FIRST_CALL_PROGRESS) ? TONE_PURITY_CALL_PROGRESS : TONE_PURITY_DUAL;
    return (fit->explained >= purity * s->energy) ? best : -1;
}

/* The better of the pair seeded from the spectrum and the one from the predictor. */
static void
tone_fit_pair(const struct tone_span* s, const double peaks[2], int peak_count, struct tone_fit* best) {
    double seeds[2][2] = {{peaks[0], peaks[1]}, {0.0, 0.0}};
    int first = (peak_count == 2) ? 0 : 1;
    int count = tone_pair_from_predictor(s, seeds[1]) ? 2 : 1;
    best->explained = -1.0;
    for (int k = first; k < count; k++) {
        struct tone_fit fit = {2, {seeds[k][0], seeds[k][1]}, {0.0, 0.0}, 0.0};
        if (fabs(fit.hz[0] - fit.hz[1]) < TONE_MIN_SEPARATION_HZ) {
            continue;
        }
        tone_refine(s, &fit);
        if (fit.explained > best->explained) {
            *best = fit;
        }
    }
}

/* mbe_tone_detect() for a span whose signal is `samples` long, the rest zero:
 * the level floor applies to those samples. */
static int
tone_detect_part(mbe_fft_plan* fft, const float span[MBE_TONE_SPAN], int samples, struct mbe_tone_detection* out) {
    struct tone_span s;
    double peaks[2] = {0.0, 0.0};
    double share = 0.0;
    double energy = 0.0;
    for (int i = 0; i < MBE_TONE_SPAN; i++) {
        energy += (double)span[i] * (double)span[i];
    }
    if (energy < samples * TONE_FLOOR_AMPLITUDE * TONE_FLOOR_AMPLITUDE / 2.0 || !tone_active_span(span, &s)) {
        return 0;
    }
    int peak_count = tone_peaks(fft, span, peaks, &share);
    if (peak_count <= 0 || share < TONE_PEAK_SHARE) {
        return (peak_count < 0) ? peak_count : 0;
    }
    struct tone_fit single = {1, {peaks[0], 0.0}, {0.0, 0.0}, 0.0};
    tone_refine(&s, &single);
    int id = tone_single_id(&s, &single);
    if (id >= 0) {
        out->id = id;
        out->amplitude = (float)single.amplitude[0];
        return 1;
    }
    struct tone_fit pair = {2, {0.0, 0.0}, {0.0, 0.0}, -1.0};
    tone_fit_pair(&s, peaks, peak_count, &pair);
    id = (pair.explained > 0.0) ? tone_dual_id(&s, &pair) : -1;
    if (id < 0) {
        return 0;
    }
    if (id >= TONE_FIRST_CALL_PROGRESS) {
        struct tone_span whole = {span, MBE_TONE_SPAN, energy};
        struct tone_fit fit = pair;
        tone_project(&whole, &fit);
        if (!(fit.explained >= TONE_PURITY_CALL_PROGRESS * energy)) {
            return 0;
        }
    }
    out->id = id;
    out->amplitude = (float)sqrt(pair.amplitude[0] * pair.amplitude[1]);
    return 1;
}

int
mbe_tone_detect(mbe_fft_plan* fft, const float span[MBE_TONE_SPAN], struct mbe_tone_detection* out) {
    return tone_detect_part(fft, span, MBE_TONE_SPAN, out);
}

/* A tone in the newest or the oldest half of span alone, as mbe_tone_detect(). */
static int
tone_detect_half(mbe_fft_plan* fft, const float span[MBE_TONE_SPAN], int newest, struct mbe_tone_detection* out) {
    float half[MBE_TONE_SPAN] = {0.0f};
    const size_t start = newest ? (size_t)(MBE_TONE_SPAN / 2) : 0;
    memcpy(half + start, span + start, (size_t)(MBE_TONE_SPAN / 2) * sizeof(float));
    return tone_detect_part(fft, half, MBE_TONE_SPAN / 2, out);
}

/* The tone in span after a frame with DTMF, KNOX or single tone previous (-1
 * for none), as mbe_tone_detect(). Where one tone changes directly to another
 * the span holds both and fits neither, and DVSI's encoders send the newer one
 * if it fills about half of the span, else the older. So a span without a tone
 * carries a DTMF or KNOX tone that fills its newest half, or else the previous
 * tone if that fills the oldest half. Not a single tone in the newest half:
 * 80 samples of voiced speech can pass for one. */
static int
tone_detect_after(mbe_fft_plan* fft, const float span[MBE_TONE_SPAN], int previous, struct mbe_tone_detection* out) {
    int status = mbe_tone_detect(fft, span, out);
    if (status != 0 || previous < 0) {
        return status;
    }
    status = tone_detect_half(fft, span, 1, out);
    if (status < 0 || (status > 0 && out->id >= TONE_FIRST_DUAL && out->id < TONE_FIRST_CALL_PROGRESS)) {
        return status;
    }
    status = tone_detect_half(fft, span, 0, out);
    return (status > 0 && out->id != previous) ? 0 : status;
}

void
mbe_tone_tracker_reset(struct mbe_tone_tracker* tracker) {
    tracker->run_id = -1;
    tracker->run = 0;
    tracker->hold = 0;
    tracker->sending = 0;
    tracker->last.id = -1;
    tracker->last.amplitude = 0.0f;
    tracker->tone_id = -1;
}

int
mbe_tone_track(struct mbe_tone_tracker* tracker, mbe_fft_plan* fft, const float span[MBE_TONE_SPAN],
               struct mbe_tone_detection* out) {
    struct mbe_tone_detection tone;
    int status = tone_detect_after(fft, span, tracker->tone_id, &tone);
    if (status < 0) {
        return status;
    }
    if (status > 0 && tone.id < TONE_FIRST_CALL_PROGRESS) {
        mbe_tone_tracker_reset(tracker);
        tracker->tone_id = tone.id;
        *out = tone;
        return 1;
    }
    tracker->tone_id = -1;
    if (status == 0) {
        tracker->run_id = -1;
        tracker->run = 0;
    } else if (tone.id != tracker->run_id) {
        tracker->run_id = tone.id;
        tracker->run = 1;
    } else if (tracker->run < MBE_TONE_CP_CONFIRM) {
        tracker->run++;
    }
    if (status > 0 && ((tracker->sending && tone.id == tracker->last.id) || tracker->run >= MBE_TONE_CP_CONFIRM)) {
        tracker->sending = 1;
        tracker->hold = MBE_TONE_CP_HOLD;
        tracker->last = tone;
        *out = tone;
        return 1;
    }
    if (tracker->sending && tracker->hold > 0) {
        tracker->hold--;
        *out = tracker->last;
        return 1;
    }
    tracker->sending = 0;
    return 0;
}
