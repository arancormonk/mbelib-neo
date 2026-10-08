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
 * @brief TIA-102.BABA speech analysis shared by the IMBE and AMBE+2 encoders.
 *
 * Ported from ham_digital_modes' float tia_102_baba modules: encoder.rs
 * (FrameAnalyzer, HighPassFilter), pitch.rs, pitch_refinement.rs, vuv.rs and
 * spectral_amplitude.rs. Equation numbers refer to TIA-102.BABA (2003).
 */

#include "mbe_frame_analysis.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Annex B: the initial pitch window w_I(n), n = 0..150 (w_I(-n) = w_I(n)). */
static const double mbe_fa_wi_half[151] = {
    0.09207659, 0.09206691, 0.09203795, 0.09198972, 0.09192220, 0.09183547, 0.09172948, 0.09160442, 0.09146029,
    0.09129713, 0.09111510, 0.09091420, 0.09069464, 0.09045644, 0.09019975, 0.08992471, 0.08963151, 0.08932017,
    0.08899096, 0.08864402, 0.08827950, 0.08789764, 0.08749853, 0.08708246, 0.08664963, 0.08620022, 0.08573447,
    0.08525267, 0.08475497, 0.08424167, 0.08371302, 0.08316927, 0.08261070, 0.08203757, 0.08145022, 0.08084887,
    0.08023385, 0.07960547, 0.07896401, 0.07830978, 0.07764313, 0.07696434, 0.07627381, 0.07557184, 0.07485873,
    0.07413481, 0.07340052, 0.07265611, 0.07190202, 0.07113854, 0.07036604, 0.06958490, 0.06879549, 0.06799815,
    0.06719328, 0.06638119, 0.06556234, 0.06473708, 0.06390573, 0.06306869, 0.06222639, 0.06137912, 0.06052732,
    0.05967136, 0.05881159, 0.05794842, 0.05708219, 0.05621329, 0.05534209, 0.05446897, 0.05359430, 0.05271840,
    0.05184171, 0.05096456, 0.05008731, 0.04921030, 0.04833391, 0.04745849, 0.04658438, 0.04571191, 0.04484143,
    0.04397328, 0.04310780, 0.04224531, 0.04138612, 0.04053057, 0.03967896, 0.03883159, 0.03798876, 0.03715080,
    0.03631795, 0.03549054, 0.03466884, 0.03385311, 0.03304362, 0.03224064, 0.03144442, 0.03065519, 0.02987323,
    0.02909875, 0.02833197, 0.02757312, 0.02682242, 0.02608007, 0.02534628, 0.02462122, 0.02390509, 0.02319806,
    0.02250030, 0.02181198, 0.02113325, 0.02046427, 0.01980515, 0.01915605, 0.01851709, 0.01788837, 0.01727001,
    0.01666212, 0.01606477, 0.01547807, 0.01490209, 0.01433691, 0.01378257, 0.01323915, 0.01270669, 0.01218523,
    0.01167480, 0.01117544, 0.01068715, 0.01020996, 0.00974387, 0.00928887, 0.00884496, 0.00841213, 0.00799034,
    0.00757957, 0.00717979, 0.00679095, 0.00641300, 0.00604590, 0.00568957, 0.00534396, 0.00500898, 0.00468457,
    0.00437064, 0.00406710, 0.00377385, 0.00349080, 0.00321783, 0.00295485, 0.00270174,
};

/* Annex D: the lowpass filter h_LPF(n), n = 0..10. */
static const double mbe_fa_lpf_half[11] = {
    0.351338, 0.278990, 0.118754, -0.015116, -0.055990, -0.026955, 0.008800, 0.016601, 0.005666, -0.002831, -0.002898,
};

/* Annex C: the refinement window w_R(n), n = 0..110. */
static const double mbe_fa_wr_half[111] = {
    1.000000, 0.999774, 0.999095, 0.997966, 0.996386, 0.994358, 0.991884, 0.988967, 0.985610, 0.981817, 0.977592,
    0.972940, 0.967866, 0.962377, 0.956477, 0.950174, 0.943474, 0.936386, 0.928916, 0.921074, 0.912868, 0.904307,
    0.895400, 0.886157, 0.876589, 0.866705, 0.856516, 0.846033, 0.835267, 0.824231, 0.812935, 0.801391, 0.789612,
    0.777610, 0.765397, 0.752986, 0.740390, 0.727620, 0.714692, 0.701616, 0.688406, 0.675076, 0.661638, 0.648105,
    0.634490, 0.620807, 0.607067, 0.593284, 0.579470, 0.565639, 0.551802, 0.537971, 0.524160, 0.510379, 0.496640,
    0.482955, 0.469336, 0.455793, 0.442337, 0.428978, 0.415727, 0.402594, 0.389588, 0.376718, 0.363994, 0.351425,
    0.339018, 0.326782, 0.314724, 0.302851, 0.291171, 0.279689, 0.268413, 0.257347, 0.246497, 0.235869, 0.225466,
    0.215294, 0.205355, 0.195653, 0.186192, 0.176974, 0.168001, 0.159276, 0.150799, 0.142572, 0.134596, 0.126872,
    0.119398, 0.112176, 0.105205, 0.098483, 0.092009, 0.085782, 0.079801, 0.074062, 0.068563, 0.063303, 0.058277,
    0.053482, 0.048915, 0.044573, 0.040451, 0.036546, 0.032852, 0.029365, 0.026081, 0.022995, 0.020102, 0.017397,
    0.014873,
};

static double
mbe_fa_wi(int n) {
    return mbe_fa_wi_half[(n < 0) ? -n : n];
}

static double
mbe_fa_lpf(int n) {
    return mbe_fa_lpf_half[(n < 0) ? -n : n];
}

static double
mbe_fa_wr(int n) {
    return mbe_fa_wr_half[(n < 0) ? -n : n];
}

static double
mbe_fa_candidate(int i) {
    return 21.0 + (0.5 * (double)i);
}

static double
mbe_fa_wr_dft_direct(int m) {
    double acc = 0.0;
    for (int n = -110; n <= 110; n++) {
        double theta = -2.0 * M_PI * (double)m * (double)n / 16384.0;
        acc += mbe_fa_wr(n) * cos(theta);
    }
    return acc;
}

/* W_R(m) (eq 30) is even in m; it is tabulated for 0 <= m <= 8192, which
 * covers every index the analysis evaluates. */
static double
mbe_fa_wr_dft(const struct mbe_fa_tables* tables, int m) {
    int am = (m < 0) ? -m : m;
    if (am <= MBE_FA_WR_DFT_HALF) {
        return tables->wr_dft[am];
    }
    return mbe_fa_wr_dft_direct(m);
}

void
mbe_fa_init_tables(struct mbe_fa_tables* tables) {
    for (int m = 0; m <= MBE_FA_WR_DFT_HALF; m++) {
        tables->wr_dft[m] = mbe_fa_wr_dft_direct(m);
    }
    for (int k = 0; k < MBE_FA_BINS; k++) {
        double theta = -2.0 * M_PI * (double)k / 256.0;
        tables->twiddle_re[k] = cos(theta);
        tables->twiddle_im[k] = sin(theta);
    }
    for (int i = 0; i < MBE_FA_CANDIDATES; i++) {
        double p = mbe_fa_candidate(i);
        double lo = 0.8 * p;
        double hi = 1.2 * p;
        int first = -1;
        int last = -1;
        for (int j = 0; j < MBE_FA_CANDIDATES; j++) {
            double q = mbe_fa_candidate(j);
            if (q >= lo && q <= hi) {
                if (first < 0) {
                    first = j;
                }
                last = j;
            }
        }
        tables->range_lo[i] = first;
        tables->range_hi[i] = last;
    }
    tables->wr_tap_sum = 0.0;
    for (int n = -110; n <= 110; n++) {
        tables->wr_tap_sum += mbe_fa_wr(n);
    }
}

void
mbe_fa_reset(struct mbe_fa_state* state) {
    memset(state, 0, sizeof(*state));
    /* "Upon initialization ... P_hat_-1 and P_hat_-2 are assumed to be equal to 100." */
    state->prev1_pitch = 100.0;
    state->prev2_pitch = 100.0;
}

void
mbe_fa_voicing_reset(struct mbe_fa_voicing_state* voicing) {
    memset(voicing, 0, sizeof(*voicing));
    voicing->xi_max = 20000.0;
}

/* Floor of v / 2^shift. Right-shifting a negative signed integer is
 * implementation-defined in C, so negative values are handled explicitly. */
static int64_t
mbe_fa_floor_shift(int64_t v, int shift) {
    if (v >= 0) {
        return v >> shift;
    }
    return -((-(v + 1)) >> shift) - 1;
}

/* 0.99 in Q30, round(0.99 * 2^30). */
#define MBE_FA_HP_POLE_Q30 ((int64_t)1063004406)

/*
 * Eq 3: H(z) = (1 - z^-1) / (1 - 0.99 z^-1), computed in the reference's
 * integer arithmetic, which its floating-point and fixed-point encoders share.
 * The feedback product is split at 2^30 so that it cannot overflow 64 bits.
 */
static double
mbe_fa_high_pass(struct mbe_fa_state* state, double sample) {
    int64_t x = (int64_t)round(sample);
    int64_t q = state->hp_prev_output_q16;
    int64_t q_hi = mbe_fa_floor_shift(q, 30);
    int64_t q_lo = q - (q_hi * ((int64_t)1 << 30));
    int64_t feedback =
        (q_hi * MBE_FA_HP_POLE_Q30) + mbe_fa_floor_shift((q_lo * MBE_FA_HP_POLE_Q30) + ((int64_t)1 << 29), 30);
    int64_t y_q16 = ((x - state->hp_prev_input) * ((int64_t)1 << 16)) + feedback;
    state->hp_prev_input = x;
    state->hp_prev_output_q16 = y_q16;
    return (double)mbe_fa_floor_shift(y_q16 + ((int64_t)1 << 15), 16);
}

/* E(P) (eq 5-8) at every candidate for the frame centred at buf[center]. */
static void
mbe_fa_error_table(const double* buf, int center, double error[MBE_FA_CANDIDATES]) {
    double s_lpf[301];
    double a[301];
    double r_table[MBE_FA_R_LAGS + 1];
    double energy = 0.0;
    double w4_sum = 0.0;

    /* Eq 9: s_LPF(n) = sum s(n - j) h_LPF(j). */
    for (int i = 0; i < 301; i++) {
        int j = i - 150;
        double acc = 0.0;
        for (int t = -10; t <= 10; t++) {
            acc += buf[center + j - t] * mbe_fa_lpf(t);
        }
        s_lpf[i] = acc;
    }
    for (int i = 0; i < 301; i++) {
        double w = mbe_fa_wi(i - 150);
        double v = s_lpf[i] * w;
        double w2 = w * w;
        energy += v * v;
        w4_sum += w2 * w2;
        a[i] = s_lpf[i] * w * w;
    }
    /* Eq 7, r(t) at every integer lag; terms with j + t outside the window vanish. */
    for (int t = 0; t <= MBE_FA_R_LAGS; t++) {
        double sum = 0.0;
        for (int i = 0; i < 301 - t; i++) {
            sum += a[i] * a[i + t];
        }
        r_table[t] = sum;
    }

    for (int c = 0; c < MBE_FA_CANDIDATES; c++) {
        double p = mbe_fa_candidate(c);
        int n_max = (int)floor(150.0 / p);
        double sum = 0.0;
        for (int n = 1; n <= n_max; n++) {
            /* Eq 8: linear interpolation of r(t). */
            double t = (double)n * p;
            double t_floor = floor(t);
            int lag = (int)t_floor;
            int lag_up = (lag + 1 < MBE_FA_R_LAGS) ? lag + 1 : MBE_FA_R_LAGS;
            sum += ((1.0 + t_floor - t) * r_table[lag]) + ((t - t_floor) * r_table[lag_up]);
        }
        double r_sum = r_table[0] + (2.0 * sum);
        double denominator = energy * (1.0 - (p * w4_sum));
        if (fabs(denominator) < 1e-12) {
            error[c] = 1.0; /* silent input: no pitch evidence */
        } else {
            /* E is in theory a ratio in [0, 1]; the interpolation of eq 8 can
             * make it slightly negative for a periodic signal, which would
             * leave the look-ahead ratio tests ill-defined. */
            error[c] = fmax((energy - (p * r_sum)) / denominator, 0.0);
        }
    }
}

static int
mbe_fa_index_of(double p) {
    int idx = (int)round((p - 21.0) / 0.5);
    return (idx < MBE_FA_CANDIDATES - 1) ? idx : MBE_FA_CANDIDATES - 1;
}

static double
mbe_fa_nearest_candidate(double p) {
    double clamped = (p < 21.0) ? 21.0 : ((p > 122.0) ? 122.0 : p);
    return (round((clamped - 21.0) / 0.5) * 0.5) + 21.0;
}

/* Look-back tracking (eq 10-12). */
static double
mbe_fa_look_back(const double e0[MBE_FA_CANDIDATES], double p_prev1, double e_prev1, double e_prev2, double* ce_b) {
    double lo = 0.8 * p_prev1;
    double hi = 1.2 * p_prev1;
    double best_p = 0.0;
    double best_e = 0.0;
    int found = 0;
    for (int i = 0; i < MBE_FA_CANDIDATES; i++) {
        double p = mbe_fa_candidate(i);
        if (p >= lo && p <= hi && (!found || e0[i] < best_e)) {
            best_p = p;
            best_e = e0[i];
            found = 1;
        }
    }
    *ce_b = best_e + e_prev1 + e_prev2;
    return best_p;
}

/* Look-ahead tracking (eq 13-20). */
static double
mbe_fa_look_ahead(const struct mbe_fa_tables* tables, const double e0[MBE_FA_CANDIDATES],
                  const double e1[MBE_FA_CANDIDATES], const double e2[MBE_FA_CANDIDATES], double* ce_f_out) {
    double best_e2[MBE_FA_CANDIDATES];
    double ce_f[MBE_FA_CANDIDATES];

    /* The inner minimum over P2 depends only on P1, so it is computed once per
     * P1 and reused for every P0. */
    for (int i = 0; i < MBE_FA_CANDIDATES; i++) {
        /* Every range holds its own candidate, so it is never empty. */
        double m = e2[tables->range_lo[i]];
        for (int j = tables->range_lo[i] + 1; j <= tables->range_hi[i]; j++) {
            m = fmin(m, e2[j]);
        }
        best_e2[i] = m;
    }
    for (int i0 = 0; i0 < MBE_FA_CANDIDATES; i0++) {
        double m = e1[tables->range_lo[i0]] + best_e2[tables->range_lo[i0]];
        for (int i1 = tables->range_lo[i0] + 1; i1 <= tables->range_hi[i0]; i1++) {
            m = fmin(m, e1[i1] + best_e2[i1]);
        }
        ce_f[i0] = e0[i0] + m;
    }

    int best = 0;
    for (int i = 1; i < MBE_FA_CANDIDATES; i++) {
        if (ce_f[i] < ce_f[best]) {
            best = i;
        }
    }
    double p_hat_0 = mbe_fa_candidate(best);
    double ce_f_p_hat_0 = ce_f[best];

    /* Sub-multiples P_hat_0 / n snapped to the candidate set, smallest first. */
    double submultiples[8];
    int count = 0;
    for (int n = 2; count < 8; n++) {
        double raw = p_hat_0 / (double)n;
        if (raw < 21.0) {
            break;
        }
        submultiples[count++] = mbe_fa_nearest_candidate(raw);
    }
    for (int i = count - 1; i >= 0; i--) {
        double candidate = submultiples[i];
        double ce = ce_f[mbe_fa_index_of(candidate)];
        /* Eq 18-19 are evaluated in multiplied form, which remains defined
         * when CE_F(P_hat_0) is zero. */
        int satisfies_18 = (ce <= 0.85) && (ce <= 1.7 * ce_f_p_hat_0);
        int satisfies_19 = (ce <= 0.4) && (ce <= 3.5 * ce_f_p_hat_0);
        int satisfies_20 = (ce <= 0.05);
        if (satisfies_18 || satisfies_19 || satisfies_20) {
            *ce_f_out = ce;
            return candidate;
        }
    }
    *ce_f_out = ce_f_p_hat_0;
    return p_hat_0;
}

/* S_w(m) (eq 29), the 256-point DFT of s(n) w_R(n), for m = -127..128. */
static void
mbe_fa_spectrum(const struct mbe_fa_tables* tables, const double* buf, int center, struct mbe_fa_spectrum* sw) {
    double windowed[221];
    for (int n = -110; n <= 110; n++) {
        windowed[n + 110] = buf[center + n] * mbe_fa_wr(n);
    }
    for (int i = 0; i < MBE_FA_BINS; i++) {
        int m = i - 127;
        double re = 0.0;
        double im = 0.0;
        for (int j = 0; j < 221; j++) {
            int n = j - 110;
            int k = (m * n) % 256;
            if (k < 0) {
                k += 256;
            }
            re += windowed[j] * tables->twiddle_re[k];
            im += windowed[j] * tables->twiddle_im[k];
        }
        sw->re[i] = re;
        sw->im[i] = im;
    }
}

static void
mbe_fa_sw_at(const struct mbe_fa_spectrum* sw, int m, double* re, double* im) {
    if (m >= -127 && m <= 128) {
        *re = sw->re[m + 127];
        *im = sw->im[m + 127];
    } else {
        *re = 0.0;
        *im = 0.0;
    }
}

/* a_l and b_l (eq 26-27 and 32-33), in DFT bins. */
static double
mbe_fa_band_lo(int l, double omega0) {
    return (256.0 / (2.0 * M_PI)) * ((double)l - 0.5) * omega0;
}

static double
mbe_fa_band_hi(int l, double omega0) {
    return (256.0 / (2.0 * M_PI)) * ((double)l + 0.5) * omega0;
}

static int
mbe_fa_wr_index(int m, int l, double omega0) {
    return (int)floor((64.0 * (double)m) - ((16384.0 / (2.0 * M_PI)) * (double)l * omega0) + 0.5);
}

/* The synthetic spectrum S_w(m, omega0) of eq 25, caching A_l (eq 28) per harmonic. */
struct mbe_fa_synthetic {
    const struct mbe_fa_tables* tables;
    const struct mbe_fa_spectrum* sw;
    double omega0;
    int max_l;
    unsigned char have[MBE_FA_MAX_HARMONICS + 2];
    double amp_re[MBE_FA_MAX_HARMONICS + 2];
    double amp_im[MBE_FA_MAX_HARMONICS + 2];
};

static void
mbe_fa_synthetic_init(struct mbe_fa_synthetic* syn, const struct mbe_fa_tables* tables,
                      const struct mbe_fa_spectrum* sw, double omega0, int max_l) {
    syn->tables = tables;
    syn->sw = sw;
    syn->omega0 = omega0;
    syn->max_l = (max_l <= MBE_FA_MAX_HARMONICS + 1) ? max_l : MBE_FA_MAX_HARMONICS + 1;
    memset(syn->have, 0, sizeof(syn->have));
    memset(syn->amp_re, 0, sizeof(syn->amp_re));
    memset(syn->amp_im, 0, sizeof(syn->amp_im));
}

/* A_l(omega0) (eq 28): least-squares harmonic amplitude from the bins of its band. */
static void
mbe_fa_harmonic_amplitude(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, int l, double omega0,
                          double* out_re, double* out_im) {
    int m_lo = (int)ceil(mbe_fa_band_lo(l, omega0));
    int m_hi = (int)ceil(mbe_fa_band_hi(l, omega0));
    double num_re = 0.0;
    double num_im = 0.0;
    double denominator = 0.0;
    for (int m = m_lo; m < m_hi; m++) {
        double wr = mbe_fa_wr_dft(tables, mbe_fa_wr_index(m, l, omega0));
        double re;
        double im;
        mbe_fa_sw_at(sw, m, &re, &im);
        num_re += re * wr;
        num_im += im * wr;
        denominator += wr * wr;
    }
    if (fabs(denominator) < 1e-12) {
        *out_re = 0.0;
        *out_im = 0.0;
    } else {
        double inv = 1.0 / denominator;
        *out_re = num_re * inv;
        *out_im = num_im * inv;
    }
}

static void
mbe_fa_synthetic_at(struct mbe_fa_synthetic* syn, int m, double* re, double* im) {
    for (int l = 0; l <= syn->max_l; l++) {
        int m_lo = (int)ceil(mbe_fa_band_lo(l, syn->omega0));
        int m_hi = (int)ceil(mbe_fa_band_hi(l, syn->omega0));
        if (m >= m_lo && m < m_hi) {
            if (!syn->have[l]) {
                mbe_fa_harmonic_amplitude(syn->tables, syn->sw, l, syn->omega0, &syn->amp_re[l], &syn->amp_im[l]);
                syn->have[l] = 1;
            }
            double wr = mbe_fa_wr_dft(syn->tables, mbe_fa_wr_index(m, l, syn->omega0));
            *re = syn->amp_re[l] * wr;
            *im = syn->amp_im[l] * wr;
            return;
        }
    }
    *re = 0.0;
    *im = 0.0;
}

/* E_R(omega0) (eq 24). */
static double
mbe_fa_refinement_error(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, double omega0) {
    double l_estimate = floor((0.9254 * M_PI / omega0) - 0.5);
    int upper_m = (int)floor(l_estimate * (256.0 / (2.0 * M_PI)) * omega0);
    int max_l = (int)fmax(l_estimate, 0.0) + 1;
    struct mbe_fa_synthetic syn;
    double sum = 0.0;

    mbe_fa_synthetic_init(&syn, tables, sw, omega0, max_l);
    for (int m = 50; m <= upper_m; m++) {
        double re;
        double im;
        double s_re;
        double s_im;
        mbe_fa_sw_at(sw, m, &re, &im);
        mbe_fa_synthetic_at(&syn, m, &s_re, &s_im);
        double d_re = re - s_re;
        double d_im = im - s_im;
        sum += (d_re * d_re) + (d_im * d_im);
    }
    return sum;
}

/* Quarter-sample refinement (5.1.5): the best of P_hat_I +- 1/8 .. 9/8. */
static double
mbe_fa_refine_pitch(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, double p_hat_i) {
    static const double offsets[10] = {-9.0 / 8.0, -7.0 / 8.0, -5.0 / 8.0, -3.0 / 8.0, -1.0 / 8.0,
                                       1.0 / 8.0,  3.0 / 8.0,  5.0 / 8.0,  7.0 / 8.0,  9.0 / 8.0};
    double best_omega0 = 0.0;
    double best_error = 0.0;
    for (int i = 0; i < 10; i++) {
        double omega0 = 2.0 * M_PI / (p_hat_i + offsets[i]);
        double e = mbe_fa_refinement_error(tables, sw, omega0);
        if (i == 0 || e < best_error) {
            best_error = e;
            best_omega0 = omega0;
        }
    }
    return best_omega0;
}

int
mbe_fa_push(const struct mbe_fa_tables* tables, struct mbe_fa_state* state, const double input[MBE_FA_FRAME],
            struct mbe_fa_result* result) {
    memmove(state->buf, state->buf + MBE_FA_FRAME, (MBE_FA_BUFFER - MBE_FA_FRAME) * sizeof(state->buf[0]));
    for (int i = 0; i < MBE_FA_FRAME; i++) {
        state->buf[MBE_FA_BUFFER - MBE_FA_FRAME + i] = mbe_fa_high_pass(state, input[i]);
    }

    /* E(P) of the newest frame whose margin is complete (k + 2). */
    memmove(state->error[0], state->error[1], 2 * sizeof(state->error[0]));
    mbe_fa_error_table(state->buf, MBE_FA_CENTER + (2 * MBE_FA_FRAME), state->error[2]);

    if (state->pushed < MBE_FA_LOOKAHEAD) {
        state->pushed++;
        return 0;
    }

    const double* e0 = state->error[0];
    double ce_b;
    double ce_f;
    double p_b = mbe_fa_look_back(e0, state->prev1_pitch, state->prev1_error, state->prev2_error, &ce_b);
    double p_f = mbe_fa_look_ahead(tables, e0, state->error[1], state->error[2], &ce_f);
    /* Eq 21-23. */
    double p_initial = ((ce_b <= 0.48) || (ce_b <= ce_f)) ? p_b : p_f;
    double e_initial = e0[mbe_fa_index_of(p_initial)];

    mbe_fa_spectrum(tables, state->buf, MBE_FA_CENTER, &result->sw);
    result->omega0_hat = mbe_fa_refine_pitch(tables, &result->sw, p_initial);
    result->initial_pitch = p_initial;
    result->initial_pitch_error = e_initial;
    memcpy(result->slot, state->buf + MBE_FA_CENTER, sizeof(result->slot));

    state->prev2_pitch = state->prev1_pitch;
    state->prev2_error = state->prev1_error;
    state->prev1_pitch = p_initial;
    state->prev1_error = e_initial;
    return 1;
}

int
mbe_fa_harmonics_count(double omega0) {
    double inner = floor((M_PI / omega0) + 0.25);
    return (int)floor(0.9254 * inner);
}

int
mbe_fa_bands_count(int harmonics) {
    if (harmonics <= 36) {
        return (harmonics + 2) / 3;
    }
    return 12;
}

/* D_k (eq 35-36): the fraction of the band's energy the harmonic model fails to explain. */
static double
mbe_fa_voicing_measure(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, int k, int harmonics,
                       double omega0, int is_highest) {
    int m_lo = (int)ceil(mbe_fa_band_lo((3 * k) - 2, omega0));
    int upper_l = is_highest ? harmonics : 3 * k;
    int m_hi = (int)ceil(mbe_fa_band_hi(upper_l, omega0));
    double error_energy = 0.0;
    double real_energy = 0.0;
    struct mbe_fa_synthetic syn;

    mbe_fa_synthetic_init(&syn, tables, sw, omega0, harmonics);
    for (int m = m_lo; m < m_hi; m++) {
        double re;
        double im;
        double s_re;
        double s_im;
        mbe_fa_sw_at(sw, m, &re, &im);
        mbe_fa_synthetic_at(&syn, m, &s_re, &s_im);
        double d_re = re - s_re;
        double d_im = im - s_im;
        error_energy += (d_re * d_re) + (d_im * d_im);
        real_energy += (re * re) + (im * im);
    }
    if (fabs(real_energy) < 1e-12) {
        return 1.0;
    }
    return error_energy / real_energy;
}

int
mbe_fa_determine_voicing(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, double omega0,
                         double initial_pitch_error, struct mbe_fa_voicing_state* voicing,
                         unsigned char bands[MBE_FA_MAX_BANDS]) {
    int harmonics = mbe_fa_harmonics_count(omega0);
    int k_hat = mbe_fa_bands_count(harmonics);
    double wr0 = mbe_fa_wr_dft(tables, 0);
    double wr0_sq = wr0 * wr0;
    double xi_lf = 0.0;
    double xi_hf = 0.0;

    if (k_hat > MBE_FA_MAX_BANDS) {
        k_hat = MBE_FA_MAX_BANDS;
    }
    if (k_hat < 0) {
        k_hat = 0;
    }

    /* Eq 38-40. */
    for (int m = 0; m <= 63; m++) {
        xi_lf += (sw->re[m + 127] * sw->re[m + 127]) + (sw->im[m + 127] * sw->im[m + 127]);
    }
    xi_lf /= wr0_sq;
    for (int m = 64; m <= 128; m++) {
        xi_hf += (sw->re[m + 127] * sw->re[m + 127]) + (sw->im[m + 127] * sw->im[m + 127]);
    }
    xi_hf /= wr0_sq;
    double xi_0 = xi_lf + xi_hf;

    /* Eq 41. */
    double xi_max;
    if (xi_0 > voicing->xi_max) {
        xi_max = (0.5 * voicing->xi_max) + (0.5 * xi_0);
    } else {
        double decayed = (0.99 * voicing->xi_max) + (0.01 * xi_0);
        xi_max = (decayed > 20000.0) ? decayed : 20000.0;
    }

    /* Eq 42. */
    double m_xi = ((0.0025 * xi_max) + xi_0) / ((0.01 * xi_max) + xi_0);
    if (!(xi_lf >= 5.0 * xi_hf)) {
        m_xi *= sqrt(xi_lf / (5.0 * xi_hf));
    }

    for (int k = 1; k <= k_hat; k++) {
        double d_k = mbe_fa_voicing_measure(tables, sw, k, harmonics, omega0, k == k_hat);
        int previous_voiced = (k - 1 < voicing->prev_count) && voicing->prev[k - 1];
        double theta;
        /* Eq 37. */
        if (initial_pitch_error > 0.5 && k >= 2) {
            theta = 0.0;
        } else if (previous_voiced) {
            theta = 0.5625 * (1.0 - (0.3096 * ((double)k - 1.0) * omega0)) * m_xi;
        } else {
            theta = 0.45 * (1.0 - (0.3096 * ((double)k - 1.0) * omega0)) * m_xi;
        }
        bands[k - 1] = (unsigned char)(d_k < theta);
    }

    voicing->xi_max = xi_max;
    voicing->prev_count = k_hat;
    memcpy(voicing->prev, bands, (size_t)k_hat);
    return k_hat;
}

/* Eq 43. */
static double
mbe_fa_voiced_amplitude(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, int l, double omega0) {
    int m_lo = (int)ceil(mbe_fa_band_lo(l, omega0));
    int m_hi = (int)ceil(mbe_fa_band_hi(l, omega0));
    double signal_energy = 0.0;
    double window_energy = 0.0;
    for (int m = m_lo; m < m_hi; m++) {
        double re;
        double im;
        mbe_fa_sw_at(sw, m, &re, &im);
        signal_energy += (re * re) + (im * im);
        double wr = mbe_fa_wr_dft(tables, mbe_fa_wr_index(m, l, omega0));
        window_energy += wr * wr;
    }
    if (fabs(window_energy) < 1e-12) {
        return 0.0;
    }
    return sqrt(signal_energy / window_energy);
}

/* Eq 44. */
static double
mbe_fa_unvoiced_amplitude(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, int l, double omega0) {
    int m_lo = (int)ceil(mbe_fa_band_lo(l, omega0));
    int m_hi = (int)ceil(mbe_fa_band_hi(l, omega0));
    double band_width = (double)(m_hi - m_lo);
    double signal_energy = 0.0;
    if (band_width <= 0.0) {
        return 0.0;
    }
    for (int m = m_lo; m < m_hi; m++) {
        double re;
        double im;
        mbe_fa_sw_at(sw, m, &re, &im);
        signal_energy += (re * re) + (im * im);
    }
    return (1.0 / tables->wr_tap_sum) * sqrt(signal_energy / band_width);
}

void
mbe_fa_spectral_amplitudes(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, int harmonics,
                           int bands_count, double omega0, const unsigned char* bands,
                           double amplitudes[MBE_FA_MAX_HARMONICS + 1]) {
    if (harmonics > MBE_FA_MAX_HARMONICS) {
        harmonics = MBE_FA_MAX_HARMONICS;
    }
    amplitudes[0] = 0.0;
    for (int l = 1; l <= harmonics; l++) {
        int k = (l + 2) / 3;
        if (k > bands_count) {
            k = bands_count;
        }
        if (k < 1) {
            k = 1;
        }
        int voiced = (k <= bands_count) && bands[k - 1];
        amplitudes[l] =
            voiced ? mbe_fa_voiced_amplitude(tables, sw, l, omega0) : mbe_fa_unvoiced_amplitude(tables, sw, l, omega0);
    }
}
