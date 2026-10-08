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
 * @brief TIA-102.BABA speech analysis shared by the IMBE and AMBE+2 encoders (private API).
 *
 * C port of the floating-point analysis in Bruce Perens' ham_digital_modes
 * (hams_open, daemons/ham_digital_modes/src/ambe/float/tia_102_baba): the
 * input high-pass filter (eq 3), the pitch error function E(P) (eq 5) with
 * look-back and look-ahead tracking (5.1.2-5.1.4), quarter-sample pitch
 * refinement (5.1.5), voiced/unvoiced determination (5.2) and spectral
 * amplitude estimation (5.3). Computation is in double precision and follows
 * the reference's order of operations, so that for the same input it produces
 * the same results as the reference.
 *
 * The analysis looks two frames ahead: the frame returned after pushing
 * input frame n is centred on the first sample of input frame n - 3. The
 * first three pushes return no analysis.
 */
#ifndef MBELIB_NEO_INTERNAL_MBE_FRAME_ANALYSIS_H
#define MBELIB_NEO_INTERNAL_MBE_FRAME_ANALYSIS_H

#include <stdint.h>

#define MBE_FA_FRAME         160
#define MBE_FA_MARGIN        160 /* E(P) reach either side of a centre: 150 window + 10 filter taps */
#define MBE_FA_BUFFER        800 /* frame k centred at 160, its look-ahead frames at 320 and 480 */
#define MBE_FA_CENTER        160
#define MBE_FA_LOOKAHEAD     3   /* pushes before the first analysed frame */
#define MBE_FA_CANDIDATES    203 /* P = 21, 21.5, ..., 122 */
#define MBE_FA_R_LAGS        151
#define MBE_FA_BINS          256  /* S_w(m), m = -127..128 */
#define MBE_FA_WR_DFT_HALF   8192 /* W_R(m) tabulated for |m| <= 8192 */
#define MBE_FA_MAX_BANDS     12
#define MBE_FA_MAX_HARMONICS 56

/** Constant tables, computed once per encoder context. */
struct mbe_fa_tables {
    double wr_dft[MBE_FA_WR_DFT_HALF + 1]; /* W_R(m), eq 30 */
    double twiddle_re[MBE_FA_BINS];        /* exp(-2 pi j k / 256) */
    double twiddle_im[MBE_FA_BINS];
    int range_lo[MBE_FA_CANDIDATES]; /* candidates within 0.8 P .. 1.2 P (eq 14 and 16) */
    int range_hi[MBE_FA_CANDIDATES];
    double wr_tap_sum; /* sum of w_R(n), eq 44 */
};

/** The windowed spectrum S_w(m) of eq 29, m = -127..128 at index m + 127. */
struct mbe_fa_spectrum {
    double re[MBE_FA_BINS];
    double im[MBE_FA_BINS];
};

/** Per-stream pitch analysis state. */
struct mbe_fa_state {
    int64_t hp_prev_input;              /* eq 3 high-pass filter, integer as in the reference */
    int64_t hp_prev_output_q16;         /* its previous output, Q16 */
    double buf[MBE_FA_BUFFER];          /* filtered input, the newest frame at the end */
    double error[3][MBE_FA_CANDIDATES]; /* E(P) of frames k, k+1, k+2 */
    double prev1_pitch, prev1_error;    /* P(-1), E(-1) */
    double prev2_pitch, prev2_error;    /* P(-2), E(-2) */
    int pushed;                         /* pushes so far, saturating at MBE_FA_LOOKAHEAD */
};

/** One analysed frame. */
struct mbe_fa_result {
    double omega0_hat;          /* refined fundamental, radians per sample */
    double initial_pitch;       /* P_hat_I, samples */
    double initial_pitch_error; /* E(P_hat_I) */
    struct mbe_fa_spectrum sw;
    double slot[MBE_FA_FRAME]; /* the frame's own 160 input samples, for tone detection */
};

/** Voicing state carried between frames: eq 41's energy tracker and the previous band decisions. */
struct mbe_fa_voicing_state {
    double xi_max;
    int prev_count;
    unsigned char prev[MBE_FA_MAX_BANDS];
};

void mbe_fa_init_tables(struct mbe_fa_tables* tables);
void mbe_fa_reset(struct mbe_fa_state* state);
void mbe_fa_voicing_reset(struct mbe_fa_voicing_state* voicing);

/*
 * High-pass filter and buffer one frame of 16-bit scale input, then analyse
 * the frame MBE_FA_LOOKAHEAD frames back. Returns 1 with result filled, or 0
 * while the look-ahead is still filling.
 */
int mbe_fa_push(const struct mbe_fa_tables* tables, struct mbe_fa_state* state, const double input[MBE_FA_FRAME],
                struct mbe_fa_result* result);

/** L_hat (eq 31) for a fundamental in radians per sample. */
int mbe_fa_harmonics_count(double omega0);

/** K_hat (eq 34). */
int mbe_fa_bands_count(int harmonics);

/*
 * V/UV decision per band (5.2) at fundamental omega0, updating the voicing
 * state. Writes mbe_fa_bands_count(mbe_fa_harmonics_count(omega0)) decisions
 * and returns that count.
 */
int mbe_fa_determine_voicing(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, double omega0,
                             double initial_pitch_error, struct mbe_fa_voicing_state* voicing,
                             unsigned char bands[MBE_FA_MAX_BANDS]);

/*
 * Spectral amplitudes M_hat_l (5.3) for l = 1..harmonics into amplitudes[l],
 * each by the voiced or unvoiced estimator of its band (bands has bands_count
 * entries).
 */
void mbe_fa_spectral_amplitudes(const struct mbe_fa_tables* tables, const struct mbe_fa_spectrum* sw, int harmonics,
                                int bands_count, double omega0, const unsigned char* bands,
                                double amplitudes[MBE_FA_MAX_HARMONICS + 1]);

#endif /* MBELIB_NEO_INTERNAL_MBE_FRAME_ANALYSIS_H */
