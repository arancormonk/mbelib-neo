// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief MBE speech analysis for the encoders (private API).
 *
 * Pitch, per-harmonic voicing measures and spectral magnitudes following the
 * method of ANSI/TIA-102.BABA chapter 5, with the voicing-independent
 * magnitude estimate of US 5,701,390. All signals are on the 16-bit sample
 * scale the standard's thresholds assume.
 */
#ifndef MBELIB_NEO_INTERNAL_MBE_SPEECH_ANALYSIS_H
#define MBELIB_NEO_INTERNAL_MBE_SPEECH_ANALYSIS_H

#include "mbe_unvoiced_fft.h"

#define MBE_ANALYSIS_FRAME       160
#define MBE_ANALYSIS_PITCH_HALF  150 /* w_I spans -150..150 */
#define MBE_ANALYSIS_WR_HALF     110 /* w_R spans -110..110 */
#define MBE_ANALYSIS_LPF_HALF    10  /* pitch low-pass spans -10..10 */
/* Past samples kept so w_I plus the low-pass fit around the analysis centre,
 * one sample before the current frame: samples -161..159. */
#define MBE_ANALYSIS_HISTORY     (MBE_ANALYSIS_PITCH_HALF + MBE_ANALYSIS_LPF_HALF + 1)
#define MBE_ANALYSIS_BUFFER      (MBE_ANALYSIS_HISTORY + MBE_ANALYSIS_FRAME)
#define MBE_ANALYSIS_WTAB_STEPS  64 /* W_R samples per FFT bin, as eq (25) indexes it */
#define MBE_ANALYSIS_WTAB_BINS   7  /* |nu| below 7 bins covers every harmonic band */
#define MBE_ANALYSIS_WTAB_LEN    ((2 * MBE_ANALYSIS_WTAB_BINS * MBE_ANALYSIS_WTAB_STEPS) + 1)
#define MBE_ANALYSIS_PITCH_COUNT 217 /* P = 20, 20.5, ..., 128 */
#define MBE_ANALYSIS_HARMONICS   56
#define MBE_ANALYSIS_COLUMNS     8 /* 500 Hz voicing columns */

/** Constant tables, computed once per encoder context. */
struct mbe_analysis_tables {
    float w_i[(2 * MBE_ANALYSIS_PITCH_HALF) + 1];
    float w_r[(2 * MBE_ANALYSIS_WR_HALF) + 1];
    float lpf[(2 * MBE_ANALYSIS_LPF_HALF) + 1];
    float w_r_dtft[MBE_ANALYSIS_WTAB_LEN]; /* W_R(nu), nu in 1/64-bin steps */
    float w_i4_sum;                        /* sum of w_I^4 */
    float w_r_sum;                         /* W_R(0) */
    float w_r_sq_sum;                      /* sum of w_R^2 */
};

/** Per-stream analysis state. */
struct mbe_analysis_state {
    float hp_in;                                      /* previous input of the eq (3) DC filter */
    float hp_out;                                     /* previous output of the eq (3) DC filter */
    float buf[MBE_ANALYSIS_BUFFER];                   /* DC-filtered samples; the newest frame at the end */
    float pitch_prev[2];                              /* look-back history P(-1), P(-2) */
    float error_prev[2];                              /* look-back history E(-1), E(-2) */
    int trusted;                                      /* previous frame gave a periodic pitch track */
    float xi_max;                                     /* eq (41) energy tracker */
    unsigned char columns_prev[MBE_ANALYSIS_COLUMNS]; /* last transmitted voicing per column */
};

/** One frame's analysis. Arrays are indexed by harmonic number (1-based). */
struct mbe_analysis_result {
    float f0;                                    /* refined fundamental, cycles per sample */
    float pitch_error;                           /* E at the tracked pitch */
    int harmonics;                               /* harmonics below Nyquist at f0 */
    float magnitude[MBE_ANALYSIS_HARMONICS + 1]; /* voicing-independent M_l */
    float fit_error[MBE_ANALYSIS_HARMONICS + 1]; /* normalized fit error D_l */
    float energy_factor;                         /* M(xi) of eq (42) */
};

void mbe_analysis_init_tables(struct mbe_analysis_tables* tables);
void mbe_analysis_reset(struct mbe_analysis_state* state);

/* DC-filter one frame (16-bit scale) into the history. */
void mbe_analysis_push(struct mbe_analysis_state* state, const float input[MBE_ANALYSIS_FRAME]);

/* Analyze the window centred one sample before the newest frame. Returns 0
 * or a negative MBE_STATUS_* value. */
int mbe_analysis_frame(const struct mbe_analysis_tables* tables, struct mbe_analysis_state* state, mbe_fft_plan* fft,
                       mbe_acf_plan* acf, struct mbe_analysis_result* result);

/* Voicing threshold Theta(k, 0.1309) of eq (37) for 1-based column k, with
 * hysteresis from the previously transmitted column voicing. */
float mbe_analysis_column_threshold(const struct mbe_analysis_state* state, const struct mbe_analysis_result* result,
                                    int column);

/* Record the transmitted voicing for hysteresis. */
void mbe_analysis_commit(struct mbe_analysis_state* state, const unsigned char columns[MBE_ANALYSIS_COLUMNS]);

#endif /* MBELIB_NEO_INTERNAL_MBE_SPEECH_ANALYSIS_H */
