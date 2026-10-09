// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Analysis front end shared by the speech encoders (private API).
 *
 * Every encoder validates and scales its input the same way, DC-filters it
 * into one mbe_speech_analysis.c history and analyses the window centred one
 * sample before the newest frame. The front end owns that state, the analysis
 * tables and the FFT and autocorrelation plans.
 */
#ifndef MBELIB_NEO_INTERNAL_MBE_ENCODER_H
#define MBELIB_NEO_INTERNAL_MBE_ENCODER_H

#include "mbe_speech_analysis.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

#define MBE_ENCODER_SAMPLES MBE_ANALYSIS_FRAME

struct mbe_encoder_frontend {
    struct mbe_analysis_state analysis;
    mbe_fft_plan* fft;
    mbe_acf_plan* acf;
    struct mbe_analysis_tables tables;
};

/* Allocate the plans and compute the tables. Returns 0, or -1 with nothing
 * left allocated. The front end must start zeroed. */
int mbe_encoder_frontend_open(struct mbe_encoder_frontend* fe);
/* Free the plans; a zeroed front end is accepted. */
void mbe_encoder_frontend_close(struct mbe_encoder_frontend* fe);
/* Restore the freshly opened analysis state (the tables are kept). */
void mbe_encoder_frontend_reset(struct mbe_encoder_frontend* fe);

/* Finite and within +-2^20, checked on the bit pattern so the test survives
 * fast-math. One bad sample would otherwise poison the DC filter and the
 * analysis history for the rest of the stream. */
int mbe_encoder_samples_valid(const float samples[MBE_ENCODER_SAMPLES]);

/* The prediction history the encoders read from the caller's prev_mp
 * (log2Ml[0..56] and gamma) is finite and within +-2^20, so the arithmetic on
 * it stays finite; a malformed state is rejected before it reaches the
 * analysis. */
int mbe_encoder_history_valid(const mbe_parms* prev_mp);

/* The model an encoder returns is finite (w0, gamma, and log2Ml and Ml over
 * its harmonics). Decoders never leave it otherwise from any state they can
 * reach; a crafted history can overflow it, and the encoders reject that. */
int mbe_encoder_model_valid(const mbe_parms* mp);

/* 16-bit PCM to the float input scale ([-1, 1)). */
void mbe_encoder_short_to_float(const short in[MBE_ENCODER_SAMPLES], float out[MBE_ENCODER_SAMPLES]);

/* DC-filter one validated frame (float scale) into the history on the 16-bit
 * scale the analysis expects, then analyse it. Returns 0 or a negative
 * MBE_STATUS_* value. */
int mbe_encoder_frontend_analyze(struct mbe_encoder_frontend* fe, const float samples[MBE_ENCODER_SAMPLES],
                                 struct mbe_analysis_result* result);

/*
 * A codec's lowest gain step puts a floor under the mean log2 amplitude it can
 * send. A frame whose own mean lies below it (a quiet tone: one strong harmonic
 * over empty bands) would decode with the excess pushed into its peak. This
 * scales the deviations of a[1..L] from their mean mean_a toward a flat
 * envelope at floor_mean until its energy matches the frame's (a flat envelope
 * is the quietest the floor can play), and returns floor_mean, the new mean.
 */
float mbe_encoder_fit_floor(float* a, int L, float mean_a, float floor_mean);

#endif /* MBELIB_NEO_INTERNAL_MBE_ENCODER_H */
