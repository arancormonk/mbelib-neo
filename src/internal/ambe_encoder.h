// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Spectral amplitude quantization and FEC shared by the AMBE encoders
 *        (private API).
 *
 * D-STAR (AMBE 3600x2400) and AMBE+2 (3600x2450) quantize spectral amplitudes
 * the same way (TIA-102.BABA-1 4.3): a gain step against half the previous
 * gain, a prediction residual against the previous frame's log2 amplitudes,
 * four block DCTs, the PRBA vector and the higher-order coefficients. Only the
 * tables, their sizes and the prediction coefficient differ, so each codec
 * passes a struct ambe_enc_tables. Both also share the 72-bit frame's FEC.
 */
#ifndef MBELIB_NEO_INTERNAL_AMBE_ENCODER_H
#define MBELIB_NEO_INTERNAL_AMBE_ENCODER_H

#include "mbelib-neo/mbelib.h"

/** One codec's quantizer tables, as its decoder dequantizes them. */
struct ambe_enc_tables {
    const float* dg; /* gain steps */
    int dg_count;
    const float (*prba24)[3];
    const float (*prba58)[4];
    const float (*hoc[4])[4]; /* b5..b8 */
    int hoc_rows[4];
    int hoc_step[4];        /* row stride: D-STAR's b8 reaches only even rows */
    const int (*lmprbl)[4]; /* block lengths by L */
    float rho;              /* prediction coefficient */
    int fit_floor_gain;     /* see ambe_enc_quantize_gain() */
};

/** Per-frame quantization workspace; no inter-frame state. */
struct ambe_enc_frame {
    const struct ambe_dct_cache* cache;
    float p[57];            /* interpolated previous log2 amplitudes */
    float a[57];            /* target log2 amplitudes, less 0.5 log2 L */
    float Tl[57], Tl_q[57]; /* prediction residual and its reconstruction */
    float Cik[5][18], Cik_q[5][18], Gm[9];
    int Ji[5], L, b[9];
    float f0q, gamma_q, mean_a, mean_p;
};

/* Mean of a[1..L] into mean_a. */
void ambe_enc_mean_amplitude(struct ambe_enc_frame* q);

/*
 * b2: the gain step nearest gamma - gamma(-1) / 2 (TIA-102.BABA-1 eqs 9-10),
 * and gamma_q. With fit_floor_gain, a frame whose target lies below the first
 * step has its envelope flattened toward the level that step can play
 * (mbe_encoder_fit_floor()).
 */
void ambe_enc_quantize_gain(struct ambe_enc_frame* q, const struct ambe_enc_tables* t, const mbe_parms* prev_mp);

/* Prediction residual Tl against the decoder's previous log2 amplitudes. */
void ambe_enc_prediction(struct ambe_enc_frame* q, const struct ambe_enc_tables* t, const mbe_parms* prev_mp);

/* b3..b8: block DCTs, the PRBA vector and the higher-order coefficients. */
void ambe_enc_quantize_residual(struct ambe_enc_frame* q, const struct ambe_enc_tables* t);

/*
 * The 72-bit frame of 49 parameter bits: C0 as (23,12) Golay plus even parity,
 * C1 as (23,12) Golay scrambled by the sequence seeded from C0's data, and C2
 * and C3 unprotected, in the plane layout the decoders read. Returns the
 * scrambled even parity of C1's codeword, which D-STAR sends in place of
 * ambe_d[24].
 */
int ambe_enc_fec(const char ambe_d[49], char ambe_fr[4][24]);

/*
 * AMBE+2 steps exposed for the tests (ambe3600x2450_enc.c). b0 is a voice
 * pitch code (0..119) and m and a are indexed by harmonic 1..L of its
 * AmbeLtable count.
 *
 * mbe_ambe2450_quantize_vuv: b1 for amplitudes m and the band decisions of
 * mbe_analysis_band_voicing() (TIA-102.BABA-1 4.2).
 * mbe_ambe2450_quantize_amplitudes: b2..b8 for log2 amplitudes a in the
 * decoder's domain (its log2Ml, eq 8 less 0.5 log2 L) against prev_mp, and
 * all 49 bits packed.
 */
int mbe_ambe2450_quantize_vuv(int b0, const float m[57], const unsigned char bands[12]);
void mbe_ambe2450_quantize_amplitudes(int b0, int b1, const float a[57], const mbe_parms* prev_mp, char ambe_d[49]);

#endif /* MBELIB_NEO_INTERNAL_AMBE_ENCODER_H */
