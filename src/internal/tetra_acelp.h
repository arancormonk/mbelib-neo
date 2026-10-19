// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Internal interface of the TETRA full-rate ACELP decoder (EN 300 395-2 clause 4).
 *
 * The public entry points live in mbelib.h. This header exposes the decoder state and the
 * individual decoding stages so the unit tests can exercise them separately.
 */

#ifndef MBELIB_NEO_INTERNAL_TETRA_ACELP_H
#define MBELIB_NEO_INTERNAL_TETRA_ACELP_H

#include <stdint.h>

#define TETRA_M                10 /* LP order */
#define TETRA_SUBFRAME         60 /* samples per subframe */
#define TETRA_SUBFRAMES        4
#define TETRA_FRAME            240 /* samples per frame */
#define TETRA_BITS             137 /* B1..B137 */
#define TETRA_PIT_MIN          20
#define TETRA_PIT_MAX          143
#define TETRA_INTERP_LEN       16 /* filter half length in samples */
/* Past excitation kept for the adaptive codebook: the largest lag (143 + 1 from a -1/3 fraction
 * on a clamped relative lag) plus the filter half length, rounded up. */
#define TETRA_EXC_HIST         (TETRA_PIT_MAX + 1 + TETRA_INTERP_LEN + 1)
/* Both interpolation filters' first 18 taps have an absolute sum <= 65535.
 * Their MAC prefixes cannot overflow int32 for any int16 samples. Later taps
 * can overflow and must saturate in order. test_tetra_simd checks this bound. */
#define TETRA_INTERP_SAFE_TAPS 18

/** Decoded parameters of one subframe (EN 300 395-2 table 3). */
typedef struct tetra_subframe_params {
    int pitch;  /* 8-bit (first subframe) or 5-bit relative pitch index */
    int pulse1; /* 5 bits */
    int pulse2; /* 3 bits */
    int pulse3; /* 3 bits */
    int pulse4; /* 3 bits */
    int sign;   /* global sign flag */
    int shift;  /* pulse position shift */
    int gain;   /* 6-bit gain index */
} tetra_subframe_params;

/** Decoded parameters of one frame. */
typedef struct tetra_frame_params {
    int lsp_index[3]; /* 8, 9 and 9 bits */
    tetra_subframe_params sub[TETRA_SUBFRAMES];
} tetra_frame_params;

/** Per-stream decoder state. */
typedef struct tetra_acelp_state {
    int16_t lsp_old[TETRA_M];                         /* previous frame's LSPs, Q15 cosine domain */
    int16_t exc[TETRA_EXC_HIST + TETRA_SUBFRAME];     /* past excitation, newest last */
    int16_t syn_mem[TETRA_M];                         /* synthesis filter memory, oldest first */
    int16_t energy_p;                                 /* quantized pitch energy, Q8 log2 */
    int16_t energy_c;                                 /* quantized code energy, Q8 log2 */
    int16_t last_pitch;                               /* integer lag of the last good subframe 4 */
    int good_seen;                                    /* a good frame has been decoded since the reset */
    tetra_subframe_params last_code[TETRA_SUBFRAMES]; /* previous frame's codebook indices */
} tetra_acelp_state;

void tetra_acelp_reset(tetra_acelp_state* st);
/** Unpack B1..B137 (one 0/1 byte each). */
void tetra_acelp_unpack(const uint8_t* bits, tetra_frame_params* fp);
/** Decode one frame into TETRA_FRAME samples; a NULL fp is a bad frame, which is concealed. */
void tetra_acelp_decode(tetra_acelp_state* st, const tetra_frame_params* fp, int16_t* pcm);

/* Stages, exposed for tests. */
void tetra_lsp_decode(const int lsp_index[3], int16_t lsp[TETRA_M]);
/** Nonzero when the cosines strictly decrease, i.e. the LSPs give a stable synthesis filter. */
int tetra_lsp_ordered(const int16_t lsp[TETRA_M]);
void tetra_lsp_interpolate(const int16_t old_lsp[TETRA_M], const int16_t new_lsp[TETRA_M],
                           int16_t out[TETRA_SUBFRAMES][TETRA_M]);
void tetra_lsp_to_lpc(const int16_t lsp[TETRA_M], int16_t a[TETRA_M + 1]);
void tetra_decode_pitch(int subframe, int index, int t1, int* t0, int* frac);
void tetra_pulse_positions(const tetra_subframe_params* sp, int pos[4]);
void tetra_syn_filt(const int16_t a[TETRA_M + 1], const int16_t* x, int16_t* y, int len, int16_t mem[TETRA_M],
                    int update);
void tetra_codevector(const int16_t a[TETRA_M + 1], int t0, const tetra_subframe_params* sp,
                      int16_t code[TETRA_SUBFRAME]);
/** Write the adaptive codebook vector over exc_cur[0..59]; exc_cur[-TETRA_EXC_HIST..-1] is the past
 * excitation. Samples within the subframe are repeated for lags below the subframe length. */
void tetra_adaptive_vector(int16_t* exc_cur, int t0, int frac);
/** AVX2 part of tetra_adaptive_vector() for a fractional lag (tetra_acelp_avx2.c): writes the leading
 * outputs eight at a time and returns how many, 0 when the lag is too short for it. Built only for x86-64
 * with AVX2, either by the target or for the run-time dispatch. */
int tetra_adaptive_avx2(int16_t* exc_cur, int t0, int first, const int16_t* taps);
void tetra_log2(int32_t x, int16_t* exponent, int16_t* fraction);
int32_t tetra_pow2(int16_t exponent, int16_t fraction);
/** Pitch and code gains (Q12, the pitch gain limited to 1.2, the code gain kept to its low 16 bits) from
 * the quantized energies (Q8) and the energies of the adaptive vector plus one, the code vector and the
 * synthesis filter's response to an impulse of 1024 (clause 4.2.2.6, equations 31-39). */
void tetra_gains(int16_t energy_p, int16_t energy_c, int32_t v_energy, int32_t c_energy, int32_t pg_energy, int16_t* gp,
                 int16_t* gc);
/** One excitation sample u = gp v + gc c (equation 40) from Q12 gains, rounded, kept to its low 16 bits. */
int16_t tetra_excitation(int16_t gp, int16_t v, int16_t gc, int16_t c);
/** Mix one subframe in place; gp is in 0..4915, gc may span int16. */
void tetra_mix_excitation(int16_t gp, int16_t gc, int16_t v[TETRA_SUBFRAME], const int16_t c[TETRA_SUBFRAME]);
/** Saturating sum of squares plus init (0 or 1), for 0..TETRA_SUBFRAME samples. */
int32_t tetra_energy(const int16_t* x, int len, int32_t init);

#endif /* MBELIB_NEO_INTERNAL_TETRA_ACELP_H */
