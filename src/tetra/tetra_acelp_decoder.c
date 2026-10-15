// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief TETRA full-rate ACELP speech decoder (EN 300 395-2 clause 4.2.3) and its public API.
 */

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mbe_fixed_point.h"
#include "mbe_result.h"
#include "mbelib-neo/mbelib.h"
#include "tetra_acelp.h"
#include "tetra_acelp_codebooks.h"

/* Initial past quantized LSP vector, EN 300 395-2 clause 4.2.2.3 (Q15 cosine domain). */
static const int16_t tetra_lsp_init[TETRA_M] = {30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000};

#define TETRA_PRED_CONST    (-768) /* -3.0 in the Q8 log2 energy domain, equation 35 */
#define TETRA_ENERGY_MAX_P  6912   /* 27.0, clause 4.2.2.6 */
#define TETRA_ENERGY_MAX_C  6400   /* 25.0 */
#define TETRA_CONCEAL_DECAY 128    /* 0.5, equation 42 */
#define TETRA_PG_IMPULSE    1024   /* impulse amplitude of the prediction-gain response */

void
tetra_acelp_reset(tetra_acelp_state* st) {
    memset(st, 0, sizeof(*st));
    memcpy(st->lsp_old, tetra_lsp_init, sizeof(st->lsp_old));
    st->last_pitch = TETRA_PIT_MIN;
}

static int
read_bits(const uint8_t* bits, int* pos, int count) {
    int value = 0;
    while (count-- > 0) {
        value = (value << 1) | bits[(*pos)++];
    }
    return value;
}

void
tetra_acelp_unpack(const uint8_t* bits, tetra_frame_params* fp) {
    int pos = 0;
    int k;
    fp->lsp_index[0] = read_bits(bits, &pos, 8);
    fp->lsp_index[1] = read_bits(bits, &pos, 9);
    fp->lsp_index[2] = read_bits(bits, &pos, 9);
    for (k = 0; k < TETRA_SUBFRAMES; k++) {
        tetra_subframe_params* sp = &fp->sub[k];
        sp->pitch = read_bits(bits, &pos, k == 0 ? 8 : 5);
        sp->pulse4 = read_bits(bits, &pos, 3);
        sp->pulse3 = read_bits(bits, &pos, 3);
        sp->pulse2 = read_bits(bits, &pos, 3);
        sp->pulse1 = read_bits(bits, &pos, 5);
        sp->sign = read_bits(bits, &pos, 1);
        sp->shift = read_bits(bits, &pos, 1);
        sp->gain = read_bits(bits, &pos, 6);
    }
}

/** Energy of x plus init, accumulated with saturation. */
static int32_t
energy(const int16_t* x, int len, int32_t init) {
    int32_t acc = init;
    int n;
    for (n = 0; n < len; n++) {
        acc = mbe_fx_L_mac0(acc, x[n], x[n]);
    }
    return acc;
}

/** Predicted energy (equation 35), floored at zero. */
static int32_t
predict_energy(int32_t own, int32_t other) {
    int32_t e = (int32_t)mbe_fx_asr64(2 * (int64_t)own + other, 2) + TETRA_PRED_CONST;
    return e > 0 ? e : 0;
}

/** Quantized pitch and code energies of one subframe (Q8). */
static void
subframe_energies(tetra_acelp_state* st, int gain_index, int bad_frame) {
    int32_t ep;
    int32_t ec;
    if (bad_frame) {
        /* Equation 42, without letting the energies go negative. */
        ep = (int32_t)st->energy_p - TETRA_CONCEAL_DECAY;
        ec = (int32_t)st->energy_c - TETRA_CONCEAL_DECAY;
        ep = ep > 0 ? ep : 0;
        ec = ec > 0 ? ec : 0;
    } else {
        ep = predict_energy(st->energy_p, st->energy_c) + tetra_gain_cb[gain_index][0];
        ec = predict_energy(st->energy_c, st->energy_p) + tetra_gain_cb[gain_index][1];
    }
    st->energy_p = mbe_fx_sat16(ep < TETRA_ENERGY_MAX_P ? ep : TETRA_ENERGY_MAX_P);
    st->energy_c = mbe_fx_sat16(ec < TETRA_ENERGY_MAX_C ? ec : TETRA_ENERGY_MAX_C);
}

/** Decode one subframe into the excitation buffer and synthesize it. */
static void
decode_subframe(tetra_acelp_state* st, const int16_t a[TETRA_M + 1], const tetra_subframe_params* sp, int t0, int frac,
                int bad_frame, int16_t* out) {
    int16_t* exc_cur = &st->exc[TETRA_EXC_HIST];
    int16_t code[TETRA_SUBFRAME];
    int16_t h[TETRA_SUBFRAME];
    int16_t x[TETRA_SUBFRAME];
    int16_t zero_mem[TETRA_M];
    int16_t s[TETRA_SUBFRAME];
    int32_t pg;
    int16_t gp;
    int16_t gc;
    int n;

    tetra_adaptive_vector(exc_cur, t0, frac);
    tetra_codevector(a, t0, sp, code);
    subframe_energies(st, sp->gain, bad_frame);

    memset(x, 0, sizeof(x));
    x[0] = TETRA_PG_IMPULSE;
    memset(zero_mem, 0, sizeof(zero_mem));
    tetra_syn_filt(a, x, h, TETRA_SUBFRAME, zero_mem, 0);
    pg = energy(h, TETRA_SUBFRAME, 0);
    tetra_gains(st->energy_p, st->energy_c, energy(exc_cur, TETRA_SUBFRAME, 1), energy(code, TETRA_SUBFRAME, 0), pg,
                &gp, &gc);

    for (n = 0; n < TETRA_SUBFRAME; n++) {
        exc_cur[n] = tetra_excitation(gp, exc_cur[n], gc, code[n]);
    }
    tetra_syn_filt(a, exc_cur, s, TETRA_SUBFRAME, st->syn_mem, 1);
    for (n = 0; n < TETRA_SUBFRAME; n++) {
        out[n] = mbe_fx_add(s[n], s[n]); /* post-processing: x2 with saturation */
    }
    memmove(st->exc, &st->exc[TETRA_SUBFRAME], TETRA_EXC_HIST * sizeof(int16_t));
}

void
tetra_acelp_decode(tetra_acelp_state* st, const tetra_frame_params* fp, int16_t* pcm) {
    int16_t lsp_new[TETRA_M];
    int16_t lsp_sub[TETRA_SUBFRAMES][TETRA_M];
    int16_t a[TETRA_M + 1];
    int t1 = 0;
    int k;

    if (!fp && !st->good_seen) {
        /* A bad frame with nothing to conceal from yet: stay silent and keep the reset state,
         * except that the first previous-frame LSP becomes 0, as in the reference decoder. The
         * first good frame then interpolates from that vector (docs/tetra-acelp.md). */
        st->lsp_old[0] = 0;
        memset(pcm, 0, TETRA_FRAME * sizeof(int16_t));
        return;
    }
    memcpy(lsp_new, st->lsp_old, sizeof(lsp_new)); /* a bad frame keeps the last good LSPs */
    if (fp) {
        int16_t decoded[TETRA_M];
        tetra_lsp_decode(fp->lsp_index, decoded);
        /* LSPs still crossed after the boundary widening are discarded the same way. */
        if (tetra_lsp_ordered(decoded)) {
            memcpy(lsp_new, decoded, sizeof(lsp_new));
        }
    }
    tetra_lsp_interpolate(st->lsp_old, lsp_new, lsp_sub);
    for (k = 0; k < TETRA_SUBFRAMES; k++) {
        const tetra_subframe_params* sp = fp ? &fp->sub[k] : &st->last_code[k];
        int t0 = st->last_pitch; /* a bad frame repeats the last good subframe-4 lag, fraction zero */
        int frac = 0;
        if (fp) {
            tetra_decode_pitch(k, sp->pitch, t1, &t0, &frac);
            if (k == 0) {
                t1 = t0;
            }
        }
        tetra_lsp_to_lpc(lsp_sub[k], a);
        decode_subframe(st, a, sp, t0, frac, fp == NULL, &pcm[(ptrdiff_t)k * TETRA_SUBFRAME]);
        if (fp && k == TETRA_SUBFRAMES - 1) {
            st->last_pitch = (int16_t)t0;
        }
    }
    memcpy(st->lsp_old, lsp_new, sizeof(st->lsp_old));
    if (fp) {
        memcpy(st->last_code, fp->sub, sizeof(st->last_code));
        st->good_seen = 1;
    }
}

/* === Public API === */

struct mbe_tetra_decoder {
    tetra_acelp_state state;
};

mbe_tetra_decoder*
mbe_tetraDecoderAlloc(void) {
    mbe_tetra_decoder* dec = (mbe_tetra_decoder*)malloc(sizeof(*dec));
    if (dec) {
        tetra_acelp_reset(&dec->state);
    }
    return dec;
}

void
mbe_tetraDecoderReset(mbe_tetra_decoder* dec) {
    if (dec) {
        tetra_acelp_reset(&dec->state);
    }
}

void
mbe_tetraDecoderFree(mbe_tetra_decoder* dec) {
    free(dec);
}

int
mbe_tetraDecodeFrame(mbe_tetra_decoder* dec, const uint8_t* bits, int bad_frame, int16_t* pcm) {
    tetra_frame_params fp;
    if (!dec || !pcm) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    if (!bad_frame) {
        /* Checked before anything is touched, so a rejected frame can be retried as bad. */
        int status = mbe_validate_bits_u8(bits, MBE_TETRA_FRAME_BITS);
        if (status != 0) {
            return status;
        }
        tetra_acelp_unpack(bits, &fp);
    }
    tetra_acelp_decode(&dec->state, bad_frame ? NULL : &fp, pcm);
    return 0;
}
