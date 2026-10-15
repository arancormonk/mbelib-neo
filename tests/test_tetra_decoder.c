// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief TETRA ACELP decoder: decoding stages, concealment, the public API contract, and
 *        regression hashes of whole streams.
 */

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"
#include "tetra_acelp.h"
#include "tetra_acelp_codebooks.h"

/* Same initial LSPs as the decoder (EN 300 395-2 clause 4.2.2.3). */
static const int16_t k_lsp_init[TETRA_M] = {30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000};

/* === Helpers === */

typedef struct lcg {
    uint32_t state;
} lcg;

static uint32_t
lcg_next(lcg* g) {
    g->state = g->state * 1664525u + 1013904223u;
    return g->state >> 8;
}

static void
put_field(uint8_t* bits, int* pos, int count, int value) {
    int i;
    for (i = count - 1; i >= 0; i--) {
        bits[(*pos)++] = (uint8_t)((value >> i) & 1);
    }
}

/** Pack frame parameters in the table 3 order (the inverse of tetra_acelp_unpack). */
static void
pack_frame(const tetra_frame_params* fp, uint8_t bits[MBE_TETRA_FRAME_BITS]) {
    int pos = 0;
    int k;
    put_field(bits, &pos, 8, fp->lsp_index[0]);
    put_field(bits, &pos, 9, fp->lsp_index[1]);
    put_field(bits, &pos, 9, fp->lsp_index[2]);
    for (k = 0; k < TETRA_SUBFRAMES; k++) {
        const tetra_subframe_params* sp = &fp->sub[k];
        put_field(bits, &pos, k == 0 ? 8 : 5, sp->pitch);
        put_field(bits, &pos, 3, sp->pulse4);
        put_field(bits, &pos, 3, sp->pulse3);
        put_field(bits, &pos, 3, sp->pulse2);
        put_field(bits, &pos, 5, sp->pulse1);
        put_field(bits, &pos, 1, sp->sign);
        put_field(bits, &pos, 1, sp->shift);
        put_field(bits, &pos, 6, sp->gain);
    }
    assert(pos == MBE_TETRA_FRAME_BITS);
}

static void
random_bits(lcg* g, uint8_t bits[MBE_TETRA_FRAME_BITS]) {
    int i;
    for (i = 0; i < MBE_TETRA_FRAME_BITS; i++) {
        bits[i] = (uint8_t)((lcg_next(g) >> 7) & 1u);
    }
}

/** A plausible frame: ordered LSPs, moderate gains, every pulse inside the subframe. */
static void
voiced_frame(lcg* g, uint8_t bits[MBE_TETRA_FRAME_BITS]) {
    tetra_frame_params fp;
    int k;
    fp.lsp_index[0] = 40 + (int)(lcg_next(g) % 40u);
    fp.lsp_index[1] = 60 + (int)(lcg_next(g) % 60u);
    fp.lsp_index[2] = 60 + (int)(lcg_next(g) % 60u);
    for (k = 0; k < TETRA_SUBFRAMES; k++) {
        tetra_subframe_params* sp = &fp.sub[k];
        sp->pitch = (k == 0) ? 60 + (int)(lcg_next(g) % 120u) : (int)(lcg_next(g) % 32u);
        sp->pulse1 = (int)(lcg_next(g) % 30u);
        sp->pulse2 = (int)(lcg_next(g) % 7u);
        sp->pulse3 = (int)(lcg_next(g) % 7u);
        sp->pulse4 = (int)(lcg_next(g) % 7u);
        sp->sign = (int)(lcg_next(g) & 1u);
        sp->shift = (int)(lcg_next(g) & 1u);
        sp->gain = 20 + (int)(lcg_next(g) % 40u);
    }
    pack_frame(&fp, bits);
}

static uint32_t
fnv1a_pcm(uint32_t h, const int16_t* pcm, size_t count) {
    size_t i;
    for (i = 0; i < count; i++) {
        uint16_t v = (uint16_t)pcm[i];
        h = (h ^ (v & 0xFFu)) * 16777619u;
        h = (h ^ (v >> 8)) * 16777619u;
    }
    return h;
}

static int
pcm_is_silent(const int16_t* pcm) {
    int n;
    for (n = 0; n < MBE_TETRA_FRAME_SAMPLES; n++) {
        if (pcm[n] != 0) {
            return 0;
        }
    }
    return 1;
}

/* === Bit unpacking (tables 3 and 4) === */

static void
test_unpack_round_trip(void) {
    lcg g = {7u};
    int trial;
    for (trial = 0; trial < 200; trial++) {
        uint8_t bits[MBE_TETRA_FRAME_BITS];
        uint8_t again[MBE_TETRA_FRAME_BITS];
        tetra_frame_params fp;
        random_bits(&g, bits);
        tetra_acelp_unpack(bits, &fp);
        pack_frame(&fp, again);
        assert(memcmp(bits, again, sizeof(bits)) == 0);
    }
}

static void
test_unpack_field_positions(void) {
    uint8_t bits[MBE_TETRA_FRAME_BITS];
    tetra_frame_params fp;
    int i;

    /* B1 is the MSB of the 8-bit LSP1-3 index; B18..B26 hold LSP7-10. */
    memset(bits, 0, sizeof(bits));
    bits[0] = 1;
    bits[25] = 1;
    tetra_acelp_unpack(bits, &fp);
    assert(fp.lsp_index[0] == 128 && fp.lsp_index[1] == 0 && fp.lsp_index[2] == 1);

    /* Subframe 3: pulse 2 is B95-B97 and pulse 1 is B98-B102 (table 4; table 3 misprints them). */
    memset(bits, 0, sizeof(bits));
    for (i = 98; i <= 102; i++) {
        bits[i - 1] = 1;
    }
    tetra_acelp_unpack(bits, &fp);
    assert(fp.sub[2].pulse1 == 31 && fp.sub[2].pulse2 == 0);
    memset(bits, 0, sizeof(bits));
    for (i = 95; i <= 97; i++) {
        bits[i - 1] = 1;
    }
    tetra_acelp_unpack(bits, &fp);
    assert(fp.sub[2].pulse2 == 7 && fp.sub[2].pulse1 == 0 && fp.sub[2].pulse3 == 0);

    /* Subframe 4's gain index is the last field, B132-B137. */
    memset(bits, 0, sizeof(bits));
    bits[131] = 1;
    bits[136] = 1;
    tetra_acelp_unpack(bits, &fp);
    assert(fp.sub[3].gain == 33);
}

/* === Pitch === */

static void
test_pitch_first_subframe(void) {
    int index;
    int prev = 0;
    for (index = 0; index < 256; index++) {
        int t0;
        int frac;
        int thirds;
        tetra_decode_pitch(0, index, 0, &t0, &frac);
        assert(frac >= -1 && frac <= 1);
        thirds = 3 * t0 + frac;
        assert(index == 0 || thirds > prev);
        if (index >= 197) {
            assert(frac == 0);
        }
        prev = thirds;
    }
    {
        int t0;
        int frac;
        tetra_decode_pitch(0, 0, 0, &t0, &frac);
        assert(3 * t0 + frac == 58); /* 19 1/3 */
        tetra_decode_pitch(0, 196, 0, &t0, &frac);
        assert(3 * t0 + frac == 254); /* 84 2/3 */
        tetra_decode_pitch(0, 197, 0, &t0, &frac);
        assert(t0 == 85 && frac == 0);
        tetra_decode_pitch(0, 255, 0, &t0, &frac);
        assert(t0 == 143 && frac == 0);
    }
}

static void
test_pitch_relative(void) {
    int t1;
    for (t1 = 19; t1 <= 143; t1++) {
        int tmin = t1 - 5 < TETRA_PIT_MIN ? TETRA_PIT_MIN : t1 - 5;
        int index;
        int prev = 0;
        if (tmin + 9 > TETRA_PIT_MAX) {
            tmin = TETRA_PIT_MAX - 9;
        }
        for (index = 0; index < 32; index++) {
            int t0;
            int frac;
            int thirds;
            tetra_decode_pitch(1, index, t1, &t0, &frac);
            thirds = 3 * t0 + frac;
            assert(frac >= -1 && frac <= 1);
            assert(index == 0 || thirds == prev + 1);
            /* 32 thirds from tmin - 2/3 to tmin + 9 + 2/3. */
            assert(thirds >= 3 * tmin - 2 && thirds <= 3 * (tmin + 9) + 2);
            prev = thirds;
        }
    }
}

/* === Algebraic codevector === */

static void
test_pulse_positions(void) {
    tetra_subframe_params sp;
    int pos[4];
    memset(&sp, 0, sizeof(sp));
    sp.pulse1 = 29;
    sp.pulse2 = 7;
    sp.pulse3 = 6;
    sp.pulse4 = 0;
    sp.shift = 1;
    tetra_pulse_positions(&sp, pos);
    assert(pos[0] == 59 && pos[1] == 59 && pos[2] == 53 && pos[3] == 7);
    sp.pulse1 = 31;
    sp.pulse3 = 7;
    sp.pulse4 = 7;
    tetra_pulse_positions(&sp, pos);
    assert(pos[0] == 63 && pos[2] == 61 && pos[3] == 63);
}

static void
test_codevector_absent_pulses(void) {
    int16_t a[TETRA_M + 1];
    int16_t code[TETRA_SUBFRAME];
    int16_t neg[TETRA_SUBFRAME];
    tetra_subframe_params sp;
    int n;
    tetra_lsp_to_lpc(k_lsp_init, a);

    /* Pulses at 60 and beyond are outside the subframe; pulse 2 is always inside. */
    memset(&sp, 0, sizeof(sp));
    sp.pulse1 = 30;
    sp.pulse2 = 7;
    sp.pulse3 = 7;
    sp.pulse4 = 7;
    sp.shift = 1; /* pulse 2 at 59 */
    tetra_codevector(a, 143, &sp, code);
    for (n = 0; n < 59; n++) {
        assert(code[n] == 0);
    }
    assert(code[59] == -4096); /* pulse 2 has amplitude -1, and F(z) starts at 1.0 (Q12) */

    sp.sign = 1;
    tetra_codevector(a, 143, &sp, neg);
    assert(neg[59] == 4096);

    /* Pulse 1 at 0 with amplitude 1.4142 shapes the whole subframe. */
    sp.sign = 0;
    sp.pulse1 = 0;
    sp.shift = 0;
    sp.pulse2 = 7; /* 58 */
    tetra_codevector(a, 143, &sp, code);
    assert(code[0] == 5792);
}

static void
test_codevector_sharpening(void) {
    int16_t a[TETRA_M + 1];
    int16_t plain[TETRA_SUBFRAME];
    int16_t sharp[TETRA_SUBFRAME];
    tetra_subframe_params sp;
    int n;
    tetra_lsp_to_lpc(k_lsp_init, a);
    memset(&sp, 0, sizeof(sp));
    sp.pulse2 = 7;
    sp.pulse3 = 7;
    sp.pulse4 = 7; /* only pulse 1 (at 0) and pulse 2 (at 58) inside */
    tetra_codevector(a, 60, &sp, plain);
    tetra_codevector(a, 30, &sp, sharp);
    for (n = 0; n < 30; n++) {
        assert(plain[n] == sharp[n]); /* f(n) += 0.8 f(n - 30) only from n = 30 */
    }
    /* At n = 30 the echo of the pulse at 0 adds about 0.8 * 5792. */
    assert(abs((sharp[30] - plain[30]) - 4634) < 64);
}

static void
test_codevector_wraps(void) {
    /* A(z) of a subframe from random bits whose shaped pulses add up past 16 bits. Like the
     * reference decoder, the codevector keeps the low 16 bits of each sum (docs/tetra-acelp.md). */
    static const int16_t a[TETRA_M + 1] = {4096, -15148, 30060, 28874, 28014, -8564, -7298, 11777, -7368, 2341, -146};
    int16_t code[TETRA_SUBFRAME];
    tetra_subframe_params sp;
    memset(&sp, 0, sizeof(sp));
    sp.pulse1 = 8;
    sp.pulse2 = 4;
    sp.pulse3 = 3;
    sp.pulse4 = 6;
    sp.sign = 1;
    sp.shift = 1;
    tetra_codevector(a, 33, &sp, code);
    assert(code[21] == -19200); /* 46336 */
    assert(code[23] == 19201);  /* -46335 */
    assert(code[33] == 13568);  /* 79104 */
    assert(code[38] == -165);   /* 65371 */
}

/* === LSPs and the LP filter === */

static void
test_lsp_decode_and_widening(void) {
    int idx[3];
    int16_t lsp[TETRA_M];
    int i;
    int i1;
    int widened = 0;
    idx[0] = 0;
    idx[1] = 0;
    idx[2] = 0;
    tetra_lsp_decode(idx, lsp);
    for (i = 0; i < 2; i++) {
        assert(lsp[i] == tetra_lsp_cb1[0][i]);
    }
    for (i = 7; i < TETRA_M; i++) {
        assert(lsp[i] == tetra_lsp_cb3[0][i - 6]);
    }
    /* Every first sub-vector against a fixed second one: the q3/q4 boundary ends up >= 915 apart. */
    idx[1] = 300;
    idx[2] = 300;
    for (i1 = 0; i1 < 256; i1++) {
        int raw_gap = tetra_lsp_cb1[i1][2] - tetra_lsp_cb2[300][0];
        idx[0] = i1;
        tetra_lsp_decode(idx, lsp);
        if (raw_gap < 916) {
            widened++;
            assert(lsp[2] - lsp[3] >= 915 && lsp[2] - lsp[3] <= 916);
            /* Both sides move by half the shortfall. */
            assert(lsp[2] - tetra_lsp_cb1[i1][2] == (916 - raw_gap) / 2);
        } else {
            assert(lsp[2] == tetra_lsp_cb1[i1][2] && lsp[3] == tetra_lsp_cb2[300][0]);
        }
    }
    assert(widened > 0);
}

static void
test_lsp_ordering(void) {
    int16_t lsp[TETRA_M];
    memcpy(lsp, k_lsp_init, sizeof(lsp));
    assert(tetra_lsp_ordered(lsp));
    lsp[4] = lsp[5] + 1;
    assert(tetra_lsp_ordered(lsp));
    lsp[4] = lsp[5]; /* equal cosines are rejected */
    assert(!tetra_lsp_ordered(lsp));
    lsp[4] = lsp[5] - 1;
    assert(!tetra_lsp_ordered(lsp));
}

static void
test_lsp_interpolation(void) {
    int16_t new_lsp[TETRA_M];
    int16_t sub[TETRA_SUBFRAMES][TETRA_M];
    int i;
    for (i = 0; i < TETRA_M; i++) {
        new_lsp[i] = (int16_t)(k_lsp_init[i] - 400 * (i + 1));
    }
    tetra_lsp_interpolate(k_lsp_init, new_lsp, sub);
    for (i = 0; i < TETRA_M; i++) {
        /* 3/4, 1/2 and 1/4 of the old vector, rounded toward minus infinity. */
        assert(sub[0][i] == (int16_t)floor(0.75 * k_lsp_init[i] + 0.25 * new_lsp[i]));
        assert(sub[1][i] == (int16_t)floor(0.5 * k_lsp_init[i] + 0.5 * new_lsp[i]));
        assert(sub[2][i] == (int16_t)floor(0.25 * k_lsp_init[i] + 0.75 * new_lsp[i]));
        assert(sub[3][i] == new_lsp[i]);
    }
}

/** A(z) from the LSPs in double precision (equations 17-18), Q12. */
static void
lsp_to_lpc_double(const int16_t lsp[TETRA_M], double a[TETRA_M + 1]) {
    double f1[12] = {1.0};
    double f2[12] = {1.0};
    int i;
    int j;
    for (i = 0; i < 5; i++) {
        double c1 = -2.0 * lsp[(ptrdiff_t)2 * i] / 32768.0;
        double c2 = -2.0 * lsp[(ptrdiff_t)2 * i + 1] / 32768.0;
        for (j = 2 * i + 2; j >= 2; j--) {
            f1[j] += c1 * f1[j - 1] + f1[j - 2];
            f2[j] += c2 * f2[j - 1] + f2[j - 2];
        }
        f1[1] += c1;
        f2[1] += c2;
    }
    for (i = 1; i <= 5; i++) {
        double p = f1[i] + f1[i - 1];
        double q = f2[i] - f2[i - 1];
        a[i] = 0.5 * (p + q) * 4096.0;
        a[TETRA_M + 1 - i] = 0.5 * (p - q) * 4096.0;
    }
    a[0] = 4096.0;
}

static void
test_lsp_to_lpc(void) {
    lcg g = {3u};
    int trial;
    for (trial = 0; trial < 300; trial++) {
        int idx[3];
        int16_t lsp[TETRA_M];
        int16_t a[TETRA_M + 1];
        double ref[TETRA_M + 1];
        int i;
        idx[0] = (int)(lcg_next(&g) % 256u);
        idx[1] = (int)(lcg_next(&g) % 512u);
        idx[2] = (int)(lcg_next(&g) % 512u);
        tetra_lsp_decode(idx, lsp);
        if (!tetra_lsp_ordered(lsp)) {
            continue;
        }
        tetra_lsp_to_lpc(lsp, a);
        lsp_to_lpc_double(lsp, ref);
        assert(a[0] == 4096);
        for (i = 1; i <= TETRA_M; i++) {
            /* A coefficient beyond +-8 does not fit Q12 in 16 bits and wraps; skip those. */
            if (fabs(ref[i]) < 32767.0) {
                assert(fabs(a[i] - ref[i]) <= 1.0);
            }
        }
    }
}

/* Interpolated LSP vectors from random codebook transitions, and the Q12 coefficients ETSI's
 * reference decoder uses for them, read back from its output (a single pulse per subframe, so the
 * output decays freely through 1/A(z)) and proven unique over the whole 16-bit range. Each vector has
 * odd LSPs, and in each the usual double-precision product (2 hi q + 2 floor(lo q / 2^15)) gives one
 * coefficient one away from these. */
static void
test_lsp_to_lpc_reference(void) {
    static const struct {
        int16_t lsp[TETRA_M];
        int16_t a[TETRA_M];
    } k_cases[] = {
        {{30292, 26654, 26057, 25141, 18702, 8508, -3942, -11835, -21330, -30701},
         {-8443, 8389, -4195, 1799, 529, -2030, 2141, -882, 113, 95}},
        {{30560, 29038, 25231, 24315, 15228, 11942, -17880, -22825, -26879, -30264},
         {-4808, -426, 4168, 1474, -5680, 2416, 3136, -1556, -2100, 2339}},
        {{30496, 27960, 19049, 18133, 1680, -9682, -15246, -18143, -24760, -28796},
         {-86, 1484, -1775, 2878, -99, 512, -1404, 1487, 816, 1378}},
        {{31778, 30872, 25132, 19156, 4453, -2313, -6228, -11296, -24776, -30248},
         {-4566, 3258, -3317, 861, 1179, -72, -919, 1064, -990, 1073}},
    };

    size_t c;
    for (c = 0; c < sizeof(k_cases) / sizeof(k_cases[0]); c++) {
        int16_t a[TETRA_M + 1];
        int i;
        tetra_lsp_to_lpc(k_cases[c].lsp, a);
        assert(a[0] == 4096);
        for (i = 1; i <= TETRA_M; i++) {
            assert(a[i] == k_cases[c].a[i - 1]);
        }
    }
}

static void
test_synthesis_filter(void) {
    int16_t a[TETRA_M + 1] = {4096, -2048, 0, 0, 0, 0, 0, 0, 0, 0, 0}; /* 1 / (1 - 0.5 z^-1) */
    int16_t x[TETRA_SUBFRAME];
    int16_t y[TETRA_SUBFRAME];
    int16_t mem[TETRA_M];
    int16_t mem_before[TETRA_M];
    memset(x, 0, sizeof(x));
    memset(mem, 0, sizeof(mem));
    x[0] = 1000;
    tetra_syn_filt(a, x, y, TETRA_SUBFRAME, mem, 1);
    assert(y[0] == 1000 && y[1] == 500 && y[2] == 250 && y[3] == 125 && y[4] == 63 && y[5] == 32);
    assert(mem[TETRA_M - 1] == y[TETRA_SUBFRAME - 1]);

    /* update = 0 leaves the memory alone; the output saturates instead of wrapping. */
    memcpy(mem_before, mem, sizeof(mem));
    x[0] = 32767;
    x[1] = 32767;
    tetra_syn_filt(a, x, y, 2, mem, 0);
    assert(memcmp(mem, mem_before, sizeof(mem)) == 0);
    assert(y[1] == 32767);
}

/* === Adaptive codebook === */

static void
test_adaptive_vector(void) {
    int16_t buf[TETRA_EXC_HIST + TETRA_SUBFRAME];
    int16_t* cur = &buf[TETRA_EXC_HIST];
    int n;

    /* Integer lag below the subframe length repeats the last T0 samples. */
    for (n = 0; n < TETRA_EXC_HIST; n++) {
        buf[n] = (int16_t)(n * 7 - 500);
    }
    tetra_adaptive_vector(cur, 25, 0);
    for (n = 0; n < TETRA_SUBFRAME; n++) {
        assert(cur[n] == cur[n - 25]);
    }

    /* A fractional lag of a constant past gives (almost) the constant back. */
    for (n = 0; n < TETRA_EXC_HIST; n++) {
        buf[n] = 1000;
    }
    tetra_adaptive_vector(cur, 143, 1);
    for (n = 0; n < TETRA_SUBFRAME; n++) {
        assert(abs(cur[n] - 1000) <= 20);
    }
    tetra_adaptive_vector(cur, 144, -1);
    for (n = 0; n < TETRA_SUBFRAME; n++) {
        assert(abs(cur[n] - 1000) <= 20);
    }
}

/* === Gains === */

static void
test_log2_pow2(void) {
    int32_t x;
    int16_t e;
    int16_t f;
    for (x = 1; x > 0 && x < INT32_MAX / 3; x = x * 3 + 1) {
        double exact = log2((double)x);
        tetra_log2(x, &e, &f);
        assert(fabs(e + f / 32768.0 - exact) < 1.0 / 2048);
    }
    tetra_log2(0, &e, &f);
    assert(e == 0 && f == 0);
    tetra_log2(-5, &e, &f);
    assert(e == 0 && f == 0);
    for (e = 0; e <= 30; e += 5) {
        for (f = 0; f < 32767 - 4096; f = (int16_t)(f + 4096)) {
            double exact = pow(2.0, e + f / 32768.0);
            int32_t got = tetra_pow2(e, f);
            assert(fabs(got - exact) <= 1.0 + exact * 1e-4);
        }
    }
}

/* Gains from the quantized energies (Q8) and the energies of the adaptive vector (plus one), the code
 * vector and the synthesis filter's impulse response. The first seven pairs are what ETSI's reference
 * decoder produced on subframes of this file's voiced and random streams whose excitation determines
 * the pair uniquely; each rules out at least one other plausible computation (separate logs of the
 * energies, a log table scaled by 32768, a normalized code energy, the log rounded instead of
 * truncated). The last three were worked out from the computation in docs/tetra-acelp.md: a negative
 * log (an adaptive vector of almost no energy), the smallest possible code energy (see
 * test_code_energy_floor) and saturated energies. */
static void
test_gains(void) {
    static const struct {
        int16_t ep;
        int16_t ec;
        int32_t v;
        int32_t c;
        int32_t pg;
        int16_t gp;
        int16_t gc;
    } k_cases[] = {
        {6912, 6132, 203887612, 95280214, 5305354, 4915, 4760}, /* pitch gain at its 1.2 limit */
        {4466, 5994, 32412758, 78162297, 3285482, 1087, 5540},
        {5315, 4373, 35456188, 88262392, 5830866, 2462, 436},
        {3716, 4344, 659049, 87242289, 6180387, 2012, 410},
        {968, 4248, 22708, 47961888, 12466714, 185, 342},
        {5338, 3827, 23679096, 101696287, 3301618, 4130, 258},
        {2671, 3857, 32060785, 21733092, 2747147, 105, 637},
        {0, 3000, 10, 1 << 24, 3000000, 4845, 217},
        {0, 2600, 1, 2876416, 1 << 20, 4915, 521},
        {6912, 6400, INT32_MAX, INT32_MAX, INT32_MAX, 143, 72},
        /* A code gain past 16 bits keeps its low 16 bits (34265 and 68436), as in the reference
         * decoder; the pitch gain is capped instead (51783 here, which wrapped would be negative). */
        {5764, 6400, 165943905, 3067269, 6680446, 1951, -31271},
        {4525, 6400, 448152011, 2929961, 1748851, 433, 2900},
        {512, 2600, 1, 2876416, 1 << 20, 4915, 521},
    };

    size_t i;
    for (i = 0; i < sizeof(k_cases) / sizeof(k_cases[0]); i++) {
        int16_t gp = -1;
        int16_t gc = -1;
        tetra_gains(k_cases[i].ep, k_cases[i].ec, k_cases[i].v, k_cases[i].c, k_cases[i].pg, &gp, &gc);
        assert(gp == k_cases[i].gp);
        assert(gc == k_cases[i].gc);
    }
}

/* The code energy's high half is used without normalization, so it must never be zero: the earliest
 * pulse sample alone is at least (1.4142 - 1) x 4096 = 1696 in magnitude, whatever the pulses. */
static void
test_code_energy_floor(void) {
    int16_t a[TETRA_M + 1];
    int64_t least = INT64_MAX;
    int combo;
    tetra_lsp_to_lpc(k_lsp_init, a);
    for (combo = 0; combo < (1 << 16); combo++) {
        tetra_subframe_params sp;
        int16_t code[TETRA_SUBFRAME];
        int64_t e = 0;
        int n;
        memset(&sp, 0, sizeof(sp));
        sp.pulse1 = combo & 31;
        sp.pulse2 = (combo >> 5) & 7;
        sp.pulse3 = (combo >> 8) & 7;
        sp.pulse4 = (combo >> 11) & 7;
        sp.sign = (combo >> 14) & 1;
        sp.shift = (combo >> 15) & 1;
        tetra_codevector(a, 19, &sp, code); /* the shortest integer lag (of 19 1/3): sharpened */
        for (n = 0; n < TETRA_SUBFRAME; n++) {
            e += (int64_t)code[n] * code[n];
        }
        if (e < least) {
            least = e;
        }
    }
    assert(least >= (int64_t)1696 * 1696);
    assert((least >> 16) >= 43);
}

static void
test_energy_limits(void) {
    tetra_acelp_state st;
    tetra_frame_params fp;
    int16_t pcm[TETRA_FRAME];
    int frame;
    int k;
    int gain;
    int best_p = 0;
    int best_c = 0;
    /* The gain indices with the largest pitch and code prediction errors. */
    for (gain = 1; gain < 64; gain++) {
        if (tetra_gain_cb[gain][0] > tetra_gain_cb[best_p][0]) {
            best_p = gain;
        }
        if (tetra_gain_cb[gain][1] > tetra_gain_cb[best_c][1]) {
            best_c = gain;
        }
    }
    tetra_acelp_reset(&st);
    memset(&fp, 0, sizeof(fp));
    fp.lsp_index[0] = 60;
    fp.lsp_index[1] = 90;
    fp.lsp_index[2] = 90;
    for (frame = 0; frame < 20; frame++) {
        for (k = 0; k < TETRA_SUBFRAMES; k++) {
            fp.sub[k].gain = (frame + k) % 2 ? best_p : best_c;
        }
        tetra_acelp_decode(&st, &fp, pcm);
        assert(st.energy_p <= 6912 && st.energy_c <= 6400); /* 27.0 and 25.0 */
    }
    assert(st.energy_p == 6912 || st.energy_c == 6400);
}

static void
test_excitation_wraps(void) {
    /* u = (gp v + gc c + 2048) >> 12, keeping the low 16 bits as the reference decoder does. */
    assert(tetra_excitation(4915, 32767, 32767, 32767) == -26233);  /* 301447 */
    assert(tetra_excitation(4915, -32768, -32768, 32767) == 26224); /* -301456 */
    assert(tetra_excitation(1, -1, 4369, 30720) == 32767);          /* 134215679 rounds to 32767 */
    assert(tetra_excitation(0, 0, 4369, 30720) == -32768);          /* 134215680 rounds to 32768 */
    assert(tetra_excitation(1, -2048, 4096, -32768) == -32768);     /* rounds to -32768 */
    assert(tetra_excitation(1, -2049, 4096, -32768) == 32767);      /* rounds to -32769 */
    assert(tetra_excitation(1, -2048, 0, 0) == 0 && tetra_excitation(1, 2048, 0, 0) == 1); /* ties up */
}

/* === Concealment === */

static void
test_bad_frames_after_reset(void) {
    tetra_acelp_state st;
    tetra_acelp_state expect;
    int16_t pcm[TETRA_FRAME];
    int frame;
    lcg g = {9u};
    int pass;
    tetra_acelp_reset(&st);
    tetra_acelp_reset(&expect);
    /* Silent, and like the reference decoder the previous frame's LSPs lose their first LSP,
     * however many bad frames come before the first good one; the same after a reset of a used
     * state. */
    expect.lsp_old[0] = 0;
    for (pass = 0; pass < 2; pass++) {
        for (frame = 0; frame < 256; frame++) {
            memset(pcm, 0x55, sizeof(pcm));
            tetra_acelp_decode(&st, NULL, pcm);
            assert(pcm_is_silent(pcm));
            assert(memcmp(&st, &expect, sizeof(st)) == 0);
        }
        for (frame = 0; frame < 5; frame++) {
            uint8_t bits[MBE_TETRA_FRAME_BITS];
            tetra_frame_params fp;
            voiced_frame(&g, bits);
            tetra_acelp_unpack(bits, &fp);
            tetra_acelp_decode(&st, &fp, pcm);
        }
        tetra_acelp_reset(&st);
    }
}

/** Decode one good frame whose LSPs are rejected (they stay crossed after widening). */
static void
decode_rejected_lsp_frame(tetra_acelp_state* st, lcg* g) {
    uint8_t bits[MBE_TETRA_FRAME_BITS];
    tetra_frame_params fp;
    int16_t pcm[TETRA_FRAME];
    voiced_frame(g, bits);
    tetra_acelp_unpack(bits, &fp);
    fp.lsp_index[0] = 0;
    fp.lsp_index[1] = 0;
    fp.lsp_index[2] = 0;
    tetra_acelp_decode(st, &fp, pcm);
}

static void
test_bad_start_transitions(void) {
    tetra_acelp_state st;
    int16_t pcm[TETRA_FRAME];
    int16_t decoded[TETRA_M];
    int32_t ep;
    int32_t ec;
    lcg g = {13u};
    static const int zero_index[3] = {0, 0, 0};
    tetra_lsp_decode(zero_index, decoded);
    assert(!tetra_lsp_ordered(decoded));

    /* A good frame ends the bad start even when its LSPs are rejected: the next bad frame is
     * concealed (energies lowered) and the reset LSPs stay as they are. */
    tetra_acelp_reset(&st);
    decode_rejected_lsp_frame(&st, &g);
    assert(st.good_seen);
    assert(memcmp(st.lsp_old, k_lsp_init, sizeof(k_lsp_init)) == 0);
    ep = st.energy_p;
    ec = st.energy_c;
    assert(ep > 0 || ec > 0);
    tetra_acelp_decode(&st, NULL, pcm);
    assert(st.energy_p == (ep - 512 > 0 ? ep - 512 : 0));
    assert(st.energy_c == (ec - 512 > 0 ? ec - 512 : 0));
    assert(memcmp(st.lsp_old, k_lsp_init, sizeof(k_lsp_init)) == 0);

    /* After a bad start, a good frame with rejected LSPs keeps the vector the bad start left, and
     * a later bad frame conceals from it. */
    tetra_acelp_reset(&st);
    tetra_acelp_decode(&st, NULL, pcm);
    decode_rejected_lsp_frame(&st, &g);
    assert(st.good_seen);
    assert(st.lsp_old[0] == 0 && memcmp(&st.lsp_old[1], &k_lsp_init[1], (TETRA_M - 1) * sizeof(int16_t)) == 0);
    ep = st.energy_p;
    ec = st.energy_c;
    assert(ep > 0 || ec > 0);
    tetra_acelp_decode(&st, NULL, pcm);
    assert(st.energy_p == (ep - 512 > 0 ? ep - 512 : 0));
    assert(st.energy_c == (ec - 512 > 0 ? ec - 512 : 0));
    assert(st.lsp_old[0] == 0);
}

static void
test_concealment_state(void) {
    tetra_acelp_state st;
    tetra_frame_params fp;
    int16_t pcm[TETRA_FRAME];
    int16_t lsp_before[TETRA_M];
    lcg g = {11u};
    uint8_t bits[MBE_TETRA_FRAME_BITS];
    int32_t ep;
    int32_t ec;
    int t0;
    int frac;
    int frame;

    tetra_acelp_reset(&st);
    voiced_frame(&g, bits);
    tetra_acelp_unpack(bits, &fp);
    tetra_acelp_decode(&st, &fp, pcm);
    memcpy(lsp_before, st.lsp_old, sizeof(lsp_before));
    tetra_decode_pitch(0, fp.sub[0].pitch, 0, &t0, &frac);
    tetra_decode_pitch(3, fp.sub[3].pitch, t0, &t0, &frac);
    assert(st.last_pitch == t0);

    ep = st.energy_p;
    ec = st.energy_c;
    tetra_acelp_decode(&st, NULL, pcm);
    /* Each subframe lowers both energies by 0.5 (128 in Q8), never below zero. */
    assert(st.energy_p == (ep - 512 > 0 ? ep - 512 : 0));
    assert(st.energy_c == (ec - 512 > 0 ? ec - 512 : 0));
    assert(memcmp(st.lsp_old, lsp_before, sizeof(lsp_before)) == 0);
    assert(memcmp(st.last_code, fp.sub, sizeof(fp.sub)) == 0);
    assert(st.last_pitch == t0);

    for (frame = 0; frame < 40; frame++) {
        tetra_acelp_decode(&st, NULL, pcm);
    }
    assert(st.energy_p == 0 && st.energy_c == 0);
}

static void
test_crossed_lsps_keep_previous(void) {
    tetra_acelp_state st;
    tetra_frame_params fp;
    int16_t pcm[TETRA_FRAME];
    int16_t decoded[TETRA_M];
    lcg g = {5u};
    /* q3 lies far below q4 here, and widening that boundary pushes q4 past q5. */
    fp.lsp_index[0] = 254;
    fp.lsp_index[1] = 73;
    fp.lsp_index[2] = 415;
    tetra_lsp_decode(fp.lsp_index, decoded);
    assert(!tetra_lsp_ordered(decoded));
    {
        uint8_t bits[MBE_TETRA_FRAME_BITS];
        tetra_frame_params good;
        voiced_frame(&g, bits);
        tetra_acelp_unpack(bits, &good);
        memcpy(fp.sub, good.sub, sizeof(fp.sub));
        tetra_acelp_reset(&st);
        tetra_acelp_decode(&st, &good, pcm);
    }
    {
        int16_t before[TETRA_M];
        memcpy(before, st.lsp_old, sizeof(before));
        tetra_acelp_decode(&st, &fp, pcm);
        assert(memcmp(st.lsp_old, before, sizeof(before)) == 0);
    }
}

/* === Public API === */

static void
test_api_arguments(void) {
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    uint8_t bits[MBE_TETRA_FRAME_BITS];
    int16_t pcm[MBE_TETRA_FRAME_SAMPLES];
    assert(dec);
    memset(bits, 0, sizeof(bits));
    assert(mbe_tetraDecodeFrame(NULL, bits, 0, pcm) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_tetraDecodeFrame(dec, bits, 0, NULL) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_tetraDecodeFrame(dec, NULL, 0, pcm) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_tetraDecodeFrame(dec, NULL, 1, pcm) == 0);
    mbe_tetraDecoderReset(NULL);
    mbe_tetraDecoderFree(NULL);
    mbe_tetraDecoderFree(dec);
}

/** Decode `frames` voiced frames from seed into one hash, frame `reject_at` first sent with bad bits. */
static uint32_t
run_with_rejection(int frames, int reject_at, uint8_t bad_value) {
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    lcg g = {21u};
    uint32_t h = 2166136261u;
    int f;
    assert(dec);
    for (f = 0; f < frames; f++) {
        uint8_t bits[MBE_TETRA_FRAME_BITS];
        int16_t pcm[MBE_TETRA_FRAME_SAMPLES];
        voiced_frame(&g, bits);
        if (f == reject_at) {
            int16_t sentinel[MBE_TETRA_FRAME_SAMPLES];
            uint8_t saved = bits[70];
            memset(pcm, 0x5A, sizeof(pcm));
            memcpy(sentinel, pcm, sizeof(pcm));
            bits[70] = bad_value;
            assert(mbe_tetraDecodeFrame(dec, bits, 0, pcm) == MBE_STATUS_INVALID_BITS);
            assert(memcmp(pcm, sentinel, sizeof(pcm)) == 0);
            bits[70] = saved;
        }
        assert(mbe_tetraDecodeFrame(dec, bits, 0, pcm) == 0);
        h = fnv1a_pcm(h, pcm, MBE_TETRA_FRAME_SAMPLES);
    }
    mbe_tetraDecoderFree(dec);
    return h;
}

static void
test_api_rejection_leaves_state(void) {
    uint32_t clean = run_with_rejection(30, -1, 0);
    assert(run_with_rejection(30, 12, 2) == clean);
    assert(run_with_rejection(30, 12, 0xFF) == clean);
    assert(run_with_rejection(30, 0, 0x80) == clean);
}

/** Decode with bad frames first..first + 2, their payloads and indicator values chosen by `variant`. */
static uint32_t
run_bad_frames(int variant, int first) {
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    lcg g = {33u};
    uint32_t h = 2166136261u;
    int f;
    assert(dec);
    for (f = 0; f < 20; f++) {
        uint8_t bits[MBE_TETRA_FRAME_BITS];
        int16_t pcm[MBE_TETRA_FRAME_SAMPLES];
        voiced_frame(&g, bits);
        if (f >= first && f < first + 3) {
            static const int indicators[3] = {1, 2, -1};
            const uint8_t* payload = bits;
            if (variant == 1) {
                payload = NULL;
            } else if (variant == 2) {
                memset(bits, 0xFF, sizeof(bits)); /* not read: no MBE_STATUS_INVALID_BITS */
            }
            assert(mbe_tetraDecodeFrame(dec, payload, variant == 3 ? indicators[f - first] : 1, pcm) == 0);
        } else {
            assert(mbe_tetraDecodeFrame(dec, bits, 0, pcm) == 0);
        }
        h = fnv1a_pcm(h, pcm, MBE_TETRA_FRAME_SAMPLES);
    }
    mbe_tetraDecoderFree(dec);
    return h;
}

static void
test_api_bad_frame_payload_ignored(void) {
    int first;
    for (first = 0; first <= 5; first += 5) { /* at the start (before any good frame) and later */
        uint32_t base = run_bad_frames(0, first);
        assert(run_bad_frames(1, first) == base);
        assert(run_bad_frames(2, first) == base);
        assert(run_bad_frames(3, first) == base);
    }
}

static void
test_api_reset_and_reentrancy(void) {
    mbe_tetra_decoder* a = mbe_tetraDecoderAlloc();
    mbe_tetra_decoder* b = mbe_tetraDecoderAlloc();
    mbe_tetra_decoder* solo = mbe_tetraDecoderAlloc();
    lcg ga = {101u};
    lcg gb = {202u};
    lcg gs = {101u};
    int f;
    assert(a && b && solo);
    /* Interleaving two streams gives each the output it has alone. */
    for (f = 0; f < 40; f++) {
        uint8_t bits_a[MBE_TETRA_FRAME_BITS];
        uint8_t bits_b[MBE_TETRA_FRAME_BITS];
        uint8_t bits_s[MBE_TETRA_FRAME_BITS];
        int16_t pa[MBE_TETRA_FRAME_SAMPLES];
        int16_t pb[MBE_TETRA_FRAME_SAMPLES];
        int16_t ps[MBE_TETRA_FRAME_SAMPLES];
        voiced_frame(&ga, bits_a);
        random_bits(&gb, bits_b);
        voiced_frame(&gs, bits_s);
        assert(mbe_tetraDecodeFrame(a, bits_a, f % 9 == 4, pa) == 0);
        assert(mbe_tetraDecodeFrame(b, bits_b, f % 5 == 1, pb) == 0);
        assert(mbe_tetraDecodeFrame(solo, bits_s, f % 9 == 4, ps) == 0);
        assert(memcmp(pa, ps, sizeof(pa)) == 0);
    }
    /* After a reset, a used context decodes like a new one. */
    mbe_tetraDecoderReset(b);
    gb.state = 101u;
    ga.state = 101u;
    mbe_tetraDecoderReset(a);
    for (f = 0; f < 10; f++) {
        uint8_t bits[MBE_TETRA_FRAME_BITS];
        int16_t pa[MBE_TETRA_FRAME_SAMPLES];
        int16_t pb[MBE_TETRA_FRAME_SAMPLES];
        voiced_frame(&ga, bits);
        assert(mbe_tetraDecodeFrame(a, bits, 0, pa) == 0);
        voiced_frame(&gb, bits);
        assert(mbe_tetraDecodeFrame(b, bits, 0, pb) == 0);
        assert(memcmp(pa, pb, sizeof(pa)) == 0);
    }
    mbe_tetraDecoderFree(a);
    mbe_tetraDecoderFree(b);
    mbe_tetraDecoderFree(solo);
}

/* === Regression hashes ===
 *
 * FNV-1a over the little-endian PCM of whole streams. The decoder is integer arithmetic only,
 * so these hold on every compiler, architecture and build type. They pin this decoder's output.
 * Checked outside the test suite: ETSI's reference decoder produces the same PCM for every stream
 * here. The random stream with losses wraps codevector samples, the overloaded stream also wraps
 * code gains and excitation samples, and the bad-start stream starts from the reference's state
 * after bad frames (docs/tetra-acelp.md).
 */

typedef enum stream_kind { STREAM_VOICED, STREAM_RANDOM, STREAM_OVERLOADED } stream_kind;

/** Random LSPs, pulses near the subframe end and high gain indices: large code gains. */
static void
overloaded_frame(lcg* g, uint8_t bits[MBE_TETRA_FRAME_BITS]) {
    tetra_frame_params fp;
    int k;
    fp.lsp_index[0] = (int)(lcg_next(g) % 256u);
    fp.lsp_index[1] = (int)(lcg_next(g) % 512u);
    fp.lsp_index[2] = (int)(lcg_next(g) % 512u);
    for (k = 0; k < TETRA_SUBFRAMES; k++) {
        tetra_subframe_params* sp = &fp.sub[k];
        sp->pitch = (int)(lcg_next(g) % (k == 0 ? 256u : 32u));
        sp->pulse4 = 6 + (int)(lcg_next(g) % 2u);
        sp->pulse3 = 6 + (int)(lcg_next(g) % 2u);
        sp->pulse2 = 6 + (int)(lcg_next(g) % 2u);
        sp->pulse1 = 26 + (int)(lcg_next(g) % 6u);
        sp->sign = (int)(lcg_next(g) % 2u);
        sp->shift = (int)(lcg_next(g) % 2u);
        sp->gain = 40 + (int)(lcg_next(g) % 24u);
    }
    pack_frame(&fp, bits);
}

/** Hash `frames` frames; frame f is bad when bad(f) says so. */
static uint32_t
stream_hash(stream_kind kind, uint32_t seed, int frames, int (*bad)(int, uint32_t)) {
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    lcg g = {seed};
    lcg loss = {seed ^ 0xA5A5A5A5u};
    uint32_t h = 2166136261u;
    int f;
    assert(dec);
    for (f = 0; f < frames; f++) {
        uint8_t bits[MBE_TETRA_FRAME_BITS];
        int16_t pcm[MBE_TETRA_FRAME_SAMPLES];
        if (kind == STREAM_VOICED) {
            voiced_frame(&g, bits);
        } else if (kind == STREAM_RANDOM) {
            random_bits(&g, bits);
        } else {
            overloaded_frame(&g, bits);
        }
        assert(mbe_tetraDecodeFrame(dec, bits, bad(f, lcg_next(&loss)), pcm) == 0);
        h = fnv1a_pcm(h, pcm, MBE_TETRA_FRAME_SAMPLES);
    }
    mbe_tetraDecoderFree(dec);
    return h;
}

static int
no_losses(int f, uint32_t r) {
    (void)f;
    (void)r;
    return 0;
}

static int
five_percent(int f, uint32_t r) {
    (void)f;
    return r % 100u < 5u;
}

static int
bad_start(int f, uint32_t r) {
    (void)r;
    return f < 5;
}

static int
long_erasure(int f, uint32_t r) {
    (void)r;
    return f >= 50 && f < 350;
}

static int
alternate(int f, uint32_t r) {
    (void)r;
    return f % 2 == 1;
}

static void
test_regression_hashes(void) {
    assert(stream_hash(STREAM_VOICED, 1u, 2000, no_losses) == 0x57F200F9u);
    assert(stream_hash(STREAM_VOICED, 2u, 2000, five_percent) == 0xC1F40AD8u);
    assert(stream_hash(STREAM_RANDOM, 3u, 2000, no_losses) == 0x1195B24Fu);
    assert(stream_hash(STREAM_RANDOM, 4u, 2000, five_percent) == 0xEE127CFCu);
    assert(stream_hash(STREAM_VOICED, 5u, 200, bad_start) == 0xCC61F145u);
    assert(stream_hash(STREAM_VOICED, 6u, 450, long_erasure) == 0x40CC6E44u);
    assert(stream_hash(STREAM_VOICED, 7u, 200, alternate) == 0x30360BE9u);
    assert(stream_hash(STREAM_OVERLOADED, 12u, 400, five_percent) == 0xB296BA68u);
}

int
main(void) {
    test_unpack_round_trip();
    test_unpack_field_positions();
    test_pitch_first_subframe();
    test_pitch_relative();
    test_pulse_positions();
    test_codevector_absent_pulses();
    test_codevector_sharpening();
    test_codevector_wraps();
    test_lsp_decode_and_widening();
    test_lsp_ordering();
    test_lsp_interpolation();
    test_lsp_to_lpc();
    test_lsp_to_lpc_reference();
    test_synthesis_filter();
    test_adaptive_vector();
    test_log2_pow2();
    test_gains();
    test_code_energy_floor();
    test_energy_limits();
    test_excitation_wraps();
    test_bad_frames_after_reset();
    test_bad_start_transitions();
    test_concealment_state();
    test_crossed_lsps_keep_previous();
    test_api_arguments();
    test_api_rejection_leaves_state();
    test_api_bad_frame_payload_ignored();
    test_api_reset_and_reentrancy();
    test_regression_hashes();
    return 0;
}
