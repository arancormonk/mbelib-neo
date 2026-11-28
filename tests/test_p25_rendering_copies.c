// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief P25's copies of the shared mid-frame and own-model code against the
 * shared code, on the behaviour they share.
 *
 * P25's rendering has its own copies of the mid-frame filler and of the column
 * helpers (p25_fill_mid_bands(), p25_column_band_state(), p25_columns_voice()
 * and p25_own_model() in src/core/mbe_subframe.c), so that the other
 * renderings' code keeps its inlining: one body serving both moved D-STAR's
 * and rate 33's samples in the 32-bit Release build. The copies differ from
 * the shared code only where P25's rules do (the phase step on unvoiced bands,
 * each frame's state as its own model plays it by the columns). This test
 * compiles that file's static functions here (the library's own object is
 * not linked: every symbol it defines is defined by this translation unit)
 * and checks, over every ordered pair of the 208 P25 pitch codes, each pair
 * with one of 64 pairs of voicing classes (every class pair on many code
 * pairs), that the copies give exactly what the shared code gives where
 * P25's rules do not apply: the column tests on P25's samples, each band's
 * mid-frame state with neither frame's own model read by the columns, the
 * mid-frame bands with no unvoiced band carrying a finite log2Ml (no step)
 * and no column rule (or one with neither frame's own model read by it), and
 * the own model with no step, the columns applied to the shared one.
 */
#include "../src/core/mbe_subframe.c"

#include <assert.h>
#include <stdio.h>

/* Whether two floats have the same bits. */
static int
same_float(const float* a, const float* b) {
    uint32_t x, y;
    memcpy(&x, a, sizeof(x));
    memcpy(&y, b, sizeof(y));
    return x == y;
}

/* Whether two subframe models are the same, floats by their bits. */
static int
same_subframe(const struct mbe_subframe_model* a, const struct mbe_subframe_model* b) {
    int same =
        same_float(&a->w0, &b->w0) && a->L == b->L && a->any_voiced == b->any_voiced && a->any_pulsed == b->any_pulsed;
    for (int l = 0; l < 57; ++l) {
        same &= a->voiced[l] == b->voiced[l] && same_float(&a->magnitude[l], &b->magnitude[l])
                && same_float(&a->phase[l], &b->phase[l]);
    }
    return same;
}

/* The model P25's decoder gives pitch code b0 (w0, L and K), its bands in
 * voicing class `kind`, magnitudes and log2 amplitudes that vary by band. */
static void
code_model(mbe_parms* mp, int b0, int kind) {
    static const unsigned b0_index[8] = {0, 1, 2, 3, 4, 5, 85, 86};
    char imbe_d[88] = {0};
    mbe_parms prev, enhanced;
    for (int i = 0; i < 8; ++i) {
        imbe_d[b0_index[i]] = (char)((b0 >> (7 - i)) & 1);
    }
    mbe_initMbeParms(mp, &prev, &enhanced);
    const int rc = mbe_decodeImbe4400Parms(imbe_d, mp, &prev);
    assert(rc == 0);
    (void)rc;
    const int L = mp->L;
    for (int l = 1; l <= 56; ++l) {
        int v;
        switch (kind) {
            case 0: v = 0; break;                                 /* none */
            case 1: v = 1; break;                                 /* all */
            case 2: v = l <= L / 2; break;                        /* a low region */
            case 3: v = ((l - 1) / 3) % 2 == 0; break;            /* alternate bands of three */
            case 4: v = l == 7; break;                            /* a lone harmonic */
            case 5: v = 3 * l > 2 * L; break;                     /* the top only */
            case 6: v = (l == 5) ? 2 : (l <= (2 * L) / 3); break; /* a band with no state */
            default: v = l == L; break;                           /* the last band alone */
        }
        mp->Vl[l] = v;
        mp->Ml[l] = (float)(1.0 + 0.5 * sin(0.37 * l + b0));
        mp->log2Ml[l] = (float)(2.0 + 1.5 * cos(0.21 * l + 0.1 * b0));
    }
}

/* A copy of a model with no unvoiced band carrying a finite log2Ml: the
 * unvoiced bands' log2Ml not a number, so that P25's phase source adds no
 * step and enters the floor, as the shared one does. */
static void
without_step(mbe_parms* mp) {
    const uint32_t nan_bits = 0x7FC00000u;
    for (int l = 1; l <= mp->L; ++l) {
        if (mp->Vl[l] == 0) {
            memcpy(&mp->log2Ml[l], &nan_bits, sizeof(nan_bits));
        }
    }
}

/* The mid-frame model's pitch and count as mbe_subframe_mid_model() sets
 * them, before its bands are filled. */
static int
mid_frame(const mbe_parms* prev, const mbe_parms* cur, struct mbe_subframe_model* model) {
    memset(model, 0, sizeof(*model));
    const int shorter = (prev->L < cur->L) ? prev->L : cur->L;
    const int longer = (prev->L < cur->L) ? cur->L : prev->L;
    const float w0 = mid_pitch(prev, cur, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    if (!positive_normal(&w0)) {
        return 0;
    }
    const float harmonics = floorf((0.9254f * (float)M_PI) / w0);
    model->w0 = w0;
    model->L = (harmonics < (float)shorter) ? shorter : ((harmonics > (float)longer) ? longer : (int)harmonics);
    return 1;
}

/* The column tests and the band states of one pair. */
static long
check_columns(const mbe_parms* prev, const mbe_parms* cur, const struct mbe_subframe_model* model) {
    struct column_rule rule;
    memset(&rule, 0, sizeof(rule));
    if (!column_rule(prev, cur, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &rule)) {
        return 0;
    }
    assert(!rule.dstar);
    for (int l = 1; l <= model->L; ++l) {
        assert(p25_columns_voice(&rule, model->w0, l) == columns_voice(&rule, model->w0, l));
        if (!p25_reads_by_columns(prev) && !p25_reads_by_columns(cur)) {
            assert(p25_column_band_state(prev, cur, &rule, model, l) == column_band_state(prev, cur, &rule, model, l));
        }
    }
    return model->L;
}

/* The mid-frame bands of one pair, with no step: no column rule, and the rule
 * where neither frame's own model is read by it. */
static int
check_fill(const mbe_parms* prev, const mbe_parms* cur, const struct mbe_subframe_model* model, int late) {
    mbe_parms a = *prev, b = *cur;
    without_step(&a);
    without_step(&b);
    struct mbe_subframe_model shared = *model, copy = *model;
    fill_mid_bands(&a, &b, &b, late, NULL, NULL, &shared);
    p25_fill_mid_bands(&a, &b, &b, late, NULL, NULL, &copy);
    assert(same_subframe(&shared, &copy));
    struct column_rule rule;
    memset(&rule, 0, sizeof(rule));
    if (column_rule(&a, &b, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &rule) && !p25_reads_by_columns(&a)
        && !p25_reads_by_columns(&b)) {
        shared = *model;
        copy = *model;
        fill_mid_bands(&a, &b, &b, late, &rule, NULL, &shared);
        p25_fill_mid_bands(&a, &b, &b, late, &rule, NULL, &copy);
        assert(same_subframe(&shared, &copy));
        return 1;
    }
    return 0;
}

/* The own model, with no step: P25's against the shared one with the columns
 * applied to it. */
static void
check_own(const mbe_parms* mp) {
    mbe_parms a = *mp;
    without_step(&a);
    struct mbe_subframe_model shared, copy;
    mbe_subframe_own_model(&a, &shared);
    p25_column_voicing(&a, &shared);
    p25_own_model(&a, &copy);
    assert(same_subframe(&shared, &copy));
}

int
main(void) {
    static mbe_parms models[8][208];
    for (int kind = 0; kind < 8; ++kind) {
        for (int b0 = 0; b0 < 208; ++b0) {
            code_model(&models[kind][b0], b0, kind);
            check_own(&models[kind][b0]);
        }
    }
    long bands = 0, pairs = 0, ruled = 0;
    int classes[64] = {0};
    for (int a = 0; a < 208; ++a) {
        for (int b = 0; b < 208; ++b) {
            const int c = ((a * 7) + b) % 64;
            const mbe_parms* prev = &models[c / 8][a];
            mbe_parms cur = models[c % 8][b];
            if ((c % 8) == 0 && ((a + b) % 2) == 0) {
                /* A frame with no voiced band lacking a pitch: the column rule
                 * applies, and with no voiced band before it neither own
                 * model is read by it. */
                mbe_set_frame_mark(&cur, MBE_FRAME_MARK_PITCH_ABSENT);
            }
            struct mbe_subframe_model model;
            if (!mid_frame(prev, &cur, &model)) {
                continue;
            }
            bands += check_columns(prev, &cur, &model);
            ruled += check_fill(prev, &cur, &model, (a + b) % 3 == 0);
            ++pairs;
            classes[c] = 1;
        }
    }
    for (int c = 0; c < 64; ++c) {
        assert(classes[c]);
    }
    assert(pairs == 208L * 208L && bands > 0 && ruled > 0);
    printf("p25 rendering copies: %ld code pairs, %ld mid-frame bands by the columns, %ld fills under the rule, "
           "alike\n",
           pairs, bands, ruled);
    return 0;
}
