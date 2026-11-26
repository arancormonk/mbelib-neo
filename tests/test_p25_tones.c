// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief P25 tones sent as voice frames: the decoder's tone test and the tone
 * model the synthesis plays for a frame that passes it.
 *
 * DVSI's P25 encoder sends a tone as voice frames at the pitch code of the
 * tone's Annex J fundamental, with one voicing word per tone and the tone's
 * harmonics large; DVSI's decoder recognises such a frame and plays a tone
 * model. The frames here are built with the IMBE encoder's quantizer
 * (mbe_imbe4400_quantize_amplitudes()) from a pitch code, a voicing word and
 * log2 amplitudes, and decoded as the library decodes P25.
 */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "imbe4400_internal.h"
#include "mbe_adaptive.h"
#include "mbe_subframe.h"
#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* The word the decoder leaves in Vl[0] for a tone frame (imbe7200x4400.c). */
#define TONE_FLAG (1 << 24)

static int
tone_word(int b0, int l1, int l2) {
    return TONE_FLAG | (b0 << 16) | (l2 << 8) | l1;
}

/* L and K of pitch code b0 as the decoder derives them (TIA-102.BABA eqs 46-48). */
static int
harmonics_of(int b0, int* K) {
    const float w0 = (float)(4 * M_PI) / (float)((float)b0 + 39.5);
    const int L = (int)(0.9254 * (int)((M_PI / w0) + 0.25));
    *K = (L < 37) ? (L + 2) / 3 : 12;
    return L;
}

/* A P25 voice frame at pitch code b0 with the K-bit voicing word `word` (band
 * 1 the most significant bit) and log2 amplitude logm[l], quantized against
 * the history `prev`. */
static void
build_frame(char imbe_d[88], int b0, unsigned word, const float logm[57], const mbe_parms* prev) {
    int K;
    (void)harmonics_of(b0, &K);
    unsigned char bands[12] = {0};
    for (int band = 1; band <= K; ++band) {
        bands[band - 1] = (unsigned char)((word >> (K - band)) & 1u);
    }
    mbe_imbe4400_quantize_amplitudes(b0, bands, logm, prev, 0, imbe_d);
}

/* log2 amplitudes: `tone` at the harmonics l1 and l2, `rest` elsewhere. */
static void
tone_amplitudes(float logm[57], int L, int l1, int l2, float tone, float rest) {
    for (int l = 0; l < 57; ++l) {
        logm[l] = (l >= 1 && l <= L && (l == l1 || l == l2)) ? tone : rest;
    }
}

/* The rms of the decoded magnitudes off the tone over that of the tone's. */
static double
others_over_tone(const mbe_parms* mp, int l1, int l2) {
    double tone = 0.0, others = 0.0;
    int tones = 0;
    for (int l = 1; l <= mp->L; ++l) {
        const double power = (double)mp->Ml[l] * (double)mp->Ml[l];
        if (l == l1 || l == l2) {
            tone += power;
            ++tones;
        } else {
            others += power;
        }
    }
    return sqrt((others / (double)(mp->L - tones)) / (tone / (double)tones));
}

/* Decode one frame as P25 from a stream's history (left unchanged). */
static mbe_parms
decode(const char imbe_d[88], const mbe_parms* history) {
    mbe_parms cur, prev = *history;
    memset(&cur, 0, sizeof(cur));
    const int rc = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
    assert(rc == 0);
    (void)rc;
    return cur;
}

static mbe_parms
fresh_history(void) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    return prev;
}

/* A frame's bits and decoded model after eight frames of the same target: the
 * prediction then carries most of the spread, so the decoded magnitudes come
 * near the target's. */
static mbe_parms
settled(char imbe_d[88], int b0, unsigned word, const float logm[57]) {
    mbe_parms history = fresh_history(), frame;
    for (int f = 0; f < 8; ++f) {
        build_frame(imbe_d, b0, word, logm, &history);
        frame = decode(imbe_d, &history);
        history = frame;
    }
    return frame;
}

/* sinc^2(pi l / 256). */
static double
droop(int l) {
    const double x = M_PI * (double)l / 256.0;
    return (sin(x) / x) * (sin(x) / x);
}

/* The tone model's pitch for b0: eq 46's fundamental on a 2^-19 grid of cycles
 * per sample, rounded to float once. */
static float
tone_pitch(int b0) {
    return (float)(2.0 * M_PI * floor(2097152.0 / (double)((2 * b0) + 79)) / 524288.0);
}

static int
same_float(float a, float b) {
    uint32_t x, y;
    memcpy(&x, &a, sizeof(x));
    memcpy(&y, &b, sizeof(y));
    return x == y;
}

/* Two subframe models alike, every float by its bits. */
static int
same_subframe_model(const struct mbe_subframe_model* a, const struct mbe_subframe_model* b) {
    if (a->L != b->L || a->any_voiced != b->any_voiced || !same_float(a->w0, b->w0)) {
        return 0;
    }
    for (int l = 1; l <= a->L; ++l) {
        if (a->voiced[l] != b->voiced[l] || !same_float(a->magnitude[l], b->magnitude[l])
            || !same_float(a->phase[l], b->phase[l])) {
            return 0;
        }
    }
    return 1;
}

/*
 * The tone test on frames built from the table's pairs: DTMF 1 (b0 53, word
 * 1110000, harmonics 4 and 7) passes; the same amplitudes with another voicing
 * word do not; with a third harmonic raised, the frame passes exactly while
 * the others' rms stays under 0.375 of the tone's (passing near 0.3, failing
 * under 0.45); a call-progress pair (162: b0 92, 0100000000, harmonics 4 and
 * 5) does not; a single tone (30: b0 12, 1100, harmonic 3) passes, and where
 * several single tones share a pair (b0 17, 1000: harmonics 1, 2 and 3) the
 * dominant harmonic's row is taken; with the other harmonics raised, a single
 * tone passes exactly while their rms stays under 0.0235 of the tone's
 * (passing just under it and failing just over it).
 */
static void
test_tone_test(void) {
    char imbe_d[88];
    float logm[57];
    int K;

    const int L = harmonics_of(53, &K);
    assert(L == 21 && K == 7);
    tone_amplitudes(logm, L, 4, 7, 11.0f, 2.0f);
    mbe_parms frame = settled(imbe_d, 53, 0x70u, logm);
    assert(others_over_tone(&frame, 4, 7) < 0.05);
    assert(frame.Vl[0] == tone_word(53, 4, 7));

    frame = settled(imbe_d, 53, 0x78u, logm);
    assert(others_over_tone(&frame, 4, 7) < 0.05);
    assert(frame.Vl[0] == 0);

    /* The tone harmonics lower, the third raised past where the share crosses. */
    double highest_pass = 0.0, lowest_fail = 1.0e9;
    tone_amplitudes(logm, L, 4, 7, 8.0f, 2.0f);
    for (int step = 0; step <= 40; ++step) {
        logm[10] = 6.0f + (0.125f * (float)step);
        frame = settled(imbe_d, 53, 0x70u, logm);
        const double share = others_over_tone(&frame, 4, 7);
        const int passes = share < 0.375;
        assert((frame.Vl[0] == tone_word(53, 4, 7)) == passes);
        highest_pass = passes ? fmax(highest_pass, share) : highest_pass;
        lowest_fail = passes ? lowest_fail : fmin(lowest_fail, share);
    }
    assert(highest_pass > 0.3 && lowest_fail < 0.45);

    const int cp_L = harmonics_of(92, &K);
    assert(K == 10);
    tone_amplitudes(logm, cp_L, 4, 5, 11.0f, 2.0f);
    frame = settled(imbe_d, 92, 0x100u, logm);
    assert(others_over_tone(&frame, 4, 5) < 0.05 && frame.Vl[0] == 0);

    const int single_L = harmonics_of(12, &K);
    assert(K == 4);
    tone_amplitudes(logm, single_L, 3, 3, 11.0f, 2.0f);
    frame = settled(imbe_d, 12, 0xCu, logm);
    assert(frame.Vl[0] == tone_word(12, 3, 3));

    const int shared_L = harmonics_of(17, &K);
    for (int l = 1; l <= 3; ++l) {
        tone_amplitudes(logm, shared_L, l, l, 11.0f, 2.0f);
        frame = settled(imbe_d, 17, 0x8u, logm);
        assert(frame.Vl[0] == tone_word(17, l, l));
    }
    /* A single tone holds the other harmonics far lower: raised together. */
    highest_pass = 0.0;
    lowest_fail = 1.0e9;
    for (int step = 0; step <= 64; ++step) {
        tone_amplitudes(logm, single_L, 3, 3, 11.0f, 3.0f + (0.0625f * (float)step));
        frame = settled(imbe_d, 12, 0xCu, logm);
        const double share = others_over_tone(&frame, 3, 3);
        const int passes = share < 0.0235;
        assert((frame.Vl[0] == tone_word(12, 3, 3)) == passes);
        highest_pass = passes ? fmax(highest_pass, share) : highest_pass;
        lowest_fail = passes ? lowest_fail : fmin(lowest_fail, share);
    }
    assert(highest_pass > 0.021 && lowest_fail < 0.026);
    puts("p25 tones: the table's pairs, the voicing word, the 0.375 and 0.0235 shares and the call-progress pairs");
}

/* A stream's three states. */
struct stream {
    mbe_parms cur, prev, enhanced;
};

static void
stream_init(struct stream* s) {
    mbe_initMbeParms(&s->cur, &s->prev, &s->enhanced);
}

/* Play one frame of parameter bits through the P25 (or, with `provoice`, the
 * ProVoice) data path. */
static void
play(float out[160], struct stream* s, const char imbe_d[88], int provoice) {
    mbe_process_result result;
    mbe_initProcessResult(&result);
    if (provoice) {
        result.flags |= MBE_PROCESS_FLAG_PROVOICE;
    }
    const int rc = mbe_processImbe4400Dataf(out, &result, imbe_d, &s->cur, &s->prev, &s->enhanced);
    assert(rc >= 0);
    (void)rc;
}

/*
 * A tone frame's synthesis model is the tone model: the 2^-19 pitch of its
 * code, the tone's harmonics voiced at one amplitude (0.988 times the rms of
 * their decoded magnitudes, times sinc^2(pi l / 256)), every other harmonic
 * unvoiced at 0, L and the decoded log2Ml kept, marked as a tone; the
 * prediction history keeps the decoded frame. The mid-frame model of a marked
 * frame is its own, and advancing it unrendered leaves the pulse train and the
 * noise generator where rendering it does. A repeat of the tone frame plays the tone model again. The
 * same bits through the ProVoice path play the decoded model.
 */
static void
test_tone_model(void) {
    struct stream s;
    char imbe_d[88];
    float logm[57], out[160];
    int K;
    stream_init(&s);
    const int L = harmonics_of(53, &K);
    tone_amplitudes(logm, L, 4, 7, 11.0f, 2.0f);
    logm[7] = 10.2f; /* the two decoded tone magnitudes unequal */
    build_frame(imbe_d, 53, 0x70u, logm, &s.prev);
    const mbe_parms decoded = decode(imbe_d, &s.prev);
    const mbe_parms last = s.enhanced;
    play(out, &s, imbe_d, 0);

    const mbe_parms* model = &s.enhanced;
    assert(mbe_frame_mark(model) == MBE_FRAME_MARK_TONE);
    assert(same_float(model->w0, tone_pitch(53)) && model->L == decoded.L && model->K == decoded.K);
    const double rms = sqrt(0.5 * (((double)decoded.Ml[4] * decoded.Ml[4]) + ((double)decoded.Ml[7] * decoded.Ml[7])));
    assert(fabs((double)decoded.Ml[4] / (double)decoded.Ml[7] - 1.0) > 0.2);
    for (int l = 1; l <= model->L; ++l) {
        const int tone = (l == 4) || (l == 7);
        assert(model->Vl[l] == tone);
        if (tone) {
            assert(fabs((double)model->Ml[l] / (0.988 * rms * droop(l)) - 1.0) < 1e-6);
        } else {
            assert(model->Ml[l] == 0.0f);
        }
        assert(same_float(model->log2Ml[l], decoded.log2Ml[l]));
    }
    /* The prediction history is the decoded frame. */
    for (int l = 1; l <= decoded.L; ++l) {
        assert(same_float(s.prev.Ml[l], decoded.Ml[l]) && same_float(s.prev.log2Ml[l], decoded.log2Ml[l]));
        assert(s.prev.Vl[l] == decoded.Vl[l]);
    }
    assert(same_float(s.prev.w0, decoded.w0));

    /* The mid-frame model of a marked frame is its own. */
    struct mbe_subframe_model mid, own;
    mbe_subframe_mid_model(&last, model, 0, MBE_VOICED_SYNTHESIS_DVSI, &mid);
    mbe_subframe_own_model(model, &own);
    assert(same_subframe_model(&mid, &own));
    /* Advancing the marked frame unrendered leaves the train rendering it leaves. */
    mbe_parms rendered = *model, advanced = *model;
    float scratch[160] = {0.0f};
    mbe_synthesize_subframes(scratch, &rendered, &last, NULL, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI);
    mbe_advance_subframes(&advanced, &last, 0, MBE_VOICED_SYNTHESIS_DVSI);
    for (int i = 0; i < 96; ++i) {
        assert(same_float(rendered.noiseOverlap[i], advanced.noiseOverlap[i]));
    }
    for (int l = 0; l <= 56; ++l) {
        assert(same_float(rendered.PSIl[l], advanced.PSIl[l]) && same_float(rendered.PHIl[l], advanced.PHIl[l]));
    }
    assert(same_float(rendered.noiseSeed, advanced.noiseSeed));

    /* A frame the repeat criteria reject repeats it: the history's tone model
     * again, not the rejected frame's, itself a valid tone frame of another
     * tone (KNOX: b0 56, word 01100000, harmonics 4 and 7). */
    const int other_L = harmonics_of(56, &K);
    tone_amplitudes(logm, other_L, 4, 7, 11.0f, 2.0f);
    char rejected[88];
    build_frame(rejected, 56, 0x60u, logm, &s.prev);
    const mbe_parms rejected_model = decode(rejected, &s.prev);
    assert(rejected_model.Vl[0] == tone_word(56, 4, 7) && rejected_model.L != decoded.L);
    mbe_process_result result;
    mbe_initProcessResult(&result);
    result.total_errors = 6; /* past the parameter-bit API's repeat criterion */
    const int rc = mbe_processImbe4400Dataf(out, &result, rejected, &s.cur, &s.prev, &s.enhanced);
    assert(rc >= 0);
    (void)rc;
    assert((result.flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
    assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_TONE && same_float(s.enhanced.w0, tone_pitch(53)));
    assert(s.enhanced.L == decoded.L && s.enhanced.Vl[0] == tone_word(53, 4, 7));
    for (int l = 1; l <= s.enhanced.L; ++l) {
        assert(s.enhanced.Vl[l] == ((l == 4) || (l == 7)));
    }

    /* ProVoice: the decoded model, at the nominal pitch. */
    struct stream pv;
    stream_init(&pv);
    play(out, &pv, imbe_d, 1);
    assert(mbe_frame_mark(&pv.enhanced) != MBE_FRAME_MARK_TONE && pv.enhanced.Vl[0] == 0);
    assert(same_float(pv.enhanced.w0, (float)(4 * M_PI) / (float)(53.0f + 39.5)));
    assert(pv.enhanced.Ml[1] > 0.0f);
    puts("p25 tones: the tone model, the history kept, the own model at mid-frame, repeats, ProVoice untouched");
}

/* A fit of two components: x[i] ~ the sum over k of cosine[k] cos(2 pi
 * cycles[k] i) + sine[k] sin(2 pi cycles[k] i). */
struct components {
    double cosine[2];
    double sine[2];
};

/* Basis function j of the fit at sample i: component j / 2, its cosine for an
 * even j and its sine for an odd one. */
static double
basis_value(const double cycles[2], int j, int i) {
    const double angle = 2.0 * M_PI * cycles[j / 2] * (double)i;
    return (j % 2 == 0) ? cos(angle) : sin(angle);
}

static double
component_amplitude(const struct components* c, int k) {
    return hypot(c->cosine[k], c->sine[k]);
}

/* The fitted components at sample i (counted as the fit counted them). */
static double
components_at(const struct components* c, const double cycles[2], int i) {
    double y = 0.0;
    for (int k = 0; k < 2; ++k) {
        const double angle = 2.0 * M_PI * cycles[k] * (double)i;
        y += (c->cosine[k] * cos(angle)) + (c->sine[k] * sin(angle));
    }
    return y;
}

/* Least-squares fit of components at `cycles[k]` cycles per sample over
 * x[0..n), by Gaussian elimination on the 4x4 normal equations; returns the
 * residual's share of the energy. */
static double
fit_components(const float* x, int n, const double cycles[2], struct components* fit) {
    double a[4][4] = {{0.0}}, b[4] = {0.0}, c[4] = {0.0};
    for (int i = 0; i < n; ++i) {
        for (int r = 0; r < 4; ++r) {
            b[r] += basis_value(cycles, r, i) * (double)x[i];
            for (int q = 0; q < 4; ++q) {
                a[r][q] += basis_value(cycles, r, i) * basis_value(cycles, q, i);
            }
        }
    }
    for (int p = 0; p < 4; ++p) {
        for (int r = p + 1; r < 4; ++r) {
            const double m = a[r][p] / a[p][p];
            for (int q = p; q < 4; ++q) {
                a[r][q] -= m * a[p][q];
            }
            b[r] -= m * b[p];
        }
    }
    for (int p = 3; p >= 0; --p) {
        double sum = b[p];
        for (int q = p + 1; q < 4; ++q) {
            sum -= a[p][q] * c[q];
        }
        c[p] = sum / a[p][p];
    }
    fit->cosine[0] = c[0];
    fit->sine[0] = c[1];
    fit->cosine[1] = c[2];
    fit->sine[1] = c[3];
    double energy = 0.0, left = 0.0;
    for (int i = 0; i < n; ++i) {
        const double y = components_at(fit, cycles, i);
        energy += (double)x[i] * (double)x[i];
        left += ((double)x[i] - y) * ((double)x[i] - y);
    }
    return left / energy;
}

/*
 * Three tone frames of one DTMF tone render, in the last frame, two steady
 * components at harmonics 4 and 7 of the 2^-19 pitch, nothing else, at one
 * model amplitude (each component over its droop and the output high-pass's
 * gain at it, within 1%). The frame after them, a quiet unvoiced one, holds
 * the tone fading out over its first 80 samples: the tone of the frame before,
 * continued, times 1 - n / 80 (within 3% of its peak), and under 3% of it after.
 */
static void
test_tone_run_renders(void) {
    struct stream s;
    char imbe_d[88];
    float logm[57], out[4][160];
    int K;
    stream_init(&s);
    const int L = harmonics_of(53, &K);
    for (int f = 0; f < 3; ++f) {
        tone_amplitudes(logm, L, 4, 7, 11.0f, 2.0f);
        logm[7] = 10.5f;
        build_frame(imbe_d, 53, 0x70u, logm, &s.prev);
        play(out[f], &s, imbe_d, 0);
        assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_TONE);
    }
    const double fundamental = floor(2097152.0 / (2.0 * 53.0 + 79.0)) / 524288.0;
    const double cycles[2] = {4.0 * fundamental, 7.0 * fundamental};
    struct components steady;
    /* The second half, where the frame's own model plays alone. */
    const double residual = fit_components(&out[2][80], 80, cycles, &steady);
    assert(residual < 1e-4);
    double model[2];
    for (int k = 0; k < 2; ++k) {
        const int l = (k == 0) ? 4 : 7;
        model[k] =
            component_amplitude(&steady, k) / (droop(l) * (double)mbe_low_band_gain((float)(cycles[k] * 8000.0)));
    }
    assert(fabs(model[0] / model[1] - 1.0) < 0.01);

    /* The frame after: unvoiced and quiet (b0 150, which marks nothing). */
    int unvoiced_K;
    const int unvoiced_L = harmonics_of(150, &unvoiced_K);
    tone_amplitudes(logm, unvoiced_L, 0, 0, 0.0f, -6.0f);
    build_frame(imbe_d, 150, 0u, logm, &s.prev);
    play(out[3], &s, imbe_d, 0);
    assert(mbe_frame_mark(&s.enhanced) != MBE_FRAME_MARK_TONE);
    /* The tone of the last tone frame, as fitted, continued: it fades out over
     * the first 80 samples, under the quiet frame's noise. */
    double peak = 0.0, worst = 0.0;
    for (int n = 0; n < 160; ++n) {
        const double tone = components_at(&steady, cycles, n + 80);
        peak = fmax(peak, fabs(tone));
        if (n < 80) {
            worst = fmax(worst, fabs((double)out[3][n] - ((1.0 - (double)n / 80.0) * tone)));
        }
    }
    assert(worst < 0.03 * peak);
    /* After them nothing of the tone is left, beyond the narrow noise the
     * mid-frame model's unvoiced bands put near it (the log-domain mean of the
     * tone's magnitudes and the quiet frame's). */
    struct components left;
    (void)fit_components(&out[3][80], 80, cycles, &left);
    for (int k = 0; k < 2; ++k) {
        assert(component_amplitude(&left, k) < 0.03 * component_amplitude(&steady, k));
    }
    puts("p25 tones: a tone run renders its two components at one level; the frame after fades it out");
}

/* Two frames' samples alike, bit for bit. */
static int
same_samples(const float* a, const float* b, int n) {
    for (int i = 0; i < n; ++i) {
        if (!same_float(a[i], b[i])) {
            return 0;
        }
    }
    return 1;
}

/* Parameter bits no pitch code reads as a voice frame (b0 255): the data path
 * repeats the last frame. */
static void
special_frame(char imbe_d[88]) {
    memset(imbe_d, 1, 88);
}

/* A stream after one DTMF 1 frame (b0 53, word 1110000, harmonics 4 and 7):
 * its history holds the decoded frame and the tone word. */
static void
stream_after_tone(struct stream* s) {
    char imbe_d[88];
    float logm[57], out[160];
    int K;
    stream_init(s);
    const int L = harmonics_of(53, &K);
    tone_amplitudes(logm, L, 4, 7, 11.0f, 2.0f);
    build_frame(imbe_d, 53, 0x70u, logm, &s->prev);
    play(out, s, imbe_d, 0);
    assert(s->prev.Vl[0] == tone_word(53, 4, 7) && mbe_frame_mark(&s->enhanced) == MBE_FRAME_MARK_TONE);
}

/* The frame a repeat of `s` plays, through the P25 or ProVoice data path. */
static void
play_repeat(float out[160], struct stream* s, int provoice) {
    char imbe_d[88];
    special_frame(imbe_d);
    play(out, s, imbe_d, provoice);
}

/*
 * A word in Vl[0] counts as a tone only as a decode leaves it: a pitch code
 * with rows, harmonics that are one of that code's rows, the model at that
 * code's harmonic count, and tone magnitudes that are positive, finite and
 * normal. A repeat of a history holding any other word plays as the same
 * history without it: a code past the table (207, with a loud harmonic 1),
 * a code inside the table without rows (11), a row of another code (harmonics
 * 7 and 9 at b0 53), the word of another code's row in this code's model (b0
 * 56's 4 and 7), the real word in the model of another pitch code (b0 60), and
 * the real word with a tone magnitude negative, zero, subnormal, infinite or
 * not a number. The real word itself plays the tone.
 */
static void
test_impossible_tone_words(void) {
    static const uint32_t magnitudes[] = {0xC4FA0000u /* -2000 */, 0x00000000u, 0x00000001u, 0x7F800000u, 0x7FC00000u};
    struct stream base;
    float with_word[160], without[160];
    stream_after_tone(&base);

    struct stream s = base;
    play_repeat(with_word, &s, 0);
    assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_TONE && same_float(s.enhanced.w0, tone_pitch(53)));

    /* The model of another pitch code (b0 60, harmonics 4 and 7 loud, a voicing word no row has). */
    struct stream other;
    char imbe_d[88];
    float logm[57], out[160];
    int K;
    stream_init(&other);
    tone_amplitudes(logm, harmonics_of(60, &K), 4, 7, 11.0f, 2.0f);
    build_frame(imbe_d, 60, 0x7Fu, logm, &other.prev);
    play(out, &other, imbe_d, 0);
    assert(other.prev.Vl[0] == 0 && other.prev.L != base.prev.L);

    const int words = 4 + (int)(sizeof(magnitudes) / sizeof(magnitudes[0])) + 1;
    for (int c = 0; c < words; ++c) {
        struct stream forged = (c == words - 1) ? other : base;
        if (c == 0) {
            forged.prev.Vl[0] = tone_word(207, 1, 1);
            forged.prev.Ml[1] = 2000.0f;
        } else if (c == 1) {
            forged.prev.Vl[0] = tone_word(11, 3, 3);
        } else if (c == 2) {
            forged.prev.Vl[0] = tone_word(53, 7, 9);
        } else if (c == 3) {
            forged.prev.Vl[0] = tone_word(56, 4, 7);
        } else if (c < words - 1) {
            memcpy(&forged.prev.Ml[4], &magnitudes[c - 4], sizeof(float));
        } else {
            forged.prev.Vl[0] = tone_word(53, 4, 7);
        }
        struct stream plain = forged;
        plain.prev.Vl[0] = 0;
        play_repeat(with_word, &forged, 0);
        play_repeat(without, &plain, 0);
        assert(mbe_frame_mark(&forged.enhanced) != MBE_FRAME_MARK_TONE);
        assert(same_samples(with_word, without, 160));
        assert(same_float(forged.enhanced.w0, plain.enhanced.w0));
    }
    puts("p25 tones: a word counts only as a decode leaves it");
}

/*
 * Only the tone model is marked MBE_FRAME_MARK_TONE, by the synthesis of the
 * frame that plays it: a history whose mark says tone but which holds no
 * tone word plays its repeat as the same history unmarked, through the P25
 * and the ProVoice data paths, and leaves no such mark. The decoders' own
 * marks reach the history untouched (a frame with no voiced band at b0 60:
 * MBE_FRAME_MARK_PITCH_ABSENT).
 */
static void
test_forged_tone_mark(void) {
    char imbe_d[88];
    float logm[57], out[160], forged_out[160], plain_out[160];
    int K;
    for (int provoice = 0; provoice < 2; ++provoice) {
        struct stream s;
        stream_init(&s);
        const int L = harmonics_of(60, &K);
        tone_amplitudes(logm, L, 0, 0, 0.0f, 6.0f);
        build_frame(imbe_d, 60, 0x7Fu, logm, &s.prev);
        play(out, &s, imbe_d, provoice);
        assert(mbe_frame_mark(&s.prev) == MBE_FRAME_MARK_NONE && s.prev.Vl[0] == 0);

        struct stream forged = s, plain = s;
        mbe_set_frame_mark(&forged.prev, MBE_FRAME_MARK_TONE);
        play_repeat(forged_out, &forged, provoice);
        play_repeat(plain_out, &plain, provoice);
        assert(same_samples(forged_out, plain_out, 160));
        assert(mbe_frame_mark(&forged.enhanced) != MBE_FRAME_MARK_TONE);
        assert(mbe_frame_mark(&forged.prev) != MBE_FRAME_MARK_TONE);

        /* The decoder's own mark: no voiced band at a code outside 128..183. */
        build_frame(imbe_d, 60, 0x00u, logm, &s.prev);
        play(out, &s, imbe_d, provoice);
        assert(mbe_frame_mark(&s.prev) == MBE_FRAME_MARK_PITCH_ABSENT);
        assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_PITCH_ABSENT);
    }
    puts("p25 tones: the tone mark is the tone model's alone");
}

/*
 * The public parameter decoder leaves 0 in Vl[0] on every frame that is not a
 * tone: a special frame (b0 255, return 1) decoded into a state that held a
 * tone word clears it there, and leaves the history it is given untouched (its
 * word included: the process API's repeat takes the tone from it).
 */
static void
test_public_decoder_clears_the_word(void) {
    char imbe_d[88];
    float logm[57];
    int K;
    mbe_parms history = fresh_history();
    tone_amplitudes(logm, harmonics_of(53, &K), 4, 7, 11.0f, 2.0f);
    build_frame(imbe_d, 53, 0x70u, logm, &history);
    mbe_parms cur = decode(imbe_d, &history);
    assert(cur.Vl[0] == tone_word(53, 4, 7));
    mbe_parms prev = cur;
    const mbe_parms before = prev;

    special_frame(imbe_d);
    const int rc = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
    assert(rc == 1);
    (void)rc;
    assert(cur.Vl[0] == 0);
    assert(prev.Vl[0] == before.Vl[0] && prev.L == before.L && same_float(prev.w0, before.w0));
    for (int l = 0; l <= 56; ++l) {
        assert(prev.Vl[l] == before.Vl[l] && same_float(prev.Ml[l], before.Ml[l])
               && same_float(prev.log2Ml[l], before.log2Ml[l]));
    }
    puts("p25 tones: the public decoder clears the word on a special frame");
}

/*
 * ProVoice keeps no tone word: Vl[0] reads 0 in all three states after its
 * frames, a repeat of a history that holds a P25 tone word included (a stream
 * that changed codec), which plays no tone model either.
 */
static void
test_provoice_keeps_no_word(void) {
    struct stream pv;
    float out[160];
    stream_after_tone(&pv);
    play_repeat(out, &pv, 1);
    assert(mbe_frame_mark(&pv.enhanced) != MBE_FRAME_MARK_TONE);
    assert(pv.cur.Vl[0] == 0 && pv.prev.Vl[0] == 0 && pv.enhanced.Vl[0] == 0);
    puts("p25 tones: ProVoice keeps no word");
}

/*
 * A tone across the repeat limit: three repeats of a tone frame play the tone
 * model (marked, at the tone pitch, not muted), the fourth mutes (the noise
 * generator still runs), and the fifth resets the model, which clears the word
 * and the mark.
 */
static void
test_tone_repeat_limit(void) {
    struct stream s;
    char imbe_d[88];
    float out[160];
    mbe_process_result result;
    stream_after_tone(&s);
    special_frame(imbe_d);
    for (int r = 1; r <= 5; ++r) {
        const float seed = s.enhanced.noiseSeed;
        mbe_initProcessResult(&result);
        const int rc = mbe_processImbe4400Dataf(out, &result, imbe_d, &s.cur, &s.prev, &s.enhanced);
        assert(rc >= 0);
        (void)rc;
        assert((result.flags & MBE_PROCESS_FLAG_REPEAT) != 0u);
        const int muted = (result.flags & MBE_PROCESS_FLAG_MUTE) != 0u;
        if (r <= 3) {
            assert(!muted && s.enhanced.Vl[0] == tone_word(53, 4, 7));
            assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_TONE && same_float(s.enhanced.w0, tone_pitch(53)));
        } else if (r == 4) {
            assert(muted && s.enhanced.repeatCount == 4);
            assert(!same_float(s.enhanced.noiseSeed, seed));
        } else {
            assert(s.cur.Vl[0] == 0 && s.prev.Vl[0] == 0 && s.enhanced.Vl[0] == 0);
            assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_NONE);
        }
    }
    puts("p25 tones: three repeats play the tone, the fourth mutes, the reset clears the word and the mark");
}

/* FNV-1a over 32-bit values, each as four bytes from the least significant:
 * a step is a bijection of the state, so a change of any one value changes
 * the result. */
static uint32_t
fnv1a_add(uint32_t hash, uint32_t value) {
    for (int byte = 0; byte < 4; ++byte) {
        hash ^= (value >> (8 * byte)) & 0xFFu;
        hash *= UINT32_C(16777619);
    }
    return hash;
}

/*
 * The P25 tone table pinned: its 116 rows (b0, b1, l1, l2), sorted by b0, each
 * with 1 <= l1 <= l2 <= L of b0 and a voicing word of b0's K bits, the rows of
 * one (b0, b1) pair all dual or all single tones; the index derived from the
 * rows (the first row of each code 0..197 and the end); and an FNV-1a hash of
 * the rows' fields in order and then of the 199 index entries (0x93A8CEBD),
 * which any one-unit change of a field or an entry moves.
 */
#define P25_TONE_TABLE_HASH UINT32_C(0x93A8CEBD)

static void
test_tone_table_pinned(void) {
    uint32_t hash = UINT32_C(2166136261);
    int rows = 0, last_b0 = -1;
    int b0, b1, l1, l2;
    while (mbe_imbe4400_tone_row(rows, &b0, &b1, &l1, &l2)) {
        int K;
        const int L = harmonics_of(b0, &K);
        assert(b0 >= last_b0 && b0 <= 197);
        assert(l1 >= 1 && l1 <= l2 && l2 <= L && b1 > 0 && b1 < (1 << K));
        for (int other = 0; other < rows; ++other) {
            int ob0, ob1, ol1, ol2;
            const int found = mbe_imbe4400_tone_row(other, &ob0, &ob1, &ol1, &ol2);
            assert(found);
            assert(ob0 != b0 || ob1 != b1 || ((ol1 == ol2) == (l1 == l2)));
            assert(ob0 != b0 || ob1 != b1 || ol1 != l1 || ol2 != l2);
        }
        hash = fnv1a_add(hash, (uint32_t)b0);
        hash = fnv1a_add(hash, (uint32_t)b1);
        hash = fnv1a_add(hash, (uint32_t)l1);
        hash = fnv1a_add(hash, (uint32_t)l2);
        last_b0 = b0;
        ++rows;
    }
    assert(rows == 116);
    for (int code = 0; code <= 198; ++code) {
        int below = 0;
        for (int i = 0; i < rows; ++i) {
            const int found = mbe_imbe4400_tone_row(i, &b0, NULL, NULL, NULL);
            assert(found);
            below += (b0 < code) ? 1 : 0;
        }
        const int entry = mbe_imbe4400_tone_first_row(code);
        assert(entry == below);
        hash = fnv1a_add(hash, (uint32_t)entry);
    }
    assert(mbe_imbe4400_tone_first_row(-1) == -1 && mbe_imbe4400_tone_first_row(199) == -1);
    printf("p25 tones: the table's 116 rows and index pinned (hash 0x%08X)\n", (unsigned)hash);
    assert(hash == P25_TONE_TABLE_HASH);
}

int
main(void) {
    test_tone_test();
    test_tone_model();
    test_tone_run_renders();
    test_impossible_tone_words();
    test_forged_tone_mark();
    test_public_decoder_clears_the_word();
    test_provoice_keeps_no_word();
    test_tone_repeat_limit();
    test_tone_table_pinned();
    return 0;
}
