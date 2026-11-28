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
 * tone's fundamental (Annex J's f0 for a dual tone, but DVSI's own for 145 and
 * 147; 31.25 ID / l for a single tone, l Annex J's harmonic number), with a
 * voicing word that is the same in every steady frame of the tone (12,555
 * steady frames compared; 112 frames at tone boundaries carry another pair)
 * and the tone's harmonics large; DVSI's decoder recognises such a frame and
 * plays a tone model. The frames here are built with the IMBE encoder's quantizer
 * (mbe_imbe4400_quantize_amplitudes()) from a pitch code, a voicing word and
 * log2 amplitudes, and decoded as the library decodes P25.
 */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "imbe4400_internal.h"
#include "mbe_adaptive.h"
#include "mbe_p25_tone_model.h"
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
 * unvoiced at 0, L and the decoded log2Ml kept, marked as a tone, the
 * low cut's tracker left as it stood; the prediction history keeps the
 * decoded frame. The mid-frame model of a marked
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
    /* The tone model leaves DVSI's low-cut tracker where a fresh stream has
     * it, 110 Hz (its mark, and two voiced bands only). */
    assert(same_float(model->noiseOverlap[MBE_LOW_CUT_TRACK_SLOT], 110.0f));
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
    mbe_subframe_mid_model(&last, model, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE, &mid);
    mbe_subframe_own_model(model, &own);
    assert(same_subframe_model(&mid, &own));
    /* Advancing the marked frame unrendered leaves the train rendering it leaves. */
    mbe_parms rendered = *model, advanced = *model;
    float scratch[160] = {0.0f};
    mbe_synthesize_subframes(scratch, &rendered, &last, NULL, 1.0f, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
    mbe_advance_subframes(&advanced, &last, 0, MBE_VOICED_SYNTHESIS_DVSI_IMBE);
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

/* A stream after one frame at pitch code b0 with voicing word `word`, log2
 * amplitude `tone` at harmonics l1 and l2 and `rest` elsewhere, played through
 * the P25 data path: its history holds the decoded frame. */
static void
stream_after_frame(struct stream* s, int b0, unsigned word, int l1, int l2, float tone, float rest) {
    char imbe_d[88];
    float logm[57], out[160];
    int K;
    stream_init(s);
    tone_amplitudes(logm, harmonics_of(b0, &K), l1, l2, tone, rest);
    build_frame(imbe_d, b0, word, logm, &s->prev);
    play(out, s, imbe_d, 0);
}

/* A stream after one DTMF 1 frame (b0 53, word 1110000, harmonics 4 and 7):
 * its history holds the decoded frame and the tone word. */
static void
stream_after_tone(struct stream* s) {
    stream_after_frame(s, 53, 0x70u, 4, 7, 11.0f, 2.0f);
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
 * A word in Vl[0] counts as a tone only as a decode leaves it (the decode's
 * test run again on the model gives it: a pitch code with rows, the model at
 * that code's pitch, harmonics that are one of that code's rows, tone
 * magnitudes positive, finite and normal). A repeat of a history holding any
 * other word plays as the same history without it: a code past the table
 * (207, with a loud harmonic 1),
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

/* A repeat of `s` with the word `word` in its history plays as the same
 * history without a word: unmarked, sample for sample alike. */
static void
assert_word_plays_as_none(const struct stream* s, int word) {
    float with_word[160], without[160];
    struct stream forged = *s, plain = *s;
    forged.prev.Vl[0] = word;
    plain.prev.Vl[0] = 0;
    play_repeat(with_word, &forged, 0);
    play_repeat(without, &plain, 0);
    assert(mbe_frame_mark(&forged.enhanced) != MBE_FRAME_MARK_TONE);
    assert(same_samples(with_word, without, 160));
    assert(same_float(forged.enhanced.w0, plain.enhanced.w0));
}

/*
 * A word counts only if the decode's own test, run again on the model it sits
 * in, gives exactly that word: the model at the pitch its code decodes to
 * (w0 by its bits, L and K), and the tone test on its voicing and decoded
 * magnitudes. Codes 1 and 2 have the same L (9), K (3) and a row at harmonic 7
 * with the voicing word 001: code 1's word in a real code-2 tone frame's
 * history plays as no word, and so does code 1's word in a real code-1 frame
 * of that voicing whose other magnitudes are too large for a single tone (no
 * word decoded). The real words play the tone.
 */
static void
test_word_from_the_model(void) {
    int K1, K2;
    const int L1 = harmonics_of(1, &K1);
    const int L2 = harmonics_of(2, &K2);
    assert(L1 == 9 && L2 == 9 && K1 == 3 && K2 == 3);
    (void)L1;
    (void)L2;

    struct stream code2;
    stream_after_frame(&code2, 2, 0x001u, 7, 7, 11.0f, 2.0f);
    assert(code2.prev.Vl[0] == tone_word(2, 7, 7) && mbe_frame_mark(&code2.enhanced) == MBE_FRAME_MARK_TONE);
    assert_word_plays_as_none(&code2, tone_word(1, 7, 7));

    struct stream speech;
    stream_after_frame(&speech, 1, 0x001u, 7, 7, 11.0f, 8.0f);
    assert(speech.prev.Vl[0] == 0 && others_over_tone(&speech.prev, 7, 7) > 0.0235);
    assert(mbe_frame_mark(&speech.enhanced) != MBE_FRAME_MARK_TONE);
    assert_word_plays_as_none(&speech, tone_word(1, 7, 7));

    /* The real words: a code-1 tone frame plays the tone, and its repeat. */
    struct stream code1;
    float out[160];
    stream_after_frame(&code1, 1, 0x001u, 7, 7, 11.0f, 2.0f);
    assert(code1.prev.Vl[0] == tone_word(1, 7, 7) && mbe_frame_mark(&code1.enhanced) == MBE_FRAME_MARK_TONE);
    play_repeat(out, &code1, 0);
    assert(mbe_frame_mark(&code1.enhanced) == MBE_FRAME_MARK_TONE && same_float(code1.enhanced.w0, tone_pitch(1)));
    play_repeat(out, &code2, 0);
    assert(mbe_frame_mark(&code2.enhanced) == MBE_FRAME_MARK_TONE && same_float(code2.enhanced.w0, tone_pitch(2)));
    puts("p25 tones: a word counts only if the model's own test gives it");
}

/*
 * Every row of the table through a real decode: ten frames of the row's pitch
 * code and voicing word with its tone harmonics loud play the tone model once
 * the prediction has settled (from the fourth), and three repeats after them
 * play it again, at the row's tone pitch and harmonics.
 */
static void
test_every_row_plays(void) {
    int i = 0, b0, b1, l1, l2;
    while (mbe_imbe4400_tone_row(i, &b0, &b1, &l1, &l2)) {
        struct stream s;
        char imbe_d[88];
        float logm[57], out[160];
        int K;
        stream_init(&s);
        tone_amplitudes(logm, harmonics_of(b0, &K), l1, l2, 11.0f, 2.0f);
        for (int f = 0; f < 13; ++f) {
            if (f < 10) {
                build_frame(imbe_d, b0, (unsigned)b1, logm, &s.prev);
                play(out, &s, imbe_d, 0);
            } else {
                play_repeat(out, &s, 0);
            }
            if (f >= 3) {
                assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_TONE);
                assert(s.prev.Vl[0] == tone_word(b0, l1, l2) && same_float(s.enhanced.w0, tone_pitch(b0)));
            }
        }
        ++i;
    }
    assert(i == 116);
    puts("p25 tones: every row plays the tone through a real decode and three repeats");
}

/*
 * Only the tone model is marked MBE_FRAME_MARK_TONE, by the synthesis of the
 * frame that plays it: a history whose mark says tone but which holds no
 * tone word plays its repeat as the same history unmarked, through the P25
 * and the ProVoice data paths. A voiced model then carries the mark a voiced
 * model rendered carries (MBE_FRAME_MARK_ONSET_COLUMNS); a model with no voiced
 * band carries none, in all three states, where its decode had left
 * MBE_FRAME_MARK_PITCH_ABSENT (a frame with no voiced band at b0 60): the
 * decoded mark the caller's state overwrote is not recovered. The decoders' own
 * marks reach the history untouched.
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
        assert(mbe_frame_mark(&s.prev) == MBE_FRAME_MARK_ONSET_COLUMNS && s.prev.Vl[0] == 0);

        struct stream forged = s, plain = s;
        mbe_set_frame_mark(&forged.prev, MBE_FRAME_MARK_TONE);
        play_repeat(forged_out, &forged, provoice);
        play_repeat(plain_out, &plain, provoice);
        assert(same_samples(forged_out, plain_out, 160));
        /* The forged mark gives way to the one a decode leaves a voiced model. */
        assert(mbe_frame_mark(&forged.enhanced) == MBE_FRAME_MARK_ONSET_COLUMNS);
        assert(mbe_frame_mark(&forged.prev) == MBE_FRAME_MARK_ONSET_COLUMNS);
        assert(mbe_frame_mark(&plain.enhanced) == MBE_FRAME_MARK_ONSET_COLUMNS);

        /* The decoder's own mark: no voiced band at a code outside 128..183. */
        build_frame(imbe_d, 60, 0x00u, logm, &s.prev);
        play(out, &s, imbe_d, provoice);
        assert(mbe_frame_mark(&s.prev) == MBE_FRAME_MARK_PITCH_ABSENT);
        assert(mbe_frame_mark(&s.enhanced) == MBE_FRAME_MARK_PITCH_ABSENT);

        /* That history marked as a tone: a repeat plays it as unmarked. */
        struct stream unvoiced_forged = s, unvoiced_plain = s;
        mbe_set_frame_mark(&unvoiced_forged.prev, MBE_FRAME_MARK_TONE);
        mbe_set_frame_mark(&unvoiced_plain.prev, MBE_FRAME_MARK_NONE);
        play_repeat(forged_out, &unvoiced_forged, provoice);
        play_repeat(plain_out, &unvoiced_plain, provoice);
        assert(same_samples(forged_out, plain_out, 160));
        assert(mbe_frame_mark(&unvoiced_forged.cur) == MBE_FRAME_MARK_NONE);
        assert(mbe_frame_mark(&unvoiced_forged.prev) == MBE_FRAME_MARK_NONE);
        assert(mbe_frame_mark(&unvoiced_forged.enhanced) == MBE_FRAME_MARK_NONE);
        /* Left alone, the repeat keeps the decoded mark. */
        struct stream decoded = s;
        play_repeat(out, &decoded, provoice);
        assert(mbe_frame_mark(&decoded.prev) == MBE_FRAME_MARK_PITCH_ABSENT);
        assert(mbe_frame_mark(&decoded.enhanced) == MBE_FRAME_MARK_PITCH_ABSENT);
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

/* Two states alike byte for byte (copied to bytes, not compared as floats). */
static int
same_state(const mbe_parms* a, const mbe_parms* b) {
    unsigned char x[sizeof(mbe_parms)], y[sizeof(mbe_parms)];
    memcpy(x, a, sizeof(x));
    memcpy(y, b, sizeof(y));
    for (size_t i = 0; i < sizeof(x); ++i) {
        if (x[i] != y[i]) {
            return 0;
        }
    }
    return 1;
}

/*
 * A status error is not a frame: the public parameter decoder given a bit that
 * is not 0 or 1 after a tone frame returns MBE_STATUS_INVALID_BITS and writes
 * nothing, so cur_mp keeps the tone word it held and both states are
 * unchanged byte for byte; given no state it returns
 * MBE_STATUS_INVALID_ARGUMENT and writes nothing either.
 */
static void
test_public_decoder_status_error(void) {
    char imbe_d[88];
    float logm[57];
    int K;
    mbe_parms history = fresh_history();
    tone_amplitudes(logm, harmonics_of(53, &K), 4, 7, 11.0f, 2.0f);
    build_frame(imbe_d, 53, 0x70u, logm, &history);
    mbe_parms cur = decode(imbe_d, &history);
    assert(cur.Vl[0] == tone_word(53, 4, 7));
    mbe_parms prev = cur;
    const mbe_parms cur_before = cur, prev_before = prev;

    imbe_d[5] = 2;
    int rc = mbe_decodeImbe4400Parms(imbe_d, &cur, &prev);
    assert(rc == MBE_STATUS_INVALID_BITS);
    assert(same_state(&cur, &cur_before) && same_state(&prev, &prev_before) && cur.Vl[0] == tone_word(53, 4, 7));
    imbe_d[5] = 0;
    rc = mbe_decodeImbe4400Parms(imbe_d, &cur, NULL);
    assert(rc == MBE_STATUS_INVALID_ARGUMENT);
    assert(same_state(&cur, &cur_before));
    (void)rc;
    puts("p25 tones: a status error of the public decoder writes nothing");
}

/*
 * ProVoice keeps no tone word: Vl[0] reads 0 in all three states after its
 * frames, a repeat of a history that holds a P25 tone word included (a stream
 * moved between the codecs' entry points, outside the stream contract), which
 * plays no tone model either. A P25 repeat after that frame plays the repeated
 * model as speech: the decoded tone frame's own voicing at its P25 pitch,
 * unmarked, no word anywhere.
 */
static void
test_provoice_keeps_no_word(void) {
    struct stream pv;
    float out[160];
    stream_after_tone(&pv);
    const mbe_parms decoded = pv.prev;
    play_repeat(out, &pv, 1);
    assert(mbe_frame_mark(&pv.enhanced) != MBE_FRAME_MARK_TONE);
    assert(pv.cur.Vl[0] == 0 && pv.prev.Vl[0] == 0 && pv.enhanced.Vl[0] == 0);

    play_repeat(out, &pv, 0);
    assert(mbe_frame_mark(&pv.enhanced) != MBE_FRAME_MARK_TONE);
    assert(pv.cur.Vl[0] == 0 && pv.prev.Vl[0] == 0 && pv.enhanced.Vl[0] == 0);
    /* Speech plays the code's pitch, eq 46 on the 2^-19 grid, which is the
     * tone model's too: the voicing and the mark tell the two apart. */
    assert(same_float(pv.enhanced.w0, decoded.w0) && same_float(pv.enhanced.w0, tone_pitch(53)));
    assert(pv.enhanced.L == decoded.L);
    for (int l = 1; l <= decoded.L; ++l) {
        assert(pv.enhanced.Vl[l] == decoded.Vl[l]);
    }
    puts("p25 tones: ProVoice keeps no word; a P25 repeat after it plays speech");
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

/* FNV-1a over 32-bit values, each as four bytes from the least significant.
 * Not every change of a value moves the result (two 32-bit values can hash
 * alike from one state); every one-unit change of any entry of the tables
 * pinned with it was checked to move it. */
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
 * the rows' fields in order and then of the 199 index entries (0x93A8CEBD):
 * every one-unit change of a field or an entry was checked to move it.
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

/* Whether two states are alike but for the tone word in Vl[0]. */
static int
same_state_but_word(const mbe_parms* a, const mbe_parms* b) {
    mbe_parms x = *a, y = *b;
    x.Vl[0] = 0;
    y.Vl[0] = 0;
    return same_state(&x, &y);
}

/* `s` with the tone harmonics' history magnitudes at `level` and every other
 * magnitude 1..L at 0. */
static void
history_at_level(struct stream* s, int l1, int l2, float level) {
    for (int l = 1; l <= s->prev.L; ++l) {
        s->prev.Ml[l] = (l == l1 || l == l2) ? level : 0.0f;
    }
}

/*
 * The levels the tone player writes are decided where the tone is taken: a
 * repeat of a real DTMF 1 history (b0 53, harmonics 4 and 7) whose two tone
 * magnitudes a caller set to FLT_MIN, the others 0, would be played at
 * subnormal levels (0.988 FLT_MIN times each droop, about 0.987 FLT_MIN),
 * which the predicate refuses: the tone is not taken, and the repeat plays
 * as the same history with the tone word cleared, sample for sample, the
 * three returned states alike but for the word itself (which a repeat
 * carries), unmarked. Over every one of the 116 rows and a ladder of levels
 * from FLT_MIN to 1e30, every model the player marks passes the predicate
 * (mbe_imbe4400_tone_model()), and every repeat not marked carries no tone
 * mark; both kinds occur. Above the renderer's magnitude limit (1e30 here
 * near it) the predicate passes and the rendering caps the magnitudes as it
 * caps any.
 */
static void
test_tone_levels_decided_once(void) {
    struct stream base;
    stream_after_tone(&base);
    struct stream forged = base;
    history_at_level(&forged, 4, 7, FLT_MIN);
    struct stream plain = forged;
    plain.prev.Vl[0] = 0;
    float with_word[160], without[160];
    play_repeat(with_word, &forged, 0);
    play_repeat(without, &plain, 0);
    assert(mbe_frame_mark(&forged.enhanced) != MBE_FRAME_MARK_TONE);
    assert(same_samples(with_word, without, 160));
    assert(same_state_but_word(&forged.cur, &plain.cur) && same_state_but_word(&forged.prev, &plain.prev)
           && same_state_but_word(&forged.enhanced, &plain.enhanced));

    int i = 0, b0, b1, l1, l2, marked = 0, refused = 0;
    while (mbe_imbe4400_tone_row(i, &b0, &b1, &l1, &l2)) {
        struct stream row;
        char imbe_d[88];
        float logm[57], out[160];
        int K;
        stream_init(&row);
        tone_amplitudes(logm, harmonics_of(b0, &K), l1, l2, 11.0f, 2.0f);
        for (int f = 0; f < 4; ++f) {
            build_frame(imbe_d, b0, (unsigned)b1, logm, &row.prev);
            play(out, &row, imbe_d, 0);
        }
        assert(row.prev.Vl[0] == tone_word(b0, l1, l2));
        for (int rung = 0; rung < 47; ++rung) {
            /* FLT_MIN times 31^rung (to about 1.5e29), then 1e30 */
            const double level = (rung < 46) ? (double)FLT_MIN * pow(31.0, (double)rung) : 1.0e30;
            for (int step = 0; step < 3; ++step) {
                /* FLT_MIN itself, and the two floats above it at each rung */
                const float at = (float)(level * (1.0 + (double)step * 1.0e-6));
                struct stream t = row;
                history_at_level(&t, l1, l2, at);
                play_repeat(out, &t, 0);
                const int plays = mbe_frame_mark(&t.enhanced) == MBE_FRAME_MARK_TONE;
                const int passes = mbe_imbe4400_tone_model(&t.enhanced);
                assert(plays == passes);
                marked += plays;
                refused += !plays;
            }
        }
        ++i;
    }
    assert(i == 116 && marked > 0 && refused > 0);
    printf("p25 tones: the levels decided once (%d played models pass the predicate, %d refused unmarked)\n", marked,
           refused);
}

int
main(void) {
    test_tone_levels_decided_once();
    test_tone_test();
    test_tone_model();
    test_tone_run_renders();
    test_impossible_tone_words();
    test_word_from_the_model();
    test_every_row_plays();
    test_forged_tone_mark();
    test_public_decoder_clears_the_word();
    test_public_decoder_status_error();
    test_provoice_keeps_no_word();
    test_tone_repeat_limit();
    test_tone_table_pinned();
    return 0;
}
