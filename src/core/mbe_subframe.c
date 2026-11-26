// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/** @file Speech synthesis in 10 ms subframes, as DVSI's decoders render it (mbe_subframe.h). */
#include "mbe_subframe.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "ambe3600x2450_internal.h"
#include "mbe_adaptive.h"
#include "mbe_math.h"
#include "mbe_validation.h"
#include "mbe_voiced.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FRAME           MBE_VOICED_FRAME
#define SUBFRAME        MBE_SUBFRAME
#define PEAK            MBE_PULSE_PEAK
/* Pulses one interval may start (a pitch period of 5 samples: no codec's
 * pitch comes near it) and pulses a frame's rendering can see. */
#define INTERVAL_PULSES 16
#define TRAIN_PULSES    (MBE_PULSE_HISTORY + 2 * INTERVAL_PULSES)
/* Band states of a model. */
#define BAND_UNVOICED   0
#define BAND_VOICED     1
#define BAND_NONE       2
/* log2 magnitude of a band with no amplitude, and the level under which an
 * interpolated magnitude counts as none. */
#define LOG2_SILENT     (-1000.0)
#define LOG2_FLOOR      (-200.0)
/* The largest magnitude a model takes. */
#define MAGNITUDE_LIMIT 1.0e30f
/* Where noiseOverlap holds the pulse phase, after the count and the pulses. */
#define PHASE_SLOT      (1 + MBE_PULSE_HISTORY)
/* The bound on a log2 amplitude the phases are regenerated from: far past
 * any codec's (D-STAR's decoder reaches about 49), where a caller's state can
 * put anything. */
#define PHASE_LOG2_MAX  1000.0f

/* A positive, finite float that is not subnormal, read from its bit pattern
 * so that fast-math cannot fold the test into accepting NaN (see
 * mbe_encoder.c). Subnormal values are left out: where the caller's
 * floating-point environment flushes them to zero, as audio hosts set it,
 * they are zero to the arithmetic, and a pitch of zero divides into positions
 * that are not numbers. */
static int
positive_normal(const float* x) {
    uint32_t bits;
    memcpy(&bits, x, sizeof(bits));
    return bits - UINT32_C(0x00800000) < UINT32_C(0x7F000000);
}

static int
any_finite(const float* x) {
    uint32_t bits;
    memcpy(&bits, x, sizeof(bits));
    return (bits & 0x7FFFFFFFu) < 0x7F800000u;
}

/* A model's magnitude, kept where sums of 56 harmonics stay in the float
 * range: no codec's amplitudes come near the limit, but a caller's state can
 * hold anything. Takes a positive finite value or +infinity. */
static float
limited_magnitude(float magnitude) {
    const float limit = MAGNITUDE_LIMIT;
    uint32_t bits, limit_bits;
    memcpy(&bits, &magnitude, sizeof(bits));
    memcpy(&limit_bits, &limit, sizeof(limit_bits));
    return (bits > limit_bits) ? limit : magnitude;
}

static int
band_state(const mbe_parms* mp, int l) {
    if (mp->Vl[l] == 0) {
        return BAND_UNVOICED;
    }
    return (mp->Vl[l] == 1) ? BAND_VOICED : BAND_NONE;
}

static void
note_voiced(struct mbe_subframe_model* model) {
    model->any_voiced = 0;
    for (int l = 1; l <= model->L; l++) {
        if (model->voiced[l] == BAND_VOICED) {
            model->any_voiced = 1;
        }
    }
}

/* A frame's decoded log2 amplitude (log2Ml), which its harmonic phases are
 * regenerated from: not a finite number, the phase kernel's floor; finite,
 * kept within +-PHASE_LOG2_MAX. Read from the bit pattern (see
 * positive_normal()). */
static double
phase_log2_amplitude(const float* x) {
    const float limit = PHASE_LOG2_MAX;
    uint32_t bits, limit_bits;
    memcpy(&bits, x, sizeof(bits));
    memcpy(&limit_bits, &limit, sizeof(limit_bits));
    const uint32_t magnitude = bits & UINT32_C(0x7FFFFFFF);
    if (magnitude >= UINT32_C(0x7F800000)) {
        return (double)MBE_PHASE_LOG2_FLOOR;
    }
    if (magnitude > limit_bits) {
        return (bits >> 31) ? -(double)limit : (double)limit;
    }
    return (double)*x;
}

/* The log2 amplitudes a frame's phases come from, [1..L]. DVSI's decoders
 * regenerate the phase of a voiced harmonic from the frame's decoded log2
 * spectral amplitudes, before spectral enhancement, the high-band shelf and
 * any scaling, whatever each band's voicing, not from the amplitudes they
 * render (as US 5,701,390 applies its kernel to the decoded magnitudes;
 * docs/dvsi-parity.md, "Phase regeneration"). */
static void
phase_log2_amplitudes(const mbe_parms* mp, double out[57]) {
    for (int l = 1; l <= mp->L; l++) {
        out[l] = phase_log2_amplitude(&mp->log2Ml[l]);
    }
}

void
mbe_subframe_own_model(const mbe_parms* mp, struct mbe_subframe_model* model) {
    double log2_amplitudes[57];
    float envelope[57] = {0};

    memset(model, 0, sizeof(*model));
    if (!mbe_harmonic_count_is_valid(mp->L)) {
        return;
    }
    model->L = mp->L;
    if (!positive_normal(&mp->w0)) {
        for (int l = 1; l <= model->L; l++) {
            model->voiced[l] = BAND_NONE;
        }
        return;
    }
    model->w0 = mp->w0;
    for (int l = 1; l <= mp->L; l++) {
        model->voiced[l] = (unsigned char)band_state(mp, l);
        model->magnitude[l] = positive_normal(&mp->Ml[l]) ? limited_magnitude(mp->Ml[l]) : 0.0f;
    }
    note_voiced(model);
    phase_log2_amplitudes(mp, log2_amplitudes);
    for (int l = 1; l <= mp->L; l++) {
        envelope[l] = (float)log2_amplitudes[l];
    }
    mbe_regenerate_phases_log2(envelope, model->L, model->phase);
}

/* log2 of each band's magnitude with the unvoiced factor taken out. */
static void
common_log2_magnitudes(const mbe_parms* mp, double out[57]) {
    const double unvoiced = log2((double)MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt((double)mp->w0));
    for (int l = 1; l <= mp->L; l++) {
        if (!positive_normal(&mp->Ml[l])) {
            out[l] = LOG2_SILENT;
        } else {
            out[l] = log2((double)mp->Ml[l]) - ((mp->Vl[l] == 0) ? unvoiced : 0.0);
        }
    }
}

/* The harmonic at or below a position among a frame's L harmonics: the whole
 * part of a position from 1 up to L, with *within set; else the first
 * harmonic below them and the last at or past it. Decided on the bit pattern,
 * so that the conversion only ever sees a value in range: a position that is
 * not a number (fast-math lets float comparisons pass one) takes an end. */
static int
harmonic_below(float position, int L, int* within) {
    const float last = (float)L;
    uint32_t bits, last_bits;
    memcpy(&bits, &position, sizeof(bits));
    memcpy(&last_bits, &last, sizeof(last_bits));
    /* Positive floats order as their bit patterns do; 1.0f is 0x3F800000. */
    *within = (bits >= UINT32_C(0x3F800000)) && (bits < last_bits);
    if (*within) {
        return (int)position;
    }
    /* Up to +infinity (0x7F800000) from L on: the last. */
    return ((bits >= last_bits) && (bits <= UINT32_C(0x7F800000))) ? L : 1;
}

/* A frame's harmonic nearest angular frequency w, 1..L. */
static int
harmonic_at(const mbe_parms* mp, float w) {
    int within;
    return harmonic_below((w / mp->w0) + 0.5f, mp->L, &within);
}

/* A frame's log2 magnitude at angular frequency w: linear between the two
 * harmonics around it, the end harmonic's beyond either end. */
static double
log2_magnitude_at(const mbe_parms* mp, const double log2_magnitude[57], float w) {
    const float position = w / mp->w0;
    int within;
    const int below = harmonic_below(position, mp->L, &within);
    if (!within) {
        return log2_magnitude[below];
    }
    const double t = (double)position - (double)below;
    return log2_magnitude[below] + (t * (log2_magnitude[below + 1] - log2_magnitude[below]));
}

static int
model_has_voiced_band(const mbe_parms* mp) {
    for (int l = 1; l <= mp->L; l++) {
        if (mp->Vl[l] == 1) {
            return 1;
        }
    }
    return 0;
}

/* Whether two pitches are the same, read from their bit patterns. */
static int
same_pitch(const float* a, const float* b) {
    uint32_t a_bits, b_bits;
    memcpy(&a_bits, a, sizeof(a_bits));
    memcpy(&b_bits, b, sizeof(b_bits));
    return a_bits == b_bits;
}

/* The half-octave wave h(P) of rate 33's mid-frame pitch at a frame's pitch
 * w0 (positive and normal), P = 2 pi / w0 its period in samples
 * (MBE_AMBE2450_MID_PITCH_WAVE_C0, mbe_subframe.h). In double with the
 * double-precision functions, as the P25 pitch law is (imbe7200x4400.c), so
 * that the pitch does not move with a libm's single-precision ones. */
static double
ambe2450_mid_pitch_wave(float w0) {
    const double t = 4.0 * M_PI * log2((2.0 * M_PI) / (double)w0);
    return MBE_AMBE2450_MID_PITCH_WAVE_C0 + (MBE_AMBE2450_MID_PITCH_WAVE_C1 * cos(t))
           + (MBE_AMBE2450_MID_PITCH_WAVE_S1 * sin(t));
}

/* Rate 33's mid-frame pitch between two voiced frames of different pitch a
 * and b (each positive and normal): their geometric mean times
 * MBE_AMBE2450_MID_PITCH_FACTOR times 1 + (h(Pa) + h(Pb)) / 2, formed in
 * double from the stored floats and rounded to float once.
 * Its range: the product of two floats is exact in double and the mean lies
 * between a and b, so it is finite, positive and normal as they are; h lies
 * within C0 -+ sqrt(C1^2 + S1^2), 4.244e-5 -+ 8.440e-4, so the two factors
 * together lie between 0.99889 and 1.00059 (cos and sin are within [-1, 1]
 * for any finite t; t is finite for every positive normal w0, within
 * 4 pi (2.65 - 128) .. 4 pi (2.65 + 126)). The pitch is then positive and
 * finite in double. A mean so near FLT_MAX that the product passes it gives
 * 0 (the conversion of a double past the float range is not defined), and one
 * so near FLT_MIN that the product falls under it rounds to a subnormal float
 * or 0: the caller rejects either as no pitch (positive_normal()), as it
 * rejects every pitch that is not positive and normal. */
static float
ambe2450_mid_pitch(float a, float b) {
    /* The factor as the float it is declared, also under x87 excess
     * precision, where its constant expression is otherwise evaluated wider. */
    const float factor = mbe_round_to_float(MBE_AMBE2450_MID_PITCH_FACTOR);
    const double mean = sqrt((double)a * (double)b);
    const double wave = 0.5 * (ambe2450_mid_pitch_wave(a) + ambe2450_mid_pitch_wave(b));
    const double pitch = mean * (double)factor * (1.0 + wave);
    return (pitch <= (double)FLT_MAX) ? (float)pitch : 0.0f;
}

/* The mid-frame pitch: between the two frames' when both are pitched speech
 * or neither is, the voiced frame's when only one has a voiced band, unless
 * the unvoiced frame before a voiced one is marked as carrying its pitch
 * (MBE_FRAME_MARK_PITCH_CARRIES).
 * Between two voiced frames of different pitch, the rate-33 rendering plays
 * the geometric mean times MBE_AMBE2450_MID_PITCH_FACTOR times the half-octave
 * wave (ambe2450_mid_pitch(); docs/dvsi-parity.md, "Rate 33: mid-frame
 * pitch"), a fit to DVSI's rate-33 decodes (mbe_subframe.h): DVSI's pulse
 * train moves against ours by about 0.05 samples over a frame whose pitch
 * code differs from the frame before's, and by the estimator's noise
 * (0.018) over one whose code is the same, and the offset follows the two
 * frames' periods with a period of half an octave. How the decoder comes to
 * it is not known. D-STAR and P25 show no such jump at a step. Both pitches
 * are positive and normal here (mbe_subframe_mid_model()). */
static float
mid_pitch(const mbe_parms* prev_mp, const mbe_parms* cur_mp, enum mbe_voiced_synthesis rendering) {
    const int prev_voiced = model_has_voiced_band(prev_mp);
    const int cur_voiced = model_has_voiced_band(cur_mp);
    if (prev_voiced && !cur_voiced) {
        return prev_mp->w0;
    }
    if (cur_voiced && !prev_voiced && mbe_frame_mark(prev_mp) != MBE_FRAME_MARK_PITCH_CARRIES) {
        return cur_mp->w0;
    }
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 && prev_voiced && cur_voiced
        && !same_pitch(&prev_mp->w0, &cur_mp->w0)) {
        return ambe2450_mid_pitch(prev_mp->w0, cur_mp->w0);
    }
    return sqrtf(prev_mp->w0 * cur_mp->w0);
}

/* Whether the mid-frame model keeps the previous frame's voicing instead of
 * taking this frame's: DVSI's decoders hold it through the first 10 ms of a
 * frame with no voiced band that the decoder marked MBE_FRAME_MARK_PITCH_ABSENT
 * when the previous frame's voicing reaches MBE_VOICING_REACH_HZ, read at its
 * harmonic at or next above that frequency. */
static int
voicing_persists(const mbe_parms* prev_mp, const mbe_parms* cur_mp) {
    if (mbe_frame_mark(cur_mp) != MBE_FRAME_MARK_PITCH_ABSENT || model_has_voiced_band(cur_mp)) {
        return 0;
    }
    const float position = ceilf((MBE_VOICING_REACH_HZ * ((2.0f * (float)M_PI) / 8000.0f)) / prev_mp->w0);
    int within;
    return prev_mp->Vl[harmonic_below(position, prev_mp->L, &within)] == 1;
}

/* ---- the harmonics a P25 onset starts with ---- */

/* The 500 Hz columns a P25 frame's voicing is read in after a frame with no
 * voiced band (MBE_FRAME_MARK_ONSET_COLUMNS; a rule fitted to DVSI's P25
 * decodes, docs/dvsi-parity.md): column 0 from 0 Hz
 * to MBE_VOICING_REACH_HZ, each next one ONSET_COLUMN_HZ wide, ONSET_COLUMNS of
 * them, up to 3505 Hz; a harmonic within ONSET_REACH_HZ below a covered column
 * sounds with it. */
#define ONSET_COLUMNS   7
#define ONSET_COLUMN_HZ 500.0
#define ONSET_REACH_HZ  250.0

struct onset_columns {
    double edge[ONSET_COLUMNS];     /* column k's upper edge, in harmonics of the frame's pitch */
    double reach[ONSET_COLUMNS];    /* the same edges ONSET_REACH_HZ lower */
    int covered[ONSET_COLUMNS + 1]; /* column k covered; [ONSET_COLUMNS], above 3505 Hz, never */
};

/* Whether a frame is voiced at its harmonic at or next above `position` (in
 * harmonics of its pitch, positive and finite): unvoiced past L. Decided on
 * the bit pattern, positive doubles ordering as their patterns do, so that the
 * conversion only ever sees a value from 0 to L. */
static int
voiced_from(const mbe_parms* mp, double position) {
    const double last = (double)mp->L;
    uint64_t bits, last_bits;
    memcpy(&bits, &position, sizeof(bits));
    memcpy(&last_bits, &last, sizeof(last_bits));
    if (bits > last_bits) {
        return 0;
    }
    const int l = (int)ceil(position);
    return mp->Vl[(l < 1) ? 1 : l] == 1;
}

/* The columns of a frame with a positive, normal pitch and 1..56 harmonics.
 * The edges' positions among its harmonics are quotients formed in double
 * from the stored pitch: for every P25 pitch code (and ProVoice's nominal
 * pitch) each lies at least 1.5e-4 of a harmonic from a whole number
 * (test_p25_onset_edges in test_subframe_synthesis.c), far past any
 * difference the arithmetic of a build (x87 excess precision, fast-math's
 * reordering) can make, so every build reads the same harmonics. */
static void
onset_columns(const mbe_parms* mp, struct onset_columns* columns) {
    const double harmonics_per_hz = M_PI / (4000.0 * (double)mp->w0);
    int lower = mp->Vl[1] == 1;
    for (int k = 0; k < ONSET_COLUMNS; k++) {
        const double hz = (double)MBE_VOICING_REACH_HZ + (ONSET_COLUMN_HZ * (double)k);
        columns->edge[k] = hz * harmonics_per_hz;
        columns->reach[k] = (hz - ONSET_REACH_HZ) * harmonics_per_hz;
        const int upper = voiced_from(mp, columns->edge[k]);
        columns->covered[k] = lower && upper;
        lower = upper;
    }
    columns->covered[ONSET_COLUMNS] = 0;
}

/* The column holding a position among the frame's harmonics, by `edges`. */
static int
column_of(const double edges[ONSET_COLUMNS], double position) {
    int k = 0;
    while (k < ONSET_COLUMNS && position >= edges[k]) {
        k++;
    }
    return k;
}

/* Whether the mid-frame model's harmonic at `position` (in harmonics of the
 * frame's pitch: a whole number at an onset, where the mid-frame pitch is the
 * frame's) sounds at an onset: in a covered column, or within ONSET_REACH_HZ
 * below one. */
static int
sounds_at_onset(const struct onset_columns* columns, double position) {
    return columns->covered[column_of(columns->edge, position)]
           || columns->covered[column_of(columns->reach, position)];
}

/* The columns of the current frame when its onset is read in them (in P25's
 * rendering, a frame marked MBE_FRAME_MARK_ONSET_COLUMNS after one with no
 * voiced band), else NULL. */
static const struct onset_columns*
onset_rule(const mbe_parms* prev_mp, const mbe_parms* cur_mp, enum mbe_voiced_synthesis rendering,
           struct onset_columns* columns) {
    if (rendering != MBE_VOICED_SYNTHESIS_DVSI_IMBE || mbe_frame_mark(cur_mp) != MBE_FRAME_MARK_ONSET_COLUMNS
        || model_has_voiced_band(prev_mp)) {
        return NULL;
    }
    onset_columns(cur_mp, columns);
    return columns;
}

/* ---- the mid-frame voicing of rate 33's rows 13 and 17 ---- */

/* The 500 Hz columns the mid-frame model of a frame so marked voices, as V/UV
 * rows (mbe_ambe2450_harmonic_voiced()). */
static const int mid_column_0[8] = {1, 0, 0, 0, 0, 0, 0, 0};
static const int mid_columns_0_1[8] = {1, 1, 0, 0, 0, 0, 0, 0};

/* The mid-frame row of the current frame in the rate-33 rendering
 * (MBE_FRAME_MARK_MID_COLUMN_0, MBE_FRAME_MARK_MID_COLUMNS_0_1), or NULL for a
 * frame that takes its voicing at mid-frame from the frames around it. */
static const int*
mid_row(const mbe_parms* cur_mp, enum mbe_voiced_synthesis rendering) {
    if (rendering != MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) {
        return NULL;
    }
    switch (mbe_frame_mark(cur_mp)) {
        case MBE_FRAME_MARK_MID_COLUMN_0: return mid_column_0;
        case MBE_FRAME_MARK_MID_COLUMNS_0_1: return mid_columns_0_1;
        default: return NULL;
    }
}

/* Whether a row voices harmonic l of the mid-frame model, at its own pitch,
 * by the rate-33 decoder's rule. A harmonic at or above Nyquist (a caller's
 * pitch can put one there, and its column past any index) lies in the top
 * column, which neither row voices; the test is on the bit pattern of the
 * positive (or infinite) frequency, so that fast-math cannot pass one by. */
static int
row_voices(const int row[8], const struct mbe_subframe_model* model, int l) {
    const float w = model->w0 * (float)l;
    const float nyquist = (float)M_PI;
    uint32_t bits, nyquist_bits;
    memcpy(&bits, &w, sizeof(bits));
    memcpy(&nyquist_bits, &nyquist, sizeof(nyquist_bits));
    if (bits >= nyquist_bits) {
        return 0;
    }
    return mbe_ambe2450_harmonic_voiced(row, l, model->w0 / (2.0f * (float)M_PI));
}

/* The mid-frame model's voicing, magnitudes and phases, its pitch and harmonic
 * count being set: each band's state from the frame whose voicing holds at
 * mid-frame, and the two frames' magnitudes at the band's frequency, met in
 * the log domain. At an onset of a frame marked MBE_FRAME_MARK_ONSET_COLUMNS
 * (`onset`), a voiced band that does not sound there (sounds_at_onset()) is
 * unvoiced; a frame with a mid-frame row (`row`, mid_row()) takes each band's
 * voicing from the row (a band with no state keeps none). Its phases come
 * from the mean of the two frames' decoded
 * log2 amplitudes (phase_log2_amplitudes()), each read at the band's
 * frequency as the magnitudes are and used as the decoder stored it: where
 * common_log2_magnitudes() takes the unvoiced factor out of an unvoiced band's
 * Ml, nothing is taken out of or put into its log2Ml, which carries no such
 * factor (an AMBE decoder builds an ordinary voice frame's unvoiced band as
 * Ml = MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt(w0) times 2^log2Ml; Ml then also
 * carries what log2Ml does not, such as rate 33's further silence scale and
 * spectral enhancement). */
static void
fill_mid_bands(const mbe_parms* prev_mp, const mbe_parms* cur_mp, int late_onset, const struct onset_columns* onset,
               const int* row, struct mbe_subframe_model* model) {
    const double unvoiced = log2((double)MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt((double)model->w0));
    const mbe_parms* voicing = voicing_persists(prev_mp, cur_mp) ? prev_mp : cur_mp;
    double prev_log2[57], cur_log2[57], prev_phase_log2[57], cur_phase_log2[57];
    float envelope[57] = {0};

    common_log2_magnitudes(prev_mp, prev_log2);
    common_log2_magnitudes(cur_mp, cur_log2);
    phase_log2_amplitudes(prev_mp, prev_phase_log2);
    phase_log2_amplitudes(cur_mp, cur_phase_log2);
    for (int l = 1; l <= model->L; l++) {
        const float w = model->w0 * (float)l;
        int state = band_state(voicing, harmonic_at(voicing, w));
        if (late_onset && state == BAND_VOICED) {
            state = BAND_UNVOICED;
        }
        if (onset != NULL && state == BAND_VOICED
            && !sounds_at_onset(onset, ((double)model->w0 * (double)l) / (double)cur_mp->w0)) {
            state = BAND_UNVOICED;
        }
        if (row != NULL && state != BAND_NONE) {
            state = row_voices(row, model, l) ? BAND_VOICED : BAND_UNVOICED;
        }
        const double mean = 0.5 * (log2_magnitude_at(prev_mp, prev_log2, w) + log2_magnitude_at(cur_mp, cur_log2, w));
        model->voiced[l] = (unsigned char)state;
        if (mean > LOG2_FLOOR) {
            /* Limited before it narrows: a caller's state can put it past what a float holds. */
            const double magnitude = exp2(mean + ((state == BAND_UNVOICED) ? unvoiced : 0.0));
            model->magnitude[l] = (magnitude < (double)MAGNITUDE_LIMIT) ? (float)magnitude : MAGNITUDE_LIMIT;
        }
        envelope[l] =
            (float)(0.5
                    * (log2_magnitude_at(prev_mp, prev_phase_log2, w) + log2_magnitude_at(cur_mp, cur_phase_log2, w)));
    }
    note_voiced(model);
    mbe_regenerate_phases_log2(envelope, model->L, model->phase);
}

void
mbe_subframe_mid_model(const mbe_parms* prev_mp, const mbe_parms* cur_mp, int late_onset,
                       enum mbe_voiced_synthesis rendering, struct mbe_subframe_model* model) {
    memset(model, 0, sizeof(*model));
    if (!mbe_harmonic_count_is_valid(prev_mp->L) || !mbe_harmonic_count_is_valid(cur_mp->L)) {
        return;
    }
    if (mbe_frame_mark(cur_mp) == MBE_FRAME_MARK_PREVIOUS_HELD) {
        /* The previous frame's model stands at mid-frame unchanged. */
        mbe_subframe_own_model(prev_mp, model);
        return;
    }
    if (mbe_frame_mark(cur_mp) == MBE_FRAME_MARK_TONE) {
        /* A P25 tone frame's own model stands at mid-frame: DVSI's decoder
         * brings a tone's amplitude from the last frame's to its own over the
         * first subframe and holds it over the second (a first tone frame is at
         * full level from sample 80; the frame after the last holds the tone
         * fading out over its first 80 samples). */
        mbe_subframe_own_model(cur_mp, model);
        return;
    }
    const int shorter = (prev_mp->L < cur_mp->L) ? prev_mp->L : cur_mp->L;
    const int longer = (prev_mp->L < cur_mp->L) ? cur_mp->L : prev_mp->L;
    model->L = shorter;
    for (int l = 1; l <= shorter; l++) {
        model->voiced[l] = BAND_NONE;
    }
    if (!positive_normal(&prev_mp->w0) || !positive_normal(&cur_mp->w0)) {
        return;
    }
    const float w0 = mid_pitch(prev_mp, cur_mp, rendering);
    if (!positive_normal(&w0)) {
        return;
    }
    const float harmonics = floorf((0.9254f * (float)M_PI) / w0);
    model->w0 = w0;
    model->L = (harmonics < (float)shorter) ? shorter : ((harmonics > (float)longer) ? longer : (int)harmonics);
    struct onset_columns columns;
    fill_mid_bands(prev_mp, cur_mp, late_onset, onset_rule(prev_mp, cur_mp, rendering, &columns),
                   mid_row(cur_mp, rendering), model);
}

/* ---- the pulse train ---- */

struct pulse_train {
    double phase;             /* cycles since the last pulse, at the point reached */
    int count;                /* pulses that can still sound */
    float time[TRAIN_PULSES]; /* in samples from the frame's start (negative: before it) */
};

/* The train a stream carries into a frame: the pulse phase reached at the last
 * frame's end and the pulses it kept. It runs on from a frame that has a
 * voiced band, and from one marked MBE_FRAME_MARK_PITCH_ABSENT or
 * MBE_FRAME_MARK_PREVIOUS_HELD, through whose first half the voicing before it
 * may have lasted; after any other frame it starts afresh, whatever a codec
 * path left in the state.
 * Values that are not what a frame stores (another rendering's noise overlap,
 * a caller's own) are dropped. */
static void
load_pulse_train(const mbe_parms* prev_mp, const struct mbe_subframe_model* prev, struct pulse_train* train) {
    const float stored = prev_mp->noiseOverlap[0];
    const float phase = prev_mp->noiseOverlap[PHASE_SLOT];
    uint32_t bits;

    train->phase = 0.0;
    train->count = 0;
    if (!mbe_subframe_state_is_marked(prev_mp)) {
        return;
    }
    const enum mbe_frame_mark mark = mbe_frame_mark(prev_mp);
    if (!prev->any_voiced && mark != MBE_FRAME_MARK_PITCH_ABSENT && mark != MBE_FRAME_MARK_PREVIOUS_HELD) {
        return;
    }
    memcpy(&bits, &stored, sizeof(bits));
    /* 32.0f is 0x42000000: a count of 0..MBE_PULSE_HISTORY. */
    if (bits > UINT32_C(0x42000000)) {
        return;
    }
    memcpy(&bits, &phase, sizeof(bits));
    /* 1.0f is 0x3F800000: a phase of 0 up to, not including, one cycle. */
    if (bits < UINT32_C(0x3F800000)) {
        train->phase = (double)phase;
    }
    const int count = (int)stored;
    for (int i = 0; i < count; i++) {
        const float before_end = prev_mp->noiseOverlap[1 + i];
        if (any_finite(&before_end) && before_end <= 0.0f && before_end >= -(float)FRAME) {
            train->time[train->count++] = before_end;
        }
    }
}

/* Pulses of one interval starting at `start`: the pitch moves linearly from
 * `from` to `to` cycles per sample across it, so the pulse phase at t into the
 * interval is phase + from t + (to - from) t^2 / (2 SUBFRAME), and a pulse
 * falls wherever that reaches a whole number. */
static void
interval_pulses(struct pulse_train* train, float start, double from, double to) {
    const double slope = (to - from) / (double)SUBFRAME;
    int cycle = (train->phase > 0.0) ? 1 : 0;

    for (int pulses = 0; pulses < INTERVAL_PULSES; pulses++, cycle++) {
        const double target = (double)cycle - train->phase;
        double t;
        if (fabs(slope) * (double)SUBFRAME < 1.0e-9 * from) {
            t = target / from;
        } else {
            /* Negative once the phase asked for lies past what a falling
             * pitch reaches: no further pulse in this interval. */
            const double radicand = (from * from) + (2.0 * slope * target);
            if (radicand < 0.0) {
                break;
            }
            t = (sqrt(radicand) - from) / slope;
        }
        if (t >= (double)SUBFRAME) {
            break;
        }
        if (train->count < TRAIN_PULSES) {
            train->time[train->count++] = start + (float)t;
        }
    }
    const double reached = train->phase + (0.5 * (double)SUBFRAME * (from + to));
    train->phase = reached - floor(reached);
}

/* Where a step between two voiced models counts as about an octave: the new
 * pitch period within 20% of twice the old one, or of half of it. */
#define OCTAVE_DOWN_LOW  1.6
#define OCTAVE_DOWN_HIGH 2.4
#define OCTAVE_UP_LOW    0.4
#define OCTAVE_UP_HIGH   0.6

/* About an octave apart the pulse train does not glide: DVSI's decoders first
 * take the old train an octave down (its phase and rate halved) when the new
 * period is 1.6 to 2.4 times the old, or up (its phase doubled, modulo a
 * cycle, and its rate doubled) at 0.4 to 0.6 times, and then cross linearly to
 * the new pitch as at any step. Measured on the tone-to-tone changes of DVSI's
 * rate-33 and D-STAR tone vectors (129 of 134 classified, each within 0.0012
 * cycle of its prediction; 5 cannot be told) and on the speech frames with
 * such steps. The inner edges are bracketed, not measured: 1.534 glides and
 * 1.618 halves; 0.611 glides and 0.584 doubles. Nothing brackets the outer
 * ones: 2.475 glides, but no step was observed between the largest halved
 * (1.72) and it, and the smallest observed to double is 0.528; 0.4 and 2.4 are
 * this library's choice, by symmetry. The train keeps its pulses' launch
 * times, not their models: a pulse launched before the step is rendered
 * through both end models of each later interval it reaches, as before; only
 * the phase and the rate the interval starts from change. */
static void
octave_step(struct pulse_train* train, double* from, double to) {
    const double ratio = *from / to; /* the new period over the old */
    if (ratio > OCTAVE_DOWN_LOW && ratio < OCTAVE_DOWN_HIGH) {
        train->phase *= 0.5;
        *from *= 0.5;
    } else if (ratio > OCTAVE_UP_LOW && ratio < OCTAVE_UP_HIGH) {
        const double doubled = 2.0 * train->phase;
        train->phase = doubled - floor(doubled);
        *from *= 2.0;
    }
}

/* One interval of the train, from model `a` to model `b`: through voiced
 * speech the pulses follow the pitch from one to the other (after an octave
 * step, octave_step()); into or out of it they keep the voiced model's; with
 * neither voiced the train starts afresh. */
static void
advance_interval(struct pulse_train* train, float start, const struct mbe_subframe_model* a,
                 const struct mbe_subframe_model* b) {
    double from = (double)a->w0 / (2.0 * M_PI);
    const double to = (double)b->w0 / (2.0 * M_PI);
    if (a->any_voiced && b->any_voiced) {
        octave_step(train, &from, to);
        interval_pulses(train, start, from, to);
    } else if (a->any_voiced) {
        interval_pulses(train, start, from, from);
    } else if (b->any_voiced) {
        interval_pulses(train, start, to, to);
    } else {
        /* A new train: the old one's pulses do not sound through the models after it. */
        train->phase = 0.0;
        train->count = 0;
    }
}

/* What the frame leaves the next: the pulse phase and the pulses of the last
 * frame's length, most recent first, in noiseOverlap; and, as the TIA-102.BABA
 * rendering leaves them, harmonic l's phase at the frame's end before the
 * regenerated phase (PSIl) and with it (PHIl). Element 0 of PSIl and PHIl is
 * not touched (the frame's mark lives there). */
static void
store_pulse_train(mbe_parms* cur_mp, const struct mbe_subframe_model* own, const struct pulse_train* train) {
    /* Under 20 pi for a pitch below Nyquist; a pitch no codec gives (a caller's
     * state) must not leave phases a float cannot hold. */
    const double advance = (2.0 * M_PI * train->phase) - ((double)PEAK * (double)own->w0);
    const float linear = (fabs(advance) < 1.0e6) ? (float)advance : 0.0f;
    /* Rounded to float before it is compared with 1: under x87 excess
     * precision the comparison could otherwise see the unrounded value, just
     * under 1, and store 1.0f (mbe_round_to_float(), mbe_math.h). */
    const float phase = mbe_round_to_float((float)train->phase);
    int kept = 0;

    for (int l = 1; l <= 56; l++) {
        cur_mp->PSIl[l] = (float)l * linear;
        cur_mp->PHIl[l] = cur_mp->PSIl[l] + own->phase[l];
    }
    memset(cur_mp->noiseOverlap, 0, sizeof(cur_mp->noiseOverlap));
    for (int i = train->count - 1; i >= 0 && kept < MBE_PULSE_HISTORY; i--) {
        if (train->time[i] >= 0.0f) {
            cur_mp->noiseOverlap[1 + kept++] = train->time[i] - (float)FRAME;
        }
    }
    cur_mp->noiseOverlap[0] = (float)kept;
    /* A phase that rounds up to a whole cycle is the start of the next one:
     * the stored phase lies in [0, 1). */
    cur_mp->noiseOverlap[PHASE_SLOT] = (phase < 1.0f) ? phase : 0.0f;
    mbe_subframe_state_mark(cur_mp);
}

/* ---- rendering the voiced part ---- */

/* A model's voiced bands, ready to render. */
struct rendered_model {
    int count;        /* voiced bands */
    float w0;         /* radians per sample */
    float period;     /* samples */
    float real[56];   /* amplitude times the cosine of the band's phase */
    float imag[56];   /* amplitude times its sine */
    float step_c[56]; /* cosine of the band's advance per sample */
    float step_s[56]; /* sine of it */
    short harmonic[56];
};

static void
prepare_render(const struct mbe_subframe_model* model, float voiced_gain, struct rendered_model* render) {
    const float hz_per_harmonic = model->w0 * (8000.0f / (2.0f * (float)M_PI));

    render->count = 0;
    render->w0 = model->w0;
    render->period = model->any_voiced ? (2.0f * (float)M_PI) / model->w0 : 0.0f;
    for (int l = 1; l <= model->L && model->any_voiced; l++) {
        const float w = model->w0 * (float)l;
        /* A harmonic at or above Nyquist would alias back into the band. */
        if (model->voiced[l] != BAND_VOICED || !(w < (float)M_PI)) {
            continue;
        }
        const float hz = hz_per_harmonic * (float)l;
        const float amplitude = 2.0f * model->magnitude[l] * voiced_gain * mbe_low_band_gain(hz);
        const int k = render->count++;
        float s, c;
        mbe_sincosf(model->phase[l] + mbe_low_band_phase(hz), &s, &c);
        render->real[k] = amplitude * c;
        render->imag[k] = amplitude * s;
        mbe_sincosf(w, &render->step_s[k], &render->step_c[k]);
        render->harmonic[k] = (short)l;
    }
}

/* The window of a response at t samples after its pulse: rising over PEAK
 * samples, whole until one pitch period after the pulse, falling over PEAK
 * more. A response lasts a period and PEAK samples, and where two of a steady
 * train overlap their windows add up to one, so the train adds up to its
 * harmonics. */
static float
response_window(float t, float period) {
    float s, c;
    if (t < (float)PEAK) {
        mbe_sincosf((float)M_PI * (t / (float)PEAK), &s, &c);
        return 0.5f * (1.0f - c);
    }
    if (t < period) {
        return 1.0f;
    }
    mbe_sincosf((float)M_PI * ((t - period) / (float)PEAK), &s, &c);
    return 0.5f * (1.0f + c);
}

/* Each band's oscillator at the response's first rendered sample: its phase
 * there is l w0 (t - PEAK) plus the band's own, t samples after the pulse. */
static void
start_oscillators(const struct rendered_model* render, float t, float* real, float* imag) {
    float turn_s, turn_c, s = 0.0f, c = 1.0f;
    int harmonic = 0;

    mbe_sincosf(render->w0 * (t - (float)PEAK), &turn_s, &turn_c);
    for (int k = 0; k < render->count; k++) {
        for (; harmonic < render->harmonic[k]; harmonic++) {
            const float next = (c * turn_c) - (s * turn_s);
            s = (s * turn_c) + (c * turn_s);
            c = next;
        }
        real[k] = (render->real[k] * c) - (render->imag[k] * s);
        imag[k] = (render->imag[k] * c) + (render->real[k] * s);
    }
}

/* Add to one interval the response a model plays for a pulse `pulse` samples
 * after the interval's start (negative: before it), under the crossfade that
 * brings the model in (`rising`) or takes it out. */
static void
render_pulse(float* out, const struct rendered_model* render, float pulse, int rising) {
    const float last = pulse + render->period + (float)PEAK;
    float real[56], imag[56];
    int first = (pulse > 0.0f) ? (int)ceilf(pulse) : 0;

    if (render->count == 0 || !(last > 0.0f) || first >= SUBFRAME) {
        return;
    }
    start_oscillators(render, (float)first - pulse, real, imag);
    for (int n = first; n < SUBFRAME && (float)n < last; n++) {
        const float fade = (float)n * (1.0f / (float)SUBFRAME);
        float sum = 0.0f;
        for (int k = 0; k < render->count; k++) {
            sum += real[k];
            const float next = (real[k] * render->step_c[k]) - (imag[k] * render->step_s[k]);
            imag[k] = (imag[k] * render->step_c[k]) + (real[k] * render->step_s[k]);
            real[k] = next;
        }
        out[n] += response_window((float)n - pulse, render->period) * (rising ? fade : (1.0f - fade)) * sum;
    }
}

static void
render_interval(float* out, float start, const struct pulse_train* train, const struct rendered_model* from,
                const struct rendered_model* to) {
    for (int i = 0; i < train->count; i++) {
        const float pulse = train->time[i] - start;
        render_pulse(out, from, pulse, 0);
        render_pulse(out, to, pulse, 1);
    }
}

/* ---- the unvoiced part ---- */

static void
unvoiced_subframe(const struct mbe_subframe_model* model, struct mbe_unvoiced_subframe* sub) {
    memset(sub, 0, sizeof(*sub));
    sub->w0 = model->w0;
    sub->L = model->L;
    for (int l = 1; l <= model->L; l++) {
        if (model->voiced[l] == BAND_UNVOICED) {
            sub->amplitude[l] = MBE_UNVOICED_SUBFRAME_SCALE * model->magnitude[l];
        }
    }
}

/* Whether the frame's voiced bands start in its own subframe: a frame so
 * marked (mbe_set_voicing_onset_delay()), in the rate-33 rendering whatever
 * the frame before (V/UV row 15), in the others only after a frame with no
 * voiced band (the mark's meaning there before the rate-33 rows' reading; in
 * P25's rendering a voiced model never reaches here with the mark, which
 * mark_rendered_model() replaces). */
static int
starts_late(const mbe_parms* cur_mp, const mbe_parms* prev_mp, enum mbe_voiced_synthesis rendering) {
    if (mbe_voicing_onset_delay(cur_mp) == 0) {
        return 0;
    }
    return rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 || !model_has_voiced_band(prev_mp);
}

/* A state with no overlap of this rendering's: what a frame that starts from
 * silence overlaps. */
static const mbe_parms no_overlap = {0};

/* A model that plays nothing and has no voiced band. */
static void
silence_model(struct mbe_subframe_model* model) {
    model->any_voiced = 0;
    for (int l = 1; l <= model->L; l++) {
        model->voiced[l] = BAND_NONE;
        model->magnitude[l] = 0.0f;
    }
}

/* Whether the synthesis takes the frame: both harmonic counts 1..56. A frame
 * it does not take renders nothing and leaves the state as it is. */
static int
frame_accepted(const mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    return mbe_harmonic_count_is_valid(cur_mp->L) && mbe_harmonic_count_is_valid(prev_mp->L);
}

/* A frame's three models (the last frame's own, mid-frame, this frame's own)
 * and the pulse train it starts from, for rendering and advancing alike.
 * Starting from silence, the last frame's own model plays nothing and counts
 * as having no voiced band, and the train it carried is dropped: a new one
 * starts at the frame's start, at phase 0. The mid-frame model is the
 * ordinary one, except for a frame marked MBE_FRAME_MARK_PREVIOUS_HELD. That
 * hold continues what was heard: a silence frame plays on the model the
 * listener heard last through its first 10 ms. A model that was not heard
 * (the one a decoder kept through a mute or a tone) is not held, so such a
 * frame from silence takes its own model at mid-frame: it rises from silence
 * to it over the first subframe and plays it over the second. DVSI's vectors
 * hold no silence frame right after a mute or a tone; this is the reading that
 * plays no speech from before them. Returns 0 when a harmonic count is not
 * 1..56. */
static int
prepare_frame(const mbe_parms* cur_mp, const mbe_parms* prev_mp, int from_silence, enum mbe_voiced_synthesis rendering,
              struct mbe_subframe_model models[3], struct pulse_train* train) {
    if (!frame_accepted(cur_mp, prev_mp)) {
        return 0;
    }
    mbe_subframe_own_model(prev_mp, &models[0]);
    mbe_subframe_mid_model(prev_mp, cur_mp, starts_late(cur_mp, prev_mp, rendering), rendering, &models[1]);
    mbe_subframe_own_model(cur_mp, &models[2]);
    if (from_silence) {
        silence_model(&models[0]);
        if (mbe_frame_mark(cur_mp) == MBE_FRAME_MARK_PREVIOUS_HELD) {
            models[1] = models[2];
        }
        train->phase = 0.0;
        train->count = 0;
    } else {
        load_pulse_train(prev_mp, &models[0], train);
    }
    return 1;
}

/* In P25's rendering, the mark of the model rendered, as it stands after
 * adaptive smoothing (which can voice harmonics of a frame decoded with none):
 * MBE_FRAME_MARK_ONSET_COLUMNS when it has a voiced band and is not a tone
 * model (MBE_FRAME_MARK_TONE), so that its onset is read in the 500 Hz columns
 * whatever mark it was decoded or repeated with (a caller's
 * MBE_FRAME_MARK_LATE_ONSET or MBE_FRAME_MARK_PREVIOUS_HELD there is
 * overwritten and has no effect); with no voiced band, not that mark (none
 * instead), any other mark kept. The prediction history the data path
 * committed before keeps the decoded model and its mark. Other renderings
 * leave the mark as it is. Called only on a frame the synthesis accepts
 * (frame_accepted()), so that a rejected call leaves the state untouched. */
static void
mark_rendered_model(mbe_parms* cur_mp, enum mbe_voiced_synthesis rendering) {
    if (rendering != MBE_VOICED_SYNTHESIS_DVSI_IMBE) {
        return;
    }
    const enum mbe_frame_mark mark = mbe_frame_mark(cur_mp);
    if (mark == MBE_FRAME_MARK_TONE) {
        return;
    }
    if (model_has_voiced_band(cur_mp)) {
        mbe_set_frame_mark(cur_mp, MBE_FRAME_MARK_ONSET_COLUMNS);
    } else if (mark == MBE_FRAME_MARK_ONSET_COLUMNS) {
        mbe_set_frame_mark(cur_mp, MBE_FRAME_MARK_NONE);
    }
}

void
mbe_synthesize_subframes(float* aout_buf, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                         float voiced_gain, int from_silence, enum mbe_voiced_synthesis rendering) {
    struct mbe_subframe_model models[3]; /* the last frame's own, mid-frame, this frame's own */
    struct rendered_model renders[3];
    struct mbe_unvoiced_subframe noise[2];
    struct pulse_train train;

    if (!frame_accepted(cur_mp, prev_mp)) {
        return;
    }
    mark_rendered_model(cur_mp, rendering);
    if (!prepare_frame(cur_mp, prev_mp, from_silence, rendering, models, &train)) {
        return;
    }
    for (int i = 0; i < 3; i++) {
        prepare_render(&models[i], voiced_gain, &renders[i]);
    }
    for (int k = 0; k < 2; k++) {
        const float start = (float)(k * SUBFRAME);
        advance_interval(&train, start, &models[k], &models[k + 1]);
        render_interval(&aout_buf[(size_t)k * SUBFRAME], start, &train, &renders[k], &renders[k + 1]);
    }
    unvoiced_subframe(&models[1], &noise[0]);
    unvoiced_subframe(&models[2], &noise[1]);
    /* Reads prev_mp's overlap and noise state before cur_mp's are written. */
    if (plan != NULL) {
        mbe_synthesizeUnvoicedSubframes(aout_buf, cur_mp, from_silence ? &no_overlap : prev_mp, plan, &noise[0],
                                        &noise[1]);
    } else {
        mbe_dvsi_noise_skip_frame(cur_mp);
    }
    store_pulse_train(cur_mp, &models[2], &train);
}

void
mbe_advance_subframes(mbe_parms* cur_mp, const mbe_parms* prev_mp, int from_silence,
                      enum mbe_voiced_synthesis rendering) {
    struct mbe_subframe_model models[3];
    struct pulse_train train;

    if (!frame_accepted(cur_mp, prev_mp)) {
        return;
    }
    mark_rendered_model(cur_mp, rendering);
    if (!prepare_frame(cur_mp, prev_mp, from_silence, rendering, models, &train)) {
        return;
    }
    advance_interval(&train, 0.0f, &models[0], &models[1]);
    advance_interval(&train, (float)SUBFRAME, &models[1], &models[2]);
    mbe_dvsi_noise_skip_frame(cur_mp);
    store_pulse_train(cur_mp, &models[2], &train);
}
