// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/** @file Speech synthesis in 10 ms subframes, as DVSI's decoders render it (mbe_subframe.h). */
#include "mbe_subframe.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

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

void
mbe_subframe_own_model(const mbe_parms* mp, struct mbe_subframe_model* model) {
    float magnitudes[57] = {0};

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
        magnitudes[l] = model->magnitude[l];
    }
    note_voiced(model);
    mbe_regenerate_phases(magnitudes, model->L, model->phase);
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

/* The mid-frame pitch: between the two frames' when both are pitched speech
 * or neither is, the voiced frame's when only one has a voiced band, unless
 * the unvoiced frame before a voiced one is marked as carrying its pitch
 * (MBE_FRAME_MARK_PITCH_CARRIES). */
static float
mid_pitch(const mbe_parms* prev_mp, const mbe_parms* cur_mp) {
    const int prev_voiced = model_has_voiced_band(prev_mp);
    const int cur_voiced = model_has_voiced_band(cur_mp);
    if (prev_voiced && !cur_voiced) {
        return prev_mp->w0;
    }
    if (cur_voiced && !prev_voiced && mbe_frame_mark(prev_mp) != MBE_FRAME_MARK_PITCH_CARRIES) {
        return cur_mp->w0;
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

/* The mid-frame model's voicing, magnitudes and phases, its pitch and harmonic
 * count being set: each band's state from the frame whose voicing holds at
 * mid-frame, and the two frames' magnitudes at the band's frequency, met in
 * the log domain. */
static void
fill_mid_bands(const mbe_parms* prev_mp, const mbe_parms* cur_mp, int late_onset, struct mbe_subframe_model* model) {
    const double unvoiced = log2((double)MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt((double)model->w0));
    const mbe_parms* voicing = voicing_persists(prev_mp, cur_mp) ? prev_mp : cur_mp;
    double prev_log2[57], cur_log2[57];
    float magnitudes[57] = {0};

    common_log2_magnitudes(prev_mp, prev_log2);
    common_log2_magnitudes(cur_mp, cur_log2);
    for (int l = 1; l <= model->L; l++) {
        const float w = model->w0 * (float)l;
        int state = band_state(voicing, harmonic_at(voicing, w));
        if (late_onset && state == BAND_VOICED) {
            state = BAND_UNVOICED;
        }
        const double mean = 0.5 * (log2_magnitude_at(prev_mp, prev_log2, w) + log2_magnitude_at(cur_mp, cur_log2, w));
        model->voiced[l] = (unsigned char)state;
        if (mean > LOG2_FLOOR) {
            /* Limited before it narrows: a caller's state can put it past what a float holds. */
            const double magnitude = exp2(mean + ((state == BAND_UNVOICED) ? unvoiced : 0.0));
            model->magnitude[l] = (magnitude < (double)MAGNITUDE_LIMIT) ? (float)magnitude : MAGNITUDE_LIMIT;
            magnitudes[l] = model->magnitude[l];
        }
    }
    note_voiced(model);
    mbe_regenerate_phases(magnitudes, model->L, model->phase);
}

void
mbe_subframe_mid_model(const mbe_parms* prev_mp, const mbe_parms* cur_mp, int late_onset,
                       struct mbe_subframe_model* model) {
    memset(model, 0, sizeof(*model));
    if (!mbe_harmonic_count_is_valid(prev_mp->L) || !mbe_harmonic_count_is_valid(cur_mp->L)) {
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
    const float w0 = mid_pitch(prev_mp, cur_mp);
    if (!positive_normal(&w0)) {
        return;
    }
    const float harmonics = floorf((0.9254f * (float)M_PI) / w0);
    model->w0 = w0;
    model->L = (harmonics < (float)shorter) ? shorter : ((harmonics > (float)longer) ? longer : (int)harmonics);
    fill_mid_bands(prev_mp, cur_mp, late_onset, model);
}

/* ---- the pulse train ---- */

struct pulse_train {
    double phase;             /* cycles since the last pulse, at the point reached */
    int count;                /* pulses that can still sound */
    float time[TRAIN_PULSES]; /* in samples from the frame's start (negative: before it) */
};

/* The train a stream carries into a frame: the pulse phase reached at the last
 * frame's end and the pulses it kept. It runs on from a frame that has a
 * voiced band, and from one marked MBE_FRAME_MARK_PITCH_ABSENT,
 * through whose first half the voicing before it may have lasted; after any
 * other frame it starts afresh, whatever a codec path left in the state.
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
    if (!prev->any_voiced && mbe_frame_mark(prev_mp) != MBE_FRAME_MARK_PITCH_ABSENT) {
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

/* One interval of the train, from model `a` to model `b`: through voiced
 * speech the pulses follow the pitch from one to the other; into or out of it
 * they keep the voiced model's; with neither voiced the train starts afresh. */
static void
advance_interval(struct pulse_train* train, float start, const struct mbe_subframe_model* a,
                 const struct mbe_subframe_model* b) {
    const double from = (double)a->w0 / (2.0 * M_PI);
    const double to = (double)b->w0 / (2.0 * M_PI);
    if (a->any_voiced && b->any_voiced) {
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
    const float phase = (float)train->phase;
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
    /* A phase that rounds up to a whole cycle is the start of the next one. */
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
 * marked (mbe_set_voicing_onset_delay()) after a frame with no voiced band. */
static int
starts_late(const mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    return mbe_voicing_onset_delay(cur_mp) != 0 && !model_has_voiced_band(prev_mp);
}

void
mbe_synthesize_subframes(float* aout_buf, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                         float voiced_gain) {
    struct mbe_subframe_model models[3]; /* the last frame's own, mid-frame, this frame's own */
    struct rendered_model renders[3];
    struct mbe_unvoiced_subframe noise[2];
    struct pulse_train train;

    if (!mbe_harmonic_count_is_valid(cur_mp->L) || !mbe_harmonic_count_is_valid(prev_mp->L)) {
        return;
    }
    mbe_subframe_own_model(prev_mp, &models[0]);
    mbe_subframe_mid_model(prev_mp, cur_mp, starts_late(cur_mp, prev_mp), &models[1]);
    mbe_subframe_own_model(cur_mp, &models[2]);
    load_pulse_train(prev_mp, &models[0], &train);
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
        mbe_synthesizeUnvoicedSubframes(aout_buf, cur_mp, prev_mp, plan, &noise[0], &noise[1]);
    } else {
        mbe_dvsi_noise_skip_frame(cur_mp);
    }
    store_pulse_train(cur_mp, &models[2], &train);
}

void
mbe_advance_subframes(mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    struct mbe_subframe_model models[3];
    struct pulse_train train;

    if (!mbe_harmonic_count_is_valid(cur_mp->L) || !mbe_harmonic_count_is_valid(prev_mp->L)) {
        return;
    }
    mbe_subframe_own_model(prev_mp, &models[0]);
    mbe_subframe_mid_model(prev_mp, cur_mp, starts_late(cur_mp, prev_mp), &models[1]);
    mbe_subframe_own_model(cur_mp, &models[2]);
    load_pulse_train(prev_mp, &models[0], &train);
    advance_interval(&train, 0.0f, &models[0], &models[1]);
    advance_interval(&train, (float)SUBFRAME, &models[1], &models[2]);
    mbe_dvsi_noise_skip_frame(cur_mp);
    store_pulse_train(cur_mp, &models[2], &train);
}
