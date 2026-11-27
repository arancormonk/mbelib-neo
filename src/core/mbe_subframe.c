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
#include "mbe_pulsed.h"
#include "mbe_validation.h"
#include "mbe_voiced.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FRAME             MBE_VOICED_FRAME
#define SUBFRAME          MBE_SUBFRAME
#define PEAK              MBE_PULSE_PEAK
/* Pulses one interval may start (a pitch period of 5 samples: no codec's
 * pitch comes near it) and pulses a frame's rendering can see. */
#define INTERVAL_PULSES   16
#define TRAIN_PULSES      (MBE_PULSE_HISTORY + 2 * INTERVAL_PULSES)
/* Band states of a model: a pulsed band (D-STAR's and P25's pulsed state,
 * pulse_bands()) plays neither the voiced pulse train nor noise. */
#define BAND_UNVOICED     0
#define BAND_VOICED       1
#define BAND_NONE         2
#define BAND_PULSED       3
/* log2 magnitude of a band with no amplitude, and the level under which an
 * interpolated magnitude counts as none. */
#define LOG2_SILENT       (-1000.0)
#define LOG2_FLOOR        (-200.0)
/* The largest magnitude a model takes. */
#define MAGNITUDE_LIMIT   1.0e30f
/* Where noiseOverlap holds the pulse phase, after the count and the pulses. */
#define PHASE_SLOT        (1 + MBE_PULSE_HISTORY)
/* Where it holds whether the frame's mid-frame model emitted its pulse
 * (pulsed_frame()): exactly PULSED_EMITTED for one that did, 0 otherwise. */
#define PULSED_SLOT       MBE_PULSED_SLOT
#define PULSED_EMITTED    MBE_PULSED_EMITTED
/* Where P25's rendering holds whether the frame's reserved code was read
 * (exactly PULSED_READ if it was) and the validity counter it left. */
#define PULSED_READ_SLOT  MBE_PULSED_READ_SLOT
#define PULSED_READ       MBE_PULSED_READ
#define PULSED_COUNT_SLOT MBE_PULSED_COUNT_SLOT
/* The bound on a log2 amplitude the phases are regenerated from: far past
 * any codec's (D-STAR's decoder reaches about 49), where a caller's state can
 * put anything. */
#define PHASE_LOG2_MAX    1000.0f

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
 * together lie between 0.998898 and 1.000587 (cos and sin are within [-1, 1]
 * for any finite t; t is finite for every positive normal w0, log2 P lying in
 * (log2(2 pi) - 128, log2(2 pi) + 126], log2(2 pi) = 2.6515). The pitch is
 * then positive and finite in double. A mean so near FLT_MAX that the product
 * passes it gives 0, no pitch, which the caller rejects (positive_normal()) as
 * it rejects every pitch that is not positive and normal (the conversion of a
 * double past the float range is not defined): the two largest floats, in
 * either order, reach about 1.0005 FLT_MAX. At the lower end no pair falls
 * under FLT_MIN: that would take both pitches within about 0.23% of FLT_MIN,
 * where the two factors together are about 1.0005. */
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

/* ---- the pulsed state of D-STAR and P25 ---- */

/*
 * DVSI's P25 and D-STAR decoders play a voice frame with no voiced band at one
 * of 56 reserved pitch codes (P25 b0 128..183, D-STAR b0 64..119, mbe_pulsed.h)
 * in a third voicing state (docs/dvsi-parity.md, "D-STAR and P25: the pulsed
 * state"): each code names, for the frame's own subframe and for its mid-frame
 * subframe, which of eight 500 Hz channels are pulsed, and a band there plays
 * no noise but the response to one pulse. The patterns are Table 2 of US
 * 7,970,606 (its P25 pitch values 128 to 183), indexed b0 - 128 for P25 and
 * b0 - 64 for D-STAR: [code][0] the own subframe's, [code][1] the mid-frame
 * subframe's, channel k (k = 0..7, 500 k to 500 (k + 1) Hz, the last open
 * above) lowest first, 'P' pulsed. That D-STAR's codes carry the same
 * patterns, and that DVSI's decoders play them so, are measurements on DVSI's
 * decodes (a DVMEGA DVstick 30's full-band curves corroborate them at one
 * D-STAR code, no more); the codes whose own subframe has a
 * pulsed channel (0..7, 24..39, 48..55) are those of D-STAR's
 * MBE_FRAME_MARK_PITCH_CARRIES.
 */
static const char pulsed_patterns[MBE_PULSED_CODES][2][9] = {
    {"UUUUUUUP", "UUUUUUUU"}, {"UUUUUUPU", "UUUUUUUU"}, {"UUUUUUPP", "UUUUUUUU"}, {"UUUUUPUU", "UUUUUUUU"}, /*  0 */
    {"UUUUUPPP", "UUUUUUUU"}, {"UUUUPPPP", "UUUUUUUU"}, {"UUUPUUUU", "UUUUUUUU"}, {"UUUPPUUU", "UUUUUUUU"}, /*  4 */
    {"UUUUUUUU", "PUUUUUUU"}, {"UUUUUUUU", "UUUPPUUU"}, {"UUUUUUUU", "UUUPUUUU"}, {"UUUUUUUU", "UUPPPPPU"}, /*  8 */
    {"UUUUUUUU", "UUPUUUUU"}, {"UUUUUUUU", "UUPPUUUU"}, {"UUUUUUUU", "UUPPPUUU"}, {"UUUUUUUU", "UUUPPPUU"}, /* 12 */
    {"UUUUUUUU", "UUUUUUUP"}, {"UUUUUUUU", "UUUUUPPP"}, {"UUUUUUUU", "UUUUUUPU"}, {"UUUUUUUU", "UUUUUPUU"}, /* 16 */
    {"UUUUUUUU", "UUUUUUPP"}, {"UUUUUUUU", "UUUUPPPP"}, {"UUUUUUUU", "UUUUPUUU"}, {"UUUUUUUU", "UUUUPPUU"}, /* 20 */
    {"PUUUUUUU", "UUUUUUUU"}, {"UUPPPUUU", "UUUUUUUU"}, {"PPPUUUUU", "UUUUUUUU"}, {"UUUPPPUU", "UUUUUUUU"}, /* 24 */
    {"PPUUUUUU", "UUUUUUUU"}, {"PPPPUUUU", "UUUUUUUU"}, {"UUUPPPPP", "UUUUUUUU"}, {"UUPPUUUU", "UUUUUUUU"}, /* 28 */
    {"PPPPPPPP", "UUUUUUUU"}, {"UUUUPUUU", "UUUUUUUU"}, {"UUPPPPPP", "UUUUUUUU"}, {"PPPPPUUU", "UUUUUUUU"}, /* 32 */
    {"UUPUUUUU", "UUUUUUUU"}, {"PPPPPPUU", "UUUUUUUU"}, {"UPPPPPPP", "UUUUUUUU"}, {"PPPPPPPU", "UUUUUUUU"}, /* 36 */
    {"UUUUUUUU", "PPPPPPPP"}, {"UUUUUUUU", "PPPPPUUU"}, {"UUUUUUUU", "PPPUUUUU"}, {"UUUUUUUU", "UPPPPPPP"}, /* 40 */
    {"UUUUUUUU", "PPUUUUUU"}, {"UUUUUUUU", "UUUPPPPP"}, {"UUUUUUUU", "UUPPPPPP"}, {"UUUUUUUU", "PPPPPPUU"}, /* 44 */
    {"PPPPPPPP", "PPPPPPPP"}, {"PPPPPPPP", "UUPPPPPP"}, {"PPPUUUUU", "PPPUUUUU"}, {"UUUPPPPP", "UUUPPPPP"}, /* 48 */
    {"UUPPPPPP", "UUPPPPPP"}, {"PPPPPPPP", "PPPUUUUU"}, {"PPUUUUUU", "PPUUUUUU"}, {"UUPUUUUU", "UUPUUUUU"}, /* 52 */
};

/* Which of a code's two patterns applies: the own subframe's, to a frame's own
 * model (at its end, and as the next frame's first model), or the mid-frame
 * subframe's, to its mid-frame model. */
#define PULSED_OWN      0
#define PULSED_MID      1
/* The second channel a band is read in lies this far above it. */
#define PULSED_REACH_HZ 250.0

/* The IMBE renderings: P25's, and ProVoice's, which is P25's without the
 * pulsed state. */
static int
imbe_rendering(enum mbe_voiced_synthesis rendering) {
    return rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE || rendering == MBE_VOICED_SYNTHESIS_DVSI_PROVOICE;
}

/* The reserved code (0..55) a model has in `rendering`, or -1: one with no
 * voiced band, as the synthesis renders it (after adaptive smoothing, which
 * can voice a band of a frame decoded with none: it then has no code), whose
 * pitch is the one the codec's decoder stores for a reserved code, by its bit
 * pattern (mbe_imbe4400_pulsed_code() for P25, mbe_ambe2400_pulsed_code() for
 * D-STAR, each with its decoder's own pitch function). Nothing but the model's
 * voicing and pitch selects a pattern: no frame mark is read, so a caller's
 * model at a reserved pitch with no voiced band is a frame with that code (it
 * holds nothing that could be forged). Only D-STAR's rendering
 * (MBE_VOICED_SYNTHESIS_DVSI) and P25's (MBE_VOICED_SYNTHESIS_DVSI_IMBE) have
 * reserved codes; every other value (rate 33's, ProVoice's, the TIA
 * rendering's, any other) has none. Whether P25's rendering reads a code a
 * model has is the validity counter's to say (read_code()). The harmonic
 * count is 1..56 here. */
static int
pulsed_code(const mbe_parms* mp, enum mbe_voiced_synthesis rendering) {
    if ((rendering != MBE_VOICED_SYNTHESIS_DVSI && rendering != MBE_VOICED_SYNTHESIS_DVSI_IMBE)
        || !positive_normal(&mp->w0) || model_has_voiced_band(mp)) {
        return -1;
    }
    return (rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE) ? mbe_imbe4400_pulsed_code(mp->w0)
                                                         : mbe_ambe2400_pulsed_code(mp->w0);
}

/* Whether a state with this rendering's subframe mark holds exactly `value` in
 * noiseOverlap[slot], by its bit pattern: anything else a caller's state holds
 * there (not a number, an infinity, a subnormal, a negative zero, a near
 * value), and any value in a state the other rendering wrote, is not it. */
static int
stored_exactly(const mbe_parms* mp, int slot, float value) {
    uint32_t bits, want;
    if (!mbe_subframe_state_is_marked(mp)) {
        return 0;
    }
    memcpy(&bits, &mp->noiseOverlap[slot], sizeof(bits));
    memcpy(&want, &value, sizeof(want));
    return bits == want;
}

/*
 * P25's validity counter as the frame before left it (PULSED_COUNT_SLOT):
 * one of the whole numbers MBE_PULSED_COUNT_MIN..MBE_PULSED_COUNT_MAX stored
 * as a float, by its bit pattern; anything else, and a state without this
 * rendering's mark (a fresh stream), is 0.
 */
static int
stored_count(const mbe_parms* prev_mp) {
    for (int count = MBE_PULSED_COUNT_MIN; count <= MBE_PULSED_COUNT_MAX; count++) {
        if (stored_exactly(prev_mp, PULSED_COUNT_SLOT, (float)count)) {
            return count;
        }
    }
    return 0;
}

/*
 * The reserved code of this frame's model that the rendering reads, or -1.
 * D-STAR's reads every code a model has (pulsed_code()). P25's reads it only
 * while the validity counter the frame before left stands at 0 or more
 * (stored_count()), as US 7,970,606 describes its decoder (FIG. 8, steps 815
 * to 850: the table is read only while ValidCount is 0 or more), a guard
 * against a stream from an encoder that does not use the codes. The counter is
 * the document's, not a measurement: on DVSI's own P25 encodings it never
 * falls below 0 (docs/dvsi-parity.md, "D-STAR and P25: the pulsed state").
 */
static int
read_code(const mbe_parms* cur_mp, const mbe_parms* prev_mp, enum mbe_voiced_synthesis rendering) {
    const int code = pulsed_code(cur_mp, rendering);
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE && stored_count(prev_mp) < 0) {
        return -1;
    }
    return code;
}

/* The reserved code of the last frame's model that its frame read, or -1: in
 * P25's rendering only where that frame's state says its code was read
 * (PULSED_READ_SLOT), in D-STAR's every code it has. */
static int
last_read_code(const mbe_parms* prev_mp, enum mbe_voiced_synthesis rendering) {
    const int code = pulsed_code(prev_mp, rendering);
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE && !stored_exactly(prev_mp, PULSED_READ_SLOT, PULSED_READ)) {
        return -1;
    }
    return code;
}

/*
 * The validity counter this frame leaves in P25's rendering (0 in every other,
 * where none is kept), from the one the frame before left. The document
 * counts received frames: the counter moves for a model at the pitch some P25
 * pitch code decodes to (b0 0..207, mbe_imbe4400_pitch_code()) and no other.
 * A model with a voiced band (a tone model included) leaves it; one with no
 * voiced band at the pitch of a code of the table, one of the 56 or the
 * default code b0 25 (whose table entry is unvoiced in both subframes:
 * nothing is read for it), raises it by one, to MBE_PULSED_COUNT_MAX at most;
 * one with no voiced band at another code's pitch lowers it by one, to
 * MBE_PULSED_COUNT_MIN at least; one at a pitch no code decodes to (the repeat
 * limit's reset model and the initial state, at eq 46's nominal pitch of b0
 * 134, or a model a caller built) leaves it. The bounds and the steps are the
 * document's; that the counter starts at 0, counts the model as rendered
 * (after adaptive smoothing; a muted or repeated frame by the model it renders
 * or advances: a repeat within the limit replays a decoded model, and so
 * counts by its code again) and lets a model at no code's pitch leave it are
 * this library's choices.
 */
static int
next_count(const mbe_parms* cur_mp, const mbe_parms* prev_mp, enum mbe_voiced_synthesis rendering) {
    if (rendering != MBE_VOICED_SYNTHESIS_DVSI_IMBE) {
        return 0;
    }
    const int count = stored_count(prev_mp);
    if (model_has_voiced_band(cur_mp) || !positive_normal(&cur_mp->w0)) {
        return count;
    }
    const int b0 = mbe_imbe4400_pitch_code(cur_mp->w0);
    if (b0 < 0) {
        return count;
    }
    if (b0 == MBE_PULSED_P25_DEFAULT || (b0 >= MBE_PULSED_P25_FIRST && b0 < MBE_PULSED_P25_FIRST + MBE_PULSED_CODES)) {
        return (count < MBE_PULSED_COUNT_MAX) ? count + 1 : MBE_PULSED_COUNT_MAX;
    }
    return (count > MBE_PULSED_COUNT_MIN) ? count - 1 : MBE_PULSED_COUNT_MIN;
}

/* Whether a code's pattern (PULSED_OWN or PULSED_MID) has a pulsed channel;
 * none for no code (-1). */
static int
pattern_pulses(int code, int pattern) {
    if (code < 0) {
        return 0;
    }
    for (int k = 0; k < 8; k++) {
        if (pulsed_patterns[code][pattern][k] == 'P') {
            return 1;
        }
    }
    return 0;
}

/* The channel (0..7) of a frequency in Hz (finite, not negative): its 500 Hz
 * column, at most 7. Decided on the bit pattern, positive doubles ordering as
 * their patterns do, so that the conversion only ever sees a value under 7. */
static int
pulsed_channel(double hz) {
    const double position = hz / 500.0;
    const double top = 7.0;
    uint64_t bits, top_bits;
    memcpy(&bits, &position, sizeof(bits));
    memcpy(&top_bits, &top, sizeof(top_bits));
    return (bits >= top_bits) ? 7 : (int)position;
}

/*
 * A model's pulsed bands under a code's pattern: in a model with no voiced
 * band, each unvoiced harmonic l, at f = l w0 4000 / pi Hz, whose channel, or
 * the channel of f + 250 Hz, is pulsed. A pulsed band takes its magnitude as a
 * voiced one (the unvoiced factor MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt(w0) of
 * that model taken out: the pulse is the voiced pulse response) and plays no
 * noise (unvoiced_subframe()). A band with no state keeps none. The channels
 * of 500 Hz from 0 Hz are the table's. The 250 Hz reach was first fitted on
 * the study's numerical model of P25 frames deep inside runs with no voiced
 * band (against 0, 125 and 375 Hz) and measured again, on DVSI's decodes, by
 * scratch builds of the library after the gate sets were declared (against 0,
 * 125, 375 and 500 Hz). The frequency is formed in double from the float
 * pitch: a harmonic just under a channel's edge stays under it (in float it
 * could round onto it). A model with a voiced band (a
 * mid-frame model a caller's MBE_FRAME_MARK_PITCH_ABSENT held, say) takes no
 * pulsed band, so that the two excitations never meet in one model.
 */
static void
pulse_bands(struct mbe_subframe_model* model, int code, int pattern) {
    if (!pattern_pulses(code, pattern) || model->any_voiced || !positive_normal(&model->w0)) {
        return;
    }
    /* As a float, also under x87 excess precision. */
    const float factor = mbe_round_to_float(MBE_UNVOICED_MAGNITUDE_FACTOR / sqrtf(model->w0));
    const char* channels = pulsed_patterns[code][pattern];
    for (int l = 1; l <= model->L; l++) {
        const double hz = (double)l * (double)model->w0 * (4000.0 / M_PI);
        if (model->voiced[l] != BAND_UNVOICED
            || (channels[pulsed_channel(hz)] != 'P' && channels[pulsed_channel(hz + PULSED_REACH_HZ)] != 'P')) {
            continue;
        }
        /* Limited before it narrows: a caller's mid-frame pitch can make the factor tiny. */
        const double magnitude = (double)model->magnitude[l] / (double)factor;
        model->voiced[l] = BAND_PULSED;
        model->magnitude[l] = (magnitude < (double)MAGNITUDE_LIMIT) ? (float)magnitude : MAGNITUDE_LIMIT;
        model->any_pulsed = 1;
    }
}

/*
 * The mid-frame pitch between two frames with no voiced band of which at least
 * one has a code that is read (read_code(), last_read_code()), as DVSI's
 * decoders play it: in P25 this
 * frame's; in D-STAR the frame before's, unless this frame's own-subframe
 * pattern has a pulsed channel, and then this frame's, or the geometric mean of
 * the two if the frame before's own-subframe pattern has one too. Measured on
 * DVSI's decodes of the development and validation vectors against the
 * geometric mean everywhere and against each codec's rule in the other. Returns
 * 0, and leaves *w0, where neither frame has a code that is read: between two
 * frames with no voiced band and no such code the geometric mean stays (DVSI's
 * default-code frames all share one pitch, so nothing there was measured).
 */
static int
pulsed_mid_pitch(const mbe_parms* prev_mp, const mbe_parms* cur_mp, enum mbe_voiced_synthesis rendering, float* w0) {
    const int cur = read_code(cur_mp, prev_mp, rendering);
    const int prev = last_read_code(prev_mp, rendering);
    if (cur < 0 && prev < 0) {
        return 0;
    }
    if (imbe_rendering(rendering)) {
        *w0 = cur_mp->w0;
    } else if (!pattern_pulses(cur, PULSED_OWN)) {
        *w0 = prev_mp->w0;
    } else {
        *w0 = pattern_pulses(prev, PULSED_OWN) ? sqrtf(prev_mp->w0 * cur_mp->w0) : cur_mp->w0;
    }
    return 1;
}

/* The mid-frame pitch: between the two frames' when both are pitched speech
 * or neither is, the voiced frame's when only one has a voiced band, unless
 * the unvoiced frame before a voiced one is marked as carrying its pitch
 * (MBE_FRAME_MARK_PITCH_CARRIES); between two frames with no voiced band of
 * which one has a reserved code that is read, in D-STAR's and P25's
 * renderings, as pulsed_mid_pitch() says.
 * Between two voiced frames of different pitch, the rate-33 rendering plays
 * the geometric mean times MBE_AMBE2450_MID_PITCH_FACTOR times the half-octave
 * wave (ambe2450_mid_pitch(); docs/dvsi-parity.md, "Rate 33: mid-frame
 * pitch"), a fit to DVSI's rate-33 decodes (mbe_subframe.h): over a frame
 * whose pitch code steps by one to five codes from the frame before's, the
 * lag between DVSI's decode and ours read at its two boundaries changes by a
 * median of 0.040 to 0.066 samples, against 0.018 over one whose code is the
 * same (the level of everything that is not the mid-frame pitch), and a wave
 * in the two frames' periods with a period of half an octave fits that
 * change. How the decoder comes to it is not known. The same measurement's
 * grouped medians show no such jump at a step in D-STAR or P25. Both pitches
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
    float pulsed;
    if (!prev_voiced && !cur_voiced && pulsed_mid_pitch(prev_mp, cur_mp, rendering, &pulsed)) {
        return pulsed;
    }
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 && prev_voiced && cur_voiced
        && !same_pitch(&prev_mp->w0, &cur_mp->w0)) {
        return ambe2450_mid_pitch(prev_mp->w0, cur_mp->w0);
    }
    return sqrtf(prev_mp->w0 * cur_mp->w0);
}

/* ---- the mid-frame voicing of D-STAR and P25: 500 Hz columns ---- */

/*
 * The mid-frame model of D-STAR's and P25's renderings is voiced by one rule
 * (docs/dvsi-parity.md, "D-STAR and P25: mid-frame voicing by columns"). It
 * was fitted to per-harmonic labels of DVSI's decodes of its development and
 * validation vectors (each label says which of two test decodes is nearer
 * DVSI's in a harmonic's band), and the DVMEGA DVstick 30's single-band step
 * curves of #157 come closer with it; it is not taken from a document, and what
 * DVSI's decoders compute is not known. A band in the same state in both
 * frames (each read at its harmonic nearest the band, nearest_harmonic())
 * keeps it. One whose state differs is voiced at mid-frame if and only if M(k)
 * or M(k') holds, k the 500 Hz column of the band's frequency and k' the
 * column of the frequency 250 Hz higher:
 *
 *     M(k) = (Sp(k) and Sc(k)) or (Sp(k) != Sc(k) and Sp(k + 1) != Sc(k + 1))
 *
 * with Sp and Sc the two frames' voicing sampled at the columns (the samples,
 * not the bands' states), the sample past the last column unvoiced in both. A
 * frame with no voiced band is unvoiced in every band and sample; a band that
 * plays nothing in either frame (BAND_NONE, a caller's state) plays nothing
 * (column_band_state()). At an onset M(k) reduces to Sc(k) and Sc(k + 1); at a
 * frame marked MBE_FRAME_MARK_PITCH_ABSENT with no voiced band, to Sp(k) and
 * Sp(k + 1).
 */

/* Columns 0..7, and the sample past the last, always unvoiced. */
#define COLUMNS          8
/* P25's columns: column 0 from 0 Hz to MBE_VOICING_REACH_HZ, each next one
 * P25_COLUMN_HZ wide, column 7 above 3505 Hz. */
#define P25_COLUMN_HZ    500.0
/* The second column the rule reads lies this far above the band. */
#define COLUMN_REACH_HZ  250.0
/* P25's last band, harmonics 34 to L of a frame with more than 36 (12 bands). */
#define P25_FULL_BANDS_L 36
#define P25_LAST_BAND    34

/* A frame's voicing sampled at the columns: [COLUMNS] is past the last. */
struct column_samples {
    int voiced[COLUMNS + 1];
};

/* Whether a frame is voiced at its harmonic at or next above `position` (in
 * harmonics of its pitch, positive and finite), harmonics past `last` (at most
 * L) reading unvoiced. Decided on the bit pattern, positive doubles ordering
 * as their patterns do, so that the conversion only ever sees a value from 0
 * to `last`. */
static int
voiced_from(const mbe_parms* mp, double position, int last) {
    const double end = (double)last;
    uint64_t bits, end_bits;
    memcpy(&bits, &position, sizeof(bits));
    memcpy(&end_bits, &end, sizeof(end_bits));
    if (bits > end_bits) {
        return 0;
    }
    const int l = (int)ceil(position);
    return mp->Vl[(l < 1) ? 1 : l] == 1;
}

/* P25's samples (the reading of PR #198's onset columns): S(0) at harmonic 1,
 * S(k) for k = 1..7 at the harmonic at or next above MBE_VOICING_REACH_HZ +
 * P25_COLUMN_HZ (k - 1), unvoiced past L. The last band of a frame with more
 * than P25_FULL_BANDS_L harmonics, voiced while harmonic P25_LAST_BAND - 1 is
 * not, reads as unvoiced (a fit to the labels of DVSI's decodes, which mostly
 * do not sound such a band at mid-frame where a column it covers would:
 * docs/dvsi-parity.md); the frame's own model still voices it. The edges'
 * positions among the frame's harmonics are quotients formed in
 * double from its stored pitch, positive and normal: for every P25 pitch code
 * (and ProVoice's nominal pitch) each lies at least 1.5e-4 of a harmonic from a
 * whole number (test_p25_onset_edges in test_subframe_synthesis.c), far past
 * any difference the arithmetic of a build (x87 excess precision, fast-math's
 * reordering) can make. */
static void
p25_samples(const mbe_parms* mp, struct column_samples* samples) {
    const double harmonics_per_hz = M_PI / (4000.0 * (double)mp->w0);
    const int isolated_last =
        (mp->L > P25_FULL_BANDS_L) && (mp->Vl[P25_LAST_BAND] == 1) && (mp->Vl[P25_LAST_BAND - 1] != 1);
    const int last = isolated_last ? P25_LAST_BAND - 1 : mp->L;
    samples->voiced[0] = mp->Vl[1] == 1;
    for (int k = 1; k < COLUMNS; k++) {
        const double hz = (double)MBE_VOICING_REACH_HZ + (P25_COLUMN_HZ * (double)(k - 1));
        samples->voiced[k] = voiced_from(mp, hz * harmonics_per_hz, last);
    }
    samples->voiced[COLUMNS] = 0;
}

/* P25's column (0..7) of a frequency in Hz (finite, not negative). */
static int
p25_column(double hz) {
    int k = 0;
    while (k < COLUMNS - 1 && hz >= (double)MBE_VOICING_REACH_HZ + (P25_COLUMN_HZ * (double)k)) {
        k++;
    }
    return k;
}

/* A D-STAR position among the 500 Hz columns as the decoder forms it to voice
 * harmonic l (mbe_ambe2400_harmonic_voiced(): l 16 f0, rounded to float, also
 * under x87 excess precision). */
static float
dstar_position(int l, float f0) {
    return mbe_round_to_float((float)l * 16.0f * f0);
}

/* The D-STAR column holding a position (positive, or zero, or infinite), at
 * most 7, as the decoder takes it: its whole part. Decided on the bit pattern,
 * so that the conversion only ever sees a value under 7: a caller's pitch can
 * put a position past any int. */
static int
dstar_column(float position) {
    const float top = 7.0f;
    uint32_t bits, top_bits;
    memcpy(&bits, &position, sizeof(bits));
    memcpy(&top_bits, &top, sizeof(top_bits));
    return (bits >= top_bits) ? 7 : (int)position;
}

/* A frame's fundamental in cycles per sample, rounded to float (x87). */
static float
dstar_f0(float w0) {
    return mbe_round_to_float(w0 / (2.0f * (float)M_PI));
}

/* D-STAR's samples: both columns of 1 kHz band j (b1's bit for it) take the
 * frame's voicing at the band's lowest harmonic, the first whose column (by
 * the decoder's arithmetic) lies in the band; a band that holds no harmonic (a
 * caller's state) reads as unvoiced. The pitch is the stored one, positive and
 * normal; at the ordinary D-STAR pitch codes, 0 to 123, the columns it gives
 * are the decoder's, each position at least 6.5e-5 of a column from a whole
 * number and 1.8e-5 from a half in every build measured
 * (test_dstar_column_edges in test_subframe_synthesis.c). Code 124, the
 * silence model (L 14, 250 Hz), puts its harmonics exactly on positions 0.5
 * and 1.0, a margin of zero; it decodes with no voiced band, and an exact edge
 * lies in the upper column (test_dstar_exact_column_edge). */
static void
dstar_samples(const mbe_parms* mp, struct column_samples* samples) {
    const float f0 = dstar_f0(mp->w0);
    int band = 0;
    memset(samples, 0, sizeof(*samples));
    for (int l = 1; l <= mp->L && band < COLUMNS / 2; l++) {
        const int own = dstar_column(dstar_position(l, f0)) / 2;
        if (own >= band) {
            const int column = own + own;
            samples->voiced[column] = mp->Vl[l] == 1;
            samples->voiced[column + 1] = mp->Vl[l] == 1;
            band = own + 1;
        }
    }
}

/* The two frames' samples and how the rendering finds a band's columns. */
struct column_rule {
    struct column_samples prev, cur;
    int dstar; /* D-STAR's columns, else P25's */
};

/* Whether the mid-frame model's bands are voiced by the columns: in D-STAR's
 * and P25's renderings, for a frame with a voiced band, or with none when it is
 * marked MBE_FRAME_MARK_PITCH_ABSENT. No mark of the previous model is read
 * (that model is caller state: a P25 tone model before this frame is read in
 * the columns as any model is). Sets the samples. */
static int
column_rule(const mbe_parms* prev_mp, const mbe_parms* cur_mp, enum mbe_voiced_synthesis rendering,
            struct column_rule* rule) {
    if (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) {
        return 0;
    }
    if (!model_has_voiced_band(cur_mp) && mbe_frame_mark(cur_mp) != MBE_FRAME_MARK_PITCH_ABSENT) {
        return 0;
    }
    rule->dstar = !imbe_rendering(rendering);
    if (rule->dstar) {
        dstar_samples(prev_mp, &rule->prev);
        dstar_samples(cur_mp, &rule->cur);
    } else {
        p25_samples(prev_mp, &rule->prev);
        p25_samples(cur_mp, &rule->cur);
    }
    return 1;
}

/* M(k): column k (0..7) voiced at mid-frame. */
static int
column_voiced(const struct column_rule* rule, int k) {
    const int* p = rule->prev.voiced;
    const int* c = rule->cur.voiced;
    return (p[k] && c[k]) || ((p[k] != c[k]) && (p[k + 1] != c[k + 1]));
}

/* Whether harmonic l of the mid-frame model, of pitch w0 (positive and
 * normal), lies in a column voiced at mid-frame, or 250 Hz below one. Where
 * the mid-frame pitch is a frame's own (an onset after a frame whose pitch does
 * not carry, a frame lacking a pitch, two frames of one pitch code) the
 * harmonics lie as far from the columns' edges as the frames' do
 * (test_p25_onset_edges, test_dstar_column_edges, for the ordinary codes);
 * between two frames of different pitch, an onset after a pitch-carrying code
 * among them, their geometric mean can put a harmonic anywhere. */
static int
columns_voice(const struct column_rule* rule, float w0, int l) {
    if (rule->dstar) {
        const float position = dstar_position(l, dstar_f0(w0));
        const int above = dstar_column(mbe_round_to_float(position + 0.5f));
        return column_voiced(rule, dstar_column(position)) || column_voiced(rule, above);
    }
    const double hz = ((double)w0 * (double)l * 4000.0) / M_PI;
    return column_voiced(rule, p25_column(hz)) || column_voiced(rule, p25_column(hz + COLUMN_REACH_HZ));
}

/* A frame's harmonic nearest harmonic l of the mid-frame model of pitch w0
 * (both pitches positive and normal), 1..L, for the columns' rule: the
 * position l w0 / mp->w0 formed in double from the stored pitches, a half
 * added and its whole part taken, so that a position exactly halfway between
 * two harmonics (a tie) takes the upper one. In double because the float
 * quotient harmonic_at() forms rounds a position just under a half (33.4999986
 * of a P25 frame's harmonics at b0 134, its mid-frame harmonic 36 after a
 * frame at b0 161) to the half itself, and so to the wrong neighbour: a band
 * voiced in both frames was taken for one that changes. The index is made on
 * the bit pattern, positive doubles ordering as their patterns do. */
static int
nearest_harmonic(const mbe_parms* mp, float w0, int l) {
    const double position = (((double)w0 * (double)l) / (double)mp->w0) + 0.5;
    const double first = 1.0;
    const double last = (double)mp->L;
    uint64_t bits, first_bits, last_bits;
    memcpy(&bits, &position, sizeof(bits));
    memcpy(&first_bits, &first, sizeof(first_bits));
    memcpy(&last_bits, &last, sizeof(last_bits));
    if (bits < first_bits) {
        return 1;
    }
    return (bits >= last_bits) ? mp->L : (int)position;
}

/* The state of the mid-frame model's band l by the columns: each frame's state
 * at its harmonic nearest the band (nearest_harmonic()); a band that plays
 * nothing in either frame (BAND_NONE: a caller's voicing that is neither 0 nor
 * 1) plays nothing; one in the same state in both keeps it; else the columns
 * decide. */
static int
column_band_state(const mbe_parms* prev_mp, const mbe_parms* cur_mp, const struct column_rule* rule,
                  const struct mbe_subframe_model* model, int l) {
    const int before = band_state(prev_mp, nearest_harmonic(prev_mp, model->w0, l));
    const int now = band_state(cur_mp, nearest_harmonic(cur_mp, model->w0, l));
    if (before == BAND_NONE || now == BAND_NONE) {
        return BAND_NONE;
    }
    if (before == now) {
        return now;
    }
    return columns_voice(rule, model->w0, l) ? BAND_VOICED : BAND_UNVOICED;
}

/* Whether the rate-33 rendering's mid-frame model keeps the previous frame's
 * voicing instead of taking this frame's: a frame with no voiced band marked
 * MBE_FRAME_MARK_PITCH_ABSENT, after a frame voiced at its harmonic at or next
 * above MBE_VOICING_REACH_HZ. No rate-33 decoder writes that mark (a caller's
 * state can hold it); this is the hold every rendering applied before the
 * columns' rule replaced it in D-STAR's and P25's, kept so that rate 33's
 * rendering stays exactly as it was. */
static int
rate33_voicing_holds(const mbe_parms* prev_mp, const mbe_parms* cur_mp) {
    if (mbe_frame_mark(cur_mp) != MBE_FRAME_MARK_PITCH_ABSENT || model_has_voiced_band(cur_mp)) {
        return 0;
    }
    const float position = ceilf((MBE_VOICING_REACH_HZ * ((2.0f * (float)M_PI) / 8000.0f)) / prev_mp->w0);
    int within;
    return prev_mp->Vl[harmonic_below(position, prev_mp->L, &within)] == 1;
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
 * count being set: each band's state by the columns (`columns`, column_rule())
 * or, where they do not apply, that of the frame whose voicing stands at
 * mid-frame (`voicing`) at the band's frequency, and the
 * two frames' magnitudes at the band's frequency, met in the log domain. A late
 * onset (`late_onset`) voices no band; a frame with a mid-frame row (`row`,
 * mid_row()) takes each band's voicing from the row (a band with no state
 * keeps none). Its phases come from the mean of the two frames' decoded
 * log2 amplitudes (phase_log2_amplitudes()), each read at the band's
 * frequency as the magnitudes are and used as the decoder stored it: where
 * common_log2_magnitudes() takes the unvoiced factor out of an unvoiced band's
 * Ml, nothing is taken out of or put into its log2Ml, which carries no such
 * factor (an AMBE decoder builds an ordinary voice frame's unvoiced band as
 * Ml = MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt(w0) times 2^log2Ml; Ml then also
 * carries what log2Ml does not, such as rate 33's further silence scale and
 * spectral enhancement). */
static void
fill_mid_bands(const mbe_parms* prev_mp, const mbe_parms* cur_mp, const mbe_parms* voicing, int late_onset,
               const struct column_rule* columns, const int* row, struct mbe_subframe_model* model) {
    const double unvoiced = log2((double)MBE_UNVOICED_MAGNITUDE_FACTOR / sqrt((double)model->w0));
    double prev_log2[57], cur_log2[57], prev_phase_log2[57], cur_phase_log2[57];
    float envelope[57] = {0};

    common_log2_magnitudes(prev_mp, prev_log2);
    common_log2_magnitudes(cur_mp, cur_log2);
    phase_log2_amplitudes(prev_mp, prev_phase_log2);
    phase_log2_amplitudes(cur_mp, cur_phase_log2);
    for (int l = 1; l <= model->L; l++) {
        const float w = model->w0 * (float)l;
        int state = (columns != NULL) ? column_band_state(prev_mp, cur_mp, columns, model, l)
                                      : band_state(voicing, harmonic_at(voicing, w));
        if (late_onset && state == BAND_VOICED) {
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

/* The mid-frame model's bands (fill_mid_bands()), its pitch and harmonic count
 * being set: by the columns where they apply (column_rule()), else by the
 * voicing of the frame that stands at mid-frame, this frame's but in the
 * rate-33 rendering's hold (rate33_voicing_holds()). */
static void
voice_mid_model(const mbe_parms* prev_mp, const mbe_parms* cur_mp, int late_onset, enum mbe_voiced_synthesis rendering,
                struct mbe_subframe_model* model) {
    struct column_rule columns;
    memset(&columns, 0, sizeof(columns));
    const int by_columns = column_rule(prev_mp, cur_mp, rendering, &columns);
    const int held = (rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450) && rate33_voicing_holds(prev_mp, cur_mp);
    fill_mid_bands(prev_mp, cur_mp, held ? prev_mp : cur_mp, late_onset, by_columns ? &columns : NULL,
                   mid_row(cur_mp, rendering), model);
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
    voice_mid_model(prev_mp, cur_mp, late_onset, rendering, model);
}

/* ---- the pulse train ---- */

struct pulse_train {
    double phase;             /* cycles since the last pulse, at the point reached */
    int count;                /* pulses that can still sound */
    float time[TRAIN_PULSES]; /* in samples from the frame's start (negative: before it) */
};

/* Whether the train runs on from the last frame (`prev_mp`, the frame as it
 * was synthesized, whose rendered model is `prev`) in `rendering`: from a frame
 * whose rendered model has a voiced band; from one marked
 * MBE_FRAME_MARK_PITCH_ABSENT or MBE_FRAME_MARK_PREVIOUS_HELD, through whose
 * first half the voicing before it may have lasted (every rendering); and, in
 * the rate-33 rendering only, from one marked MBE_FRAME_MARK_MID_COLUMN_0 (V/UV
 * row 17), whose mid-frame model voices column 0 (mid_row()), so that the
 * train ran through its first half and out of the mid-frame model over its
 * second (advance_interval()). That is the train store_pulse_train() left.
 * The new case is a row-17 frame rendered with no voiced band: one whose end
 * the adaptive smoothing of channel errors voices carried the train before
 * the row-17 train change as after. Whether the carried pulses sound depends
 * on the frame after: they do when its mid-frame model is voiced (a frame with
 * a voiced band of its own there, or another row-17 frame), and a frame whose
 * first subframe has no voiced band at either end starts afresh. From the
 * mid-frame voicing change, which left row 17's own model with no voiced
 * band, to the row-17 train change, a frame voiced at mid-frame after it
 * started a new train at its first sample. Of the 12 row-17 voice frames of
 * DVSI's development vectors (searched for on development and validation;
 * validation holds none) rendered with no voiced band between two frames with
 * a voiced band, frames 2 to n - 5 of each file, the voiced run after the
 * frame (its first 12 frames at most) comes nearer DVSI's decode by more than
 * 0.02 of correlation with the train running on in 6 and farther in none
 * (docs/dvsi-parity.md, "Rate 33: rows 13, 15 and 17"). A lost-frame repeat
 * does not carry the mark (ambe2450_load_repeat()). These marks of the
 * previous model are read as the caller's state holds them, with no check of
 * where they came from; each decides only whether the train that state stored
 * is loaded, and load_pulse_train() bounds what is loaded. (No mark of the
 * previous model chooses which bands are voiced: column_rule().) In D-STAR's
 * and P25's renderings MBE_FRAME_MARK_MID_COLUMN_0 on the previous model plays
 * as no mark does: the samples are those of the same state with no mark. */
static int
train_runs_on(const mbe_parms* prev_mp, const struct mbe_subframe_model* prev, enum mbe_voiced_synthesis rendering) {
    if (prev->any_voiced) {
        return 1;
    }
    const enum mbe_frame_mark mark = mbe_frame_mark(prev_mp);
    if (mark == MBE_FRAME_MARK_PITCH_ABSENT || mark == MBE_FRAME_MARK_PREVIOUS_HELD) {
        return 1;
    }
    return rendering == MBE_VOICED_SYNTHESIS_DVSI_AMBE2450 && mark == MBE_FRAME_MARK_MID_COLUMN_0;
}

/* The train a stream carries into a frame: the pulse phase reached at the last
 * frame's end and the pulses it kept, when it runs on from that frame
 * (train_runs_on()); after any other frame it starts afresh, whatever a codec
 * path left in the state.
 * What is read is bounded, not checked to be what a frame stores: a state
 * without this rendering's subframe mark (another rendering's noise overlap)
 * loads nothing; a pulse count whose bit pattern lies above 32.0f's (a
 * negative or not-a-number one included) loads nothing, and a count inside the
 * range is truncated to a whole number (a caller's 1.5 plays as 1); a phase
 * outside [0, 1) by its bit pattern reads as 0; a pulse time that is not
 * finite, is positive, or lies more than a frame back is skipped. */
static void
load_pulse_train(const mbe_parms* prev_mp, const struct mbe_subframe_model* prev, enum mbe_voiced_synthesis rendering,
                 struct pulse_train* train) {
    const float stored = prev_mp->noiseOverlap[0];
    const float phase = prev_mp->noiseOverlap[PHASE_SLOT];
    uint32_t bits;

    train->phase = 0.0;
    train->count = 0;
    if (!mbe_subframe_state_is_marked(prev_mp) || !train_runs_on(prev_mp, prev, rendering)) {
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

/* What a frame's pulsed state leaves the next (pulsed_frame()): which of the
 * first two models emit, whether this frame's code was read, and P25's
 * validity counter (0 in the other renderings). */
struct pulsed_reading {
    int emits[2];
    int read;
    int count;
};

/* What the frame leaves the next: the pulse phase and the pulses of the last
 * frame's length, most recent first, in noiseOverlap, with its pulsed state
 * (`pulsed`, pulsed_frame()): whether its mid-frame model emitted a pulse
 * (PULSED_SLOT) and, in P25's rendering, whether its code was read
 * (PULSED_READ_SLOT) and the validity counter (PULSED_COUNT_SLOT); and, as
 * the TIA-102.BABA rendering leaves
 * them, harmonic l's phase at the frame's end before the regenerated phase
 * (PSIl) and with it (PHIl). Element 0 of PSIl and PHIl is not touched (the
 * frame's mark lives there). */
static void
store_pulse_train(mbe_parms* cur_mp, const struct mbe_subframe_model* own, const struct pulse_train* train,
                  const struct pulsed_reading* pulsed) {
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
    cur_mp->noiseOverlap[PULSED_SLOT] = pulsed->emits[1] ? PULSED_EMITTED : 0.0f;
    cur_mp->noiseOverlap[PULSED_READ_SLOT] = pulsed->read ? PULSED_READ : 0.0f;
    cur_mp->noiseOverlap[PULSED_COUNT_SLOT] = (float)pulsed->count;
    mbe_subframe_state_mark(cur_mp);
}

/* ---- rendering the voiced part ---- */

/* A model's voiced bands, or its pulsed ones, ready to render. */
struct rendered_model {
    int count;        /* bands */
    float w0;         /* radians per sample */
    float period;     /* samples */
    float real[56];   /* amplitude times the cosine of the band's phase */
    float imag[56];   /* amplitude times its sine */
    float step_c[56]; /* cosine of the band's advance per sample */
    float step_s[56]; /* sine of it */
    short harmonic[56];
};

/* A model's bands in `state` (BAND_VOICED, played by the voiced pulse train,
 * or BAND_PULSED, played by the pulse of the pulsed state), at the voiced
 * render gain: the same response for both. */
static void
prepare_render(const struct mbe_subframe_model* model, float voiced_gain, int state, struct rendered_model* render) {
    const float hz_per_harmonic = model->w0 * (8000.0f / (2.0f * (float)M_PI));
    const int any = (state == BAND_PULSED) ? model->any_pulsed : model->any_voiced;

    render->count = 0;
    render->w0 = model->w0;
    render->period = any ? (2.0f * (float)M_PI) / model->w0 : 0.0f;
    for (int l = 1; l <= model->L && any; l++) {
        const float w = model->w0 * (float)l;
        /* A harmonic at or above Nyquist would alias back into the band. */
        if (model->voiced[l] != state || !(w < (float)M_PI)) {
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

/* One interval of the voiced pulse train: every pulse plays the voiced bands
 * of the model at the interval's start going out and of the one at its end
 * coming in. */
static void
render_interval(float* out, float start, const struct pulse_train* train, const struct rendered_model* from,
                const struct rendered_model* to) {
    for (int i = 0; i < train->count; i++) {
        const float pulse = train->time[i] - start;
        render_pulse(out, from, pulse, 0);
        render_pulse(out, to, pulse, 1);
    }
}

/* The pulse a model of the pulsed state emits, MBE_PULSED_DELAY samples into
 * the interval it starts (pulsed_frame()): through its own pulsed bands under
 * its falling crossfade and through the next model's pulsed bands under its
 * rising one, the voiced pulse response (render_pulse()), over this interval
 * only: it does not ring on past the next model's point. The voiced train's
 * pulses play no pulsed band and this pulse no voiced band. */
static void
render_emitted_pulse(float* out, const struct rendered_model* from, const struct rendered_model* to) {
    render_pulse(out, from, (float)MBE_PULSED_DELAY, 0);
    render_pulse(out, to, (float)MBE_PULSED_DELAY, 1);
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

/* Whether the mid-frame model of the frame a state holds emitted its pulse of
 * the pulsed state (PULSED_SLOT): only a state with this rendering's subframe
 * mark holding exactly PULSED_EMITTED there, by its bit pattern (0x3F800000);
 * anything else a caller's state holds (not a number, an infinity, a
 * subnormal, a negative value, 0.5, 2.0) is "not emitted", and so is a state
 * the other rendering wrote. */
static int
mid_model_emitted(const mbe_parms* prev_mp) {
    return stored_exactly(prev_mp, PULSED_SLOT, PULSED_EMITTED);
}

/*
 * The pulsed state of a frame's three models (docs/dvsi-parity.md, "D-STAR and
 * P25: the pulsed state"), in D-STAR's and P25's renderings: the last frame's
 * own model under its code's own-subframe pattern, the mid-frame model under
 * this frame's mid-frame pattern and this frame's own model under its
 * own-subframe pattern (pulse_bands(); the codes from the models themselves,
 * pulsed_code(), where the rendering reads them: in P25's only while the
 * validity counter allows, read_code(), and the last frame's code only where
 * that frame read it, last_read_code(); a code that is not read is no code,
 * and the frame renders as it would at another pitch with no voiced band). The
 * mid-frame model that a frame marked
 * MBE_FRAME_MARK_PREVIOUS_HELD or MBE_FRAME_MARK_TONE puts in place (the
 * previous model held, or a tone frame's own model) takes no pattern. A model
 * with a pulsed band emits one pulse MBE_PULSED_DELAY samples after its point
 * unless the model before it emitted one: emits[0] for the last frame's own
 * model (the one before it is that frame's mid-frame model, whose emission the
 * state carries, mid_model_emitted()), emits[1] for the mid-frame model; this
 * frame's own model emits, if at all, in the next frame. A voiced model never
 * counts as having emitted (it has no pulsed band), nor does one with no
 * pulsed band. Starting from silence, the last frame's own model plays
 * nothing, so it has no pulsed band and emits nothing. Measured on DVSI's
 * decodes by a scratch build of the library: the alternation, against a model
 * before a voiced one emitting nothing and against a voiced model counting as
 * having emitted. The 24 samples were first fitted on the study's numerical
 * model of P25 frames deep inside runs with no voiced band (against 23 and 25)
 * and measured again, on DVSI's decodes, by scratch builds of the library
 * after the gate sets were declared (against 22, 23, 25 and 26).
 */
static void
pulsed_frame(const mbe_parms* cur_mp, const mbe_parms* prev_mp, enum mbe_voiced_synthesis rendering,
             struct mbe_subframe_model models[3], struct pulsed_reading* pulsed) {
    const int code = read_code(cur_mp, prev_mp, rendering);
    const enum mbe_frame_mark mark = mbe_frame_mark(cur_mp);
    pulse_bands(&models[0], last_read_code(prev_mp, rendering), PULSED_OWN);
    if (mark != MBE_FRAME_MARK_PREVIOUS_HELD && mark != MBE_FRAME_MARK_TONE) {
        pulse_bands(&models[1], code, PULSED_MID);
    }
    pulse_bands(&models[2], code, PULSED_OWN);
    pulsed->emits[0] = models[0].any_pulsed && !mid_model_emitted(prev_mp);
    pulsed->emits[1] = models[1].any_pulsed && !pulsed->emits[0];
    pulsed->read = (rendering == MBE_VOICED_SYNTHESIS_DVSI_IMBE) && code >= 0;
    pulsed->count = next_count(cur_mp, prev_mp, rendering);
}

/* A frame's three models (the last frame's own, mid-frame, this frame's own),
 * their pulsed state (pulsed_frame(): which of the first two emit, what the
 * frame leaves the next) and the
 * pulse train it starts from, for rendering and advancing alike.
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
              struct mbe_subframe_model models[3], struct pulse_train* train, struct pulsed_reading* pulsed) {
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
        load_pulse_train(prev_mp, &models[0], rendering, train);
    }
    pulsed_frame(cur_mp, prev_mp, rendering, models, pulsed);
    return 1;
}

/* In P25's rendering, the mark of the model rendered, as it stands after
 * adaptive smoothing (which can voice harmonics of a frame decoded with none):
 * MBE_FRAME_MARK_ONSET_COLUMNS when it has a voiced band and is not a tone
 * model (MBE_FRAME_MARK_TONE), whatever mark it was decoded or repeated with,
 * so that a caller's MBE_FRAME_MARK_LATE_ONSET or MBE_FRAME_MARK_PREVIOUS_HELD
 * there is overwritten and has no effect; with no voiced band, not that mark
 * (none instead), any other mark kept. The rendering reads its voiced models'
 * bands in the 500 Hz columns whatever their mark (column_rule()), so
 * MBE_FRAME_MARK_ONSET_COLUMNS no longer chooses bands; this function still
 * reads it (to clear it from a model with no voiced band) and writes it, and
 * the marking still decides that no other mark acts on a voiced model. The prediction history the data path
 * committed before keeps the decoded model and its mark. Other renderings
 * leave the mark as it is. Called only on a frame the synthesis accepts
 * (frame_accepted()), so that a rejected call leaves the state untouched. */
static void
mark_rendered_model(mbe_parms* cur_mp, enum mbe_voiced_synthesis rendering) {
    if (!imbe_rendering(rendering)) {
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
    struct rendered_model renders[3], pulsed[3];
    struct mbe_unvoiced_subframe noise[2];
    struct pulse_train train;
    struct pulsed_reading reading;

    if (!frame_accepted(cur_mp, prev_mp)) {
        return;
    }
    mark_rendered_model(cur_mp, rendering);
    if (!prepare_frame(cur_mp, prev_mp, from_silence, rendering, models, &train, &reading)) {
        return;
    }
    for (int i = 0; i < 3; i++) {
        prepare_render(&models[i], voiced_gain, BAND_VOICED, &renders[i]);
        prepare_render(&models[i], voiced_gain, BAND_PULSED, &pulsed[i]);
    }
    for (int k = 0; k < 2; k++) {
        const float start = (float)(k * SUBFRAME);
        advance_interval(&train, start, &models[k], &models[k + 1]);
        render_interval(&aout_buf[(size_t)k * SUBFRAME], start, &train, &renders[k], &renders[k + 1]);
        if (reading.emits[k]) {
            render_emitted_pulse(&aout_buf[(size_t)k * SUBFRAME], &pulsed[k], &pulsed[k + 1]);
        }
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
    store_pulse_train(cur_mp, &models[2], &train, &reading);
}

void
mbe_advance_subframes(mbe_parms* cur_mp, const mbe_parms* prev_mp, int from_silence,
                      enum mbe_voiced_synthesis rendering) {
    struct mbe_subframe_model models[3];
    struct pulse_train train;
    struct pulsed_reading reading;

    if (!frame_accepted(cur_mp, prev_mp)) {
        return;
    }
    mark_rendered_model(cur_mp, rendering);
    if (!prepare_frame(cur_mp, prev_mp, from_silence, rendering, models, &train, &reading)) {
        return;
    }
    advance_interval(&train, 0.0f, &models[0], &models[1]);
    advance_interval(&train, (float)SUBFRAME, &models[1], &models[2]);
    mbe_dvsi_noise_skip_frame(cur_mp);
    /* Where rendering the frame leaves it, the pulsed state too (the emission,
     * the code read, P25's validity counter). */
    store_pulse_train(cur_mp, &models[2], &train, &reading);
}
