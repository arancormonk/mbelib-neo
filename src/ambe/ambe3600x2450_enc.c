// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief AMBE+2 3600x2450 (DMR, NXDN, YSF, P25 Phase 2) speech encoder.
 *
 * Quantizes the speech model of mbe_speech_analysis.c as TIA-102.BABA-1
 * clause 4 describes (clause numbers below refer to that addendum):
 *
 *  - b0: the Annex A fundamental nearest the analysed one (4.1).
 *  - b1: the V/UV vector among the first 17 of Annex B nearest the band
 *    decisions of TIA-102.BABA 5.2, weighted by the squared spectral
 *    amplitudes (4.2, eq 4).
 *  - Spectral amplitudes by the estimator of each band's decision
 *    (TIA-102.BABA 5.3, eqs 43 and 44), taken to the log domain with the
 *    transmitted voicing (eq 8), then the gain, prediction residual, PRBA and
 *    higher-order coefficients of 4.3 against the decoder's previous frame
 *    (ambe_encoder.c).
 *  - DTMF, KNOX, call-progress and single tones are sent as tone frames (7.2).
 *
 * A frame whose gain falls below the first step is flattened toward the level
 * that step can play (mbe_encoder_fit_floor()), as D-STAR's are; the standard
 * leaves such frames to the nearest step, which plays a quiet tone far too
 * loud. cur_mp is what mbe_decodeAmbe2450Parms() decodes from the emitted bits
 * given prev_mp.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ambe3600x2400_internal.h"
#include "ambe3600x2450_const.h"
#include "ambe3600x2450_internal.h"
#include "ambe_encoder.h"
#include "mbe_encoder.h"
#include "mbe_result.h"
#include "mbe_speech_analysis.h"
#include "mbe_tone.h"
#include "mbe_tone_detect.h"
#include "mbelib-neo/mbelib.h"

#define AMBE2450_VUV_CODES     17 /* 4.2: only the first 17 vectors are searched */

/*
 * Level: AMBE2450_ENC_MAG_SCALE scales the analysis amplitudes so the decoded
 * level matches DVSI's encoding of the same input, both decoded by this
 * library: the median, over DVSI's rate-33 development vectors, of the level
 * that puts our encoding at DVSI's (as AMBE2400_ENC_MAG_SCALE for D-STAR).
 */
#define AMBE2450_ENC_MAG_SCALE 1.08f

struct mbe_ambe2450_encoder {
    struct mbe_encoder_frontend fe;
    unsigned char bands_prev[MBE_ANALYSIS_BANDS]; /* last voice frame's band decisions, eq 37 */
    struct mbe_tone_tracker tones;
};

static const struct ambe_enc_tables ambe2450_enc_tables = {
    .dg = AmbeDg,
    .dg_count = 32,
    .dg_scale = 1.0f,
    .dg_offset = 0.0f,
    .dg_history = 0.5f,
    .prba24 = AmbePRBA24,
    .prba58 = AmbePRBA58,
    .hoc = {AmbeHOCb5, AmbeHOCb6, AmbeHOCb7, AmbeHOCb8},
    .hoc_rows = {32, 16, 16, 8},
    .hoc_step = {1, 1, 1, 1},
    .lmprbl = AmbeLmprbl,
    .rho = 0.65f,
    .fit_floor_gain = 1, /* as D-STAR: see mbe_encoder_fit_floor() */
};

void
mbe_ambe2450EncoderReset(mbe_ambe2450_encoder* enc) {
    if (enc == NULL) {
        return;
    }
    mbe_encoder_frontend_reset(&enc->fe);
    memset(enc->bands_prev, 0, sizeof(enc->bands_prev));
    mbe_tone_tracker_reset(&enc->tones);
}

mbe_ambe2450_encoder*
mbe_ambe2450EncoderAlloc(void) {
    mbe_ambe2450_encoder* enc = calloc(1, sizeof(*enc));
    if (enc == NULL) {
        return NULL;
    }
    if (mbe_encoder_frontend_open(&enc->fe) != 0) {
        free(enc);
        return NULL;
    }
    mbe_ambe2450EncoderReset(enc);
    return enc;
}

void
mbe_ambe2450EncoderFree(mbe_ambe2450_encoder* enc) {
    if (enc != NULL) {
        mbe_encoder_frontend_close(&enc->fe);
        free(enc);
    }
}

int
mbe_ambe2450_pitch_code(float f0) {
    if (!(f0 > 0.0f)) {
        return MBE_AMBE2450_PITCH_CODES - 1;
    }
    long b0 = lroundf((MBE_AMBE2450_LOG2_F0_OFFSET - log2f(f0)) / MBE_AMBE2450_LOG2_F0_SLOPE);
    return (int)((b0 < 0) ? 0 : ((b0 >= MBE_AMBE2450_PITCH_CODES) ? MBE_AMBE2450_PITCH_CODES - 1 : b0));
}

/* eq 5: the V/UV vector element of harmonic l at fundamental f0q. */
static int
ambe2450_enc_column(float f0q, int l) {
    int jl = (int)((float)l * 16.0f * f0q);
    return (jl > 7) ? 7 : jl;
}

int
mbe_ambe2450_quantize_vuv(int b0, const float m[57], const unsigned char bands[MBE_ANALYSIS_BANDS]) {
    /* 4.2 eq 4: the vector minimizing the squared-amplitude weighted
     * disagreement with the band decisions; ties keep the lower index. */
    const float f0q = mbe_ambe2450_f0(b0);
    const int L = mbe_ambe2450_harmonic_count(b0);
    int best = 0;
    float best_distance = 0.0f;
    for (int n = 0; n < AMBE2450_VUV_CODES; n++) {
        float distance = 0.0f;
        for (int l = 1; l <= L; l++) {
            int voiced = bands[mbe_analysis_band_of(l) - 1];
            if (AmbeVuv[n][ambe2450_enc_column(f0q, l)] != voiced) {
                distance += m[l] * m[l];
            }
        }
        if (n == 0 || distance < best_distance) {
            best_distance = distance;
            best = n;
        }
    }
    return best;
}

/* eq 8 less its 0.5 log2 L, which the gain adds: an unvoiced amplitude is
 * raised by the inverse of the decoder's unvoiced scale 0.2046 / sqrt(w0). */
static void
ambe2450_enc_log_amplitudes(int b0, int b1, const float m[57], float a[57]) {
    const float f0q = mbe_ambe2450_f0(b0);
    const float unvoiced_offset = (0.5f * log2f(f0q * 2.0f * (float)M_PI)) + 2.289f;
    for (int l = 1; l <= mbe_ambe2450_harmonic_count(b0); l++) {
        int voiced = AmbeVuv[b1][ambe2450_enc_column(f0q, l)];
        a[l] = log2f(m[l] + 1e-12f) + (voiced ? 0.0f : unvoiced_offset);
    }
}

void
mbe_ambe2450_quantize_amplitudes(int b0, int b1, const float a[57], const mbe_parms* prev_mp, char ambe_d[49]) {
    struct ambe_enc_frame q = {0};
    q.cache = mbe_ambe2400_get_dct_cache();
    q.b[0] = b0;
    q.b[1] = b1;
    q.f0q = mbe_ambe2450_f0(b0);
    q.L = mbe_ambe2450_harmonic_count(b0);
    memcpy(q.a, a, sizeof(q.a));
    ambe_enc_mean_amplitude(&q);
    ambe_enc_quantize_gain(&q, &ambe2450_enc_tables, prev_mp);
    ambe_enc_prediction(&q, &ambe2450_enc_tables, prev_mp);
    ambe_enc_quantize_residual(&q, &ambe2450_enc_tables);
    mbe_ambe2450_pack(q.b, ambe_d);
}

/* Quantize and pack a voice frame; cur_mp is the decoder's model of it. */
static int
ambe2450_enc_voice(mbe_ambe2450_encoder* enc, const struct mbe_analysis_result* res, char ambe_d[49], mbe_parms* cur_mp,
                   const mbe_parms* prev_mp) {
    unsigned char bands[MBE_ANALYSIS_BANDS];
    float m[57] = {0};
    float a[57] = {0};

    const int b0 = mbe_ambe2450_pitch_code(res->f0);
    const int L = mbe_ambe2450_harmonic_count(b0);
    mbe_analysis_band_voicing(res, L, enc->bands_prev, bands);
    for (int l = 1; l <= L; l++) {
        int voiced = bands[mbe_analysis_band_of(l) - 1];
        m[l] = AMBE2450_ENC_MAG_SCALE * (voiced ? res->voiced_magnitude[l] : res->noise_magnitude[l]);
    }
    const int b1 = mbe_ambe2450_quantize_vuv(b0, m, bands);
    ambe2450_enc_log_amplitudes(b0, b1, m, a);
    mbe_ambe2450_quantize_amplitudes(b0, b1, a, prev_mp, ambe_d);

    /* The decoder's update, on a private copy of the caller's history. */
    mbe_parms history = *prev_mp;
    int kind = mbe_decodeAmbe2450Parms(ambe_d, cur_mp, &history);
    if (kind != MBE_AMBE2450_FRAME_VOICE) {
        return (kind < 0) ? kind : MBE_STATUS_INVALID_ARGUMENT;
    }
    if (!mbe_encoder_model_valid(cur_mp)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    memcpy(enc->bands_prev, bands, sizeof(enc->bands_prev));
    return MBE_AMBE2450_FRAME_VOICE;
}

/*
 * 7.2 Table 10: u0 = 63 then AD(6..1); the tone index four times from bit 12
 * (u1, u2 and u3 hold it whole or split, which also fills b0's low bits); AD(0)
 * at bit 44; the last four bits 0.
 */
static void
ambe2450_enc_pack_tone(int tone_id, int ad, char ambe_d[49]) {
    memset(ambe_d, 0, 49);
    for (int i = 0; i < 6; i++) {
        ambe_d[i] = 1;
        ambe_d[6 + i] = (char)((ad >> (6 - i)) & 1);
    }
    for (int i = 0; i < 32; i++) {
        ambe_d[12 + i] = (char)((tone_id >> (7 - (i % 8))) & 1);
    }
    ambe_d[44] = (char)(ad & 1);
}

/* AD for a component amplitude (peak, 16-bit scale) from the inverse of the
 * level law the decoder plays tone frames with (mbe_tone_ambe2450_level_db).
 * DVSI's encoder sends half a step more: on its rate-33 tone vectors this
 * rounding gives DVSI's AD on about 72% of tone frames and one step off on
 * nearly all others. */
static int
ambe2450_enc_tone_level(float amplitude) {
    double rms_db = 20.0 * log10(fmax((double)amplitude, 1e-3) / (M_SQRT2 * 32768.0));
    long ad = lround(127.5 + ((rms_db - MBE_TONE_AMBE2450_DB_AT_AD127) / MBE_TONE_AMBE2450_DB_PER_STEP));
    return (int)((ad < 0) ? 0 : ((ad > 127) ? 127 : ad));
}

static int
ambe2450_encode(mbe_ambe2450_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                const mbe_parms* prev_mp) {
    struct mbe_analysis_result res;
    int status = mbe_encoder_frontend_analyze(&enc->fe, samples, &res);
    if (status < 0) {
        return status;
    }
    struct mbe_tone_detection tone;
    status = mbe_tone_track(&enc->tones, enc->fe.fft,
                            mbe_analysis_span(&enc->fe.analysis, MBE_TONE_SPAN, MBE_TONE_OFFSET), &tone);
    if (status < 0) {
        return status;
    }
    if (status > 0) {
        /* A tone frame carries no speech model and leaves the prediction
         * history alone (4.3); cur_mp mirrors that. */
        ambe2450_enc_pack_tone(tone.id, ambe2450_enc_tone_level(tone.amplitude), ambe_d);
        *cur_mp = *prev_mp;
        return mbe_encoder_model_valid(cur_mp) ? MBE_AMBE2450_FRAME_TONE : MBE_STATUS_INVALID_ARGUMENT;
    }
    return ambe2450_enc_voice(enc, &res, ambe_d, cur_mp, prev_mp);
}

int
mbe_encodeAmbe2450Parms(mbe_ambe2450_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                        const mbe_parms* prev_mp) {
    if (enc == NULL || samples == NULL || ambe_d == NULL || cur_mp == NULL || prev_mp == NULL || cur_mp == prev_mp
        || !mbe_encoder_samples_valid(samples) || !mbe_encoder_history_valid(prev_mp)) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    /* A frame that fails (a history whose model would overflow) leaves no trace. */
    const struct mbe_analysis_state saved = enc->fe.analysis;
    const struct mbe_tone_tracker saved_tones = enc->tones;
    mbe_parms out = *cur_mp;
    int status = ambe2450_encode(enc, samples, ambe_d, &out, prev_mp);
    if (status < 0) {
        enc->fe.analysis = saved;
        enc->tones = saved_tones;
        return status;
    }
    *cur_mp = out;
    return status;
}

int
mbe_encodeAmbe2450ParmsShort(mbe_ambe2450_encoder* enc, const short* samples, char ambe_d[49], mbe_parms* cur_mp,
                             const mbe_parms* prev_mp) {
    float float_buf[MBE_ENCODER_SAMPLES];

    if (samples == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    mbe_encoder_short_to_float(samples, float_buf);
    return mbe_encodeAmbe2450Parms(enc, float_buf, ambe_d, cur_mp, prev_mp);
}

int
mbe_encodeAmbe3600x2450Frame(const char ambe_d[49], char ambe_fr[4][24]) {
    if (ambe_d == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    int bit_status = mbe_validate_bits(ambe_d, 49u);
    if (bit_status < 0) {
        return bit_status;
    }
    (void)ambe_enc_fec(ambe_d, ambe_fr);
    return 0;
}
