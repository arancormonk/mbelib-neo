// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 *
 * Copyright (C) 2010 mbelib Author
 * GPG Key ID: 0xEA5EFE2C (9E7A 5527 9CDC EBF7 BF1B  D772 4F98 E863 EA5E FE2C)
 *
 * Portions were originally under the ISC license; this mbelib-neo
 * distribution is provided under GPL-2.0-or-later. See LICENSE for details.
 */

/**
 * @file
 * @brief AMBE 3600x2400 (D-STAR DV) speech encoder.
 *
 * Performs the inverse of mbe_decodeAmbe2400Parms():
 *
 *  - 20 ms of 8 kHz PCM is windowed and analysed (pitch, voicing,
 *    harmonic spectral amplitudes) to produce the same parameters the
 *    decoder consumes.
 *  - Every parameter is quantized against the exact tables the decoder
 *    dequantizes from (AmbePlusLtable/AmbePlusVuv/AmbePlusDg/
 *    AmbePlusPRBA24/AmbePlusPRBA58/AmbePlusHOCb5..b8), so the
 *    reconstructed frame is bit-compatible with this library's
 *    mbe_decodeAmbe2400Parms()/mbe_processAmbe3600x2400*() path and follows
 *    the D-STAR AMBE bit layout (interleave, scrambler and Golay parity
 *    cross-checked against the MMDVM tables). Interoperability with DVSI
 *    hardware has not been verified.
 *  - The prediction state (log2Ml, gamma) is advanced with the
 *    QUANTIZED values, mirroring the decoder state so the encoder and
 *    decoder never drift apart.
 *
 * Tone (DTMF) encoding is not supported; silent input produces the
 * standard AMBE silence frame (b0 == 127, tone index 128).
 */

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "ambe3600x2400_const.h"
#include "ambe3600x2400_internal.h"
#include "mbe_ecc.h"
#include "mbe_unvoiced_fft.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"

#define AMBE2400_ENC_FFT_SIZE    256
#define AMBE2400_ENC_SAMPLES     160
#define AMBE2400_ENC_SILENCE_RMS 0.0015f

/* f0 = exp2(-4.311767578125 - 2.1336e-2 * (b0 + 0.5)) */
#define AMBE2400_ENC_F0_OFFSET   (-4.311767578125f)
#define AMBE2400_ENC_F0_STEP     (-0.021336f)

/*
 * Calibration: the AMBE 2400 gain field (AmbePlusDg) is non-negative and
 * spans about 5.35 log2 units (32 dB), so the encoder keeps a fixed nominal
 * level. An input AGC normalizes every frame to AMBE2400_ENC_AGC_TARGET
 * RMS and this scale places that nominal level in the middle of the
 * codec's magnitude range. 77.0 was tuned empirically so that encoded
 * frames decode at a level comparable to off-air reference recordings,
 * keeping the gain field within the range its codebook spans.
 */
#define AMBE2400_ENC_MAG_SCALE   77.0f
#define AMBE2400_ENC_AGC_TARGET  0.08f
#define AMBE2400_ENC_AGC_ALPHA   0.9f
#define AMBE2400_ENC_AGC_MIN     0.1f
#define AMBE2400_ENC_AGC_MAX     10.0f

/*
 * Caller-owned analysis state: one context per stream.
 *
 * history holds the previous 160 samples so the encoder can analyse
 * a 256 sample (32 ms) window centred on the current frame start without
 * blocking for look-ahead.
 */
struct mbe_ambe2400_encoder {
    float history[AMBE2400_ENC_SAMPLES];
    float agc_gain;
    float agc_rms;
    float hist_gain;
    int silence_run;
    float prev_lag;
    float max_energy;
    mbe_fft_plan* fft;
};

void
mbe_ambe2400EncoderReset(mbe_ambe2400_encoder* enc) {
    if (enc == NULL) {
        return;
    }
    mbe_fft_plan* fft = enc->fft;
    memset(enc, 0, sizeof(*enc));
    enc->fft = fft;
    enc->agc_gain = 1.0f;
    enc->agc_rms = AMBE2400_ENC_AGC_TARGET;
    enc->hist_gain = 1.0f;
    enc->max_energy = 1e-6f;
}

mbe_ambe2400_encoder*
mbe_ambe2400EncoderAlloc(void) {
    mbe_ambe2400_encoder* enc = calloc(1, sizeof(*enc));
    if (enc == NULL) {
        return NULL;
    }
    enc->fft = mbe_fft_plan_alloc();
    if (enc->fft == NULL) {
        free(enc);
        return NULL;
    }
    mbe_ambe2400EncoderReset(enc);
    return enc;
}

void
mbe_ambe2400EncoderFree(mbe_ambe2400_encoder* enc) {
    if (enc != NULL) {
        mbe_fft_plan_free(enc->fft);
        free(enc);
    }
}

/*
 * Compute log2 spectral magnitudes and raw per-harmonic voicing from a
 * pitch estimate over the windowed analysis buffer.
 */
static int
ambe2400_enc_spectrum(mbe_ambe2400_encoder* enc, const float* windowed, float f0q, int L, float mag[57],
                      int vl_ana[57]) {
    float fft_out[AMBE2400_ENC_FFT_SIZE];
    const float f0_bin = f0q * (float)AMBE2400_ENC_FFT_SIZE;
    float max_band_energy = 0.0f;

    int status = mbe_fft_forward_real(enc->fft, windowed, fft_out);
    if (status < 0) {
        return status;
    }

    /* Ordered real FFT: [DC, Nyquist, re1, im1, ..., re127, im127]. */
    for (int l = 1; l <= L; l++) {
        float lo = ((float)l - 0.5f) * f0_bin;
        float hi = ((float)l + 0.5f) * f0_bin;
        int lo_bin = (int)lo;
        int hi_bin = (int)hi + 1;
        float band_energy = 0.0f;
        float band_peak = 0.0f;
        int n_bins = 0;

        if (lo_bin < 1) {
            lo_bin = 1;
        }
        if (hi_bin > (AMBE2400_ENC_FFT_SIZE / 2) - 1) {
            hi_bin = (AMBE2400_ENC_FFT_SIZE / 2) - 1;
        }

        for (int b = lo_bin; b <= hi_bin; b++) {
            float re = fft_out[2 * (size_t)b];
            float im = fft_out[(2 * b) + 1];
            float e = (re * re) + (im * im);
            band_energy += e;
            n_bins++;
            if (e > band_peak) {
                band_peak = e;
            }
        }

        mag[l] = sqrtf(band_energy + 1e-12f);
        if (band_energy > max_band_energy) {
            max_band_energy = band_energy;
        }
        /* Voiced when the dominant bin stands clear of the in-band
         * background (window leakage raises the background uniformly).
         * Thresholds are deliberately loose: too strict a criterion makes
         * soft/sloped frames drop to "unvoiced", which sounds like the voice
         * cutting out. The VUV codebook + hysteresis smooth the rest. */
        {
            float background = (n_bins > 1) ? ((band_energy - band_peak) / (float)(n_bins - 1)) : band_energy;
            vl_ana[l] = (band_peak > (1.35f * background)) && (band_energy > (0.0008f * max_band_energy + 1e-12f));
        }
    }
    return 0;
}

static float
ambe2400_enc_pitch_strength(const float* buf, int n, float lag_f) {
    float strength;
    /* Voicing strength: normalized autocorrelation at the chosen lag.
     * ~1.0 for clean periodic speech, ~0.0 for noise. */
    int lag = (int)(lag_f + 0.5f);
    if (lag < 20) {
        lag = 20;
    }
    if (lag > 127) {
        lag = 127;
    }
    float num = 0.0f;
    float den_a = 0.0f;
    float den_b = 0.0f;
    for (int i = 0; i + lag < n; i++) {
        num += buf[i] * buf[i + lag];
        den_a += buf[i] * buf[i];
        den_b += buf[i + lag] * buf[i + lag];
    }
    float den = sqrtf(den_a * den_b);
    float r = (den > 1e-12f) ? (num / den) : 0.0f;
    strength = r;
    if (strength < 0.0f) {
        strength = 0.0f;
    }
    if (strength > 1.0f) {
        strength = 1.0f;
    }
    return strength;
}

static float
ambe2400_enc_refine_lag(const float amdf[128], int lag) {
    /* Parabolic refinement */
    float denom = 0.0f;
    if (lag > 20 && lag < 127) {
        denom = amdf[lag - 1] + amdf[lag + 1] - (2.0f * amdf[lag]);
    }
    float lag_f = (float)lag;
    if (fabsf(denom) > 1e-12f) {
        float shift = 0.5f * (amdf[lag - 1] - amdf[lag + 1]) / denom;
        if (shift > -1.0f && shift < 1.0f) {
            lag_f += shift;
        }
    }

    return lag_f;
}

static float
ambe2400_enc_pitch_candidate(const mbe_ambe2400_encoder* enc, const float amdf[128], float global_min) {
    const float tol = global_min * 1.4f;
    int have_prev = (enc->prev_lag >= (float)20);

    float best_lag_f = 20.0f;
    float best_score = 1e30f;

    for (int lag = 20; lag <= 127; lag++) {
        if ((lag > 20 && amdf[lag] >= amdf[lag - 1]) || (lag < 127 && amdf[lag] > amdf[lag + 1])) {
            continue;
        }
        if (amdf[lag] > tol) {
            continue;
        }

        float lag_f = ambe2400_enc_refine_lag(amdf, lag);

        float score;
        if (have_prev) {
            /* AMDF quality + continuity (octave-aware). */
            score = (amdf[lag] / (global_min + 1e-12f)) + (0.8f * fabsf(log2f(lag_f / enc->prev_lag)));
        } else {
            /* Cold start: among the tied minima, the shortest period is
             * the true fundamental (kills the octave-down error). */
            score = (amdf[lag] / (global_min + 1e-12f)) + (0.25f * (lag_f / 64.0f));
        }

        if (score < best_score) {
            best_score = score;
            best_lag_f = lag_f;
        }
    }

    return best_lag_f;
}

/*
 * Pitch search: AMDF over lags 20..127 (63..400 Hz) with octave-error
 * protection. For a periodic signal the AMDF has near-equal minima at the
 * true period and its multiples; a plain global minimum therefore often
 * latches onto double the true period (half the pitch). We collect every
 * local minimum within a tolerance of the global minimum and pick among
 * them by continuity (and, on cold start, the shortest period).
 */
static float
ambe2400_enc_pitch(mbe_ambe2400_encoder* enc, const float* buf, int n, float* strength) {
    const int min_lag = 20;
    const int max_lag = 127;
    float amdf[128];
    float global_min = 1e30f;

    for (int lag = min_lag; lag <= max_lag; lag++) {
        float acc = 0.0f;
        int cnt = 0;
        for (int i = 0; i + lag < n; i++) {
            float d = buf[i] - buf[i + lag];
            acc += fabsf(d);
            cnt++;
        }
        amdf[lag] = acc / (float)cnt;
        if (amdf[lag] < global_min) {
            global_min = amdf[lag];
        }
    }

    float best_lag_f = ambe2400_enc_pitch_candidate(enc, amdf, global_min);
    int have_prev = (enc->prev_lag >= (float)min_lag);

    /* Final guard: limit frame-to-frame jumps to 25%. */
    if (have_prev) {
        float max_jump = 0.25f * enc->prev_lag;
        float diff = best_lag_f - enc->prev_lag;
        if (fabsf(diff) > max_jump) {
            best_lag_f = enc->prev_lag + (diff > 0.0f ? max_jump : -max_jump);
        }
    }

    enc->prev_lag = best_lag_f;

    if (strength != NULL) {
        *strength = ambe2400_enc_pitch_strength(buf, n, best_lag_f);
    }
    return best_lag_f;
}

/* Temporary analysis and quantization workspace; no inter-frame state. */
struct ambe2400_enc_frame {
    const struct ambe_dct_cache* cache;
    float buf[AMBE2400_ENC_FFT_SIZE];
    float windowed[AMBE2400_ENC_FFT_SIZE];
    float p[57];
    float mag[57], a[57], Tl[57], Tl_q[57];
    float Cik[5][18], Cik_q[5][18], Gm[9];
    int Vl_ana[57], Ji[5], L, b[9];
    float f0q, gamma_q, mean_a, mean_p, strength;
};

static void
ambe2400_enc_window(const mbe_ambe2400_encoder* enc, struct ambe2400_enc_frame* q, const float* pcm) {
    /* Build the centred analysis window: last 128 of history + first 128 of pcm.
     * The AGC gain is applied per frame; ramp it smoothly across the window so
     * the frame boundary does not introduce a step discontinuity (which smears
     * the spectrum and destroys the voicing decision). */
    for (int i = 0; i < 128; i++) {
        q->buf[i] = enc->history[AMBE2400_ENC_SAMPLES - 128 + i] / enc->hist_gain;
    }
    for (int i = 0; i < 128; i++) {
        q->buf[128 + i] = pcm[i] / enc->agc_gain;
    }
    for (int i = 0; i < AMBE2400_ENC_FFT_SIZE; i++) {
        float t = (float)i / (float)(AMBE2400_ENC_FFT_SIZE - 1);
        float g = enc->hist_gain + (enc->agc_gain - enc->hist_gain) * t;
        q->buf[i] *= g;
    }

    /* DC removal */
    float dc = 0.0f;
    for (int i = 0; i < AMBE2400_ENC_FFT_SIZE; i++) {
        dc += q->buf[i];
    }
    dc /= (float)AMBE2400_ENC_FFT_SIZE;

    for (int i = 0; i < AMBE2400_ENC_FFT_SIZE; i++) {
        float w = (float)(0.54 - (0.46 * cos((2.0 * M_PI * (double)i) / (double)(AMBE2400_ENC_FFT_SIZE - 1))));
        q->buf[i] -= dc;
        q->windowed[i] = q->buf[i] * w;
    }
}

static void
ambe2400_enc_quantize_pitch(mbe_ambe2400_encoder* enc, struct ambe2400_enc_frame* q) {
    /* Pitch */
    q->strength = 0.0f;
    float lag_f = ambe2400_enc_pitch(enc, q->buf, AMBE2400_ENC_FFT_SIZE, &q->strength);
    float f0 = 1.0f / lag_f;

    q->b[0] = (int)lroundf((log2f(f0) - AMBE2400_ENC_F0_OFFSET) / AMBE2400_ENC_F0_STEP - 0.5f);
    if (q->b[0] < 0) {
        q->b[0] = 0;
    }
    if (q->b[0] > 125) {
        q->b[0] = 125;
    }

    q->L = mbe_clamp_harmonic_count((int)AmbePlusLtable[q->b[0]]);
    q->f0q = exp2f(AMBE2400_ENC_F0_OFFSET + (AMBE2400_ENC_F0_STEP * ((float)q->b[0] + 0.5f)));
}

static void
ambe2400_enc_voicing(mbe_ambe2400_encoder* enc, struct ambe2400_enc_frame* q) {
    /* Voicing: energy-aware. Track the frame's energy against a
     * running max and bias weak, non-periodic frames toward unvoiced. A weak
     * frame the AMDF still judges periodic is kept voiced; weak + noisy
     * frames (consonants, silences) are forced unvoiced so the decoder
     * renders them as soft noise rather than a wrong periodic "splatter".
     * Otherwise force voiced harmonics up to a cutoff that rises with the
     * voicing strength, since the peak/background test under-voices high
     * harmonics of gliding/soft voiced speech. */
    {
        float frame_e = 0.0f;
        for (int i = 0; i < AMBE2400_ENC_FFT_SIZE; i++) {
            frame_e += q->buf[i] * q->buf[i];
        }
        frame_e /= (float)AMBE2400_ENC_FFT_SIZE;

        if (frame_e > enc->max_energy) {
            enc->max_energy = frame_e;
        } else {
            enc->max_energy = (0.99f * enc->max_energy) + (0.01f * frame_e);
        }

        float rel = frame_e / (enc->max_energy + 1e-12f);

        if (rel < 0.03f && q->strength < 0.65f) {
            for (int l = 1; l <= q->L; l++) {
                q->Vl_ana[l] = 0;
            }
        } else {
            int max_jl = -1;
            if (q->strength > 0.60f) {
                max_jl = 7; /* ~3.5 kHz: fully voiced */
            } else if (q->strength > 0.40f) {
                max_jl = 5; /* ~2.5 kHz */
            }
            if (max_jl >= 0) {
                for (int l = 1; l <= q->L; l++) {
                    int jl = (int)((float)l * 16.0f * q->f0q);
                    if (jl <= max_jl) {
                        q->Vl_ana[l] = 1;
                    }
                }
            }
        }
    }
}

static void
ambe2400_enc_quantize_vuv(struct ambe2400_enc_frame* q, const mbe_parms* prev_mp) {
    /* V/UV quantization */
    {
        int best_row = 0;
        float best_dist = 1e30f;

        for (int row = 0; row < 16; row++) {
            float dist = 0.0f;
            for (int l = 1; l <= q->L; l++) {
                int jl = (int)((float)l * 16.0f * q->f0q);
                if (jl > 7) {
                    jl = 7;
                }
                float d = (float)(q->Vl_ana[l] - AmbePlusVuv[row][jl]);
                /* hysteresis: a band voiced last frame that is borderline now
                 * stays voiced (halve the cost of retaining voicing) */
                if (q->Vl_ana[l] == 0 && prev_mp->Vl[l] == 1) {
                    d *= 0.5f;
                }
                dist += fabsf(d);
            }
            if (dist < best_dist) {
                best_dist = dist;
                best_row = row;
            }
        }
        q->b[1] = best_row;
    }
}

static void
ambe2400_enc_quantize_gain(struct ambe2400_enc_frame* q, const mbe_parms* prev_mp) {
    q->mean_a = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        q->a[l] = log2f((q->mag[l] * AMBE2400_ENC_MAG_SCALE) + 1e-12f);
        q->mean_a += q->a[l];
    }
    q->mean_a /= (float)q->L;

    /* Gain quantization */
    {
        float gamma_raw = q->mean_a + (0.5f * log2f((float)q->L));
        float target = gamma_raw - (0.5f * prev_mp->gamma);
        int best = 0;
        float best_err = 1e30f;

        for (int c = 0; c < 64; c++) {
            float err = fabsf(AmbePlusDg[c] - target);
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        q->b[2] = best;
        q->gamma_q = AmbePlusDg[q->b[2]] + (0.5f * prev_mp->gamma);
    }
}

static void
ambe2400_enc_prediction(struct ambe2400_enc_frame* q, const mbe_parms* prev_mp) {
    int prev_L;
    /* Spectral prediction and residual */
    prev_L = mbe_clamp_harmonic_count(prev_mp->L);
    q->L = mbe_clamp_harmonic_count(q->L);
    float prev_log2Ml[57];
    memcpy(prev_log2Ml, prev_mp->log2Ml, sizeof(prev_log2Ml));
    for (int l = prev_L + 1; l <= q->L; l++) {
        prev_log2Ml[l] = prev_log2Ml[prev_L];
    }
    prev_log2Ml[0] = prev_log2Ml[1];

    q->mean_p = 0.0f;
    for (int l = 1; l <= q->L; l++) {
        float flokl = ambe2400_prediction_position(prev_L, q->L, l);
        int intkl = (int)flokl;
        float deltal = flokl - (float)intkl;
        int upper = intkl + 1;
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }
        q->p[l] = ambe2400_interpolate_prediction(deltal, prev_log2Ml[intkl], prev_log2Ml[upper]);
        q->mean_p += q->p[l];
    }
    q->mean_p /= (float)q->L;

    for (int l = 1; l <= q->L; l++) {
        q->Tl[l] = q->a[l] - q->mean_a - (0.65f * (q->p[l] - q->mean_p));
    }
}

static void
ambe2400_enc_block_dct(struct ambe2400_enc_frame* q) {
    /* Block DCTs */
    q->Ji[1] = AmbePlusLmprbl[q->L][0];
    q->Ji[2] = AmbePlusLmprbl[q->L][1];
    q->Ji[3] = AmbePlusLmprbl[q->L][2];
    q->Ji[4] = AmbePlusLmprbl[q->L][3];

    {
        int l = 1;
        for (int blk = 1; blk <= 4; blk++) {
            int ji = q->Ji[blk];
            for (int k = 1; k <= ji; k++) {
                float sum = 0.0f;
                for (int j = 1; j <= ji; j++) {
                    sum += q->Tl[l + j - 1] * q->cache->idct_cos[ji][j][k];
                }
                /* Decoder IDCT is Tl[j] = sum_k a_k Cik[k] cos(theta_kj),
                 * a_1=1, a_k=2. Exact inverse: Cik[k] = (1/ji) sum_j Tl[j] cos. */
                q->Cik[blk][k] = (1.0f / (float)ji) * sum;
            }
            l += ji;
        }
    }
}

static void
ambe2400_enc_prba_dct(struct ambe2400_enc_frame* q) {
    float Ri[9];
    /* Ri from block DC-ish terms */
    const float sqrt2 = 1.41421356237f;
    Ri[1] = q->Cik[1][1] + (sqrt2 * q->Cik[1][2]);
    Ri[2] = q->Cik[1][1] - (sqrt2 * q->Cik[1][2]);
    Ri[3] = q->Cik[2][1] + (sqrt2 * q->Cik[2][2]);
    Ri[4] = q->Cik[2][1] - (sqrt2 * q->Cik[2][2]);
    Ri[5] = q->Cik[3][1] + (sqrt2 * q->Cik[3][2]);
    Ri[6] = q->Cik[3][1] - (sqrt2 * q->Cik[3][2]);
    Ri[7] = q->Cik[4][1] + (sqrt2 * q->Cik[4][2]);
    Ri[8] = q->Cik[4][1] - (sqrt2 * q->Cik[4][2]);

    /* Gm: 8-point DCT of Ri (Gm[1] is the discarded DC term) */
    for (int m = 2; m <= 8; m++) {
        float sum = 0.0f;
        for (int i = 1; i <= 8; i++) {
            sum += Ri[i] * q->cache->ri_cos[m][i];
        }
        q->Gm[m] = sum / 8.0f;
    }
}

static void
ambe2400_enc_quantize_prba(struct ambe2400_enc_frame* q) {
    /* PRBA24 (b3): Gm[2..4] */
    {
        int best = 0;
        float best_err = 1e30f;
        for (int c = 0; c < 512; c++) {
            float err = 0.0f;
            for (int m = 0; m < 3; m++) {
                float d = q->Gm[2 + m] - AmbePlusPRBA24[c][m];
                err += d * d;
            }
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        q->b[3] = best;
    }

    /* PRBA58 (b4): Gm[5..8] */
    {
        int best = 0;
        float best_err = 1e30f;
        for (int c = 0; c < 128; c++) {
            float err = 0.0f;
            for (int m = 0; m < 4; m++) {
                float d = q->Gm[5 + m] - AmbePlusPRBA58[c][m];
                err += d * d;
            }
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        q->b[4] = best;
    }
}

static const float (*const ambe2400_hoc_tables[4])[4] = {AmbePlusHOCb5, AmbePlusHOCb6, AmbePlusHOCb7, AmbePlusHOCb8};

static void
ambe2400_enc_quantize_hoc(struct ambe2400_enc_frame* q) {
    /* HOC blocks (b5..b8): Cik[blk][3..min(Ji,6)] */
    {
        int codes[4];

        for (int blk = 0; blk < 4; blk++) {
            int ji = q->Ji[blk + 1];
            int kmax = (ji < 6) ? ji : 6;
            int best = 0;
            float best_err = 1e30f;

            /* Block 4 carries only bits 3..1, so only even rows are reachable. */
            for (int c = 0; c < 16; c += (blk == 3) ? 2 : 1) {
                float err = 0.0f;
                for (int k = 3; k <= kmax; k++) {
                    float d = q->Cik[blk + 1][k] - ambe2400_hoc_tables[blk][c][k - 3];
                    err += d * d;
                }
                if (err < best_err) {
                    best_err = err;
                    best = c;
                }
            }
            codes[blk] = best;
        }
        q->b[5] = codes[0];
        q->b[6] = codes[1];
        q->b[7] = codes[2];
        q->b[8] = codes[3];
    }
}

static void
ambe2400_enc_reconstruct_coefficients(struct ambe2400_enc_frame* q) {
    float Ri_q[9];
    mbe_ambe2400_reconstruct_prba(q->b[3], q->b[4], Ri_q);
    mbe_ambe2400_reconstruct_cik(Ri_q, q->b + 5, q->Ji, q->Cik_q);
}

static void
ambe2400_enc_pack(const struct ambe2400_enc_frame* q, char ambe_d[49]) {
    /* Pack the 49-bit AMBE word */
    memset(ambe_d, 0, 49);
    for (int i = 0; i < 6; i++) {
        ambe_d[i] = (char)((q->b[0] >> (6 - i)) & 1);
    }
    ambe_d[48] = (char)(q->b[0] & 1);

    ambe_d[38] = (char)((q->b[1] >> 3) & 1);
    ambe_d[39] = (char)((q->b[1] >> 2) & 1);
    ambe_d[40] = (char)((q->b[1] >> 1) & 1);
    ambe_d[41] = (char)(q->b[1] & 1);

    ambe_d[6] = (char)((q->b[2] >> 5) & 1);
    ambe_d[7] = (char)((q->b[2] >> 4) & 1);
    ambe_d[8] = (char)((q->b[2] >> 3) & 1);
    ambe_d[9] = (char)((q->b[2] >> 2) & 1);
    ambe_d[42] = (char)((q->b[2] >> 1) & 1);
    ambe_d[43] = (char)(q->b[2] & 1);

    ambe_d[10] = (char)((q->b[3] >> 8) & 1);
    ambe_d[11] = (char)((q->b[3] >> 7) & 1);
    ambe_d[12] = (char)((q->b[3] >> 6) & 1);
    ambe_d[13] = (char)((q->b[3] >> 5) & 1);
    ambe_d[14] = (char)((q->b[3] >> 4) & 1);
    ambe_d[15] = (char)((q->b[3] >> 3) & 1);
    ambe_d[16] = (char)((q->b[3] >> 2) & 1);
    ambe_d[44] = (char)((q->b[3] >> 1) & 1);
    ambe_d[45] = (char)(q->b[3] & 1);

    ambe_d[17] = (char)((q->b[4] >> 6) & 1);
    ambe_d[18] = (char)((q->b[4] >> 5) & 1);
    ambe_d[19] = (char)((q->b[4] >> 4) & 1);
    ambe_d[20] = (char)((q->b[4] >> 3) & 1);
    ambe_d[21] = (char)((q->b[4] >> 2) & 1);
    ambe_d[46] = (char)((q->b[4] >> 1) & 1);
    ambe_d[47] = (char)(q->b[4] & 1);

    ambe_d[22] = (char)((q->b[5] >> 3) & 1);
    ambe_d[23] = (char)((q->b[5] >> 2) & 1);
    ambe_d[25] = (char)((q->b[5] >> 1) & 1);
    ambe_d[26] = (char)(q->b[5] & 1);

    ambe_d[27] = (char)((q->b[6] >> 3) & 1);
    ambe_d[28] = (char)((q->b[6] >> 2) & 1);
    ambe_d[29] = (char)((q->b[6] >> 1) & 1);
    ambe_d[30] = (char)(q->b[6] & 1);

    ambe_d[31] = (char)((q->b[7] >> 3) & 1);
    ambe_d[32] = (char)((q->b[7] >> 2) & 1);
    ambe_d[33] = (char)((q->b[7] >> 1) & 1);
    ambe_d[34] = (char)(q->b[7] & 1);

    ambe_d[35] = (char)((q->b[8] >> 3) & 1);
    ambe_d[36] = (char)((q->b[8] >> 2) & 1);
    ambe_d[37] = (char)((q->b[8] >> 1) & 1);

    /* ambe_d[24] is the spare bit; filled by the FEC layer below. */
}

static void
ambe2400_enc_fill_parms(const struct ambe2400_enc_frame* q, mbe_parms* cur_mp) {
    /* Write quantized parameters to cur_mp (decoder-equivalent state) */
    cur_mp->w0 = q->f0q * (float)2 * M_PI;
    cur_mp->L = q->L;
    cur_mp->K = 0;
    cur_mp->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
    for (int l = 1; l <= q->L; l++) {
        int jl = (int)((float)l * 16.0f * q->f0q);
        if (jl > 7) {
            jl = 7;
        }
        cur_mp->Vl[l] = AmbePlusVuv[q->b[1]][jl];
        if (cur_mp->Vl[l] == 1) {
            cur_mp->K++;
        }
    }
    cur_mp->gamma = q->gamma_q;
}

/* Quantize a voice frame, reconstruct its predictor state, and pack 49 bits. */
static int
ambe2400_encode_voice(mbe_ambe2400_encoder* enc, const float* pcm, char ambe_d[49], mbe_parms* cur_mp,
                      const mbe_parms* prev_mp) {
    struct ambe2400_enc_frame q = {0};
    q.cache = mbe_ambe2400_get_dct_cache();
    ambe2400_enc_window(enc, &q, pcm);
    ambe2400_enc_quantize_pitch(enc, &q);
    int status = ambe2400_enc_spectrum(enc, q.windowed, q.f0q, q.L, q.mag, q.Vl_ana);
    if (status < 0) {
        return status;
    }
    ambe2400_enc_voicing(enc, &q);
    ambe2400_enc_quantize_vuv(&q, prev_mp);
    ambe2400_enc_quantize_gain(&q, prev_mp);
    ambe2400_enc_prediction(&q, prev_mp);
    ambe2400_enc_block_dct(&q);
    ambe2400_enc_prba_dct(&q);
    ambe2400_enc_quantize_prba(&q);
    ambe2400_enc_quantize_hoc(&q);
    ambe2400_enc_reconstruct_coefficients(&q);
    mbe_ambe2400_inverse_dct_tl(q.Cik_q, q.Ji, q.Tl_q);
    ambe2400_enc_pack(&q, ambe_d);
    ambe2400_enc_fill_parms(&q, cur_mp);
    /* Share the decoder update without changing the caller's predictor. */
    mbe_parms prediction_prev = *prev_mp;
    mbe_ambe2400_update_spectral_amplitudes(cur_mp, &prediction_prev, q.Tl_q, 0.2046f / sqrtf(cur_mp->w0));
    return 0;
}

/*
 * Silence frame: b0 = 127, tone index = 128 (the standard AMBE silence
 * pattern). The decoder resets its parameter state on this frame, so the
 * encoder resets its prediction state as well.
 */
static void
ambe2400_encode_silence(char ambe_d[49], mbe_parms* cur_mp) {
    memset(ambe_d, 0, 49);
    for (int i = 0; i < 6; i++) {
        ambe_d[i] = 1;
    }
    ambe_d[48] = 1;

    /* tone index 128: t7t6t5 = 100, low bits 0 */
    ambe_d[6] = 0;
    ambe_d[7] = 0;
    ambe_d[8] = 0;
    ambe_d[9] = 0;
    ambe_d[42] = 0;
    ambe_d[43] = 0;
    ambe_d[10] = 0;
    ambe_d[11] = 0;

    cur_mp->w0 = (float)((M_PI / 32.0) * (2.0 * M_PI));
    cur_mp->L = 15;
    cur_mp->K = 0;
    cur_mp->gamma = 0.0f;
    cur_mp->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
    for (int l = 0; l <= 56; l++) {
        cur_mp->Ml[l] = 1.0f;
        cur_mp->Vl[l] = 0;
        cur_mp->log2Ml[l] = 0.0f;
    }
}

static bool
ambe2400_enc_silence_gate(mbe_ambe2400_encoder* enc, float rms) {
    /* Silence gate with hang-over: a frame is only "silence" after a few
     * consecutive quiet frames (squelch close delay); speech opens the gate
     * immediately. This avoids the squelch-like flapping between silence and
     * speech on boundary frames. Noise within speech is left to caller-side
     * preprocessing; the gate only decides silence frames. */
    bool is_silence = false;
    if (rms < AMBE2400_ENC_SILENCE_RMS) {
        if (enc->silence_run < 5) {
            enc->silence_run++;
        }
        is_silence = (enc->silence_run >= 5);
    } else {
        enc->silence_run = 0;
    }

    return is_silence;
}

static void
ambe2400_enc_apply_agc(mbe_ambe2400_encoder* enc, const float* samples, float rms, bool is_silence, float* agc_buf) {
    /* Input AGC: track level only on speech frames, normalize to the
     * nominal target. Noise frames would drag the level estimate down and
     * over-drive the speech. */
    if (!is_silence && rms > 1e-6f) {
        enc->agc_rms = (AMBE2400_ENC_AGC_ALPHA * enc->agc_rms) + ((1.0f - AMBE2400_ENC_AGC_ALPHA) * rms);
    }
    enc->agc_gain = AMBE2400_ENC_AGC_TARGET / (enc->agc_rms + 1e-9f);
    if (enc->agc_gain < AMBE2400_ENC_AGC_MIN) {
        enc->agc_gain = AMBE2400_ENC_AGC_MIN;
    }
    if (enc->agc_gain > AMBE2400_ENC_AGC_MAX) {
        enc->agc_gain = AMBE2400_ENC_AGC_MAX;
    }
    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        agc_buf[i] = samples[i] * enc->agc_gain;
    }
}

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of PCM into AMBE 2400
 *        parameter bits.
 *
 * @param samples Input PCM floats (160), nominal range [-1, 1].
 * @param ambe_d  Output parameter bits (49).
 * @param cur_mp  Output: quantized (decoder-equivalent) parameters.
 * @param prev_mp Input: previous frame state (see mbe_initMbeParms()).
 * @return 0 for a voice frame, 1 for a silence frame, negative on error.
 */
int
mbe_encodeAmbe2400Parms(mbe_ambe2400_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                        const mbe_parms* prev_mp) {
    float agc_buf[AMBE2400_ENC_SAMPLES];
    float rms = 0.0f;

    if (enc == NULL || samples == NULL || ambe_d == NULL || cur_mp == NULL || prev_mp == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        rms += samples[i] * samples[i];
    }
    rms = sqrtf(rms / (float)AMBE2400_ENC_SAMPLES);

    bool is_silence = ambe2400_enc_silence_gate(enc, rms);

    ambe2400_enc_apply_agc(enc, samples, rms, is_silence, agc_buf);

    if (is_silence) {
        ambe2400_encode_silence(ambe_d, cur_mp);
    } else {
        int ret = ambe2400_encode_voice(enc, agc_buf, ambe_d, cur_mp, prev_mp);
        if (ret < 0) {
            return ret;
        }
    }

    memcpy(enc->history, agc_buf, (size_t)AMBE2400_ENC_SAMPLES * sizeof(float));
    enc->hist_gain = enc->agc_gain;

    return (int)is_silence;
}

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of 16-bit PCM into AMBE 2400
 *        parameter bits.
 *
 * @see mbe_encodeAmbe2400Parms for details.
 */
int
mbe_encodeAmbe2400ParmsShort(mbe_ambe2400_encoder* enc, const short* samples, char ambe_d[49], mbe_parms* cur_mp,
                             const mbe_parms* prev_mp) {
    float float_buf[AMBE2400_ENC_SAMPLES];

    if (samples == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        float_buf[i] = (float)samples[i] / 32768.0f;
    }
    return mbe_encodeAmbe2400Parms(enc, float_buf, ambe_d, cur_mp, prev_mp);
}

static void
ambe2400_enc_c0(const char* ambe_d, char ambe_fr[4][24]) {
    char cw[23], data12[12];
    /* C0: 12 data bits (ambe_d[0..11]) -> (23,12) Golay + even parity */
    for (int i = 0; i < 12; i++) {
        data12[i] = ambe_d[i];
    }
    mbe_golay2312_encode(data12, cw);
    for (int j = 0; j < 23; j++) {
        ambe_fr[0][j + 1] = cw[j];
    }
    int ones = 0;
    for (int j = 0; j < 24; j++) {
        ones += (ambe_fr[0][j] & 1);
    }
    ambe_fr[0][0] = (char)(ones & 1);
}

static void
ambe2400_enc_pr(const char ambe_fr[4][24], unsigned short pr[25]) {
    unsigned short foo = 0;
    /* PR scramble sequence from the C0 data word */
    for (int i = 23; i >= 12; i--) {
        foo <<= 1;
        foo |= (unsigned short)(ambe_fr[0][i] & 1);
    }
    pr[0] = (unsigned short)(16 * foo);
    for (int i = 1; i < 25; i++) {
        pr[i] = (unsigned short)((173 * pr[i - 1]) + 13849 - (65536 * (((173 * pr[i - 1]) + 13849) / 65536)));
    }
}

static void
ambe2400_enc_c1(const char* ambe_d, char ambe_fr[4][24], const unsigned short pr[25]) {
    char cw[23], data12[12];
    /* C1: 12 data bits (ambe_d[12..23]) -> Golay, XOR p bits 23..1 */
    for (int i = 0; i < 12; i++) {
        data12[i] = ambe_d[12 + i];
    }
    mbe_golay2312_encode(data12, cw);
    for (int j = 0; j < 23; j++) {
        /* fr[1][j] = b word bit (j+1), scrambled with p bit (j+1) = pr[24-(j+1)] >> 15 */
        int prbit = pr[23 - j] / 32768;
        ambe_fr[1][j] = (char)((cw[j] & 1) ^ prbit);
    }

    /* Spare bit (ambe_d[24] = fr[2][10] on air): b word bit 0 = parity ^ p bit 0 */
    int ones = 0;
    for (int j = 0; j < 23; j++) {
        ones += (cw[j] & 1);
    }
    ambe_fr[2][10] = (char)(((ones & 1) ^ (pr[24] / 32768)) & 1);
}

static void
ambe2400_enc_raw(const char* ambe_d, char ambe_fr[4][24]) {
    /* C2/C3: raw bits (decoder reads them MSB-first) */
    for (int j = 0; j < 11; j++) {
        ambe_fr[2][j] = ambe_d[24 + (10 - j)];
    }
    for (int j = 0; j < 14; j++) {
        ambe_fr[3][j] = ambe_d[35 + (13 - j)];
    }
}

/**
 * @brief Encode 49 AMBE 2400 parameter bits into a 72-bit D-STAR DV
 *        data frame (FEC + interleave), in the same plane layout the
 *        decoder (mbe_decodeAmbe3600x2400Frame) consumes.
 *
 * ambe_d[24] is the spare bit; on output it carries the scrambled even
 * parity of the second Golay codeword. All other input bits round-trip
 * exactly.
 *
 * @param ambe_d  Input parameter bits (49).
 * @param ambe_fr Output frame as 4x24 bitplanes.
 * @return 0 on success.
 */
int
mbe_encodeAmbe3600x2400Frame(const char ambe_d[49], char ambe_fr[4][24]) {
    unsigned short pr[25];

    if (ambe_d == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 49; i++) {
        if (ambe_d[i] != 0 && ambe_d[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }

    memset(ambe_fr, 0, 4 * sizeof(ambe_fr[0]));

    ambe2400_enc_c0(ambe_d, ambe_fr);
    ambe2400_enc_pr((const char (*)[24])ambe_fr, pr);
    ambe2400_enc_raw(ambe_d, ambe_fr);
    ambe2400_enc_c1(ambe_d, ambe_fr, pr);

    return 0;
}
