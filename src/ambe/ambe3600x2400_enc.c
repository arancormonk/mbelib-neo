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
 *    reconstructed frame is bit-compatible with this library's decoder
 *    and (best effort) with DVSI AMBE-3000 based receivers.
 *  - The prediction state (log2Ml, gamma) is advanced with the
 *    QUANTIZED values, mirroring the decoder state so the encoder and
 *    decoder never drift apart.
 *
 * Tone (DTMF) encoding is not supported; silent input produces the
 * standard AMBE silence frame (b0 == 127, tone index 128).
 */

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "ambe3600x2400_const.h"
#include "ecc_const.h"
#include "mbe_compiler.h"
#include "mbe_ecc.h"
#include "mbe_validation.h"
#include "mbelib-neo/mbelib.h"
#include "pffft.h"

#define AMBE2400_ENC_FFT_SIZE    256
#define AMBE2400_ENC_SAMPLES     160
#define AMBE2400_ENC_SILENCE_RMS 0.0015f

/* f0 = exp2(-4.311767578125 - 2.1336e-2 * (b0 + 0.5)) */
#define AMBE2400_ENC_F0_OFFSET   (-4.311767578125f)
#define AMBE2400_ENC_F0_STEP     (-0.021336f)

/*
 * Calibration: the AMBE 2400 gain field (AmbePlusDg) is non-negative and
 * spans ~5.3 dB, so the encoder must keep the input at a fixed nominal
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
 * Thread-local analysis state: one encode session per thread.
 *
 * pcm_history holds the previous 160 samples so the encoder can analyse
 * a 256 sample (32 ms) window centred on the current frame without
 * blocking for look-ahead.
 */
static MBE_THREAD_LOCAL float ambe2400_enc_history[AMBE2400_ENC_SAMPLES];
static MBE_THREAD_LOCAL float ambe2400_enc_agc_gain = 1.0f;
static MBE_THREAD_LOCAL float ambe2400_enc_agc_rms = AMBE2400_ENC_AGC_TARGET;
static MBE_THREAD_LOCAL float ambe2400_enc_hist_gain = 1.0f;
static MBE_THREAD_LOCAL int ambe2400_enc_silence_run = 0;
static MBE_THREAD_LOCAL float ambe2400_enc_prev_lag = 0.0f;
static MBE_THREAD_LOCAL float ambe2400_enc_max_energy = 1e-6f;

static MBE_THREAD_LOCAL PFFFT_Setup* ambe2400_enc_fft;
static MBE_THREAD_LOCAL float* ambe2400_enc_work;

struct ambe2400_dct_cache {
    int inited;
    float blk_cos[18][18][18]; /* [ji][j][k], ji=1..17, j,k=1..ji */
    float prba_cos[9][9];      /* [m][i] 8-point DCT basis */
};

static MBE_THREAD_LOCAL struct ambe2400_dct_cache ambe2400_enc_cache;

static struct ambe2400_dct_cache*
ambe2400_enc_get_dct_cache(void) {
    struct ambe2400_dct_cache* cache = &ambe2400_enc_cache;
    if (cache->inited) {
        return cache;
    }

    /* blk_cos[ji][j][k] = cos(pi*(k-1)*(j-0.5)/ji) */
    for (int ji = 1; ji <= 17; ji++) {
        for (int j = 1; j <= ji; j++) {
            for (int k = 1; k <= ji; k++) {
                cache->blk_cos[ji][j][k] = cosf((M_PI * (float)(k - 1) * ((float)j - 0.5f)) / (float)ji);
            }
        }
    }

    /* prba_cos[m][i] = cos(pi*(m-1)*(i-0.5)/8), m,i = 1..8 */
    for (int m = 1; m <= 8; m++) {
        for (int i = 1; i <= 8; i++) {
            cache->prba_cos[m][i] = cosf((M_PI * (float)(m - 1) * ((float)i - 0.5f)) / 8.0f);
        }
    }

    cache->inited = 1;
    return cache;
}

static int
ambe2400_enc_fft_init(void) {
    if (ambe2400_enc_fft == NULL) {
        ambe2400_enc_fft = pffft_new_setup(AMBE2400_ENC_FFT_SIZE, PFFFT_REAL);
        if (ambe2400_enc_fft == NULL) {
            return -1;
        }
        ambe2400_enc_work = (float*)pffft_aligned_malloc((size_t)AMBE2400_ENC_FFT_SIZE * sizeof(float));
        if (ambe2400_enc_work == NULL) {
            pffft_destroy_setup(ambe2400_enc_fft);
            ambe2400_enc_fft = NULL;
            return -1;
        }
    }
    return 0;
}

/*
 * Compute log2 spectral magnitudes and raw per-harmonic voicing from a
 * pitch estimate over the windowed analysis buffer.
 */
static void
ambe2400_enc_spectrum(const float* windowed, float f0q, int L, float mag[57], int vl_ana[57]) {
    /* pffft requires 16-byte aligned buffers; stack arrays may not be. */
    static MBE_THREAD_LOCAL float* fft_in;
    static MBE_THREAD_LOCAL float* fft_out;
    const float f0_bin = f0q * (float)AMBE2400_ENC_FFT_SIZE;
    float max_band_energy = 0.0f;

    if (fft_in == NULL) {
        fft_in = (float*)pffft_aligned_malloc((size_t)AMBE2400_ENC_FFT_SIZE * sizeof(float));
        fft_out = (float*)pffft_aligned_malloc((size_t)AMBE2400_ENC_FFT_SIZE * sizeof(float));
        if (fft_in == NULL || fft_out == NULL) {
            return;
        }
    }

    memcpy(fft_in, windowed, (size_t)AMBE2400_ENC_FFT_SIZE * sizeof(float));
    pffft_transform_ordered(ambe2400_enc_fft, fft_in, fft_out, ambe2400_enc_work, PFFFT_FORWARD);

    /* Ordered real FFT layout: [DC, re1, im1, re2, im2, ..., reN/2] */
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

        for (int b = lo_bin; b <= hi_bin && b <= (AMBE2400_ENC_FFT_SIZE / 2) - 1; b++) {
            float re = fft_out[2 * b];
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
ambe2400_enc_pitch(const float* buf, int n, float* strength) {
    const int min_lag = 20;
    const int max_lag = 127;
    float amdf[128];
    float global_min = 1e30f;
    int global_min_lag = 0;

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
            global_min_lag = lag;
        }
    }

    const float tol = global_min * 1.4f;
    int have_prev = (ambe2400_enc_prev_lag >= (float)min_lag);

    float best_lag_f = (float)global_min_lag;
    float best_score = 1e30f;

    for (int lag = min_lag; lag <= max_lag; lag++) {
        if ((lag > min_lag && amdf[lag] >= amdf[lag - 1]) || (lag < max_lag && amdf[lag] > amdf[lag + 1])) {
            continue;
        }
        if (amdf[lag] > tol) {
            continue;
        }

        /* Parabolic refinement */
        float denom = 0.0f;
        if (lag > min_lag && lag < max_lag) {
            denom = amdf[lag - 1] + amdf[lag + 1] - (2.0f * amdf[lag]);
        }
        float lag_f = (float)lag;
        if (fabsf(denom) > 1e-12f) {
            float shift = 0.5f * (amdf[lag - 1] - amdf[lag + 1]) / denom;
            if (shift > -1.0f && shift < 1.0f) {
                lag_f += shift;
            }
        }

        float score;
        if (have_prev) {
            /* AMDF quality + continuity (octave-aware). */
            score = (amdf[lag] / (global_min + 1e-12f)) + (0.8f * fabsf(log2f(lag_f / ambe2400_enc_prev_lag)));
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

    /* Final guard: limit frame-to-frame jumps to 25%. */
    if (have_prev) {
        float max_jump = 0.25f * ambe2400_enc_prev_lag;
        float diff = best_lag_f - ambe2400_enc_prev_lag;
        if (fabsf(diff) > max_jump) {
            best_lag_f = ambe2400_enc_prev_lag + (diff > 0.0f ? max_jump : -max_jump);
        }
    }

    ambe2400_enc_prev_lag = best_lag_f;

    if (strength != NULL) {
        /* Voicing strength: normalized autocorrelation at the chosen lag.
         * ~1.0 for clean periodic speech, ~0.0 for noise. */
        int lag = (int)(best_lag_f + 0.5f);
        if (lag < min_lag) {
            lag = min_lag;
        }
        if (lag > max_lag) {
            lag = max_lag;
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
        *strength = r;
        if (*strength < 0.0f) {
            *strength = 0.0f;
        }
        if (*strength > 1.0f) {
            *strength = 1.0f;
        }
    }
    return best_lag_f;
}

/*
 * Quantize all parameters, pack the 49-bit AMBE word and advance the
 * prediction state. Returns 0 for a voice frame, 1 for a silence frame.
 */
static int
ambe2400_encode_voice(const float* pcm, char ambe_d[49], mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    const struct ambe2400_dct_cache* cache = ambe2400_enc_get_dct_cache();
    float buf[AMBE2400_ENC_FFT_SIZE];
    float windowed[AMBE2400_ENC_FFT_SIZE];
    float mag[57] = {0.0f};
    float a[57] = {0.0f};
    float p[57] = {0.0f};
    float Tl[57] = {0.0f};
    float Tl_q[57] = {0.0f};
    float Cik[5][18];
    float Cik_q[5][18];
    float Ri[9];
    float Ri_q[9];
    float Gm[9];
    float log2Ml_q[57] = {0.0f};
    int Vl_ana[57] = {0};
    int Ji[5];
    int L;
    int b0;
    int b1;
    int b2;
    int b3;
    int b4;
    int b5;
    int b6;
    int b7;
    int b8;
    float f0q;
    float gamma_q;
    int prev_L;

    memset(Cik, 0, sizeof(Cik));
    memset(Cik_q, 0, sizeof(Cik_q));

    /* Build the centred analysis window: last 128 of history + first 128 of pcm.
     * The AGC gain is applied per frame; ramp it smoothly across the window so
     * the frame boundary does not introduce a step discontinuity (which smears
     * the spectrum and destroys the voicing decision). */
    for (int i = 0; i < 128; i++) {
        buf[i] = ambe2400_enc_history[AMBE2400_ENC_SAMPLES - 128 + i] / ambe2400_enc_hist_gain;
    }
    for (int i = 0; i < 128; i++) {
        buf[128 + i] = pcm[i] / ambe2400_enc_agc_gain;
    }
    for (int i = 0; i < AMBE2400_ENC_FFT_SIZE; i++) {
        float t = (float)i / (float)(AMBE2400_ENC_FFT_SIZE - 1);
        float g = ambe2400_enc_hist_gain + (ambe2400_enc_agc_gain - ambe2400_enc_hist_gain) * t;
        buf[i] *= g;
    }

    /* DC removal */
    float dc = 0.0f;
    for (int i = 0; i < AMBE2400_ENC_FFT_SIZE; i++) {
        dc += buf[i];
    }
    dc /= (float)AMBE2400_ENC_FFT_SIZE;

    for (int i = 0; i < AMBE2400_ENC_FFT_SIZE; i++) {
        float w = (float)(0.54 - (0.46 * cos((2.0 * M_PI * (double)i) / (double)(AMBE2400_ENC_FFT_SIZE - 1))));
        buf[i] -= dc;
        windowed[i] = buf[i] * w;
    }
    /* Pitch */
    float strength = 0.0f;
    float lag_f = ambe2400_enc_pitch(buf, AMBE2400_ENC_FFT_SIZE, &strength);
    float f0 = 1.0f / lag_f;

    b0 = (int)lroundf((log2f(f0) - AMBE2400_ENC_F0_OFFSET) / AMBE2400_ENC_F0_STEP - 0.5f);
    if (b0 < 0) {
        b0 = 0;
    }
    if (b0 > 125) {
        b0 = 125;
    }

    L = mbe_clamp_harmonic_count((int)AmbePlusLtable[b0]);
    f0q = exp2f(AMBE2400_ENC_F0_OFFSET + (AMBE2400_ENC_F0_STEP * ((float)b0 + 0.5f)));

    /* Spectral analysis */
    ambe2400_enc_spectrum(windowed, f0q, L, mag, Vl_ana);

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
            frame_e += buf[i] * buf[i];
        }
        frame_e /= (float)AMBE2400_ENC_FFT_SIZE;

        if (frame_e > ambe2400_enc_max_energy) {
            ambe2400_enc_max_energy = frame_e;
        } else {
            ambe2400_enc_max_energy = (0.99f * ambe2400_enc_max_energy) + (0.01f * frame_e);
        }

        float rel = frame_e / (ambe2400_enc_max_energy + 1e-12f);

        if (rel < 0.03f && strength < 0.65f) {
            for (int l = 1; l <= L; l++) {
                Vl_ana[l] = 0;
            }
        } else {
            int max_jl = -1;
            if (strength > 0.60f) {
                max_jl = 7; /* ~3.5 kHz: fully voiced */
            } else if (strength > 0.40f) {
                max_jl = 5; /* ~2.5 kHz */
            }
            if (max_jl >= 0) {
                for (int l = 1; l <= L; l++) {
                    int jl = (int)((float)l * 16.0f * f0q);
                    if (jl <= max_jl) {
                        Vl_ana[l] = 1;
                    }
                }
            }
        }
    }

    float mean_a = 0.0f;
    for (int l = 1; l <= L; l++) {
        a[l] = log2f((mag[l] * AMBE2400_ENC_MAG_SCALE) + 1e-12f);
        mean_a += a[l];
    }
    mean_a /= (float)L;

    /* V/UV quantization */
    {
        int best_row = 0;
        float best_dist = 1e30f;

        for (int row = 0; row < 16; row++) {
            float dist = 0.0f;
            for (int l = 1; l <= L; l++) {
                int jl = (int)((float)l * 16.0f * f0q);
                if (jl > 7) {
                    jl = 7;
                }
                float d = (float)(Vl_ana[l] - AmbePlusVuv[row][jl]);
                /* hysteresis: a band voiced last frame that is borderline now
                 * stays voiced (halve the cost of retaining voicing) */
                if (Vl_ana[l] == 0 && prev_mp->Vl[l] == 1) {
                    d *= 0.5f;
                }
                dist += fabsf(d);
            }
            if (dist < best_dist) {
                best_dist = dist;
                best_row = row;
            }
        }
        b1 = best_row;
    }

    /* Gain quantization */
    {
        float gamma_raw = mean_a + (0.5f * log2f((float)L));
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
        b2 = best;
        gamma_q = AmbePlusDg[b2] + (0.5f * prev_mp->gamma);
    }

    /* Spectral prediction and residual */
    prev_L = mbe_clamp_harmonic_count(prev_mp->L);
    float prev_log2Ml[57];
    memcpy(prev_log2Ml, prev_mp->log2Ml, sizeof(prev_log2Ml));
    for (int l = prev_L + 1; l <= L; l++) {
        prev_log2Ml[l] = prev_log2Ml[prev_L];
    }
    prev_log2Ml[0] = prev_log2Ml[1];

    float mean_p = 0.0f;
    for (int l = 1; l <= L; l++) {
        float flokl = ((float)prev_L / (float)L) * (float)l;
        int intkl = (int)flokl;
        float deltal = flokl - (float)intkl;
        int upper = intkl + 1;
        float v_lo;
        float v_hi;

        if (intkl > MBE_MAX_HARMONIC_BANDS) {
            intkl = MBE_MAX_HARMONIC_BANDS;
        }
        if (upper > MBE_MAX_HARMONIC_BANDS) {
            upper = MBE_MAX_HARMONIC_BANDS;
        }

        v_lo = prev_log2Ml[intkl];
        v_hi = prev_log2Ml[upper];

        p[l] = ((1.0f - deltal) * v_lo) + (deltal * v_hi);
        mean_p += p[l];
    }
    mean_p /= (float)L;

    for (int l = 1; l <= L; l++) {
        Tl[l] = a[l] - mean_a - (0.65f * (p[l] - mean_p));
    }

    /* Block DCTs */
    Ji[1] = AmbePlusLmprbl[L][0];
    Ji[2] = AmbePlusLmprbl[L][1];
    Ji[3] = AmbePlusLmprbl[L][2];
    Ji[4] = AmbePlusLmprbl[L][3];

    {
        int l = 1;
        for (int blk = 1; blk <= 4; blk++) {
            int ji = Ji[blk];
            for (int k = 1; k <= ji; k++) {
                float sum = 0.0f;
                for (int j = 1; j <= ji; j++) {
                    sum += Tl[l + j - 1] * cache->blk_cos[ji][j][k];
                }
                /* Decoder IDCT is Tl[j] = sum_k a_k Cik[k] cos(theta_kj),
                 * a_1=1, a_k=2. Exact inverse: Cik[k] = (1/ji) sum_j Tl[j] cos. */
                Cik[blk][k] = (1.0f / (float)ji) * sum;
            }
            l += ji;
        }
    }

    /* Ri from block DC-ish terms */
    const float sqrt2 = 1.41421356237f;
    Ri[1] = Cik[1][1] + (sqrt2 * Cik[1][2]);
    Ri[2] = Cik[1][1] - (sqrt2 * Cik[1][2]);
    Ri[3] = Cik[2][1] + (sqrt2 * Cik[2][2]);
    Ri[4] = Cik[2][1] - (sqrt2 * Cik[2][2]);
    Ri[5] = Cik[3][1] + (sqrt2 * Cik[3][2]);
    Ri[6] = Cik[3][1] - (sqrt2 * Cik[3][2]);
    Ri[7] = Cik[4][1] + (sqrt2 * Cik[4][2]);
    Ri[8] = Cik[4][1] - (sqrt2 * Cik[4][2]);

    /* Gm: 8-point DCT of Ri (Gm[1] is the discarded DC term) */
    for (int m = 2; m <= 8; m++) {
        float sum = 0.0f;
        for (int i = 1; i <= 8; i++) {
            sum += Ri[i] * cache->prba_cos[m][i];
        }
        Gm[m] = sum / 8.0f;
    }

    /* PRBA24 (b3): Gm[2..4] */
    {
        int best = 0;
        float best_err = 1e30f;
        for (int c = 0; c < 512; c++) {
            float err = 0.0f;
            for (int m = 0; m < 3; m++) {
                float d = Gm[2 + m] - AmbePlusPRBA24[c][m];
                err += d * d;
            }
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        b3 = best;
    }

    /* PRBA58 (b4): Gm[5..8] */
    {
        int best = 0;
        float best_err = 1e30f;
        for (int c = 0; c < 128; c++) {
            float err = 0.0f;
            for (int m = 0; m < 4; m++) {
                float d = Gm[5 + m] - AmbePlusPRBA58[c][m];
                err += d * d;
            }
            if (err < best_err) {
                best_err = err;
                best = c;
            }
        }
        b4 = best;
    }

    /* HOC blocks (b5..b8): Cik[blk][3..min(Ji,6)] */
    {
        int codes[4];
        const float (*tables[4])[4] = {AmbePlusHOCb5, AmbePlusHOCb6, AmbePlusHOCb7, AmbePlusHOCb8};

        for (int blk = 0; blk < 4; blk++) {
            int ji = Ji[blk + 1];
            int kmax = (ji < 6) ? ji : 6;
            int best = 0;
            float best_err = 1e30f;

            /* Block 4 carries only bits 3..1, so only even rows are reachable. */
            for (int c = 0; c < 16; c += (blk == 3) ? 2 : 1) {
                float err = 0.0f;
                for (int k = 3; k <= kmax; k++) {
                    float d = Cik[blk + 1][k] - tables[blk][c][k - 3];
                    err += d * d;
                }
                if (err < best_err) {
                    best_err = err;
                    best = c;
                }
            }
            codes[blk] = best;
        }
        b5 = codes[0];
        b6 = codes[1];
        b7 = codes[2];
        b8 = codes[3];
    }

    /*
     * Reconstruct the quantized frame exactly as the decoder would, and
     * use it to advance the prediction state.
     */
    {
        float Gm_q[9];
        memset(Gm_q, 0, sizeof(Gm_q));
        Gm_q[1] = 0.0f;
        for (int m = 0; m < 3; m++) {
            Gm_q[2 + m] = AmbePlusPRBA24[b3][m];
        }
        for (int m = 0; m < 4; m++) {
            Gm_q[5 + m] = AmbePlusPRBA58[b4][m];
        }

        for (int i = 1; i <= 8; i++) {
            float sum = 0.0f;
            for (int m = 1; m <= 8; m++) {
                int am = (m == 1) ? 1 : 2;
                sum += (float)am * Gm_q[m] * cache->prba_cos[m][i];
            }
            Ri_q[i] = sum;
        }

        const float rconst = 1.0f / (2.0f * sqrt2);
        Cik_q[1][1] = 0.5f * (Ri_q[1] + Ri_q[2]);
        Cik_q[1][2] = rconst * (Ri_q[1] - Ri_q[2]);
        Cik_q[2][1] = 0.5f * (Ri_q[3] + Ri_q[4]);
        Cik_q[2][2] = rconst * (Ri_q[3] - Ri_q[4]);
        Cik_q[3][1] = 0.5f * (Ri_q[5] + Ri_q[6]);
        Cik_q[3][2] = rconst * (Ri_q[5] - Ri_q[6]);
        Cik_q[4][1] = 0.5f * (Ri_q[7] + Ri_q[8]);
        Cik_q[4][2] = rconst * (Ri_q[7] - Ri_q[8]);

        {
            const float (*tables[4])[4] = {AmbePlusHOCb5, AmbePlusHOCb6, AmbePlusHOCb7, AmbePlusHOCb8};
            int codes[4] = {b5, b6, b7, b8};
            for (int blk = 0; blk < 4; blk++) {
                int ji = Ji[blk + 1];
                int kmax = (ji < 6) ? ji : 6;
                for (int k = 3; k <= kmax; k++) {
                    Cik_q[blk + 1][k] = tables[blk][codes[blk]][k - 3];
                }
            }
        }

        /* Block IDCT back to Tl_q */
        {
            int l = 1;
            for (int blk = 1; blk <= 4; blk++) {
                int ji = Ji[blk];
                for (int j = 1; j <= ji; j++) {
                    float sum = 0.0f;
                    for (int k = 1; k <= ji; k++) {
                        int ak = (k == 1) ? 1 : 2;
                        sum += (float)ak * Cik_q[blk][k] * cache->blk_cos[ji][j][k];
                    }
                    Tl_q[l] = sum;
                    l++;
                }
            }
        }

        /* Mirror of ambe2400_update_spectral_amplitudes() */
        float mean_tlq = 0.0f;
        for (int l = 1; l <= L; l++) {
            mean_tlq += Tl_q[l];
        }
        mean_tlq /= (float)L;

        for (int l = 1; l <= L; l++) {
            log2Ml_q[l] = Tl_q[l] - mean_tlq + (0.65f * (p[l] - mean_p)) + gamma_q - (0.5f * log2f((float)L));
        }
    }

    /* Pack the 49-bit AMBE word */
    memset(ambe_d, 0, 49);
    for (int i = 0; i < 6; i++) {
        ambe_d[i] = (char)((b0 >> (6 - i)) & 1);
    }
    ambe_d[48] = (char)(b0 & 1);

    ambe_d[38] = (char)((b1 >> 3) & 1);
    ambe_d[39] = (char)((b1 >> 2) & 1);
    ambe_d[40] = (char)((b1 >> 1) & 1);
    ambe_d[41] = (char)(b1 & 1);

    ambe_d[6] = (char)((b2 >> 5) & 1);
    ambe_d[7] = (char)((b2 >> 4) & 1);
    ambe_d[8] = (char)((b2 >> 3) & 1);
    ambe_d[9] = (char)((b2 >> 2) & 1);
    ambe_d[42] = (char)((b2 >> 1) & 1);
    ambe_d[43] = (char)(b2 & 1);

    ambe_d[10] = (char)((b3 >> 8) & 1);
    ambe_d[11] = (char)((b3 >> 7) & 1);
    ambe_d[12] = (char)((b3 >> 6) & 1);
    ambe_d[13] = (char)((b3 >> 5) & 1);
    ambe_d[14] = (char)((b3 >> 4) & 1);
    ambe_d[15] = (char)((b3 >> 3) & 1);
    ambe_d[16] = (char)((b3 >> 2) & 1);
    ambe_d[44] = (char)((b3 >> 1) & 1);
    ambe_d[45] = (char)(b3 & 1);

    ambe_d[17] = (char)((b4 >> 6) & 1);
    ambe_d[18] = (char)((b4 >> 5) & 1);
    ambe_d[19] = (char)((b4 >> 4) & 1);
    ambe_d[20] = (char)((b4 >> 3) & 1);
    ambe_d[21] = (char)((b4 >> 2) & 1);
    ambe_d[46] = (char)((b4 >> 1) & 1);
    ambe_d[47] = (char)(b4 & 1);

    ambe_d[22] = (char)((b5 >> 3) & 1);
    ambe_d[23] = (char)((b5 >> 2) & 1);
    ambe_d[25] = (char)((b5 >> 1) & 1);
    ambe_d[26] = (char)(b5 & 1);

    ambe_d[27] = (char)((b6 >> 3) & 1);
    ambe_d[28] = (char)((b6 >> 2) & 1);
    ambe_d[29] = (char)((b6 >> 1) & 1);
    ambe_d[30] = (char)(b6 & 1);

    ambe_d[31] = (char)((b7 >> 3) & 1);
    ambe_d[32] = (char)((b7 >> 2) & 1);
    ambe_d[33] = (char)((b7 >> 1) & 1);
    ambe_d[34] = (char)(b7 & 1);

    ambe_d[35] = (char)((b8 >> 3) & 1);
    ambe_d[36] = (char)((b8 >> 2) & 1);
    ambe_d[37] = (char)((b8 >> 1) & 1);

    /* ambe_d[24] is the spare bit; filled by the FEC layer below. */

    /* Write quantized parameters to cur_mp (decoder-equivalent state) */
    cur_mp->w0 = f0q * (float)(2.0 * M_PI);
    cur_mp->L = L;
    cur_mp->K = 0;
    cur_mp->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;
    for (int l = 1; l <= L; l++) {
        int jl = (int)((float)l * 16.0f * f0q);
        if (jl > 7) {
            jl = 7;
        }
        cur_mp->Vl[l] = AmbePlusVuv[b1][jl];
        if (cur_mp->Vl[l] == 1) {
            cur_mp->K++;
        }
        cur_mp->log2Ml[l] = log2Ml_q[l];
        cur_mp->Ml[l] = exp2f(log2Ml_q[l]);
        if (cur_mp->Vl[l] == 0) {
            cur_mp->Ml[l] *= 0.2046f / sqrtf(cur_mp->w0);
        }
    }
    cur_mp->gamma = gamma_q;

    return 0;
}

/*
 * Silence frame: b0 = 127, tone index = 128 (the standard AMBE silence
 * pattern). The decoder resets its parameter state on this frame, so the
 * encoder resets its prediction state as well.
 */
static void
ambe2400_encode_silence(char ambe_d[49], mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    (void)prev_mp;
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
mbe_encodeAmbe2400Parms(const float* samples, char ambe_d[49], mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    float agc_buf[AMBE2400_ENC_SAMPLES];
    float rms = 0.0f;
    int ret;

    if (samples == NULL || ambe_d == NULL || cur_mp == NULL || prev_mp == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    if (ambe2400_enc_fft_init() < 0) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        rms += samples[i] * samples[i];
    }
    rms = sqrtf(rms / (float)AMBE2400_ENC_SAMPLES);

    /* Silence gate with hang-over: a frame is only "silence" after a few
     * consecutive quiet frames (squelch close delay); speech opens the gate
     * immediately. This avoids the squelch-like flapping between silence and
     * speech on boundary frames. In-speech noise is handled by the optional
     * libspecbleach front-end in the caller, not here. */
    bool is_silence = false;
    if (rms < AMBE2400_ENC_SILENCE_RMS) {
        if (ambe2400_enc_silence_run < 5) {
            ambe2400_enc_silence_run++;
        }
        is_silence = (ambe2400_enc_silence_run >= 5);
    } else {
        ambe2400_enc_silence_run = 0;
    }

    /* Input AGC: track level only on speech frames, normalize to the
     * nominal target. Noise frames would drag the level estimate down and
     * over-drive the speech. */
    if (!is_silence && rms > 1e-6f) {
        ambe2400_enc_agc_rms =
            (AMBE2400_ENC_AGC_ALPHA * ambe2400_enc_agc_rms) + ((1.0f - AMBE2400_ENC_AGC_ALPHA) * rms);
    }
    ambe2400_enc_agc_gain = AMBE2400_ENC_AGC_TARGET / (ambe2400_enc_agc_rms + 1e-9f);
    if (ambe2400_enc_agc_gain < AMBE2400_ENC_AGC_MIN) {
        ambe2400_enc_agc_gain = AMBE2400_ENC_AGC_MIN;
    }
    if (ambe2400_enc_agc_gain > AMBE2400_ENC_AGC_MAX) {
        ambe2400_enc_agc_gain = AMBE2400_ENC_AGC_MAX;
    }
    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        agc_buf[i] = samples[i] * ambe2400_enc_agc_gain;
    }

    if (is_silence) {
        ambe2400_encode_silence(ambe_d, cur_mp, prev_mp);
    } else {
        ret = ambe2400_encode_voice(agc_buf, ambe_d, cur_mp, prev_mp);
        if (ret < 0) {
            return ret;
        }
    }

    memcpy(ambe2400_enc_history, agc_buf, (size_t)AMBE2400_ENC_SAMPLES * sizeof(float));
    ambe2400_enc_hist_gain = ambe2400_enc_agc_gain;

    return is_silence ? 1 : 0;
}

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of 16-bit PCM into AMBE 2400
 *        parameter bits.
 *
 * @see mbe_encodeAmbe2400Parms for details.
 */
int
mbe_encodeAmbe2400ParmsShort(const short* samples, char ambe_d[49], mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    float float_buf[AMBE2400_ENC_SAMPLES];

    if (samples == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < AMBE2400_ENC_SAMPLES; i++) {
        float_buf[i] = (float)samples[i] / 32768.0f;
    }
    return mbe_encodeAmbe2400Parms(float_buf, ambe_d, cur_mp, prev_mp);
}

/*
 * D-STAR DV frame interleave tables (air bit -> plane/index), from the
 * D-STAR specification. Air bit i (0..71, after the 24-bit sync word)
 * carries ambe_fr[dW[i]][dX[i]].
 */
static const int ambe2400_enc_dW[72] = {
    0, 0, 3, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, 1, 3, 2, 1, 1, 0, 0, 3, 2, 0, 0, 3, 2, 1, 1, 0, 0, 1, 1, 0, 0,
    3, 2, 1, 1, 3, 2, 1, 1, 0, 0, 3, 2, 0, 0, 3, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, 1, 3, 3, 2, 1, 0, 0, 3, 3,
};

static const int ambe2400_enc_dX[72] = {
    10, 22, 11, 9, 10, 22, 11, 23, 8, 20, 9, 21, 10, 8, 9, 21, 8, 6,  7,  19, 8, 20, 9, 7,
    6,  18, 7,  5, 6,  18, 7,  19, 4, 16, 5, 17, 6,  4, 5, 17, 4, 2,  3,  15, 4, 16, 5, 3,
    2,  14, 3,  1, 2,  14, 3,  15, 0, 12, 1, 13, 2,  0, 1, 13, 0, 12, 10, 11, 0, 12, 1, 13,
};

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
    char cw[23];
    char data12[12];
    unsigned short pr[25];
    unsigned short foo = 0;
    int ones;

    if (ambe_d == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 49; i++) {
        if (ambe_d[i] != 0 && ambe_d[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }

    memset(ambe_fr, 0, 4 * 24 * sizeof(char));

    /* C0: 12 data bits (ambe_d[0..11]) -> (23,12) Golay + even parity */
    for (int i = 0; i < 12; i++) {
        data12[i] = ambe_d[i];
    }
    mbe_golay2312_encode(data12, cw);
    for (int j = 0; j < 23; j++) {
        ambe_fr[0][j + 1] = cw[j];
    }
    ones = 0;
    for (int j = 0; j < 24; j++) {
        ones += (ambe_fr[0][j] & 1);
    }
    ambe_fr[0][0] = (char)(ones & 1);

    /* PR scramble sequence from the C0 data word */
    for (int i = 23; i >= 12; i--) {
        foo <<= 1;
        foo |= (unsigned short)(ambe_fr[0][i] & 1);
    }
    pr[0] = (unsigned short)(16 * foo);
    for (int i = 1; i < 25; i++) {
        pr[i] = (unsigned short)((173 * pr[i - 1]) + 13849 - (65536 * (((173 * pr[i - 1]) + 13849) / 65536)));
    }

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

    /* C2/C3: raw bits (decoder reads them MSB-first) */
    for (int j = 0; j < 11; j++) {
        ambe_fr[2][j] = ambe_d[24 + (10 - j)];
    }
    for (int j = 0; j < 14; j++) {
        ambe_fr[3][j] = ambe_d[35 + (13 - j)];
    }

    /* Spare bit (ambe_d[24] = fr[2][10] on air): b word bit 0 = parity ^ p bit 0 */
    ones = 0;
    for (int j = 0; j < 23; j++) {
        ones += (cw[j] & 1);
    }
    ambe_fr[2][10] = (char)(((ones & 1) ^ (pr[24] / 32768)) & 1);

    return 0;
}

/**
 * @brief Serialize a 72-bit AMBE 3600x2400 frame into the 9 data bytes
 *        of a D-STAR DV frame (the 24-bit sync word is NOT included).
 *
 * Bytes are packed in air order (LSB first within each byte, matching
 * the GMSK modulator). This is the exact inverse of
 * mbe_decodeDStarDVData().
 *
 * @param ambe_fr  Input frame as 4x24 bitplanes.
 * @param bytes9   Output 9 bytes (72 bits).
 * @return 0 on success.
 */
int
mbe_encodeDStarDVData(const char ambe_fr[4][24], unsigned char bytes9[9]) {
    if (ambe_fr == NULL || bytes9 == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    memset(bytes9, 0, 9);
    for (int i = 0; i < 72; i++) {
        if (ambe_fr[ambe2400_enc_dW[i]][ambe2400_enc_dX[i]] & 1) {
            bytes9[i >> 3] |= (unsigned char)(1u << (i & 7));
        }
    }
    return 0;
}

/**
 * @brief Extract a 72-bit AMBE 3600x2400 frame from the 9 data bytes of
 *        a D-STAR DV frame.
 *
 * @see mbe_encodeDStarDVData for the byte/bit convention.
 *
 * @param bytes9  Input 9 bytes (72 bits).
 * @param ambe_fr Output frame as 4x24 bitplanes.
 * @return 0 on success.
 */
int
mbe_decodeDStarDVData(const unsigned char bytes9[9], char ambe_fr[4][24]) {
    if (bytes9 == NULL || ambe_fr == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }

    memset(ambe_fr, 0, 4 * 24 * sizeof(char));
    for (int i = 0; i < 72; i++) {
        ambe_fr[ambe2400_enc_dW[i]][ambe2400_enc_dX[i]] = (char)((bytes9[i >> 3] >> (i & 7)) & 1u);
    }
    return 0;
}
