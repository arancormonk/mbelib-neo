// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 *
 * Copyright (C) 2010 mbelib Author
 * GPG Key ID: 0xEA5EFE2C (9E7A 5527 9CDC EBF7 BF1B  D772 4F98 E863 EA5E FE2C)
 *
 * Portions were originally under the ISC license; this mbelib-neo
 * distribution is provided under GPL-2.0-or-later. See LICENSE for details.
 */

/** @file Shared AMBE 2400 DCT and spectral reconstruction (private API). */
#ifndef MBELIB_NEO_INTERNAL_AMBE3600X2400_H
#define MBELIB_NEO_INTERNAL_AMBE3600X2400_H

#include <math.h>

#include "mbe_math.h"
#include "mbelib-neo/mbelib.h"

/*
 * Weight of the previous frame's log2 spectral amplitudes in the D-STAR
 * prediction (log2Ml = Tl + rho * interpolated previous - mean terms).
 * mbelib took 0.65 from AMBE+2 (TIA-102.BABA-1). Decoding DVSI's own
 * AMBE-3000 D-STAR test vectors with 0.65 compresses the spectral envelope to
 * about 60% in the log domain and plays the 2-4 kHz bands 8-10 dB hotter than
 * DVSI's decoder relative to 0-1 kHz; 0.80 matches DVSI's output within about
 * 1 dB per band. The encoder predicts with the same value so that its frames
 * reconstruct the intended envelope in DVSI decoders.
 */
#define MBE_AMBE2400_PREDICTION_RHO       0.8f

/*
 * D-STAR's gain is absolute: gamma = SCALE * AmbePlusDg[b2] + OFFSET, with no
 * term from the previous frame. mbelib carried over rate 33's differential law
 * (TIA-102.BABA-1 eq 26, Dg[b2] + 0.5 gamma(-1)), whose steady state is the
 * same 2 Dg, so levels matched on average while onsets after silence came out
 * up to 2.6 log2 (16 dB) low. In DVSI's re-quantizations of its decoded
 * D-STAR frames into rate 33, the gain implied for fully voiced frames sits a
 * median 0.037 log2 from this law per index (0.35 from the differential one),
 * and the offset is where the rate-33 gain indices agree best.
 * The encoder quantizes the gain against the same levels.
 */
#define MBE_AMBE2400_GAIN_SCALE           2.0f
#define MBE_AMBE2400_GAIN_OFFSET          (-0.09f)

/*
 * D-STAR mutes a hard-decision stream while its error rate (0.95 eps_R +
 * 0.001064 eps_T, as rate 33 keeps it, with D-STAR's own counts) is above
 * this. DVSI's AMBE-3000 mutes D-STAR at sustained channel error rates near
 * 10% and never at 5% or below; on its 1-10% test vectors this threshold
 * agrees with DVSI's muting on 100% of the loud frames at 1-5% and 90% at
 * 10%, where it starts muting earlier than DVSI does (0.7 s into a constant
 * 10% channel, against DVSI's 3.3 s). It approximates DVSI's rule,
 * which integrates the errors more slowly than this rate does. Soft decisions
 * do not mute on it: at 10% DVSI mutes 53% of the loud frames there, against
 * 88% with hard decisions, and no rule tried reproduced that. A muted frame
 * is synthesized as usual (model, phase, overlap and noise state all advance)
 * and played as mute noise at DVSI's level (MBE_AMBE2400_MUTE_NOISE_AMPLITUDE).
 */
#define MBE_AMBE2400_MUTE_ERROR_RATE      0.042f

/* Amplitude (int16 scale) of the uniform noise D-STAR mutes to. The int16
 * conversion truncates, so it is 0 on two thirds of the samples and +-1 on the
 * rest, as DVSI's muted D-STAR output is (dam_e10_hd: 67% 0, 16.5% each +-1). */
#define MBE_AMBE2400_MUTE_NOISE_AMPLITUDE 1.5f

/*
 * D-STAR pitch law: f0 = 2^(OFFSET - SLOPE * (b0 + 0.5)) cycles per sample.
 * mbelib's -4.311767578125 - 0.021336 * (b0 + 0.5) plays D-STAR flat, by 4%
 * at b0 13 down to 1% at b0 120. These constants fit the pitch of DVSI's own
 * decoded AMBE-3000 D-STAR test vectors (pulse trains and speech, b0 13..123)
 * within 1 cent; the pitch DVSI carries into its D-STAR to P25 rate conversion
 * and the pitch of the original input speech agree with them. The encoder
 * quantizes with the same law.
 */
#define MBE_AMBE2400_LOG2_F0_OFFSET       (-4.24738f)
#define MBE_AMBE2400_LOG2_F0_SLOPE        0.0217705f

/* Fundamental (cycles per sample) of pitch code b0. */
static inline float
mbe_ambe2400_f0(int b0) {
    return exp2f(MBE_AMBE2400_LOG2_F0_OFFSET - (MBE_AMBE2400_LOG2_F0_SLOPE * ((float)b0 + 0.5f)));
}

/*
 * Harmonic count for fundamental f0: every harmonic at or below 0.9254 of
 * Nyquist. The AMBE+2 L table is exactly floor(0.9254 * pi / w0) of its own w0
 * table; D-STAR's law gives a higher f0 for the same b0, and DVSI's D-STAR
 * decoder synthesizes harmonics only up to that same 3.70 kHz edge. Reusing the
 * AMBE+2 table put up to four harmonics too many into D-STAR's spectral blocks.
 * The block tables start at 9 harmonics, which caps f0 above 412 Hz (b0 0..1).
 */
static inline int
mbe_ambe2400_harmonic_count(float f0) {
    int L = (int)(0.9254f * 0.5f / f0);
    return (L < 9) ? 9 : ((L > 56) ? 56 : L);
}

/*
 * Tone frame (b0 126) index: bits 6..8 select its three most significant bits
 * through this table, and its five low bits sit at 9, 42, 43, 10 and 11. The
 * decoder reads and the encoder writes the index with these.
 */
static inline int
mbe_ambe2400_tone_high_bits(int selector) {
    static const unsigned char high[8] = {4, 0, 1, 2, 3, 7, 6, 5};
    return high[selector & 7];
}

static inline int
mbe_ambe2400_tone_index(const char* ambe_d) {
    int selector = (ambe_d[6] << 2) | (ambe_d[7] << 1) | ambe_d[8];
    return (mbe_ambe2400_tone_high_bits(selector) << 5) | (ambe_d[9] << 4) | (ambe_d[42] << 3) | (ambe_d[43] << 2)
           | (ambe_d[10] << 1) | ambe_d[11];
}

/*
 * What mbe_decodeAmbe2400Parms() returns besides 0 (voice) and a playable
 * tone index (5..122, 128..147). On DVSI's AMBE-3000, b0 125, a b0 126 tone
 * frame whose index has no tone, and each b0 127 frame measured play sample
 * for sample as a lost frame does, which repeats the last frame:
 * MBE_AMBE2400_FRAME_LOST. The index-128 fields of one of them play as a tone
 * under b0 126, so b0 127 is not read as a tone frame, and every b0 127 frame
 * is taken as lost. Among them is the all-ones DV frame (FF x 9), which DVSI's
 * host software writes into a bit file for a frame it does not send; only its
 * decoding side plays that frame as comfort noise, by flagging it to the chip
 * (docs/dvsi-parity.md).
 */
#define MBE_AMBE2400_FRAME_LOST   2

/* b0 124, the silence frame (the D-STAR null frame among them). */
#define MBE_AMBE2400_SILENCE_B0   124

/* The highest pitch code the encoder sends; DVSI's encoder sends none above it. */
#define MBE_AMBE2400_MAX_VOICE_B0 123

/* Write index 0..255 into a tone frame; every selector maps to a distinct
 * value, so each index has exactly one. */
static inline void
mbe_ambe2400_set_tone_index(char* ambe_d, int index) {
    int selector = 0;
    while (selector < 7 && mbe_ambe2400_tone_high_bits(selector) != ((index >> 5) & 7)) {
        selector++;
    }
    ambe_d[6] = (char)((selector >> 2) & 1);
    ambe_d[7] = (char)((selector >> 1) & 1);
    ambe_d[8] = (char)(selector & 1);
    ambe_d[9] = (char)((index >> 4) & 1);
    ambe_d[42] = (char)((index >> 3) & 1);
    ambe_d[43] = (char)((index >> 2) & 1);
    ambe_d[10] = (char)((index >> 1) & 1);
    ambe_d[11] = (char)(index & 1);
}

struct ambe_dct_cache {
    int inited;
    float ri_cos[9][9];         /* [m][i] for m=1..8, i=1..8 (index 0 unused) */
    float idct_cos[18][18][18]; /* [ji][j][k] for ji=1..17, j=1..ji, k=1..ji */
};

/* Not exported from the shared library; present as ordinary globals in the
 * static archive, hence the mbe_ prefix.
 *
 * Arrays use one-based harmonic/block indices. L and prev_L must be clamped
 * to 1..56; Ji comes from AmbePlusLmprbl and codes contains b5..b8. */
const struct ambe_dct_cache* mbe_ambe2400_get_dct_cache(void);
void mbe_ambe2400_reconstruct_prba(int b3, int b4, float Ri[9]);
void mbe_ambe2400_reconstruct_ri(const float Gm[9], float Ri[9]);
void mbe_ambe2400_reconstruct_cik(const float Ri[9], const int codes[4], const int Ji[5], float Cik[5][18]);
void mbe_ambe2400_inverse_dct_tl(float Cik[5][18], const int Ji[5], float Tl[57]);

/* Preserve the decoder's full update loop, including its previous-spectrum
 * fix-ups and Ml scaling. The encoder passes a private copy of prev_mp. */
void mbe_ambe2400_update_spectral_amplitudes(mbe_parms* cur_mp, mbe_parms* prev_mp, const float Tl[57], float unvc);

/* Small arithmetic helpers stay inline to preserve the decoder's evaluation
 * order before fast-math/LTO transforms its surrounding loops. */
static inline float
ambe2400_prediction_position(int prev_L, int L, int l) {
    return mbe_prediction_position(prev_L, L, l);
}

static inline float
ambe2400_interpolate_prediction(float deltal, float lower, float upper) {
    return (((float)1 - deltal) * lower) + (deltal * upper);
}

/* c1 and c2 are separately weighted prediction taps; BigGamma is
 * gamma - 0.5*log2(L) - mean(Tl). */
static inline float
ambe2400_reconstruct_log2Ml(float Tl, float c1, float c2, float Sum43, float BigGamma) {
    return Tl + c1 + c2 - Sum43 + BigGamma;
}

#endif /* MBELIB_NEO_INTERNAL_AMBE3600X2400_H */
