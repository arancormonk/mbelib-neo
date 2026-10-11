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
#define MBE_AMBE2400_PREDICTION_RHO 0.8f

/*
 * D-STAR pitch law: f0 = 2^(OFFSET - SLOPE * (b0 + 0.5)) cycles per sample.
 * mbelib's -4.311767578125 - 0.021336 * (b0 + 0.5) plays D-STAR flat, by 4%
 * at b0 13 down to 1% at b0 120. These constants fit the pitch of DVSI's own
 * decoded AMBE-3000 D-STAR test vectors (pulse trains and speech, b0 13..123)
 * within 1 cent; the pitch DVSI carries into its D-STAR to P25 rate conversion
 * and the pitch of the original input speech agree with them. The encoder
 * quantizes with the same law.
 */
#define MBE_AMBE2400_LOG2_F0_OFFSET (-4.24738f)
#define MBE_AMBE2400_LOG2_F0_SLOPE  0.0217705f

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
