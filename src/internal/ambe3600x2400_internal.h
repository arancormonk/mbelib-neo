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
    return ((float)prev_L / (float)L) * (float)l;
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
