// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Shared AMBE 2400 DCT and spectral reconstruction (private API). */
#ifndef MBELIB_NEO_INTERNAL_AMBE3600X2400_H
#define MBELIB_NEO_INTERNAL_AMBE3600X2400_H

#include "mbelib-neo/mbelib.h"

struct ambe_dct_cache {
    int inited;
    float ri_cos[9][9];         /* [m][i] for m=1..8, i=1..8 (index 0 unused) */
    float idct_cos[18][18][18]; /* [ji][j][k] for ji=1..17, j=1..ji, k=1..ji */
};

/* Hidden by the library's default visibility, like ambe_common.h helpers.
 * Arrays use one-based harmonic/block indices. L and prev_L must be clamped
 * to 1..56; Ji comes from AmbePlusLmprbl and codes contains b5..b8. */
const struct ambe_dct_cache* ambe2400_get_dct_cache(void);
void ambe2400_reconstruct_prba(int b3, int b4, float Ri[9]);
void ambe2400_reconstruct_ri(const float Gm[9], float Ri[9]);
void ambe2400_reconstruct_cik(const float Ri[9], const int codes[4], const int Ji[5], float Cik[5][18]);
void ambe2400_inverse_dct_tl(float Cik[5][18], const int Ji[5], float Tl[57]);

/* Preserve the decoder's full update loop, including its previous-spectrum
 * fix-ups and Ml scaling. The encoder passes a private copy of prev_mp. */
void ambe2400_update_spectral_amplitudes(mbe_parms* cur_mp, mbe_parms* prev_mp, const float Tl[57], float unvc);

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
