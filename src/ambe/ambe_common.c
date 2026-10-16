// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Internal helpers shared by AMBE 3600x2400 and 3600x2450 paths.
 *
 * Provides common ECC correction for C0, demodulation of C1 using a
 * pseudo-random sequence derived from C0, and packing of 49 AMBE
 * parameter bits from the four bitplanes.
 */

#include "ambe_common.h"
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "mbe_adaptive.h"
#include "mbe_ecc.h"
#include "mbelib-neo/mbelib.h"

int
mbe_eccAmbe3600C0_common(char fr[4][24]) {
    int j, errs;
    char in[23], out[23];
    for (j = 0; j < 23; j++) {
        in[j] = fr[0][j + 1];
    }
    errs = mbe_golay2312(in, out);
    for (j = 0; j < 23; j++) {
        fr[0][j + 1] = out[j];
    }
    /* JMBE C0 uses Golay24: when the protected 23-bit codeword has no
     * syndrome, enforce even parity by correcting the extra parity bit. */
    if (errs == 0) {
        int ones = 0;
        for (j = 0; j < 24; j++) {
            ones += (fr[0][j] & 1);
        }
        if ((ones & 1) != 0) {
            fr[0][0] ^= 1;
            errs = 1;
        }
    }
    return errs;
}

/* The Golay(23,12) codeword of the data bits a decode left in out[11..22]
 * (MSB at 22). */
static void
ambe_golay23_codeword(const char out[23], char codeword[23]) {
    char data[12];
    for (int i = 0; i < 12; i++) {
        data[i] = out[22 - i];
    }
    mbe_golay2312_encode(data, codeword);
}

/* That codeword, and how many of the received bits in[] differ from it. */
static int
ambe_golay23_errors(const char in[23], const char out[23], char codeword[23]) {
    ambe_golay23_codeword(out, codeword);
    int errs = 0;
    for (int i = 0; i < 23; i++) {
        errs += (codeword[i] & 1) != (in[i] & 1);
    }
    return errs;
}

int
mbe_eccAmbe3600C0Counted_common(char fr[4][24]) {
    char in[23], out[23], codeword[23];
    for (int j = 0; j < 23; j++) {
        in[j] = fr[0][j + 1];
    }
    (void)mbe_golay2312(in, out);
    int errs = ambe_golay23_errors(in, out, codeword);
    int parity = 0;
    for (int j = 0; j < 23; j++) {
        fr[0][j + 1] = codeword[j];
        parity ^= codeword[j] & 1;
    }
    if ((fr[0][0] & 1) != parity) {
        fr[0][0] = (char)parity;
        errs++;
    }
    return errs;
}

/* Rounded reliability cost of the received soft bits in[] that disagree with
 * codeword[], in units of a fully reliable bit. */
static int
ambe_soft_cost(const mbe_soft_bit* in, const char* codeword, int count, int extra) {
    int cost = extra;
    for (int i = 0; i < count; i++) {
        if ((in[i].bit & 1u) != (unsigned)(codeword[i] & 1)) {
            cost += in[i].reliability;
        }
    }
    return (cost + 127) / 255;
}

int
mbe_eccAmbe3600C0SoftCounted_common(mbe_soft_bit fr[4][24]) {
    mbe_soft_bit in[23];
    char out[23], codeword[23];
    for (int j = 0; j < 23; j++) {
        in[j] = fr[0][j + 1];
    }
    (void)mbe_golay2312Soft(in, out);
    ambe_golay23_codeword(out, codeword);
    int parity = 0;
    for (int j = 0; j < 23; j++) {
        fr[0][j + 1].bit = (uint8_t)(codeword[j] & 1);
        parity ^= codeword[j] & 1;
    }
    int extra = 0;
    if ((fr[0][0].bit & 1u) != (unsigned)parity) {
        extra = fr[0][0].reliability;
        fr[0][0].bit = (uint8_t)parity;
    }
    return ambe_soft_cost(in, codeword, 23, extra);
}

int
mbe_eccAmbe3600C0Soft_common(mbe_soft_bit fr[4][24]) {
    int j, errs;
    mbe_soft_bit in[23];
    char out[23];

    for (j = 0; j < 23; j++) {
        in[j] = fr[0][j + 1];
    }
    errs = mbe_golay2312Soft(in, out);
    for (j = 0; j < 23; j++) {
        fr[0][j + 1].bit = (uint8_t)(out[j] & 1);
    }

    if (errs == 0) {
        int ones = 0;
        for (j = 0; j < 24; j++) {
            ones += (int)(fr[0][j].bit & 1u);
        }
        if ((ones & 1) != 0) {
            fr[0][0].bit ^= 1u;
            errs = 1;
        }
    }
    return errs;
}

void
mbe_demodulateAmbe3600Data_common(char fr[4][24]) {
    int i, j, k;
    unsigned short pr[115];
    unsigned short foo = 0;

    /* create pseudo-random modulator */
    for (i = 23; i >= 12; i--) {
        foo <<= 1;
        foo |= fr[0][i];
    }
    pr[0] = (unsigned short)(16 * foo);
    for (i = 1; i < 24; i++) {
        pr[i] = (unsigned short)((173 * pr[i - 1]) + 13849 - (65536 * (((173 * pr[i - 1]) + 13849) / 65536)));
    }
    for (i = 1; i < 24; i++) {
        pr[i] = (unsigned short)(pr[i] / 32768);
    }

    /* demodulate fr with pr */
    k = 1;
    for (j = 22; j >= 0; j--) {
        fr[1][j] = (char)(fr[1][j] ^ pr[k]);
        k++;
    }
}

void
mbe_demodulateAmbe3600DataSoft_common(mbe_soft_bit fr[4][24]) {
    int i, j, k;
    unsigned short pr[115];
    unsigned short foo = 0;

    for (i = 23; i >= 12; i--) {
        foo <<= 1;
        foo |= (unsigned short)(fr[0][i].bit & 1u);
    }
    pr[0] = (unsigned short)(16 * foo);
    for (i = 1; i < 24; i++) {
        pr[i] = (unsigned short)((173 * pr[i - 1]) + 13849 - (65536 * (((173 * pr[i - 1]) + 13849) / 65536)));
    }
    for (i = 1; i < 24; i++) {
        pr[i] = (unsigned short)(pr[i] / 32768);
    }

    k = 1;
    for (j = 22; j >= 0; j--) {
        fr[1][j].bit = (uint8_t)((fr[1][j].bit & 1u) ^ (uint8_t)pr[k]);
        k++;
    }
}

int
mbe_eccAmbe3600Data_common(char fr[4][24], char* out49) {
    int j, errs;
    char *ambe, gin[24], gout[24];
    ambe = out49;
    /* just copy C0 */
    for (j = 23; j > 11; j--) {
        *ambe = fr[0][j];
        ambe++;
    }
    /* ecc and copy C1 */
    for (j = 0; j < 23; j++) {
        gin[j] = fr[1][j];
    }
    errs = mbe_golay2312(gin, gout);
    for (j = 22; j > 10; j--) {
        *ambe = gout[j];
        ambe++;
    }
    /* just copy C2 */
    for (j = 10; j >= 0; j--) {
        *ambe = fr[2][j];
        ambe++;
    }
    /* just copy C3 */
    for (j = 13; j >= 0; j--) {
        *ambe = fr[3][j];
        ambe++;
    }
    return errs;
}

int
mbe_eccAmbe3600DataCounted_common(char fr[4][24], char* out49) {
    char in[23], codeword[23];
    /* out49 may start at fr: packing C2 overwrites the received C1 parity. */
    memcpy(in, fr[1], sizeof(in));
    (void)mbe_eccAmbe3600Data_common(fr, out49);
    /* C1's decoded data is already in out49[12..23], MSB first. Re-encode
     * it to count parity errors too, without decoding the same word twice. */
    mbe_golay2312_encode(out49 + 12, codeword);
    int errs = 0;
    for (int j = 0; j < 23; j++) {
        errs += (codeword[j] & 1) != (in[j] & 1);
    }
    return errs;
}

int
mbe_eccAmbe3600DataSoftCounted_common(mbe_soft_bit fr[4][24], char* out49) {
    mbe_soft_bit in[23];
    char codeword[23];
    memcpy(in, fr[1], sizeof(in));
    (void)mbe_eccAmbe3600DataSoft_common(fr, out49);
    mbe_golay2312_encode(out49 + 12, codeword);
    return ambe_soft_cost(in, codeword, 23, 0);
}

int
mbe_eccAmbe3600DataSoft_common(mbe_soft_bit fr[4][24], char* out49) {
    int j, errs;
    char *ambe, gout[23];
    mbe_soft_bit gin[23];

    ambe = out49;
    for (j = 23; j > 11; j--) {
        *ambe = (char)(fr[0][j].bit & 1u);
        ambe++;
    }

    for (j = 0; j < 23; j++) {
        gin[j] = fr[1][j];
    }
    errs = mbe_golay2312Soft(gin, gout);
    for (j = 22; j > 10; j--) {
        *ambe = gout[j];
        ambe++;
    }

    for (j = 10; j >= 0; j--) {
        *ambe = (char)(fr[2][j].bit & 1u);
        ambe++;
    }
    for (j = 13; j >= 0; j--) {
        *ambe = (char)(fr[3][j].bit & 1u);
        ambe++;
    }
    return errs;
}

void
mbe_setAmbeDefaultModel_common(mbe_parms* mp) {
    if (!mp) {
        return;
    }

    mp->w0 = MBE_AMBE_INIT_W0;
    mp->L = MBE_AMBE_INIT_L;
    mp->K = 0;
    mp->gamma = 0.0f;
    for (int l = 0; l <= 56; l++) {
        mp->Ml[l] = 1.0f;
        mp->Vl[l] = 0;
        mp->log2Ml[l] = 0.0f; /* log2(1.0) */
    }
}

void
mbe_initAmbeParms_common(mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    if (!cur_mp || !prev_mp || !prev_mp_enhanced) {
        return;
    }

    prev_mp->swn = 0;
    prev_mp->tonePhase = 0;
    mbe_setAmbeDefaultModel_common(prev_mp);
    for (int l = 0; l <= 56; l++) {
        prev_mp->PHIl[l] = 0.0f;
        prev_mp->PSIl[l] = 0.0f;
    }

    prev_mp->localEnergy = MBE_DEFAULT_LOCAL_ENERGY;
    prev_mp->amplitudeThreshold = MBE_DEFAULT_AMPLITUDE_THRESHOLD;
    prev_mp->errorRate = 0.0f;
    prev_mp->errorCountTotal = 0;
    prev_mp->errorCount4 = 0;

    prev_mp->repeatCount = 0;
    prev_mp->mutingThreshold = MBE_MUTING_THRESHOLD_AMBE;

    prev_mp->noiseSeed = -1.0f;
    memset(prev_mp->noiseOverlap, 0, sizeof(prev_mp->noiseOverlap));
    memset(prev_mp->previousUw, 0, sizeof(prev_mp->previousUw));

    *cur_mp = *prev_mp;
    *prev_mp_enhanced = *prev_mp;
}

void
mbe_setAmbeErasureParms_common(mbe_parms* mp, const mbe_parms* state_src) {
    if (!mp) {
        return;
    }

    const mbe_parms* continuity = state_src ? state_src : mp;

    mp->swn = 0;
    mp->tonePhase = 0;
    mp->w0 = 0.0f; /* JMBE AMBEFundamentalFrequency.W120..W123 */
    mp->L = 9;
    mp->K = 0;
    mp->gamma = 0.0f;

    for (int l = 0; l <= 56; l++) {
        mp->Ml[l] = 1.0f;
        mp->Vl[l] = 0;
        mp->log2Ml[l] = 0.0f; /* log2(1.0) */
        mp->PHIl[l] = continuity->PHIl[l];
        mp->PSIl[l] = continuity->PSIl[l];
    }

    mp->localEnergy = MBE_DEFAULT_LOCAL_ENERGY;
    mp->amplitudeThreshold = MBE_DEFAULT_AMPLITUDE_THRESHOLD;

    mp->noiseSeed = continuity->noiseSeed;
    memcpy(mp->noiseOverlap, continuity->noiseOverlap, sizeof(mp->noiseOverlap));
    memcpy(mp->previousUw, continuity->previousUw, sizeof(mp->previousUw));
}

void
mbe_ensureAmbeDefaults_common(mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced) {
    if (!prev_mp) {
        return;
    }

    if (fabsf(prev_mp->mutingThreshold - MBE_MUTING_THRESHOLD_AMBE) > 1e-6f) {
        mbe_initAmbeParms_common(cur_mp, prev_mp, prev_mp_enhanced);
    }
}
