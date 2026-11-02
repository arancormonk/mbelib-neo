// SPDX-License-Identifier: GPL-2.0-or-later
/* Offline reference formulas from the codec initializers. Not linked into the library.
 * Analysis-window formulas: Copyright (C) 2026 by arancormonk.
 * See docs/encoder-tables.md for build and regeneration instructions. */
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#define REFERENCE_PI            3.14159265358979323846
#define MBE_ANALYSIS_PITCH_HALF 150
#define MBE_ANALYSIS_WR_HALF    110
#define MBE_ANALYSIS_LPF_HALF   10
#define MBE_ANALYSIS_WTAB_BINS  7
#define MBE_ANALYSIS_WTAB_STEPS 64
#define ANALYSIS_FFT            256
#define W_I_BETA                5.25
#define W_R_BETA                6.0
#define LPF_CUTOFF              0.1757
#define LPF_BETA                4.18

struct mbe_analysis_tables {
    float w_i[301], w_r[221], lpf[21], w_r_dtft[897];
    float w_i4_sum, w_r_sum, w_r_sq_sum;
};

static void
emit(const char* group, float value) {
    uint32_t bits;
    memcpy(&bits, &value, sizeof bits);
    printf("%s %08" PRIx32 "\n", group, bits);
}

static double
bessel_i0(double x) {
    double sum = 1.0;
    double term = 1.0;
    for (int k = 1; k < 50; ++k) {
        term *= (x / (2.0 * (double)k)) * (x / (2.0 * (double)k));
        sum += term;
        if (term < sum * 1e-17) {
            break;
        }
    }
    return sum;
}

static double
kaiser(int n, int half, double beta) {
    double r = (double)n / (double)half;
    return bessel_i0(beta * sqrt(fmax(0.0, 1.0 - (r * r)))) / bessel_i0(beta);
}

static void
init_windows(struct mbe_analysis_tables* t) {
    double sum_sq = 0.0;
    for (int n = -MBE_ANALYSIS_PITCH_HALF; n <= MBE_ANALYSIS_PITCH_HALF; ++n) {
        double w = kaiser(n, MBE_ANALYSIS_PITCH_HALF, W_I_BETA);
        sum_sq += w * w;
    }
    double w_i4 = 0.0;
    for (int n = -MBE_ANALYSIS_PITCH_HALF; n <= MBE_ANALYSIS_PITCH_HALF; ++n) {
        double w = kaiser(n, MBE_ANALYSIS_PITCH_HALF, W_I_BETA) / sqrt(sum_sq); /* eq (6): sum w_I^2 = 1 */
        t->w_i[n + MBE_ANALYSIS_PITCH_HALF] = (float)w;
        w_i4 += w * w * w * w;
    }
    t->w_i4_sum = (float)w_i4;

    double w_r_sum = 0.0;
    double w_r_sq = 0.0;
    for (int n = -MBE_ANALYSIS_WR_HALF; n <= MBE_ANALYSIS_WR_HALF; ++n) {
        double w = kaiser(n, MBE_ANALYSIS_WR_HALF, W_R_BETA);
        t->w_r[n + MBE_ANALYSIS_WR_HALF] = (float)w;
        w_r_sum += w;
        w_r_sq += w * w;
    }
    t->w_r_sum = (float)w_r_sum;
    t->w_r_sq_sum = (float)w_r_sq;
}

static void
init_lowpass(struct mbe_analysis_tables* t) {
    double taps[(2 * MBE_ANALYSIS_LPF_HALF) + 1];
    double sum = 0.0;
    for (int n = -MBE_ANALYSIS_LPF_HALF; n <= MBE_ANALYSIS_LPF_HALF; ++n) {
        double x = 2.0 * REFERENCE_PI * LPF_CUTOFF * (double)n;
        double sinc = (n == 0) ? 1.0 : sin(x) / x;
        taps[n + MBE_ANALYSIS_LPF_HALF] = 2.0 * LPF_CUTOFF * sinc * kaiser(n, MBE_ANALYSIS_LPF_HALF, LPF_BETA);
        sum += taps[n + MBE_ANALYSIS_LPF_HALF];
    }
    for (int i = 0; i < (2 * MBE_ANALYSIS_LPF_HALF) + 1; ++i) {
        t->lpf[i] = (float)(taps[i] / sum); /* unit gain at DC */
    }
}

/* W_R(nu) of eq (30) for nu in bins of the 256-point DFT: real, since w_R is
 * symmetric about n = 0. Zero beyond the table, where bands never reach. */
static void
init_window_transform(struct mbe_analysis_tables* t) {
    const int centre = MBE_ANALYSIS_WTAB_BINS * MBE_ANALYSIS_WTAB_STEPS;
    /* Both w_R and its transform are even: fold the sum and mirror the table. */
    for (int i = 0; i <= centre; ++i) {
        double nu = (double)i / (double)MBE_ANALYSIS_WTAB_STEPS;
        double sum = (double)t->w_r[MBE_ANALYSIS_WR_HALF];
        for (int n = 1; n <= MBE_ANALYSIS_WR_HALF; ++n) {
            sum += 2.0 * (double)t->w_r[n + MBE_ANALYSIS_WR_HALF]
                   * cos(2.0 * REFERENCE_PI * nu * (double)n / ANALYSIS_FFT);
        }
        t->w_r_dtft[centre + i] = (float)sum;
        t->w_r_dtft[centre - i] = (float)sum;
    }
}

int
main(void) {
    if (sizeof(float) != sizeof(uint32_t) || FLT_RADIX != 2 || FLT_MANT_DIG != 24 || FLT_MAX_EXP != 128) {
        fputs("The reference requires IEEE binary32 float.\n", stderr);
        return 1;
    }
    struct mbe_analysis_tables t;
    init_windows(&t);
    init_lowpass(&t);
    init_window_transform(&t);
    for (int i = 0; i < 301; i++) {
        emit("analysis", t.w_i[i]);
    }
    for (int i = 0; i < 221; i++) {
        emit("analysis", t.w_r[i]);
    }
    for (int i = 0; i < 21; i++) {
        emit("analysis", t.lpf[i]);
    }
    for (int i = 0; i < 897; i++) {
        emit("analysis", t.w_r_dtft[i]);
    }
    emit("analysis", t.w_i4_sum);
    emit("analysis", t.w_r_sum);
    emit("analysis", t.w_r_sq_sum);
    for (int m = 1; m <= 8; m++) {
        for (int i = 1; i <= 8; i++) {
            emit("ambe", cosf((REFERENCE_PI * (float)(m - 1) * ((float)i - 0.5f)) / 8.0f));
        }
    }
    for (int ji = 1; ji <= 17; ji++) {
        for (int j = 1; j <= ji; j++) {
            for (int k = 1; k <= ji; k++) {
                emit("ambe", cosf((REFERENCE_PI * (float)(k - 1) * ((float)j - 0.5f)) / (float)ji));
            }
        }
    }
    for (int codec = 0; codec < 2; codec++) {
        for (int b0 = 0; b0 < 124; b0++) {
            float f0 = codec ? exp2f(-4.322605f - (0.0220427f * (float)b0))
                             : exp2f(-4.24738f - (0.0217705f * ((float)b0 + 0.5f)));
            emit("pitch", f0 * (float)2 * REFERENCE_PI);
        }
    }
    for (int b0 = 0; b0 < 208; b0++) {
        float w0 = (float)(4 * REFERENCE_PI) / ((float)b0 + 39.5f);
        printf("harmonics %u\n", (unsigned)(unsigned char)(0.9254 * (int)((REFERENCE_PI / w0) + 0.25)));
    }
    for (int i = 0; i < 160; i++) {
        emit("hann", 0.5f - (0.5f * cosf(2.0f * (float)REFERENCE_PI * (float)(i + 1) / 161.0f)));
    }
    for (int J = 1; J <= 10; J++) {
        for (int k = 1; k <= J; k++) {
            for (int j = 1; j <= J; j++) {
                emit("imbe", cosf((float)REFERENCE_PI * (float)(k - 1) * ((float)j - 0.5f) / (float)J));
            }
        }
    }
    return ferror(stdout) ? 1 : 0;
}
