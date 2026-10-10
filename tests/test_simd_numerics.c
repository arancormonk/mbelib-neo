// SPDX-License-Identifier: GPL-2.0-or-later
/* Numerical oracles for private SIMD kernels, linked to the selected library. */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "mbe_tone_fit.h"
#include "mbe_voiced.h" // IWYU pragma: keep (MBE_VOICED_KERNELS selects the voiced checks)

#define PI 3.14159265358979323846

/* Independent Gaussian elimination with pivoting; production uses Cholesky. */
static void
solve_reference(double a[4][5], int n, double coef[4]) {
    assert(n > 0 && n <= 4);
    for (int k = 0; k < n; ++k) {
        int pivot = k;
        for (int j = k + 1; j < n; ++j) {
            if (fabs(a[j][k]) > fabs(a[pivot][k])) {
                pivot = j;
            }
        }
        assert(fabs(a[pivot][k]) > 1e-12);
        for (int j = 0; j <= n; ++j) {
            double swap = a[k][j];
            a[k][j] = a[pivot][j];
            a[pivot][j] = swap;
        }
        for (int i = k + 1; i < n; ++i) {
            double factor = a[i][k] / a[k][k];
            for (int j = k; j <= n; ++j) {
                a[i][j] -= factor * a[k][j];
            }
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = a[i][n];
        for (int j = i + 1; j < n; ++j) {
            sum -= a[i][j] * coef[j];
        }
        coef[i] = sum / a[i][i];
    }
}

static double
reference_tone(const float* samples, int n, const double hz[2], int count, double amp[2]) {
    assert(count == 1 || count == 2);
    double a[4][5] = {{0}}, rhs[4] = {0}, coef[4] = {0};
    int dim = 2 * count;
    for (int i = 0; i < n; ++i) {
        double v[4] = {0};
        for (int k = 0; k < count; ++k) {
            double phase = 2.0 * PI * hz[k] * i / 8000.0;
            v[2u * (size_t)k] = cos(phase);
            v[2u * (size_t)k + 1u] = sin(phase);
        }
        for (int r = 0; r < dim; ++r) {
            rhs[r] += v[r] * samples[i];
            for (int q = 0; q < dim; ++q) {
                a[r][q] += v[r] * v[q];
            }
        }
    }
    for (int r = 0; r < dim; ++r) {
        a[r][dim] = rhs[r];
    }
    solve_reference(a, dim, coef);
    double energy = 0;
    for (int r = 0; r < dim; ++r) {
        energy += rhs[r] * coef[r];
    }
    amp[0] = hypot(coef[0], coef[1]);
    if (count == 2) {
        amp[1] = hypot(coef[2], coef[3]);
    }
    return energy;
}

static void
test_tone_fit(void) {
    const double pairs[][2] = {{60.1, 75.2}, {350, 490}, {440, 480}, {697, 1209}, {1900, 2300}, {3910, 3949}};
    uint32_t rng = 0x135a19u;
    for (int count = 1; count <= 2; ++count) {
        for (unsigned p = 0; p < sizeof(pairs) / sizeof(pairs[0]); ++p) {
            for (int n = 72; n <= 160; ++n) {
                float samples[161];
                for (int i = 0; i <= n; ++i) {
                    rng = rng * 1664525u + 1013904223u;
                    samples[i] =
                        (float)(1000.0 * cos(2.0 * PI * pairs[p][0] * i / 8000.0 + 0.77)
                                + 500.0 * sin(2.0 * PI * pairs[p][1] * i / 8000.0) + 0.1 * (double)(rng >> 24));
                }
                double expected_amp[2], actual_amp[2];
                /* Offset tests unaligned loads and every possible scalar tail. */
                double expected = reference_tone(samples + 1, n, pairs[p], count, expected_amp);
                double actual = mbe_tone_fit(samples + 1, n, pairs[p], count, actual_amp);
                assert(fabs(actual - expected) <= 1e-9 * fmax(1.0, fabs(expected)));
                for (int k = 0; k < count; ++k) {
                    assert(fabs(actual_amp[k] - expected_amp[k]) <= 1e-9 * fmax(1.0, expected_amp[k]));
                }
            }
        }
    }
}

#if defined(MBE_VOICED_KERNELS)
static void
test_voiced_kernels(void) {
    float window[320];
    for (int n = 0; n < 320; ++n) {
        window[n] = n < 160 ? (float)n / 160.0f : (float)(320 - n) / 160.0f;
    }
    const float frequencies[] = {0.05f, 0.1f, 0.7f, 1.5f, 3.13f};
    for (int trial = 0; trial < 80; ++trial) {
        float out[160] = {0};
        struct mbe_voiced_component p = {(float)trial * 5.13f, frequencies[trial % 5], 1.1f, trial % 3 != 0};
        struct mbe_voiced_component c = {(float)trial * -2.41f, frequencies[(trial + 1) % 5], 0.7f, trial % 3 != 1};
        mbe_voiced_windowed(out, window, &p, &c);
        for (int n = 0; n < 160; ++n) {
            double expected = p.active ? p.gain * window[n + 160] * cos((double)p.phase + (double)p.omega * n) : 0;
            expected += c.active ? c.gain * window[n] * cos((double)c.phase + (double)c.omega * n) : 0;
            assert(fabs(out[n] - expected) <= 3e-5);
        }
        for (int harmonic = 1; harmonic < 8; ++harmonic) {
            float phase = (float)trial * 6.15f;
            float delta = trial & 1 ? 0.009f : -0.009f;
            if (trial % 3 == 0) {
                delta = 0.0f;
            }
            float frequency = 0.1f * (float)harmonic;
            for (int n = 0; n < 160; ++n) {
                out[n] = 0.25f;
            }
            mbe_voiced_interpolated(out, phase, 0.5f, 1.25f, frequency, delta, harmonic);
            for (int n = 0; n < 160; ++n) {
                /* The historical float phase expression, including its rounding. */
                float theta = phase + frequency * (float)n + delta * (float)(harmonic * n * n) / 320.0f;
                float amp = 0.5f + ((float)n / 160.0f) * 0.75f;
                double expected = 0.25 + 2.0 * amp * cos((double)theta);
                assert(fabs(out[n] - expected) <= 2e-4);
            }
        }
    }
}

/* Caller-supplied phases far beyond the decoder's range, up to near the
 * largest float. Each lane's seed rotates the phase's own cosine and sine by
 * its offset, so the output follows the exact phase; the oracle does the same
 * in double, since phase + offset would round the offset away. A fast-math
 * library build need not reduce huge cos/sin arguments exactly, so there
 * phases beyond LARGE_PHASE_EXACT only have to give a finite waveform within
 * its amplitude; any departure from the exact phase is still reported. */
#if defined(MBELIB_TEST_LIBRARY_FAST_MATH)
#define LARGE_PHASE_EXACT 65536.0
#else
#define LARGE_PHASE_EXACT HUGE_VAL
#endif

static void
test_large_phases(void) {
    const float phases[] = {-2048.0f, 511.9f, 512.1f, 65536.0f, 0x1p52f, -3.0e38f};
    for (unsigned k = 0; k < sizeof(phases) / sizeof(phases[0]); ++k) {
        for (int chirp = 0; chirp < 2; ++chirp) {
            const float delta = chirp ? 0.009f : 0.0f;
            const int harmonic = 3;
            const int exact = fabs((double)phases[k]) <= LARGE_PHASE_EXACT;
            float out[160] = {0};
            mbe_voiced_interpolated(out, phases[k], 0.5f, 1.0f, 0.125f, delta, harmonic);
            const double pc = cos((double)phases[k]), ps = sin((double)phases[k]);
            double worst = 0.0;
            int worst_n = 0;
            for (int n = 0; n < 160; ++n) {
                double offset = 0.125 * n + (double)delta * harmonic * n * n / 320.0;
                double amp = 0.5 + (double)n / 320.0;
                double expected = 2.0 * amp * (pc * cos(offset) - ps * sin(offset));
                double error = fabs(out[n] - expected);
                if (!(error <= worst)) {
                    worst = error;
                    worst_n = n;
                }
                /* Finite and within the waveform's amplitude: fabs(NaN) fails. */
                assert(fabs(out[n]) <= 2.0 * amp + 1e-3);
            }
            if (!(worst <= 2e-4)) {
                /* stderr is unbuffered, so the report survives the assertion's abort. */
                fprintf(stderr, "large phase %g, pitch delta %g: largest error %g at sample %d%s\n", (double)phases[k],
                        (double)delta, worst, worst_n, exact ? "" : " (fast-math library: not required)");
            }
            assert(!exact || worst <= 2e-4);
        }
    }
}
#endif

int
main(void) {
    test_tone_fit();
#if defined(MBE_VOICED_KERNELS)
    test_voiced_kernels();
    test_large_phases();
#endif
    puts("SIMD numerical oracles passed");
    return 0;
}
