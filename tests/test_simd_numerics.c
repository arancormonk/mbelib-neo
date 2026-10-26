// SPDX-License-Identifier: GPL-2.0-or-later
/* Numerical oracles for private SIMD kernels, linked to the selected library.
 * The kernels that a run-time AVX2 dispatch selects run under every kernel
 * ISA this build and CPU can select. */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "ambe_encoder.h"
#include "mbe_cpu.h"
#include "mbe_tone_fit.h"
#include "mbe_unvoiced_fft.h"
#include "mbe_voiced.h" // IWYU pragma: keep (MBE_VOICED_KERNELS selects the voiced checks)

#define PI 3.14159265358979323846

static void
assert_float_bits_equal(const float* actual, const float* expected, int count) {
    for (int i = 0; i < count; ++i) {
        uint32_t actual_bits, expected_bits;
        memcpy(&actual_bits, &actual[i], sizeof(actual_bits));
        memcpy(&expected_bits, &expected[i], sizeof(expected_bits));
        assert(actual_bits == expected_bits);
    }
}

static void
test_noise_sequence(void) {
    /* Every possible LCG state, with three frames to exercise both overlap
     * and returned state. The oracle advances one sample at a time, never
     * using the production batch constants. Offset buffers exercise callers
     * that do not provide SIMD alignment. */
    for (uint32_t initial = 0; initial < 53125u; ++initial) {
        float seed = (float)initial;
        uint32_t state = initial;
        float storage[MBE_FFT_SIZE + 1], overlap_storage[MBE_NOISE_OVERLAP + 1];
        float* actual = storage + 1;
        float* overlap = overlap_storage + 1;
        float expected[MBE_FFT_SIZE], expected_overlap[MBE_NOISE_OVERLAP];
        for (int i = 0; i < MBE_NOISE_OVERLAP; ++i) {
            overlap[i] = expected_overlap[i] = (float)(i * 37);
        }
        for (int frame = 0; frame < 3; ++frame) {
            memcpy(expected, expected_overlap, sizeof(expected_overlap));
            for (int i = MBE_NOISE_OVERLAP; i < MBE_FFT_SIZE; ++i) {
                expected[i] = (float)state;
                state = (171u * state + 11213u) % 53125u;
            }
            memcpy(expected_overlap, expected + MBE_FFT_SIZE - MBE_NOISE_OVERLAP, sizeof(expected_overlap));
            float expected_seed = (float)state;
            mbe_generate_noise_with_overlap(actual, &seed, overlap);
            assert_float_bits_equal(actual, expected, MBE_FFT_SIZE);
            assert_float_bits_equal(overlap, expected_overlap, MBE_NOISE_OVERLAP);
            assert_float_bits_equal(&seed, &expected_seed, 1);
        }
    }
}

static int
nearest_reference(const float* table, int width, int rows, int step, const float* target, int used) {
    int best = 0;
    double best_error = 1e30;
    for (int row = 0; row < rows; row += step) {
        double error = 0.0;
        for (int m = 0; m < used; ++m) {
            double d = (double)target[m] - table[row * width + m];
            error += d * d;
        }
        if (error < best_error) {
            best = row;
            best_error = error;
        }
    }
    return best;
}

static void
test_nearest_codebook(void) {
    float storage[33 * 4 + 1], target_storage[5];
    float* table = storage + 1;
    float* target = target_storage + 1;
    uint32_t rng = 0x793afb21u;
    for (int width = 3; width <= 4; ++width) {
        for (int i = 0; i < 33 * width; ++i) {
            rng = rng * 1664525u + 1013904223u;
            /* Small integers make float and double distances exact, even
             * under fast-math, so this oracle can check the chosen index. */
            table[i] = (float)((int)(rng >> 24) - 128);
        }
        for (int m = 0; m < width; ++m) {
            target[m] = table[17 * width + m];
            table[4 * width + m] = target[m]; /* tied minima in different lanes */
        }
        for (int rows = 1; rows <= 33; ++rows) {
            for (int step = 1; step <= 3; ++step) {
                for (int used = 0; used <= width; ++used) {
                    assert(ambe_enc_nearest(table, width, rows, step, target, used)
                           == nearest_reference(table, width, rows, step, target, used));
                }
            }
        }
    }
}

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
    float window_storage[321];
    float* window = window_storage + 1;
    for (int n = 0; n < 320; ++n) {
        window[n] = n < 160 ? (float)n / 160.0f : (float)(320 - n) / 160.0f;
    }
    const float frequencies[] = {0.05f, 0.1f, 0.7f, 1.5f, 3.13f};
    for (int trial = 0; trial < 80; ++trial) {
        float out_storage[161] = {0};
        float* out = out_storage + 1;
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
            /* Odd trials raise the moving value to the fourth power (DVSI's law). */
            const int fourth_power = trial & 1;
            mbe_voiced_interpolated(out, phase, 0.5f, 1.25f, frequency, delta, harmonic, fourth_power);
            for (int n = 0; n < 160; ++n) {
                /* The historical float phase expression, including its rounding. */
                float theta = phase + frequency * (float)n + delta * (float)(harmonic * n * n) / 320.0f;
                float amp = 0.5f + ((float)n / 160.0f) * 0.75f;
                if (fourth_power) {
                    amp = amp * amp * amp * amp;
                }
                double expected = 0.25 + 2.0 * amp * cos((double)theta);
                assert(fabs(out[n] - expected) <= 2e-4 * (fourth_power ? 2.5 : 1.0));
            }
        }
    }
}

/* Caller-supplied phases far beyond the decoder's range, up to near the
 * largest float. Each lane's seed rotates the phase's own cosine and sine by
 * its offset, so the output follows the exact phase; the oracle does the same
 * in double, since phase + offset would round the offset away. The library
 * compiles these kernels with IEEE semantics in fast-math builds too, so the
 * check holds there as well. */
static void
test_large_phases(void) {
    const float phases[] = {-2048.0f, 511.9f, 512.1f, 65536.0f, 0x1p52f, -3.0e38f};
    for (unsigned k = 0; k < sizeof(phases) / sizeof(phases[0]); ++k) {
        for (int chirp = 0; chirp < 2; ++chirp) {
            const float delta = chirp ? 0.009f : 0.0f;
            const int harmonic = 3;
            float out[160] = {0};
            mbe_voiced_interpolated(out, phases[k], 0.5f, 1.0f, 0.125f, delta, harmonic, 0);
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
                fprintf(stderr, "large phase %g, pitch delta %g: largest error %g at sample %d\n", (double)phases[k],
                        (double)delta, worst, worst_n);
            }
            assert(worst <= 2e-4);
        }
    }
}
#endif

static void
test_dispatched_kernels(void) {
    test_tone_fit();
#if defined(MBE_VOICED_KERNELS)
    test_voiced_kernels();
    test_large_phases();
#endif
}

int
main(void) {
    test_noise_sequence();
    test_nearest_codebook();
    test_dispatched_kernels();
    const int forced[] = {MBE_CPU_ISA_SSE2, MBE_CPU_ISA_AVX2};
    for (unsigned i = 0; i < sizeof(forced) / sizeof(forced[0]); ++i) {
        if (mbe_cpu_select_isa(forced[i]) == 0) {
            printf("kernel ISA %d\n", mbe_cpu_isa());
            test_dispatched_kernels();
        }
    }
    assert(mbe_cpu_select_isa(MBE_CPU_ISA_AUTO) == 0);
    puts("SIMD numerical oracles passed");
    return 0;
}
