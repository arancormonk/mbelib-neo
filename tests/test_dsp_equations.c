// SPDX-License-Identifier: GPL-2.0-or-later
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "dsp_equation_helpers.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

#include <limits.h>

struct frequency_bin {
    double real;
    double imag;
};

static void
direct_dft(const float* input, struct frequency_bin* spectrum) {
    for (int k = 0; k <= 128; ++k) {
        double real = 0.0;
        double imag = 0.0;
        for (int n = 0; n < 256; ++n) {
            const double angle = (2.0 * TEST_PI / 256.0) * (double)(k * n);
            real += input[n] * cos(angle);
            imag -= input[n] * sin(angle);
        }
        spectrum[k].real = real;
        spectrum[k].imag = imag;
    }
}

static void
test_dft(mbe_fft_plan* plan) {
    float storage[257];
    float ordered_storage[257];
    float* input = storage + 1;
    float* ordered = ordered_storage + 1;
    for (int n = 0; n < 256; ++n) {
        input[n] =
            (float)(0.75 + 1.2 * sin(2.0 * TEST_PI * 3.0 * n / 256.0) + 0.33 * cos(2.0 * TEST_PI * 37.0 * n / 256.0)
                    + (n % 2 == 0 ? 0.2 : -0.2) + (n % 17) / 31.0);
    }
    struct frequency_bin expected[129];
    direct_dft(input, expected);
    assert(mbe_fft_forward_real(plan, input, ordered) == 0);
    assert_close(ordered[0], expected[0].real, 0.0001, 2.0e-6);
    assert_close(ordered[1], expected[128].real, 0.0001, 2.0e-6);
    for (int k = 1; k < 128; ++k) {
        assert_close(ordered[2u * (size_t)k], expected[k].real, 0.0001, 2.0e-6);
        assert_close(ordered[2u * (size_t)k + 1u], expected[k].imag, 0.0001, 2.0e-6);
    }
    assert(mbe_fft_forward_real(plan, input, input) == 0);
    assert(compare_bytes(input, ordered, 256 * sizeof(float)) == 0);
    assert(mbe_fft_forward_real(NULL, input, ordered) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_fft_forward_real(plan, NULL, ordered) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_fft_forward_real(plan, input, NULL) == MBE_STATUS_INVALID_ARGUMENT);
}

static void
check_acf(mbe_acf_plan* plan, int len, int lags) {
    float x_storage[513];
    float r_storage[513];
    float* x = x_storage + 1;
    float* r = r_storage + 1;
    for (int n = 0; n < 512; ++n) {
        x[n] = (float)(0.01 * (n % 13) + sin(0.07 * n) - cos(0.17 * n));
    }
    assert(mbe_acf_compute(plan, x, len, r, lags) == 0);
    for (int lag = 0; lag <= lags; ++lag) {
        double expected = 0.0;
        for (int n = 0; n + lag < len; ++n) {
            expected += (double)x[n] * x[n + lag];
        }
        assert_close(r[lag], expected, 0.0001, 3.0e-6);
    }
    assert(mbe_acf_compute(plan, x, len, x, lags) == 0);
    assert(compare_bytes(x, r, (size_t)(lags + 1) * sizeof(float)) == 0);
}

static void
test_acf(void) {
    mbe_acf_plan* plan = mbe_acf_plan_alloc();
    assert(plan != NULL);
    check_acf(plan, 0, 0);
    check_acf(plan, 0, 511);
    check_acf(plan, 1, 511);
    check_acf(plan, 160, 100);
    check_acf(plan, 257, 255);
    check_acf(plan, 400, 112);
    check_acf(plan, 512, 0);
    float values[512] = {0};
    assert(mbe_acf_compute(NULL, values, 1, values, 1) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, NULL, 1, values, 1) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, values, 1, NULL, 1) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, values, -1, values, 1) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, values, 1, values, -1) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, values, 513, values, 0) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, values, 0, values, 512) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, values, 512, values, 1) == MBE_STATUS_INVALID_ARGUMENT);
    assert(mbe_acf_compute(plan, values, INT_MAX, values, INT_MAX) == MBE_STATUS_INVALID_ARGUMENT);
    mbe_acf_plan_free(plan);
    mbe_acf_plan_free(NULL);
}

static void
test_window(void) {
    for (int n = -300; n <= 300; ++n) {
        const float expected = (float)window_equation(n);
        assert(float_bits(mbe_synthesisWindow(n)) == float_bits(expected));
    }
    assert(float_bits(mbe_synthesisWindow(INT_MIN)) == 0);
    assert(float_bits(mbe_synthesisWindow(INT_MAX)) == 0);
}

static double
wola_equation(const float* prior, const double* current, int n, int length) {
    const double old_weight = window_equation(n);
    const double new_weight = window_equation(n - length);
    double numerator = 0.0;
    if (old_weight > 0.0) {
        numerator += old_weight * prior[n + 128];
    }
    if (new_weight > 0.0) {
        numerator += new_weight * current[n - length + 128];
    }
    const double denominator = old_weight * old_weight + new_weight * new_weight;
    return denominator > 0.0 ? numerator / denominator : 0.0;
}

static void
check_wola(int length) {
    float previous[256];
    float current[256];
    double reference_current[256];
    float output[256];
    for (int n = 0; n < 256; ++n) {
        previous[n] = (float)(sin(n * 0.17) * 23.0);
        current[n] = (float)(cos(n * 0.13) * 31.0);
        reference_current[n] = current[n];
        output[n] = 3.0f;
    }
    mbe_wola_combine(output, previous, current, length);
    for (int n = 0; n < length; ++n) {
        const double expected = 3.0 + wola_equation(previous, reference_current, n, length);
        assert_close(output[n], expected, 0.00002, 2.0e-6);
    }
    for (int n = length; n < 256; ++n) {
        assert(float_bits(output[n]) == float_bits(3.0f));
    }
    mbe_wola_combine(NULL, previous, current, length);
    mbe_wola_combine(output, NULL, current, length);
    mbe_wola_combine(output, previous, NULL, length);
    mbe_wola_combine(output, previous, current, -1);
    mbe_wola_combine(output, previous, current, 257);
}

static void
shape_reference(struct frequency_bin* spectrum, const mbe_parms* model) {
    struct frequency_bin shaped[129] = {{0}};
    const float spacing = MBE_256_OVER_2PI * model->w0;
    for (int l = 1; l <= model->L; ++l) {
        const int lower = (int)ceilf(((float)l - 0.5f) * spacing);
        const int upper = (int)ceilf(((float)l + 0.5f) * spacing);
        double sum = 0.0;
        for (int k = lower; k < upper; ++k) {
            sum += spectrum[k].real * spectrum[k].real + spectrum[k].imag * spectrum[k].imag;
        }
        if (model->Vl[l] != 0 || upper <= lower || sum <= 1.0e-10) {
            continue;
        }
        const double factor = 146.17696 * model->Ml[l] / sqrt(sum / (upper - lower));
        for (int k = lower; k < upper; ++k) {
            shaped[k].real = spectrum[k].real * factor;
            shaped[k].imag = spectrum[k].imag * factor;
        }
    }
    memcpy(spectrum, shaped, sizeof(shaped));
}

static void
reference_unvoiced(const float* noise, const mbe_parms* model, double* current) {
    float windowed[256];
    for (int n = 0; n < 256; ++n) {
        windowed[n] = (float)(noise[n] * window_equation(n - 128));
    }
    struct frequency_bin spectrum[129];
    direct_dft(windowed, spectrum);
    shape_reference(spectrum, model);
    for (int n = 0; n < 256; ++n) {
        double sum = spectrum[0].real;
        for (int k = 1; k < 128; ++k) {
            const double angle = (2.0 * TEST_PI / 256.0) * (double)(k * n);
            sum += 2.0 * (spectrum[k].real * cos(angle) - spectrum[k].imag * sin(angle));
        }
        current[n] = sum / 256.0;
    }
}

static void
test_unvoiced(mbe_fft_plan* plan) {
    mbe_parms current;
    new_model(&current);
    current.L = 20;
    current.w0 = 0.13f;
    for (int l = 1; l <= current.L; ++l) {
        current.Ml[l] = 10.0f + 0.2f * (float)l;
        current.Vl[l] = l % 3 == 1;
    }
    for (int n = 0; n < 256; ++n) {
        current.previousUw[n] = (float)(cos(n * 0.13) * 3.0);
    }
    const mbe_parms previous = current;
    float noise[256];
    uint32_t random = 3147;
    for (int n = 0; n < 256; ++n) {
        noise[n] = (float)random;
        random = (random * 171u + 11213u) % 53125u;
    }
    double expected_current[256];
    reference_unvoiced(noise, &current, expected_current);
    float output[160];
    for (int n = 0; n < 160; ++n) {
        output[n] = 1.5f;
    }
    mbe_synthesizeUnvoicedFFTWithNoise(output, &current, &current, plan, noise);
    for (int n = 0; n < 256; ++n) {
        assert_close(current.previousUw[n], expected_current[n], 0.00015, 5.0e-6);
    }
    for (int n = 0; n < 160; ++n) {
        assert_close(output[n], 1.5 + wola_equation(previous.previousUw, expected_current, n, 160), 0.0003, 6.0e-6);
    }
}

static void
test_empty_bands(mbe_fft_plan* plan) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    float output[160] = {0};
    float noise[256] = {0};
    mbe_synthesizeUnvoicedFFTWithNoise(output, &current, &previous, plan, noise);
    for (int n = 0; n < 256; ++n) {
        assert((float_bits(current.previousUw[n]) & UINT32_C(0x7fffffff)) == 0);
    }
    current.w0 = 0.0f;
    for (int n = 0; n < 256; ++n) {
        noise[n] = (float)n;
    }
    mbe_synthesizeUnvoicedFFTWithNoise(output, &current, &previous, plan, noise);
    for (int n = 0; n < 160; ++n) {
        assert((float_bits(output[n]) & UINT32_C(0x7fffffff)) == 0);
    }
    mbe_synthesizeUnvoicedFFT(NULL, &current, &previous, plan);
    mbe_synthesizeUnvoicedFFT(output, NULL, &previous, plan);
    mbe_synthesizeUnvoicedFFT(output, &current, NULL, plan);
    mbe_synthesizeUnvoicedFFT(output, &current, &previous, NULL);
}

static void
check_unvoiced_harmonics(mbe_fft_plan* plan, int current_count, int previous_count) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    current.L = current_count;
    previous.L = previous_count;
    current.w0 = 0.05f;
    current.noiseSeed = 3147.0f;
    current.previousUw[0] = 77.0f;
    mbe_parms generated = current;
    float noise[256] = {0};
    uint32_t state = 3147;
    for (int n = 96; n < 256; ++n) {
        noise[n] = (float)state;
        state = (171u * state + 11213u) % 53125u;
    }
    double expected[256];
    reference_unvoiced(noise, &current, expected);
    float output[160] = {0};
    float generated_output[160] = {0};
    mbe_synthesizeUnvoicedFFTWithNoise(output, &current, &previous, plan, noise);
    mbe_synthesizeUnvoicedFFT(generated_output, &generated, &previous, plan);
    assert(float_bits(generated.noiseSeed) == float_bits((float)state));
    for (int n = 0; n < 256; ++n) {
        assert_close(current.previousUw[n], expected[n], 0.00015, 5.0e-6);
        assert_close(generated.previousUw[n], expected[n], 0.00015, 5.0e-6);
    }
    for (int n = 0; n < 160; ++n) {
        const double sample = wola_equation(previous.previousUw, expected, n, 160);
        assert_close(output[n], sample, 0.0003, 6.0e-6);
        assert_close(generated_output[n], sample, 0.0003, 6.0e-6);
    }
}

static void
test_unvoiced_harmonic_interval(mbe_fft_plan* plan) {
    const int counts[] = {1, 8, 56};
    for (unsigned i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i) {
        for (unsigned j = 0; j < sizeof(counts) / sizeof(counts[0]); ++j) {
            check_unvoiced_harmonics(plan, counts[i], counts[j]);
        }
    }
}

static void
check_invalid_unvoiced_harmonics(mbe_fft_plan* plan, int current_count, int previous_count) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    current.L = current_count;
    previous.L = previous_count;
    current.previousUw[0] = 77.0f;
    const mbe_parms saved = current;
    float output[160];
    float noise[256] = {0};
    for (int n = 0; n < 160; ++n) {
        output[n] = 3.0f;
    }
    mbe_synthesizeUnvoicedFFTWithNoise(output, &current, &previous, plan, noise);
    mbe_synthesizeUnvoicedFFT(output, &current, &previous, plan);
    assert(compare_bytes(&current, &saved, sizeof(current)) == 0);
    for (int n = 0; n < 160; ++n) {
        assert(float_bits(output[n]) == float_bits(3.0f));
    }
}

static void
test_invalid_unvoiced_harmonics(mbe_fft_plan* plan) {
    check_invalid_unvoiced_harmonics(plan, 0, 1);
    check_invalid_unvoiced_harmonics(plan, 57, 1);
    check_invalid_unvoiced_harmonics(plan, 1, 0);
    check_invalid_unvoiced_harmonics(plan, 1, 57);
}

/* Runtime bit-pattern inputs also exercise guards in fast-math library builds. */
static void
test_unvoiced_pitch_guard(mbe_fft_plan* plan) {
    const uint32_t patterns[] = {0,
                                 UINT32_C(0x80000000),
                                 UINT32_C(0xbf800000),
                                 UINT32_C(0x7f800000),
                                 UINT32_C(0xff800000),
                                 UINT32_C(0x7fc00001)};
    for (size_t k = 0; k < sizeof(patterns) / sizeof(patterns[0]); ++k) {
        mbe_parms current;
        mbe_parms previous;
        new_model(&current);
        new_model(&previous);
        memcpy(&current.w0, &patterns[k], sizeof(current.w0));
        current.previousUw[0] = 77.0f;
        previous.previousUw[128] = 1.0f;
        const mbe_parms saved = current;
        float output[160] = {0};
        float noise[256] = {0};
        mbe_synthesizeUnvoicedFFTWithNoise(output, &current, &previous, plan, noise);
        if (k < 2) {
            assert_close(output[0], 1.0, 1.0e-7, 0.0);
            assert((float_bits(current.previousUw[0]) & UINT32_C(0x7fffffff)) == 0);
        } else {
            assert(compare_bytes(&current, &saved, sizeof(current)) == 0);
            assert((float_bits(output[0]) & UINT32_C(0x7fffffff)) == 0);
        }
    }
}

int
main(void) {
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    test_dft(plan);
    test_acf();
    test_window();
    check_wola(80);
    check_wola(160);
    check_wola(256);
    test_unvoiced(plan);
    test_empty_bands(plan);
    test_unvoiced_harmonic_interval(plan);
    test_invalid_unvoiced_harmonics(plan);
    test_unvoiced_pitch_guard(plan);
    mbe_fft_plan_free(plan);
    mbe_fft_plan_free(NULL);
    puts("independent DFT/ACF/WOLA/unvoiced equations: passed");
    return 0;
}
