// SPDX-License-Identifier: GPL-2.0-or-later
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "dsp_equation_helpers.h"
#include "mbe_adaptive.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

#include <float.h>
#include <limits.h>

static void
test_lcg(void) {
    float samples[320];
    float seed = 3147.0f;
    uint32_t state = 3147;
    mbe_generate_noise_lcg(samples, 320, &seed);
    for (int n = 0; n < 320; ++n) {
        assert(float_bits(samples[n]) == float_bits((float)state));
        state = (171u * state + 11213u) % 53125u;
    }
    assert(float_bits(seed) == float_bits((float)state));
    const float saved = seed;
    mbe_generate_noise_lcg(NULL, 2, &seed);
    mbe_generate_noise_lcg(samples, -1, &seed);
    mbe_generate_noise_lcg(samples, 2, NULL);
    assert(float_bits(seed) == float_bits(saved));
}

static void
test_seed_reduction(void) {
    const uint32_t seeds[] = {0, 53125, 53126, 10000000, UINT32_C(4294967040)};
    for (unsigned i = 0; i < sizeof(seeds) / sizeof(seeds[0]); ++i) {
        float seed = (float)seeds[i];
        uint32_t state = (uint32_t)seed % 53125u;
        float samples[2];
        mbe_generate_noise_lcg(samples, 2, &seed);
        for (int n = 0; n < 2; ++n) {
            assert(float_bits(samples[n]) == float_bits((float)state));
            state = (171u * state + 11213u) % 53125u;
        }
        assert(float_bits(seed) == float_bits((float)state));
    }
}

static void
test_overlap(void) {
    float seed = -1.0f;
    float overlap[96];
    float first[256];
    float second[256];
    memset(overlap, 0x55, sizeof(overlap));
    mbe_generate_noise_with_overlap(first, &seed, overlap);
    for (int n = 0; n < 256; ++n) {
        assert(float_bits(first[n]) == 0);
    }
    assert(float_bits(seed) == float_bits(3147.0f));
    mbe_generate_noise_with_overlap(first, &seed, overlap);
    uint32_t state = 3147;
    for (int n = 0; n < 96; ++n) {
        assert(float_bits(first[n]) == 0);
    }
    for (int n = 96; n < 256; ++n) {
        assert(float_bits(first[n]) == float_bits((float)state));
        state = (171u * state + 11213u) % 53125u;
    }
    mbe_generate_noise_with_overlap(second, &seed, overlap);
    assert(compare_bytes(first + 160, second, 96 * sizeof(float)) == 0);
    assert(compare_bytes(second + 160, overlap, sizeof(overlap)) == 0);
    for (int n = 96; n < 256; ++n) {
        assert(float_bits(second[n]) == float_bits((float)state));
        state = (171u * state + 11213u) % 53125u;
    }
    assert(float_bits(seed) == float_bits((float)state));
}

static void
test_seed_override(void) {
    float samples[256];
    float overlap[96] = {0};
    float seed = 123.0f;
    mbe_seedUnvoicedNoiseLcg(0);
    mbe_generate_noise_with_overlap(samples, &seed, overlap);
    seed = -1.0f;
    mbe_generate_noise_with_overlap(samples, &seed, overlap);
    assert(float_bits(seed) == float_bits((float)(UINT32_C(0x6d25357b) % 53125u)));
    seed = -1.0f;
    mbe_generate_noise_with_overlap(samples, &seed, overlap);
    assert(float_bits(seed) == float_bits(3147.0f));
    mbe_seedUnvoicedNoiseLcg(UINT32_MAX);
    seed = -1.0f;
    mbe_generate_noise_with_overlap(samples, &seed, overlap);
    assert(float_bits(seed) == float_bits((float)(UINT32_MAX % 53125u)));
}

static void
check_comfort(uint32_t seed, float amplitude, int initialize) {
    float samples[160];
    if (initialize) {
        mbe_seedComfortNoiseRng(seed);
    }
    const uint32_t mapped = seed == 0 ? UINT32_C(0x6d25357b) : seed;
    uint64_t state = (uint64_t)mapped ^ UINT64_C(0x5deece66d);
    mbe_synthesizeUniformNoisef(samples, amplitude);
    for (int n = 0; n < 160; ++n) {
        state = (state * UINT64_C(0x5deece66d) + 11) % (UINT64_C(1) << 48);
        const double expected = (((double)(state >> 24) / 16777216.0) * 2.0 - 1.0) * amplitude;
        assert_close(samples[n], expected, 2.0e-6, 1.0e-7);
    }
}

static void
test_comfort(void) {
    check_comfort(UINT32_C(0x12345678), 1.0f, 0);
    check_comfort(0, 3.0f, 1);
    check_comfort(UINT32_MAX, -0.5f, 1);
    float floats[160];
    short expected[160];
    short actual[160];
    mbe_seedComfortNoiseRng(37);
    mbe_synthesizeComfortNoisef(floats);
    mbe_floattoshort(floats, expected);
    mbe_seedComfortNoiseRng(37);
    mbe_synthesizeComfortNoise(actual);
    assert(compare_bytes(actual, expected, sizeof(actual)) == 0);
    mbe_synthesizeUniformNoisef(NULL, 1.0f);
    mbe_synthesizeComfortNoise(NULL);
    mbe_synthesizeComfortNoisef(NULL);
}

static double
reference_threshold(const mbe_parms* current, double energy) {
    if (current->errorRate <= 0.005f && current->errorCountTotal <= 4) {
        return DBL_MAX;
    }
    const double power = pow(energy, 3.0 / 8.0);
    if (current->errorRate <= 0.0125f && current->errorCount4 == 0) {
        return 45.255 * power / exp(277.26 * current->errorRate);
    }
    return 1.414 * power;
}

static void
check_smoothing(float rate, int errors, int block4) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    current.L = 9;
    current.errorRate = rate;
    current.errorCountTotal = errors;
    current.errorCount4 = block4;
    previous.localEnergy = 80000.0f;
    previous.amplitudeThreshold = 3000;
    const double energy = 0.95 * 80000.0 + 0.05 * 35000.0;
    const double threshold = reference_threshold(&current, energy);
    for (int l = 1; l <= 9; ++l) {
        current.Ml[l] = 130.0f * (float)l;
    }
    const int amplitude = rate <= 0.005f && errors <= 6 ? 20480 : 9000 - 300 * errors;
    const double sum = 130.0 * 45.0;
    const double scale = amplitude > sum ? 1.0 : fmax(0.0, (double)amplitude / sum);
    mbe_applyAdaptiveSmoothingWithRm0(&current, &previous, 35000.0f);
    assert_close(current.localEnergy, energy, 0.01, 0.0);
    assert(current.amplitudeThreshold == amplitude);
    for (int l = 1; l <= 9; ++l) {
        assert(current.Vl[l] == (130.0 * l > threshold));
        assert_close(current.Ml[l], 130.0 * l * scale, 0.0002, 1.0e-6);
    }
}

static void
test_smoothing_branches(void) {
    const float rates[] = {0.0f, 0.005f, 0.0125f, 0.02f};
    const int errors[] = {4, 5, 6, 7, 40};
    for (unsigned r = 0; r < sizeof(rates) / sizeof(rates[0]); ++r) {
        for (unsigned e = 0; e < sizeof(errors) / sizeof(errors[0]); ++e) {
            check_smoothing(rates[r], errors[e], 0);
            check_smoothing(rates[r], errors[e], 1);
        }
    }
    check_smoothing(nextafterf(0.005f, 1.0f), 4, 0);
    check_smoothing(nextafterf(0.0125f, 1.0f), 4, 0);
}

static void
test_threshold_boundaries(void) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    previous.localEnergy = 65536.0f;
    current.errorRate = 0.02f;
    const float boundary = 1.414f * 64.0f;
    current.Ml[1] = nextafterf(boundary, 0.0f);
    current.Ml[2] = boundary;
    current.Ml[3] = nextafterf(boundary, FLT_MAX);
    mbe_applyAdaptiveSmoothingWithRm0(&current, &previous, 65536.0f);
    assert(current.Vl[1] == 0);
    assert(current.Vl[2] == 0);
    assert(current.Vl[3] == 1);

    /* Equation 112's other branch, at error rate 0 (expf(0) is 1, so the
     * boundary is exact on every target): five errors, none in block 4. */
    mbe_parms gained;
    new_model(&gained);
    gained.errorRate = 0.0f;
    gained.errorCountTotal = 5;
    gained.errorCount4 = 0;
    const float gained_boundary = 45.255f * 64.0f;
    gained.Ml[1] = nextafterf(gained_boundary, 0.0f);
    gained.Ml[2] = gained_boundary;
    gained.Ml[3] = nextafterf(gained_boundary, FLT_MAX);
    mbe_applyAdaptiveSmoothingWithRm0(&gained, &previous, 65536.0f);
    assert(gained.Vl[1] == 0);
    assert(gained.Vl[2] == 0);
    assert(gained.Vl[3] == 1);

    current.errorRate = 0.0125f;
    current.errorCountTotal = 4;
    assert(!mbe_requiresAdaptiveSmoothing(&current));
    current.errorCountTotal = 5;
    assert(mbe_requiresAdaptiveSmoothing(&current));
    current.errorCountTotal = 4;
    current.errorRate = nextafterf(0.0125f, 1.0f);
    assert(mbe_requiresAdaptiveSmoothing(&current));
    current.errorRate = current.mutingThreshold;
    assert(!mbe_requiresMuting(&current));
    current.errorRate = nextafterf(current.mutingThreshold, 1.0f);
    assert(mbe_requiresMuting(&current));
    current.repeatCount = MBE_MAX_FRAME_REPEATS - 1;
    assert(!mbe_isMaxFrameRepeat(&current));
    ++current.repeatCount;
    assert(mbe_isMaxFrameRepeat(&current));
}

static void
test_damaged_and_alias_state(void) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    previous.localEnergy = 0.0f;
    previous.amplitudeThreshold = 0;
    current.errorRate = 0.03f;
    current.errorCountTotal = 7;
    mbe_applyAdaptiveSmoothingWithRm0(&current, &previous, 0.0f);
    assert_close(current.localEnergy, 71250.0, 0.01, 0.0);
    assert(current.amplitudeThreshold == 20480 + 6000 - 2100);
    previous.localEnergy = 10000.0f;
    current.errorCountTotal = INT_MAX;
    mbe_applyAdaptiveSmoothingWithRm0(&current, &previous, 0.0f);
    assert_close(current.localEnergy, 10000.0, 0.001, 0.0);
    for (int l = 1; l <= current.L; ++l) {
        assert(float_bits(current.Ml[l]) == 0);
    }
    new_model(&current);
    previous = current;
    mbe_parms expected = current;
    mbe_applyAdaptiveSmoothing(&expected, &previous);
    mbe_applyAdaptiveSmoothing(&current, &current);
    assert(compare_bytes(&current, &expected, sizeof(current)) == 0);
}

static void
check_harmonic_smoothing(int current_count, int previous_count) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    current.L = current_count;
    previous.L = previous_count;
    for (int l = 1; l <= current.L; ++l) {
        current.Ml[l] = 2.0f;
    }
    const mbe_parms initial = current;
    mbe_applyAdaptiveSmoothing(&current, &previous);
    assert_close(current.localEnergy, 0.95 * 75000.0 + 0.05 * 4.0 * current_count, 0.01, 0.0);
    current = initial;
    mbe_applyAdaptiveSmoothingWithRm0(&current, &previous, 20000.0f);
    assert_close(current.localEnergy, 72250.0, 0.01, 0.0);
}

static void
test_harmonic_interval(void) {
    const int counts[] = {1, 8, 56};
    for (unsigned i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i) {
        for (unsigned j = 0; j < sizeof(counts) / sizeof(counts[0]); ++j) {
            check_harmonic_smoothing(counts[i], counts[j]);
        }
    }
}

static void
check_invalid_smoothing(int current_count, int previous_count) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    current.L = current_count;
    previous.L = previous_count;
    const mbe_parms saved = current;
    mbe_applyAdaptiveSmoothing(&current, &previous);
    mbe_applyAdaptiveSmoothingWithRm0(&current, &previous, 5.0f);
    assert(compare_bytes(&current, &saved, sizeof(current)) == 0);
}

static void
test_invalid_smoothing(void) {
    mbe_parms current;
    mbe_parms previous;
    new_model(&current);
    new_model(&previous);
    const mbe_parms saved = current;
    mbe_applyAdaptiveSmoothing(NULL, &previous);
    mbe_applyAdaptiveSmoothing(&current, NULL);
    assert(compare_bytes(&current, &saved, sizeof(current)) == 0);
    check_invalid_smoothing(0, 1);
    check_invalid_smoothing(57, 1);
    check_invalid_smoothing(1, 0);
    check_invalid_smoothing(1, 57);
    assert(!mbe_requiresAdaptiveSmoothing(NULL));
    assert(!mbe_requiresMuting(NULL));
    assert(!mbe_isMaxFrameRepeat(NULL));
}

int
main(void) {
    test_lcg();
    test_seed_reduction();
    test_overlap();
    test_seed_override();
    test_comfort();
    test_smoothing_branches();
    test_threshold_boundaries();
    test_damaged_and_alias_state();
    test_harmonic_interval();
    test_invalid_smoothing();
    puts("independent adaptive/noise equations: passed");
    return 0;
}
