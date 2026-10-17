// SPDX-License-Identifier: GPL-2.0-or-later
/* Independently implemented from TIA-102.BABA equations 117--126 and the policies in legal/implementation-provenance.md. */
#include "mbe_unvoiced_fft.h"
#include "mbe_compiler.h"

#include "pffft.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct transform_storage {
    PFFFT_Setup* setup;
    float* allocation;
    float* time;
    float* spectrum;
    float* shaped;
    float* work;
};

struct mbe_fft_plan {
    struct transform_storage fft;
};

struct mbe_acf_plan {
    struct transform_storage fft;
};

static MBE_THREAD_LOCAL uint32_t pending_noise_seed;
static MBE_THREAD_LOCAL int noise_seed_pending;

static void
release_transform(struct transform_storage* fft) {
    if (fft->setup != NULL) {
        pffft_destroy_setup(fft->setup);
    }
    if (fft->allocation != NULL) {
        pffft_aligned_free(fft->allocation);
    }
}

static int
prepare_transform(struct transform_storage* fft, int size) {
    fft->setup = pffft_new_setup(size, PFFFT_REAL);
    if (fft->setup == NULL) {
        return 0;
    }
    fft->allocation = pffft_aligned_malloc(4u * (size_t)size * sizeof(float));
    if (fft->allocation == NULL) {
        return 0;
    }
    fft->time = fft->allocation;
    fft->spectrum = fft->time + size;
    fft->shaped = fft->spectrum + size;
    fft->work = fft->shaped + size;
    return 1;
}

mbe_fft_plan*
mbe_fft_plan_alloc(void) {
    mbe_fft_plan* plan = calloc(1, sizeof(*plan));
    if (plan != NULL && !prepare_transform(&plan->fft, MBE_FFT_SIZE)) {
        mbe_fft_plan_free(plan);
        return NULL;
    }
    return plan;
}

void
mbe_fft_plan_free(mbe_fft_plan* plan) {
    if (plan != NULL) {
        release_transform(&plan->fft);
        free(plan);
    }
}

mbe_acf_plan*
mbe_acf_plan_alloc(void) {
    mbe_acf_plan* plan = calloc(1, sizeof(*plan));
    if (plan != NULL && !prepare_transform(&plan->fft, MBE_ACF_SIZE)) {
        mbe_acf_plan_free(plan);
        return NULL;
    }
    return plan;
}

void
mbe_acf_plan_free(mbe_acf_plan* plan) {
    if (plan != NULL) {
        release_transform(&plan->fft);
        free(plan);
    }
}

int
mbe_fft_forward_real(mbe_fft_plan* plan, const float input[MBE_FFT_SIZE], float output[MBE_FFT_SIZE]) {
    if (plan == NULL || input == NULL || output == NULL) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    struct transform_storage* fft = &plan->fft;
    memcpy(fft->time, input, MBE_FFT_SIZE * sizeof(float));
    pffft_transform_ordered(fft->setup, fft->time, fft->spectrum, fft->work, PFFFT_FORWARD);
    memcpy(output, fft->spectrum, MBE_FFT_SIZE * sizeof(float));
    return 0;
}

int
mbe_acf_compute(mbe_acf_plan* plan, const float* x, int len, float* r, int lags) {
    if (plan == NULL || x == NULL || r == NULL || len < 0 || lags < 0 || len > MBE_ACF_SIZE || lags >= MBE_ACF_SIZE
        || len > MBE_ACF_SIZE - lags) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    struct transform_storage* fft = &plan->fft;
    memset(fft->time, 0, MBE_ACF_SIZE * sizeof(float));
    memcpy(fft->time, x, (size_t)len * sizeof(float));
    pffft_transform_ordered(fft->setup, fft->time, fft->spectrum, fft->work, PFFFT_FORWARD);
    fft->spectrum[0] *= fft->spectrum[0];
    fft->spectrum[1] *= fft->spectrum[1];
    for (int k = 1; k < MBE_ACF_SIZE / 2; ++k) {
        const float real = fft->spectrum[2u * (size_t)k];
        const float imag = fft->spectrum[2u * (size_t)k + 1u];
        fft->spectrum[2u * (size_t)k] = real * real + imag * imag;
        fft->spectrum[2u * (size_t)k + 1u] = 0.0f;
    }
    pffft_transform_ordered(fft->setup, fft->spectrum, fft->time, fft->work, PFFFT_BACKWARD);
    for (int lag = 0; lag <= lags; ++lag) {
        r[lag] = fft->time[lag] * (1.0f / MBE_ACF_SIZE);
    }
    return 0;
}

static uint32_t
usable_noise_seed(float seed) {
    uint32_t bits;
    memcpy(&bits, &seed, sizeof(bits));
    if (bits == UINT32_C(0x80000000)) {
        return 0;
    }
    /* Largest binary32 value below 2^32; reject before integer conversion. */
    if (bits > UINT32_C(0x4f7fffff)) {
        return (uint32_t)MBE_LCG_DEFAULT_SEED;
    }
    return (uint32_t)seed % 53125u;
}

void
mbe_generate_noise_lcg(float* buffer, int count, float* seed) {
    if (buffer == NULL || seed == NULL || count <= 0) {
        return;
    }
    uint32_t state = usable_noise_seed(*seed);
    for (int i = 0; i < count; ++i) {
        buffer[i] = (float)state;
        state = (171u * state + 11213u) % 53125u;
    }
    *seed = (float)state;
}

void
mbe_seedUnvoicedNoiseLcg(uint32_t seed) {
    if (seed == 0) {
        seed = UINT32_C(0x6d25357b);
    }
    pending_noise_seed = seed % 53125u;
    noise_seed_pending = 1;
}

void
mbe_generate_noise_with_overlap(float* buffer, float* seed, float* overlap) {
    if (buffer == NULL || seed == NULL || overlap == NULL) {
        return;
    }
    if (*seed < 0.0f) {
        memset(buffer, 0, MBE_FFT_SIZE * sizeof(float));
        memset(overlap, 0, MBE_NOISE_OVERLAP * sizeof(float));
        *seed = noise_seed_pending ? (float)pending_noise_seed : MBE_LCG_DEFAULT_SEED;
        noise_seed_pending = 0;
        return;
    }
    memmove(buffer, overlap, MBE_NOISE_OVERLAP * sizeof(float));
    mbe_generate_noise_lcg(buffer + MBE_NOISE_OVERLAP, MBE_FFT_SIZE - MBE_NOISE_OVERLAP, seed);
    memmove(overlap, buffer + MBE_FFT_SIZE - MBE_NOISE_OVERLAP, MBE_NOISE_OVERLAP * sizeof(float));
}

float
mbe_synthesisWindow(int n) {
    /* Annex I's taper samples, rounded separately to binary32. */
    static const float taper[51] = {0.00f, 0.02f, 0.04f, 0.06f, 0.08f, 0.10f, 0.12f, 0.14f, 0.16f, 0.18f, 0.20f,
                                    0.22f, 0.24f, 0.26f, 0.28f, 0.30f, 0.32f, 0.34f, 0.36f, 0.38f, 0.40f, 0.42f,
                                    0.44f, 0.46f, 0.48f, 0.50f, 0.52f, 0.54f, 0.56f, 0.58f, 0.60f, 0.62f, 0.64f,
                                    0.66f, 0.68f, 0.70f, 0.72f, 0.74f, 0.76f, 0.78f, 0.80f, 0.82f, 0.84f, 0.86f,
                                    0.88f, 0.90f, 0.92f, 0.94f, 0.96f, 0.98f, 1.00f};
    if (n <= -105 || n >= 105) {
        return 0.0f;
    }
    if (n < 0) {
        n = -n;
    }
    if (n <= 55) {
        return 1.0f;
    }
    return taper[105 - n];
}

void
mbe_wola_combine(float* output, const float* prevUw, const float* currUw, int N) {
    if (output == NULL || prevUw == NULL || currUw == NULL || N <= 0 || N > MBE_FFT_SIZE) {
        return;
    }
    for (int n = 0; n < N; ++n) {
        const float previous_window = mbe_synthesisWindow(n);
        const float current_window = mbe_synthesisWindow(n - N);
        const float weight = previous_window * previous_window + current_window * current_window;
        float sample = 0.0f;
        if (previous_window > 0.0f) {
            sample += previous_window * prevUw[128 + n];
        }
        if (current_window > 0.0f) {
            sample += current_window * currUw[128 + n - N];
        }
        if (weight > 0.0f) {
            output[n] += sample / weight;
        }
    }
}

static int
valid_unvoiced_harmonics(const mbe_parms* mp) {
    return mp != NULL && mp->L >= 1 && mp->L <= 56;
}

static int
valid_unvoiced_model(const mbe_parms* mp) {
    if (!valid_unvoiced_harmonics(mp)) {
        return 0;
    }
    uint32_t bits;
    memcpy(&bits, &mp->w0, sizeof(bits));
    return bits <= UINT32_C(0x7f7fffff) || bits == UINT32_C(0x80000000);
}

static int
band_edge(float frequency) {
    if (frequency >= 128.0f) {
        return 128;
    }
    return (int)ceilf(frequency);
}

static float
band_power(const float* spectrum, int lower, int upper) {
    float power = 0.0f;
    for (int k = lower; k < upper; ++k) {
        const float real = spectrum[2u * (size_t)k];
        const float imag = k == 0 ? 0.0f : spectrum[2u * (size_t)k + 1u];
        power += real * real + imag * imag;
    }
    return power;
}

static void
shape_unvoiced(struct transform_storage* fft, const mbe_parms* mp) {
    memset(fft->shaped, 0, MBE_FFT_SIZE * sizeof(float));
    const float spacing = MBE_256_OVER_2PI * mp->w0;
    for (int l = 1; l <= mp->L; ++l) {
        if (mp->Vl[l] != 0) {
            continue;
        }
        const int lower = band_edge(((float)l - 0.5f) * spacing);
        const int upper = band_edge(((float)l + 0.5f) * spacing);
        const float power = band_power(fft->spectrum, lower, upper);
        if (upper <= lower || power <= 1.0e-10f) {
            continue;
        }
        const float scale = (MBE_UNVOICED_SCALE_COEFF * mp->Ml[l]) / sqrtf(power / (float)(upper - lower));
        for (int k = lower; k < upper; ++k) {
            fft->shaped[2u * (size_t)k] = fft->spectrum[2u * (size_t)k] * scale;
            if (k != 0) {
                fft->shaped[2u * (size_t)k + 1u] = fft->spectrum[2u * (size_t)k + 1u] * scale;
            }
        }
    }
}

void
mbe_synthesizeUnvoicedFFTWithNoise(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                                   const float* noise_buffer) {
    if (output == NULL || !valid_unvoiced_model(cur_mp) || !valid_unvoiced_harmonics(prev_mp) || plan == NULL
        || noise_buffer == NULL) {
        return;
    }
    struct transform_storage* fft = &plan->fft;
    for (int i = 0; i < MBE_FFT_SIZE; ++i) {
        fft->time[i] = noise_buffer[i] * mbe_synthesisWindow(i - MBE_FFT_SIZE / 2);
    }
    pffft_transform_ordered(fft->setup, fft->time, fft->spectrum, fft->work, PFFFT_FORWARD);
    shape_unvoiced(fft, cur_mp);
    pffft_transform_ordered(fft->setup, fft->shaped, fft->time, fft->work, PFFFT_BACKWARD);
    for (int i = 0; i < MBE_FFT_SIZE; ++i) {
        fft->time[i] *= 1.0f / MBE_FFT_SIZE;
    }
    mbe_wola_combine(output, prev_mp->previousUw, fft->time, 160);
    memcpy(cur_mp->previousUw, fft->time, sizeof(cur_mp->previousUw));
}

void
mbe_synthesizeUnvoicedFFT(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan) {
    if (output == NULL || !valid_unvoiced_model(cur_mp) || !valid_unvoiced_harmonics(prev_mp) || plan == NULL) {
        return;
    }
    float noise[MBE_FFT_SIZE];
    mbe_generate_noise_with_overlap(noise, &cur_mp->noiseSeed, cur_mp->noiseOverlap);
    mbe_synthesizeUnvoicedFFTWithNoise(output, cur_mp, prev_mp, plan, noise);
}
