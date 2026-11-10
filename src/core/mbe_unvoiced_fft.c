// SPDX-License-Identifier: GPL-2.0-or-later
/* Independently implemented from TIA-102.BABA equations 117--126 and the policies in legal/implementation-provenance.md. */
#include "mbe_unvoiced_fft.h"
#include "mbe_adaptive.h"
#include "mbe_compiler.h"
#include "mbe_voiced.h"

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
    /* DVSI's unvoiced subframes: the square-root triangular window and the
     * output high-pass at each bin (real, imaginary). */
    float subframe_window[MBE_SUBFRAME_WINDOW];
    float output_filter[2 * (MBE_FFT_SIZE / 2 + 1)];
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

/* The subframe window, sqrt of a triangle over two subframes, and the output
 * high-pass the synthesis renders (mbe_low_band_gain(), mbe_low_band_phase())
 * at the transform's bins. */
static void
prepare_subframe_tables(mbe_fft_plan* plan) {
    for (int n = 0; n < MBE_SUBFRAME_WINDOW; ++n) {
        const int distance = (n < MBE_SUBFRAME) ? (MBE_SUBFRAME - n) : (n - MBE_SUBFRAME);
        plan->subframe_window[n] = sqrtf((float)(MBE_SUBFRAME - distance) / (float)MBE_SUBFRAME);
    }
    for (int k = 0; k <= MBE_FFT_SIZE / 2; ++k) {
        const float hz = (float)k * (8000.0f / (float)MBE_FFT_SIZE);
        const float gain = mbe_low_band_gain(hz);
        const double phase = (double)mbe_low_band_phase(hz);
        plan->output_filter[2u * (size_t)k] = gain * (float)cos(phase);
        plan->output_filter[2u * (size_t)k + 1u] = gain * (float)sin(phase);
    }
}

mbe_fft_plan*
mbe_fft_plan_alloc(void) {
    mbe_fft_plan* plan = calloc(1, sizeof(*plan));
    if (plan != NULL && !prepare_transform(&plan->fft, MBE_FFT_SIZE)) {
        mbe_fft_plan_free(plan);
        return NULL;
    }
    if (plan != NULL) {
        prepare_subframe_tables(plan);
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

/* MBE_SUBFRAME_STATE_MARK (-1.0f, 0xBF800000), by its bits. */
static int
subframe_mark_is_set(const float* slot) {
    uint32_t bits;
    memcpy(&bits, slot, sizeof(bits));
    return bits == UINT32_C(0xBF800000);
}

void
mbe_subframe_state_mark(mbe_parms* mp) {
    if (mp != NULL) {
        mp->noiseOverlap[MBE_SUBFRAME_STATE_SLOT] = MBE_SUBFRAME_STATE_MARK;
    }
}

int
mbe_subframe_state_is_marked(const mbe_parms* mp) {
    return mp != NULL && subframe_mark_is_set(&mp->noiseOverlap[MBE_SUBFRAME_STATE_SLOT]);
}

void
mbe_subframe_state_adopt(mbe_parms* mp) {
    if (mp == NULL || mbe_subframe_state_is_marked(mp)) {
        return;
    }
    /* The other rendering's noise overlap is no pulse train. */
    memset(mp->noiseOverlap, 0, sizeof(mp->noiseOverlap));
    mbe_subframe_state_mark(mp);
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
    /* What the subframe rendering left there is its pulse train, not noise. */
    if (subframe_mark_is_set(&overlap[MBE_SUBFRAME_STATE_SLOT])) {
        memset(overlap, 0, MBE_NOISE_OVERLAP * sizeof(float));
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

uint32_t
mbe_dvsi_noise_next(uint32_t state) {
    return (MBE_DVSI_NOISE_MULTIPLIER * (state & 0xFFFFu) + MBE_DVSI_NOISE_INCREMENT) & 0xFFFFu;
}

uint32_t
mbe_dvsi_noise_skip(uint32_t state, int samples) {
    state &= 0xFFFFu;
    for (int i = 0; i < samples; ++i) {
        state = mbe_dvsi_noise_next(state);
    }
    return state;
}

/* The generator state a stream's noiseSeed holds. A negative value is a fresh
 * stream, which starts from 0 as every DVSI stream does (the thread seed of
 * mbe_setThreadRngSeed() is not used: DVSI's generator has none); anything
 * else that is not a 16-bit state restarts it likewise. The bit pattern is
 * read so that fast-math cannot fold the test. */
static uint32_t
dvsi_noise_state(const float* seed) {
    uint32_t bits;
    memcpy(&bits, seed, sizeof(bits));
    /* 65535.0f is 0x477FFF00: the largest state a float holds exactly. */
    return (bits <= UINT32_C(0x477FFF00)) ? (uint32_t)*seed : 0u;
}

static void
clear_subframe_tail(mbe_parms* mp) {
    memset(mp->previousUw, 0, sizeof(mp->previousUw));
}

void
mbe_dvsi_noise_skip_frame(mbe_parms* mp) {
    if (mp == NULL) {
        return;
    }
    mp->noiseSeed = (float)mbe_dvsi_noise_skip(dvsi_noise_state(&mp->noiseSeed), MBE_VOICED_FRAME);
    clear_subframe_tail(mp);
    mbe_subframe_state_adopt(mp);
}

static int
valid_subframe(const struct mbe_unvoiced_subframe* sub) {
    uint32_t bits;
    if (sub == NULL || sub->L < 1 || sub->L > 56) {
        return 0;
    }
    /* A positive, finite fundamental that is not subnormal (zero to the
     * arithmetic where the caller's environment flushes denormals). */
    memcpy(&bits, &sub->w0, sizeof(bits));
    return bits - UINT32_C(0x00800000) < UINT32_C(0x7F000000);
}

static int
subframe_band_edge(float frequency) {
    const int top = MBE_FFT_SIZE / 2;
    const float rounded = frequency + 0.5f;
    if (!(rounded >= 1.0f)) {
        return 0;
    }
    if (rounded >= (float)top) {
        return top;
    }
    /* 1 to 127; the bounds are restated for the analyzers, which do not
     * follow the float comparisons above into the conversion. */
    const int edge = (int)rounded;
    if (edge < 0) {
        return 0;
    }
    return (edge > top) ? top : edge;
}

/* Scale each band that carries noise against its own noise power and put it
 * through the output high-pass; a band's edges round to the nearest bin. */
static void
shape_subframe(mbe_fft_plan* plan, const struct mbe_unvoiced_subframe* sub) {
    struct transform_storage* fft = &plan->fft;
    const float spacing = MBE_256_OVER_2PI * sub->w0;
    memset(fft->shaped, 0, MBE_FFT_SIZE * sizeof(float));
    for (int l = 1; l <= sub->L; ++l) {
        const int lower = subframe_band_edge(((float)l - 0.5f) * spacing);
        const int upper = subframe_band_edge(((float)l + 0.5f) * spacing);
        if (!(sub->amplitude[l] > 0.0f) || upper <= lower) {
            continue;
        }
        const float power = band_power(fft->spectrum, lower, upper);
        if (power <= 1.0e-10f) {
            continue;
        }
        const float scale = sub->amplitude[l] / sqrtf(power / (float)(upper - lower));
        /* Bin 0 is real and the high-pass removes it. */
        for (int k = (lower == 0) ? 1 : lower; k < upper; ++k) {
            const float real = fft->spectrum[2u * (size_t)k] * scale;
            const float imag = fft->spectrum[2u * (size_t)k + 1u] * scale;
            const float filter_real = plan->output_filter[2u * (size_t)k];
            const float filter_imag = plan->output_filter[2u * (size_t)k + 1u];
            fft->shaped[2u * (size_t)k] = (real * filter_real) - (imag * filter_imag);
            fft->shaped[2u * (size_t)k + 1u] = (real * filter_imag) + (imag * filter_real);
        }
    }
}

/* One subframe: MBE_SUBFRAME_WINDOW noise samples in, as many windowed output samples out. */
static void
unvoiced_subframe(mbe_fft_plan* plan, const float* noise, const struct mbe_unvoiced_subframe* sub, float* out) {
    struct transform_storage* fft = &plan->fft;
    if (!valid_subframe(sub)) {
        memset(out, 0, MBE_SUBFRAME_WINDOW * sizeof(float));
        return;
    }
    for (int i = 0; i < MBE_SUBFRAME_WINDOW; ++i) {
        fft->time[i] = noise[i] * plan->subframe_window[i];
    }
    memset(fft->time + MBE_SUBFRAME_WINDOW, 0, (MBE_FFT_SIZE - MBE_SUBFRAME_WINDOW) * sizeof(float));
    pffft_transform_ordered(fft->setup, fft->time, fft->spectrum, fft->work, PFFFT_FORWARD);
    shape_subframe(plan, sub);
    pffft_transform_ordered(fft->setup, fft->shaped, fft->time, fft->work, PFFFT_BACKWARD);
    for (int i = 0; i < MBE_SUBFRAME_WINDOW; ++i) {
        out[i] = fft->time[i] * (1.0f / MBE_FFT_SIZE) * plan->subframe_window[i];
    }
}

void
mbe_synthesizeUnvoicedSubframes(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                                const struct mbe_unvoiced_subframe* mid, const struct mbe_unvoiced_subframe* own) {
    float noise[MBE_VOICED_FRAME + MBE_SUBFRAME];
    float subframe[MBE_SUBFRAME_WINDOW];

    if (output == NULL || cur_mp == NULL || prev_mp == NULL || plan == NULL || mid == NULL || own == NULL) {
        return;
    }
    /* Output sample n of the frame is the state n + 1 steps on. */
    uint32_t state = dvsi_noise_state(&cur_mp->noiseSeed);
    for (int i = 0; i < MBE_VOICED_FRAME + MBE_SUBFRAME; ++i) {
        state = mbe_dvsi_noise_next(state);
        noise[i] = (float)(int16_t)(uint16_t)state;
        if (i == MBE_VOICED_FRAME - 1) {
            cur_mp->noiseSeed = (float)state;
        }
    }
    /* The last frame's own subframe plays out over the first half (read
     * before cur_mp's is written: a caller may pass one set as both), unless
     * the other rendering wrote that state: its overlap is not a subframe. */
    if (mbe_subframe_state_is_marked(prev_mp)) {
        for (int n = 0; n < MBE_SUBFRAME; ++n) {
            output[n] += prev_mp->previousUw[MBE_UNVOICED_TAIL_OFFSET + n];
        }
    }
    /* Mid-frame subframe, centred on sample 80; then the frame's own, centred
     * on its end, whose second half the next frame plays. */
    unvoiced_subframe(plan, noise, mid, subframe);
    for (int n = 0; n < MBE_SUBFRAME_WINDOW; ++n) {
        output[n] += subframe[n];
    }
    unvoiced_subframe(plan, noise + MBE_SUBFRAME, own, subframe);
    for (int n = 0; n < MBE_SUBFRAME; ++n) {
        output[MBE_SUBFRAME + n] += subframe[n];
    }
    clear_subframe_tail(cur_mp);
    memcpy(cur_mp->previousUw + MBE_UNVOICED_TAIL_OFFSET, subframe + MBE_SUBFRAME, MBE_SUBFRAME * sizeof(float));
    mbe_subframe_state_adopt(cur_mp);
}

/* `subframe_overlap`: prev_mp's previousUw is the subframe rendering's, which
 * is one subframe's windowed half, not a transform output to window. */
static void
unvoiced_fft_with_noise(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                        const float* noise_buffer, int subframe_overlap) {
    static const float no_overlap[MBE_FFT_SIZE] = {0.0f};
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
    mbe_wola_combine(output, subframe_overlap ? no_overlap : prev_mp->previousUw, fft->time, 160);
    memcpy(cur_mp->previousUw, fft->time, sizeof(cur_mp->previousUw));
}

int
mbe_spec_frame_noise(float* noise, mbe_parms* cur_mp, const mbe_parms* prev_mp) {
    if (noise == NULL || cur_mp == NULL || prev_mp == NULL) {
        return 0;
    }
    /* Read before the draw, which replaces the mark where one state serves
     * as both frames. prev_mp itself is left as it is: a caller may render a
     * frame of its own from a codec path's state and go on with that path. */
    const int subframe_overlap = mbe_subframe_state_is_marked(prev_mp);
    /* What the subframe rendering left in cur_mp goes with its mark: the
     * pulse train in the draw, the overlap here, since a frame can be
     * rejected before it writes its own. */
    if (mbe_subframe_state_is_marked(cur_mp)) {
        memset(cur_mp->previousUw, 0, sizeof(cur_mp->previousUw));
    }
    mbe_generate_noise_with_overlap(noise, &cur_mp->noiseSeed, cur_mp->noiseOverlap);
    return subframe_overlap;
}

void
mbe_synthesizeUnvoicedFFTFrom(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                              const float* noise_buffer, int subframe_overlap) {
    if (output == NULL || !valid_unvoiced_model(cur_mp) || !valid_unvoiced_harmonics(prev_mp) || plan == NULL
        || noise_buffer == NULL) {
        return;
    }
    unvoiced_fft_with_noise(output, cur_mp, prev_mp, plan, noise_buffer, subframe_overlap);
}

void
mbe_synthesizeUnvoicedFFTWithNoise(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                                   const float* noise_buffer) {
    mbe_synthesizeUnvoicedFFTFrom(output, cur_mp, prev_mp, plan, noise_buffer, mbe_subframe_state_is_marked(prev_mp));
}

void
mbe_synthesizeUnvoicedFFT(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan) {
    if (output == NULL || !valid_unvoiced_model(cur_mp) || !valid_unvoiced_harmonics(prev_mp) || plan == NULL) {
        return;
    }
    float noise[MBE_FFT_SIZE];
    const int subframe_overlap = mbe_spec_frame_noise(noise, cur_mp, prev_mp);
    unvoiced_fft_with_noise(output, cur_mp, prev_mp, plan, noise, subframe_overlap);
}
