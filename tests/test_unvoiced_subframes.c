// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief DVSI's noise generator and 10 ms unvoiced subframes against direct equations.
 */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "dsp_equation_helpers.h"
#include "mbe_adaptive.h"
#include "mbe_unvoiced_fft.h"
#include "mbelib-neo/mbelib.h"

static void
test_generator(void) {
    static unsigned char seen[65536];
    uint32_t state = 0;

    assert(mbe_dvsi_noise_next(0u) == 13849u);
    assert(mbe_dvsi_noise_next(13849u) == 50430u);
    /* Every state once, then back to the start. */
    for (int i = 0; i < 65536; ++i) {
        assert(state < 65536u && !seen[state]);
        seen[state] = 1;
        state = mbe_dvsi_noise_next(state);
    }
    assert(state == 0u);

    uint32_t stepped = 4711u;
    for (int i = 0; i < 160; ++i) {
        stepped = mbe_dvsi_noise_next(stepped);
    }
    assert(mbe_dvsi_noise_skip(4711u, 160) == stepped);
    assert(mbe_dvsi_noise_skip(4711u, 0) == 4711u);
    assert(mbe_dvsi_noise_skip(4711u, -5) == 4711u);
    assert(mbe_dvsi_noise_skip(4711u, 65536) == 4711u);
    /* Only the low 16 bits are state. */
    assert(mbe_dvsi_noise_next(0x12340000u | 77u) == mbe_dvsi_noise_next(77u));
}

/* Whether two runs of samples are the same, bit for bit. */
static int
same_samples(const float* a, const float* b, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        uint32_t a_bits, b_bits;
        memcpy(&a_bits, &a[i], sizeof(a_bits));
        memcpy(&b_bits, &b[i], sizeof(b_bits));
        if (a_bits != b_bits) {
            return 0;
        }
    }
    return 1;
}

/* Whether a stream's noise state is the generator state `state` (a float holds every one exactly). */
static int
noise_state_is(const mbe_parms* mp, uint32_t state) {
    const float want = (float)state;
    uint32_t have_bits, want_bits;
    memcpy(&have_bits, &mp->noiseSeed, sizeof(have_bits));
    memcpy(&want_bits, &want, sizeof(want_bits));
    return have_bits == want_bits;
}

/* The samples a frame draws: sample i of the frame is the state i + 1 steps on. */
static void
frame_noise(uint32_t state, double noise[240]) {
    for (int i = 0; i < 240; ++i) {
        state = mbe_dvsi_noise_next(state);
        noise[i] = (double)(int16_t)(uint16_t)state;
    }
}

static double
subframe_window(int n) {
    const double t = (double)(n - 80);
    return sqrt(t < 0.0 ? (t + 80.0) / 80.0 : (80.0 - t) / 80.0);
}

/* The output high-pass at bin k of the 256-point transform: s^2 / (s^2 + sqrt2 s + 1), s = j hz / 100. */
static void
output_filter(int k, double* real, double* imag) {
    const double x = ((double)k * 8000.0 / 256.0) / 100.0;
    if (k == 0) {
        *real = 0.0;
        *imag = 0.0;
        return;
    }
    const double gain = (x * x) / sqrt((x * x * x * x) + 1.0);
    const double phase = TEST_PI - atan2(sqrt(2.0) * x, 1.0 - (x * x));
    *real = gain * cos(phase);
    *imag = gain * sin(phase);
}

/* One subframe by the equations: window, 256-point DFT, bands scaled against
 * their own power between edges rounded to the nearest bin, the output
 * high-pass, inverse DFT, the window again. */
static void
reference_subframe(const double* noise, const struct mbe_unvoiced_subframe* sub, double out[160]) {
    double time[256] = {0};
    double real[129], imag[129], shaped_real[129] = {0}, shaped_imag[129] = {0};

    for (int n = 0; n < 160; ++n) {
        time[n] = noise[n] * subframe_window(n);
    }
    for (int k = 0; k <= 128; ++k) {
        real[k] = 0.0;
        imag[k] = 0.0;
        for (int n = 0; n < 256; ++n) {
            const double angle = (2.0 * TEST_PI / 256.0) * (double)(k * n);
            real[k] += time[n] * cos(angle);
            imag[k] -= time[n] * sin(angle);
        }
    }
    const double spacing = 256.0 * (double)sub->w0 / (2.0 * TEST_PI);
    for (int l = 1; l <= sub->L; ++l) {
        int lower = (int)floor(((double)l - 0.5) * spacing + 0.5);
        int upper = (int)floor(((double)l + 0.5) * spacing + 0.5);
        lower = lower > 128 ? 128 : lower;
        upper = upper > 128 ? 128 : upper;
        double power = 0.0;
        for (int k = lower; k < upper; ++k) {
            power += real[k] * real[k] + (k == 0 ? 0.0 : imag[k] * imag[k]);
        }
        if (!(sub->amplitude[l] > 0.0f) || upper <= lower || power <= 0.0) {
            continue;
        }
        const double scale = (double)sub->amplitude[l] / sqrt(power / (double)(upper - lower));
        for (int k = lower; k < upper; ++k) {
            double filter_real, filter_imag;
            output_filter(k, &filter_real, &filter_imag);
            const double bin_imag = (k == 0) ? 0.0 : imag[k];
            shaped_real[k] = scale * (real[k] * filter_real - bin_imag * filter_imag);
            shaped_imag[k] = scale * (real[k] * filter_imag + bin_imag * filter_real);
        }
    }
    for (int n = 0; n < 160; ++n) {
        double sample = shaped_real[0];
        for (int k = 1; k < 128; ++k) {
            const double angle = (2.0 * TEST_PI / 256.0) * (double)(k * n);
            sample += 2.0 * (shaped_real[k] * cos(angle) - shaped_imag[k] * sin(angle));
        }
        out[n] = (sample / 256.0) * subframe_window(n);
    }
}

static void
fill_subframe(struct mbe_unvoiced_subframe* sub, float w0, int L, float level, int gap) {
    memset(sub, 0, sizeof(*sub));
    sub->w0 = w0;
    sub->L = L;
    for (int l = 1; l <= L; ++l) {
        sub->amplitude[l] = (gap != 0 && l % gap == 0) ? 0.0f : level * (1.0f + 0.03f * (float)l);
    }
}

static void
test_subframes_match_equations(mbe_fft_plan* plan) {
    mbe_parms current, previous;
    struct mbe_unvoiced_subframe mid, own;
    double noise[240], mid_out[160], own_out[160];
    float output[160];
    const uint32_t state = 40001u;

    new_model(&current);
    new_model(&previous);
    current.noiseSeed = (float)state;
    for (int n = 0; n < 256; ++n) {
        previous.previousUw[n] = (float)(40.0 * sin(0.37 * n));
        current.previousUw[n] = 99.0f;
    }
    mbe_subframe_state_mark(&previous); /* a subframe's second half, as this rendering leaves it */
    fill_subframe(&mid, 0.0613f, 47, 30.0f, 4);
    fill_subframe(&own, 0.1960f, 14, 90.0f, 0);
    for (int n = 0; n < 160; ++n) {
        output[n] = 2.5f;
    }
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &mid, &own);

    frame_noise(state, noise);
    reference_subframe(noise, &mid, mid_out);
    reference_subframe(noise + 80, &own, own_out);
    for (int n = 0; n < 160; ++n) {
        double expected = 2.5 + mid_out[n];
        if (n < 80) {
            expected += (double)previous.previousUw[MBE_UNVOICED_TAIL_OFFSET + n];
        } else {
            expected += own_out[n - 80];
        }
        assert_close(output[n], expected, 0.01, 2.0e-4);
    }
    /* The second half of the own subframe is kept for the next frame; nothing else is. */
    for (int n = 0; n < 256; ++n) {
        const int tail = n - MBE_UNVOICED_TAIL_OFFSET;
        const double expected = (tail >= 0 && tail < 80) ? own_out[80 + tail] : 0.0;
        assert_close(current.previousUw[n], expected, 0.01, 2.0e-4);
    }
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(state, 160)));
    assert(mbe_subframe_state_is_marked(&current));

    /* A state the other rendering wrote becomes this one's whole: its noise
     * overlap is given up, not marked as a pulse train. */
    mbe_parms adopted = previous;
    for (int n = 0; n < 96; ++n) {
        adopted.noiseOverlap[n] = 1000.0f + (float)n;
    }
    adopted.noiseSeed = (float)state;
    float scratch[160] = {0.0f};
    mbe_synthesizeUnvoicedSubframes(scratch, &adopted, &previous, plan, &mid, &own);
    assert(mbe_subframe_state_is_marked(&adopted));
    for (int n = 0; n < MBE_SUBFRAME_STATE_SLOT; ++n) {
        assert(adopted.noiseOverlap[n] == 0.0f);
    }
    for (int n = 0; n < 96; ++n) {
        adopted.noiseOverlap[n] = 1000.0f + (float)n;
    }
    mbe_dvsi_noise_skip_frame(&adopted);
    assert(mbe_subframe_state_is_marked(&adopted));
    for (int n = 0; n < MBE_SUBFRAME_STATE_SLOT; ++n) {
        assert(adopted.noiseOverlap[n] == 0.0f);
    }

    /* An overlap this rendering did not leave is not played. */
    float marked[160] = {0.0f}, unmarked[160] = {0.0f}, none[160] = {0.0f};
    mbe_parms state_a = previous, state_b = previous, state_c = previous;
    mbe_parms next_a = previous, next_b = previous, next_c = previous;
    next_a.noiseSeed = next_b.noiseSeed = next_c.noiseSeed = (float)state;
    state_b.noiseOverlap[MBE_SUBFRAME_STATE_SLOT] = 12345.0f; /* a noise sample, as the other rendering leaves */
    memset(state_c.previousUw, 0, sizeof(state_c.previousUw));
    mbe_synthesizeUnvoicedSubframes(marked, &next_a, &state_a, plan, &mid, &own);
    mbe_synthesizeUnvoicedSubframes(unmarked, &next_b, &state_b, plan, &mid, &own);
    mbe_synthesizeUnvoicedSubframes(none, &next_c, &state_c, plan, &mid, &own);
    assert(same_samples(unmarked, none, 160));
    assert(!same_samples(marked, none, 160));
}

/* A band's noise at its amplitude: level set by MBE_UNVOICED_SUBFRAME_SCALE's convention. */
static void
test_window_squares_sum_to_one(mbe_fft_plan* plan) {
    mbe_parms current, previous;
    struct mbe_unvoiced_subframe model;
    float output[160];
    double energy[2] = {0.0, 0.0};

    new_model(&current);
    new_model(&previous);
    current.noiseSeed = 0.0f;
    fill_subframe(&model, 0.0981748f, 29, 100.0f, 0);
    /* Stationary model: after the first frame every block of the frame holds the same power. */
    for (int frame = 0; frame < 400; ++frame) {
        memset(output, 0, sizeof(output));
        mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &model, &model);
        previous = current; /* the state a stream carries to its next frame */
        if (frame == 0) {
            continue;
        }
        for (int n = 0; n < 160; ++n) {
            energy[(n >= 40 && n < 120) ? 0 : 1] += (double)output[n] * output[n];
        }
    }
    /* Around the subframe centres against around the crossovers. */
    assert_close(energy[0] / energy[1], 1.0, 0.0, 0.03);
}

static void
test_cold_start_and_skip(mbe_fft_plan* plan) {
    mbe_parms current, previous;
    struct mbe_unvoiced_subframe model;
    float output[160] = {0};

    new_model(&current);
    new_model(&previous);
    assert(current.noiseSeed < 0.0f);
    fill_subframe(&model, 0.1960f, 14, 50.0f, 0);
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &model, &model);
    /* A fresh stream starts DVSI's generator from 0. */
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 160)));

    new_model(&current);
    for (int n = 0; n < 256; ++n) {
        current.previousUw[n] = 7.0f;
    }
    mbe_dvsi_noise_skip_frame(&current);
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 160)));
    for (int n = 0; n < 256; ++n) {
        assert(current.previousUw[n] == 0.0f);
    }
    mbe_dvsi_noise_skip_frame(&current);
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 320)));
    mbe_dvsi_noise_skip_frame(NULL);
}

static void
test_invalid_arguments(mbe_fft_plan* plan) {
    mbe_parms current, previous;
    struct mbe_unvoiced_subframe model, bad;
    float output[160];
    const float not_a_number = nanf("");

    new_model(&current);
    new_model(&previous);
    fill_subframe(&model, 0.1960f, 14, 50.0f, 0);
    for (int n = 0; n < 160; ++n) {
        output[n] = 1.0f;
    }
    /* A model no band can be read from draws nothing and still moves the generator. */
    bad = model;
    bad.L = 57;
    current.noiseSeed = 5.0f;
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &bad, &bad);
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(5u, 160)));
    bad = model;
    bad.w0 = not_a_number;
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &bad, &bad);
    bad.w0 = -0.2f;
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &bad, &bad);
    bad.w0 = 0.0f;
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &bad, &bad);
    for (int n = 0; n < 160; ++n) {
        assert(output[n] == 1.0f);
    }
    /* A noise state that is not a generator state restarts it. */
    current.noiseSeed = not_a_number;
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &model, &model);
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 160)));
    current.noiseSeed = 1.0e30f;
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &model, &model);
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 160)));
    for (int n = 0; n < 160; ++n) {
        assert(isfinite(output[n]));
    }
    mbe_synthesizeUnvoicedSubframes(NULL, &current, &previous, plan, &model, &model);
    mbe_synthesizeUnvoicedSubframes(output, NULL, &previous, plan, &model, &model);
    mbe_synthesizeUnvoicedSubframes(output, &current, NULL, plan, &model, &model);
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, NULL, &model);
    mbe_synthesizeUnvoicedSubframes(output, &current, &previous, plan, &model, NULL);
}

/* ---- through the synthesis core ---- */

static void
unvoiced_model(mbe_parms* mp, float level) {
    mp->w0 = (float)(2.0 * TEST_PI / 32.0);
    mp->L = 14;
    for (int l = 1; l <= 56; ++l) {
        mp->Vl[l] = 0;
        mp->Ml[l] = (l <= mp->L) ? level : 0.0f;
    }
}

static double
block_energy(const float* samples, int from, int to) {
    double energy = 0.0;
    for (int n = from; n < to; ++n) {
        energy += (double)samples[n] * samples[n];
    }
    return energy / (double)(to - from);
}

/* A level step in unvoiced speech crosses the whole frame and passes the
 * decibel midpoint at mid-frame, as DVSI's does. */
static void
test_core_level_step(void) {
    double before = 0.0, middle = 0.0, after = 0.0, early = 0.0, late = 0.0;

    for (int trial = 0; trial < 300; ++trial) {
        mbe_parms current, previous, enhanced;
        float output[160];
        mbe_initMbeParms(&current, &previous, &enhanced);
        unvoiced_model(&previous, 400.0f);
        previous.noiseSeed = (float)((trial * 977) & 0xFFFF);
        /* Three steady frames, the step, two more. */
        for (int frame = 0; frame < 6; ++frame) {
            current = previous;
            unvoiced_model(&current, frame < 3 ? 400.0f : 50.0f);
            mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
            if (frame == 2) {
                before += block_energy(output, 0, 160);
            } else if (frame == 3) {
                early += block_energy(output, 16, 48);
                middle += block_energy(output, 64, 96);
                late += block_energy(output, 112, 144);
            } else if (frame == 5) {
                after += block_energy(output, 0, 160);
            }
            previous = current;
        }
    }
    const double step_db = 10.0 * log10(before / after);
    assert_close(step_db, 20.0 * log10(400.0 / 50.0), 0.5, 0.0);
    /* Mid-frame sits half the step down; the quarters lie on either side of it. */
    assert_close(10.0 * log10(before / middle), 0.5 * step_db, 1.0, 0.0);
    assert(10.0 * log10(before / early) > 0.5 && 10.0 * log10(before / early) < 0.5 * step_db - 1.0);
    assert(10.0 * log10(before / late) > 0.5 * step_db + 1.0 && 10.0 * log10(before / late) < step_db - 0.5);
}

static void
test_core_generator_and_silent_history(void) {
    mbe_parms current, previous, enhanced;
    float output[160];

    mbe_initMbeParms(&current, &previous, &enhanced);
    /* A silent last frame: the mid-frame subframe is silent, the own one is not. */
    unvoiced_model(&previous, 0.0f);
    unvoiced_model(&current, 300.0f);
    mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    for (int n = 0; n < 80; ++n) {
        assert(output[n] == 0.0f);
    }
    assert(block_energy(output, 100, 160) > 1.0);
    for (int n = 0; n < 160; ++n) {
        assert(isfinite(output[n]));
    }
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 160)));
    previous = current;
    mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 320)));
}

/* The public synthesis keeps TIA-102.BABA's generator and its silent first buffer. */
static void
test_public_synthesis_keeps_specification(void) {
    mbe_parms current, previous, enhanced;
    float output[160];

    mbe_initMbeParms(&current, &previous, &enhanced);
    unvoiced_model(&previous, 300.0f);
    current = previous;
    mbe_synthesizeSpeechf(output, &current, &previous);
    for (int n = 0; n < 160; ++n) {
        assert(output[n] == 0.0f);
    }
    assert(current.noiseSeed == 3147.0f);
    previous = current;
    mbe_synthesizeSpeechf(output, &current, &previous);
    uint32_t state = 3147u;
    for (int i = 0; i < 160; ++i) {
        state = (171u * state + 11213u) % 53125u;
    }
    assert(current.noiseSeed == (float)state);
}

/* One state handed back and forth between the two renderings stays finite and bounded. */
static void
test_mixed_renderings(void) {
    mbe_parms current, previous, enhanced;
    float output[160];

    mbe_initMbeParms(&current, &previous, &enhanced);
    unvoiced_model(&previous, 300.0f);
    for (int frame = 0; frame < 40; ++frame) {
        current = previous;
        if ((frame / 3) % 2 == 0) {
            mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
        } else {
            mbe_synthesizeSpeechf(output, &current, &previous);
        }
        for (int n = 0; n < 160; ++n) {
            assert(isfinite(output[n]) && fabsf(output[n]) < 20000.0f);
        }
        previous = current;
    }
}

/*
 * The two renderings keep different things in previousUw and noiseOverlap.
 * A stream moved from one to the other starts the new one's unvoiced part
 * from nothing, where playing the other's buffers would put a cut-off
 * transform output, or pulse times read as noise, into the audio.
 */
static void
test_renderings_do_not_play_each_others_state(void) {
    mbe_parms current, previous, enhanced, kept;
    float output[160], reference[160];

    /* TIA-102.BABA frames, then a subframe frame: as from a state with no overlap. */
    mbe_initMbeParms(&current, &previous, &enhanced);
    unvoiced_model(&previous, 300.0f);
    for (int frame = 0; frame < 3; ++frame) {
        current = previous;
        mbe_synthesizeSpeechf(output, &current, &previous);
        previous = current;
    }
    assert(!mbe_subframe_state_is_marked(&previous));
    assert(block_energy(previous.previousUw, 0, 256) > 0.0);
    kept = previous;
    current = previous;
    mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    assert(mbe_subframe_state_is_marked(&current));
    previous = kept;
    memset(previous.previousUw, 0, sizeof(previous.previousUw));
    current = previous;
    mbe_synthesizeRepeatedSpeechf(reference, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    assert(same_samples(output, reference, 160));
    /* It fades in: the first samples carry no more than the rising window lets through. */
    assert(block_energy(output, 0, 8) < 0.02 * block_energy(output, 80, 88));

    /* Subframe frames, then a TIA-102.BABA frame: as from no overlap and no noise overlap. */
    mbe_initMbeParms(&current, &previous, &enhanced);
    unvoiced_model(&previous, 300.0f);
    for (int frame = 0; frame < 3; ++frame) {
        current = previous;
        mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
        previous = current;
    }
    assert(mbe_subframe_state_is_marked(&previous));
    assert(block_energy(previous.previousUw, 0, 256) > 0.0);
    kept = previous;
    current = previous;
    mbe_synthesizeSpeechf(output, &current, &previous);
    assert(!mbe_subframe_state_is_marked(&current));
    previous = kept;
    memset(previous.previousUw, 0, sizeof(previous.previousUw));
    memset(previous.noiseOverlap, 0, sizeof(previous.noiseOverlap));
    current = previous;
    mbe_synthesizeSpeechf(reference, &current, &previous);
    assert(same_samples(output, reference, 160));
    /* The same with one state as both frames, where drawing the noise replaces the mark. */
    current = kept;
    mbe_synthesizeSpeechf(output, &current, &current);
    assert(same_samples(output, reference, 160));

    /* A frame the caller renders on the side leaves the codec path's state as
     * it was: the frame that path plays next is the one it plays without it. */
    mbe_parms scratch = kept, borrowed = kept;
    mbe_synthesizeSpeechf(output, &scratch, &borrowed);
    assert(mbe_subframe_state_is_marked(&borrowed));
    assert(same_samples(borrowed.previousUw, kept.previousUw, 256));
    assert(same_samples(borrowed.noiseOverlap, kept.noiseOverlap, 96));
    assert(same_samples(&borrowed.noiseSeed, &kept.noiseSeed, 1));
    current = kept;
    previous = kept;
    mbe_synthesizeRepeatedSpeechf(reference, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    current = borrowed;
    mbe_synthesizeRepeatedSpeechf(output, &current, &borrowed, MBE_VOICED_SYNTHESIS_DVSI);
    assert(same_samples(output, reference, 160));

    /* A frame this rendering rejects after drawing its noise (a pitch that is
     * none) leaves no subframe overlap behind under an unmarked state: the
     * frame after it plays as from a state with none. */
    mbe_parms rejected = kept, valid_next, clean_next;
    previous = kept;
    rejected.w0 = -0.2f;
    mbe_synthesizeSpeechf(output, &rejected, &previous);
    assert(mbe_subframe_state_is_marked(&rejected) || block_energy(rejected.previousUw, 0, 256) <= 0.0);
    rejected.w0 = kept.w0;
    valid_next = rejected;
    mbe_synthesizeSpeechf(output, &valid_next, &rejected);
    clean_next = rejected;
    memset(rejected.previousUw, 0, sizeof(rejected.previousUw));
    memset(clean_next.previousUw, 0, sizeof(clean_next.previousUw));
    mbe_synthesizeSpeechf(reference, &clean_next, &rejected);
    assert(same_samples(output, reference, 160));

    /* Within one rendering the overlap is played: without it the frame differs. */
    for (int subframes = 0; subframes < 2; ++subframes) {
        mbe_initMbeParms(&current, &previous, &enhanced);
        unvoiced_model(&previous, 300.0f);
        for (int frame = 0; frame < 4; ++frame) {
            current = previous;
            if (frame == 3) {
                kept = previous;
            }
            if (subframes) {
                mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
            } else {
                mbe_synthesizeSpeechf(output, &current, &previous);
            }
            previous = current;
        }
        previous = kept;
        memset(previous.previousUw, 0, sizeof(previous.previousUw));
        current = previous;
        if (subframes) {
            mbe_synthesizeRepeatedSpeechf(reference, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
        } else {
            mbe_synthesizeSpeechf(reference, &current, &previous);
        }
        assert(!same_samples(output, reference, 160));
    }
}

/* The unvoiced synthesis of the TIA-102.BABA rendering called on its own, with
 * one state the subframe rendering wrote as both frames: its overlap is not
 * windowed in, though drawing the noise replaces the state's mark. */
static void
test_spec_unvoiced_with_one_subframe_state(mbe_fft_plan* plan) {
    mbe_parms current, previous, enhanced, kept, reference_state;
    float output[160], reference[160];

    mbe_initMbeParms(&current, &previous, &enhanced);
    unvoiced_model(&previous, 300.0f);
    for (int frame = 0; frame < 3; ++frame) {
        current = previous;
        mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
        previous = current;
    }
    kept = previous;
    assert(mbe_subframe_state_is_marked(&kept) && block_energy(kept.previousUw, 0, 256) > 0.0);

    memset(output, 0, sizeof(output));
    current = kept;
    mbe_synthesizeUnvoicedFFT(output, &current, &current, plan);
    assert(!mbe_subframe_state_is_marked(&current));

    memset(reference, 0, sizeof(reference));
    reference_state = kept;
    memset(reference_state.previousUw, 0, sizeof(reference_state.previousUw));
    memset(reference_state.noiseOverlap, 0, sizeof(reference_state.noiseOverlap));
    current = reference_state;
    mbe_synthesizeUnvoicedFFT(reference, &current, &reference_state, plan);
    assert(block_energy(reference, 0, 160) > 0.0);
    assert(same_samples(output, reference, 160));
}

/* A mixed model: some bands voiced, level and pitch set by the caller. */
static void
mixed_model(mbe_parms* mp, float w0, int L, float level, int voiced_every) {
    mp->w0 = w0;
    mp->L = L;
    for (int l = 1; l <= 56; ++l) {
        mp->Vl[l] = (voiced_every != 0 && l % voiced_every == 0) ? 1 : 0;
        mp->Ml[l] = (l <= L) ? level * (1.0f + 0.01f * (float)l) : 0.0f;
    }
}

/*
 * The mid-frame subframe reads both models: pitch, amplitudes and harmonic
 * counts a caller may have filled with anything. Nothing may index outside a
 * model or the transform, the output stays finite where the inputs are, and
 * ordinary frames after it play as if it had not happened (apart from the
 * generator, which keeps running).
 */
static void
test_core_invalid_models(void) {
    const float bad_pitch[] = {0.0f, -0.0f, -0.3f,    1.0e-42f,        1.0e-30f,        1.0e30f,
                               3.2f, 50.0f, nanf(""), (float)INFINITY, -(float)INFINITY};
    const float bad_level[] = {0.0f, -5.0f, 1.0e-42f, 1.0e30f, nanf(""), (float)INFINITY};
    const int counts[] = {1, 2, 9, 55, 56};
    mbe_parms current, previous, enhanced, clean;
    float output[160], reference[160];

    mbe_initMbeParms(&clean, &previous, &enhanced);
    mixed_model(&clean, 0.21f, 13, 250.0f, 3);

    for (size_t which = 0; which < 2; ++which) {
        for (size_t i = 0; i < sizeof(bad_pitch) / sizeof(bad_pitch[0]); ++i) {
            for (size_t c = 0; c < sizeof(counts) / sizeof(counts[0]); ++c) {
                previous = clean;
                current = clean;
                mbe_parms* broken = which ? &current : &previous;
                broken->w0 = bad_pitch[i];
                broken->L = counts[c];
                mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
                /* An unvoiced frame with a pitch that is not one draws no noise of its own. */
                assert(noise_state_is(&current, mbe_dvsi_noise_skip(0u, 160)));
                for (int n = 0; n < 256; ++n) {
                    assert(isfinite(current.previousUw[n]));
                }
            }
        }
        for (size_t i = 0; i < sizeof(bad_level) / sizeof(bad_level[0]); ++i) {
            previous = clean;
            current = clean;
            mbe_parms* broken = which ? &current : &previous;
            for (int l = 1; l <= 56; ++l) {
                broken->Ml[l] = bad_level[i];
            }
            /* All unvoiced, so only the unvoiced path reads the amplitudes. */
            mixed_model(which ? &previous : &current, 0.21f, 13, 250.0f, 0);
            memset(broken->Vl, 0, sizeof(broken->Vl));
            mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
            for (int n = 0; n < 160; ++n) {
                assert(isfinite(output[n]));
            }
            for (int n = 0; n < 256; ++n) {
                assert(isfinite(current.previousUw[n]));
            }
        }
    }

    /* Harmonic counts out of range are rejected before anything is read. */
    previous = clean;
    current = clean;
    current.L = 57;
    mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    current = clean;
    previous.L = 0;
    mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);

    /* Recovery, unvoiced speech: the frame after a previous model with no
     * usable pitch has a silent mid-frame subframe and its own subframe as
     * usual, and the frame after that plays exactly as in a stream whose
     * broken frame was simply silent. */
    mixed_model(&clean, 0.21f, 13, 250.0f, 0);
    previous = clean;
    current = clean;
    previous.w0 = nanf("");
    mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    for (int n = 0; n < 160; ++n) {
        assert(isfinite(output[n]));
    }
    mbe_parms silent = clean, expected = clean;
    for (int l = 1; l <= 56; ++l) {
        silent.Ml[l] = 0.0f;
    }
    mbe_synthesizeRepeatedSpeechf(reference, &expected, &silent, MBE_VOICED_SYNTHESIS_DVSI);
    assert(compare_bytes(output, reference, sizeof(output)) == 0);
    previous = current;
    mbe_synthesizeRepeatedSpeechf(output, &current, &previous, MBE_VOICED_SYNTHESIS_DVSI);
    silent = expected;
    mbe_synthesizeRepeatedSpeechf(reference, &expected, &silent, MBE_VOICED_SYNTHESIS_DVSI);
    assert(compare_bytes(output, reference, sizeof(output)) == 0);
}

/* A stream's unvoiced noise is its own state: replayed from a copy, or
 * interleaved with another stream, it plays the same samples. */
static void
test_core_streams_replay(void) {
    mbe_parms a, a_prev, b, b_prev, enhanced;
    float expected[8][160], actual[160], discarded[160];

    mbe_initMbeParms(&a, &a_prev, &enhanced);
    mbe_initMbeParms(&b, &b_prev, &enhanced);
    mixed_model(&a_prev, 0.10f, 24, 120.0f, 3);
    mixed_model(&b_prev, 0.125f, 20, 90.0f, 0);
    for (int frame = 0; frame < 3; ++frame) {
        a = a_prev;
        mbe_synthesizeRepeatedSpeechf(discarded, &a, &a_prev, MBE_VOICED_SYNTHESIS_DVSI);
        a_prev = a;
        b = b_prev;
        mbe_synthesizeRepeatedSpeechf(discarded, &b, &b_prev, MBE_VOICED_SYNTHESIS_DVSI);
        b_prev = b;
    }
    mbe_parms replay_prev = a_prev, interleaved_prev = a_prev;
    for (int frame = 0; frame < 8; ++frame) {
        a = a_prev;
        mbe_synthesizeRepeatedSpeechf(expected[frame], &a, &a_prev, MBE_VOICED_SYNTHESIS_DVSI);
        a_prev = a;
    }
    assert(compare_bytes(expected[0], expected[1], sizeof(expected[0])) != 0);
    for (int frame = 0; frame < 8; ++frame) {
        mbe_parms replay = replay_prev;
        mbe_synthesizeRepeatedSpeechf(actual, &replay, &replay_prev, MBE_VOICED_SYNTHESIS_DVSI);
        replay_prev = replay;
        assert(compare_bytes(expected[frame], actual, sizeof(actual)) == 0);
    }
    for (int frame = 0; frame < 8; ++frame) {
        b = b_prev;
        mbe_synthesizeRepeatedSpeechf(discarded, &b, &b_prev, MBE_VOICED_SYNTHESIS_DVSI);
        b_prev = b;
        mbe_parms interleaved = interleaved_prev;
        mbe_synthesizeRepeatedSpeechf(actual, &interleaved, &interleaved_prev, MBE_VOICED_SYNTHESIS_DVSI);
        interleaved_prev = interleaved;
        assert(compare_bytes(expected[frame], actual, sizeof(actual)) == 0);
    }
    /* The thread seed does not reach this generator. */
    mbe_setThreadRngSeed(0xC0FFEEu);
    mbe_initMbeParms(&a, &a_prev, &enhanced);
    mixed_model(&a_prev, 0.10f, 24, 120.0f, 0);
    a = a_prev;
    mbe_synthesizeRepeatedSpeechf(actual, &a, &a_prev, MBE_VOICED_SYNTHESIS_DVSI);
    assert(a.noiseSeed == (float)mbe_dvsi_noise_skip(0u, 160));
}

int
main(void) {
    mbe_fft_plan* plan = mbe_fft_plan_alloc();
    assert(plan != NULL);
    test_generator();
    test_subframes_match_equations(plan);
    test_window_squares_sum_to_one(plan);
    test_cold_start_and_skip(plan);
    test_invalid_arguments(plan);
    test_core_level_step();
    test_core_generator_and_silent_history();
    test_public_synthesis_keeps_specification();
    test_mixed_renderings();
    test_renderings_do_not_play_each_others_state();
    test_spec_unvoiced_with_one_subframe_state(plan);
    test_core_invalid_models();
    test_core_streams_replay();
    mbe_fft_plan_free(plan);
    return 0;
}
