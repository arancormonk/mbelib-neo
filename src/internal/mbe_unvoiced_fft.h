// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 *
 * Internal header for FFT-based unvoiced synthesis.
 * Implements TIA-102.BABA equations 117-126 for high-quality unvoiced audio.
 */

#ifndef MBEINT_MBE_UNVOICED_FFT_H
#define MBEINT_MBE_UNVOICED_FFT_H

#include <stdint.h>

#include "mbelib-neo/mbelib.h"

/* FFT size for unvoiced synthesis */
#define MBE_FFT_SIZE 256
/* The SSE2/NEON synthesis loops process four floats per step without a
 * scalar tail, so the size must stay a multiple of four. */
#if (MBE_FFT_SIZE % 4) != 0
#error "MBE_FFT_SIZE must be a multiple of four for SIMD synthesis"
#endif

/* Scaling coefficient for unvoiced synthesis (Algorithm #120) */
#define MBE_UNVOICED_SCALE_COEFF 146.17696f

/* 256 / (2 * PI) constant for bin frequency calculation */
#define MBE_256_OVER_2PI         (256.0f / (2.0f * 3.14159265358979323846f))

/* LCG constants from TIA-102.BABA equation 117 */
#define MBE_LCG_A                171.0f
#define MBE_LCG_B                11213.0f
#define MBE_LCG_M                53125.0f

/* Default LCG seed from TIA-102.BABA Annex A */
#define MBE_LCG_DEFAULT_SEED     3147.0f

/* Noise buffer overlap size for continuity between frames */
#define MBE_NOISE_OVERLAP        96

/**
 * @brief Opaque FFT plan handle for unvoiced synthesis.
 */
typedef struct mbe_fft_plan mbe_fft_plan;

/**
 * @brief Allocate an FFT plan for unvoiced synthesis.
 *
 * Creates a reusable plan for 256-point real FFT operations.
 * Call mbe_fft_plan_free() when done.
 *
 * @return Allocated plan, or NULL on failure.
 */
mbe_fft_plan* mbe_fft_plan_alloc(void);

/**
 * @brief Free an FFT plan.
 *
 * @param plan Plan to free.
 */
void mbe_fft_plan_free(mbe_fft_plan* plan);

/**
 * @brief Forward 256-point real FFT using the plan's aligned scratch buffers.
 *
 * Input and output may be unaligned. Output contains DC at [0], Nyquist at
 * [1], and interleaved real/imaginary bins 1..127 at [2..255]. The plan must
 * not be used concurrently. No allocations occur during the transform.
 * @return 0 on success, MBE_STATUS_INVALID_ARGUMENT for NULL arguments.
 */
int mbe_fft_forward_real(mbe_fft_plan* plan, const float input[MBE_FFT_SIZE], float output[MBE_FFT_SIZE]);

/** Transform size of the autocorrelation plan: len + lags must not exceed it. */
#define MBE_ACF_SIZE 512

/**
 * @brief Opaque plan for FFT autocorrelation of short real sequences.
 */
typedef struct mbe_acf_plan mbe_acf_plan;

/**
 * @brief Allocate an autocorrelation plan (512-point real FFT and buffers).
 * @return Allocated plan, or NULL on failure.
 */
mbe_acf_plan* mbe_acf_plan_alloc(void);

/**
 * @brief Free an autocorrelation plan.
 * @param plan Plan to free (NULL is ignored).
 */
void mbe_acf_plan_free(mbe_acf_plan* plan);

/**
 * @brief Autocorrelation r[lag] = sum_n x[n] x[n + lag] for lag = 0..lags.
 *
 * Computed through a zero-padded MBE_ACF_SIZE-point FFT, so it matches the
 * direct sum up to floating-point rounding. No allocations occur.
 * @return 0, or MBE_STATUS_INVALID_ARGUMENT for NULL pointers, a negative
 *         count, lags of MBE_ACF_SIZE or more, or len + lags beyond
 *         MBE_ACF_SIZE.
 */
int mbe_acf_compute(mbe_acf_plan* plan, const float* x, int len, float* r, int lags);

/**
 * @brief Generate LCG noise samples for unvoiced synthesis.
 *
 * Implements the TIA-102.BABA equation 117 linear congruential generator.
 * Uses x' = (171*x + 11213) mod 53125.
 *
 * @param buffer Output buffer for noise samples.
 * @param count Number of samples to generate.
 * @param seed In/out LCG state.
 */
void mbe_generate_noise_lcg(float* buffer, int count, float* seed);

/**
 * @brief Seed the thread-local default LCG state for unvoiced noise.
 *
 * Applies to the next cold-start sequence in mbe_generate_noise_with_overlap().
 *
 * @param seed Thread-local seed value (0 maps to a fixed non-zero default).
 */
void mbe_seedUnvoicedNoiseLcg(uint32_t seed);

/**
 * @brief Generate the next overlapping 256-sample noise buffer.
 *
 * State is encoded in existing public fields to preserve ABI:
 * - overlap[0..95] stores current buffer samples [0..95]
 * - seed stores the current tail generator state for samples [96..255]
 *
 * Library cold-start policy: seed < 0 emits an all-zero buffer,
 * then primes state for the next call with MBE_LCG_DEFAULT_SEED.
 *
 * @param buffer Output buffer of 256 samples.
 * @param seed In/out LCG state.
 * @param overlap In/out 96-sample head/tail state buffer.
 */
void mbe_generate_noise_with_overlap(float* buffer, float* seed, float* overlap);

/**
 * @brief Synthesize unvoiced speech using FFT method.
 *
 * Implements TIA-102.BABA equations 117-126:
 * - Generate noise buffer (with overlap continuity)
 * - Apply synthesis window
 * - FFT -> scale unvoiced bins -> IFFT
 * - WOLA combine with previous frame
 *
 * @param output Output buffer of 160 samples (added to existing content).
 * @param cur_mp Current frame parameters.
 * @param prev_mp Previous frame parameters (provides previousUw for WOLA).
 * @param plan FFT plan (reusable).
 */
void mbe_synthesizeUnvoicedFFT(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan);

/**
 * @brief Synthesize unvoiced speech using a pre-generated noise buffer.
 *
 * Same as mbe_synthesizeUnvoicedFFT but uses the provided noise buffer
 * instead of generating new noise. The caller controls per-frame noise
 * consumption independently of the other synthesis stages.
 *
 * @param output Output buffer of 160 samples (added to existing content).
 * @param cur_mp Current frame parameters.
 * @param prev_mp Previous frame parameters (provides previousUw for WOLA).
 * @param plan FFT plan (reusable).
 * @param noise_buffer Pre-generated 256-sample noise buffer.
 */
void mbe_synthesizeUnvoicedFFTWithNoise(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                                        const float* noise_buffer);

/**
 * @brief Draw a frame's noise for the TIA-102.BABA rendering (eq 117).
 *
 * As mbe_generate_noise_with_overlap() on cur_mp's noise state. Returns
 * whether prev_mp's previousUw is the subframe rendering's (see
 * MBE_SUBFRAME_STATE_MARK), for mbe_synthesizeUnvoicedFFTFrom(): read before
 * the draw, which replaces the mark where one state serves as both frames.
 * prev_mp is not written. What the subframe rendering left in cur_mp, which
 * this frame replaces, is given up with its mark.
 *
 * @param noise   Receives MBE_FFT_SIZE noise samples.
 * @param cur_mp  Current frame: noise state in and out.
 * @param prev_mp Previous frame.
 * @return Nonzero when prev_mp's overlap is not a transform output to window in.
 */
int mbe_spec_frame_noise(float* noise, mbe_parms* cur_mp, const mbe_parms* prev_mp);

/**
 * @brief mbe_synthesizeUnvoicedFFTWithNoise() for noise drawn by
 *        mbe_spec_frame_noise(), with what that returned.
 */
void mbe_synthesizeUnvoicedFFTFrom(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                                   const float* noise_buffer, int subframe_overlap);

/*
 * DVSI's unvoiced synthesis (docs/dvsi-parity.md, "Unvoiced synthesis"), which
 * the codec paths render with MBE_VOICED_SYNTHESIS_DVSI. DVSI's decoders, in
 * every mode, synthesize a 10 ms subframe at a time: the noise under a
 * square-root triangular window two subframes long, a 256-point transform,
 * each unvoiced band scaled against its own noise power, the inverse
 * transform, the same window again and overlap-add. A frame carries two
 * subframes: one mid-frame, whose model lies between the previous frame's and
 * this one's, and the frame's own at its end.
 */
#define MBE_SUBFRAME                80
#define MBE_SUBFRAME_WINDOW         160 /* two subframes */

/* DVSI's noise generator, identical in every mode and stream:
 * x <- 173 x + 13849 mod 2^16 from x = 0, read as a signed 16-bit sample, one
 * step per output sample of every frame that is not a tone frame. */
#define MBE_DVSI_NOISE_MULTIPLIER   173u
#define MBE_DVSI_NOISE_INCREMENT    13849u

/* Where previousUw keeps the second half of a frame's own subframe. */
#define MBE_UNVOICED_TAIL_OFFSET    128

/* The two renderings keep different things in the same fields: the subframe
 * rendering a subframe's windowed second half in previousUw and the pulse
 * train in noiseOverlap, the TIA-102.BABA rendering a whole transform output
 * and noise samples. A state the subframe rendering wrote holds this value in
 * the last element of noiseOverlap, where the other keeps a noise sample
 * (never negative), and neither rendering plays what the other left: a caller
 * that moves a stream from one to the other gets a fade-in, not the other's
 * buffer as audio. */
#define MBE_SUBFRAME_STATE_SLOT     (MBE_NOISE_OVERLAP - 1)
#define MBE_SUBFRAME_STATE_MARK     (-1.0f)

/* An unvoiced band's amplitude against its normalised noise:
 * MBE_UNVOICED_SCALE_COEFF * sqrt(MBE_SUBFRAME_WINDOW / MBE_FFT_SIZE). */
#define MBE_UNVOICED_SUBFRAME_SCALE (MBE_UNVOICED_SCALE_COEFF * 0.790569415f)

/** What one subframe's noise is shaped by. */
struct mbe_unvoiced_subframe {
    float w0;            /**< Fundamental, radians per sample. */
    int L;               /**< Harmonic count. */
    float amplitude[57]; /**< [1..L]: scaled amplitude of each band's noise; 0 for a band with none. */
};

/** @brief Mark a state's overlap buffers as the subframe rendering's. */
void mbe_subframe_state_mark(mbe_parms* mp);

/** @brief Whether a state's overlap buffers are the subframe rendering's. */
int mbe_subframe_state_is_marked(const mbe_parms* mp);

/**
 * @brief Make a state the subframe rendering's before a subframe overlap is
 *        written into it: one that is not marked loses the other rendering's
 *        noise overlap (which is no pulse train) and is marked.
 */
void mbe_subframe_state_adopt(mbe_parms* mp);

/** @brief DVSI's noise generator: the state one sample on. */
uint32_t mbe_dvsi_noise_next(uint32_t state);

/** @brief DVSI's noise generator: the state @p samples samples on. */
uint32_t mbe_dvsi_noise_skip(uint32_t state, int samples);

/**
 * @brief Move a stream's DVSI noise state one frame on without drawing from it.
 *
 * For a frame that plays no synthesized speech but is not a tone frame (a
 * mute): DVSI's generator runs through those and stands still through tones.
 * The second half of the last subframe is dropped, since nothing it could
 * overlap was played.
 *
 * @param mp Parameter set holding the stream's noise state (noiseSeed, previousUw).
 */
void mbe_dvsi_noise_skip_frame(mbe_parms* mp);

/**
 * @brief Add a frame's unvoiced part, as DVSI's decoders synthesize it.
 *
 * Draws MBE_VOICED_FRAME + MBE_SUBFRAME samples of DVSI's noise from
 * cur_mp->noiseSeed and moves the state one frame on. A negative noiseSeed
 * (a fresh stream) starts the generator from 0, as every DVSI stream does;
 * mbe_setThreadRngSeed() does not seed it.
 *
 * @param output  160 samples; the unvoiced part is added.
 * @param cur_mp  Current frame: noise state in and out, previousUw out.
 * @param prev_mp Previous frame: previousUw holds the last subframe's second half.
 * @param plan    Transform plan.
 * @param mid     The mid-frame subframe's model.
 * @param own     The frame's own subframe's model.
 */
void mbe_synthesizeUnvoicedSubframes(float* output, mbe_parms* cur_mp, const mbe_parms* prev_mp, mbe_fft_plan* plan,
                                     const struct mbe_unvoiced_subframe* mid, const struct mbe_unvoiced_subframe* own);

/**
 * @brief Look up the shared synthesis window at an offset from its center.
 *
 * The window is defined for indices -105 to +105.
 * Returns 0 for out-of-range indices.
 *
 * @param n Sample index relative to center.
 * @return Window value at index n.
 */
float mbe_synthesisWindow(int n);

/**
 * @brief Apply WOLA (Weighted Overlap-Add) to combine frames.
 *
 * Implements Algorithm #126 to smoothly combine previous and current
 * frame outputs.
 *
 * @param output Output buffer of N samples.
 * @param prevUw Previous frame's inverse FFT output (256 samples).
 * @param currUw Current frame's inverse FFT output (256 samples).
 * @param N Frame length (typically 160).
 */
void mbe_wola_combine(float* output, const float* prevUw, const float* currUw, int N);

#endif /* MBEINT_MBE_UNVOICED_FFT_H */
