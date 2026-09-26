// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Shared golden PCM fixtures: the single-frame synthesis model and
 *        error-free, voice-only AMBE frame sequences.
 *
 * Included by tests/test_golden_pcm.c and tools/gen_golden.c so the pinned
 * hashes and the generator always hash identical input.
 */

#ifndef MBELIB_TESTS_GOLDEN_SEQUENCES_H
#define MBELIB_TESTS_GOLDEN_SEQUENCES_H

#include <stddef.h>
#include <stdint.h>

#include "mbelib-neo/mbelib.h"

#define GOLDEN_AMBE_FRAMES    24
#define GOLDEN_AMBE2450_SEED  0x2450A5A5u
#define GOLDEN_AMBE2400_SEED  0x2400C3C3u
#define GOLDEN_FNV1A_OFFSET   2166136261u
#define GOLDEN_SYNTH_RNG_SEED 0xC0FFEEu

typedef int (*golden_ambe_process_fn)(float* aout_buf, mbe_process_result* result, const char ambe_d[49],
                                      mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/** FNV-1a hashes of a whole decoded sequence (float and int16 PCM) plus its energy. */
struct golden_hashes {
    uint32_t f32;
    uint32_t s16;
    double sumsq; /**< Sum of squared float samples, for sanity bounds. */
};

static uint32_t
golden_fnv1a32_update(uint32_t h, const void* data, size_t len) {
    const uint8_t* p = (const uint8_t*)data;
    for (size_t i = 0; i < len; ++i) {
        h ^= p[i];
        h *= 16777619u;
    }
    return h;
}

static uint32_t
golden_xorshift32(uint32_t* state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

/*
 * Bit 0 carries the b0 MSB for both AMBE 3600x2400 and 3600x2450. Clearing it
 * keeps b0 < 64, which rules out tone, erasure and silence frame types, so the
 * sequence exercises only the voice decode and synthesis path.
 */
static void
golden_fill_ambe_voice_frame(char ambe_d[49], uint32_t* state) {
    for (int i = 0; i < 49; ++i) {
        ambe_d[i] = (char)(golden_xorshift32(state) & 1u);
    }
    ambe_d[0] = 0;
}

/**
 * @brief Deterministic single-frame synthesis model for the float/int16 golden hashes.
 *
 * w0 = 0.105 rad with L = 27, the harmonic count the IMBE rule
 * floor(0.9254 * floor(pi / w0 + 0.25)) gives, so every harmonic lies below
 * Nyquist (27 * 0.105 < pi) like the models the decoders produce.
 *
 * @param cur  Output current parameter set.
 * @param prev Output previous parameter set (copy of current).
 */
static void
golden_fill_single_frame(mbe_parms* cur, mbe_parms* prev) {
    mbe_parms enh;
    mbe_initMbeParms(cur, prev, &enh);
    cur->w0 = 0.105f;
    cur->L = 27;
    for (int l = 1; l <= cur->L; ++l) {
        cur->Vl[l] = (l % 4) ? 1 : 0;
        cur->Ml[l] = 0.035f + 0.0015f * (float)l;
        cur->PHIl[l] = (float)l * 0.03f;
        cur->PSIl[l] = (float)l * 0.02f;
    }
    *prev = *cur;
}

/**
 * @brief Decode an error-free voice-only AMBE sequence and hash its PCM.
 * @return 0 on success, or the negative status returned by @p process.
 */
static int
golden_hash_ambe_voice_sequence(golden_ambe_process_fn process, uint32_t seed, struct golden_hashes* out) {
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms prev_enh;
    float pcm_f[160];
    short pcm_s[160];
    char ambe_d[49];
    uint32_t state = seed;

    out->f32 = GOLDEN_FNV1A_OFFSET;
    out->s16 = GOLDEN_FNV1A_OFFSET;
    out->sumsq = 0.0;
    mbe_setThreadRngSeed(GOLDEN_SYNTH_RNG_SEED);
    mbe_initMbeParms(&cur, &prev, &prev_enh);
    for (int f = 0; f < GOLDEN_AMBE_FRAMES; ++f) {
        golden_fill_ambe_voice_frame(ambe_d, &state);
        int rc = process(pcm_f, NULL, ambe_d, &cur, &prev, &prev_enh);
        if (rc < 0) {
            return rc;
        }
        mbe_floattoshort(pcm_f, pcm_s);
        for (int i = 0; i < 160; ++i) {
            out->sumsq += (double)pcm_f[i] * (double)pcm_f[i];
        }
        out->f32 = golden_fnv1a32_update(out->f32, pcm_f, sizeof(pcm_f));
        out->s16 = golden_fnv1a32_update(out->s16, pcm_s, sizeof(pcm_s));
    }
    return 0;
}

#endif /* MBELIB_TESTS_GOLDEN_SEQUENCES_H */
