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

#define GOLDEN_AMBE_FRAMES     24
#define GOLDEN_AMBE2450_SEED   0x2450A5A5u
#define GOLDEN_AMBE2400_SEED   0x2400C3C3u
#define GOLDEN_IMBE_FRAMES     24
#define GOLDEN_IMBE_LOUD_SEED  0x4400A5A5u
#define GOLDEN_IMBE_QUIET_SEED 0x44005A5Au
#define GOLDEN_FNV1A_OFFSET    2166136261u
#define GOLDEN_SYNTH_RNG_SEED  0xC0FFEEu

typedef int (*golden_ambe_process_fn)(float* aout_buf, mbe_process_result* result, const char ambe_d[49],
                                      mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
typedef int (*golden_ambe_process_short_fn)(short* aout_buf, mbe_process_result* result, const char ambe_d[49],
                                            mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/** A codec's float and int16 parameter-bit entry points. */
struct golden_ambe_codec {
    golden_ambe_process_fn process_float;
    golden_ambe_process_short_fn process_short;
};

/** FNV-1a hashes of a whole decoded sequence (float and int16 PCM) plus its energy. */
struct golden_hashes {
    uint32_t f32;
    uint32_t s16;
    double sumsq; /**< Sum of squared float samples, for sanity bounds. */
};

/** Continue a 32-bit FNV-1a hash over @p len bytes of @p data. */
static inline uint32_t
golden_fnv1a32_update(uint32_t h, const void* data, size_t len) {
    const uint8_t* p = (const uint8_t*)data;
    for (size_t i = 0; i < len; ++i) {
        h ^= p[i];
        h *= 16777619u;
    }
    return h;
}

/** 32-bit FNV-1a hash of @p len bytes of @p data. */
static inline uint32_t
golden_fnv1a32(const void* data, size_t len) {
    return golden_fnv1a32_update(GOLDEN_FNV1A_OFFSET, data, len);
}

static inline uint32_t
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
static inline void
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
static inline void
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
 *
 * The float hash covers the float entry point's output, the int16 hash what
 * the codec's own int16 entry point returns (which rounds to nearest on the
 * AMBE paths, unlike mbe_floattoshort()). Each pass starts from a freshly
 * initialized and seeded state and decodes the same frames.
 * @return 0 on success, or the negative status returned by an entry point.
 */
static inline int
golden_hash_ambe_voice_sequence(const struct golden_ambe_codec* codec, uint32_t seed, struct golden_hashes* out) {
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms prev_enh;
    float pcm_f[160];
    short pcm_s[160];
    char ambe_d[49];

    out->f32 = GOLDEN_FNV1A_OFFSET;
    out->s16 = GOLDEN_FNV1A_OFFSET;
    out->sumsq = 0.0;
    for (int pass = 0; pass < 2; ++pass) {
        uint32_t state = seed;
        mbe_setThreadRngSeed(GOLDEN_SYNTH_RNG_SEED);
        mbe_initMbeParms(&cur, &prev, &prev_enh);
        for (int f = 0; f < GOLDEN_AMBE_FRAMES; ++f) {
            golden_fill_ambe_voice_frame(ambe_d, &state);
            int rc = (pass == 0) ? codec->process_float(pcm_f, NULL, ambe_d, &cur, &prev, &prev_enh)
                                 : codec->process_short(pcm_s, NULL, ambe_d, &cur, &prev, &prev_enh);
            if (rc < 0) {
                return rc;
            }
            if (pass == 0) {
                for (int i = 0; i < 160; ++i) {
                    out->sumsq += (double)pcm_f[i] * (double)pcm_f[i];
                }
                out->f32 = golden_fnv1a32_update(out->f32, pcm_f, sizeof(pcm_f));
            } else {
                out->s16 = golden_fnv1a32_update(out->s16, pcm_s, sizeof(pcm_s));
            }
        }
    }
    return 0;
}

/*
 * P25 (IMBE 4400) voice frames of random bits at a pitch code under 208 (b0's
 * bits 0..5, 85 and 86), drawn again until the frame's decoded mean log2
 * magnitude (the public parameter decoder on the sequence's own history) lies
 * on the wanted side: quiet under -1.85, where P25's quiet-frame gain plays it
 * lower, loud over -1.5, where it does not. The bits are a function of the
 * seed alone: the history is the decoder's, which no synthesis touches.
 */
static inline int
golden_fill_imbe_voice_frame(char imbe_d[88], uint32_t* state, int quiet, mbe_parms* history) {
    for (int attempt = 0; attempt < 1000; ++attempt) {
        for (int i = 0; i < 88; ++i) {
            imbe_d[i] = (char)(golden_xorshift32(state) & 1u);
        }
        /* b0 = bits 0..5, 85, 86, MSB first: clear bit 0 and bit 1 so b0 < 64, a voice pitch. */
        imbe_d[0] = 0;
        imbe_d[1] = 0;
        mbe_parms cur = *history, prev = *history;
        if (mbe_decodeImbe4400Parms(imbe_d, &cur, &prev) != 0 || cur.L < 1 || cur.L > 56) {
            continue;
        }
        double sum = 0.0;
        for (int l = 1; l <= cur.L; ++l) {
            sum += (double)cur.log2Ml[l];
        }
        const double x = sum / (double)cur.L;
        if (quiet ? (x < -1.85) : (x > -1.5)) {
            *history = cur;
            return 0;
        }
    }
    return -1;
}

/**
 * @brief Decode an error-free P25 voice sequence, quiet or loud, through the
 *        IMBE 4400 data path and hash its PCM (float and int16 entry points).
 *
 * With provoice nonzero the frames go through the same entry points tagged
 * MBE_PROCESS_FLAG_PROVOICE, as ProVoice's frame paths hand them on.
 * @return 0 on success, -1 if no frame could be drawn, or the negative status
 *         an entry point returned.
 */
static inline int
golden_hash_imbe_voice_sequence(uint32_t seed, int quiet, int provoice, struct golden_hashes* out) {
    mbe_parms cur;
    mbe_parms prev;
    mbe_parms prev_enh;
    float pcm_f[160];
    short pcm_s[160];
    char imbe_d[88];

    out->f32 = GOLDEN_FNV1A_OFFSET;
    out->s16 = GOLDEN_FNV1A_OFFSET;
    out->sumsq = 0.0;
    for (int pass = 0; pass < 2; ++pass) {
        uint32_t state = seed;
        mbe_parms history, scratch;
        mbe_initMbeParms(&scratch, &history, &prev_enh);
        mbe_setThreadRngSeed(GOLDEN_SYNTH_RNG_SEED);
        mbe_initMbeParms(&cur, &prev, &prev_enh);
        for (int f = 0; f < GOLDEN_IMBE_FRAMES; ++f) {
            if (golden_fill_imbe_voice_frame(imbe_d, &state, quiet, &history) != 0) {
                return -1;
            }
            mbe_process_result result;
            mbe_initProcessResult(&result);
            if (provoice) {
                result.flags |= MBE_PROCESS_FLAG_PROVOICE;
            }
            int rc = (pass == 0) ? mbe_processImbe4400Dataf(pcm_f, &result, imbe_d, &cur, &prev, &prev_enh)
                                 : mbe_processImbe4400Data(pcm_s, &result, imbe_d, &cur, &prev, &prev_enh);
            if (rc < 0) {
                return rc;
            }
            if (pass == 0) {
                for (int i = 0; i < 160; ++i) {
                    out->sumsq += (double)pcm_f[i] * (double)pcm_f[i];
                }
                out->f32 = golden_fnv1a32_update(out->f32, pcm_f, sizeof(pcm_f));
            } else {
                out->s16 = golden_fnv1a32_update(out->s16, pcm_s, sizeof(pcm_s));
            }
        }
    }
    return 0;
}

/** The AMBE 3600x2450 and D-STAR entry points the golden sequences decode with. */
static const struct golden_ambe_codec golden_ambe2450_codec = {mbe_processAmbe2450Dataf, mbe_processAmbe2450Data};
static const struct golden_ambe_codec golden_ambe2400_codec = {mbe_processAmbe2400Dataf, mbe_processAmbe2400Data};

#endif /* MBELIB_TESTS_GOLDEN_SEQUENCES_H */
