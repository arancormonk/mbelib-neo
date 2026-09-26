// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Utility: generate golden FNV-1a hashes for deterministic synthesis.
 */

#include <stdint.h>
#include <stdio.h>

#include "golden_sequences.h"
#include "mbelib-neo/mbelib.h"

/**
 * @brief Program entry: prints float and int16 golden hashes to stdout.
 */
int
main(void) {
    float out_f[160];
    short out_s[160];
    mbe_parms cur, prev;

    mbe_setThreadRngSeed(0xC0FFEEu);
    golden_fill_single_frame(&cur, &prev);
    mbe_synthesizeSpeechf(out_f, &cur, &prev);

    // Hash float bytes and also short bytes after conversion
    uint32_t hf = golden_fnv1a32(out_f, sizeof(out_f));
    mbe_floattoshort(out_f, out_s);
    uint32_t hs = golden_fnv1a32(out_s, sizeof(out_s));

    printf("GOLDEN_F32_FNV1A=0x%08X\n", (unsigned)hf);
    printf("GOLDEN_S16_FNV1A=0x%08X\n", (unsigned)hs);

    // Error-free voice-only AMBE sequences (see tests/golden_sequences.h)
    struct golden_hashes ambe2450;
    struct golden_hashes ambe2400;
    if (golden_hash_ambe_voice_sequence(mbe_processAmbe2450Dataf, GOLDEN_AMBE2450_SEED, &ambe2450) != 0
        || golden_hash_ambe_voice_sequence(mbe_processAmbe2400Dataf, GOLDEN_AMBE2400_SEED, &ambe2400) != 0) {
        fprintf(stderr, "AMBE golden sequence decode failed\n");
        return 1;
    }
    printf("GOLDEN_AMBE2450_F32_FNV1A=0x%08X\n", (unsigned)ambe2450.f32);
    printf("GOLDEN_AMBE2450_S16_FNV1A=0x%08X\n", (unsigned)ambe2450.s16);
    printf("GOLDEN_AMBE2400_F32_FNV1A=0x%08X\n", (unsigned)ambe2400.f32);
    printf("GOLDEN_AMBE2400_S16_FNV1A=0x%08X\n", (unsigned)ambe2400.s16);
    return 0;
}
