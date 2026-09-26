// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include "mbelib-neo/mbelib.h"

static void
fill_bits(char* bits, std::size_t bit_count, const std::uint8_t* data, std::size_t size, std::size_t offset,
          bool raw_bits) {
    std::memset(bits, 0, bit_count);
    if (size <= offset) {
        return;
    }

    const std::uint8_t* payload = data + offset;
    const std::size_t payload_size = size - offset;
    for (std::size_t i = 0; i < bit_count; ++i) {
        const std::uint8_t byte = payload[(i / 8U) % payload_size];
        bits[i] = raw_bits ? static_cast<char>(byte) : static_cast<char>((byte >> (i % 8U)) & 1U);
    }
}

/*
 * AMBE 3600x2450 stream on one parameter triplet, so repeats, mutes, the
 * silence history freeze and mute recovery see real cross-frame state. Each
 * 8-byte chunk carries C0/protected error counts (consistent with the total,
 * C0 context valid) and 7 bytes of parameter bits (MSB first).
 */
static void
process_ambe2450_stream(const std::uint8_t* data, std::size_t size) {
    mbe_parms cur = {};
    mbe_parms prev = {};
    mbe_parms prev_enh = {};
    mbe_initMbeParms(&cur, &prev, &prev_enh);

    const std::size_t max_frames = 64U;
    for (std::size_t f = 0; f < max_frames && (f + 1U) * 8U <= size; ++f) {
        const std::uint8_t* chunk = data + f * 8U;
        char ambe_d[49];
        for (std::size_t i = 0; i < sizeof(ambe_d); ++i) {
            ambe_d[i] = static_cast<char>((chunk[1U + i / 8U] >> (7U - i % 8U)) & 1U);
        }
        mbe_process_result result;
        mbe_initProcessResult(&result);
        result.flags = MBE_PROCESS_FLAG_C0_VALID;
        result.c0_errors = static_cast<int>(chunk[0] & 0x07U);
        result.protected_errors = static_cast<int>((chunk[0] >> 3U) & 0x0FU);
        result.total_errors = result.c0_errors + result.protected_errors;

        float out[160] = {};
        if (mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh) < 0) {
            __builtin_trap();
        }
        for (float sample : out) {
            if (!std::isfinite(sample)) {
                __builtin_trap();
            }
        }
        if (prev_enh.L < 1 || prev_enh.L > 56 || prev.L < 1 || prev.L > 56) {
            __builtin_trap();
        }
    }
}

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    if (data == nullptr || size == 0U) {
        return 0;
    }

    mbe_parms cur = {};
    mbe_parms prev = {};
    mbe_parms prev_enh = {};
    mbe_initMbeParms(&cur, &prev, &prev_enh);

    float out[160] = {};
    mbe_process_result result = {};
    result.total_errors = (size > 1U) ? static_cast<int>(data[1] & 0x0FU) : 0;
    const bool raw_bits = (data[0] & 0x80U) != 0U;

    if ((data[0] & 0x03U) == 0x03U) {
        process_ambe2450_stream(data + 1, size - 1U);
        return 0;
    }

    switch (data[0] % 3U) {
        case 0: {
            char imbe_d[88];
            fill_bits(imbe_d, sizeof(imbe_d), data, size, 3U, raw_bits);
            (void)mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &prev_enh);
            break;
        }
        case 1: {
            char ambe_d[49];
            fill_bits(ambe_d, sizeof(ambe_d), data, size, 3U, raw_bits);
            (void)mbe_processAmbe2400Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh);
            break;
        }
        default: {
            char ambe_d[49];
            fill_bits(ambe_d, sizeof(ambe_d), data, size, 3U, raw_bits);
            (void)mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh);
            break;
        }
    }

    return 0;
}
