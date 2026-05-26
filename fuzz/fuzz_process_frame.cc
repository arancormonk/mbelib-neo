// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

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
    const int uvquality = (size > 2U) ? static_cast<int>((data[2] % 8U) + 1U) : 8;
    const bool raw_bits = (data[0] & 0x80U) != 0U;

    switch (data[0] % 3U) {
        case 0: {
            char imbe_d[88];
            fill_bits(imbe_d, sizeof(imbe_d), data, size, 3U, raw_bits);
            (void)mbe_processImbe4400Dataf(out, &result, imbe_d, &cur, &prev, &prev_enh, uvquality);
            break;
        }
        case 1: {
            char ambe_d[49];
            fill_bits(ambe_d, sizeof(ambe_d), data, size, 3U, raw_bits);
            (void)mbe_processAmbe2400Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh, uvquality);
            break;
        }
        default: {
            char ambe_d[49];
            fill_bits(ambe_d, sizeof(ambe_d), data, size, 3U, raw_bits);
            (void)mbe_processAmbe2450Dataf(out, &result, ambe_d, &cur, &prev, &prev_enh, uvquality);
            break;
        }
    }

    return 0;
}
