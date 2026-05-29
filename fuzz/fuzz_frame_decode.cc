// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

#include <cstddef>
#include <cstdint>

#include "mbelib-neo/mbelib.h"

static std::uint8_t
byte_at(const std::uint8_t* data, std::size_t size, std::size_t index) {
    return data[index % size];
}

template <std::size_t Rows, std::size_t Cols>
static void
fill_hard_frame(char (&frame)[Rows][Cols], const std::uint8_t* data, std::size_t size, std::size_t offset,
                bool raw_bits) {
    for (std::size_t row = 0; row < Rows; ++row) {
        for (std::size_t col = 0; col < Cols; ++col) {
            const std::size_t bit_index = row * Cols + col;
            const std::uint8_t value = byte_at(data, size, offset + (bit_index / 8U));
            frame[row][col] = raw_bits ? static_cast<char>(value) : static_cast<char>((value >> (bit_index % 8U)) & 1U);
        }
    }
}

template <std::size_t Rows, std::size_t Cols>
static void
fill_soft_frame(mbe_soft_bit (&frame)[Rows][Cols], const std::uint8_t* data, std::size_t size, std::size_t offset,
                bool raw_bits) {
    for (std::size_t row = 0; row < Rows; ++row) {
        for (std::size_t col = 0; col < Cols; ++col) {
            const std::size_t bit_index = row * Cols + col;
            const std::uint8_t value = byte_at(data, size, offset + (bit_index / 8U));
            const std::uint8_t reliability = byte_at(data, size, offset + 32U + bit_index);
            frame[row][col].bit = raw_bits ? static_cast<std::uint8_t>(value & 3U)
                                           : static_cast<std::uint8_t>((value >> (bit_index % 8U)) & 1U);
            frame[row][col].reliability = reliability;
        }
    }
}

static void
init_state(mbe_parms* cur, mbe_parms* prev, mbe_parms* prev_enh) {
    mbe_initMbeParms(cur, prev, prev_enh);
}

struct FuzzState {
    mbe_parms cur = {};
    mbe_parms prev = {};
    mbe_parms prev_enh = {};
    float out[160] = {};
    mbe_process_result result = {};
    bool raw_bits = false;
};

static FuzzState
make_state(const std::uint8_t* data) {
    FuzzState state = {};
    init_state(&state.cur, &state.prev, &state.prev_enh);
    state.raw_bits = (data[0] & 0x80U) != 0U;
    return state;
}

static void
run_ambe2400_case(std::uint8_t mode, FuzzState& state, const std::uint8_t* data, std::size_t size) {
    switch (mode) {
        case 0: {
            char frame[4][24];
            char bits[49] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeAmbe3600x2400Frame(frame, bits, &state.result);
            break;
        }
        case 1: {
            mbe_soft_bit frame[4][24];
            char bits[49] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeAmbe3600x2400SoftFrame(frame, bits, &state.result);
            break;
        }
        case 2: {
            char frame[4][24];
            char bits[49] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processAmbe3600x2400Framef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                 &state.prev_enh);
            break;
        }
        default: {
            mbe_soft_bit frame[4][24];
            char bits[49] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processAmbe3600x2400SoftFramef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                     &state.prev_enh);
            break;
        }
    }
}

static void
run_ambe2450_case(std::uint8_t mode, FuzzState& state, const std::uint8_t* data, std::size_t size) {
    switch (mode) {
        case 0: {
            char frame[4][24];
            char bits[49] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeAmbe3600x2450Frame(frame, bits, &state.result);
            break;
        }
        case 1: {
            mbe_soft_bit frame[4][24];
            char bits[49] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeAmbe3600x2450SoftFrame(frame, bits, &state.result);
            break;
        }
        case 2: {
            char frame[4][24];
            char bits[49] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processAmbe3600x2450Framef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                 &state.prev_enh);
            break;
        }
        default: {
            mbe_soft_bit frame[4][24];
            char bits[49] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processAmbe3600x2450SoftFramef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                     &state.prev_enh);
            break;
        }
    }
}

static void
run_imbe7200_case(std::uint8_t mode, FuzzState& state, const std::uint8_t* data, std::size_t size) {
    switch (mode) {
        case 0: {
            char frame[8][23];
            char bits[88] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeImbe7200x4400Frame(frame, bits, &state.result);
            break;
        }
        case 1: {
            mbe_soft_bit frame[8][23];
            char bits[88] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeImbe7200x4400SoftFrame(frame, bits, &state.result);
            break;
        }
        case 2: {
            char frame[8][23];
            char bits[88] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processImbe7200x4400Framef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                 &state.prev_enh);
            break;
        }
        default: {
            mbe_soft_bit frame[8][23];
            char bits[88] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processImbe7200x4400SoftFramef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                     &state.prev_enh);
            break;
        }
    }
}

static void
run_imbe7100_case(std::uint8_t mode, FuzzState& state, const std::uint8_t* data, std::size_t size) {
    switch (mode) {
        case 0: {
            char frame[7][24];
            char bits[88] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeImbe7100x4400Frame(frame, bits, &state.result);
            break;
        }
        case 1: {
            mbe_soft_bit frame[7][24];
            char bits[88] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_decodeImbe7100x4400SoftFrame(frame, bits, &state.result);
            break;
        }
        case 2: {
            char frame[7][24];
            char bits[88] = {};
            fill_hard_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processImbe7100x4400Framef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                 &state.prev_enh);
            break;
        }
        default: {
            mbe_soft_bit frame[7][24];
            char bits[88] = {};
            fill_soft_frame(frame, data, size, 2U, state.raw_bits);
            (void)mbe_processImbe7100x4400SoftFramef(state.out, &state.result, frame, bits, &state.cur, &state.prev,
                                                     &state.prev_enh);
            break;
        }
    }
}

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    if (data == nullptr || size == 0U) {
        return 0;
    }

    FuzzState state = make_state(data);
    const std::uint8_t choice = data[0] % 16U;
    const std::uint8_t mode = choice % 4U;

    switch (choice / 4U) {
        case 0: run_ambe2400_case(mode, state, data, size); break;
        case 1: run_ambe2450_case(mode, state, data, size); break;
        case 2: run_imbe7200_case(mode, state, data, size); break;
        default: run_imbe7100_case(mode, state, data, size); break;
    }

    return 0;
}
