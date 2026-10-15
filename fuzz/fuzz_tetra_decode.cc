// SPDX-License-Identifier: GPL-2.0-or-later
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "mbelib-neo/mbelib.h"

// Each frame takes one control byte (bad-frame flag) and its payload: 18 packed bytes in the
// default mode, or 137 raw bytes (mostly invalid bits) when the first input byte is odd.
static constexpr std::size_t kPackedBytes = (MBE_TETRA_FRAME_BITS + 7) / 8;
static constexpr int kMaxFrames = 48;

static bool
bits_valid(const std::uint8_t* bits) {
    for (int i = 0; i < MBE_TETRA_FRAME_BITS; i++) {
        if (bits[i] > 1U) {
            return false;
        }
    }
    return true;
}

// Decode one frame and check the API contract; returns false when the decoder misbehaves.
static bool
decode_checked(mbe_tetra_decoder* dec, const std::uint8_t* bits, int bad_frame, std::int16_t* pcm) {
    std::int16_t sentinel[MBE_TETRA_FRAME_SAMPLES];
    std::memset(pcm, 0x5A, sizeof(sentinel));
    std::memcpy(sentinel, pcm, sizeof(sentinel));
    int status = mbe_tetraDecodeFrame(dec, bits, bad_frame, pcm);
    if (bad_frame != 0 || bits_valid(bits)) {
        return status == 0;
    }
    // A rejected frame leaves the output alone and can be retried as a bad frame.
    if (status != MBE_STATUS_INVALID_BITS || std::memcmp(pcm, sentinel, sizeof(sentinel)) != 0) {
        return false;
    }
    return mbe_tetraDecodeFrame(dec, bits, 1, pcm) == 0;
}

// The frame's bits: raw bytes, or the packed payload unpacked MSB first.
static void
read_bits(const std::uint8_t* payload, bool raw, std::uint8_t* bits) {
    for (int i = 0; i < MBE_TETRA_FRAME_BITS; i++) {
        const auto index = static_cast<std::size_t>(i);
        bits[i] = raw ? payload[index] : static_cast<std::uint8_t>((payload[index / 8U] >> (7 - i % 8)) & 1U);
    }
}

// A reset context decodes the stream's first frame exactly as the new context did.
static bool
reset_replays(mbe_tetra_decoder* dec, const std::uint8_t* bits, int bad_frame, const std::int16_t* expected) {
    std::int16_t pcm[MBE_TETRA_FRAME_SAMPLES];
    mbe_tetraDecoderReset(dec);
    return mbe_tetraDecodeFrame(dec, bad_frame != 0 ? nullptr : bits, bad_frame, pcm) == 0
           && std::memcmp(pcm, expected, sizeof(pcm)) == 0;
}

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    if (size < 1U) {
        return 0;
    }
    const bool raw = (data[0] & 1U) != 0U;
    const std::size_t step = 1U + (raw ? static_cast<std::size_t>(MBE_TETRA_FRAME_BITS) : kPackedBytes);
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    if (dec == nullptr) {
        return 0;
    }
    std::uint8_t first_bits[MBE_TETRA_FRAME_BITS] = {};
    std::int16_t first_pcm[MBE_TETRA_FRAME_SAMPLES] = {};
    int first_bad = 1;
    int frames = 0;
    for (std::size_t offset = 1U; offset + step <= size && frames < kMaxFrames; offset += step, frames++) {
        std::uint8_t bits[MBE_TETRA_FRAME_BITS];
        std::int16_t pcm[MBE_TETRA_FRAME_SAMPLES];
        const int bad_frame = data[offset] < 16U ? ((data[offset] & 1U) != 0U ? -1 : 1) : 0;
        read_bits(data + offset + 1U, raw, bits);
        if (!decode_checked(dec, bits, bad_frame, pcm)) {
            std::abort();
        }
        if (frames == 0) {
            std::memcpy(first_bits, bits, sizeof(bits));
            std::memcpy(first_pcm, pcm, sizeof(pcm));
            first_bad = (bad_frame != 0 || !bits_valid(bits)) ? 1 : 0;
        }
    }
    if (frames > 0 && !reset_replays(dec, first_bits, first_bad, first_pcm)) {
        std::abort();
    }
    mbe_tetraDecoderFree(dec);
    return 0;
}
