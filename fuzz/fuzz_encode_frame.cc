// SPDX-License-Identifier: GPL-2.0-or-later
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "mbelib-neo/mbelib.h"

static void
check_status(int status) {
    if (status < 0) {
        std::abort();
    }
}

static void
fuzz_wire_frame(const std::uint8_t* data) {
    char frame[4][24], encoded[4][24], bits[49];
    unsigned char bytes[9];
    check_status(mbe_decodeDStarDVData(data, frame));
    check_status(mbe_encodeDStarDVData(frame, bytes));
    if (std::memcmp(data, bytes, sizeof(bytes)) != 0) {
        std::abort();
    }
    check_status(mbe_decodeAmbe3600x2400Frame(frame, bits, nullptr));
    check_status(mbe_encodeAmbe3600x2400Frame(bits, encoded));
    check_status(mbe_encodeDStarDVData(encoded, bytes));
}

// Parameter bits from the input, one per byte's low bit.
static void
fuzz_bits(const std::uint8_t* data, std::size_t size, char* bits, std::size_t count) {
    for (std::size_t i = 0; i < count; i++) {
        bits[i] = static_cast<char>(i < size ? data[i] & 1U : 0U);
    }
}

// Parameter bits -> FEC frame -> parameter bits is the identity for every bit.
static void
fuzz_param_frames(const std::uint8_t* data, std::size_t size) {
    char ambe[49], ambe_back[49], ambe_fr[4][24];
    fuzz_bits(data, size, ambe, sizeof(ambe));
    check_status(mbe_encodeAmbe3600x2450Frame(ambe, ambe_fr));
    if (mbe_decodeAmbe3600x2450Frame(ambe_fr, ambe_back, nullptr) != 0 || std::memcmp(ambe, ambe_back, 49) != 0) {
        std::abort();
    }
    char imbe[88], imbe_back[88], imbe_fr[8][23];
    fuzz_bits(data, size, imbe, sizeof(imbe));
    check_status(mbe_encodeImbe7200x4400Frame(imbe, imbe_fr));
    if (mbe_decodeImbe7200x4400Frame(imbe_fr, imbe_back, nullptr) != 0 || std::memcmp(imbe, imbe_back, 88) != 0) {
        std::abort();
    }
}

static void
read_pcm(const std::uint8_t* data, std::size_t offset, std::size_t limit, short pcm[160]) {
    for (std::size_t i = 0; i < 160U; i++) {
        pcm[i] = 0;
        if (offset + 2U * i + 1U < limit) {
            int value = data[offset + 2U * i] | (static_cast<int>(data[offset + 2U * i + 1U]) << 8);
            pcm[i] = static_cast<short>(value >= 32768 ? value - 65536 : value);
        }
    }
}

static void
fuzz_pcm(const std::uint8_t* data, std::size_t size) {
    mbe_ambe2400_encoder* dstar = mbe_ambe2400EncoderAlloc();
    mbe_ambe2450_encoder* ambe = mbe_ambe2450EncoderAlloc();
    mbe_imbe4400_encoder* imbe = mbe_imbe4400EncoderAlloc();
    if (dstar != nullptr && ambe != nullptr && imbe != nullptr) {
        mbe_parms state[3][3] = {};
        for (auto& s : state) {
            mbe_initMbeParms(&s[0], &s[1], &s[2]);
        }
        // Bound the work per input while exercising state transitions and partial PCM.
        std::size_t limit = size < 3200U ? size : 3200U;
        for (std::size_t offset = 0; offset < limit; offset += 320U) {
            short pcm[160];
            char bits[49], imbe_bits[88], frame[4][24], imbe_frame[8][23];
            unsigned char bytes[9];
            read_pcm(data, offset, limit, pcm);
            check_status(mbe_encodeAmbe2400ParmsShort(dstar, pcm, bits, &state[0][0], &state[0][1]));
            check_status(mbe_encodeAmbe3600x2400Frame(bits, frame));
            check_status(mbe_encodeDStarDVData(frame, bytes));
            check_status(mbe_encodeAmbe2450ParmsShort(ambe, pcm, bits, &state[1][0], &state[1][1]));
            check_status(mbe_encodeAmbe3600x2450Frame(bits, frame));
            check_status(mbe_encodeImbe4400ParmsShort(imbe, pcm, imbe_bits, &state[2][0], &state[2][1]));
            check_status(mbe_encodeImbe7200x4400Frame(imbe_bits, imbe_frame));
            for (auto& s : state) {
                mbe_moveMbeParms(&s[0], &s[1]);
            }
        }
    }
    mbe_ambe2400EncoderFree(dstar);
    mbe_ambe2450EncoderFree(ambe);
    mbe_imbe4400EncoderFree(imbe);
}

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    if (size >= 9U) {
        fuzz_wire_frame(data);
    }
    if (size != 0U) {
        fuzz_param_frames(data, size);
        fuzz_pcm(data, size);
    }
    return 0;
}
