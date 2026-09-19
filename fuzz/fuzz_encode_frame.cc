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

static void
fuzz_pcm(const std::uint8_t* data, std::size_t size) {
    mbe_ambe2400_encoder* enc = mbe_ambe2400EncoderAlloc();
    if (enc == nullptr) {
        return;
    }
    mbe_parms cur = {}, prev = {}, enhanced = {};
    mbe_initMbeParms(&cur, &prev, &enhanced);
    // Bound the work per input while exercising state transitions and partial PCM.
    std::size_t limit = size < 3200U ? size : 3200U;
    for (std::size_t offset = 0; offset < limit; offset += 320U) {
        short pcm[160] = {};
        for (std::size_t i = 0; i < 160U && offset + 2U * i + 1U < limit; i++) {
            int value = data[offset + 2U * i] | (static_cast<int>(data[offset + 2U * i + 1U]) << 8);
            pcm[i] = static_cast<short>(value >= 32768 ? value - 65536 : value);
        }
        char bits[49], frame[4][24];
        unsigned char bytes[9];
        check_status(mbe_encodeAmbe2400ParmsShort(enc, pcm, bits, &cur, &prev));
        check_status(mbe_encodeAmbe3600x2400Frame(bits, frame));
        check_status(mbe_encodeDStarDVData(frame, bytes));
        mbe_moveMbeParms(&cur, &prev);
    }
    mbe_ambe2400EncoderFree(enc);
}

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    if (size >= 9U) {
        fuzz_wire_frame(data);
    }
    if (size != 0U) {
        fuzz_pcm(data, size);
    }
    return 0;
}
