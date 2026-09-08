// SPDX-License-Identifier: GPL-2.0-or-later
/**
 * @file
 * @brief Developer-only ECC/framing fixtures from canonical parameter bits.
 */
#include "mbe_quality_frames.h"

#include <math.h>
#include <string.h>

#include "ecc_const.h"
#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static unsigned int
read_bits(const char* data, int count) {
    unsigned int value = 0;
    for (int i = 0; i < count; ++i) {
        value = (value << 1) | (unsigned int)data[i];
    }
    return value;
}

static unsigned int
encode_golay(unsigned int data) {
    unsigned int parity = 0;
    for (int i = 0; i < 12; ++i) {
        if ((data >> (11 - i)) & 1u) {
            parity ^= (unsigned int)golayGenerator[i];
        }
    }
    return (data << 11) | parity;
}

static void
put_code(char* row, unsigned int code, int count) {
    for (int i = 0; i < count; ++i) {
        row[i] = (char)((code >> i) & 1u);
    }
}

static int
encode_hamming(const char* data, const int generator[4], char* row) {
    unsigned int block = read_bits(data, 11) << 4;
    /* The low four code bits are parity in both Hamming variants. Solve the
     * same zero-syndrome equations as test_ecc rather than assume a variant's
     * parity order matches the other one.
     */
    for (unsigned int parity = 0; parity < 16; ++parity) {
        unsigned int code = block | parity;
        unsigned int syndrome = 0;
        for (int r = 0; r < 4; ++r) {
            unsigned int masked = code & (unsigned int)generator[r];
            unsigned int bit = 0;
            for (int i = 0; i < 15; ++i) {
                bit ^= (masked >> i) & 1u;
            }
            syndrome |= bit << r;
        }
        if (syndrome == 0) {
            put_code(row, code, 15);
            return 0;
        }
    }
    return MBE_STATUS_INVALID_ARGUMENT;
}

static int
pack_ambe(const char* data, char frame[4][24], int is2400) {
    put_code(frame[0] + 1, encode_golay(read_bits(data, 12)), 23);
    /* C0 is extended Golay: column zero makes all 24 bits even parity. */
    for (int i = 1; i < 24; ++i) {
        frame[0][0] ^= frame[0][i];
    }
    put_code(frame[1], encode_golay(read_bits(data + 12, 12)), 23);
    for (int i = 0; i < 11; ++i) {
        frame[2][10 - i] = data[24 + i];
    }
    for (int i = 0; i < 14; ++i) {
        frame[3][13 - i] = data[35 + i];
    }
    return is2400 ? mbe_demodulateAmbe3600x2400Data(frame) : mbe_demodulateAmbe3600x2450Data(frame);
}

static int
pack_imbe7200(const char* data, char frame[8][23]) {
    for (size_t row = 0; row < 4; ++row) {
        put_code(frame[row], encode_golay(read_bits(data + 12 * row, 12)), 23);
    }
    for (size_t row = 4; row < 7; ++row) {
        int ret = encode_hamming(data + 48 + 11 * (row - 4), hammingGenerator, frame[row]);
        if (ret < 0) {
            return ret;
        }
    }
    for (int i = 0; i < 7; ++i) {
        frame[7][6 - i] = data[81 + i];
    }
    return mbe_demodulateImbe7200x4400Data(frame);
}

static void
imbe7100_from_canonical(const char* y, char x[88]) {
    static const int b0_indices[] = {0, 1, 2, 3, 4, 5, 85, 86};
    int b0 = 0;
    for (int i = 0; i < 8; ++i) {
        b0 = (b0 << 1) | y[b0_indices[i]];
    }
    /* Keep the converter's float casts and double operations exactly: K
     * depends on the rounding of w0, including at fundamental boundaries.
     */
    float w0 = ((float)(4 * M_PI) / (float)((float)b0 + 39.5));
    int L = (int)(0.9254 * (int)((M_PI / w0) + 0.25));
    int K = L < 37 ? (int)((float)(L + 2) / (float)3) : 12;

    x[0] = y[87]; /* Preserve the status bit, including nonzero fixtures. */
    x[42] = y[48 + K];
    x[43] = y[49 + K];
    for (int i = 0; i < K; ++i) {
        x[44 + i] = y[48 + i];
    }
    int src = 0;
    int dst = 1;
    while (src < 87) {
        x[dst] = y[src];
        if (++src == 48) {
            src += K + 2;
        }
        if (++dst == 42) {
            dst += K + 2;
        }
    }
}

static int
pack_imbe7100(const char* data, char frame[7][24]) {
    char x[88];
    imbe7100_from_canonical(data, x);
    /* Shortened C0: the high five of the twelve Golay data bits are zero. */
    put_code(frame[0] + 1, encode_golay(read_bits(x, 7)), 18);
    put_code(frame[1] + 1, encode_golay(read_bits(x + 7, 12)), 23);
    put_code(frame[2], encode_golay(read_bits(x + 19, 12)), 23);
    put_code(frame[3], encode_golay(read_bits(x + 31, 12)), 23);
    for (size_t row = 4; row < 6; ++row) {
        int ret = encode_hamming(x + 43 + 11 * (row - 4), imbe7100x4400hammingGenerator, frame[row]);
        if (ret < 0) {
            return ret;
        }
    }
    for (int i = 0; i < 23; ++i) {
        frame[6][22 - i] = x[65 + i];
    }
    /* XOR is its own inverse. The public demodulator also consumes the PR
     * bit for ignored row1 column0; skipping that position shifts later rows.
     */
    return mbe_demodulateImbe7100x4400Data(frame);
}

int
mbe_quality_frame_from_data(const char* codec, const char* data, size_t count, char* frame, size_t capacity) {
    if (!codec || !data || !frame) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    int mode;
    size_t data_count;
    int frame_count;
    if (strcmp(codec, "imbe7200") == 0) {
        mode = 0;
        data_count = 88;
        frame_count = 184;
    } else if (strcmp(codec, "imbe7100") == 0) {
        mode = 1;
        data_count = 88;
        frame_count = 168;
    } else if (strcmp(codec, "ambe2450") == 0 || strcmp(codec, "ambe2400") == 0) {
        mode = strcmp(codec, "ambe2400") == 0 ? 3 : 2;
        data_count = 49;
        frame_count = 96;
    } else {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    if (count != data_count || capacity < (size_t)frame_count) {
        return MBE_STATUS_INVALID_ARGUMENT;
    }
    for (size_t i = 0; i < count; ++i) {
        if (data[i] != 0 && data[i] != 1) {
            return MBE_STATUS_INVALID_BITS;
        }
    }

    /* Stage the rectangular frame so errors never mutate caller output and
     * overlapping input/output remains safe. Unused cells stay numeric zero.
     */
    union {
        char imbe7200[8][23];
        char imbe7100[7][24];
        char ambe[4][24];
    } packed = {0};

    int ret;
    if (mode == 0) {
        ret = pack_imbe7200(data, packed.imbe7200);
    } else if (mode == 1) {
        ret = pack_imbe7100(data, packed.imbe7100);
    } else {
        ret = pack_ambe(data, packed.ambe, mode == 3);
    }
    if (ret < 0) {
        return ret;
    }
    memcpy(frame, &packed, (size_t)frame_count);
    return frame_count;
}
