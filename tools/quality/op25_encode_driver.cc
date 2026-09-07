// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Standalone OP25 encoder bridge: 8 kHz s16le PCM to decoder frame bits.
 *
 * Only OP25's headers and bundled mbelib belong in this executable. Linking
 * mbelib-neo here would mix incompatible mbe_parms layouts and function ABIs.
 */

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <optional>
#include <sys/stat.h>

// OP25's ambe_encoder.h requires the preceding declarations.
// clang-format off
#include "imbe_vocoder/imbe_vocoder.h"
#include "mbelib.h"
#include "p25p2_vf.h"
#include "ambe_encoder.h"
// clang-format on

static void
usage(const char* program) {
    std::fprintf(stderr,
                 "Usage: %s --mode imbe7200|ambe2450|ambe2400 --in <s16le-8k.raw> --out <frames.txt> "
                 "[--gain-adjust <log2 attenuation>]\n",
                 program);
}

static void
imbe_bits(imbe_vocoder& encoder, int16_t pcm[160], char line[97]) {
    int16_t words[8];
    encoder.imbe_encode(words, pcm);
    int pos = 0;
    // Pre-FEC u0..u3: 12 bits, u4..u6: 11 bits, u7: 7 bits, MSB first.
    for (int word = 0; word < 8; ++word) {
        const int width = word < 4 ? 12 : (word < 7 ? 11 : 7);
        for (int bit = width - 1; bit >= 0; --bit) {
            line[pos++] = static_cast<char>('0' + ((static_cast<uint16_t>(words[word]) >> bit) & 1));
        }
    }
}

static void
dstar_bits(const uint8_t codeword[72], char line[97]) {
    // Inverse of OP25 encode_dstar's result[d_list[i]] = pre_buf[i].
    static const int d_list[72] = {7,  1,  11, 21, 31, 25, 35, 45, 55, 49, 59, 69, 6,  0,  10, 20, 30, 24,
                                   34, 44, 54, 48, 58, 68, 5,  15, 9,  19, 29, 39, 33, 43, 53, 63, 57, 67,
                                   4,  14, 8,  18, 28, 38, 32, 42, 52, 62, 56, 66, 3,  13, 23, 17, 27, 37,
                                   47, 41, 51, 61, 71, 65, 2,  12, 22, 16, 26, 36, 46, 40, 50, 60, 70, 64};
    uint8_t pre[72];
    for (int i = 0; i < 72; ++i) {
        pre[i] = codeword[d_list[i]];
    }
    std::memset(line, '0', 96);
    for (int i = 0; i < 24; ++i) {
        line[23 - i] = static_cast<char>('0' + pre[i]);
    }
    // Leave the C1 PN modulation intact: the decoder frame API removes it.
    for (int i = 0; i < 23; ++i) {
        line[24 + 22 - i] = static_cast<char>('0' + pre[24 + i]);
    }
    // pre[47] is C1's extra Golay parity bit; the 2400 parameter map skips it.
    for (int i = 0; i < 11; ++i) {
        line[48 + 10 - i] = static_cast<char>('0' + pre[47 + i]);
    }
    for (int i = 0; i < 14; ++i) {
        line[72 + 13 - i] = static_cast<char>('0' + pre[58 + i]);
    }
}

int
main(int argc, char** argv) {
    const char* mode = nullptr;
    const char* input_path = nullptr;
    const char* output_path = nullptr;
    const char* gain_option = nullptr;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 == argc) {
            usage(argv[0]);
            return 2;
        }
        const char** value = nullptr;
        if (std::strcmp(argv[i], "--mode") == 0) {
            value = &mode;
        } else if (std::strcmp(argv[i], "--in") == 0) {
            value = &input_path;
        } else if (std::strcmp(argv[i], "--out") == 0) {
            value = &output_path;
        } else if (std::strcmp(argv[i], "--gain-adjust") == 0) {
            value = &gain_option;
        }
        if (value == nullptr || *value != nullptr || argv[i + 1][0] == '\0') {
            usage(argv[0]);
            return 2;
        }
        *value = argv[i + 1];
    }
    if (mode == nullptr || input_path == nullptr || output_path == nullptr) {
        usage(argv[0]);
        return 2;
    }
    const bool imbe = std::strcmp(mode, "imbe7200") == 0;
    const bool dstar = std::strcmp(mode, "ambe2400") == 0;
    if (!imbe && !dstar && std::strcmp(mode, "ambe2450") != 0) {
        usage(argv[0]);
        return 2;
    }

    struct stat input_stat, output_stat;
    const bool have_input_stat = stat(input_path, &input_stat) == 0;
    if (have_input_stat && stat(output_path, &output_stat) == 0 && input_stat.st_dev == output_stat.st_dev
        && input_stat.st_ino == output_stat.st_ino) {
        std::fprintf(stderr, "Input and output must be different files.\n");
        return 2;
    }
    if (have_input_stat && S_ISREG(input_stat.st_mode) && (input_stat.st_size % 2) != 0) {
        std::fprintf(stderr, "Input contains an odd number of bytes: %s\n", input_path);
        return 2;
    }
    std::ifstream input(input_path, std::ios::binary);
    if (!input) {
        std::fprintf(stderr, "Cannot open input: %s\n", input_path);
        return 2;
    }
    std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
    if (!output) {
        std::fprintf(stderr, "Cannot open output: %s\n", output_path);
        return 2;
    }

    std::optional<imbe_vocoder> imbe_encoder;
    std::optional<ambe_encoder> ambe;
    // OP25 analyzes int16-scale PCM; mbelib-neo's conversion applies gain 7.
    // Compensate in the log2 gain quantizer, without degrading PCM analysis.
    float gain_adjust = 2.807354922057604f; // log2(7); runner calibrates against baseline
    if (gain_option) {
        char* end;
        errno = 0;
        gain_adjust = std::strtof(gain_option, &end);
        if (errno || *end || !std::isfinite(gain_adjust) || std::fabs(gain_adjust) > 16) {
            std::fprintf(stderr, "Gain adjustment must be finite and within [-16, 16].\n");
            return 2;
        }
    }
    if (imbe) {
        imbe_encoder.emplace();
        imbe_encoder->set_gain_adjust(gain_adjust);
    } else {
        ambe.emplace();
        ambe->set_gain_adjust(gain_adjust);
        if (dstar) {
            ambe->set_dstar_mode();
            ambe->set_alt_dstar_interleave(false);
        } else {
            ambe->set_49bit_mode();
        }
    }
    const int bits = imbe ? 88 : (dstar ? 96 : 49);
    for (;;) {
        unsigned char raw[320];
        input.read(reinterpret_cast<char*>(raw), sizeof(raw));
        const std::streamsize bytes = input.gcount();
        if (input.bad() || (input.fail() && !input.eof())) {
            std::fprintf(stderr, "Error reading input: %s\n", input_path);
            return 2;
        }
        if ((bytes % 2) != 0) {
            std::fprintf(stderr, "Input contains an odd number of bytes: %s\n", input_path);
            return 2;
        }
        if (bytes == 0) {
            break;
        }
        int16_t pcm[160] = {};
        for (std::streamsize i = 0; i < bytes / 2; ++i) {
            const int sample = raw[2 * i] | (static_cast<int>(raw[2 * i + 1]) << 8);
            pcm[i] = static_cast<int16_t>(sample < 32768 ? sample : sample - 65536);
        }
        char line[97];
        if (imbe) {
            imbe_bits(*imbe_encoder, pcm, line);
        } else {
            uint8_t codeword[72];
            ambe->encode(pcm, codeword);
            if (dstar) {
                dstar_bits(codeword, line);
            } else {
                for (int i = 0; i < 49; ++i) {
                    line[i] = static_cast<char>('0' + codeword[i]);
                }
            }
        }
        line[bits] = '\n';
        output.write(line, bits + 1);
        if (!output) {
            std::fprintf(stderr, "Error writing output: %s\n", output_path);
            return 2;
        }
        if (input.eof()) {
            break;
        }
    }
    output.close();
    if (!output) {
        std::fprintf(stderr, "Error closing output: %s\n", output_path);
        return 2;
    }
    return 0;
}
