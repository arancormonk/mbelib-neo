// SPDX-License-Identifier: GPL-2.0-or-later
/**
 * @file
 * @brief Reframe canonical parameter data for developer quality fixtures.
 */
#include <stdio.h>
#include <string.h>

#include "mbe_quality_dvsi.h"
#include "mbe_quality_frames.h"
#include "mbe_quality_fs.h"

static void
usage(const char* program) {
    fprintf(stderr,
            "Usage: %s --codec imbe7200|imbe7100|ambe2450|ambe2400 --in DATA_FILE --out FRAME_FILE\n"
            "       %s --codec imbe7200|ambe2450|ambe2400 --from-dvsi VECTOR.bit --out FRAME_FILE\n",
            program, program);
}

static int
write_frame_row(FILE* staged, char* frame, int bits) {
    for (int i = 0; i < bits; ++i) {
        frame[i] = (char)('0' + frame[i]);
    }
    frame[bits] = '\n';
    if (fwrite(frame, 1, (size_t)bits + 1, staged) != (size_t)bits + 1) {
        fprintf(stderr, "Cannot stage reframed output.\n");
        return 2;
    }
    return 0;
}

/* DVSI hard-decision vectors: whole binary frames only, channel bits kept. */
static int
reframe_dvsi(FILE* input, FILE* staged, const char* codec) {
    const int frame_bytes = mbe_quality_dvsi_frame_bytes(codec);
    if (frame_bytes < 0) {
        fprintf(stderr, "DVSI import supports imbe7200, ambe2450 and ambe2400.\n");
        return 2;
    }
    unsigned char bytes[18];
    size_t frame_index = 0;
    size_t count;
    while ((count = fread(bytes, 1, (size_t)frame_bytes, input)) == (size_t)frame_bytes) {
        char frame[185];
        int bits = mbe_quality_frame_from_dvsi(codec, bytes, count, frame, sizeof(frame));
        if (bits < 0) {
            fprintf(stderr, "frame %zu: DVSI import failed (%d)\n", frame_index, bits);
            return 2;
        }
        if (write_frame_row(staged, frame, bits) != 0) {
            return 2;
        }
        ++frame_index;
    }
    if (ferror(input)) {
        fprintf(stderr, "Error reading DVSI input.\n");
        return 2;
    }
    if (count != 0) {
        fprintf(stderr, "frame %zu: truncated DVSI frame (%zu of %d bytes)\n", frame_index, count, frame_bytes);
        return 2;
    }
    if (frame_index == 0) {
        fprintf(stderr, "Input contains no DVSI frames.\n");
        return 2;
    }
    return 0;
}

static int
reframe(FILE* input, FILE* staged, const char* codec, size_t width) {
    char line[100];
    size_t frame_index = 0;
    int ch;
    while ((ch = fgetc(input)) != EOF) {
        size_t count = 0;
        do {
            if (count == sizeof(line)) {
                fprintf(stderr, "frame %zu: frame line is too long\n", frame_index);
                return 2;
            }
            line[count++] = (char)ch;
        } while (ch != '\n' && (ch = fgetc(input)) != EOF);
        /* Match the evaluator: accept LF/CRLF or a final unterminated row,
         * but never skip blank rows or whitespace within parameter bits.
         */
        while (count && (line[count - 1] == '\r' || line[count - 1] == '\n')) {
            --count;
        }
        if (count != width) {
            fprintf(stderr, "frame %zu: expected %zu parameter bits for %s, got %zu\n", frame_index, width, codec,
                    count);
            return 2;
        }
        char data[88] = {0};
        for (size_t i = 0; i < count; ++i) {
            if (line[i] != '0' && line[i] != '1') {
                fprintf(stderr, "frame %zu: expected literal binary digits\n", frame_index);
                return 2;
            }
            data[i] = (char)(line[i] - '0');
        }
        char frame[185];
        int bits = mbe_quality_frame_from_data(codec, data, count, frame, sizeof(frame));
        if (bits < 0) {
            fprintf(stderr, "frame %zu: reframing failed (%d)\n", frame_index, bits);
            return 2;
        }
        if (write_frame_row(staged, frame, bits) != 0) {
            return 2;
        }
        ++frame_index;
    }
    if (ferror(input)) {
        fprintf(stderr, "Error reading parameter input.\n");
        return 2;
    }
    if (frame_index == 0) {
        fprintf(stderr, "Input contains no parameter frames.\n");
        return 2;
    }
    return 0;
}

struct reframe_options {
    const char* codec;
    const char* input_path; /* parameter rows, or the DVSI vector when dvsi is set */
    const char* output_path;
    int dvsi;
    size_t width;
};

/* Canonical parameter bits per row, or 0 for an unknown codec. */
static size_t
parameter_width(const char* codec) {
    if (strcmp(codec, "imbe7200") == 0 || strcmp(codec, "imbe7100") == 0) {
        return 88;
    }
    if (strcmp(codec, "ambe2450") == 0 || strcmp(codec, "ambe2400") == 0) {
        return 49;
    }
    return 0;
}

/* DVSI names its 4-bit soft-decision vectors *_sd.bit. Their size is a whole
 * number of hard-decision frames, so only the name tells them apart.
 */
static int
is_soft_decision_vector(const char* path) {
    static const char suffix[] = "_sd.bit";
    const size_t length = strlen(path);
    return length >= sizeof(suffix) - 1 && strcmp(path + length - (sizeof(suffix) - 1), suffix) == 0;
}

/* Parse "--name value" pairs; returns 0 or prints usage and returns 2. */
static int
parse_options(int argc, char** argv, struct reframe_options* options) {
    const char* rows_path = NULL;
    const char* dvsi_path = NULL;
    memset(options, 0, sizeof(*options));
    for (int i = 1; i < argc; i += 2) {
        const char** value = NULL;
        if (i + 1 == argc) {
            value = NULL;
        } else if (strcmp(argv[i], "--codec") == 0) {
            value = &options->codec;
        } else if (strcmp(argv[i], "--in") == 0) {
            value = &rows_path;
        } else if (strcmp(argv[i], "--from-dvsi") == 0) {
            value = &dvsi_path;
        } else if (strcmp(argv[i], "--out") == 0) {
            value = &options->output_path;
        }
        if (!value || *value || argv[i + 1][0] == '\0') {
            usage(argv[0]);
            return 2;
        }
        *value = argv[i + 1];
    }
    const char* codec = options->codec;
    if (!codec || !options->output_path || (rows_path == NULL) == (dvsi_path == NULL)) {
        usage(argv[0]);
        return 2;
    }
    if (dvsi_path && is_soft_decision_vector(dvsi_path)) {
        fprintf(stderr, "DVSI soft-decision vectors (*_sd.bit) are not supported; use the *_hd.bit file.\n");
        return 2;
    }
    options->dvsi = dvsi_path != NULL;
    options->input_path = options->dvsi ? dvsi_path : rows_path;
    options->width = parameter_width(codec);
    if (options->width == 0) {
        usage(argv[0]);
        return 2;
    }
    return 0;
}

/* Copy the fully staged frames to the destination, rechecking aliases. */
static int
publish_staged(FILE* staged, const char* input_path, const char* output_path) {
    if (mbe_quality_paths_equal(input_path, output_path) || mbe_quality_same_file(input_path, output_path)) {
        fprintf(stderr, "Input and output must be different files.\n");
        return 2;
    }
    FILE* output = mbe_quality_open_output(output_path);
    if (!output) {
        fprintf(stderr, "Cannot open output: %s\n", output_path);
        return 2;
    }
    int ret = 0;
    char buffer[4096];
    size_t count;
    while ((count = fread(buffer, 1, sizeof(buffer), staged)) != 0) {
        if (fwrite(buffer, 1, count, output) != count) {
            fprintf(stderr, "Error writing output: %s\n", output_path);
            ret = 2;
            break;
        }
    }
    if (ferror(staged)) {
        fprintf(stderr, "Error reading temporary frame stream.\n");
        ret = 2;
    }
    if (fclose(output) != 0) {
        fprintf(stderr, "Error closing output: %s\n", output_path);
        ret = 2;
    }
    return ret;
}

int
main(int argc, char** argv) {
    struct reframe_options options;
    if (parse_options(argc, argv, &options) != 0) {
        return 2;
    }
    const char* codec = options.codec;
    const char* input_path = options.input_path;
    const char* output_path = options.output_path;
    const int dvsi = options.dvsi;
    const size_t width = options.width;
    if (mbe_quality_paths_equal(input_path, output_path) || mbe_quality_same_file(input_path, output_path)) {
        fprintf(stderr, "Input and output must be different files.\n");
        return 2;
    }
    FILE* input = fopen(input_path, "rb");
    if (!input) {
        fprintf(stderr, "Cannot open input: %s\n", input_path);
        return 2;
    }
    /* Stage every frame before opening the destination: even a malformed
     * final row must leave an existing output untouched. A temporary stream
     * keeps memory bounded and also permits nonseekable parameter inputs.
     */
    FILE* staged = tmpfile();
    if (!staged) {
        fprintf(stderr, "Cannot create temporary frame stream.\n");
        fclose(input);
        return 2;
    }
    int ret = dvsi ? reframe_dvsi(input, staged, codec) : reframe(input, staged, codec, width);
    if (fclose(input) != 0) {
        fprintf(stderr, "Error closing input: %s\n", input_path);
        ret = 2;
    }
    if (ret == 0 && (fflush(staged) != 0 || fseek(staged, 0, SEEK_SET) != 0)) {
        fprintf(stderr, "Cannot rewind temporary frame stream.\n");
        ret = 2;
    }
    if (ret != 0) {
        fclose(staged);
        return ret;
    }
    ret = publish_staged(staged, input_path, output_path);
    if (fclose(staged) != 0) {
        fprintf(stderr, "Error closing temporary frame stream.\n");
        ret = 2;
    }
    return ret;
}
