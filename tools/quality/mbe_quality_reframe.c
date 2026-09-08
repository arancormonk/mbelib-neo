// SPDX-License-Identifier: GPL-2.0-or-later
/**
 * @file
 * @brief Reframe canonical parameter data for developer quality fixtures.
 */
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "mbe_quality_frames.h"

static void
usage(const char* program) {
    fprintf(stderr, "Usage: %s --codec imbe7200|imbe7100|ambe2450|ambe2400 --in DATA_FILE --out FRAME_FILE\n", program);
}

static int
same_file(const char* input, const char* output) {
    struct stat input_stat, output_stat;
    return strcmp(input, output) == 0
           || (stat(input, &input_stat) == 0 && stat(output, &output_stat) == 0
               && input_stat.st_dev == output_stat.st_dev && input_stat.st_ino == output_stat.st_ino);
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
        char data[88];
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
        for (int i = 0; i < bits; ++i) {
            frame[i] = (char)('0' + frame[i]);
        }
        frame[bits] = '\n';
        if (fwrite(frame, 1, (size_t)bits + 1, staged) != (size_t)bits + 1) {
            fprintf(stderr, "Cannot stage reframed output.\n");
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

int
main(int argc, char** argv) {
    const char* codec = NULL;
    const char* input_path = NULL;
    const char* output_path = NULL;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 == argc) {
            usage(argv[0]);
            return 2;
        }
        const char** value = NULL;
        if (strcmp(argv[i], "--codec") == 0) {
            value = &codec;
        } else if (strcmp(argv[i], "--in") == 0) {
            value = &input_path;
        } else if (strcmp(argv[i], "--out") == 0) {
            value = &output_path;
        }
        if (!value || *value || argv[i + 1][0] == '\0') {
            usage(argv[0]);
            return 2;
        }
        *value = argv[i + 1];
    }
    if (!codec || !input_path || !output_path) {
        usage(argv[0]);
        return 2;
    }
    size_t width;
    if (strcmp(codec, "imbe7200") == 0 || strcmp(codec, "imbe7100") == 0) {
        width = 88;
    } else if (strcmp(codec, "ambe2450") == 0 || strcmp(codec, "ambe2400") == 0) {
        width = 49;
    } else {
        usage(argv[0]);
        return 2;
    }
    if (same_file(input_path, output_path)) {
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
    int ret = reframe(input, staged, codec, width);
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
    if (same_file(input_path, output_path)) {
        fprintf(stderr, "Input and output must be different files.\n");
        fclose(staged);
        return 2;
    }
    FILE* output = fopen(output_path, "wb");
    if (!output) {
        fprintf(stderr, "Cannot open output: %s\n", output_path);
        fclose(staged);
        return 2;
    }
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
    if (fclose(staged) != 0) {
        fprintf(stderr, "Error closing temporary frame stream.\n");
        ret = 2;
    }
    return ret;
}
