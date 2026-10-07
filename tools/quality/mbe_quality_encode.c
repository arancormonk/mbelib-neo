// SPDX-License-Identifier: GPL-2.0-or-later
/**
 * @file
 * @brief Encode 8 kHz speech with the AMBE 3600x2400 (D-STAR) encoder for
 *        quality measurement.
 *
 * Writes one row of 49 literal 0/1 parameter bits per 20 ms frame, the
 * Dataf input mbe_quality_eval --codec ambe2400 accepts, so encoder output
 * can be decoded and scored against the input with the existing evaluator.
 */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mbe_quality_fs.h"
#include "mbelib-neo/mbelib.h"

#define FRAME_SAMPLES 160

struct encode_options {
    const char* input_path;
    const char* output_path;
    long flush_frames;
};

static void
usage(const char* program) {
    fprintf(stderr, "Usage: %s --codec ambe2400 --in SPEECH.raw --out ROWS.txt [--flush-frames 0..50]\n", program);
}

static int
parse_count(const char* text, long* value) {
    char* end = NULL;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed < 0 || parsed > 50) {
        return -1;
    }
    *value = parsed;
    return 0;
}

struct option_slots {
    const char* codec;
    const char* flush;
};

/* Destination for a "--name" option, or NULL if the name is unknown. */
static const char**
option_slot(const char* name, struct encode_options* options, struct option_slots* slots) {
    if (strcmp(name, "--codec") == 0) {
        return &slots->codec;
    }
    if (strcmp(name, "--in") == 0) {
        return &options->input_path;
    }
    if (strcmp(name, "--out") == 0) {
        return &options->output_path;
    }
    if (strcmp(name, "--flush-frames") == 0) {
        return &slots->flush;
    }
    return NULL;
}

static int
parse_options(int argc, char** argv, struct encode_options* options) {
    struct option_slots slots = {NULL, NULL};
    memset(options, 0, sizeof(*options));
    options->flush_frames = 1; /* covers the encoder's documented 10 ms analysis delay */
    for (int i = 1; i < argc; i += 2) {
        const char** value = (i + 1 < argc) ? option_slot(argv[i], options, &slots) : NULL;
        if (!value || *value || argv[i + 1][0] == '\0') {
            usage(argv[0]);
            return 2;
        }
        *value = argv[i + 1];
    }
    if (!slots.codec || strcmp(slots.codec, "ambe2400") != 0 || !options->input_path || !options->output_path
        || (slots.flush && parse_count(slots.flush, &options->flush_frames) != 0)) {
        usage(argv[0]);
        return 2;
    }
    return 0;
}

static int
write_row(FILE* staged, const char bits[49]) {
    char row[50];
    for (int i = 0; i < 49; ++i) {
        row[i] = (char)('0' + bits[i]);
    }
    row[49] = '\n';
    if (fwrite(row, 1, sizeof(row), staged) != sizeof(row)) {
        fprintf(stderr, "Cannot stage encoded output.\n");
        return 2;
    }
    return 0;
}

static int
encode_frame(mbe_ambe2400_encoder* enc, const short pcm[FRAME_SAMPLES], mbe_parms* cur, mbe_parms* prev, FILE* staged) {
    char bits[49];
    if (mbe_encodeAmbe2400ParmsShort(enc, pcm, bits, cur, prev) < 0) {
        fprintf(stderr, "Encoder rejected a frame.\n");
        return 2;
    }
    mbe_moveMbeParms(cur, prev);
    return write_row(staged, bits);
}

/* Encode every whole or zero-padded final frame, then the flush frames. */
static int
encode_stream(FILE* input, FILE* staged, long flush_frames) {
    mbe_ambe2400_encoder* enc = mbe_ambe2400EncoderAlloc();
    if (!enc) {
        fprintf(stderr, "Cannot allocate the encoder.\n");
        return 2;
    }
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    unsigned char raw[FRAME_SAMPLES * 2];
    size_t count;
    long frames = 0;
    int ret = 0;
    while (ret == 0 && (count = fread(raw, 1, sizeof(raw), input)) != 0) {
        short pcm[FRAME_SAMPLES] = {0};
        if (count % 2 != 0) {
            fprintf(stderr, "Input is not whole 16-bit samples.\n");
            ret = 2;
            break;
        }
        for (size_t i = 0; i < count / 2; ++i) {
            pcm[i] = (short)(unsigned short)(raw[2 * i] | (raw[(2 * i) + 1] << 8));
        }
        ret = encode_frame(enc, pcm, &cur, &prev, staged);
        ++frames;
    }
    if (ret == 0 && ferror(input)) {
        fprintf(stderr, "Error reading input.\n");
        ret = 2;
    }
    if (ret == 0 && frames == 0) {
        fprintf(stderr, "Input contains no samples.\n");
        ret = 2;
    }
    for (long i = 0; ret == 0 && i < flush_frames; ++i) {
        static const short silence[FRAME_SAMPLES];
        ret = encode_frame(enc, silence, &cur, &prev, staged);
    }
    mbe_ambe2400EncoderFree(enc);
    return ret;
}

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
        fprintf(stderr, "Error reading temporary output stream.\n");
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
    struct encode_options options;
    if (parse_options(argc, argv, &options) != 0) {
        return 2;
    }
    if (mbe_quality_paths_equal(options.input_path, options.output_path)
        || mbe_quality_same_file(options.input_path, options.output_path)) {
        fprintf(stderr, "Input and output must be different files.\n");
        return 2;
    }
    FILE* input = fopen(options.input_path, "rb");
    if (!input) {
        fprintf(stderr, "Cannot open input: %s\n", options.input_path);
        return 2;
    }
    /* Stage all rows first so a failure leaves an existing output untouched. */
    FILE* staged = tmpfile();
    if (!staged) {
        fprintf(stderr, "Cannot create temporary output stream.\n");
        fclose(input);
        return 2;
    }
    int ret = encode_stream(input, staged, options.flush_frames);
    if (fclose(input) != 0) {
        fprintf(stderr, "Error closing input: %s\n", options.input_path);
        ret = 2;
    }
    if (ret == 0 && (fflush(staged) != 0 || fseek(staged, 0, SEEK_SET) != 0)) {
        fprintf(stderr, "Cannot rewind temporary output stream.\n");
        ret = 2;
    }
    if (ret == 0) {
        ret = publish_staged(staged, options.input_path, options.output_path);
    }
    if (fclose(staged) != 0) {
        fprintf(stderr, "Error closing temporary output stream.\n");
        ret = 2;
    }
    return ret;
}
