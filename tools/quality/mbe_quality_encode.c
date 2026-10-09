// SPDX-License-Identifier: GPL-2.0-or-later
/**
 * @file
 * @brief Encode 8 kHz speech with the AMBE 3600x2400 (D-STAR), AMBE+2
 *        3600x2450 or IMBE 7200x4400 encoder for quality measurement.
 *
 * Writes one row of literal 0/1 parameter bits per 20 ms frame (49 for the
 * AMBE codecs, 88 for IMBE), the Dataf input mbe_quality_eval --codec
 * ambe2400|ambe2450|imbe7200 accepts, so encoder output can be decoded and
 * scored against the input with the existing evaluator.
 */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mbe_quality_fs.h"
#include "mbelib-neo/mbelib.h"

#define FRAME_SAMPLES 160

enum encode_codec { CODEC_AMBE2400, CODEC_AMBE2450, CODEC_IMBE7200 };

struct encode_options {
    const char* input_path;
    const char* output_path;
    long flush_frames;
    enum encode_codec codec;
};

static void
usage(const char* program) {
    fprintf(stderr,
            "Usage: %s --codec ambe2400|ambe2450|imbe7200 --in SPEECH.raw --out ROWS.txt [--flush-frames 0..50]\n",
            program);
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
    for (int i = 1; i < argc; i += 2) {
        const char** value = (i + 1 < argc) ? option_slot(argv[i], options, &slots) : NULL;
        if (!value || *value || argv[i + 1][0] == '\0') {
            usage(argv[0]);
            return 2;
        }
        *value = argv[i + 1];
    }
    options->flush_frames = 1; /* covers the encoders' documented 10 ms analysis delay */
    if (slots.codec && strcmp(slots.codec, "ambe2400") == 0) {
        options->codec = CODEC_AMBE2400;
    } else if (slots.codec && strcmp(slots.codec, "ambe2450") == 0) {
        options->codec = CODEC_AMBE2450;
    } else if (slots.codec && strcmp(slots.codec, "imbe7200") == 0) {
        options->codec = CODEC_IMBE7200;
    } else {
        usage(argv[0]);
        return 2;
    }
    if (!options->input_path || !options->output_path
        || (slots.flush && parse_count(slots.flush, &options->flush_frames) != 0)) {
        usage(argv[0]);
        return 2;
    }
    return 0;
}

static int
write_row(FILE* staged, const char* bits, int count) {
    char row[89];
    for (int i = 0; i < count; ++i) {
        row[i] = (char)('0' + bits[i]);
    }
    row[count] = '\n';
    if (fwrite(row, 1, (size_t)count + 1u, staged) != (size_t)count + 1u) {
        fprintf(stderr, "Cannot stage encoded output.\n");
        return 2;
    }
    return 0;
}

/* One encoder of the selected codec and the caller-owned parameter state all three take. */
struct encoder {
    enum encode_codec codec;
    mbe_ambe2400_encoder* ambe2400;
    mbe_ambe2450_encoder* ambe2450;
    mbe_imbe4400_encoder* imbe;
    mbe_parms cur, prev, enhanced;
};

static int
encoder_open(struct encoder* enc, enum encode_codec codec) {
    memset(enc, 0, sizeof(*enc));
    enc->codec = codec;
    mbe_initMbeParms(&enc->cur, &enc->prev, &enc->enhanced);
    if (codec == CODEC_AMBE2400) {
        enc->ambe2400 = mbe_ambe2400EncoderAlloc();
        return enc->ambe2400 ? 0 : -1;
    }
    if (codec == CODEC_AMBE2450) {
        enc->ambe2450 = mbe_ambe2450EncoderAlloc();
        return enc->ambe2450 ? 0 : -1;
    }
    enc->imbe = mbe_imbe4400EncoderAlloc();
    return enc->imbe ? 0 : -1;
}

static void
encoder_close(struct encoder* enc) {
    mbe_ambe2400EncoderFree(enc->ambe2400);
    mbe_ambe2450EncoderFree(enc->ambe2450);
    mbe_imbe4400EncoderFree(enc->imbe);
}

static int
encode_frame(struct encoder* enc, const short pcm[FRAME_SAMPLES], FILE* staged) {
    char bits[88];
    int count = 49;
    int ret;
    if (enc->codec == CODEC_AMBE2400) {
        ret = mbe_encodeAmbe2400ParmsShort(enc->ambe2400, pcm, bits, &enc->cur, &enc->prev);
    } else if (enc->codec == CODEC_AMBE2450) {
        ret = mbe_encodeAmbe2450ParmsShort(enc->ambe2450, pcm, bits, &enc->cur, &enc->prev);
    } else {
        ret = mbe_encodeImbe4400ParmsShort(enc->imbe, pcm, bits, &enc->cur, &enc->prev);
        count = 88;
    }
    mbe_moveMbeParms(&enc->cur, &enc->prev);
    if (ret < 0) {
        fprintf(stderr, "Encoder rejected a frame.\n");
        return 2;
    }
    return write_row(staged, bits, count);
}

/* Encode every whole or zero-padded final frame, then the flush frames. */
static int
encode_stream(FILE* input, FILE* staged, enum encode_codec codec, long flush_frames) {
    struct encoder enc;
    if (encoder_open(&enc, codec) != 0) {
        encoder_close(&enc);
        fprintf(stderr, "Cannot allocate the encoder.\n");
        return 2;
    }
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
        ret = encode_frame(&enc, pcm, staged);
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
    static const short silence[FRAME_SAMPLES] = {0};
    for (long i = 0; ret == 0 && i < flush_frames; ++i) {
        ret = encode_frame(&enc, silence, staged);
    }
    encoder_close(&enc);
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
    int ret = encode_stream(input, staged, options.codec, options.flush_frames);
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
