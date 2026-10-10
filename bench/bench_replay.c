// SPDX-License-Identifier: GPL-2.0-or-later
/** @file In-memory codec replay. Input loading, allocation and warm-up are untimed. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if defined(_WIN32)
#include <fcntl.h>
#include <io.h>
#endif
#include "bench_parse.h"
#include "mbelib-neo/mbelib.h"

#define MAX_FRAMES 65536u
#define MAX_RUNS   31

enum codec { AMBE2400, AMBE2450, IMBE7200, IMBE7100 };

enum mode { ENCODE, HARD, SOFT };

struct options {
    enum codec codec;
    enum mode mode;
    int repeats, runs, mixed;
};

struct frame {
    short pcm[160];
    char hard[184];
    mbe_soft_bit soft[184];
};

struct stream {
    mbe_ambe2400_encoder* dstar;
    mbe_ambe2450_encoder* ambe;
    mbe_imbe4400_encoder* imbe;
    mbe_parms cur, prev, enhanced;
};

static int
parse_options(int argc, char** argv, struct options* o) {
    static const char* const codecs[] = {"ambe2400", "ambe2450", "imbe7200", "imbe7100"};
    static const char* const modes[] = {"encode", "hard", "soft"};
    if (argc < 3 || argc > 6) {
        return -1;
    }
    int c = 0, m = 0;
    while (c < 4 && strcmp(argv[2], codecs[c]) != 0) {
        ++c;
    }
    while (m < 3 && strcmp(argv[1], modes[m]) != 0) {
        ++m;
    }
    if (c == 4 || m == 3 || (c == IMBE7100 && m == ENCODE)) {
        return -1;
    }
    o->codec = (enum codec)c;
    o->mode = (enum mode)m;
    o->repeats = 1;
    o->runs = 7;
    o->mixed = argc > 5 && strcmp(argv[5], "mixed") == 0;
    if ((argc > 3 && bench_parse_int_arg(argv[3], 1, 10000, &o->repeats))
        || (argc > 4 && bench_parse_int_arg(argv[4], 1, MAX_RUNS, &o->runs))
        || (argc > 5 && !o->mixed && strcmp(argv[5], "uniform") != 0)) {
        return -1;
    }
    return 0;
}

static int
read_pcm(struct frame* f) {
    unsigned char raw[320] = {0};
    size_t n = fread(raw, 1, sizeof(raw), stdin);
    if (ferror(stdin) || n % 2u != 0u) {
        return -1;
    }
    for (size_t i = 0; i < n / 2u; ++i) {
        unsigned v = (unsigned)raw[2u * i] | ((unsigned)raw[2u * i + 1u] << 8);
        f->pcm[i] = (short)(v < 32768u ? (int)v : (int)v - 65536);
    }
    return n != 0u; /* A partial final frame is zero padded. */
}

static int
read_bits(struct frame* f, int width, int mixed, size_t frame_index) {
    char row[512];
    if (!fgets(row, sizeof(row), stdin)) {
        return ferror(stdin) ? -1 : 0;
    }
    size_t n = strcspn(row, "\r\n");
    if (n != (size_t)width || row[n + strspn(row + n, "\r\n")] != '\0') {
        return -1;
    }
    for (int i = 0; i < width; ++i) {
        if (row[i] != '0' && row[i] != '1') {
            return -1;
        }
        f->hard[i] = (char)(row[i] - '0');
        f->soft[i].bit = (uint8_t)f->hard[i];
        /* Includes zero confidence, full confidence and intermediate weights. */
        f->soft[i].reliability = mixed ? (uint8_t)((frame_index * 17u + (size_t)i * 73u) & 255u) : 255u;
    }
    return 1;
}

static int
load_frames(const struct options* o, struct frame** frames, size_t* count) {
    size_t capacity = 0;
    *frames = NULL;
    *count = 0;
    int width = o->codec == IMBE7200 ? 184 : (o->codec == IMBE7100 ? 168 : 96);
    for (;;) {
        struct frame f = {0};
        int status = o->mode == ENCODE ? read_pcm(&f) : read_bits(&f, width, o->mixed, *count);
        if (status <= 0) {
            return status == 0 && *count > 0u ? 0 : -1;
        }
        if (*count == MAX_FRAMES) {
            return -1;
        }
        if (*count == capacity) {
            size_t next = capacity ? capacity * 2u : 256u;
            struct frame* p = realloc(*frames, next * sizeof(**frames));
            if (!p) {
                return -1;
            }
            *frames = p;
            capacity = next;
        }
        (*frames)[(*count)++] = f;
    }
}

static int
stream_alloc(struct stream* s, const struct options* o) {
    memset(s, 0, sizeof(*s));
    if (o->mode != ENCODE) {
        return 0;
    }
    switch (o->codec) {
        case AMBE2400: s->dstar = mbe_ambe2400EncoderAlloc(); return s->dstar ? 0 : -1;
        case AMBE2450: s->ambe = mbe_ambe2450EncoderAlloc(); return s->ambe ? 0 : -1;
        default: s->imbe = mbe_imbe4400EncoderAlloc(); return s->imbe ? 0 : -1;
    }
}

static void
stream_reset(struct stream* s) {
    mbe_setThreadRngSeed(0x123456u);
    mbe_initMbeParms(&s->cur, &s->prev, &s->enhanced);
    if (s->dstar) {
        mbe_ambe2400EncoderReset(s->dstar);
    }
    if (s->ambe) {
        mbe_ambe2450EncoderReset(s->ambe);
    }
    if (s->imbe) {
        mbe_imbe4400EncoderReset(s->imbe);
    }
}

static int
encode_frame(struct stream* s, const struct frame* f, char bits[88]) {
    int status;
    if (s->dstar) {
        status = mbe_encodeAmbe2400ParmsShort(s->dstar, f->pcm, bits, &s->cur, &s->prev);
    } else if (s->ambe) {
        status = mbe_encodeAmbe2450ParmsShort(s->ambe, f->pcm, bits, &s->cur, &s->prev);
    } else {
        status = mbe_encodeImbe4400ParmsShort(s->imbe, f->pcm, bits, &s->cur, &s->prev);
    }
    mbe_moveMbeParms(&s->cur, &s->prev);
    return status;
}

static int
hard_frame(struct stream* s, enum codec codec, const struct frame* f, float pcm[160], char bits[88]) {
    switch (codec) {
        case AMBE2400:
            return mbe_processAmbe3600x2400Framef(pcm, NULL, (const char (*)[24])f->hard, bits, &s->cur, &s->prev,
                                                  &s->enhanced);
        case AMBE2450:
            return mbe_processAmbe3600x2450Framef(pcm, NULL, (const char (*)[24])f->hard, bits, &s->cur, &s->prev,
                                                  &s->enhanced);
        case IMBE7200:
            return mbe_processImbe7200x4400Framef(pcm, NULL, (const char (*)[23])f->hard, bits, &s->cur, &s->prev,
                                                  &s->enhanced);
        default:
            return mbe_processImbe7100x4400Framef(pcm, NULL, (const char (*)[24])f->hard, bits, &s->cur, &s->prev,
                                                  &s->enhanced);
    }
}

static int
soft_frame(struct stream* s, enum codec codec, const struct frame* f, float pcm[160], char bits[88]) {
    switch (codec) {
        case AMBE2400:
            return mbe_processAmbe3600x2400SoftFramef(pcm, NULL, (const mbe_soft_bit(*)[24])f->soft, bits, &s->cur,
                                                      &s->prev, &s->enhanced);
        case AMBE2450:
            return mbe_processAmbe3600x2450SoftFramef(pcm, NULL, (const mbe_soft_bit(*)[24])f->soft, bits, &s->cur,
                                                      &s->prev, &s->enhanced);
        case IMBE7200:
            return mbe_processImbe7200x4400SoftFramef(pcm, NULL, (const mbe_soft_bit(*)[23])f->soft, bits, &s->cur,
                                                      &s->prev, &s->enhanced);
        default:
            return mbe_processImbe7100x4400SoftFramef(pcm, NULL, (const mbe_soft_bit(*)[24])f->soft, bits, &s->cur,
                                                      &s->prev, &s->enhanced);
    }
}

static int
replay(struct stream* s, const struct options* o, const struct frame* frames, size_t count, uint32_t* checksum) {
    for (int p = 0; p < o->repeats; ++p) {
        for (size_t i = 0; i < count; ++i) {
            char bits[88];
            float pcm[160];
            int status;
            if (o->mode == ENCODE) {
                status = encode_frame(s, frames + i, bits);
            } else if (o->mode == HARD) {
                status = hard_frame(s, o->codec, frames + i, pcm, bits);
            } else {
                status = soft_frame(s, o->codec, frames + i, pcm, bits);
            }
            if (status < 0) {
                return -1;
            }
            int n = o->codec == AMBE2400 || o->codec == AMBE2450 ? 49 : 88;
            for (int b = 0; b < n; ++b) {
                *checksum = *checksum * 31u + (unsigned)bits[b];
            }
            if (o->mode != ENCODE) {
                uint32_t sample;
                memcpy(&sample, pcm + (i % 160u), sizeof(sample));
                *checksum = *checksum * 31u + sample + (unsigned)status;
            }
        }
    }
    return 0;
}

static int
compare_time(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

int
main(int argc, char** argv) {
    struct options o;
    if (parse_options(argc, argv, &o)) {
        fprintf(stderr,
                "Usage: %s encode|hard|soft ambe2400|ambe2450|imbe7200|imbe7100 "
                "[repeats=1] [runs=7] [uniform|mixed] < input\n",
                argv[0]);
        return 2;
    }
#if defined(_WIN32)
    if (_setmode(_fileno(stdin), _O_BINARY) < 0) {
        return 1;
    }
#endif
    if (o.repeats < 1 || o.runs < 1 || o.runs > MAX_RUNS) {
        return 2;
    }
    struct frame* frames;
    size_t count;
    if (load_frames(&o, &frames, &count)) {
        fprintf(stderr, "Invalid/empty input, allocation failure or more than %u frames.\n", MAX_FRAMES);
        free(frames);
        return 1;
    }
    if (count == 0u) {
        free(frames);
        return 1;
    }
    double frames_per_run = (double)count * (double)o.repeats;
    if (frames_per_run <= 0.0) {
        free(frames);
        return 2;
    }
    struct stream s;
    int status = stream_alloc(&s, &o);
    double times[MAX_RUNS];
    uint32_t checksum = 0;
    for (int r = -1; r < o.runs && status == 0; ++r) {
        stream_reset(&s);
        uint32_t run_checksum = 0;
        clock_t start = clock();
        status = replay(&s, &o, frames, count, &run_checksum);
        double elapsed = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
        if (r >= 0) {
            /* count is 1..65536 and repeats is 1..10000. Their double product
             * is exact and positive; the analyzer does not model float ranges. */
            // NOLINTNEXTLINE(clang-analyzer-optin.taint.TaintedDiv)
            times[r] = elapsed * 1e6 / frames_per_run;
            if (run_checksum != checksum) {
                status = -1;
            }
        }
        checksum = run_checksum;
    }
    mbe_ambe2400EncoderFree(s.dstar);
    mbe_ambe2450EncoderFree(s.ambe);
    mbe_imbe4400EncoderFree(s.imbe);
    free(frames);
    if (status) {
        fprintf(stderr, "Replay failed or output was not deterministic.\n");
        return 1;
    }
    qsort(times, (size_t)o.runs, sizeof(times[0]), compare_time);
    static const char* const codec_names[] = {"ambe2400", "ambe2450", "imbe7200", "imbe7100"};
    static const char* const mode_names[] = {"encode", "hard", "soft"};
    double median = 0.5 * (times[(o.runs - 1) / 2] + times[o.runs / 2]);
    printf("%s %s %s: median %.3f us/frame (min %.3f, max %.3f), %d runs, %zu frames x %d, checksum %08x\n",
           mode_names[o.mode], codec_names[o.codec], o.mixed ? "mixed" : "uniform", median, times[0], times[o.runs - 1],
           o.runs, count, o.repeats, (unsigned)checksum);
    return 0;
}
