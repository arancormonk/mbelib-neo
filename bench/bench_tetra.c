// SPDX-License-Identifier: GPL-2.0-or-later
/** @file In-memory TETRA replay; input, allocation, reset and warm-up are untimed. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#if defined(_WIN32)
#include <fcntl.h>
#include <io.h>
#endif
#include "bench_parse.h"
#include "mbelib-neo/mbelib.h"

#define MAX_FRAMES 65536u
#define MAX_RUNS   31

struct frame {
    uint8_t bits[MBE_TETRA_FRAME_BITS];
    int bad;
};

/* CHD: a bad-frame word followed by 137 bit words, all little-endian uint16. */
static int
load_frames(struct frame* frames, size_t* count) {
    unsigned char raw[2 * (MBE_TETRA_FRAME_BITS + 1)];
    *count = 0;
    for (;;) {
        size_t got = fread(raw, 1, sizeof(raw), stdin);
        if (got != sizeof(raw)) {
            return got == 0 && !ferror(stdin) && *count != 0 ? 0 : -1;
        }
        if (*count == MAX_FRAMES) {
            return -1;
        }
        struct frame* f = &frames[(*count)++];
        f->bad = raw[0] != 0 || raw[1] != 0;
        for (size_t b = 0; b < MBE_TETRA_FRAME_BITS; ++b) {
            unsigned v = raw[2 * (b + 1)] | ((unsigned)raw[2 * (b + 1) + 1] << 8);
            if (!f->bad && v > 1) {
                return -1;
            }
            f->bits[b] = (uint8_t)(v & 1u);
        }
    }
}

static int
replay(mbe_tetra_decoder* dec, const struct frame* frames, size_t count, int repeats, uint32_t* checksum) {
    for (int p = 0; p < repeats; ++p) {
        for (size_t i = 0; i < count; ++i) {
            int16_t pcm[MBE_TETRA_FRAME_SAMPLES];
            if (mbe_tetraDecodeFrame(dec, frames[i].bits, frames[i].bad, pcm) != 0) {
                return -1;
            }
            /* Sample every frame; full PCM identity is checked outside the timing harness. */
            *checksum = *checksum * 31u + (uint16_t)pcm[i % MBE_TETRA_FRAME_SAMPLES];
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
    int repeats = 10;
    int runs = 7;
    if (argc > 3 || (argc > 1 && bench_parse_int_arg(argv[1], 1, 10000, &repeats))
        || (argc > 2 && bench_parse_int_arg(argv[2], 1, MAX_RUNS, &runs))) {
        fprintf(stderr, "Usage: %s [repeats=10] [runs=7] < frames.chd\n", argv[0]);
        return 2;
    }
#if defined(_WIN32)
    if (_setmode(_fileno(stdin), _O_BINARY) < 0) {
        return 1;
    }
#endif
    struct frame* frames = (struct frame*)malloc(MAX_FRAMES * sizeof(*frames));
    size_t count = 0;
    if (!frames || load_frames(frames, &count)) {
        fprintf(stderr, "Invalid/empty CHD input, allocation failure or more than %u frames.\n", MAX_FRAMES);
        free(frames);
        return 1;
    }
    mbe_tetra_decoder* dec = mbe_tetraDecoderAlloc();
    int status = dec ? 0 : -1;
    double times[MAX_RUNS];
    uint32_t checksum = 0;
    for (int r = -1; r < runs && status == 0; ++r) {
        mbe_tetraDecoderReset(dec);
        uint32_t run_checksum = 0;
        clock_t start = clock();
        status = replay(dec, frames, count, repeats, &run_checksum);
        double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (r >= 0) {
            /* load_frames rejects empty input; repeats is 1..10000. */
            // NOLINTNEXTLINE(clang-analyzer-optin.taint.TaintedDiv)
            times[r] = elapsed * 1e6 / ((double)count * repeats);
            if (run_checksum != checksum) {
                status = -1;
            }
        }
        checksum = run_checksum;
    }
    mbe_tetraDecoderFree(dec);
    free(frames);
    if (status) {
        fprintf(stderr, "Replay failed or output was not deterministic.\n");
        return 1;
    }
    qsort(times, (size_t)runs, sizeof(times[0]), compare_time);
    double median = 0.5 * (times[(runs - 1) / 2] + times[runs / 2]);
    printf("tetra: median %.3f us/frame (min %.3f, max %.3f), %d runs, %zu frames x %d, checksum %08x\n", median,
           times[0], times[runs - 1], runs, count, repeats, (unsigned)checksum);
    return 0;
}
