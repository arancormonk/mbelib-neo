// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/** @file Independent, reference-based measurements of public-API synthesis. */
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "mbelib-neo/mbelib.h"

#define FFT_N 256
#define PI    3.14159265358979323846

typedef struct {
    double* samples;
    size_t count;
    size_t capacity;
} Signal;

typedef struct {
    const char* key;
    double value;
    int valid;
} Metric;

static Metric metrics[64];
static size_t metric_count;

typedef struct {
    int lag;
    double corr;
    int valid;
} Alignment;

typedef struct {
    const char* input_api;
    size_t nonfinite, clips;
    double peak;
    size_t tone, erasure, repeat, mute;
    size_t c0, protected_errors, c4, total;
} DecodeStats;

static const char* reference_state;

static void
fail(const char* message) {
    fprintf(stderr, "%s\n", message);
    exit(2);
}

static FILE*
open_file(const char* path, const char* mode) {
    FILE* f = fopen(path, mode);
    if (!f) {
        fprintf(stderr, "%s: %s\n", path, strerror(errno));
        exit(2);
    }
    return f;
}

/* Resolve symlinks in existing paths/parents as well as lexical "."/"..".
 * Comparing inode identities additionally catches hard-linked inputs. */
static char*
normalized_path(const char* path) {
    char* resolved = realpath(path, NULL);
    if (resolved) {
        return resolved;
    }
    if (errno != ENOENT) {
        fail("cannot resolve input/output path");
    }
    struct stat unresolved;
    if (lstat(path, &unresolved) == 0) {
        fail("cannot safely resolve a dangling input/output symlink");
    }
    char* copy = strdup(path);
    if (!copy) {
        fail("out of memory");
    }
    char* slash = strrchr(copy, '/');
    const char* leaf = slash ? slash + 1 : copy;
    if (!*leaf) {
        fail("output path has an empty filename");
    }
    if (slash) {
        *slash = '\0';
    }
    char* parent = realpath(slash ? (*copy ? copy : "/") : ".", NULL);
    if (!parent) {
        fail("cannot resolve output parent directory");
    }
    size_t size = strlen(parent) + strlen(leaf) + 2;
    resolved = malloc(size);
    if (!resolved) {
        fail("out of memory");
    }
    snprintf(resolved, size, "%s/%s", parent, leaf);
    free(parent);
    free(copy);
    return resolved;
}

static void
reject_alias(const char* a, const char* b) {
    if (!a || !b) {
        return;
    }
    struct stat sa, sb;
    if (stat(a, &sa) == 0 && stat(b, &sb) == 0 && sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino) {
        fail("outputs must be different files from inputs and each other");
    }
    char* na = normalized_path(a);
    char* nb = normalized_path(b);
    int same = !strcmp(na, nb);
    free(na);
    free(nb);
    if (same) {
        fail("outputs must be different paths from inputs and each other");
    }
}

static void
append(Signal* s, double value) {
    if (s->count == s->capacity) {
        size_t cap = s->capacity ? s->capacity * 2 : 4096;
        if (cap < s->capacity || cap > SIZE_MAX / sizeof(double)) {
            fail("signal too large");
        }
        double* p = realloc(s->samples, cap * sizeof(double));
        if (!p) {
            fail("out of memory");
        }
        s->samples = p;
        s->capacity = cap;
    }
    s->samples[s->count++] = value;
}

static uint32_t
le32(const unsigned char* p) {
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

static unsigned
le16(const unsigned char* p) {
    return (unsigned)p[0] | (unsigned)p[1] << 8;
}

static void
put_le(FILE* f, uint32_t value, unsigned bytes) {
    for (unsigned i = 0; i < bytes; ++i) {
        if (fputc((int)((value >> (8 * i)) & 255u), f) == EOF) {
            fail("output write failed");
        }
    }
}

static void
read_exact(FILE* f, void* data, size_t n) {
    if (fread(data, 1, n, f) != n) {
        fail("truncated or unreadable reference");
    }
}

static void
skip_bytes(FILE* f, uint32_t n) {
    unsigned char buf[4096];
    while (n) {
        size_t part = n < sizeof(buf) ? n : sizeof(buf);
        read_exact(f, buf, part);
        n -= (uint32_t)part;
    }
}

static void
read_pcm(FILE* f, Signal* s, uint32_t bytes) {
    if (bytes & 1u) {
        fail("reference has an incomplete s16le sample");
    }
    for (uint32_t i = 0; i < bytes / 2; ++i) {
        unsigned char b[2];
        read_exact(f, b, 2);
        unsigned v = le16(b);
        append(s, v < 32768u ? (double)v : (double)v - 65536.0);
    }
}

static Signal
read_reference(const char* path) {
    Signal s = {0};
    FILE* f = open_file(path, "rb");
    const char* ext = strrchr(path, '.');
    if (ext && strcmp(ext, ".raw") == 0) {
        int a;
        while ((a = fgetc(f)) != EOF) {
            int b = fgetc(f);
            if (b == EOF) {
                fail("reference has an incomplete s16le sample");
            }
            unsigned v = (unsigned)a | (unsigned)b << 8;
            append(&s, v < 32768u ? (double)v : (double)v - 65536.0);
        }
        if (ferror(f)) {
            fail("reference read failed");
        }
    } else if (ext && strcmp(ext, ".wav") == 0) {
        unsigned char h[12];
        read_exact(f, h, 12);
        if (memcmp(h, "RIFF", 4) != 0 || memcmp(h + 8, "WAVE", 4) != 0 || le32(h + 4) < 4) {
            fail("reference is not RIFF/WAVE");
        }
        uint32_t remaining = le32(h + 4) - 4;
        int have_fmt = 0, have_data = 0;
        while (remaining) {
            if (remaining < 8) {
                fail("invalid WAV chunk header");
            }
            read_exact(f, h, 8);
            uint32_t n = le32(h + 4);
            remaining -= 8;
            if (n > remaining || (n & 1u) > remaining - n) {
                fail("WAV chunk exceeds RIFF size");
            }
            if (!memcmp(h, "fmt ", 4)) {
                unsigned char fmt[16];
                if (have_fmt || n < 16) {
                    fail("invalid WAV fmt chunk");
                }
                read_exact(f, fmt, 16);
                if (le16(fmt) != 1 || le16(fmt + 2) != 1 || le32(fmt + 4) != 8000 || le32(fmt + 8) != 16000
                    || le16(fmt + 12) != 2 || le16(fmt + 14) != 16) {
                    fail("reference must be PCM s16le, mono, 8000 Hz");
                }
                skip_bytes(f, n - 16);
                have_fmt = 1;
            } else if (!memcmp(h, "data", 4)) {
                if (have_data) {
                    fail("multiple WAV data chunks are unsupported");
                }
                read_pcm(f, &s, n);
                have_data = 1;
            } else {
                skip_bytes(f, n);
            }
            if (n & 1u) {
                skip_bytes(f, 1);
            }
            remaining -= n + (n & 1u);
        }
        if (!have_fmt || !have_data) {
            fail("WAV needs fmt and data chunks");
        }
    } else {
        fail("reference extension must be .raw or .wav");
    }
    fclose(f);
    if (!s.count) {
        fail("reference is empty");
    }
    return s;
}

static uint32_t
write_wav(const char* path, const Signal* s) {
    if (s->count > (UINT32_MAX - 36u) / 2u) {
        fail("decoded audio exceeds WAV size limit");
    }
    FILE* f = open_file(path, "wb");
    uint32_t bytes = (uint32_t)(s->count * 2);
    fputs("RIFF", f);
    put_le(f, 36 + bytes, 4);
    fputs("WAVEfmt ", f);
    put_le(f, 16, 4);
    put_le(f, 1, 2);
    put_le(f, 1, 2);
    put_le(f, 8000, 4);
    put_le(f, 16000, 4);
    put_le(f, 2, 2);
    put_le(f, 16, 2);
    fputs("data", f);
    put_le(f, bytes, 4);
    uint32_t hash = 2166136261u;
    for (size_t i = 0; i < s->count; ++i) {
        uint16_t v = (uint16_t)(int16_t)s->samples[i];
        put_le(f, v, 2);
        hash = (hash ^ (v & 255u)) * 16777619u;
        hash = (hash ^ (v >> 8)) * 16777619u;
    }
    if (ferror(f) || fclose(f)) {
        fail("WAV write failed");
    }
    return hash;
}

static Signal
decode(const char* codec, const char* path, uint32_t seed, DecodeStats* stats) {
    FILE* f = open_file(path, "rb");
    Signal s = {0};
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    mbe_setThreadRngSeed(seed);
    const int imbe7200 = !strcmp(codec, "imbe7200"), imbe7100 = !strcmp(codec, "imbe7100");
    const int ambe2450 = !strcmp(codec, "ambe2450"), ambe2400 = !strcmp(codec, "ambe2400");
    char line[187];
    size_t frame = 0;
    int ch;
    while ((ch = fgetc(f)) != EOF) {
        size_t n = 0;
        do {
            if (n == sizeof(line)) {
                fail("frame line is too long");
            }
            line[n++] = (char)ch;
        } while (ch != '\n' && (ch = fgetc(f)) != EOF);
        while (n && (line[n - 1] == '\r' || line[n - 1] == '\n')) {
            --n;
        }
        if (!((imbe7200 && (n == 88 || n == 184)) || (imbe7100 && n == 168)
              || ((ambe2450 || ambe2400) && (n == 49 || n == 96)))) {
            fprintf(stderr, "frame %zu: length/codec mismatch\n", frame);
            exit(2);
        }
        char bits[184];
        for (size_t i = 0; i < n; ++i) {
            if (line[i] != '0' && line[i] != '1') {
                fprintf(stderr, "frame %zu: expected literal binary digits\n", frame);
                exit(2);
            }
            bits[i] = (char)(line[i] - '0');
        }
        float audio[160];
        short pcm[160];
        mbe_process_result result;
        mbe_initProcessResult(&result);
        int status;
        const char* api;
        if (imbe7200 && n == 88) {
            api = "mbe_processImbe4400Dataf";
            status = mbe_processImbe4400Dataf(audio, &result, bits, &cur, &prev, &enhanced);
        } else if (imbe7200 && n == 184) {
            char data[88];
            api = "mbe_processImbe7200x4400Framef";
            status =
                mbe_processImbe7200x4400Framef(audio, &result, (const char (*)[23])bits, data, &cur, &prev, &enhanced);
        } else if (imbe7100 && n == 168) {
            char data[88];
            api = "mbe_processImbe7100x4400Framef";
            status =
                mbe_processImbe7100x4400Framef(audio, &result, (const char (*)[24])bits, data, &cur, &prev, &enhanced);
        } else if (ambe2450 && n == 96) {
            char data[49];
            api = "mbe_processAmbe3600x2450Framef";
            status =
                mbe_processAmbe3600x2450Framef(audio, &result, (const char (*)[24])bits, data, &cur, &prev, &enhanced);
        } else if (ambe2400 && n == 96) {
            char data[49];
            api = "mbe_processAmbe3600x2400Framef";
            status =
                mbe_processAmbe3600x2400Framef(audio, &result, (const char (*)[24])bits, data, &cur, &prev, &enhanced);
        } else if (ambe2450 && n == 49) {
            api = "mbe_processAmbe2450Dataf";
            status = mbe_processAmbe2450Dataf(audio, &result, bits, &cur, &prev, &enhanced);
        } else {
            api = "mbe_processAmbe2400Dataf";
            status = mbe_processAmbe2400Dataf(audio, &result, bits, &cur, &prev, &enhanced);
        }
        if (status < 0) {
            fprintf(stderr, "frame %zu: process status %d\n", frame, status);
            exit(3);
        }
        if (stats->input_api && strcmp(stats->input_api, api) != 0) {
            fail("frame file must use one public input API");
        }
        stats->input_api = api;
        stats->tone += (result.flags & MBE_PROCESS_FLAG_TONE) != 0;
        stats->erasure += (result.flags & MBE_PROCESS_FLAG_ERASURE) != 0;
        stats->repeat += (result.flags & MBE_PROCESS_FLAG_REPEAT) != 0;
        stats->mute += (result.flags & MBE_PROCESS_FLAG_MUTE) != 0;
        stats->c0 += (size_t)result.c0_errors;
        stats->protected_errors += (size_t)result.protected_errors;
        stats->c4 += (size_t)result.c4_errors;
        stats->total += (size_t)result.total_errors;
        /* Inspect object bits, not isfinite(): Release may assume finite math.
         * These are post-library-limiter, pre-int16 observations. */
        for (size_t i = 0; i < 160; ++i) {
            uint32_t bits_value;
            memcpy(&bits_value, &audio[i], sizeof(bits_value));
            if ((bits_value & UINT32_C(0x7f800000)) == UINT32_C(0x7f800000)) {
                ++stats->nonfinite;
                continue;
            }
            double magnitude = fabs((double)audio[i]);
            stats->peak = fmax(stats->peak, magnitude);
            stats->clips += magnitude >= (double)((32767.0f * 0.95f) / 7.0f);
        }
        if (stats->nonfinite) {
            fprintf(stderr, "frame %zu: nonfinite float PCM (%zu samples)\n", frame, stats->nonfinite);
            exit(3);
        }
        mbe_floattoshort(audio, pcm);
        for (size_t i = 0; i < 160; ++i) {
            append(&s, pcm[i]);
        }
        ++frame;
    }
    if (ferror(f)) {
        fail("frame file read failed");
    }
    fclose(f);
    if (!frame) {
        fail("frame file is empty");
    }
    return s;
}

static double
energy(const double* x, size_t n) {
    double e = 0;
    for (size_t i = 0; i < n; ++i) {
        e += x[i] * x[i];
    }
    return e;
}

/* Centered 20 ms RMS envelopes at a 1 ms hop, with truncated edge windows. */
static Signal
alignment_envelope(const Signal* s) {
    Signal e = {0};
    double maximum = 0;
    for (size_t i = 0; i < s->count; i += 8) {
        size_t start = i > 80 ? i - 80 : 0;
        size_t end = i + 80 < s->count ? i + 80 : s->count;
        double v = sqrt(energy(s->samples + start, end - start) / (double)(end - start));
        append(&e, v);
        maximum = fmax(maximum, v);
    }
    double floor_value = fmax(maximum * 1e-4, 1e-12);
    for (size_t i = 0; i < e.count; ++i) {
        e.samples[i] = log10(fmax(e.samples[i], floor_value));
    }
    return e;
}

static int
correlation(const double* a, const double* b, size_t n, double* value) {
    double ax = 0, bx = 0, aa = 0, bb = 0, ab = 0;
    *value = 0;
    if (!n) {
        return 0;
    }
    for (size_t i = 0; i < n; ++i) {
        ax += a[i];
        bx += b[i];
    }
    ax /= (double)n;
    bx /= (double)n;
    for (size_t i = 0; i < n; ++i) {
        double x = a[i] - ax, y = b[i] - bx;
        aa += x * x;
        bb += y * y;
        ab += x * y;
    }
    if (aa <= 1e-20 || bb <= 1e-20) {
        return 0;
    }
    *value = ab / sqrt(aa * bb);
    return 1;
}

static Alignment
align(const Signal* ref, const Signal* dec) {
    Signal a = alignment_envelope(ref), b = alignment_envelope(dec);
    Alignment best = {0};
    for (int lag = -20; lag <= 100; ++lag) {
        size_t ra = lag < 0 ? (size_t)-lag : 0;
        size_t rb = lag > 0 ? (size_t)lag : 0;
        if (ra >= a.count || rb >= b.count) {
            continue;
        }
        size_t n = a.count - ra < b.count - rb ? a.count - ra : b.count - rb;
        if (n < 20) {
            continue;
        }
        double c;
        if (correlation(a.samples + ra, b.samples + rb, n, &c)
            && (!best.valid || c > best.corr || (c == best.corr && abs(lag * 8) < abs(best.lag)))) {
            best = (Alignment){lag * 8, c, 1};
        }
    }
    free(a.samples);
    free(b.samples);
    return best;
}

static void
metric_valid(const char* key, double value, int valid) {
    if (metric_count >= sizeof(metrics) / sizeof(metrics[0])) {
        fail("Metric capacity exceeded.");
    }
    metrics[metric_count++] = (Metric){key, value, valid};
}

static void
metric(const char* key, double value) {
    metric_valid(key, value, 1);
}

static double
db_ratio(double numerator, double denominator) {
    return 10 * log10(fmax(numerator, 1e-20) / fmax(denominator, 1e-20));
}

static double
boundary_index(const Signal* s) {
    double boundary = 0, middle = 0;
    for (size_t b = 160; b + 83 < s->count; b += 160) {
        for (size_t n = b - 4; n <= b + 3; ++n) {
            double d = s->samples[n] - s->samples[n - 1];
            double m = s->samples[n + 80] - s->samples[n + 79];
            boundary += d * d;
            middle += m * m;
        }
    }
    return db_ratio(boundary, middle);
}

/* Localized join contrast: the across-frame derivative versus its eight
 * neighbors. Keep the actual synthesis grid after alignment, not a new grid
 * starting at the trimmed buffer. Reference and decoded use the same mask. */
static double
join_contrast(const double* x, size_t n, size_t decoded_offset, const unsigned char* active, size_t nf, size_t* count) {
    double at_join = 0, neighbors = 0;
    *count = 0;
    for (size_t b = 160 - decoded_offset % 160; b + 4 < n; b += 160) {
        if (b <= 4 || b / 160 >= nf || !active[b / 160]) {
            continue;
        }
        ++*count;
        for (size_t i = b - 4; i <= b + 4; ++i) {
            double delta = x[i] - x[i - 1];
            if (i == b) {
                at_join += delta * delta;
            } else {
                neighbors += delta * delta;
            }
        }
    }
    return db_ratio(at_join, neighbors / 8);
}

/* Deliberately independent of the vocoder FFT. Unnormalized, symmetric Hann. */
static void
spectrum(const double* x, double power[129]) {
    double re[FFT_N], im[FFT_N] = {0};
    for (unsigned i = 0, j = 0; i < FFT_N; ++i) {
        re[j] = x[i] * (0.5 - 0.5 * cos(2 * PI * (double)i / (FFT_N - 1)));
        unsigned bit = FFT_N / 2;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
    }
    for (unsigned len = 2; len <= FFT_N; len *= 2) {
        double wr = cos(-2 * PI / len), wi = sin(-2 * PI / len);
        for (unsigned i = 0; i < FFT_N; i += len) {
            double ur = 1, ui = 0;
            for (unsigned j = 0; j < len / 2; ++j) {
                unsigned a = i + j, b = a + len / 2;
                double vr = re[b] * ur - im[b] * ui, vi = re[b] * ui + im[b] * ur;
                re[b] = re[a] - vr;
                im[b] = im[a] - vi;
                re[a] += vr;
                im[a] += vi;
                double next = ur * wr - ui * wi;
                ui = ur * wi + ui * wr;
                ur = next;
            }
        }
    }
    for (size_t i = 0; i <= FFT_N / 2; ++i) {
        power[i] = re[i] * re[i] + im[i] * im[i];
    }
}

static double
crest(const double* x) {
    double peak = 0;
    for (size_t i = 0; i < 160; ++i) {
        peak = fmax(peak, fabs(x[i]));
    }
    return db_ratio(peak * peak, energy(x, 160) / 160);
}

static void
measure(const Signal* dec, const Signal* reference, int forced_lag, int have_lag) {
    double pcm_peak = 0;
    size_t rails = 0, run = 0, max_run = 0;
    for (size_t i = 0; i < dec->count; ++i) {
        double magnitude = fabs(dec->samples[i]);
        pcm_peak = fmax(pcm_peak, magnitude);
        if (magnitude >= 31128) {
            ++rails;
            ++run;
            if (run > max_run) {
                max_run = run;
            }
        } else {
            run = 0;
        }
    }
    metric("pcm_peak", pcm_peak);
    metric("pcm_rail_samples", (double)rails);
    metric("pcm_max_rail_run", (double)max_run);
    metric_valid("boundary_index_db", boundary_index(dec), dec->count > 243 && pcm_peak > 0);
    Alignment automatic = reference ? align(reference, dec) : (Alignment){0};
    int lag = have_lag ? forced_lag : automatic.lag;
    metric("auto_lag_samples", automatic.lag);
    metric("lag_samples", lag);
    metric_valid("alignment_corr", automatic.corr, automatic.valid);
    metric("alignment_at_limit", automatic.valid && (automatic.lag == -160 || automatic.lag == 800));
    size_t r0 = reference && lag < 0 ? (size_t)-lag : 0, d0 = lag > 0 ? (size_t)lag : 0;
    if (d0 >= dec->count || (reference && r0 >= reference->count)) {
        fail("need at least 256 aligned samples for spectral metrics");
    }
    size_t n = dec->count - d0;
    if (reference && reference->count - r0 < n) {
        n = reference->count - r0;
    }
    if (n < FFT_N) {
        fail("need at least 256 aligned samples for spectral metrics");
    }
    metric("reference_samples", reference ? (double)reference->count : 0);
    metric("decoded_samples", (double)dec->count);
    metric("aligned_samples", (double)n);
    metric("reference_trim_head", (double)r0);
    metric("reference_trim_tail", reference ? (double)(reference->count - r0 - n) : 0);
    metric("decoded_trim_head", (double)d0);
    metric("decoded_trim_tail", (double)(dec->count - d0 - n));
    const double* d = dec->samples + d0;
    const double* r = reference ? reference->samples + r0 : d;
    size_t nf = n / 160;
    unsigned char* active = calloc(nf, 1);
    if (!active) {
        fail("out of memory");
    }
    double maximum = 0;
    for (size_t f = 0; f < nf; ++f) {
        maximum = fmax(maximum, energy(r + f * 160, 160));
    }
    double er = 0, ed = 0;
    size_t count = 0;
    for (size_t f = 0; f < nf; ++f) {
        double e = energy(r + f * 160, 160);
        if (e >= maximum * 1e-4 && e > 0) {
            active[f] = 1;
            er += e;
            ed += energy(d + f * 160, 160);
            ++count;
        }
    }
    int supported = er > 0 && ed > 0;
    int speech = reference && supported;
    reference_state = !reference ? "no_reference"
                      : er <= 0  ? "silent_reference"
                      : ed <= 0  ? "silent_decoded"
                                 : "speech";
    metric("active_frames", (double)count);
    double gain = speech ? sqrt(er / ed) : 1;
    metric_valid("level_offset_db", -20 * log10(gain), speech);
    double* normalized = NULL;
    if (speech && gain != 1) {
        normalized = malloc(n * sizeof(*normalized));
        if (!normalized) {
            fail("out of memory");
        }
        for (size_t i = 0; i < n; ++i) {
            normalized[i] = d[i] * gain;
        }
        d = normalized;
    }
    size_t joins;
    double join_dec = join_contrast(d, n, d0, active, nf, &joins);
    metric("join_count", (double)joins);
    metric_valid("join_dec_db", join_dec, supported && joins > 0);
    double join_ref = reference ? join_contrast(r, n, d0, active, nf, &joins) : 0;
    metric_valid("join_ref_db", join_ref, speech && joins > 0);
    metric_valid("join_excess_db", join_dec - join_ref, speech && joins > 0);
    double cr = 0, cd = 0;
    Signal envr = {0}, envd = {0};
    double maxr = 0, maxd = 0;
    for (size_t f = 0; f < nf; ++f) {
        if (!active[f]) {
            continue;
        }
        cr += crest(r + f * 160);
        cd += crest(d + f * 160);
        for (size_t j = 0; j < 160; j += 40) {
            double vr = sqrt(energy(r + f * 160 + j, 40) / 40);
            double vd = sqrt(energy(d + f * 160 + j, 40) / 40);
            append(&envr, vr);
            append(&envd, vd);
            maxr = fmax(maxr, vr);
            maxd = fmax(maxd, vd);
        }
    }
    for (size_t i = 0; i < envr.count; ++i) {
        envr.samples[i] = 20 * log10(fmax(envr.samples[i], fmax(maxr * 1e-3, 1e-12)));
        envd.samples[i] = 20 * log10(fmax(envd.samples[i], fmax(maxd * 1e-3, 1e-12)));
    }
    double env_corr;
    int env_valid = correlation(envr.samples, envd.samples, envr.count, &env_corr);
    metric_valid("env_corr", env_corr, speech && env_valid);
    metric_valid("crest_ref_db", count ? cr / (double)count : 0, speech);
    metric_valid("crest_delta_db", count ? (cd - cr) / (double)count : 0, speech);
    metric_valid("crest_dec_db", count ? cd / (double)count : 0, supported);
    free(envr.samples);
    free(envd.samples);

    double max_power = 0, pr[129], pd[129];
    for (size_t i = 0; i + FFT_N <= n; i += 80) {
        spectrum(r + i, pr);
        for (size_t k = 0; k <= 128; ++k) {
            max_power = fmax(max_power, pr[k]);
        }
    }
    double epsilon = fmax(max_power * 1e-6, 1e-20);
    double br[5] = {0}, bd[5] = {0}, lsd = 0;
    size_t spectra = 0;
    static const size_t edges[] = {0, 16, 32, 64, 96, 129};
    for (size_t i = 0; i + FFT_N <= n; i += 80) {
        size_t frame = (i + FFT_N / 2) / 160;
        if (frame >= nf || !active[frame]) {
            continue;
        }
        spectrum(r + i, pr);
        spectrum(d + i, pd);
        double error = 0;
        for (size_t k = 2; k <= 118; ++k) {
            double delta = db_ratio(pr[k] + epsilon, pd[k] + epsilon);
            error += delta * delta;
        }
        lsd += sqrt(error / 117);
        for (size_t band = 0; band < 5; ++band) {
            for (size_t k = edges[band]; k < edges[band + 1]; ++k) {
                br[band] += pr[k];
                bd[band] += pd[k];
            }
        }
        ++spectra;
    }
    metric("spectral_frames", (double)spectra);
    metric_valid("lsd_db", spectra ? lsd / (double)spectra : 0, speech && spectra > 0);
    static const char* delta_names[] = {"band_delta_db_0_500", "band_delta_db_500_1000", "band_delta_db_1000_2000",
                                        "band_delta_db_2000_3000", "band_delta_db_3000_4000"};
    static const char* absolute_names[] = {"band_db_0_500", "band_db_500_1000", "band_db_1000_2000",
                                           "band_db_2000_3000", "band_db_3000_4000"};
    for (size_t band = 0; band < 5; ++band) {
        /* Absolute levels: mean STFT band power in unnormalized int16 units. */
        metric_valid(delta_names[band], db_ratio(bd[band], br[band]), speech && spectra > 0);
        if (!reference) {
            metric_valid(absolute_names[band], db_ratio(bd[band], (double)(spectra ? spectra : 1)),
                         supported && spectra > 0);
        }
    }
    free(active);
    free(normalized);
}

int
main(int argc, char** argv) {
    const char *codec = NULL, *frames = NULL, *out = NULL, *ref = NULL, *json = NULL, *pcm = NULL;
    uint32_t seed = 0x12345678u;
    int have_seed = 0, have_lag = 0, lag = 0;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 >= argc || !*argv[i + 1]) {
            fail("options require nonempty --key value pairs");
        }
        for (int j = 1; j < i; j += 2) {
            if (!strcmp(argv[i], argv[j])) {
                fail("duplicate option");
            }
        }
        const char* v = argv[i + 1];
        if (!strcmp(argv[i], "--codec")) {
            codec = v;
        } else if (!strcmp(argv[i], "--frames")) {
            frames = v;
        } else if (!strcmp(argv[i], "--out")) {
            out = v;
        } else if (!strcmp(argv[i], "--ref")) {
            ref = v;
        } else if (!strcmp(argv[i], "--json")) {
            json = v;
        } else if (!strcmp(argv[i], "--decoded")) {
            pcm = v;
        } else if (!strcmp(argv[i], "--seed")) {
            char* end;
            errno = 0;
            unsigned long value = strtoul(v, &end, 0);
            if (errno || (*v != '+' && (*v < '0' || *v > '9')) || *end || value > UINT32_MAX) {
                fail("seed must be a uint32");
            }
            seed = (uint32_t)value;
            have_seed = 1;
        } else if (!strcmp(argv[i], "--lag-samples")) {
            char* end;
            errno = 0;
            long value = strtol(v, &end, 10);
            if (errno || (*v != '-' && *v != '+' && (*v < '0' || *v > '9')) || *end || value < -160 || value > 800) {
                fail("lag-samples must be an integer in -160..800");
            }
            lag = (int)value;
            have_lag = 1;
        } else {
            fail("unknown option");
        }
    }
    if (pcm) {
        if (codec || frames || out || have_seed) {
            fail("--decoded rejects --codec, --frames, --out, and --seed");
        }
    } else if (!codec || !frames || !out
               || (strcmp(codec, "imbe7200") != 0 && strcmp(codec, "imbe7100") != 0 && strcmp(codec, "ambe2400") != 0
                   && strcmp(codec, "ambe2450") != 0)) {
        fail("usage: mbe_quality_eval --codec imbe7200|imbe7100|ambe2400|ambe2450 --frames file --out decoded.wav "
             "[--ref speech.raw|speech.wav] [--seed uint32] [--lag-samples -160..800] [--json metrics.json]\n"
             "   or: mbe_quality_eval --decoded input.raw|input.wav [--ref speech.raw|speech.wav] "
             "[--lag-samples -160..800] [--json metrics.json]");
    }
    reject_alias(out, frames);
    reject_alias(out, ref);
    reject_alias(json, frames);
    reject_alias(json, ref);
    reject_alias(json, pcm);
    reject_alias(out, json);
    DecodeStats stats = {0};
    Signal decoded = pcm ? read_reference(pcm) : decode(codec, frames, seed, &stats);
    Signal reference = ref ? read_reference(ref) : (Signal){0};
    measure(&decoded, ref ? &reference : NULL, lag, have_lag);
    metric_valid("float_nonfinite_samples", (double)stats.nonfinite, !pcm);
    metric_valid("float_peak", stats.peak, !pcm);
    metric_valid("float_clip_samples", (double)stats.clips, !pcm);
    metric_valid("tone_frames", (double)stats.tone, !pcm);
    metric_valid("erasure_frames", (double)stats.erasure, !pcm);
    metric_valid("repeat_frames", (double)stats.repeat, !pcm);
    metric_valid("mute_frames", (double)stats.mute, !pcm);
    metric_valid("c0_errors", (double)stats.c0, !pcm);
    metric_valid("protected_errors", (double)stats.protected_errors, !pcm);
    metric_valid("c4_errors", (double)stats.c4, !pcm);
    metric_valid("total_errors", (double)stats.total, !pcm);
    /* Measurement never mutates decoded; WAV and hash retain the exact PCM. */
    uint32_t hash = 2166136261u;
    if (out) {
        hash = write_wav(out, &decoded);
    } else {
        for (size_t i = 0; i < decoded.count; ++i) {
            uint16_t v = (uint16_t)(int16_t)decoded.samples[i];
            hash = (hash ^ (v & 255u)) * 16777619u;
            hash = (hash ^ (v >> 8)) * 16777619u;
        }
    }
    const char* input_api = pcm ? "pcm" : stats.input_api;
    printf("schema_version=2\nreference_state=%s\ninput_api=%s\nframes=%zu pcm_fnv1a=0x%08X\n", reference_state,
           input_api, decoded.count / 160, (unsigned)hash);
    FILE* jf = json ? open_file(json, "w") : NULL;
    if (jf) {
        fprintf(jf,
                "{\n  \"schema_version\": 2,\n  \"reference_state\": \"%s\",\n  \"input_api\": \"%s\","
                "\n  \"frames\": %zu,\n  \"pcm_fnv1a\": \"0x%08X\"",
                reference_state, input_api, decoded.count / 160, (unsigned)hash);
    }
    for (size_t i = 0; i < metric_count; ++i) {
        if (metrics[i].valid) {
            printf("%s=%.17g\n", metrics[i].key, metrics[i].value);
        } else {
            printf("%s=null\n", metrics[i].key);
        }
        if (jf) {
            fprintf(jf, ",\n  \"%s\": ", metrics[i].key);
            if (!metrics[i].valid) {
                fputs("null", jf);
            } else if (!strcmp(metrics[i].key, "alignment_at_limit")) {
                fputs(metrics[i].value ? "true" : "false", jf);
            } else {
                fprintf(jf, "%.17g", metrics[i].value);
            }
        }
    }
    if (jf) {
        fputs("\n}\n", jf);
        if (ferror(jf) || fclose(jf)) {
            fail("JSON write failed");
        }
    }
    free(decoded.samples);
    free(reference.samples);
    return 0;
}
