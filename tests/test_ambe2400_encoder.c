// SPDX-License-Identifier: GPL-2.0-or-later
/* AMBE 3600x2400 encoder tests: frame FEC, DV byte packing, Golay, and
 * encoder->decoder parameter-state parity. */
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include "mbelib-neo/mbelib.h"

#ifdef MBE_ENCODER_TEST_OOM
#include <string.h>

/* GNU link wrapping covers eleven allocation points: the context, the
 * project-owned synthesis FFT plan and its buffers, and the autocorrelation
 * plan and its buffers. These failures return NULL. Allocation failures inside
 * the vendored pffft setup are not recoverable; the same limitation applies to
 * decoder plan allocation. PFFFT's same-object internal aligned calls are not
 * wrapped. */
#define ENCODER_ALLOCATIONS 11
void* encoder_real_alloc(size_t size) __asm__("__real_pffft_aligned_malloc");
void encoder_real_free(void* ptr) __asm__("__real_pffft_aligned_free");
void* encoder_real_calloc(size_t count, size_t size) __asm__("__real_calloc");
void* encoder_real_malloc(size_t size) __asm__("__real_malloc");
void encoder_real_heap_free(void* ptr) __asm__("__real_free");
extern void* encoder_fault_alloc(size_t size) __asm__("__wrap_pffft_aligned_malloc");
extern void encoder_track_free(void* ptr) __asm__("__wrap_pffft_aligned_free");
extern void* encoder_fault_calloc(size_t count, size_t size) __asm__("__wrap_calloc");
extern void* encoder_track_malloc(size_t size) __asm__("__wrap_malloc");
extern void encoder_track_heap_free(void* ptr) __asm__("__wrap_free");
static int allocation_count;
static int fail_at;
static int live_allocations;
static int heap_allocation_count;
static void* live_callocs[2];

extern void*
encoder_track_malloc(size_t size) {
    heap_allocation_count++;
    return encoder_real_malloc(size);
}

extern void*
encoder_fault_calloc(size_t count, size_t size) {
    allocation_count++;
    if (allocation_count == fail_at) {
        return NULL;
    }
    void* ptr = encoder_real_calloc(count, size);
    for (size_t i = 0; ptr != NULL && i < sizeof(live_callocs) / sizeof(live_callocs[0]); i++) {
        if (live_callocs[i] == NULL) {
            live_callocs[i] = ptr;
            break;
        }
    }
    return ptr;
}

static int
live_calloc_count(void) {
    int count = 0;
    for (size_t i = 0; i < sizeof(live_callocs) / sizeof(live_callocs[0]); i++) {
        count += live_callocs[i] != NULL;
    }
    return count;
}

extern void
encoder_track_heap_free(void* ptr) {
    for (size_t i = 0; ptr != NULL && i < sizeof(live_callocs) / sizeof(live_callocs[0]); i++) {
        if (live_callocs[i] == ptr) {
            live_callocs[i] = NULL;
        }
    }
    encoder_real_heap_free(ptr);
}

extern void*
encoder_fault_alloc(size_t size) {
    allocation_count++;
    if (allocation_count == fail_at) {
        return NULL;
    }
    void* ptr = encoder_real_alloc(size);
    if (ptr != NULL) {
        live_allocations++;
    }
    return ptr;
}

extern void
encoder_track_free(void* ptr) {
    if (ptr != NULL) {
        live_allocations--;
    }
    encoder_real_free(ptr);
}

/* The three encoders share the analysis front end and its allocations. */
struct oom_codec {
    const char* name;
    void* (*alloc)(void);
    void (*reset)(void* enc);
    void (*release)(void* enc);
    int (*encode)(void* enc, const float* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev);
    int (*encode_short)(void* enc, const short* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev);
};

static void*
dstar_alloc(void) {
    return mbe_ambe2400EncoderAlloc();
}

static void
dstar_reset(void* enc) {
    mbe_ambe2400EncoderReset(enc);
}

static void
dstar_release(void* enc) {
    mbe_ambe2400EncoderFree(enc);
}

static int
dstar_encode(void* enc, const float* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2400Parms(enc, pcm, bits, cur, prev);
}

static int
dstar_encode_short(void* enc, const short* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2400ParmsShort(enc, pcm, bits, cur, prev);
}

static void*
ambe2450_alloc(void) {
    return mbe_ambe2450EncoderAlloc();
}

static void
ambe2450_reset(void* enc) {
    mbe_ambe2450EncoderReset(enc);
}

static void
ambe2450_release(void* enc) {
    mbe_ambe2450EncoderFree(enc);
}

static int
ambe2450_encode(void* enc, const float* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2450Parms(enc, pcm, bits, cur, prev);
}

static int
ambe2450_encode_short(void* enc, const short* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeAmbe2450ParmsShort(enc, pcm, bits, cur, prev);
}

static void*
imbe_alloc(void) {
    return mbe_imbe4400EncoderAlloc();
}

static void
imbe_reset(void* enc) {
    mbe_imbe4400EncoderReset(enc);
}

static void
imbe_release(void* enc) {
    mbe_imbe4400EncoderFree(enc);
}

static int
imbe_encode(void* enc, const float* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeImbe4400Parms(enc, pcm, bits, cur, prev);
}

static int
imbe_encode_short(void* enc, const short* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev) {
    return mbe_encodeImbe4400ParmsShort(enc, pcm, bits, cur, prev);
}

static const struct oom_codec oom_codecs[] = {
    {"ambe2400", dstar_alloc, dstar_reset, dstar_release, dstar_encode, dstar_encode_short},
    {"ambe2450", ambe2450_alloc, ambe2450_reset, ambe2450_release, ambe2450_encode, ambe2450_encode_short},
    {"imbe4400", imbe_alloc, imbe_reset, imbe_release, imbe_encode, imbe_encode_short},
};

static int
check_no_encode_allocations(const struct oom_codec* codec, void* enc) {
    mbe_parms cur, prev, enhanced;
    float pcm[160];
    char bits[88];
    mbe_initMbeParms(&cur, &prev, &enhanced);
    int before = allocation_count;
    int heap_before = heap_allocation_count;
    fail_at = before + 1;
    for (int frame = 0; frame < 12; frame++) {
        for (int i = 0; i < 160; i++) {
            pcm[i] = frame < 6 ? 0.01f * sinf(0.1f * (float)i) : 0.0f;
            pcm[i] += (frame >= 8 && frame < 10) ? 0.1f * sinf(0.55f * (float)i) : 0.0f; /* a tone frame */
        }
        if (codec->encode(enc, pcm, bits, &cur, &prev) < 0 || allocation_count != before
            || heap_allocation_count != heap_before) {
            return 1;
        }
        mbe_moveMbeParms(&cur, &prev);
    }
    codec->reset(enc);
    short shorts[160] = {0};
    return codec->encode_short(enc, shorts, bits, &cur, &prev) < 0 || allocation_count != before
           || heap_allocation_count != heap_before;
}

/* Usage: test_ambe2400_encoder_oom ALLOCATION [ambe2400|ambe2450|imbe4400] */
int
main(int argc, char** argv) {
    const struct oom_codec* codec = &oom_codecs[0];
    fail_at = 0;
    for (const char* p = argc >= 2 ? argv[1] : ""; *p != '\0' && fail_at <= ENCODER_ALLOCATIONS; p++) {
        fail_at = *p >= '0' && *p <= '9' ? (fail_at * 10) + (*p - '0') : ENCODER_ALLOCATIONS + 1;
    }
    for (size_t i = 0; argc == 3 && i < sizeof(oom_codecs) / sizeof(oom_codecs[0]); i++) {
        codec = strcmp(argv[2], oom_codecs[i].name) == 0 ? &oom_codecs[i] : codec;
    }
    if (fail_at < 1 || fail_at > ENCODER_ALLOCATIONS || argc > 3 || (argc == 3 && strcmp(argv[2], codec->name) != 0)) {
        return 1;
    }
    void* enc = codec->alloc();
    if (enc != NULL || allocation_count < fail_at || live_allocations != 0 || live_calloc_count() != 0) {
        (void)fprintf(stderr, "%s: failed allocation %d: calls=%d live=%d callocs=%d\n", codec->name, fail_at,
                      allocation_count, live_allocations, live_calloc_count());
        codec->release(enc);
        return 1;
    }
    fail_at = 0;
    allocation_count = 0;
    enc = codec->alloc();
    if (enc == NULL || allocation_count != ENCODER_ALLOCATIONS) {
        (void)fprintf(stderr, "%s: encoder allocation points: %d, want %d\n", codec->name, allocation_count,
                      ENCODER_ALLOCATIONS);
        codec->release(enc);
        return 1;
    }
    int failed = check_no_encode_allocations(codec, enc);
    codec->release(enc);
    if (failed || live_allocations != 0 || live_calloc_count() != 0) {
        return 1;
    }
    (void)printf("%s allocation failure: NULL, no leaks; retry, reset and encoding: no allocations\n", codec->name);
    return 0;
}
#else

#include <stdint.h>
#include <string.h>

#include "ambe_common.h"
#include "mbe_ecc.h"

static int
float_bits_equal(float a, float b) {
    uint32_t a_bits, b_bits;
    memcpy(&a_bits, &a, sizeof(a_bits));
    memcpy(&b_bits, &b, sizeof(b_bits));
    return a_bits == b_bits;
}

/* Bitwise identity of two parameter sets, compared as bytes. */
static int
parms_identical(const mbe_parms* a, const mbe_parms* b) {
    unsigned char x[sizeof(mbe_parms)];
    unsigned char y[sizeof(mbe_parms)];
    memcpy(x, a, sizeof(x));
    memcpy(y, b, sizeof(y));
    return memcmp(x, y, sizeof(x)) == 0;
}

/* At most 2^127 in magnitude: a range test on the bit pattern read from
 * memory, which unlike an exponent-mask test or a by-value float survives
 * fast-math. */
static int
float_finite(const float* x) {
    uint32_t bits;
    memcpy(&bits, x, sizeof(bits));
    return (bits & 0x7FFFFFFFu) <= 0x7F000000u;
}

/* w0, gamma, and log2Ml and Ml over the model's harmonics are finite. */
static int
model_finite(const mbe_parms* mp) {
    int ok = float_finite(&mp->w0) && float_finite(&mp->gamma) && mp->L >= 1 && mp->L <= 56;
    for (int l = 1; ok && l <= mp->L; l++) {
        ok = float_finite(&mp->log2Ml[l]) && float_finite(&mp->Ml[l]);
    }
    return ok;
}

static uint32_t rng = 0xC0FFEE;

static uint32_t
rnd(void) {
    rng = rng * 1664525u + 1013904223u;
    return rng >> 8;
}

static int
test_frame_roundtrip(void) {
    int fails = 0;
    for (int it = 0; it < 4096; it++) {
        char d[49], fr[4][24], out[49];
        mbe_process_result res;
        for (int i = 0; i < 49; i++) {
            d[i] = (char)(rnd() & 1);
        }
        if (mbe_encodeAmbe3600x2400Frame(d, fr) != 0) {
            fails++;
            continue;
        }
        int errs = mbe_decodeAmbe3600x2400Frame((const char (*)[24])fr, out, &res);
        int mism = 0;
        for (int i = 0; i < 49; i++) {
            if (i != 24 && d[i] != out[i]) {
                mism++;
            }
        }
        if (errs != 0 || mism) {
            fails++;
            if (fails < 5) {
                printf("  frame rt: errs=%d mism=%d\n", errs, mism);
            }
        }
    }
    printf("frame FEC roundtrip: %s (%d fails / 4096)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

static int
test_frame_single_bit_errors(void) {
    for (int it = 0; it < 128; it++) {
        char d[49], fr[4][24], damaged[4][24], out[49];
        for (int i = 0; i < 49; i++) {
            d[i] = (char)(rnd() & 1);
        }
        if (mbe_encodeAmbe3600x2400Frame(d, fr) != 0) {
            return 1;
        }
        for (int plane = 0; plane < 2; plane++) {
            for (int bit = 0; bit < 12; bit++) {
                memcpy(damaged, fr, sizeof(fr));
                damaged[plane][12 - plane + bit] ^= 1;
                if (mbe_decodeAmbe3600x2400Frame((const char (*)[24])damaged, out, NULL) != 1) {
                    return 1;
                }
                for (int i = 0; i < 49; i++) {
                    if (i != 24 && d[i] != out[i]) {
                        return 1;
                    }
                }
            }
        }
    }
    puts("single-bit data error correction: ok (all C0/C1 data positions)");
    return 0;
}

static int
test_dv_bytes(void) {
    int fails = 0;
    for (int it = 0; it < 4096; it++) {
        char fr[4][24], back[4][24];
        unsigned char b[9], b2[9];
        for (int p = 0; p < 4; p++) {
            for (int i = 0; i < 24; i++) {
                fr[p][i] = (char)(rnd() & 1);
            }
        }
        mbe_encodeDStarDVData((const char (*)[24])fr, b);
        mbe_decodeDStarDVData(b, back);
        mbe_encodeDStarDVData((const char (*)[24])back, b2);
        if (memcmp(b, b2, 9) != 0) {
            fails++;
        }
        const int carried[4] = {24, 23, 11, 14};
        for (int p = 0; p < 4; p++) {
            for (int i = 0; i < carried[p]; i++) {
                if (fr[p][i] != back[p][i]) {
                    fails++;
                }
            }
        }
    }
    /* permutation sanity: every (plane,idx) hit exactly once */
    {
        char fr[4][24];
        unsigned char b[9];
        int hits = 0;
        int used_ok = 1;
        for (int p = 0; p < 4; p++) {
            for (int i = 0; i < 24; i++) {
                memset(fr, 0, sizeof fr);
                fr[p][i] = 1;
                mbe_encodeDStarDVData((const char (*)[24])fr, b);
                int ones = 0;
                for (int k = 0; k < 9; k++) {
                    for (int q = 0; q < 8; q++) {
                        ones += (b[k] >> q) & 1;
                    }
                }
                if (ones == 1) {
                    hits++;
                }
                int expect = (p == 0) || (p == 1 && i <= 22) || (p == 2 && i <= 10) || (p == 3 && i <= 13);
                if (ones != expect) {
                    used_ok = 0;
                    printf("  dv: plane %d idx %d reaches air %d times (expected %d)\n", p, i, ones, expect);
                }
            }
        }
        printf("  dv: %d of 96 plane slots carried (72 expected), used-set matches decoder layout: %s\n", hits,
               used_ok ? "yes" : "NO");
        if (hits != 72 || !used_ok) {
            fails++;
        }
    }
    printf("DV byte pack/unpack inverse: %s (%d fails)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

static int
test_golay(void) {
    int fails = 0;
    for (int it = 0; it < 4096; it++) {
        char in[12], cw[23], out[23];
        for (int i = 0; i < 12; i++) {
            in[i] = (char)((it >> (11 - i)) & 1);
        }
        mbe_golay2312_encode(in, cw);
        int errs = mbe_golay2312(cw, out);
        int mism = 0;
        for (int i = 0; i < 12; i++) {
            if (out[22 - i] != in[i]) {
                mism++;
            }
        }
        if (errs != 0 || mism || memcmp(cw, out, 23) != 0) {
            fails++;
        }
    }
    printf("golay encode->decode: %s (%d fails / 4096)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

/* Synthetic speech-like test signal: harmonic tone with pitch glide, bursts
 * of noise, and silences. */
static void
gen_signal(short* pcm, int n) {
    double ph = 0;
    for (int i = 0; i < n; i++) {
        int seg = (i / 1600) % 6; /* 200 ms segments */
        double t = (double)i / 8000.0;
        double v = 0;
        if (seg == 0 || seg == 1 || seg == 3) {
            double f0 = 120.0 + 60.0 * sin(2 * M_PI * 0.7 * t);
            ph += 2 * M_PI * f0 / 8000.0;
            for (int h = 1; h <= 20; h++) {
                v += sin(h * ph) / (h * 0.8 + 0.5);
            }
            v *= 0.25 * (0.6 + 0.4 * sin(2 * M_PI * 3.0 * t));
        } else if (seg == 2) {
            v = 0.05 * (((double)(rnd() & 0xffff) / 32768.0) - 1.0);
        } else if (seg == 4) {
            v = 0.0; /* silence */
        } else {
            double f0 = 200.0;
            ph += 2 * M_PI * f0 / 8000.0;
            for (int h = 1; h <= 10; h++) {
                v += sin(h * ph) / h;
            }
            v *= 0.15;
        }
        if (v > 0.99) {
            v = 0.99;
        }
        if (v < -0.99) {
            v = -0.99;
        }
        pcm[i] = (short)(v * 32767.0);
    }
}

static int
test_state_parity(mbe_ambe2400_encoder* enc) {
    mbe_ambe2400EncoderReset(enc);

    enum { FRAMES = 600 };

    static short pcm[FRAMES * 160];
    gen_signal(pcm, FRAMES * 160);

    mbe_parms e_cur, e_prev, e_enh; /* encoder chain */
    mbe_parms d_cur, d_prev, d_enh; /* independent decoder chain */
    mbe_initMbeParms(&e_cur, &e_prev, &e_enh);
    mbe_initMbeParms(&d_cur, &d_prev, &d_enh);

    int bad_frames = 0, voice = 0, first_bad = -1;
    float worst = 0;
    int worst_frame = -1, worst_l = -1;
    int L_mism = 0, Vl_mism = 0, gamma_mism = 0, w0_mism = 0, log2_mism = 0, Ml_mism = 0;
    for (int f = 0; f < FRAMES; f++) {
        char d[49];
        unsigned char saved_prev[sizeof(e_prev)];
        memcpy(saved_prev, &e_prev, sizeof(e_prev));
        int r = mbe_encodeAmbe2400ParmsShort(enc, pcm + (size_t)f * 160, d, &e_cur, &e_prev);
        unsigned char after_prev[sizeof(e_prev)];
        memcpy(after_prev, &e_prev, sizeof(e_prev));
        if (memcmp(saved_prev, after_prev, sizeof(e_prev)) != 0) {
            return 1;
        }
        if (r != 0) {
            printf("  encode returned %d at frame %d\n", r, f);
            return 1;
        }
        int dr = mbe_decodeAmbe2400Parms(d, &d_cur, &d_prev);
        voice++;
        if (dr != 0) {
            printf("  frame %d: decoder returned %d for a voice frame\n", f, dr);
            bad_frames++;
        }
        int frame_bad = 0;
        unsigned char decoder_w0[sizeof(d_cur.w0)], encoder_w0[sizeof(e_cur.w0)];
        memcpy(decoder_w0, &d_cur.w0, sizeof(decoder_w0));
        memcpy(encoder_w0, &e_cur.w0, sizeof(encoder_w0));
        if (memcmp(decoder_w0, encoder_w0, sizeof(decoder_w0)) != 0) {
            w0_mism++;
            frame_bad = 1;
        }
        if (d_cur.L != e_cur.L) {
            L_mism++;
            frame_bad = 1;
        }
        if (!float_bits_equal(d_cur.gamma, e_cur.gamma)) {
            gamma_mism++;
            frame_bad = 1;
        }
        for (int l = 1; l <= d_cur.L && l <= 56; l++) {
            if (d_cur.Vl[l] != e_cur.Vl[l]) {
                Vl_mism++;
                frame_bad = 1;
            }
            float dl = fabsf(d_cur.log2Ml[l] - e_cur.log2Ml[l]);
            if (dl > worst) {
                worst = dl;
                worst_frame = f;
                worst_l = l;
            }
            if (!float_bits_equal(d_cur.log2Ml[l], e_cur.log2Ml[l])) {
                log2_mism++;
                frame_bad = 1;
            }
            if (!float_bits_equal(d_cur.Ml[l], e_cur.Ml[l])) {
                Ml_mism++;
                frame_bad = 1;
            }
        }
        if (frame_bad) {
            bad_frames++;
            if (first_bad < 0) {
                first_bad = f;
            }
        }
        mbe_moveMbeParms(&d_cur, &d_prev);
        mbe_moveMbeParms(&e_cur, &e_prev);
    }
    printf("encoder/decoder state parity over %d frames (%d voice):\n", FRAMES, voice);
    printf("  frames with any mismatch: %d (first at %d)\n", bad_frames, first_bad);
    printf("  L mismatches: %d, gamma mismatches: %d, Vl mismatches: %d\n", L_mism, gamma_mism, Vl_mism);
    printf("  bitwise w0 mismatches: %d\n", w0_mism);
    printf("  bitwise log2Ml mismatches: %d, bitwise Ml mismatches: %d\n", log2_mism, Ml_mism);
    printf("  worst |log2Ml| diff: %g (frame %d, l=%d)\n", worst, worst_frame, worst_l);
    return bad_frames ? 1 : 0;
}

static int
test_c1_parity(void) {
    for (unsigned word = 0; word < 4096; word++) {
        char d[49] = {0}, fr[4][24], cw[23];
        for (int i = 0; i < 12; i++) {
            d[i] = (char)((word >> (11 - i)) & 1u);
        }
        for (int i = 12; i < 49; i++) {
            d[i] = (char)(rnd() & 1u);
        }
        if (mbe_encodeAmbe3600x2400Frame(d, fr) != 0) {
            printf("C1 parity: frame encode failed for seed %u\n", word);
            return 1;
        }
        mbe_golay2312_encode(d + 12, cw);
        unsigned parity = 0;
        for (int i = 0; i < 23; i++) {
            parity ^= (unsigned)cw[i];
        }
        uint32_t pr = word * 16u;
        for (int i = 1; i <= 24; i++) {
            pr = (173u * pr + 13849u) & 65535u;
        }
        if (fr[2][10] != (char)(parity ^ (pr >> 15))) {
            printf("C1 parity: scrambled parity mismatch for seed %u\n", word);
            return 1;
        }
    }
    puts("C1 scrambled even parity: ok (4096 seeds)");
    return 0;
}

static int
test_invalid_arguments(mbe_ambe2400_encoder* enc) {
    mbe_ambe2400EncoderReset(enc);
    float pcm[160] = {0};
    short shorts[160] = {0};
    char d[49] = {0}, fr[4][24] = {{0}};
    unsigned char bytes[9] = {0};
    mbe_parms c, p, h;
    mbe_initMbeParms(&c, &p, &h);
    const int invalid = MBE_STATUS_INVALID_ARGUMENT;

    const struct {
        const char* name;
        int status;
    } cases[] = {
        {"float PCM: null context", mbe_encodeAmbe2400Parms(NULL, pcm, d, &c, &p)},
        {"short PCM: null context", mbe_encodeAmbe2400ParmsShort(NULL, shorts, d, &c, &p)},
        {"float PCM: null samples", mbe_encodeAmbe2400Parms(enc, NULL, d, &c, &p)},
        {"float PCM: null bits", mbe_encodeAmbe2400Parms(enc, pcm, NULL, &c, &p)},
        {"float PCM: null current parameters", mbe_encodeAmbe2400Parms(enc, pcm, d, NULL, &p)},
        {"float PCM: null previous parameters", mbe_encodeAmbe2400Parms(enc, pcm, d, &c, NULL)},
        {"short PCM: null samples", mbe_encodeAmbe2400ParmsShort(enc, NULL, d, &c, &p)},
        {"short PCM: null bits", mbe_encodeAmbe2400ParmsShort(enc, shorts, NULL, &c, &p)},
        {"short PCM: null current parameters", mbe_encodeAmbe2400ParmsShort(enc, shorts, d, NULL, &p)},
        {"short PCM: null previous parameters", mbe_encodeAmbe2400ParmsShort(enc, shorts, d, &c, NULL)},
        {"frame encode: null bits", mbe_encodeAmbe3600x2400Frame(NULL, fr)},
        {"frame encode: null frame", mbe_encodeAmbe3600x2400Frame(d, NULL)},
        {"DV encode: null frame", mbe_encodeDStarDVData(NULL, bytes)},
        {"DV encode: null bytes", mbe_encodeDStarDVData((const char (*)[24])fr, NULL)},
        {"DV decode: null bytes", mbe_decodeDStarDVData(NULL, fr)},
        {"DV decode: null frame", mbe_decodeDStarDVData(bytes, NULL)},
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        if (cases[i].status != invalid) {
            printf("invalid arguments: %s returned %d (expected %d)\n", cases[i].name, cases[i].status, invalid);
            return 1;
        }
    }
    for (int i = 0; i < 49; i++) {
        d[i] = 2;
        if (mbe_encodeAmbe3600x2400Frame(d, fr) != MBE_STATUS_INVALID_BITS) {
            printf("invalid arguments: frame encode accepted invalid bit at index %d\n", i);
            return 1;
        }
        d[i] = 0;
    }
    return 0;
}

/* A rejected frame (non-finite or absurdly large samples) leaves the stream
 * exactly as if the frame had never been offered. */
static int
test_invalid_samples(mbe_ambe2400_encoder* enc) {
    const float bad_values[] = {NAN, INFINITY, -INFINITY, 3.40282347e38f /* FLT_MAX */, 2097152.0f};
    for (size_t b = 0; b < sizeof(bad_values) / sizeof(bad_values[0]); b++) {
        char reference[30][49] = {{0}};
        char observed[30][49] = {{0}};
        for (int pass = 0; pass < 2; pass++) {
            mbe_parms cur, prev, enhanced;
            mbe_ambe2400EncoderReset(enc);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            for (int f = 0; f < 30; f++) {
                float pcm[160];
                for (int i = 0; i < 160; i++) {
                    pcm[i] = 0.2f * sinf((float)(2.0 * M_PI * 200.0 * (f * 160 + i) / 8000.0));
                }
                if (pass == 1 && f == 10) {
                    float bad[160];
                    memcpy(bad, pcm, sizeof(bad));
                    bad[37] = bad_values[b];
                    char unused[49];
                    mbe_parms before = cur;
                    if (mbe_encodeAmbe2400Parms(enc, bad, unused, &cur, &prev) != MBE_STATUS_INVALID_ARGUMENT
                        || memcmp(&before, &cur, sizeof(cur)) != 0) {
                        printf("invalid samples: value %zu was not rejected cleanly\n", b);
                        return 1;
                    }
                }
                char* bits = (pass == 0) ? reference[f] : observed[f];
                if (mbe_encodeAmbe2400Parms(enc, pcm, bits, &cur, &prev) < 0) {
                    return 1;
                }
                mbe_moveMbeParms(&cur, &prev);
            }
        }
        if (memcmp(reference, observed, sizeof(reference)) != 0) {
            printf("invalid samples: value %zu changed the stream\n", b);
            return 1;
        }
    }
    puts("non-finite and out-of-range samples are rejected without touching state");
    return 0;
}

/* A prediction history with a non-finite or out-of-range log2Ml or gamma, or
 * one from which the model would overflow, is rejected with cur_mp untouched. */
static int
test_bad_history(mbe_ambe2400_encoder* enc) {
    float pcm[160] = {0};
    char d[49];
    for (int k = 0; k < 5; k++) {
        mbe_parms cur, prev, enhanced;
        mbe_ambe2400EncoderReset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        const mbe_parms before = cur;
        if (k == 0) {
            prev.log2Ml[3] = NAN;
        } else if (k == 1) {
            prev.log2Ml[0] = INFINITY;
        } else if (k == 2) {
            prev.gamma = -3e9f;
        } else if (k == 3) {
            /* Finite, but the prediction overflows: log2 +1000 over the lower
             * half of the harmonics and -1000 over the rest. */
            for (int l = 0; l <= 56; l++) {
                prev.log2Ml[l] = (2 * l <= prev.L) ? 1000.0f : -1000.0f;
            }
        } else {
            prev.gamma = 1000.0f; /* likewise through the gain */
        }
        if (mbe_encodeAmbe2400Parms(enc, pcm, d, &cur, &prev) != MBE_STATUS_INVALID_ARGUMENT) {
            printf("bad history case %d was not rejected\n", k);
            return 1;
        }
        if (!parms_identical(&cur, &before)) {
            printf("bad history case %d changed cur_mp\n", k);
            return 1;
        }
    }
    puts("malformed and overflowing prediction histories are rejected");
    return 0;
}

/* The decoder's own highest states encode. The frame b0..b8 = {118, 0, 63, 511,
 * 123, 13, 0, 0, 0} takes log2Ml[1] past 48 by its 19th repeat and to 48.6 by
 * its 200th; a bound on the history any tighter than the overflow itself
 * would refuse it. */
static int
test_reachable_history(mbe_ambe2400_encoder* enc) {
    static const char frame[] = "1110111111111111111110110010000000000000001111110";
    const int repeats[] = {19, 200};
    for (size_t r = 0; r < sizeof(repeats) / sizeof(repeats[0]); r++) {
        mbe_parms cur, prev, enhanced;
        char d[49];
        for (int i = 0; i < 49; i++) {
            d[i] = (char)(frame[i] - '0');
        }
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int f = 0; f < repeats[r]; f++) {
            if (mbe_decodeAmbe2400Parms(d, &cur, &prev) != 0) {
                return 1;
            }
            mbe_moveMbeParms(&cur, &prev);
        }
        if (prev.log2Ml[1] <= 48.0f) {
            printf("reachable history: log2Ml[1] %.4f after %d repeats\n", (double)prev.log2Ml[1], repeats[r]);
            return 1;
        }
        float pcm[160];
        for (int i = 0; i < 160; i++) {
            pcm[i] = 0.9f * sinf((float)(2.0 * M_PI * 150.0 * i / 8000.0));
        }
        mbe_ambe2400EncoderReset(enc);
        if (mbe_encodeAmbe2400Parms(enc, pcm, d, &cur, &prev) != 0 || !model_finite(&cur)) {
            printf("reachable history after %d repeats was rejected\n", repeats[r]);
            return 1;
        }
    }
    puts("the decoder's highest states (log2Ml past 48) encode");
    return 0;
}

static int
test_prediction_boundaries(mbe_ambe2400_encoder* enc) {
    mbe_ambe2400EncoderReset(enc);
    const int counts[] = {-7, 1, 10, 56, 99};
    for (size_t test = 0; test < sizeof(counts) / sizeof(counts[0]); test++) {
        mbe_parms cur, prev, enhanced, decoded, decode_prev;
        char bits[49];
        short pcm[160];
        mbe_initMbeParms(&cur, &prev, &enhanced);
        prev.L = counts[test];
        for (int l = 0; l <= 56; l++) {
            prev.log2Ml[l] = (float)l * 0.125f;
            prev.Ml[l] = exp2f(prev.log2Ml[l]);
        }
        /* Deliberately poison the zero tap; the decoder replaces it with [1]. */
        prev.log2Ml[0] = -12.0f;
        decode_prev = prev;
        decoded = cur;
        for (int i = 0; i < 160; i++) {
            pcm[i] = (short)(3000.0f * sinf(0.08f * (float)i));
        }
        if (mbe_encodeAmbe2400ParmsShort(enc, pcm, bits, &cur, &prev) != 0
            || mbe_decodeAmbe2400Parms(bits, &decoded, &decode_prev) != 0) {
            return 1;
        }
        if (cur.L != decoded.L || !float_bits_equal(cur.gamma, decoded.gamma)) {
            return 1;
        }
        for (int l = 1; l <= cur.L; l++) {
            if (!float_bits_equal(cur.log2Ml[l], decoded.log2Ml[l])) {
                return 1;
            }
        }
        if (prev.L != counts[test] || prev.log2Ml[0] != -12.0f) {
            return 1;
        }
    }
    puts("prediction boundary taps and harmonic-count clamps: ok");
    return 0;
}

/* D-STAR pitch law calibrated on DVSI's AMBE-3000 vectors (cycles per sample). */
static double
dstar_f0(int b0) {
    return exp2(-4.24738 - 0.0217705 * (b0 + 0.5));
}

/* The decoder reconstructs w0 from b0 with the calibrated law across the voice
 * range, and sizes L like the AMBE+2 table: every harmonic at or below
 * 0.9254 x 4 kHz, the top harmonic DVSI's D-STAR decoder synthesizes. */
static int
test_pitch_law(void) {
    for (int code = 0; code <= 125; code++) {
        mbe_parms cur, prev, enhanced;
        char d[49] = {0};
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int bit = 0; bit < 6; bit++) {
            d[bit] = (char)((code >> (6 - bit)) & 1);
        }
        d[48] = (char)(code & 1);
        if (mbe_decodeAmbe2400Parms(d, &cur, &prev) != 0) {
            continue;
        }
        int want = (int)(0.9254 * 0.5 / dstar_f0(code));
        want = want < 9 ? 9 : (want > 56 ? 56 : want);
        if (cur.L != want) {
            printf("pitch law: b0 %d gives L %d, want %d\n", code, cur.L, want);
            return 1;
        }
    }
    static const int codes[] = {0, 13, 40, 76, 100, 123, 125};
    for (size_t i = 0; i < sizeof(codes) / sizeof(codes[0]); i++) {
        mbe_parms cur, prev, enhanced;
        char d[49] = {0};
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int bit = 0; bit < 6; bit++) {
            d[bit] = (char)((codes[i] >> (6 - bit)) & 1);
        }
        d[48] = (char)(codes[i] & 1);
        if (mbe_decodeAmbe2400Parms(d, &cur, &prev) != 0) {
            printf("pitch law: b0 %d did not decode as voice\n", codes[i]);
            return 1;
        }
        double want = 2.0 * M_PI * dstar_f0(codes[i]);
        if (fabs((double)cur.w0 / want - 1.0) > 1e-5) {
            printf("pitch law: b0 %d gives w0 %.7f, want %.7f\n", codes[i], (double)cur.w0, want);
            return 1;
        }
    }
    puts("pitch law: b0 0..125 follow the calibrated D-STAR law and harmonic count");
    return 0;
}

static int
test_pitch_endpoint(mbe_ambe2400_encoder* enc, int period, int expected_b0) {
    mbe_ambe2400EncoderReset(enc);
    mbe_parms c, p, h;
    char d[49];
    float pcm[160];
    mbe_initMbeParms(&c, &p, &h);
    /* A long steady run, so period 20 exposes the old interior-minimum fallback. */
    for (int frame = 0; frame < 200; frame++) {
        for (int i = 0; i < 160; i++) {
            int phase = (frame * 160 + i) % period;
            pcm[i] = 0.1f * sinf((float)(2.0 * M_PI * phase / period));
        }
        if (mbe_encodeAmbe2400Parms(enc, pcm, d, &c, &p) != 0) {
            return 1;
        }
        mbe_moveMbeParms(&c, &p);
    }
    int b0 = (unsigned char)d[48];
    for (int i = 0; i < 6; i++) {
        b0 |= (int)d[i] << (6 - i);
    }
    printf("pitch period %d: b0=%d (expected %d)\n", period, b0, expected_b0);
    return b0 != expected_b0;
}

static int
test_dc_noise(mbe_ambe2400_encoder* enc) {
    mbe_ambe2400EncoderReset(enc);
    mbe_parms c, p, h;
    char d[49];
    float pcm[160];
    int voiced = 0, total = 0;
    mbe_initMbeParms(&c, &p, &h);
    for (int frame = 0; frame < 40; frame++) {
        for (int i = 0; i < 160; i++) {
            pcm[i] = 0.1f + 0.01f * ((float)(rnd() & 65535u) / 32768.0f - 1.0f);
        }
        if (mbe_encodeAmbe2400Parms(enc, pcm, d, &c, &p) != 0) {
            return 1;
        }
        if (frame >= 10) {
            total += c.L;
            for (int l = 1; l <= c.L; l++) {
                voiced += c.Vl[l];
            }
        }
        mbe_moveMbeParms(&c, &p);
    }
    printf("DC-offset noise: %d/%d voiced harmonics\n", voiced, total);
    return voiced * 20 > total; /* at most 5% */
}

static double
gauss_noise(void) {
    double sum = 0.0;
    for (int i = 0; i < 12; i++) {
        sum += (double)(rnd() & 0xffffu) / 65536.0;
    }
    return sum - 6.0;
}

/* Formant-like envelope: -6 dB/oct tilt with peaks near 500, 1500, 2500 Hz. */
static double
vowel_envelope(double hz) {
    double peaks = 1.0 + (2.0 * exp(-pow((hz - 500.0) / 150.0, 2))) + (1.5 * exp(-pow((hz - 1500.0) / 200.0, 2)))
                   + exp(-pow((hz - 2500.0) / 250.0, 2));
    return peaks / (1.0 + (hz / 300.0));
}

enum fixture_kind {
    FIX_NOISE,
    FIX_COLORED_NOISE,
    FIX_VOWEL,
    FIX_MIXED,
    FIX_MISSING_F0,
    FIX_STRONG_H2,
    FIX_GLIDE,
    FIX_ONSET
};

struct fixture {
    enum fixture_kind kind;
    double f0;    /* Hz; the glide starts here and doubles over the run */
    double level; /* noise RMS or harmonic scale, full scale = 1 */
    double phase[64];
    double lowpass;
};

static double
fixture_f0(const struct fixture* fx, int frame) {
    return (fx->kind == FIX_GLIDE) ? fx->f0 * pow(2.0, (double)frame / 150.0) : fx->f0;
}

static double
fixture_harmonics(struct fixture* fx, double f0) {
    double v = 0.0;
    for (int h = 1; h < 64 && h * f0 < 3800.0; h++) {
        double hz = h * f0;
        double a = vowel_envelope(hz);
        if ((fx->kind == FIX_MIXED && hz >= 2000.0) || (fx->kind == FIX_MISSING_F0 && h == 1)) {
            a = 0.0;
        } else if (fx->kind == FIX_STRONG_H2 && h == 2) {
            a *= 3.0;
        }
        fx->phase[h] += 2.0 * M_PI * hz / 8000.0;
        v += a * sin(fx->phase[h]);
    }
    return v * fx->level;
}

static void
fixture_frame(struct fixture* fx, int frame, float pcm[160]) {
    const double f0 = fixture_f0(fx, frame);
    for (int i = 0; i < 160; i++) {
        double v;
        if (fx->kind == FIX_NOISE || (fx->kind == FIX_ONSET && frame < 40)) {
            v = fx->level * gauss_noise();
        } else if (fx->kind == FIX_COLORED_NOISE) {
            fx->lowpass = (0.8 * fx->lowpass) + (fx->level * gauss_noise()); /* -6 dB/oct above ~300 Hz */
            v = fx->lowpass;
        } else {
            v = fixture_harmonics(fx, f0);
            if (fx->kind == FIX_MIXED) {
                v += 0.05 * fx->level * gauss_noise(); /* about 30 dB below the harmonics */
            }
        }
        pcm[i] = (float)v;
    }
}

struct fixture_stats {
    int frames;
    int voiced_harmonics;
    int harmonics;
    int band_voiced[4];
    int mixed_rows;
    int b0_close;
    int octave_errors;
};

static int
ideal_b0(double hz) {
    return (int)lround(((log2(hz / 8000.0) + 4.24738) / -0.0217705) - 0.5);
}

/* Encode a fixture and count decisions on settled voice frames. */
static int
run_fixture(mbe_ambe2400_encoder* enc, struct fixture* fx, int frames, int settle, struct fixture_stats* st) {
    mbe_parms cur, prev, enhanced;
    memset(st, 0, sizeof(*st));
    mbe_ambe2400EncoderReset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int f = 0; f < frames; f++) {
        float pcm[160];
        char d[49];
        fixture_frame(fx, f, pcm);
        int r = mbe_encodeAmbe2400Parms(enc, pcm, d, &cur, &prev);
        if (r < 0) {
            return 1;
        }
        if (r == 0 && f >= settle) {
            int b0 = d[48];
            for (int i = 0; i < 6; i++) {
                b0 |= d[i] << (6 - i);
            }
            int b1 = (d[38] << 3) | (d[39] << 2) | (d[40] << 1) | d[41];
            int want = ideal_b0(fixture_f0(fx, f));
            st->frames++;
            for (int l = 1; l <= cur.L; l++) {
                st->voiced_harmonics += cur.Vl[l];
                st->harmonics++;
            }
            for (int band = 0; band < 4; band++) {
                st->band_voiced[band] += (b1 >> (3 - band)) & 1;
            }
            st->mixed_rows += (b1 == 0x0c);
            st->b0_close += (b0 >= want - 1 && b0 <= want + 1);
            int distance = (b0 > want) ? b0 - want : want - b0;
            st->octave_errors += (distance > 30); /* 2^(30 * 0.021336) is about 1.56 */
        }
        mbe_moveMbeParms(&cur, &prev);
    }
    return 0;
}

/* Noise of any colour and level is unvoiced (the old encoder voiced all of it). */
static int
test_noise_unvoiced(mbe_ambe2400_encoder* enc) {
    const struct {
        enum fixture_kind kind;
        double level;
    } cases[] = {{FIX_NOISE, 0.1}, {FIX_NOISE, 0.0056}, {FIX_COLORED_NOISE, 0.03}};

    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
        for (int seed = 0; seed < 3; seed++) {
            struct fixture fx = {.kind = cases[c].kind, .level = cases[c].level};
            struct fixture_stats st;
            rng = 0x1234u + (uint32_t)(seed * 7919);
            if (run_fixture(enc, &fx, 200, 50, &st) != 0 || st.harmonics == 0) {
                return 1;
            }
            if (st.voiced_harmonics * 20 > st.harmonics) {
                printf("noise case %zu seed %d: %d/%d voiced harmonics\n", c, seed, st.voiced_harmonics, st.harmonics);
                return 1;
            }
        }
    }
    puts("noise (white -20/-45 dBFS, coloured): at most 5% voiced harmonics");
    return 0;
}

/* Steady vowels across the pitch range are voiced and on pitch. */
static int
test_vowels_voiced(mbe_ambe2400_encoder* enc) {
    const double f0s[] = {70.0, 95.0, 120.0, 150.0, 200.0, 240.0, 310.0};
    for (size_t i = 0; i < sizeof(f0s) / sizeof(f0s[0]); i++) {
        struct fixture fx = {.kind = FIX_VOWEL, .f0 = f0s[i], .level = 0.05};
        struct fixture_stats st;
        if (run_fixture(enc, &fx, 150, 50, &st) != 0 || st.frames == 0) {
            return 1;
        }
        for (int band = 0; band < 4; band++) {
            if (st.band_voiced[band] * 20 < st.frames * 19) {
                printf("vowel %.0f Hz: band %d voiced in %d/%d frames\n", f0s[i], band, st.band_voiced[band],
                       st.frames);
                return 1;
            }
        }
        if (st.b0_close * 20 < st.frames * 19) {
            printf("vowel %.0f Hz: b0 within 1 in %d/%d frames\n", f0s[i], st.b0_close, st.frames);
            return 1;
        }
    }
    puts("vowels 70-310 Hz: voiced in all bands, b0 within 1");
    return 0;
}

/* Harmonics below 2 kHz with noise above are voiced only below 2 kHz. */
static int
test_mixed_bands(mbe_ambe2400_encoder* enc) {
    struct fixture fx = {.kind = FIX_MIXED, .f0 = 200.0, .level = 0.05};
    struct fixture_stats st;
    rng = 0xBEEFu;
    if (run_fixture(enc, &fx, 150, 50, &st) != 0 || st.frames == 0) {
        return 1;
    }
    printf("harmonics below 2 kHz, noise above: b1 = 0x0c in %d/%d frames\n", st.mixed_rows, st.frames);
    return st.mixed_rows * 10 < st.frames * 9;
}

/* Pitch cases that trip simple trackers: no octave errors once settled. */
static int
test_pitch_cases(mbe_ambe2400_encoder* enc) {
    const struct {
        enum fixture_kind kind;
        double f0;
        int settle;
        const char* name;
    } cases[] = {{FIX_VOWEL, 240.0, 3, "240 Hz"},
                 {FIX_MISSING_F0, 120.0, 3, "missing fundamental"},
                 {FIX_STRONG_H2, 150.0, 3, "dominant 2nd harmonic"},
                 {FIX_GLIDE, 100.0, 3, "100-200 Hz glide"},
                 {FIX_ONSET, 150.0, 43, "onset after noise"}};

    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
        struct fixture fx = {.kind = cases[c].kind, .f0 = cases[c].f0, .level = 0.05};
        struct fixture_stats st;
        rng = 0xC0DEu;
        if (run_fixture(enc, &fx, 150, cases[c].settle, &st) != 0 || st.frames == 0) {
            return 1;
        }
        if (st.octave_errors != 0 || st.b0_close * 20 < st.frames * 19) {
            printf("%s: %d octave errors, b0 within 1 in %d/%d frames\n", cases[c].name, st.octave_errors, st.b0_close,
                   st.frames);
            return 1;
        }
    }
    puts("pitch: 240 Hz, missing fundamental, strong 2nd harmonic, glide, onset: on pitch");
    return 0;
}

/* Level of a steady 140 Hz voice, amplitude * sum sin(h x) / h over the given
 * number of harmonics, encoded and then decoded by the process path, in dB re
 * an RMS of 32768 on the int16 scale (float output is int16 / 7). The first 20
 * frames are skipped. */
static double
decoded_level_db(mbe_ambe2400_encoder* enc, double amplitude, int harmonics) {
    mbe_ambe2400EncoderReset(enc);
    mbe_parms ec, ep, eh, dc, dp, dh;
    mbe_initMbeParms(&ec, &ep, &eh);
    mbe_initMbeParms(&dc, &dp, &dh);
    double phase = 0.0, sum = 0.0;
    int n = 0;
    for (int frame = 0; frame < 80; frame++) {
        float pcm[160], out[160];
        char bits[49];
        for (int i = 0; i < 160; i++) {
            phase += 2.0 * M_PI * 140.0 / 8000.0;
            double v = 0.0;
            for (int h = 1; h <= harmonics; h++) {
                v += sin(h * phase) / h;
            }
            pcm[i] = (float)(amplitude * v);
        }
        if (mbe_encodeAmbe2400Parms(enc, pcm, bits, &ec, &ep) != 0
            || mbe_processAmbe2400Dataf(out, NULL, bits, &dc, &dp, &dh) < 0) {
            return 1e9;
        }
        mbe_moveMbeParms(&ec, &ep);
        for (int i = 0; frame >= 20 && i < 160; i++) {
            double x = 7.0 * (double)out[i];
            sum += x * x;
            n++;
        }
    }
    return 10.0 * log10((sum / n) + 1e-30) - (20.0 * log10(32768.0));
}

/* The encoder keeps the input level, as DVSI's does: the decoded level follows
 * the input down 20 and 40 dB, and a -21 dBFS voice decodes near its own level. */
static int
test_level_follows_input(mbe_ambe2400_encoder* enc) {
    const double input_db = 20.0 * log10(0.1 * sqrt(1.5962 / 2.0)); /* RMS of 0.1 sum sin(h x)/h, h 1..20 */
    double loud = decoded_level_db(enc, 0.1, 20);
    double quiet = decoded_level_db(enc, 0.01, 20);
    double faint = decoded_level_db(enc, 0.001, 20);
    printf("level: input %.1f dB decodes at %.1f, -20 dB at %.1f, -40 dB at %.1f\n", input_db, loud, quiet - loud,
           faint - loud);
    return fabs(loud - input_db) > 2.0 || fabs((quiet - loud) + 20.0) > 1.0 || fabs((faint - loud) + 40.0) > 2.0;
}

/* A quiet pure tone has a far lower mean log magnitude than the lowest gain
 * the codec can send. The envelope is flattened just enough that the frame
 * still decodes at the input level instead of putting the excess in the peak. */
static int
test_quiet_tone_level(mbe_ambe2400_encoder* enc) {
    int fails = 0;
    const double amplitudes[] = {0.03, 0.003, 0.0003};
    for (size_t k = 0; k < sizeof(amplitudes) / sizeof(amplitudes[0]); k++) {
        double input_db = 20.0 * log10(amplitudes[k] / sqrt(2.0));
        double decoded = decoded_level_db(enc, amplitudes[k], 1);
        printf("quiet tone: %.1f dBFS decodes at %.1f\n", input_db, decoded);
        fails += fabs(decoded - input_db) > 3.0;
    }
    return fails != 0;
}

/* Quiet and silent input is coded as voice frames, as DVSI's encoder codes it,
 * never as tone or silence frames, and decodes near silence. */
static int
test_quiet_input(mbe_ambe2400_encoder* enc) {
    mbe_ambe2400EncoderReset(enc);
    mbe_parms ec, ep, eh, dc, dp, dh;
    mbe_initMbeParms(&ec, &ep, &eh);
    mbe_initMbeParms(&dc, &dp, &dh);
    double peak = 0.0;
    for (int frame = 0; frame < 60; frame++) {
        float pcm[160], out[160];
        char bits[49];
        for (int i = 0; i < 160; i++) {
            pcm[i] = frame < 30 ? 0.0f : 3e-4f * (float)(((rnd() & 0xffff) / 32768.0) - 1.0); /* -70 dBFS */
        }
        mbe_process_result result;
        mbe_initProcessResult(&result);
        if (mbe_encodeAmbe2400Parms(enc, pcm, bits, &ec, &ep) != 0
            || mbe_processAmbe2400Dataf(out, &result, bits, &dc, &dp, &dh) < 0
            || (result.flags & MBE_PROCESS_FLAG_TONE) != 0) {
            printf("quiet input: frame %d was not coded as voice\n", frame);
            return 1;
        }
        mbe_moveMbeParms(&ec, &ep);
        for (int i = 0; i < 160; i++) {
            peak = fmax(peak, fabs(7.0 * (double)out[i]));
        }
    }
    printf("quiet input: voice frames, decoded peak %.1f\n", peak);
    return peak > 64.0;
}

/* Compare interleaved streams with standalone replays, then replay after reset.
 * The sequence exercises pitch, voicing, PCM history and silent input. */
static int
test_context_stream(mbe_ambe2400_encoder* stream, mbe_ambe2400_encoder* neighbor) {
    mbe_parms cur, prev, enhanced, other_cur, other_prev, other_enhanced;
    char expected[24][49];
    mbe_initMbeParms(&cur, &prev, &enhanced);
    mbe_initMbeParms(&other_cur, &other_prev, &other_enhanced);
    mbe_ambe2400EncoderReset(stream);
    mbe_ambe2400EncoderReset(neighbor);
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            mbe_ambe2400EncoderReset(stream);
            mbe_initMbeParms(&cur, &prev, &enhanced);
        }
        for (int frame = 0; frame < 24; frame++) {
            float pcm[160], noise[160];
            char bits[49], other_bits[49];
            for (int i = 0; i < 160; i++) {
                pcm[i] = frame % 12 < 6 ? 0.01f * sinf(0.1f * (float)(frame * 160 + i)) : 0.0f;
                noise[i] = 0.2f * cosf(0.23f * (float)(frame * 160 + i));
            }
            if (mbe_encodeAmbe2400Parms(stream, pcm, bits, &cur, &prev) < 0) {
                return 1;
            }
            if (pass == 0) {
                memcpy(expected[frame], bits, sizeof(bits));
                if (mbe_encodeAmbe2400Parms(neighbor, noise, other_bits, &other_cur, &other_prev) < 0) {
                    return 1;
                }
                mbe_moveMbeParms(&other_cur, &other_prev);
            } else if (memcmp(expected[frame], bits, sizeof(bits)) != 0) {
                return 1;
            }
            mbe_moveMbeParms(&cur, &prev);
        }
    }
    return 0;
}

static int
test_contexts(mbe_ambe2400_encoder* enc) {
    mbe_ambe2400_encoder* other = mbe_ambe2400EncoderAlloc();
    if (other == NULL) {
        return 1;
    }
    int fails = test_context_stream(enc, other);
    fails += test_context_stream(other, enc);
    mbe_ambe2400EncoderFree(other);
    mbe_ambe2400EncoderFree(NULL);
    mbe_ambe2400EncoderReset(NULL);
    printf("independent contexts and reset replay: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

int
main(void) {
    mbe_ambe2400_encoder* enc = mbe_ambe2400EncoderAlloc();
    if (enc == NULL) {
        return 1;
    }
    int fails = 0;
    fails += test_invalid_arguments(enc);
    fails += test_invalid_samples(enc);
    fails += test_bad_history(enc);
    fails += test_reachable_history(enc);
    fails += test_c1_parity();
    fails += test_golay();
    fails += test_frame_roundtrip();
    fails += test_frame_single_bit_errors();
    fails += test_dv_bytes();
    fails += test_state_parity(enc);
    fails += test_prediction_boundaries(enc);
    fails += test_pitch_law();
    /* 400 Hz, the shortest analysed period, is b0 3; b0 0..2 reach 402..418 Hz. */
    fails += test_pitch_endpoint(enc, 20, 3);
    fails += test_pitch_endpoint(enc, 127, 125);
    fails += test_dc_noise(enc);
    fails += test_noise_unvoiced(enc);
    fails += test_vowels_voiced(enc);
    fails += test_mixed_bands(enc);
    fails += test_pitch_cases(enc);
    fails += test_contexts(enc);
    fails += test_level_follows_input(enc);
    fails += test_quiet_tone_level(enc);
    fails += test_quiet_input(enc);
    mbe_ambe2400EncoderFree(enc);
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}

#endif
