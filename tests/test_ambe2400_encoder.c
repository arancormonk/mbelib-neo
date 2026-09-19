// SPDX-License-Identifier: GPL-2.0-or-later
/* Round-trip harness for PR #84: frame FEC, DV byte packing, Golay, and
 * encoder->decoder parameter-state parity. */
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <string.h>
#include "mbelib-neo/mbelib.h"

#ifdef MBE_ENCODER_TEST_OOM
/* GNU link wrapping covers the context and each project-owned aligned
 * allocation. PFFFT's same-object internal aligned calls are not wrapped. */
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
static void* live_context;

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
    live_context = encoder_real_calloc(count, size);
    return live_context;
}

extern void
encoder_track_heap_free(void* ptr) {
    if (ptr == live_context) {
        live_context = NULL;
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

static int
check_no_encode_allocations(mbe_ambe2400_encoder* enc) {
    mbe_parms cur, prev, enhanced;
    float pcm[160];
    char bits[49];
    mbe_initMbeParms(&cur, &prev, &enhanced);
    int before = allocation_count;
    int heap_before = heap_allocation_count;
    fail_at = before + 1;
    for (int frame = 0; frame < 12; frame++) {
        for (int i = 0; i < 160; i++) {
            pcm[i] = frame < 6 ? 0.01f * sinf(0.1f * (float)i) : 0.0f;
        }
        if (mbe_encodeAmbe2400Parms(enc, pcm, bits, &cur, &prev) < 0 || allocation_count != before
            || heap_allocation_count != heap_before) {
            return 1;
        }
        mbe_moveMbeParms(&cur, &prev);
    }
    mbe_ambe2400EncoderReset(enc);
    short shorts[160] = {0};
    return mbe_encodeAmbe2400ParmsShort(enc, shorts, bits, &cur, &prev) < 0 || allocation_count != before
           || heap_allocation_count != heap_before;
}

int
main(int argc, char** argv) {
    if (argc != 2 || strlen(argv[1]) != 1 || argv[1][0] < '1' || argv[1][0] > '6') {
        return 1;
    }
    fail_at = argv[1][0] - '0';
    mbe_ambe2400_encoder* enc = mbe_ambe2400EncoderAlloc();
    if (enc != NULL || live_allocations != 0 || live_context != NULL) {
        (void)fprintf(stderr, "failed allocation %d: live=%d context=%p\n", fail_at, live_allocations, live_context);
        mbe_ambe2400EncoderFree(enc);
        return 1;
    }
    fail_at = 0;
    enc = mbe_ambe2400EncoderAlloc();
    if (enc == NULL) {
        return 1;
    }
    int failed = check_no_encode_allocations(enc);
    mbe_ambe2400EncoderFree(enc);
    if (failed || live_allocations != 0 || live_context != NULL) {
        return 1;
    }
    (void)puts("allocation failure: NULL, no leaks; retry, reset and encoding: no allocations");
    return 0;
}
#else

#include <stdint.h>

#include "mbe_ecc.h"

static int
float_bits_equal(float a, float b) {
    uint32_t a_bits, b_bits;
    memcpy(&a_bits, &a, sizeof(a_bits));
    memcpy(&b_bits, &b, sizeof(b_bits));
    return a_bits == b_bits;
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

    int bad_frames = 0, silence = 0, voice = 0, first_bad = -1;
    float worst = 0;
    int worst_frame = -1, worst_l = -1;
    int L_mism = 0, Vl_mism = 0, gamma_mism = 0, w0_mism = 0, log2_mism = 0;
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
        if (r < 0) {
            printf("  encode error %d at frame %d\n", r, f);
            return 1;
        }
        int dr = mbe_decodeAmbe2400Parms(d, &d_cur, &d_prev);
        if (r == 1) {
            silence++;
            if (dr != 3) {
                printf("  frame %d: encoder said silence but decoder returned %d\n", f, dr);
                bad_frames++;
            }
            /* process path re-inits decoder state on silence */
            mbe_initMbeParms(&d_cur, &d_prev, &d_enh);
        } else {
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
            if (fabsf(d_cur.gamma - e_cur.gamma) > 1e-4f) {
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
                float scale = fmaxf(1.0f, fabsf(d_cur.Ml[l]));
                if (fabsf(d_cur.Ml[l] - e_cur.Ml[l]) > 1e-4f * scale) {
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
        }
        mbe_moveMbeParms(&e_cur, &e_prev);
    }
    printf("encoder/decoder state parity over %d frames (%d voice, %d silence):\n", FRAMES, voice, silence);
    printf("  frames with any mismatch: %d (first at %d)\n", bad_frames, first_bad);
    printf("  L mismatches: %d, gamma mismatches: %d, Vl mismatches: %d\n", L_mism, gamma_mism, Vl_mism);
    printf("  bitwise w0 mismatches: %d\n", w0_mism);
    printf("  bitwise log2Ml mismatches: %d\n", log2_mism);
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
        if (cur.L != decoded.L || fabsf(cur.gamma - decoded.gamma) > 1e-4f) {
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

static int
test_vuv_hysteresis(mbe_ambe2400_encoder* enc) {
    mbe_ambe2400EncoderReset(enc);
    int retained = 0;
    for (int fixture = 0; fixture < 64; fixture++) {
        mbe_parms cur, prev, enhanced;
        float pcm[160];
        char bits[49];
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int i = 0; i < 160; i++) {
            pcm[i] = 0.0001f * ((float)(rnd() & 65535u) / 32768.0f - 1.0f)
                     + 0.03f * sinf((float)(2.0 * M_PI * (fixture % 8 + 1) * i / 160));
            if (i == 0) {
                pcm[i] += 0.01f * (float)(fixture + 1);
            }
        }
        /* Repetition makes the analysis window, pitch and AGC settle. */
        for (int frame = 0; frame < 200; frame++) {
            if (mbe_encodeAmbe2400Parms(enc, pcm, bits, &cur, &prev) != 0) {
                return 1;
            }
        }
        for (int l = 1; l <= 56; l++) {
            prev.Vl[l] = 0;
        }
        if (mbe_encodeAmbe2400Parms(enc, pcm, bits, &cur, &prev) != 0) {
            return 1;
        }
        int without_history = cur.K;
        int harmonics = cur.L;
        for (int l = 1; l <= 56; l++) {
            prev.Vl[l] = 1;
        }
        if (mbe_encodeAmbe2400Parms(enc, pcm, bits, &cur, &prev) != 0) {
            return 1;
        }
        if (cur.L != harmonics || cur.K < without_history) {
            return 1;
        }
        if (cur.K > without_history) {
            retained++;
        }
    }
    printf("voiced history retains additional bands in %d fixtures\n", retained);
    return retained == 0;
}

static int
test_pitch_endpoint(mbe_ambe2400_encoder* enc, int period, int expected_b0) {
    mbe_ambe2400EncoderReset(enc);
    mbe_parms c, p, h;
    char d[49];
    float pcm[160];
    mbe_initMbeParms(&c, &p, &h);
    /* Let the AGC converge so period 20 exposes the old interior-minimum fallback. */
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
    return voiced > total / 2;
}

/* Compare interleaved streams with standalone replays, then replay after reset.
 * The sequence exercises pitch, voicing, AGC, PCM history and the silence gate. */
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
    fails += test_c1_parity();
    fails += test_golay();
    fails += test_frame_roundtrip();
    fails += test_frame_single_bit_errors();
    fails += test_dv_bytes();
    fails += test_state_parity(enc);
    fails += test_prediction_boundaries(enc);
    fails += test_pitch_endpoint(enc, 20, 0);
    fails += test_pitch_endpoint(enc, 127, 125);
    fails += test_dc_noise(enc);
    fails += test_vuv_hysteresis(enc);
    fails += test_contexts(enc);
    mbe_ambe2400EncoderFree(enc);
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}

#endif
