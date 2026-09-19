// SPDX-License-Identifier: GPL-2.0-or-later
/* Round-trip harness for PR #84: frame FEC, DV byte packing, Golay, and
 * encoder->decoder parameter-state parity. */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "mbelib-neo/mbelib.h"

#ifdef MBE_ENCODER_TEST_OOM
/* GNU link wrapping injects failure at each project-owned aligned allocation.
 * PFFFT's own internal calls are defined in its object, so are not wrapped. */
void* encoder_real_alloc(size_t size) __asm__("__real_pffft_aligned_malloc");
void encoder_real_free(void* ptr) __asm__("__real_pffft_aligned_free");
extern void* encoder_fault_alloc(size_t size) __asm__("__wrap_pffft_aligned_malloc");
extern void encoder_track_free(void* ptr) __asm__("__wrap_pffft_aligned_free");
static int allocation_count;
static int fail_at;
static int live_allocations;

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

int
main(int argc, char** argv) {
    if (argc != 2 || strlen(argv[1]) != 1 || argv[1][0] < '1' || argv[1][0] > '5') {
        return 1;
    }
    fail_at = argv[1][0] - '0';
    mbe_parms cur, prev, enhanced;
    float pcm[160];
    char bits[49];
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int i = 0; i < 160; i++) {
        pcm[i] = 0.1f * sinf(0.1f * (float)i);
    }
    int status = mbe_encodeAmbe2400Parms(pcm, bits, &cur, &prev);
    if (status >= 0 || live_allocations != 0) {
        (void)fprintf(stderr, "failed allocation %d: status=%d live=%d\n", fail_at, status, live_allocations);
        return 1;
    }
    fail_at = 0;
    if (mbe_encodeAmbe2400Parms(pcm, bits, &cur, &prev) != 0) {
        return 1;
    }
    (void)puts("allocation failure: propagated, cleaned up, retry succeeded");
    return 0;
}
#else

#include <stdint.h>

#include "mbe_ecc.h"

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
            }
        }
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
test_state_parity(void) {
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
    int L_mism = 0, Vl_mism = 0, gamma_mism = 0;
    for (int f = 0; f < FRAMES; f++) {
        char d[49];
        unsigned char saved_prev[sizeof(e_prev)];
        memcpy(saved_prev, &e_prev, sizeof(e_prev));
        int r = mbe_encodeAmbe2400ParmsShort(pcm + (size_t)f * 160, d, &e_cur, &e_prev);
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
                if (dl > 1e-3f) {
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
            return 1;
        }
    }
    puts("C1 scrambled even parity: ok (4096 seeds)");
    return 0;
}

static int
test_invalid_arguments(void) {
    float pcm[160] = {0};
    short shorts[160] = {0};
    char d[49] = {0}, fr[4][24] = {{0}};
    unsigned char bytes[9] = {0};
    mbe_parms c, p, h;
    mbe_initMbeParms(&c, &p, &h);
    const int invalid = MBE_STATUS_INVALID_ARGUMENT;
    if (mbe_encodeAmbe2400Parms(NULL, d, &c, &p) != invalid || mbe_encodeAmbe2400Parms(pcm, NULL, &c, &p) != invalid
        || mbe_encodeAmbe2400Parms(pcm, d, NULL, &p) != invalid
        || mbe_encodeAmbe2400Parms(pcm, d, &c, NULL) != invalid) {
        return 1;
    }
    if (mbe_encodeAmbe2400ParmsShort(NULL, d, &c, &p) != invalid
        || mbe_encodeAmbe2400ParmsShort(shorts, NULL, &c, &p) != invalid
        || mbe_encodeAmbe2400ParmsShort(shorts, d, NULL, &p) != invalid
        || mbe_encodeAmbe2400ParmsShort(shorts, d, &c, NULL) != invalid) {
        return 1;
    }
    if (mbe_encodeAmbe3600x2400Frame(NULL, fr) != invalid || mbe_encodeAmbe3600x2400Frame(d, NULL) != invalid
        || mbe_encodeDStarDVData(NULL, bytes) != invalid
        || mbe_encodeDStarDVData((const char (*)[24])fr, NULL) != invalid || mbe_decodeDStarDVData(NULL, fr) != invalid
        || mbe_decodeDStarDVData(bytes, NULL) != invalid) {
        return 1;
    }
    for (int i = 0; i < 49; i++) {
        d[i] = 2;
        if (mbe_encodeAmbe3600x2400Frame(d, fr) != MBE_STATUS_INVALID_BITS) {
            return 1;
        }
        d[i] = 0;
    }
    return 0;
}

static int
test_prediction_boundaries(void) {
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
        if (mbe_encodeAmbe2400ParmsShort(pcm, bits, &cur, &prev) != 0
            || mbe_decodeAmbe2400Parms(bits, &decoded, &decode_prev) != 0) {
            return 1;
        }
        if (cur.L != decoded.L || fabsf(cur.gamma - decoded.gamma) > 1e-4f) {
            return 1;
        }
        for (int l = 1; l <= cur.L; l++) {
            if (fabsf(cur.log2Ml[l] - decoded.log2Ml[l]) > 1e-3f) {
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
test_vuv_hysteresis(void) {
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
            if (mbe_encodeAmbe2400Parms(pcm, bits, &cur, &prev) != 0) {
                return 1;
            }
        }
        for (int l = 1; l <= 56; l++) {
            prev.Vl[l] = 0;
        }
        if (mbe_encodeAmbe2400Parms(pcm, bits, &cur, &prev) != 0) {
            return 1;
        }
        int without_history = cur.K;
        int harmonics = cur.L;
        for (int l = 1; l <= 56; l++) {
            prev.Vl[l] = 1;
        }
        if (mbe_encodeAmbe2400Parms(pcm, bits, &cur, &prev) != 0) {
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

/* Each analysis case runs in its own process because analysis state is TLS. */
static int
test_pitch_endpoint(int period, int expected_b0) {
    mbe_parms c, p, h;
    char d[49];
    float pcm[160];
    mbe_initMbeParms(&c, &p, &h);
    for (int frame = 0; frame < 16; frame++) {
        for (int i = 0; i < 160; i++) {
            int phase = (frame * 160 + i) % period;
            pcm[i] = 0.1f * sinf((float)(2.0 * M_PI * phase / period));
        }
        if (mbe_encodeAmbe2400Parms(pcm, d, &c, &p) != 0) {
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
test_dc_noise(void) {
    mbe_parms c, p, h;
    char d[49];
    float pcm[160];
    int voiced = 0, total = 0;
    mbe_initMbeParms(&c, &p, &h);
    for (int frame = 0; frame < 40; frame++) {
        for (int i = 0; i < 160; i++) {
            pcm[i] = 0.1f + 0.01f * ((float)(rnd() & 65535u) / 32768.0f - 1.0f);
        }
        if (mbe_encodeAmbe2400Parms(pcm, d, &c, &p) != 0) {
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

int
main(int argc, char** argv) {
    if (argc == 2 && strcmp(argv[1], "pitch20") == 0) {
        return test_pitch_endpoint(20, 0);
    }
    if (argc == 2 && strcmp(argv[1], "pitch127") == 0) {
        return test_pitch_endpoint(127, 125);
    }
    if (argc == 2 && strcmp(argv[1], "hysteresis") == 0) {
        return test_vuv_hysteresis();
    }
    if (argc == 2 && strcmp(argv[1], "dc_noise") == 0) {
        return test_dc_noise();
    }
    int fails = 0;
    fails += test_invalid_arguments();
    fails += test_c1_parity();
    fails += test_golay();
    fails += test_frame_roundtrip();
    fails += test_frame_single_bit_errors();
    fails += test_dv_bytes();
    fails += test_state_parity();
    fails += test_prediction_boundaries();
    printf("%s\n", fails ? "SOME TESTS FAILED" : "ALL OK");
    return fails ? 1 : 0;
}

#endif
