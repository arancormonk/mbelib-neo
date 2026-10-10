// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/* Codec-neutral checks for the AMBE+2 and IMBE encoder tests: each test
 * describes its codec with a struct enc_codec and runs these against it. */
#ifndef MBELIB_TESTS_ENCODER_TEST_SUPPORT_H
#define MBELIB_TESTS_ENCODER_TEST_SUPPORT_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct enc_codec {
    const char* name;
    int bits;      /* parameter bits per frame */
    int tone_kind; /* encoder return value for a tone frame, or -1 */
    void* (*alloc)(void);
    void (*reset)(void* enc);
    void (*release)(void* enc);
    int (*encode)(void* enc, const float* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev);
    int (*encode_short)(void* enc, const short* pcm, char* bits, mbe_parms* cur, const mbe_parms* prev);
    int (*decode)(const char* bits, mbe_parms* cur, mbe_parms* prev);
    int (*process)(float* out, mbe_process_result* result, const char* bits, mbe_parms* cur, mbe_parms* prev,
                   mbe_parms* enhanced);
    const char* peak_frame;           /* parameter bits ('0'/'1') whose repetition takes the decoder's log2Ml highest */
    int (*tone_id)(const char* bits); /* the tone index a tone frame carries, else -1; NULL without tone frames */
};

static uint32_t enc_rng = 0x2450u;

static inline uint32_t
enc_rnd(void) {
    enc_rng = (enc_rng * 1664525u) + 1013904223u;
    return enc_rng >> 8;
}

static inline double
enc_gauss(void) {
    double sum = 0.0;
    for (int i = 0; i < 12; i++) {
        sum += (double)(enc_rnd() & 0xffffu) / 65536.0;
    }
    return sum - 6.0;
}

static inline int
enc_float_bits_equal(float a, float b) {
    uint32_t x, y;
    memcpy(&x, &a, sizeof(x));
    memcpy(&y, &b, sizeof(y));
    return x == y;
}

/* Bitwise identity of two parameter sets (compared as bytes, as exact
 * reproducibility is what the checks want). */
static inline int
enc_parms_identical(const mbe_parms* a, const mbe_parms* b) {
    unsigned char x[sizeof(mbe_parms)];
    unsigned char y[sizeof(mbe_parms)];
    memcpy(x, a, sizeof(x));
    memcpy(y, b, sizeof(y));
    return memcmp(x, y, sizeof(x)) == 0;
}

/* At most 2^127 in magnitude: a range test on the bit pattern read from
 * memory, which unlike an exponent-mask test or a by-value float survives
 * fast-math. */
static inline int
enc_float_finite(const float* x) {
    uint32_t bits;
    memcpy(&bits, x, sizeof(bits));
    return (bits & 0x7FFFFFFFu) <= 0x7F000000u;
}

/* w0, gamma, and log2Ml and Ml over the model's harmonics are finite. */
static inline int
enc_model_finite(const mbe_parms* mp) {
    int ok = enc_float_finite(&mp->w0) && enc_float_finite(&mp->gamma);
    for (int l = 1; ok && l <= mp->L && l <= 56; l++) {
        ok = enc_float_finite(&mp->log2Ml[l]) && enc_float_finite(&mp->Ml[l]);
    }
    return ok;
}

/* Decoded fundamental in Hz. */
static inline double
enc_f0_hz(const mbe_parms* mp) {
    return (double)mp->w0 * 8000.0 / (2.0 * M_PI);
}

/* A speech-like test signal: a gliding harmonic voice, noise, silence and a
 * steady 200 Hz voice in 200 ms segments. */
static inline void
enc_speech_like(short* pcm, int n) {
    double phase = 0.0;
    for (int i = 0; i < n; i++) {
        int segment = (i / 1600) % 6;
        double t = (double)i / 8000.0;
        double v = 0.0;
        if (segment == 0 || segment == 1 || segment == 3) {
            phase += 2.0 * M_PI * (120.0 + (60.0 * sin(2.0 * M_PI * 0.7 * t))) / 8000.0;
            for (int h = 1; h <= 20; h++) {
                v += sin(h * phase) / ((h * 0.8) + 0.5);
            }
            v *= 0.25 * (0.6 + (0.4 * sin(2.0 * M_PI * 3.0 * t)));
        } else if (segment == 2) {
            v = 0.05 * (((double)(enc_rnd() & 0xffff) / 32768.0) - 1.0);
        } else if (segment == 5) {
            phase += 2.0 * M_PI * 200.0 / 8000.0;
            for (int h = 1; h <= 10; h++) {
                v += sin(h * phase) / h;
            }
            v *= 0.15;
        }
        v = (v > 0.99) ? 0.99 : ((v < -0.99) ? -0.99 : v);
        pcm[i] = (short)(v * 32767.0);
    }
}

/* Formant-like envelope: -6 dB/oct tilt with peaks near 500, 1500, 2500 Hz. */
static inline double
enc_vowel_envelope(double hz) {
    double peaks = 1.0 + (2.0 * exp(-pow((hz - 500.0) / 150.0, 2))) + (1.5 * exp(-pow((hz - 1500.0) / 200.0, 2)))
                   + exp(-pow((hz - 2500.0) / 250.0, 2));
    return peaks / (1.0 + (hz / 300.0));
}

enum enc_fixture_kind { ENC_NOISE, ENC_COLORED_NOISE, ENC_VOWEL, ENC_MISSING_F0, ENC_STRONG_H2, ENC_GLIDE, ENC_SERIES };

struct enc_fixture {
    enum enc_fixture_kind kind;
    double f0;    /* Hz; a glide starts here and doubles over 150 frames */
    double level; /* noise RMS or harmonic scale, full scale 1 */
    double phase[64];
    double lowpass;
};

static inline double
enc_fixture_f0(const struct enc_fixture* fx, int frame) {
    return (fx->kind == ENC_GLIDE) ? fx->f0 * pow(2.0, (double)frame / 150.0) : fx->f0;
}

static inline double
enc_fixture_harmonics(struct enc_fixture* fx, double f0) {
    double v = 0.0;
    for (int h = 1; h < 64 && h * f0 < 3800.0; h++) {
        double hz = h * f0;
        double a = enc_vowel_envelope(hz);
        if (fx->kind == ENC_SERIES) {
            a = (h <= 8) ? 1.0 / h : 0.0; /* eight harmonics, amplitude 1/h */
        } else if (fx->kind == ENC_MISSING_F0 && h == 1) {
            a = 0.0;
        } else if (fx->kind == ENC_STRONG_H2 && h == 2) {
            a *= 3.0;
        }
        fx->phase[h] += 2.0 * M_PI * hz / 8000.0;
        v += a * sin(fx->phase[h]);
    }
    return v * fx->level;
}

static inline void
enc_fixture_frame(struct enc_fixture* fx, int frame, float pcm[160]) {
    const double f0 = enc_fixture_f0(fx, frame);
    for (int i = 0; i < 160; i++) {
        double v;
        if (fx->kind == ENC_NOISE) {
            v = fx->level * enc_gauss();
        } else if (fx->kind == ENC_COLORED_NOISE) {
            fx->lowpass = (0.8 * fx->lowpass) + (fx->level * enc_gauss());
            v = fx->lowpass;
        } else {
            v = enc_fixture_harmonics(fx, f0);
        }
        pcm[i] = (float)v;
    }
}

struct enc_fixture_stats {
    int frames;           /* settled voice frames */
    int other_frames;     /* settled frames that are not voice frames */
    int on_pitch;         /* decoded f0 within 3% of the input's */
    int octave_errors;    /* decoded f0 off by more than 30% */
    int voiced_harmonics; /* below 3 kHz */
    int harmonics;        /* below 3 kHz */
};

/* Encode a fixture through the codec and judge the settled frames by the
 * model the decoder plays. */
static inline int
enc_run_fixture(const struct enc_codec* c, void* enc, struct enc_fixture* fx, int frames, int settle,
                struct enc_fixture_stats* st) {
    mbe_parms cur, prev, enhanced;
    memset(st, 0, sizeof(*st));
    c->reset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int f = 0; f < frames; f++) {
        float pcm[160];
        char bits[88];
        enc_fixture_frame(fx, f, pcm);
        int r = c->encode(enc, pcm, bits, &cur, &prev);
        if (r < 0) {
            return 1;
        }
        st->other_frames += (r != 0 && f >= settle);
        if (r == 0 && f >= settle) {
            const double want = enc_fixture_f0(fx, f);
            st->frames++;
            if (want > 0.0) { /* noise has no pitch to judge */
                double ratio = enc_f0_hz(&cur) / want;
                st->on_pitch += fabs(ratio - 1.0) <= 0.03;
                st->octave_errors += ratio < 0.7 || ratio > 1.4;
            }
            for (int l = 1; l <= cur.L && (double)l * enc_f0_hz(&cur) < 3000.0; l++) {
                st->voiced_harmonics += cur.Vl[l];
                st->harmonics++;
            }
        }
        mbe_moveMbeParms(&cur, &prev);
    }
    return 0;
}

/* The encoder's cur_mp is the decoder's model of the emitted bits; prev_mp is
 * never touched; tone frames leave cur_mp a copy of prev_mp. */
static inline int
enc_test_state_parity(const struct enc_codec* c, void* enc) {
    enum { FRAMES = 600 };

    static short pcm[FRAMES * 160];
    mbe_parms ec, ep, eh, dc, dp, dh;
    int mismatches = 0;
    enc_speech_like(pcm, FRAMES * 160);
    c->reset(enc);
    mbe_initMbeParms(&ec, &ep, &eh);
    mbe_initMbeParms(&dc, &dp, &dh);
    for (int f = 0; f < FRAMES; f++) {
        char bits[88];
        mbe_parms before = ep;
        int r = c->encode_short(enc, pcm + ((size_t)f * 160), bits, &ec, &ep);
        if (r < 0 || !enc_parms_identical(&before, &ep)) {
            printf("%s state parity: frame %d returned %d or changed prev_mp\n", c->name, f, r);
            return 1;
        }
        if (r == c->tone_kind) {
            mismatches += !enc_parms_identical(&ec, &ep);
            continue;
        }
        if (c->decode(bits, &dc, &dp) != 0) {
            return 1;
        }
        int bad = dc.L != ec.L || !enc_float_bits_equal(dc.w0, ec.w0) || !enc_float_bits_equal(dc.gamma, ec.gamma);
        for (int l = 1; l <= dc.L && l <= 56; l++) {
            bad |= dc.Vl[l] != ec.Vl[l] || !enc_float_bits_equal(dc.log2Ml[l], ec.log2Ml[l])
                   || !enc_float_bits_equal(dc.Ml[l], ec.Ml[l]);
        }
        mismatches += bad;
        mbe_moveMbeParms(&dc, &dp);
        mbe_moveMbeParms(&ec, &ep);
    }
    printf("%s encoder/decoder state parity over %d frames: %d mismatching frames\n", c->name, FRAMES, mismatches);
    return mismatches != 0;
}

/* A rejected frame leaves the stream exactly as if it had never been offered. */
static inline int
enc_test_invalid_samples(const struct enc_codec* c, void* enc) {
    const float bad_values[] = {NAN, INFINITY, -INFINITY, 3.40282347e38f, 2097152.0f};
    for (size_t b = 0; b < sizeof(bad_values) / sizeof(bad_values[0]); b++) {
        char reference[30][88];
        char observed[30][88];
        for (int pass = 0; pass < 2; pass++) {
            mbe_parms cur, prev, enhanced;
            c->reset(enc);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            for (int f = 0; f < 30; f++) {
                float pcm[160];
                for (int i = 0; i < 160; i++) {
                    pcm[i] = 0.2f * sinf((float)(2.0 * M_PI * 200.0 * (f * 160 + i) / 8000.0));
                    pcm[i] += 0.05f * sinf((float)(2.0 * M_PI * 400.0 * (f * 160 + i) / 8000.0));
                }
                if (pass == 1 && f == 10) {
                    float bad[160];
                    char unused[88];
                    mbe_parms before = cur;
                    memcpy(bad, pcm, sizeof(bad));
                    bad[37] = bad_values[b];
                    if (c->encode(enc, bad, unused, &cur, &prev) != MBE_STATUS_INVALID_ARGUMENT
                        || !enc_parms_identical(&before, &cur)) {
                        printf("%s invalid samples: value %zu was not rejected cleanly\n", c->name, b);
                        return 1;
                    }
                }
                if (c->encode(enc, pcm, (pass == 0) ? reference[f] : observed[f], &cur, &prev) < 0) {
                    return 1;
                }
                mbe_moveMbeParms(&cur, &prev);
            }
        }
        for (int f = 0; f < 30; f++) {
            if (memcmp(reference[f], observed[f], (size_t)c->bits) != 0) {
                printf("%s invalid samples: value %zu changed the stream\n", c->name, b);
                return 1;
            }
        }
    }
    printf("%s: non-finite and out-of-range samples are rejected without touching state\n", c->name);
    return 0;
}

/* Prediction against any previous harmonic count and edge values matches the
 * decoder's, without modifying the caller's prev_mp. */
static inline int
enc_test_prediction_boundaries(const struct enc_codec* c, void* enc) {
    const int counts[] = {-7, 1, 10, 56, 99};
    for (size_t test = 0; test < sizeof(counts) / sizeof(counts[0]); test++) {
        mbe_parms cur, prev, enhanced, decoded, decode_prev;
        char bits[88];
        short pcm[160];
        c->reset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        prev.L = counts[test];
        for (int l = 0; l <= 56; l++) {
            prev.log2Ml[l] = (float)l * 0.125f;
            prev.Ml[l] = exp2f(prev.log2Ml[l]);
        }
        prev.log2Ml[0] = -12.0f; /* the decoders read index 1 there */
        decode_prev = prev;
        decoded = cur;
        for (int i = 0; i < 160; i++) {
            pcm[i] = (short)(3000.0f * sinf(0.08f * (float)i) + 900.0f * sinf(0.21f * (float)i));
        }
        if (c->encode_short(enc, pcm, bits, &cur, &prev) != 0 || c->decode(bits, &decoded, &decode_prev) != 0) {
            return 1;
        }
        if (cur.L != decoded.L || !enc_float_bits_equal(cur.gamma, decoded.gamma) || prev.L != counts[test]
            || !enc_float_bits_equal(prev.log2Ml[0], -12.0f)) {
            return 1;
        }
        for (int l = 1; l <= cur.L; l++) {
            if (!enc_float_bits_equal(cur.log2Ml[l], decoded.log2Ml[l])) {
                return 1;
            }
        }
    }
    printf("%s: prediction boundary taps and harmonic-count clamps match the decoder\n", c->name);
    return 0;
}

/* The same frames as float and as 16-bit PCM give the same bits. */
static inline int
enc_test_short_float(const struct enc_codec* c, void* enc) {
    enum { FRAMES = 60 };

    static short pcm[FRAMES * 160];
    char from_short[FRAMES][88];
    mbe_parms cur, prev, enhanced;
    enc_speech_like(pcm, FRAMES * 160);
    for (int pass = 0; pass < 2; pass++) {
        c->reset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int f = 0; f < FRAMES; f++) {
            float x[160];
            char bits[88];
            const short* frame = pcm + ((size_t)f * 160);
            for (int i = 0; i < 160; i++) {
                x[i] = (float)frame[i] / 32768.0f;
            }
            int r = (pass == 0) ? c->encode_short(enc, frame, from_short[f], &cur, &prev)
                                : c->encode(enc, x, bits, &cur, &prev);
            if (r < 0 || (pass == 1 && memcmp(bits, from_short[f], (size_t)c->bits) != 0)) {
                printf("%s short/float parity: frame %d differs\n", c->name, f);
                return 1;
            }
            mbe_moveMbeParms(&cur, &prev);
        }
    }
    printf("%s: 16-bit and float input encode alike\n", c->name);
    return 0;
}

/* Interleaved streams match standalone runs; a reset replays the stream. */
static inline int
enc_test_context_stream(const struct enc_codec* c, void* stream, void* neighbor) {
    mbe_parms cur, prev, enhanced, other_cur, other_prev, other_enhanced;
    char expected[24][88];
    mbe_initMbeParms(&cur, &prev, &enhanced);
    mbe_initMbeParms(&other_cur, &other_prev, &other_enhanced);
    c->reset(stream);
    c->reset(neighbor);
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            c->reset(stream);
            mbe_initMbeParms(&cur, &prev, &enhanced);
        }
        for (int frame = 0; frame < 24; frame++) {
            float pcm[160], noise[160];
            char bits[88], other_bits[88];
            for (int i = 0; i < 160; i++) {
                pcm[i] = frame % 12 < 6 ? 0.01f * sinf(0.1f * (float)(frame * 160 + i)) : 0.0f;
                pcm[i] += frame % 12 < 6 ? 0.004f * sinf(0.3f * (float)(frame * 160 + i)) : 0.0f;
                noise[i] = 0.2f * cosf(0.23f * (float)(frame * 160 + i)) + 0.05f * (float)enc_gauss();
            }
            if (c->encode(stream, pcm, bits, &cur, &prev) < 0) {
                return 1;
            }
            if (pass == 0) {
                memcpy(expected[frame], bits, (size_t)c->bits);
                if (c->encode(neighbor, noise, other_bits, &other_cur, &other_prev) < 0) {
                    return 1;
                }
                mbe_moveMbeParms(&other_cur, &other_prev);
            } else if (memcmp(expected[frame], bits, (size_t)c->bits) != 0) {
                return 1;
            }
            mbe_moveMbeParms(&cur, &prev);
        }
    }
    return 0;
}

static inline int
enc_test_contexts(const struct enc_codec* c, void* enc) {
    void* other = c->alloc();
    if (other == NULL) {
        return 1;
    }
    int fails = enc_test_context_stream(c, enc, other);
    fails += enc_test_context_stream(c, other, enc);
    c->release(other);
    c->release(NULL);
    c->reset(NULL);
    printf("%s independent contexts and reset replay: %s\n", c->name, fails ? "FAIL" : "ok");
    return fails;
}

/* Noise of any colour and level stays essentially unvoiced. The bound is 10%
 * rather than the D-STAR encoder's 5%: these codecs decide voicing per
 * three-harmonic band (TIA-102.BABA 5.2), and low-pass noise leaves a few
 * low bands voiced; an encoder that voices noise voices most of it. */
static inline int
enc_test_noise_unvoiced(const struct enc_codec* c, void* enc) {
    const struct {
        enum enc_fixture_kind kind;
        double level;
    } cases[] = {{ENC_NOISE, 0.1}, {ENC_NOISE, 0.0056}, {ENC_COLORED_NOISE, 0.03}};

    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        struct enc_fixture fx = {.kind = cases[k].kind, .level = cases[k].level};
        struct enc_fixture_stats st;
        enc_rng = 0x1234u + (uint32_t)k;
        if (enc_run_fixture(c, enc, &fx, 200, 50, &st) != 0 || st.harmonics == 0 || st.other_frames != 0) {
            return 1;
        }
        if (st.voiced_harmonics * 10 > st.harmonics) {
            printf("%s noise case %zu: %d/%d voiced harmonics\n", c->name, k, st.voiced_harmonics, st.harmonics);
            return 1;
        }
    }
    printf("%s: white and coloured noise at most 10%% voiced\n", c->name);
    return 0;
}

/* Steady vowels across the pitch range are voiced and on pitch. */
static inline int
enc_test_vowels(const struct enc_codec* c, void* enc) {
    const double f0s[] = {70.0, 95.0, 120.0, 150.0, 200.0, 240.0, 310.0};
    for (size_t i = 0; i < sizeof(f0s) / sizeof(f0s[0]); i++) {
        struct enc_fixture fx = {.kind = ENC_VOWEL, .f0 = f0s[i], .level = 0.05};
        struct enc_fixture_stats st;
        if (enc_run_fixture(c, enc, &fx, 150, 50, &st) != 0 || st.frames != 100 || st.other_frames != 0) {
            return 1;
        }
        if (st.on_pitch * 20 < st.frames * 19 || st.voiced_harmonics * 10 < st.harmonics * 9) {
            printf("%s vowel %.0f Hz: on pitch %d/%d, voiced %d/%d harmonics\n", c->name, f0s[i], st.on_pitch,
                   st.frames, st.voiced_harmonics, st.harmonics);
            return 1;
        }
    }
    printf("%s: vowels 70-310 Hz voiced and on pitch\n", c->name);
    return 0;
}

/* Signals that trip simple pitch trackers: no octave errors once settled. */
static inline int
enc_test_pitch_cases(const struct enc_codec* c, void* enc) {
    const struct {
        enum enc_fixture_kind kind;
        double f0;
        const char* name;
    } cases[] = {{ENC_VOWEL, 240.0, "240 Hz vowel"},         {ENC_MISSING_F0, 120.0, "missing fundamental"},
                 {ENC_STRONG_H2, 150.0, "dominant 2nd"},     {ENC_GLIDE, 100.0, "100-200 Hz glide"},
                 {ENC_SERIES, 400.0 / 3.0, "133 Hz series"}, {ENC_SERIES, 200.0, "200 Hz series"},
                 {ENC_SERIES, 390.0, "390 Hz series"}};

    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        struct enc_fixture fx = {.kind = cases[k].kind, .f0 = cases[k].f0, .level = 0.05};
        struct enc_fixture_stats st;
        if (enc_run_fixture(c, enc, &fx, 150, 5, &st) != 0 || st.frames != 145 || st.other_frames != 0) {
            return 1;
        }
        if (st.octave_errors != 0 || st.on_pitch * 20 < st.frames * 19) {
            printf("%s %s: %d octave errors, on pitch in %d/%d frames\n", c->name, cases[k].name, st.octave_errors,
                   st.on_pitch, st.frames);
            return 1;
        }
    }
    printf("%s: 240 Hz, missing fundamental, strong 2nd, glide, 133/200/390 Hz series on pitch\n", c->name);
    return 0;
}

/* Decoded level (process path, first 20 frames skipped) of a steady 140 Hz
 * voice amplitude * sum sin(h x) / h over the given number of harmonics, in dB
 * re an RMS of 32768 on the int16 scale (float output is int16 / 7). */
static inline double
enc_decoded_level_db(const struct enc_codec* c, void* enc, double amplitude, int harmonics) {
    mbe_parms ec, ep, eh, dc, dp, dh;
    double phase = 0.0, sum = 0.0;
    int n = 0;
    c->reset(enc);
    mbe_initMbeParms(&ec, &ep, &eh);
    mbe_initMbeParms(&dc, &dp, &dh);
    for (int frame = 0; frame < 80; frame++) {
        float pcm[160], out[160];
        char bits[88];
        for (int i = 0; i < 160; i++) {
            double v = 0.0;
            phase += 2.0 * M_PI * 140.0 / 8000.0;
            for (int h = 1; h <= harmonics; h++) {
                v += sin(h * phase) / h;
            }
            pcm[i] = (float)(amplitude * v);
        }
        if (c->encode(enc, pcm, bits, &ec, &ep) != 0 || c->process(out, NULL, bits, &dc, &dp, &dh) < 0) {
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

/* The decoded level follows the input down 20 and 40 dB, and a -21 dBFS voice
 * decodes near its own level. */
static inline int
enc_test_level(const struct enc_codec* c, void* enc) {
    const double input_db = 20.0 * log10(0.1 * sqrt(1.5962 / 2.0)); /* RMS of 0.1 sum sin(h x)/h, h 1..20 */
    double loud = enc_decoded_level_db(c, enc, 0.1, 20);
    double quiet = enc_decoded_level_db(c, enc, 0.01, 20);
    double faint = enc_decoded_level_db(c, enc, 0.001, 20);
    printf("%s level: input %.1f dB decodes at %.1f, -20 dB at %.1f, -40 dB at %.1f\n", c->name, input_db, loud,
           quiet - loud, faint - loud);
    return fabs(loud - input_db) > 2.0 || fabs((quiet - loud) + 20.0) > 1.0 || fabs((faint - loud) + 40.0) > 2.0;
}

/* A quiet pure tone has a far lower mean log magnitude than the lowest gain
 * the codec can send; it still decodes at the input level instead of putting
 * the excess in its peak. */
static inline int
enc_test_quiet_tone_level(const struct enc_codec* c, void* enc) {
    int fails = 0;
    const double amplitudes[] = {0.03, 0.003, 0.0003};
    for (size_t k = 0; k < sizeof(amplitudes) / sizeof(amplitudes[0]); k++) {
        double input_db = 20.0 * log10(amplitudes[k] / sqrt(2.0));
        double decoded = enc_decoded_level_db(c, enc, amplitudes[k], 1);
        printf("%s quiet tone: %.1f dBFS decodes at %.1f\n", c->name, input_db, decoded);
        fails += fabs(decoded - input_db) > 3.0;
    }
    return fails != 0;
}

/* A prediction history with a non-finite or out-of-range value, or one from
 * which the model would overflow (log2 +1000 over the lower half of the
 * harmonics and -1000 over the rest, which no mean removal cancels), is rejected
 * without touching state: cur_mp keeps its contents and the stream continues
 * as if the frame had never been offered. */
static inline int
enc_test_bad_history(const struct enc_codec* c, void* enc) {
    static const struct {
        int harmonic; /* the log2Ml index, -1 for gamma, or -2 for the +-value step */
        float value;
    } cases[] = {{3, NAN}, {0, INFINITY}, {1, 3e9f}, {-1, NAN}, {-1, -INFINITY}, {-1, -3e9f}, {-2, 1000.0f}};

    for (size_t b = 0; b < sizeof(cases) / sizeof(cases[0]); b++) {
        char reference[12][88];
        char observed[12][88];
        for (int pass = 0; pass < 2; pass++) {
            mbe_parms cur, prev, enhanced;
            c->reset(enc);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            for (int f = 0; f < 12; f++) {
                float pcm[160];
                for (int i = 0; i < 160; i++) {
                    /* A voice, not a tone: a tone frame never reads the history. */
                    float phase = (float)(2.0 * M_PI * 150.0 * (f * 160 + i) / 8000.0);
                    pcm[i] = (0.1f * sinf(phase)) + (0.05f * sinf(2.0f * phase)) + (0.03f * sinf(3.0f * phase));
                }
                if (pass == 1 && f == 6) {
                    mbe_parms bad = prev;
                    const mbe_parms before = cur;
                    char unused[88];
                    if (cases[b].harmonic == -2) {
                        for (int l = 0; l <= 56; l++) {
                            bad.log2Ml[l] = (2 * l <= bad.L) ? cases[b].value : -cases[b].value;
                        }
                    } else if (cases[b].harmonic < 0) {
                        bad.gamma = cases[b].value;
                    } else {
                        bad.log2Ml[cases[b].harmonic] = cases[b].value;
                    }
                    if (c->encode(enc, pcm, unused, &cur, &bad) != MBE_STATUS_INVALID_ARGUMENT) {
                        printf("%s bad history %zu was not rejected\n", c->name, b);
                        return 1;
                    }
                    if (!enc_parms_identical(&cur, &before)) {
                        printf("%s bad history %zu changed cur_mp\n", c->name, b);
                        return 1;
                    }
                }
                if (c->encode(enc, pcm, (pass == 0) ? reference[f] : observed[f], &cur, &prev) < 0) {
                    return 1;
                }
                mbe_moveMbeParms(&cur, &prev);
            }
        }
        for (int f = 0; f < 12; f++) {
            if (memcmp(reference[f], observed[f], (size_t)c->bits) != 0) {
                printf("%s bad history %zu changed the stream\n", c->name, b);
                return 1;
            }
        }
    }
    printf("%s: malformed and overflowing prediction histories are rejected without touching state\n", c->name);
    return 0;
}

/* One harmonic at +48 over the rest at -48 with gamma +-48, at the edge of
 * what any decoder reaches, encodes to a finite model; histories far beyond
 * it (+-1000, +-2^20) encode to a finite model or are rejected with cur_mp
 * untouched. */
static inline int
enc_test_extreme_history(const struct enc_codec* c, void* enc) {
    const float magnitudes[] = {48.0f, 1000.0f, 1048576.0f};
    int rejected = 0;
    for (size_t m = 0; m < sizeof(magnitudes) / sizeof(magnitudes[0]); m++) {
        for (int k = 0; k < 4; k++) {
            mbe_parms cur, prev, enhanced;
            float pcm[160];
            char bits[88];
            c->reset(enc);
            mbe_initMbeParms(&cur, &prev, &enhanced);
            const mbe_parms before = cur;
            prev.L = (k & 1) ? 9 : 56;
            for (int l = 0; l <= 56; l++) {
                prev.log2Ml[l] = (l == 1) ? magnitudes[m] : -magnitudes[m];
            }
            prev.gamma = (k & 2) ? -magnitudes[m] : magnitudes[m];
            for (int i = 0; i < 160; i++) {
                pcm[i] = (k & 1) ? 0.9f * sinf(0.3f * (float)i) : 0.0f;
            }
            int status = c->encode(enc, pcm, bits, &cur, &prev);
            if (status == MBE_STATUS_INVALID_ARGUMENT && m > 0 && enc_parms_identical(&cur, &before)) {
                rejected++;
                continue;
            }
            if (status < 0 || !enc_model_finite(&cur)) {
                printf("%s extreme history %g/%d: status %d, finite %d\n", c->name, (double)magnitudes[m], k, status,
                       enc_model_finite(&cur));
                return 1;
            }
        }
    }
    printf("%s: extreme histories give a finite model (%d of 8 far beyond the edge rejected)\n", c->name, rejected);
    return 0;
}

/* The decoder's own highest state (its peak frame repeated) encodes: no state a
 * decoder reaches is mistaken for a malformed one. */
static inline int
enc_test_reachable_history(const struct enc_codec* c, void* enc) {
    mbe_parms cur, prev, enhanced;
    char bits[88];
    double peak = -1e9;
    if (strlen(c->peak_frame) != (size_t)c->bits) {
        return 1;
    }
    for (int i = 0; i < c->bits; i++) {
        bits[i] = (char)(c->peak_frame[i] - '0');
    }
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int f = 0; f < 40; f++) {
        if (c->decode(bits, &cur, &prev) != 0) {
            printf("%s peak frame did not decode as voice\n", c->name);
            return 1;
        }
        mbe_moveMbeParms(&cur, &prev);
    }
    for (int l = 1; l <= prev.L; l++) {
        peak = fmax(peak, (double)prev.log2Ml[l]);
    }
    c->reset(enc);
    for (int f = 0; f < 2; f++) {
        float pcm[160];
        for (int i = 0; i < 160; i++) {
            pcm[i] = (f == 0) ? 0.0f : 0.9f * sinf((float)(2.0 * M_PI * 150.0 * i / 8000.0));
        }
        int status = c->encode(enc, pcm, bits, &cur, &prev);
        if (status < 0 || !enc_model_finite(&cur)) {
            printf("%s: the decoder's highest state was rejected (status %d)\n", c->name, status);
            return 1;
        }
    }
    printf("%s: the decoder's highest state (log2 amplitude %.1f) encodes\n", c->name, peak);
    return peak < 30.0;
}

/* Silent and -70 dBFS input is coded as voice frames that decode near silence. */
static inline int
enc_test_quiet_input(const struct enc_codec* c, void* enc) {
    mbe_parms ec, ep, eh, dc, dp, dh;
    double peak = 0.0;
    c->reset(enc);
    mbe_initMbeParms(&ec, &ep, &eh);
    mbe_initMbeParms(&dc, &dp, &dh);
    for (int frame = 0; frame < 60; frame++) {
        float pcm[160], out[160];
        char bits[88];
        mbe_process_result result;
        for (int i = 0; i < 160; i++) {
            pcm[i] = frame < 30 ? 0.0f : 3e-4f * (float)(((enc_rnd() & 0xffff) / 32768.0) - 1.0);
        }
        mbe_initProcessResult(&result);
        if (c->encode(enc, pcm, bits, &ec, &ep) != 0 || c->process(out, &result, bits, &dc, &dp, &dh) < 0
            || (result.flags & (MBE_PROCESS_FLAG_TONE | MBE_PROCESS_FLAG_SILENCE)) != 0) {
            printf("%s quiet input: frame %d was not coded as voice\n", c->name, frame);
            return 1;
        }
        mbe_moveMbeParms(&ec, &ep);
        for (int i = 0; i < 160; i++) {
            peak = fmax(peak, fabs(7.0 * (double)out[i]));
        }
    }
    printf("%s quiet input: voice frames, decoded peak %.1f\n", c->name, peak);
    return peak > 64.0;
}

/* One frame of zeros flushes the tail: a burst in the last 40 samples of the
 * final input frame reaches the model by the flush frame. */
static inline int
enc_test_flush(const struct enc_codec* c, void* enc) {
    mbe_parms cur, prev, enhanced;
    double energy[12] = {0};
    c->reset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int frame = 0; frame < 12; frame++) {
        float pcm[160] = {0};
        char bits[88];
        for (int i = 120; frame == 10 && i < 160; i++) {
            pcm[i] = 0.3f * sinf((float)(2.0 * M_PI * 1000.0 * i / 8000.0)) + 0.2f * sinf((float)(0.9 * i));
        }
        if (c->encode(enc, pcm, bits, &cur, &prev) < 0) {
            return 1;
        }
        for (int l = 1; l <= cur.L; l++) {
            energy[frame] += (double)cur.Ml[l] * (double)cur.Ml[l];
        }
        mbe_moveMbeParms(&cur, &prev);
    }
    printf("%s flush: model energy before the burst %.3g, after one zero frame %.3g\n", c->name, energy[9], energy[11]);
    return energy[11] < 100.0 * (energy[9] + 1e-6);
}

/* A frame of a dial tone (350 + 440 Hz) of the given peak per component over
 * stream samples [on, off), silent elsewhere; then from samples at and after
 * dtmf_on, DTMF 5 (770 + 1336 Hz) at the same level. */
static inline void
enc_dial_frame(float pcm[160], int frame, int on, int off, int dtmf_on, double peak) {
    for (int i = 0; i < 160; i++) {
        int n = (frame * 160) + i;
        double t = (double)n / 8000.0;
        double v = 0.0;
        if (n >= on && n < off) {
            v = sin(2.0 * M_PI * 350.0 * t) + sin(2.0 * M_PI * 440.0 * t);
        } else if (dtmf_on >= 0 && n >= dtmf_on) {
            v = sin(2.0 * M_PI * 770.0 * t) + sin(2.0 * M_PI * 1336.0 * t);
        }
        pcm[i] = (float)(peak * v);
    }
}

/* Encode `frames` frames of enc_dial_frame() input; ids[f] is the tone index
 * frame f carries, or -1. If bad_at >= 0, a louder copy of that frame is first
 * offered with a finite history whose model overflows (log2 +1000 over the
 * lower half of the harmonics, -1000 above), which voice and tone frames alike
 * reject after their analysis. */
static inline int
enc_dial_run(const struct enc_codec* c, void* enc, int on, int off, int dtmf_on, int frames, int bad_at, int ids[],
             char bits[][88]) {
    mbe_parms cur, prev, enhanced;
    c->reset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int f = 0; f < frames; f++) {
        float pcm[160];
        if (f == bad_at) {
            mbe_parms bad = prev;
            const mbe_parms before = cur;
            enc_dial_frame(pcm, f, on, off, dtmf_on, 0.11);
            for (int l = 0; l <= 56; l++) {
                bad.log2Ml[l] = (2 * l <= bad.L) ? 1000.0f : -1000.0f;
            }
            if (c->encode(enc, pcm, bits[f], &cur, &bad) != MBE_STATUS_INVALID_ARGUMENT
                || !enc_parms_identical(&cur, &before)) {
                printf("%s call progress: the bad history at frame %d was not rejected cleanly\n", c->name, f);
                return 1;
            }
        }
        enc_dial_frame(pcm, f, on, off, dtmf_on, 0.1);
        if (c->encode(enc, pcm, bits[f], &cur, &prev) < 0) {
            return 1;
        }
        ids[f] = c->tone_id(bits[f]);
        mbe_moveMbeParms(&cur, &prev);
    }
    return 0;
}

/*
 * Call-progress timing (mbe_tone_track()): a dial tone (Table 9 index 160,
 * D-STAR 144, passed as dial) is sent once it has filled three spans and held
 * two frames after. 45 ms bursts never are, at any alignment; an 80 ms tone is
 * sent for frames 7..9; a frame rejected during confirmation leaves the
 * stream unchanged; a 20 ms dropout is bridged; DTMF 5 (index 133 in both
 * codecs) that follows is sent from its first detection.
 */
static inline int
enc_test_call_progress(const struct enc_codec* c, void* enc, int dial) {
    enum { FRAMES = 40 };

    int ids[FRAMES];
    int again[FRAMES];
    static char bits[FRAMES][88];
    static char replay[FRAMES][88];
    for (int shift = 0; shift < 160; shift += 20) {
        if (enc_dial_run(c, enc, 800 + shift, 1160 + shift, -1, 14, -1, ids, bits) != 0) {
            return 1;
        }
        for (int f = 0; f < 14; f++) {
            if (ids[f] >= 0) {
                printf("%s call progress: a 45 ms burst (shift %d) sent tone %d at frame %d\n", c->name, shift, ids[f],
                       f);
                return 1;
            }
        }
    }
    if (enc_dial_run(c, enc, 640, 1280, -1, 14, -1, ids, bits) != 0
        || enc_dial_run(c, enc, 640, 1280, -1, 14, 6, again, replay) != 0) {
        return 1;
    }
    for (int f = 0; f < 14; f++) {
        if (ids[f] != ((f >= 7 && f <= 9) ? dial : -1) || memcmp(bits[f], replay[f], (size_t)c->bits) != 0) {
            printf("%s call progress: 80 ms tone frame %d carries %d, replay %s\n", c->name, f, ids[f],
                   memcmp(bits[f], replay[f], (size_t)c->bits) ? "differs" : "matches");
            return 1;
        }
    }
    {
        /* A dial tone over frames 0..29 with frame 15 (2400..2559) silent. */
        int dropout[34];
        static char dropout_bits[34][88];
        mbe_parms cur, prev, enhanced;
        c->reset(enc);
        mbe_initMbeParms(&cur, &prev, &enhanced);
        for (int f = 0; f < 34; f++) {
            float pcm[160];
            enc_dial_frame(pcm, f, 0, 4800, -1, 0.1);
            for (int i = 0; i < 160; i++) {
                int n = (f * 160) + i;
                pcm[i] = (n >= 2400 && n < 2560) ? 0.0f : pcm[i];
            }
            if (c->encode(enc, pcm, dropout_bits[f], &cur, &prev) < 0) {
                return 1;
            }
            dropout[f] = c->tone_id(dropout_bits[f]);
            mbe_moveMbeParms(&cur, &prev);
            if (dropout[f] != ((f >= 3 && f <= 31) ? dial : -1)) {
                printf("%s call progress: with a 20 ms dropout, frame %d carries %d\n", c->name, f, dropout[f]);
                return 1;
            }
        }
    }
    if (enc_dial_run(c, enc, 0, 3200, 3200, FRAMES, -1, ids, bits) != 0) {
        return 1;
    }
    int first_dtmf = -1;
    for (int f = 0; f < FRAMES; f++) {
        if (first_dtmf < 0 && ids[f] == 133) {
            first_dtmf = f;
        }
        if (first_dtmf >= 0 && ids[f] != 133) {
            printf("%s call progress: frame %d after DTMF carries %d\n", c->name, f, ids[f]);
            return 1;
        }
    }
    if (first_dtmf < 20 || first_dtmf > 21) {
        printf("%s call progress: DTMF after the dial tone first sent at frame %d\n", c->name, first_dtmf);
        return 1;
    }
    printf("%s call progress: sent after three spans and held two frames; 45 ms bursts are voice\n", c->name);
    return 0;
}

/* Encode `frames` frames of tone a over stream samples [on, change), then tone
 * b over [change, off), silent elsewhere; a tone is {low Hz or 0 for a single
 * tone, high Hz, peak per component}. ids[f] is the tone index frame f
 * carries, or -1. */
static inline int
enc_tone_run(const struct enc_codec* c, void* enc, const double a[3], const double b[3], int on, int change, int off,
             int frames, int ids[]) {
    mbe_parms cur, prev, enhanced;
    c->reset(enc);
    mbe_initMbeParms(&cur, &prev, &enhanced);
    for (int f = 0; f < frames; f++) {
        float pcm[160];
        char bits[88];
        for (int i = 0; i < 160; i++) {
            int n = (f * 160) + i;
            const double* tone = (n < change) ? a : b;
            double t = (double)n / 8000.0;
            double v = ((tone[0] > 0.0) ? sin(2.0 * M_PI * tone[0] * t) : 0.0) + sin(2.0 * M_PI * tone[1] * t);
            pcm[i] = (n >= on && n < off) ? (float)(tone[2] * v) : 0.0f;
        }
        if (c->encode(enc, pcm, bits, &cur, &prev) < 0) {
            return 1;
        }
        ids[f] = c->tone_id(bits);
        mbe_moveMbeParms(&cur, &prev);
    }
    return 0;
}

/* Whether ids are voice, then DTMF 5 (index 133), then `next`, then voice,
 * with none of them skipped and nothing else. */
static inline int
enc_tone_sequence_ok(const int ids[], int frames, int next) {
    int seen = 0;
    for (int f = 0; f < frames; f++) {
        int phase = (ids[f] == 133) ? 1 : ((ids[f] == next) ? 2 : ((ids[f] < 0 && seen > 0) ? 3 : 0));
        if ((ids[f] >= 0 && phase == 0) || (phase != seen && phase != seen + 1)) {
            return 0;
        }
        seen = phase;
    }
    return seen == 3;
}

/*
 * Tone edges (mbe_tone_track()): a 60 ms DTMF 5 (index 133) or 1 kHz tone
 * (index 32) is sent for three or four frames, those whose span it fills at
 * least 45%, wherever it starts within the detector's 20-sample blocks, and is
 * not held after it ends. DTMF 5 changing directly to DTMF 9 (index digit9)
 * leaves no voice frame between the two wherever the change falls, at a level
 * well above the detector's floor and at one just above it (a 34 peak per
 * component on the 16-bit scale, below the floor for half a span).
 */
static inline int
enc_test_tone_edges(const struct enc_codec* c, void* enc, int digit9) {
    enum { FRAMES = 24 };

    static const double tones[2][3] = {{770.0, 1336.0, 0.1}, {0.0, 1000.0, 0.1}};
    static const int tone_ids[2] = {133, 32};
    int ids[FRAMES];
    for (int t = 0; t < 2; t++) {
        for (int start = 640; start < 800; start += 5) {
            if (enc_tone_run(c, enc, tones[t], tones[t], start, start + 480, start + 480, FRAMES, ids) != 0) {
                return 1;
            }
            int sent = 0;
            int first = -1;
            int last = -1;
            for (int f = 0; f < FRAMES; f++) {
                sent += (ids[f] >= 0);
                first = (first < 0 && ids[f] >= 0) ? f : first;
                last = (ids[f] >= 0) ? f : last;
            }
            int ok = sent >= 3 && sent <= 4 && last - first + 1 == sent;
            for (int f = first; ok && f <= last; f++) {
                ok = ids[f] == tone_ids[t];
            }
            if (!ok) {
                printf("%s tone edges: a 60 ms tone %d from sample %d sent as %d frames %d..%d\n", c->name, tone_ids[t],
                       start, sent, first, last);
                return 1;
            }
        }
    }
    for (int quiet = 0; quiet < 2; quiet++) {
        const double peak = quiet ? 34.0 / 32768.0 : 0.1;
        const double five[3] = {770.0, 1336.0, peak};
        const double nine[3] = {852.0, 1477.0, peak};
        for (int change = 1600; change < 1760; change += 10) {
            if (enc_tone_run(c, enc, five, nine, 640, change, change + 960, FRAMES, ids) != 0) {
                return 1;
            }
            if (!enc_tone_sequence_ok(ids, FRAMES, digit9)) {
                printf("%s tone edges: DTMF 5 changing to DTMF 9 at sample %d (peak %.4g): frame tones", c->name,
                       change, peak);
                for (int f = 0; f < FRAMES; f++) {
                    printf(" %d", ids[f]);
                }
                printf("\n");
                return 1;
            }
        }
    }
    printf("%s tone edges: 60 ms tones sent for 3-4 frames at any alignment; no voice frame at a change\n", c->name);
    return 0;
}

static inline int
enc_run_common(const struct enc_codec* c, void* enc) {
    int fails = 0;
    fails += enc_test_state_parity(c, enc);
    fails += enc_test_invalid_samples(c, enc);
    fails += enc_test_prediction_boundaries(c, enc);
    fails += enc_test_short_float(c, enc);
    fails += enc_test_contexts(c, enc);
    fails += enc_test_noise_unvoiced(c, enc);
    fails += enc_test_vowels(c, enc);
    fails += enc_test_pitch_cases(c, enc);
    fails += enc_test_level(c, enc);
    fails += enc_test_quiet_tone_level(c, enc);
    fails += enc_test_bad_history(c, enc);
    fails += enc_test_extreme_history(c, enc);
    fails += enc_test_reachable_history(c, enc);
    fails += enc_test_quiet_input(c, enc);
    fails += enc_test_flush(c, enc);
    return fails;
}

#endif /* MBELIB_TESTS_ENCODER_TEST_SUPPORT_H */
