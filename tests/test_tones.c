// SPDX-License-Identifier: GPL-2.0-or-later
/**
 * @file
 * @brief Tone frames against levels and frequencies measured on DVSI's decoders.
 *
 * Each tone component's RMS level, in dB relative to an RMS of 32768, follows
 * the law fitted to DVSI's AMBE-3000 decoded tone vectors:
 *  - AMBE+2: TIA-102.BABA-1 7.2, 0.711 dB per AD step from +3.17 dBm0 at AD 127,
 *    where +3.17 dBm0 decodes at -2.52 dB per component;
 *  - D-STAR: 3.54 + 0.355 (volume - 255) dB from the 8-bit tone volume.
 * Dual tones carry both components at that level. NOTONES builds output silence.
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mbelib-neo/mbelib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FRAMES 12

static void
put_bits(char* d, int start, int value, int count) {
    for (int i = 0; i < count; ++i) {
        d[start + i] = (char)((value >> (count - 1 - i)) & 1);
    }
}

/* TIA-102.BABA-1 Table 10 tone frame with redundant ID copies and amplitude AD. */
static void
pack_ambe2450_tone(char d[49], int tone_id, int ad) {
    memset(d, 0, 49);
    put_bits(d, 0, 63, 6);
    put_bits(d, 6, ad >> 1, 6);
    put_bits(d, 12, tone_id, 8);
    put_bits(d, 20, tone_id >> 4, 4);
    put_bits(d, 24, tone_id, 4);
    put_bits(d, 28, tone_id >> 1, 7);
    put_bits(d, 35, tone_id, 1);
    put_bits(d, 36, tone_id, 8);
    put_bits(d, 44, ad, 1);
}

/* D-STAR tone frame: b0 126, tone index in bits 6..11, 42, 43 (bits 6..8 carry
 * the index's top three bits through a permutation) and the 8-bit volume in
 * bits 12..16, 44, 45, 17. */
static void
pack_dstar_tone(char d[49], int tone_id, int volume) {
    static const int selector[8] = {1, 2, 3, 4, 0, 7, 6, 5}; /* index bits 7..5 -> bits 6..8 */
    memset(d, 0, 49);
    put_bits(d, 0, 63, 6);
    put_bits(d, 6, selector[tone_id >> 5], 3);
    d[9] = (char)((tone_id >> 4) & 1);
    d[42] = (char)((tone_id >> 3) & 1);
    d[43] = (char)((tone_id >> 2) & 1);
    d[10] = (char)((tone_id >> 1) & 1);
    d[11] = (char)(tone_id & 1);
    put_bits(d, 12, volume >> 3, 5);
    d[44] = (char)((volume >> 2) & 1);
    d[45] = (char)((volume >> 1) & 1);
    d[17] = (char)(volume & 1);
}

/* Amplitude of the component at `hz`, in dB re an RMS of 32768 on the int16 scale. */
static double
component_db(const float* pcm, int n, double hz) {
    double re = 0.0, im = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = 7.0 * (double)pcm[i]; /* float output is int16 / 7 */
        re += x * cos(2.0 * M_PI * hz * i / 8000.0);
        im += x * sin(2.0 * M_PI * hz * i / 8000.0);
    }
    double amplitude = 2.0 * sqrt(re * re + im * im) / n;
    return 20.0 * log10(amplitude / sqrt(2.0) / 32768.0 + 1e-30);
}

static double
energy(const float* pcm, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += (double)pcm[i] * pcm[i];
    }
    return sum;
}

typedef int (*process_fn)(float*, mbe_process_result*, const char[49], mbe_parms*, mbe_parms*, mbe_parms*);

/* Play the same frame FRAMES times from a fresh stream; skip the first frame's onset. */
static void
play(process_fn process, const char d[49], float* out, unsigned* flags) {
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    float frame[160];
    *flags = 0;
    for (int f = 0; f < FRAMES; ++f) {
        mbe_process_result result;
        mbe_initProcessResult(&result);
        assert(process(frame, &result, d, &cur, &prev, &enhanced) >= 0);
        *flags |= result.flags;
        if (f > 0) {
            memcpy(out + ((size_t)(f - 1) * 160u), frame, sizeof(frame));
        }
    }
}

static void
expect_tone(const char* name, process_fn process, const char d[49], double f1, double f2, double level_db) {
    float pcm[(FRAMES - 1) * 160];
    unsigned flags;
    const int n = (FRAMES - 1) * 160;
    play(process, d, pcm, &flags);
    assert(flags & MBE_PROCESS_FLAG_TONE);
#ifdef MBELIB_TEST_NOTONES
    (void)f1;
    (void)f2;
    (void)level_db;
    assert(energy(pcm, n) == 0.0);
    printf("%s: silent (NOTONES)\n", name);
#else
    double a = component_db(pcm, n, f1);
    double b = f2 > 0.0 ? component_db(pcm, n, f2) : level_db;
    printf("%s: %.1f Hz %.2f dB, %.1f Hz %.2f dB (want %.2f)\n", name, f1, a, f2, b, level_db);
    assert(fabs(a - level_db) < 0.1);
    assert(fabs(b - level_db) < 0.1);
    /* Everything else is at least 40 dB down: the components are the whole signal. */
    double components = pow(10.0, a / 10.0) + (f2 > 0.0 ? pow(10.0, b / 10.0) : 0.0);
    double total = 10.0 * log10(49.0 * energy(pcm, n) / n / (32768.0 * 32768.0));
    assert(fabs(10.0 * log10(components) - total) < 0.05);
#endif
}

static double
ambe2450_level(int ad) {
    return -2.523 + (90.30 / 127.0) * (ad - 127);
}

static double
dstar_level(int volume) {
    return 3.540 + 0.35495 * (volume - 255);
}

static void
test_ambe2450_tones(void) {
    char d[49];
    pack_ambe2450_tone(d, 30, 111); /* 937.5 Hz */
    expect_tone("AMBE+2 single AD 111", mbe_processAmbe2450Dataf, d, 937.5, 0.0, ambe2450_level(111));
    pack_ambe2450_tone(d, 129, 99); /* DTMF 1 */
    expect_tone("AMBE+2 DTMF 1 AD 99", mbe_processAmbe2450Dataf, d, 697.0, 1209.0, ambe2450_level(99));
    pack_ambe2450_tone(d, 162, 105); /* busy */
    expect_tone("AMBE+2 busy AD 105", mbe_processAmbe2450Dataf, d, 480.0, 620.0, ambe2450_level(105));
}

static void
test_dstar_tones(void) {
    char d[49];
    pack_dstar_tone(d, 30, 204);
    expect_tone("D-STAR single vol 204", mbe_processAmbe2400Dataf, d, 937.5, 0.0, dstar_level(204));
    /* DTMF: index 128 + 4 * column + row over 697/770/852/941 by 1209/1336/1477/1633 Hz. */
    pack_dstar_tone(d, 128, 206);
    expect_tone("D-STAR DTMF 1 vol 206", mbe_processAmbe2400Dataf, d, 697.0, 1209.0, dstar_level(206));
    pack_dstar_tone(d, 128 + 4 * 2 + 3, 181);
    expect_tone("D-STAR DTMF # vol 181", mbe_processAmbe2400Dataf, d, 941.0, 1477.0, dstar_level(181));
    /* Call progress: 144..147 as AMBE+2's 160..163. */
    pack_dstar_tone(d, 146, 203);
    expect_tone("D-STAR busy vol 203", mbe_processAmbe2400Dataf, d, 480.0, 620.0, dstar_level(203));
    pack_dstar_tone(d, 147, 203);
    expect_tone("D-STAR 350+490 vol 203", mbe_processAmbe2400Dataf, d, 350.0, 490.0, dstar_level(203));
}

/* The encoder's silence frame (b0 127, index 128, volume 0) stays silence and
 * resets the decoder, as the encoder resets itself: the voice frame after it
 * decodes exactly as from a fresh stream. */
static void
test_dstar_silence_frame(void) {
    char voice[49] = {0}, silence[49] = {0};
    put_bits(voice, 0, 40, 6); /* b0 80 */
    put_bits(voice, 6, 0x15, 6);
    put_bits(voice, 12, 0x2a, 6);
    put_bits(silence, 0, 63, 6);
    silence[48] = 1;
    float out[160], fresh[160], after[160];
    mbe_parms cur, prev, enhanced;
    mbe_initMbeParms(&cur, &prev, &enhanced);
    mbe_setThreadRngSeed(1);
    for (int f = 0; f < 4; ++f) {
        assert(mbe_processAmbe2400Dataf(out, NULL, voice, &cur, &prev, &enhanced) >= 0);
    }
    mbe_process_result result;
    mbe_initProcessResult(&result);
    assert(mbe_processAmbe2400Dataf(out, &result, silence, &cur, &prev, &enhanced) >= 0);
    assert(component_db(out, 160, 697.0) < -60.0 && component_db(out, 160, 1209.0) < -60.0);
    mbe_setThreadRngSeed(2);
    assert(mbe_processAmbe2400Dataf(after, NULL, voice, &cur, &prev, &enhanced) >= 0);

    mbe_parms fcur, fprev, fenhanced;
    mbe_initMbeParms(&fcur, &fprev, &fenhanced);
    mbe_setThreadRngSeed(2);
    assert(mbe_processAmbe2400Dataf(fresh, NULL, voice, &fcur, &fprev, &fenhanced) >= 0);
    for (int i = 0; i < 160; ++i) {
        uint32_t a, b;
        memcpy(&a, &after[i], sizeof(a));
        memcpy(&b, &fresh[i], sizeof(b));
        assert(a == b);
    }
    puts("D-STAR encoder silence frame: no tone, decoder reset");
}

int
main(void) {
    test_ambe2450_tones();
    test_dstar_tones();
    test_dstar_silence_frame();
    puts("tones: ok");
    return 0;
}
