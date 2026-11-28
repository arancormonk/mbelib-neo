// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Golden hash tests for deterministic synthesis and conversion.
 */

#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "golden_sequences.h"
#include "mbe_compiler.h"
#include "mbe_cpu.h"
#include "mbelib-neo/mbelib.h"

/*
 * Exact float hashes assume each float operation rounds to float, as SSE
 * arithmetic does. 32-bit x86 builds with x87 math evaluate in extended
 * precision (FLT_EVAL_METHOD 2) and cannot reproduce them; they still check
 * the int16 hash, determinism and the sanity bounds.
 */
#if defined(MBELIB_TEST_STRICT_FLOAT) && !defined(_MSC_VER) && defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD == 0
#define GOLDEN_FLOAT_EXACT 1
#endif

/*
 * Error-free voice-only AMBE sequences (tests/golden_sequences.h). Exact hashes
 * are enforced only for x86-64 Debug builds outside MSVC: rounding differences
 * accumulate over 24 frames, so other targets check determinism and energy.
 * The float hashes depend on the kernel ISA the library runs (mbe_cpu_isa():
 * scalar, SSE2 or AVX2 vector width), whether the library was compiled for it
 * or picked it at run time; so can the int16 hashes (AMBE+2's differ between
 * SSE2 and AVX2).
 * Regenerate with gen_golden (scalar: dev-debug; SIMD: Debug with
 * -DMBELIB_ENABLE_SIMD=ON, `gen_golden sse2` and `gen_golden avx2`).
 * PFFFT 0ac3d39 uses more accurate twiddle factors. These float hashes were
 * regenerated on Ubuntu 24.04 for scalar, SSE2 and AVX2; int16 hashes held.
 * The scale-free spectral weights (#158) moved every AMBE float hash and
 * AMBE+2's int16 hashes by rounding.
 * The int16 hashes cover what the codecs' own int16 entry points return, which
 * round to nearest (mbe_floattoshort() truncates); with that, two AMBE+2
 * samples that lie within 3e-4 of a half fall on different sides of it in the
 * scalar and SIMD builds, so its scalar int16 hash differs from the SIMD one.
 * The band grid of rate-33 V/UV rows 18-31 (w0 = 2.81 / L), on 9 of the
 * AMBE+2 sequence's 24 frames, moved its float and int16 hashes.
 * Voiced phases regenerated from the decoded log2 amplitudes (log2Ml) moved
 * every AMBE hash; D-STAR's scalar int16 hash now differs from the SIMD one
 * too. The public mbe_synthesizeSpeechf() hashes below did not move.
 * The mid-frame voicing of V/UV rows 13, 15 and 17 (row 17's own voicing all
 * unvoiced, column 0 voiced at mid-frame), on the AMBE+2 sequence's one row-17
 * frame (frame 11, of 24; it has no row 13 or 15), moved its float and int16
 * hashes; the D-STAR and public hashes did not move, and P25's onset columns
 * reach no golden sequence.
 * The half-octave wave of rate 33's mid-frame pitch, between two voiced
 * frames of different pitch, moved the AMBE+2 float and int16 hashes (scalar,
 * SSE2 and AVX2; the same under glibc 2.39 and 2.44); the D-STAR and public
 * hashes did not move.
 * The mid-frame voicing of D-STAR and P25 by 500 Hz columns moved the D-STAR
 * float and int16 hashes (scalar, SSE2 and AVX2; the same under glibc 2.39 and
 * 2.44); D-STAR's scalar int16 hash equals the SIMD one again. The AMBE+2 and
 * public hashes did not move, and no golden sequence is P25's.
 * The pulse train carried through a row-17 frame in rate 33's rendering, from
 * the AMBE+2 sequence's row-17 frame (frame 11) into the frame after it, moved
 * the AMBE+2 float and int16 hashes (scalar, SSE2 and AVX2; the same under
 * glibc 2.39 and 2.44); the D-STAR and public hashes did not move.
 * Rate 33's pulsed rows (V/UV rows 18-31 played as rate 33's form of the
 * pulsed state: the row's pulsed channels, a pulse at the time the pitch code
 * gives) moved the AMBE+2 float and int16 hashes (scalar, SSE2 and AVX2): its
 * frames of those rows, 1 (row 27), 2 (20), 3 (22), 4 (28), 6 (26), 8 (27), 10
 * (28), 16 (20) and 20 (22), and the frames after them; the D-STAR and public
 * hashes did not move.
 * DVSI's low band (a third-order high-pass at 75 Hz in place of the
 * second-order one at 100 Hz, and the adaptive low cut, in D-STAR's, P25's and
 * rate 33's renderings and their unvoiced synthesis) moved the AMBE+2 and
 * D-STAR float and int16 hashes (scalar, SSE2 and AVX2): every voiced and
 * unvoiced frame of both sequences renders through the new filter. The scalar
 * int16 hashes equal the SIMD ones again in both. The public hashes did not
 * move: mbe_synthesizeSpeechf() keeps the second-order filter and no cut.
 * P25's quiet frames (every magnitude a frame of a low decoded mean log2
 * magnitude renders times 2^clamp(1.84 (x + 1.73), -2, 0), and P25's 16-bit
 * output rounded to nearest) came with three P25 sequences through the IMBE
 * 4400 data path (golden_sequences.h): a loud one, whose float hashes are those
 * of the build before the change (no frame under the knee) and whose int16
 * hash moved with the rounding; a quiet one, whose float and int16 hashes moved
 * (every frame takes the gain); and the quiet one tagged as ProVoice, whose
 * hashes are the build before's (ProVoice takes no gain and keeps truncating).
 * Their int16 hashes are the same for every kernel ISA and are checked in the
 * 32-bit Debug build too; their float hashes are SSE2's and AVX2's alike.
 * Merged onto the four rules of DVSI's P25 decoder (fix/p25-decoder, 0724a02),
 * the P25 sequences' hashes moved again: the loud one's float hashes are now
 * 0724a02's (the rules move every P25 frame; the gain none of these), its
 * int16 hash the rounding of those; the quiet one's both; the ProVoice one's
 * stayed (neither change reaches ProVoice). No other hash moved. In the 32-bit
 * Debug build (x87) the loud sequence's int16 hash is its own (0xA705BBCB):
 * the merged rendering moves 16-bit samples there; the quiet and ProVoice
 * sequences' are the 64-bit builds'.
 */
#if defined(MBE_ARCH_X86_64) && defined(GOLDEN_FLOAT_EXACT) && defined(MBELIB_TEST_STRICT_INT16)
#define GOLDEN_AMBE_STRICT 1
#else
#define GOLDEN_AMBE_STRICT 0
#endif

/* The P25 sequences' int16 hashes are checked wherever the int16 hashes are strict on x86, the 32-bit x87 Debug
 * build included (its float hashes are not: x87 evaluates otherwise). */
#if (defined(MBE_ARCH_X86_64) || defined(MBE_ARCH_X86_32)) && defined(MBELIB_TEST_STRICT_INT16)
#define GOLDEN_IMBE_STRICT_INT16 1
#else
#define GOLDEN_IMBE_STRICT_INT16 0
#endif

struct golden_ambe_expect {
    uint32_t ambe2450_f32, ambe2450_s16, ambe2400_f32, ambe2400_s16;
    uint32_t imbe_loud_f32, imbe_loud_s16, imbe_quiet_f32, imbe_quiet_s16, provoice_quiet_f32, provoice_quiet_s16;
};

static struct golden_ambe_expect
golden_ambe_expected(int isa) {
    struct golden_ambe_expect e = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};
    /* The codec paths use no kernel that differs between SSE2 and AVX2. */
    if (isa == MBE_CPU_ISA_AVX2 || isa == MBE_CPU_ISA_SSE2) {
        e.ambe2450_f32 = 0xAE5328B0u;
        e.ambe2450_s16 = 0xD2EB6902u;
        e.ambe2400_f32 = 0xECAEA43Bu;
        e.ambe2400_s16 = 0x57511F5Au;
        e.imbe_loud_f32 = 0x0AEBE54Fu;
        e.imbe_quiet_f32 = 0x0EB8EE64u;
        e.provoice_quiet_f32 = 0xE8B7AA45u;
    } else {
        e.ambe2450_f32 = 0x77F99B90u;
        e.ambe2450_s16 = 0xD2EB6902u;
        e.ambe2400_f32 = 0x5BE6D192u;
        e.ambe2400_s16 = 0x57511F5Au;
        e.imbe_loud_f32 = 0x5CCDC192u;
        e.imbe_quiet_f32 = 0x28546228u;
        e.provoice_quiet_f32 = 0xF356E69Cu;
    }
    e.imbe_loud_s16 = 0xA081D020u;
#if defined(MBE_ARCH_X86_32)
    /* x87 (the 32-bit Debug build): the four rules of DVSI's P25 decoder round some of the loud sequence's
     * intermediates otherwise in extended precision, which moves 16-bit samples; the quiet and ProVoice
     * sequences' int16 hashes are the 64-bit builds'. */
    e.imbe_loud_s16 = 0xA705BBCBu;
#endif
    e.imbe_quiet_s16 = 0x6D04AF5Du;
    e.provoice_quiet_s16 = 0xDBC888BCu;
    return e;
}

static int
float_bits_equal(float a, float b) {
    uint32_t a_bits;
    uint32_t b_bits;
    memcpy(&a_bits, &a, sizeof(a_bits));
    memcpy(&b_bits, &b, sizeof(b_bits));
    return a_bits == b_bits;
}

/**
 * @brief Check one voice-only AMBE sequence: determinism, energy, and (strict builds) exact hashes.
 * @return 0 on success, 1 on failure.
 */
static int
check_ambe_sequence(const char* name, const struct golden_ambe_codec* codec, uint32_t seed, uint32_t exp_f32,
                    uint32_t exp_s16) {
    struct golden_hashes first;
    struct golden_hashes second;

    if (golden_hash_ambe_voice_sequence(codec, seed, &first) != 0
        || golden_hash_ambe_voice_sequence(codec, seed, &second) != 0) {
        fprintf(stderr, "%s: sequence decode failed\n", name);
        return 1;
    }
    if (first.f32 != second.f32 || first.s16 != second.s16) {
        fprintf(stderr, "%s: determinism failure across runs\n", name);
        return 1;
    }
    double mean_sq = first.sumsq / (double)(GOLDEN_AMBE_FRAMES * 160);
    if (!(mean_sq > 1e-6 && mean_sq < 1e8)) {
        fprintf(stderr, "%s: sanity: mean square out of range: %g\n", name, mean_sq);
        return 1;
    }
    if (GOLDEN_AMBE_STRICT && (first.f32 != exp_f32 || first.s16 != exp_s16)) {
        fprintf(stderr, "%s: hash mismatch: f32 got=0x%08X exp=0x%08X, s16 got=0x%08X exp=0x%08X\n", name,
                (unsigned)first.f32, (unsigned)exp_f32, (unsigned)first.s16, (unsigned)exp_s16);
        return 1;
    }
    return 0;
}

/**
 * @brief Check one P25 voice sequence (golden_hash_imbe_voice_sequence()): determinism, energy, and (strict builds)
 *        exact hashes.
 * @return 0 on success, 1 on failure.
 */
static int
check_imbe_sequence(const char* name, uint32_t seed, int quiet, int provoice, uint32_t exp_f32, uint32_t exp_s16) {
    struct golden_hashes first;
    struct golden_hashes second;

    if (golden_hash_imbe_voice_sequence(seed, quiet, provoice, &first) != 0
        || golden_hash_imbe_voice_sequence(seed, quiet, provoice, &second) != 0) {
        fprintf(stderr, "%s: sequence decode failed\n", name);
        return 1;
    }
    if (first.f32 != second.f32 || first.s16 != second.s16) {
        fprintf(stderr, "%s: determinism failure across runs\n", name);
        return 1;
    }
    double mean_sq = first.sumsq / (double)(GOLDEN_IMBE_FRAMES * 160);
    if (!(mean_sq > 1e-9 && mean_sq < 1e8)) {
        fprintf(stderr, "%s: sanity: mean square out of range: %g\n", name, mean_sq);
        return 1;
    }
    if ((GOLDEN_AMBE_STRICT && first.f32 != exp_f32) || (GOLDEN_IMBE_STRICT_INT16 && first.s16 != exp_s16)) {
        fprintf(stderr, "%s: hash mismatch: f32 got=0x%08X exp=0x%08X, s16 got=0x%08X exp=0x%08X\n", name,
                (unsigned)first.f32, (unsigned)exp_f32, (unsigned)first.s16, (unsigned)exp_s16);
        return 1;
    }
    return 0;
}

/** Exit status of a test whose ISA this build or CPU cannot run (CTest SKIP_RETURN_CODE). */
#define GOLDEN_SKIP 77

/**
 * @brief Select the kernel ISA from argv: [auto|sse2|avx2] [--require].
 * @return 0 when selected, GOLDEN_SKIP when unavailable (1 with --require), 2 on a usage error.
 */
static int
select_isa(int argc, char** argv) {
    int isa = MBE_CPU_ISA_AUTO;
    int require = 0;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--require") == 0) {
            require = 1;
        } else if (strcmp(argv[i], "sse2") == 0) {
            isa = MBE_CPU_ISA_SSE2;
        } else if (strcmp(argv[i], "avx2") == 0) {
            isa = MBE_CPU_ISA_AVX2;
        } else if (strcmp(argv[i], "auto") != 0) {
            fprintf(stderr, "usage: %s [auto|sse2|avx2] [--require]\n", argv[0]);
            return 2;
        }
    }
    if (mbe_cpu_select_isa(isa) != 0) {
        fprintf(stderr, "requested kernel ISA unavailable (dispatch %s, CPU AVX2 %s, active ISA %d)\n",
                mbe_cpu_dispatch_enabled() ? "on" : "off", mbe_cpu_has_avx2() ? "yes" : "no", mbe_cpu_isa());
        return require ? 1 : GOLDEN_SKIP;
    }
    return 0;
}

/**
 * @brief Test entry: determinism + arch-specific exactness or sanity checks.
 */
int
main(int argc, char** argv) {
    int status = select_isa(argc, argv);
    if (status != 0) {
        return status;
    }
    const int isa = mbe_cpu_isa();
    printf("kernel ISA %d (dispatch %s)\n", isa, mbe_cpu_dispatch_enabled() ? "on" : "off");
    const struct golden_ambe_expect ambe = golden_ambe_expected(isa);
    if (check_ambe_sequence("ambe2450", &golden_ambe2450_codec, GOLDEN_AMBE2450_SEED, ambe.ambe2450_f32,
                            ambe.ambe2450_s16)
        || check_ambe_sequence("ambe2400", &golden_ambe2400_codec, GOLDEN_AMBE2400_SEED, ambe.ambe2400_f32,
                               ambe.ambe2400_s16)
        || check_imbe_sequence("imbe-loud", GOLDEN_IMBE_LOUD_SEED, 0, 0, ambe.imbe_loud_f32, ambe.imbe_loud_s16)
        || check_imbe_sequence("imbe-quiet", GOLDEN_IMBE_QUIET_SEED, 1, 0, ambe.imbe_quiet_f32, ambe.imbe_quiet_s16)
        || check_imbe_sequence("provoice-quiet", GOLDEN_IMBE_QUIET_SEED, 1, 1, ambe.provoice_quiet_f32,
                               ambe.provoice_quiet_s16)) {
        return 1;
    }

    /*
     * On x86/x64, keep exact golden hashes in strict mode. Float-domain
     * synthesis is deterministic per build mode, but SIMD and scalar can differ
     * due to floating-point accumulation order.
     * On other arches (e.g., AArch64/NEON), rounding order can differ; we use
     * determinism and sanity checks so CI remains green.
     */
    /* Regenerated harmonic phase, the shared spec-exact WOLA window, the
     * output high-pass phase and voiced/unvoiced render gains, and an in-band
     * fixture (golden_fill_single_frame: L * w0 < pi). */
#if (defined(MBE_ARCH_X86_64) || defined(MBE_ARCH_X86_32)) && defined(GOLDEN_FLOAT_EXACT)
    uint32_t x86_expected_f32 = 0xDDD283CAu; /* scalar */
    if (isa == MBE_CPU_ISA_AVX2) {
        x86_expected_f32 = 0x7F0B1440u;
    } else if (isa == MBE_CPU_ISA_SSE2) {
        x86_expected_f32 = 0x00280DA8u;
    }
#endif
#if (defined(MBE_ARCH_X86_64) || defined(MBE_ARCH_X86_32)) && defined(MBELIB_TEST_STRICT_INT16)
    const uint32_t X86_S16_FNV1A = 0x3EC7146Du;
#endif

    float out_f[160];
    short out_s[160];
    mbe_parms cur, prev;

    mbe_setThreadRngSeed(0xC0FFEEu);
    golden_fill_single_frame(&cur, &prev);
    /* First run */
    mbe_synthesizeSpeechf(out_f, &cur, &prev);
    uint32_t hf1 = golden_fnv1a32(out_f, sizeof(out_f));
    mbe_floattoshort(out_f, out_s);
    uint32_t hs1 = golden_fnv1a32(out_s, sizeof(out_s));
    (void)hf1; /* Hash computed for potential debug/regression use */
    (void)hs1;

    /* Determinism check: same seed/params produce identical output */
    float out_f2[160];
    short out_s2[160];
    mbe_setThreadRngSeed(0xC0FFEEu);
    golden_fill_single_frame(&cur, &prev);
    mbe_synthesizeSpeechf(out_f2, &cur, &prev);
    for (int i = 0; i < 160; ++i) {
        if (!float_bits_equal(out_f[i], out_f2[i])) {
            fprintf(stderr, "determinism failure: float sample %d differs across runs\n", i);
            return 1;
        }
    }
    mbe_floattoshort(out_f2, out_s2);
    if (memcmp(out_s, out_s2, sizeof(out_s)) != 0) {
        fprintf(stderr, "determinism failure: int16 buffers differ across runs\n");
        return 1;
    }

    /* Arch-specific exactness or sanity checks */
#if defined(MBE_ARCH_X86_64) || defined(MBE_ARCH_X86_32)
    /*
     * On x86, exact float hash can differ under aggressive FP opts. In Debug
     * or when explicitly requested, enforce exact float hash. Otherwise, use
     * sanity bounds for float but always enforce exact int16 hash.
     */
#if defined(GOLDEN_FLOAT_EXACT)
    if (hf1 != x86_expected_f32) {
        fprintf(stderr, "x86 float hash mismatch: got=0x%08X exp=0x%08X\n", (unsigned)hf1, (unsigned)x86_expected_f32);
        return 1;
    }
#else
    /* Float sanity bounds on x86 when strict disabled */
    double sumsq_f = 0.0, maxabs_f = 0.0;
    for (int i = 0; i < 160; ++i) {
        double v = out_f[i];
        double a = v < 0 ? -v : v;
        sumsq_f += v * v;
        if (a > maxabs_f) {
            maxabs_f = a;
        }
    }
    double rms_f = sumsq_f / 160.0;
    if (!(rms_f > 1e-6 && rms_f < 1e+2)) {
        fprintf(stderr, "sanity: float rms out of range: %g\n", rms_f);
        return 1;
    }
    if (!(maxabs_f < 2e4)) {
        fprintf(stderr, "sanity: float maxabs too large: %g\n", maxabs_f);
        return 1;
    }
#endif
    /* int16: strict in Debug or when requested; else sanity */
#ifdef MBELIB_TEST_STRICT_INT16
    if (hs1 != X86_S16_FNV1A) {
        fprintf(stderr, "x86 int16 hash mismatch: got=0x%08X exp=0x%08X\n", (unsigned)hs1, (unsigned)X86_S16_FNV1A);
        return 1;
    }
#else
    int maxabs_s = 0;
    for (int i = 0; i < 160; ++i) {
        int v = out_s[i];
        int a = v < 0 ? -v : v;
        if (a > maxabs_s) {
            maxabs_s = a;
        }
    }
    if (!(maxabs_s < 32000)) {
        fprintf(stderr, "sanity: int16 maxabs too large: %d\n", maxabs_s);
        return 1;
    }
#endif
#else
    /* Sanity bounds for non-x86: ensure energy and dynamic range look right */
    double sumsq_f = 0.0, maxabs_f = 0.0;
    for (int i = 0; i < 160; ++i) {
        double v = out_f[i];
        double a = v < 0 ? -v : v;
        sumsq_f += v * v;
        if (a > maxabs_f) {
            maxabs_f = a;
        }
    }
    double rms_f = sumsq_f / 160.0;
    if (!(rms_f > 1e-6 && rms_f < 1e+2)) { /* absurd values guard */
        fprintf(stderr, "sanity: float rms out of range: %g\n", rms_f);
        return 1;
    }
    if (!(maxabs_f < 2e4)) { /* must be far below clipped 32760/again */
        fprintf(stderr, "sanity: float maxabs too large: %g\n", maxabs_f);
        return 1;
    }
    /* int16 bounds */
    int maxabs_s = 0;
    for (int i = 0; i < 160; ++i) {
        int v = out_s[i];
        int a = v < 0 ? -v : v;
        if (a > maxabs_s) {
            maxabs_s = a;
        }
    }
    if (!(maxabs_s < 32000)) {
        fprintf(stderr, "sanity: int16 maxabs too large: %d\n", maxabs_s);
        return 1;
    }
#endif
    return 0;
}
