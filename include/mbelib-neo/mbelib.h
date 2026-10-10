// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Public API for mbelib-neo vocoder primitives and helpers.
 *
 * This header exposes the installed `mbe_` API used to decode IMBE/AMBE frames
 * and synthesize 8 kHz PCM audio. The 2.x API uses result-returning processing
 * calls and `mbe_process_result` status reporting; it is ABI-breaking relative
 * to 1.x, while minor releases within a major version are intended to remain
 * ABI compatible.
 *
 * @note The `*f` (float) PCM APIs return samples in mbelib's historical float
 *       scale (roughly int16/7), not normalized `[-1, +1]`. Use
 *       mbe_floattoshort() to convert to `int16_t` PCM, or scale by
 *       `(7.0f / 32768.0f)` for normalized floats (approximately `[-0.95, +0.95]`
 *       after soft clipping).
 *
 * @note Process and decode APIs validate public hard and soft bit arrays
 *       strictly. Hard bits must be exactly 0 or 1, and soft bit `bit` fields
 *       must be exactly 0 or 1. Invalid arguments return
 *       `MBE_STATUS_INVALID_ARGUMENT`; invalid bit values return
 *       `MBE_STATUS_INVALID_BITS`.
 *
 * @note Threading: processing APIs are reentrant when each stream has its own
 *       `mbe_parms` state. `mbe_setThreadRngSeed()` affects only the calling
 *       thread's synthesis RNG state.
 */

#ifndef MBELIB_NEO_PUBLIC_MBEBELIB_H
#define MBELIB_NEO_PUBLIC_MBEBELIB_H

#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

// Expose project version macro. In normal builds this header is
// generated into the build include dir; provide a fallback for linting.
#if defined(__has_include)
#if __has_include("mbelib-neo/version.h")
#include "mbelib-neo/version.h"
#endif
#endif
#ifndef MBELIB_VERSION
#define MBELIB_VERSION "0.0.0-dev"
#endif

#if !defined(MBE_API)
#if defined(_WIN32) || defined(__CYGWIN__)
/*
 * On Windows, exporting from the DLL uses dllexport when building the
 * shared library and dllimport when consuming it. For static libraries,
 * dllimport is incorrect and can cause unresolved externals; in that case
 * consumers should see an empty decoration. We propagate MBE_STATIC from
 * the CMake target for static consumption.
 */
#if defined(MBE_STATIC)
#define MBE_API
#elif defined(MBE_BUILDING)
#define MBE_API __declspec(dllexport)
#else
#define MBE_API __declspec(dllimport)
#endif
#else
#define MBE_API __attribute__((visibility("default")))
#endif
#endif

#if !defined(MBE_DEPRECATED)
#if defined(_MSC_VER)
#define MBE_DEPRECATED(msg) __declspec(deprecated(msg))
#elif defined(__GNUC__) || defined(__clang__)
#define MBE_DEPRECATED(msg) __attribute__((deprecated(msg)))
#else
#define MBE_DEPRECATED(msg)
#endif
#endif

#if !defined(MBE_DEPRECATED_FOR)
#define MBE_DEPRECATED_FOR(newsym) MBE_DEPRECATED("Use " #newsym)
#endif

struct mbe_parameters {
    /** Fundamental radian frequency (w0). */
    float w0;
    /** Number of harmonic bands (L). */
    int L;
    /** Number of voiced bands (K). */
    int K;
    /** Voiced/unvoiced flags per band (1..56). */
    int Vl[57];
    /** Magnitude per band (1..56). */
    float Ml[57];
    /** Base-2 log magnitude per band (1..56). */
    float log2Ml[57];
    /** Absolute phase per band (1..56). */
    float PHIl[57];
    /** Smoothed phase per band (1..56). */
    float PSIl[57];
    /** Spectral amplitude enhancement scale. */
    float gamma;
    /** Tone synthesis phase accumulator. */
    uint32_t tonePhase;
    /** Sine wave increment for tone synthesis. */
    int swn;

    /* === Adaptive smoothing state (Algorithms #111-116) === */
    /** Local energy tracking with IIR smoothing (Algorithm #111). */
    float localEnergy;
    /** Amplitude threshold for scaling (Algorithm #115). */
    int amplitudeThreshold;
    /** Bit error rate for current frame (0.0 to 1.0). */
    float errorRate;
    /** Total bit errors detected/corrected in this frame. */
    int errorCountTotal;
    /** Coset 4 error count (IMBE-specific). */
    int errorCount4;

    /* === Frame repeat/muting state === */
    /** Consecutive repeat count (0 to MAX_FRAME_REPEATS). */
    int repeatCount;
    /** Muting threshold for this codec (IMBE: 0.0875, AMBE: 0.096). */
    float mutingThreshold;

    /* === FFT-based unvoiced synthesis state === */
    /** Previous frame inverse FFT output for WOLA (256 samples). */
    float previousUw[256];
    /** LCG noise generator state (seed). */
    float noiseSeed;
    /** Noise buffer overlap for continuity (96 samples). */
    float noiseOverlap[96];
};

typedef struct mbe_parameters mbe_parms;

/**
 * @brief Soft-decision input bit.
 *
 * `bit` carries the caller's hard decision (0 or 1). `reliability` is the
 * confidence in that hard decision, where 0 means unknown/erasure-like and
 * 255 means highly reliable.
 */
typedef struct mbe_soft_bit {
    uint8_t bit;
    uint8_t reliability;
} mbe_soft_bit;

/** Processing/result flag: frame input used soft decisions. */
#define MBE_PROCESS_FLAG_SOFT_INPUT 0x0001u
/** Processing/result flag: C0 error count is available to the synthesis path. */
#define MBE_PROCESS_FLAG_C0_VALID   0x0002u
/** Processing/result flag: IMBE C4 error count is available. */
#define MBE_PROCESS_FLAG_C4_VALID   0x0004u
/** Processing/result flag: frame was classified as tone. */
#define MBE_PROCESS_FLAG_TONE       0x0010u
/** Processing/result flag: frame was classified as erasure. */
#define MBE_PROCESS_FLAG_ERASURE    0x0020u
/** Processing/result flag: previous parameters were repeated. */
#define MBE_PROCESS_FLAG_REPEAT     0x0040u
/** Processing/result flag: output was muted/comfort-noise substituted. */
#define MBE_PROCESS_FLAG_MUTE       0x0080u
/**
 * Processing/result flag: an AMBE 3600x2450 silence frame (b0 124/125) was
 * decoded and synthesized. Silence frames do not update the prediction
 * history (TIA-102.BABA-1 4.3). Not rendered by mbe_formatProcessResult().
 */
#define MBE_PROCESS_FLAG_SILENCE    0x0100u
/**
 * Processing/result flag (context): the parameter bits come from an IMBE
 * 7100x4400 (ProVoice) frame decode. Pass that result on to the IMBE 4400
 * data API, as for the C0/C4 context, and a muted frame gets JMBE's comfort
 * noise instead of the TIA-102.BABA 7.8 level, as in the ProVoice frame APIs.
 */
#define MBE_PROCESS_FLAG_PROVOICE   0x0200u

/** Status code: invalid pointer, invalid status counters, or otherwise unusable arguments. */
#define MBE_STATUS_INVALID_ARGUMENT (-1)
/** Status code: hard or soft input bit arrays contained values other than 0 or 1. */
#define MBE_STATUS_INVALID_BITS     (-2)

/**
 * @brief Frame decode and synthesis status.
 *
 * Decode helpers initialize and populate this structure when a non-NULL
 * pointer is supplied. Synthesis calls consume available decode context,
 * update synthesis status flags, and return `total_errors`.
 */
typedef struct mbe_process_result {
    /** Corrected errors in the C0/protected header field, when `MBE_PROCESS_FLAG_C0_VALID` is set. */
    int c0_errors;
    /** Corrected errors in protected parameter fields, excluding `c0_errors`. */
    int protected_errors;
    /** Corrected IMBE C4/Hamming errors, when `MBE_PROCESS_FLAG_C4_VALID` is set. */
    int c4_errors;
    /** Total error count used for repeat/muting decisions and status formatting. */
    int total_errors;
    /** Bitwise OR of `MBE_PROCESS_FLAG_*` values. */
    unsigned flags;
} mbe_process_result;

/** @brief Reset a process result to all-zero/default values. */
MBE_API void mbe_initProcessResult(mbe_process_result* result);
/**
 * @brief Format a process result as a compact status trace.
 *
 * Writes `'='` repeated `result->total_errors` times, then any set `E`, `T`,
 * `R`, and `M` status flags in that order, followed by NUL. Output is
 * truncated to `size`.
 */
MBE_API void mbe_formatProcessResult(char* str, size_t size, const mbe_process_result* result);
/**
 * @brief Build one soft bit from a hard bit and reliability.
 * @param bit Hard decision; any non-zero value maps to 1.
 * @param reliability Confidence in the hard decision (`0` unknown/erasure-like, `255` highly reliable).
 */
MBE_API mbe_soft_bit mbe_softBitFromHard(int bit, uint8_t reliability);
/**
 * @brief Build one soft bit from a signed LLR.
 * @param llr Signed log-likelihood-like value; positive maps to bit 1, non-positive maps to bit 0.
 * @return Soft bit with reliability equal to `abs(llr)` clamped to `0..255`.
 */
MBE_API mbe_soft_bit mbe_softBitFromLlr(int16_t llr);
/**
 * @brief Convert hard 0/1 bits to soft bits with a fixed reliability.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_softBitsFromHard(const char* bits, mbe_soft_bit* soft, size_t count, uint8_t reliability);
/**
 * @brief Convert signed LLRs to soft bits.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_softBitsFromLlr(const int16_t* llr, mbe_soft_bit* soft, size_t count);

/**
 * @brief Correct a (23,12) Golay encoded block in-place and extract data.
 * @param block Pointer to packed 23-bit block (upper bits ignored). On return, contains 12-bit data.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_checkGolayBlock(long int* block);
/**
 * @brief Decode a (23,12) Golay codeword.
 * @param in  Input bits, LSB at index 0, length 23.
 * @param out Output bits, corrected, LSB at index 0, length 23.
 * @return Number of corrected bit errors in the protected portion.
 */
MBE_API int mbe_golay2312(const char* in, char* out);
/**
 * @brief Soft-decision Golay(23,12) decode.
 * @param in  Input soft bits, LSB at index 0, length 23.
 * @param out Output hard bits, length 23.
 * @return Number of hard-decision data-bit differences between `in` and the selected codeword.
 * @note Matches mbe_golay2312 output semantics: data bits are corrected, parity bits preserve the hard input.
 */
MBE_API int mbe_golay2312Soft(const mbe_soft_bit* in, char* out);
/**
 * @brief Decode a (15,11) Hamming codeword (IMBE/AMBE common use).
 * @param in  Input bits, LSB at index 0, length 15.
 * @param out Output bits, corrected, LSB at index 0, length 15.
 * @return Number of corrected bit errors (0 or 1).
 */
MBE_API int mbe_hamming1511(const char* in, char* out);
/**
 * @brief Soft-decision Hamming(15,11) decode.
 * @param in  Input soft bits, LSB at index 0, length 15.
 * @param out Output corrected hard bits, LSB at index 0, length 15.
 * @return Number of hard-decision bit differences between `in` and the selected codeword.
 */
MBE_API int mbe_hamming1511Soft(const mbe_soft_bit* in, char* out);
/**
 * @brief Decode a (15,11) Hamming codeword with IMBE 7100x4400 mapping.
 * @param in  Input bits, LSB at index 0, length 15.
 * @param out Output bits, corrected, LSB at index 0, length 15.
 * @return Number of corrected bit errors (0 or 1).
 */
MBE_API int mbe_7100x4400hamming1511(const char* in, char* out);
/**
 * @brief Soft-decision Hamming(15,11) decode with IMBE 7100x4400 mapping.
 * @param in  Input soft bits, LSB at index 0, length 15.
 * @param out Output corrected hard bits, LSB at index 0, length 15.
 * @return Number of hard-decision bit differences between `in` and the selected codeword.
 */
MBE_API int mbe_7100x4400hamming1511Soft(const mbe_soft_bit* in, char* out);

/* Prototypes from ambe3600x2400.c */
/** @brief Print AMBE 2400 parameter bits to stderr (debug). */
MBE_API void mbe_dumpAmbe2400Data(const char* ambe_d);
/** @brief Print a raw AMBE 3600x2400 frame to stderr (debug). */
MBE_API void mbe_dumpAmbe3600x2400Frame(const char ambe_fr[4][24]);
/**
 * @brief Apply ECC to AMBE 3600x2400 C0 and update in-place.
 * @param ambe_fr AMBE frame as 4x24 bitplanes.
 * @return Number of corrected errors in C0.
 */
MBE_API int mbe_eccAmbe3600x2400C0(char ambe_fr[4][24]);
/**
 * @brief Apply ECC to AMBE 3600x2400 data and pack parameters.
 * @param ambe_fr AMBE frame as 4x24 bitplanes.
 * @param ambe_d  Output parameter bits (49).
 * @return Number of corrected errors in protected fields.
 */
MBE_API int mbe_eccAmbe3600x2400Data(char ambe_fr[4][24], char* ambe_d);
/**
 * @brief Decode AMBE 2400 parameters from demodulated bits.
 * @param ambe_d  Demodulated AMBE parameter bits (49).
 * @param cur_mp  Output: current frame parameters.
 * @param prev_mp Input: previous frame parameters (for prediction).
 * @return Tone index or 0 for voice; implementation-specific non-zero for tone frames.
 */
MBE_API int mbe_decodeAmbe2400Parms(const char* ambe_d, mbe_parms* cur_mp, mbe_parms* prev_mp);
/**
 * @brief Demodulate interleaved AMBE 3600x2400 data in-place.
 * @param ambe_fr AMBE frame as 4x24 bitplanes, updated in-place.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_demodulateAmbe3600x2400Data(char ambe_fr[4][24]);
/**
 * @brief Decode a hard AMBE 3600x2400 frame to parameter bits without synthesis.
 * @param ambe_fr Input frame as 4x24 bitplanes; not modified.
 * @param ambe_d  Output parameter bits (49).
 * @param result  Optional output status; receives C0/protected/total errors and `MBE_PROCESS_FLAG_C0_VALID`.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeAmbe3600x2400Frame(const char ambe_fr[4][24], char ambe_d[49], mbe_process_result* result);
/**
 * @brief Decode a soft AMBE 3600x2400 frame to parameter bits without synthesis.
 * @param ambe_fr Input soft frame as 4x24 bitplanes; not modified.
 * @param ambe_d  Output hard parameter bits (49).
 * @param result  Optional output status; also sets `MBE_PROCESS_FLAG_SOFT_INPUT`.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeAmbe3600x2400SoftFrame(const mbe_soft_bit ambe_fr[4][24], char ambe_d[49],
                                             mbe_process_result* result);
/**
 * @brief Process AMBE 2400 parameters into 8 kHz float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional in/out status context. C0 context is used when `MBE_PROCESS_FLAG_C0_VALID` is set.
 * @param ambe_d   Demodulated parameter bits (49).
 * @param cur_mp   In/out: current frame parameters (may be enhanced).
 * @param prev_mp  In/out: previous frame parameters.
 * @param prev_mp_enhanced In/out: enhanced previous parameters for continuity.
 * @return Total error count on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_processAmbe2400Dataf(float* aout_buf, mbe_process_result* result, const char ambe_d[49],
                                     mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/**
 * @brief Process AMBE 2400 parameters into 8 kHz 16-bit PCM.
 * @see mbe_processAmbe2400Dataf for details.
 */
MBE_API int mbe_processAmbe2400Data(short* aout_buf, mbe_process_result* result, const char ambe_d[49],
                                    mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/* === AMBE 3600x2400 (D-STAR) encoding === */

/**
 * @brief Caller-owned AMBE 2400 analysis state.
 *
 * Use one context per stream. A context is not thread-safe for concurrent use;
 * any number of independent contexts may be used in one thread.
 */
typedef struct mbe_ambe2400_encoder mbe_ambe2400_encoder;

/**
 * @brief Allocate a fresh encoder context, including its FFT plan.
 *
 * Encoding never allocates. Allocation failures inside the vendored pffft
 * setup are not recoverable; the same limitation applies to decoder plan
 * allocation.
 *
 * @return Owned context, or NULL when the context or its FFT plan buffers
 * cannot be allocated.
 * @see mbe_ambe2400EncoderFree
 */
MBE_API mbe_ambe2400_encoder* mbe_ambe2400EncoderAlloc(void);

/**
 * @brief Restore freshly allocated analysis state, retaining the FFT plan.
 * @param enc Context to reset; NULL is accepted.
 * @see mbe_encodeAmbe2400Parms for resetting the caller's prediction state.
 */
MBE_API void mbe_ambe2400EncoderReset(mbe_ambe2400_encoder* enc);

/**
 * @brief Free an encoder context and its FFT plan; NULL is accepted.
 */
MBE_API void mbe_ambe2400EncoderFree(mbe_ambe2400_encoder* enc);

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of float PCM into AMBE 2400
 *        parameter bits.
 *
 * Bit-compatible with this library's mbe_decodeAmbe2400Parms()/
 * mbe_processAmbe3600x2400*() path and following the D-STAR AMBE bit layout
 * (interleave, scrambler and Golay parity cross-checked against the MMDVM
 * tables). The spectral reconstruction it targets matches DVSI's AMBE-3000
 * D-STAR test vectors; on-air interoperability has not been verified.
 * As from DVSI's encoder, quiet and silent input is coded as low-level voice,
 * and the decoded level follows the input level (there is no AGC). Earlier
 * versions normalized the level and sent the AMBE silence frame (b0 127, tone
 * index 128) for quiet input; this library's decoders play that frame as
 * comfort noise and reset. DTMF, call-progress and single tones are sent as
 * D-STAR tone frames (b0 126, with the tone index and 8-bit volume the decoder
 * plays), where DVSI's encoder sends them. A call-progress tone is sent once
 * it has filled three consecutive 20 ms analysis spans (an 80 ms tone always
 * does; DVSI starts one frame earlier) and for up to two frames after it ends,
 * which a stream that ends with the tone only gets if it encodes them. A DTMF
 * or single tone is sent from its first detection and for one frame after
 * its last, as DVSI's encoder sends the frame in which such a tone ends or
 * changes. KNOX tones, which D-STAR cannot carry, are coded as voice.
 *
 * Initialize with mbe_ambe2400EncoderAlloc() and mbe_initMbeParms(), then
 * advance prediction state with mbe_moveMbeParms(cur_mp, prev_mp) between
 * frames. prev_mp is read-only. To restart a stream, call
 * mbe_ambe2400EncoderReset() and mbe_initMbeParms().
 * State equivalence holds with both the mbe_processAmbe2400* path and a bare
 * mbe_decodeAmbe2400Parms() chain. A tone frame leaves the decoder's
 * prediction history unchanged, so cur_mp is then a copy of prev_mp.
 *
 * The analysis is centred one sample before the frame: parameters lag audio
 * by about 10 ms. The pitch analysis spans the whole frame, while the spectral
 * analysis reaches its first 110 samples and covers the rest with the next
 * call. Feed a final frame of zeros to flush the tail.
 *
 * @param enc     Caller-owned context; NULL returns MBE_STATUS_INVALID_ARGUMENT.
 * @param samples Input PCM floats (160), nominal range [-1, 1]. A frame with a
 *                non-finite sample or one beyond +-2^20 returns
 *                MBE_STATUS_INVALID_ARGUMENT and leaves the context unchanged.
 * @param ambe_d  Output parameter bits (49). ambe_d[24] is the spare bit.
 * @param cur_mp  Output: quantized (decoder-equivalent) parameters.
 * @param prev_mp Input: previous quantized frame state; never modified. A
 *                state with a log2Ml or gamma that is not finite or lies
 *                beyond +-2^20, or one from which the model would approach
 *                overflow (a log2Ml above 120; decoders stay below 50),
 *                returns MBE_STATUS_INVALID_ARGUMENT and leaves the context
 *                and cur_mp unchanged.
 * @return 0, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_encodeAmbe2400Parms(mbe_ambe2400_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                                    const mbe_parms* prev_mp);
/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of 16-bit PCM into AMBE 2400
 *        parameter bits.
 * @see mbe_encodeAmbe2400Parms for details.
 */
MBE_API int mbe_encodeAmbe2400ParmsShort(mbe_ambe2400_encoder* enc, const short* samples, char ambe_d[49],
                                         mbe_parms* cur_mp, const mbe_parms* prev_mp);
/**
 * @brief Encode 49 AMBE 2400 parameter bits into a 72-bit D-STAR DV data
 *        frame (FEC + interleave), in the decoder's plane layout.
 *
 * ambe_d[24] is the spare bit; on output it carries the scrambled even
 * parity of the second Golay codeword. All other input bits round-trip
 * exactly through
 * mbe_decodeAmbe3600x2400Frame().
 *
 * @param ambe_d  Input parameter bits (49).
 * @param ambe_fr Output frame as 4x24 bitplanes.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_encodeAmbe3600x2400Frame(const char ambe_d[49], char ambe_fr[4][24]);

/* === D-STAR DV framing === */

/**
 * @brief Serialize a 72-bit AMBE 3600x2400 frame into the 9 data bytes
 *        of a D-STAR DV frame (sync word not included).
 *
 * Bytes are packed in air order (LSB first within each byte, matching
 * the GMSK modulator). Inverse of mbe_decodeDStarDVData().
 *
 * @param ambe_fr Input frame as 4x24 bitplanes.
 * @param bytes9  Output 9 bytes (72 bits).
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_encodeDStarDVData(const char ambe_fr[4][24], unsigned char bytes9[9]);
/**
 * @brief Extract a 72-bit AMBE 3600x2400 frame from the 9 data bytes of
 *        a D-STAR DV frame (sync word not included).
 * @see mbe_encodeDStarDVData for the byte/bit convention.
 *
 * @param bytes9  Input 9 bytes (72 bits).
 * @param ambe_fr Output frame as 4x24 bitplanes.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_decodeDStarDVData(const unsigned char bytes9[9], char ambe_fr[4][24]);

/* === AMBE 3600x2400 frame processing === */

/**
 * @brief Process a complete AMBE 3600x2400 frame into 8 kHz float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional output status populated by decode and synthesis.
 * @param ambe_fr  Input frame as 4x24 bitplanes.
 * @param ambe_d   Scratch/output parameter bits (49).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 * @return Total error count on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_processAmbe3600x2400Framef(float* aout_buf, mbe_process_result* result, const char ambe_fr[4][24],
                                           char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                           mbe_parms* prev_mp_enhanced);
/**
 * @brief Process a complete AMBE 3600x2400 frame into 8 kHz 16-bit PCM.
 * @see mbe_processAmbe3600x2400Framef for details.
 */
MBE_API int mbe_processAmbe3600x2400Frame(short* aout_buf, mbe_process_result* result, const char ambe_fr[4][24],
                                          char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                          mbe_parms* prev_mp_enhanced);
/**
 * @brief Process a soft AMBE 3600x2400 frame into float PCM.
 * @param result Optional output status populated by decode and synthesis.
 * @return Final total error count used by the wrapper.
 */
MBE_API int mbe_processAmbe3600x2400SoftFramef(float* aout_buf, mbe_process_result* result,
                                               const mbe_soft_bit ambe_fr[4][24], char ambe_d[49], mbe_parms* cur_mp,
                                               mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/** @brief Process a soft AMBE 3600x2400 frame into int16 PCM; see float variant for result semantics. */
MBE_API int mbe_processAmbe3600x2400SoftFrame(short* aout_buf, mbe_process_result* result,
                                              const mbe_soft_bit ambe_fr[4][24], char ambe_d[49], mbe_parms* cur_mp,
                                              mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/* Prototypes from ambe3600x2450.c */
/** @brief Print AMBE 2450 parameter bits to stderr (debug). */
MBE_API void mbe_dumpAmbe2450Data(const char* ambe_d);
/** @brief Print a raw AMBE 3600x2450 frame to stderr (debug). */
MBE_API void mbe_dumpAmbe3600x2450Frame(const char ambe_fr[4][24]);
/** @brief ECC correction for AMBE 3600x2450 C0. */
MBE_API int mbe_eccAmbe3600x2450C0(char ambe_fr[4][24]);
/** @brief ECC and parameter packing for AMBE 3600x2450. */
MBE_API int mbe_eccAmbe3600x2450Data(char ambe_fr[4][24], char* ambe_d);
/** AMBE 2450 frame type: voice frame; the only type that may be committed as prediction history. */
#define MBE_AMBE2450_FRAME_VOICE   0
/** AMBE 2450 frame type: silence frame (b0 124/125); decoded but never committed to prev_mp. */
#define MBE_AMBE2450_FRAME_SILENCE 1
/**
 * AMBE 2450 frame type: erasure (b0 120-123, or a tone frame whose tone index is
 * invalid or whose redundant fields disagree, TIA-102.BABA-1 7.3); a frame repeat.
 */
#define MBE_AMBE2450_FRAME_ERASURE 2
/** AMBE 2450 frame type: tone frame with a usable tone index (TIA-102.BABA-1 7, 7.3). */
#define MBE_AMBE2450_FRAME_TONE    7
/**
 * @brief Classify AMBE 2450 parameter bits without decoding them.
 *
 * A tone frame is recognised by the first six bits of u0 equal to 63
 * (TIA-102.BABA-1 7), before b0 is read; otherwise b0 gives the type (4.1).
 *
 * @param ambe_d Demodulated parameter bits (49).
 * @return An `MBE_AMBE2450_FRAME_*` value, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_classifyAmbe2450Frame(const char ambe_d[49]);
/**
 * @brief Decode AMBE 2450 model parameters (no repeat, mute or synthesis handling).
 *
 * Predicts spectral amplitudes from `prev_mp`, which must hold the last valid
 * voice frame (TIA-102.BABA-1 4.4.1 eq 26, 4.4.3 eq 43). Voice and silence
 * frames both decode into `cur_mp` and return `MBE_AMBE2450_FRAME_VOICE`, as
 * in 2.1. A caller managing its own state should copy `cur_mp` into `prev_mp`
 * only when `mbe_classifyAmbe2450Frame()` reports `MBE_AMBE2450_FRAME_VOICE`;
 * silence, erasure, tone and repeated frames must not replace the prediction
 * history. The decode writes the eq 44-45 edge values (index 0 and above
 * `prev_mp->L`) into `prev_mp`, which does not change its history.
 *
 * @param ambe_d  Demodulated parameter bits (49).
 * @param cur_mp  Output: current frame parameters (voice and silence frames).
 * @param prev_mp In/out: last valid voice frame (prediction history).
 * @return `MBE_AMBE2450_FRAME_VOICE` (voice or silence), `MBE_AMBE2450_FRAME_ERASURE`
 *         or `MBE_AMBE2450_FRAME_TONE`, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_decodeAmbe2450Parms(const char* ambe_d, mbe_parms* cur_mp, mbe_parms* prev_mp);
/** @brief Demodulate AMBE 3600x2450 interleaved data. */
MBE_API int mbe_demodulateAmbe3600x2450Data(char ambe_fr[4][24]);
/**
 * @brief Decode a hard AMBE 3600x2450 frame to parameter bits without synthesis.
 * @param ambe_fr Input frame as 4x24 bitplanes; not modified.
 * @param ambe_d  Output parameter bits (49).
 * @param result  Optional output status; receives C0/protected/total errors and `MBE_PROCESS_FLAG_C0_VALID`.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeAmbe3600x2450Frame(const char ambe_fr[4][24], char ambe_d[49], mbe_process_result* result);
/**
 * @brief Decode a soft AMBE 3600x2450 frame to parameter bits without synthesis.
 * @param ambe_fr Input soft frame as 4x24 bitplanes; not modified.
 * @param ambe_d  Output hard parameter bits (49).
 * @param result  Optional output status; also sets `MBE_PROCESS_FLAG_SOFT_INPUT`.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeAmbe3600x2450SoftFrame(const mbe_soft_bit ambe_fr[4][24], char ambe_d[49],
                                             mbe_process_result* result);
/**
 * @brief Process AMBE 2450 parameters into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional in/out status context. C0 context is used when `MBE_PROCESS_FLAG_C0_VALID` is set.
 * @param ambe_d   Demodulated parameter bits (49).
 * @param cur_mp   In/out: current frame parameters (may be enhanced).
 * @param prev_mp  In/out: last valid voice frame (prediction history) plus per-frame error/repeat state.
 * @param prev_mp_enhanced In/out: last synthesized frame; also the frame-repeat source.
 * @return Total error count on success, or a negative `MBE_STATUS_*` code.
 *
 * Frame types follow TIA-102.BABA-1: silence frames are synthesized but do not
 * update the prediction history (`MBE_PROCESS_FLAG_SILENCE`); erasures (a tone
 * frame with an unusable index also sets `MBE_PROCESS_FLAG_TONE`) and frames
 * meeting the 5.6 repeat criteria repeat the last synthesized frame
 * (`MBE_PROCESS_FLAG_REPEAT`); the error rate above 0.096
 * or a 4th consecutive invalid frame mutes (`MBE_PROCESS_FLAG_MUTE`). The three
 * parameter sets must be distinct objects; aliased sets return
 * `MBE_STATUS_INVALID_ARGUMENT`.
 */
MBE_API int mbe_processAmbe2450Dataf(float* aout_buf, mbe_process_result* result, const char ambe_d[49],
                                     mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/** @brief Process AMBE 2450 parameters into 16-bit PCM. */
MBE_API int mbe_processAmbe2450Data(short* aout_buf, mbe_process_result* result, const char ambe_d[49],
                                    mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/**
 * @brief Process AMBE 3600x2450 frame into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional output status populated by decode and synthesis.
 * @param ambe_fr  Input frame as 4x24 bitplanes.
 * @param ambe_d   Scratch/output parameter bits (49).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 * @return Total error count on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_processAmbe3600x2450Framef(float* aout_buf, mbe_process_result* result, const char ambe_fr[4][24],
                                           char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                           mbe_parms* prev_mp_enhanced);
/** @brief Process AMBE 3600x2450 frame into 16-bit PCM. */
MBE_API int mbe_processAmbe3600x2450Frame(short* aout_buf, mbe_process_result* result, const char ambe_fr[4][24],
                                          char ambe_d[49], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                          mbe_parms* prev_mp_enhanced);
/**
 * @brief Process a soft AMBE 3600x2450 frame into float PCM.
 * @param result Optional output status populated by decode and synthesis.
 * @return Final total error count used by the wrapper.
 */
MBE_API int mbe_processAmbe3600x2450SoftFramef(float* aout_buf, mbe_process_result* result,
                                               const mbe_soft_bit ambe_fr[4][24], char ambe_d[49], mbe_parms* cur_mp,
                                               mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/** @brief Process a soft AMBE 3600x2450 frame into int16 PCM; see float variant for result semantics. */
MBE_API int mbe_processAmbe3600x2450SoftFrame(short* aout_buf, mbe_process_result* result,
                                              const mbe_soft_bit ambe_fr[4][24], char ambe_d[49], mbe_parms* cur_mp,
                                              mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/* === AMBE+2 3600x2450 (DMR, NXDN, YSF, P25 Phase 2) encoding === */

/**
 * @brief Caller-owned AMBE+2 2450 analysis state.
 *
 * Use one context per stream. A context is not thread-safe for concurrent use;
 * any number of independent contexts may be used in one thread.
 */
typedef struct mbe_ambe2450_encoder mbe_ambe2450_encoder;

/**
 * @brief Allocate a fresh encoder context, including its FFT plans.
 *
 * Encoding never allocates. As for mbe_ambe2400EncoderAlloc(), allocation
 * failures inside the vendored pffft setup are not recoverable.
 *
 * @return Owned context, or NULL when the context or its plan buffers cannot
 * be allocated.
 * @see mbe_ambe2450EncoderFree
 */
MBE_API mbe_ambe2450_encoder* mbe_ambe2450EncoderAlloc(void);

/**
 * @brief Restore freshly allocated analysis state, retaining the FFT plans.
 * @param enc Context to reset; NULL is accepted.
 * @see mbe_encodeAmbe2450Parms for resetting the caller's prediction state.
 */
MBE_API void mbe_ambe2450EncoderReset(mbe_ambe2450_encoder* enc);

/** @brief Free an encoder context and its FFT plans; NULL is accepted. */
MBE_API void mbe_ambe2450EncoderFree(mbe_ambe2450_encoder* enc);

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of float PCM into AMBE+2 2450
 *        parameter bits.
 *
 * Quantizes as TIA-102.BABA-1 clause 4 describes, against the tables
 * mbe_decodeAmbe2450Parms() dequantizes with; as in the D-STAR encoder, a
 * frame quieter than the lowest gain step can play is flattened toward that
 * step's level rather than played louder. DTMF, KNOX, call-progress and
 * single tones (TIA-102.BABA-1 Table 9) are sent as tone frames (7.2), timed
 * as in the D-STAR encoder: a call-progress tone once it has filled three
 * consecutive analysis spans and for up to two frames after it ends, and a
 * DTMF, KNOX or single tone from its first detection and for one frame after
 * its last. Like
 * DVSI's encoder it never sends silence or erasure frames: quiet input is
 * coded as low-level voice.
 *
 * The state contract is the D-STAR encoder's: initialize with
 * mbe_ambe2450EncoderAlloc() and mbe_initMbeParms(), and after every frame
 * call mbe_moveMbeParms(cur_mp, prev_mp). For a voice frame cur_mp is what
 * mbe_decodeAmbe2450Parms() decodes from ambe_d given prev_mp. A tone frame
 * leaves the decoder's prediction history unchanged (TIA-102.BABA-1 4.3), so
 * cur_mp is then a copy of prev_mp: a snapshot of the history, not a model of
 * the tone. To hear what was sent, run the emitted bits through
 * mbe_processAmbe2450Data() with a separate decoder state. To restart a
 * stream, call mbe_ambe2450EncoderReset() and mbe_initMbeParms().
 *
 * The timing is the D-STAR encoder's: parameters lag audio by about 10 ms,
 * and a final frame of zeros flushes the tail.
 *
 * @param enc     Caller-owned context; NULL returns MBE_STATUS_INVALID_ARGUMENT.
 * @param samples Input PCM floats (160), nominal range [-1, 1]. A frame with a
 *                non-finite sample or one beyond +-2^20 returns
 *                MBE_STATUS_INVALID_ARGUMENT and leaves the context unchanged.
 * @param ambe_d  Output parameter bits (49).
 * @param cur_mp  Output: the decoder's parameters for this frame (see above);
 *                must not be prev_mp.
 * @param prev_mp Input: previous frame state; never modified. A state with a
 *                log2Ml or gamma that is not finite or lies beyond +-2^20, or
 *                one from which the model would approach overflow (a log2Ml
 *                above 120, or a non-finite w0 or amplitude; decoders stay
 *                below 50), returns MBE_STATUS_INVALID_ARGUMENT and leaves the
 *                context and cur_mp unchanged.
 * @return MBE_AMBE2450_FRAME_VOICE or MBE_AMBE2450_FRAME_TONE, or a negative
 *         `MBE_STATUS_*` code.
 */
MBE_API int mbe_encodeAmbe2450Parms(mbe_ambe2450_encoder* enc, const float* samples, char ambe_d[49], mbe_parms* cur_mp,
                                    const mbe_parms* prev_mp);
/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of 16-bit PCM into AMBE+2 2450
 *        parameter bits.
 * @see mbe_encodeAmbe2450Parms for details.
 */
MBE_API int mbe_encodeAmbe2450ParmsShort(mbe_ambe2450_encoder* enc, const short* samples, char ambe_d[49],
                                         mbe_parms* cur_mp, const mbe_parms* prev_mp);
/**
 * @brief Encode 49 AMBE+2 2450 parameter bits into a 72-bit AMBE 3600x2450
 *        frame (TIA-102.BABA-1 5.2-5.3: Golay codes and the modulation of
 *        C1), in the plane layout mbe_decodeAmbe3600x2450Frame() reads.
 *
 * All 49 bits round-trip through mbe_decodeAmbe3600x2450Frame().
 *
 * @param ambe_d  Input parameter bits (49).
 * @param ambe_fr Output frame as 4x24 bitplanes.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_encodeAmbe3600x2450Frame(const char ambe_d[49], char ambe_fr[4][24]);

/* Prototypes from imbe7200x4400.c */
/** @brief Print IMBE 4400 parameter bits to stderr (debug). */
MBE_API void mbe_dumpImbe4400Data(const char* imbe_d);
/** @brief Print IMBE 7200x4400 parameter bits to stderr (debug). */
MBE_API void mbe_dumpImbe7200x4400Data(const char* imbe_d);
/** @brief Print a raw IMBE 7200x4400 frame to stderr (debug). */
MBE_API void mbe_dumpImbe7200x4400Frame(const char imbe_fr[8][23]);
/** @brief ECC correction for IMBE 7200x4400 C0. */
MBE_API int mbe_eccImbe7200x4400C0(char imbe_fr[8][23]);
/** @brief ECC and parameter packing for IMBE 7200x4400. */
MBE_API int mbe_eccImbe7200x4400Data(char imbe_fr[8][23], char* imbe_d);
/** @brief Decode IMBE 4400 parameters. */
MBE_API int mbe_decodeImbe4400Parms(const char* imbe_d, mbe_parms* cur_mp, mbe_parms* prev_mp);
/** @brief Demodulate IMBE 7200x4400 interleaved data. */
MBE_API int mbe_demodulateImbe7200x4400Data(char imbe[8][23]);
/**
 * @brief Decode a hard IMBE 7200x4400 frame to parameter bits without synthesis.
 * @param imbe_fr Input frame as 8x23 bitplanes; not modified.
 * @param imbe_d  Output parameter bits (88).
 * @param result  Optional output status; receives C0/protected/C4/total errors and valid-context flags.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeImbe7200x4400Frame(const char imbe_fr[8][23], char imbe_d[88], mbe_process_result* result);
/**
 * @brief Decode a soft IMBE 7200x4400 frame to parameter bits without synthesis.
 * @param imbe_fr Input soft frame as 8x23 bitplanes; not modified.
 * @param imbe_d  Output hard parameter bits (88).
 * @param result  Optional output status; also sets `MBE_PROCESS_FLAG_SOFT_INPUT`.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeImbe7200x4400SoftFrame(const mbe_soft_bit imbe_fr[8][23], char imbe_d[88],
                                             mbe_process_result* result);
/**
 * @brief Process IMBE 4400 parameters into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional in/out status context. C0/C4 context is used when the matching valid flags are set.
 * @param imbe_d   Demodulated parameter bits (88).
 * @param cur_mp   In/out: current frame parameters (may be enhanced).
 * @param prev_mp  In/out: previous frame parameters.
 * @param prev_mp_enhanced In/out: enhanced previous parameters for continuity.
 * @return Total error count on success, or a negative `MBE_STATUS_*` code.
 *
 * A muted frame outputs the TIA-102.BABA 7.8 noise (P25), or JMBE's comfort
 * noise when `result` carries `MBE_PROCESS_FLAG_PROVOICE` from a 7100x4400 decode.
 */
MBE_API int mbe_processImbe4400Dataf(float* aout_buf, mbe_process_result* result, const char imbe_d[88],
                                     mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/** @brief Process IMBE 4400 parameters into 16-bit PCM. */
MBE_API int mbe_processImbe4400Data(short* aout_buf, mbe_process_result* result, const char imbe_d[88],
                                    mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/**
 * @brief Process IMBE 7200x4400 frame into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional output status populated by decode and synthesis.
 * @param imbe_fr  Input frame as 8x23 bitplanes.
 * @param imbe_d   Scratch/output parameter bits (88).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 * @return Total error count on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_processImbe7200x4400Framef(float* aout_buf, mbe_process_result* result, const char imbe_fr[8][23],
                                           char imbe_d[88], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                           mbe_parms* prev_mp_enhanced);
/** @brief Process IMBE 7200x4400 frame into 16-bit PCM. */
MBE_API int mbe_processImbe7200x4400Frame(short* aout_buf, mbe_process_result* result, const char imbe_fr[8][23],
                                          char imbe_d[88], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                          mbe_parms* prev_mp_enhanced);
/**
 * @brief Process a soft IMBE 7200x4400 frame into float PCM.
 * @param result Optional output status populated by decode and synthesis.
 * @return Final total error count used by the wrapper.
 */
MBE_API int mbe_processImbe7200x4400SoftFramef(float* aout_buf, mbe_process_result* result,
                                               const mbe_soft_bit imbe_fr[8][23], char imbe_d[88], mbe_parms* cur_mp,
                                               mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/** @brief Process a soft IMBE 7200x4400 frame into int16 PCM; see float variant for result semantics. */
MBE_API int mbe_processImbe7200x4400SoftFrame(short* aout_buf, mbe_process_result* result,
                                              const mbe_soft_bit imbe_fr[8][23], char imbe_d[88], mbe_parms* cur_mp,
                                              mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/* === IMBE 7200x4400 (P25 Phase 1 full rate) encoding === */

/**
 * @brief Caller-owned IMBE 4400 analysis state.
 *
 * Use one context per stream. A context is not thread-safe for concurrent use;
 * any number of independent contexts may be used in one thread.
 */
typedef struct mbe_imbe4400_encoder mbe_imbe4400_encoder;

/**
 * @brief Allocate a fresh encoder context, including its FFT plans.
 *
 * Encoding never allocates. As for mbe_ambe2400EncoderAlloc(), allocation
 * failures inside the vendored pffft setup are not recoverable.
 *
 * @return Owned context, or NULL when the context or its plan buffers cannot
 * be allocated.
 * @see mbe_imbe4400EncoderFree
 */
MBE_API mbe_imbe4400_encoder* mbe_imbe4400EncoderAlloc(void);

/**
 * @brief Restore freshly allocated analysis state, retaining the FFT plans;
 *        the next frame's synchronization bit is 0 again.
 * @param enc Context to reset; NULL is accepted.
 * @see mbe_encodeImbe4400Parms for resetting the caller's prediction state.
 */
MBE_API void mbe_imbe4400EncoderReset(mbe_imbe4400_encoder* enc);

/** @brief Free an encoder context and its FFT plans; NULL is accepted. */
MBE_API void mbe_imbe4400EncoderFree(mbe_imbe4400_encoder* enc);

/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of float PCM into IMBE 4400
 *        parameter bits.
 *
 * Quantizes as TIA-102.BABA chapter 6 describes, against the tables
 * mbe_decodeImbe4400Parms() dequantizes with; as in the D-STAR encoder, a
 * frame quieter than the lowest gain level can play is flattened toward that
 * level rather than played louder. imbe_d[87] is the
 * synchronization bit of 6.5: 0 in the first frame after allocation or reset,
 * then alternating. Every frame is a voice frame.
 *
 * The state contract is the D-STAR encoder's: initialize with
 * mbe_imbe4400EncoderAlloc() and mbe_initMbeParms(), and after every frame
 * call mbe_moveMbeParms(cur_mp, prev_mp); cur_mp is what
 * mbe_decodeImbe4400Parms() decodes from imbe_d given prev_mp. To restart a
 * stream, call mbe_imbe4400EncoderReset() and mbe_initMbeParms().
 *
 * The timing is the D-STAR encoder's: parameters lag audio by about 10 ms,
 * and a final frame of zeros flushes the tail.
 *
 * @param enc     Caller-owned context; NULL returns MBE_STATUS_INVALID_ARGUMENT.
 * @param samples Input PCM floats (160), nominal range [-1, 1]. A frame with a
 *                non-finite sample or one beyond +-2^20 returns
 *                MBE_STATUS_INVALID_ARGUMENT and leaves the context unchanged.
 * @param imbe_d  Output parameter bits (88).
 * @param cur_mp  Output: the decoder's parameters for this frame; must not be
 *                prev_mp.
 * @param prev_mp Input: previous frame state; never modified. A state with a
 *                log2Ml or gamma that is not finite or lies beyond +-2^20, or
 *                one from which the model would approach overflow (a log2Ml
 *                above 120; decoders stay below 50), returns
 *                MBE_STATUS_INVALID_ARGUMENT and leaves the context and cur_mp
 *                unchanged.
 * @return 0, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_encodeImbe4400Parms(mbe_imbe4400_encoder* enc, const float* samples, char imbe_d[88], mbe_parms* cur_mp,
                                    const mbe_parms* prev_mp);
/**
 * @brief Encode 160 samples (20 ms, 8 kHz) of 16-bit PCM into IMBE 4400
 *        parameter bits.
 * @see mbe_encodeImbe4400Parms for details.
 */
MBE_API int mbe_encodeImbe4400ParmsShort(mbe_imbe4400_encoder* enc, const short* samples, char imbe_d[88],
                                         mbe_parms* cur_mp, const mbe_parms* prev_mp);
/**
 * @brief Encode 88 IMBE 4400 parameter bits into a 144-bit IMBE 7200x4400
 *        frame (TIA-102.BABA chapter 7: four (23,12) Golay and three (15,11)
 *        Hamming codes, seven unprotected bits, and the modulation of
 *        vectors 1-6), in the plane layout mbe_decodeImbe7200x4400Frame()
 *        reads. Cells past each vector's length are 0.
 *
 * All 88 bits round-trip through mbe_decodeImbe7200x4400Frame().
 *
 * @param imbe_d  Input parameter bits (88).
 * @param imbe_fr Output frame as 8x23 bitplanes.
 * @return 0 on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_encodeImbe7200x4400Frame(const char imbe_d[88], char imbe_fr[8][23]);

/* Prototypes from imbe7100x4400.c */
/** @brief Print IMBE 7100x4400 parameter bits to stderr (debug). */
MBE_API void mbe_dumpImbe7100x4400Data(const char* imbe_d);
/** @brief Print IMBE 7100x4400 frame to stderr (debug). */
MBE_API void mbe_dumpImbe7100x4400Frame(const char imbe_fr[7][24]);
/** @brief ECC correction for IMBE 7100x4400 C0. */
MBE_API int mbe_eccImbe7100x4400C0(char imbe_fr[7][24]);
/** @brief ECC and parameter packing for IMBE 7100x4400. */
MBE_API int mbe_eccImbe7100x4400Data(char imbe_fr[7][24], char* imbe_d);
/** @brief Demodulate IMBE 7100x4400 interleaved data. */
MBE_API int mbe_demodulateImbe7100x4400Data(char imbe[7][24]);
/** @brief Convert IMBE 7100x4400 parameter set into 7200x4400 layout. */
MBE_API int mbe_convertImbe7100to7200(char* imbe_d);
/**
 * @brief Decode a hard IMBE 7100x4400 frame to converted IMBE 4400 parameter bits without synthesis.
 * @param imbe_fr Input frame as 7x24 bitplanes; not modified.
 * @param imbe_d  Output parameter bits (88), converted to the 7200x4400/IMBE 4400 layout.
 * @param result  Optional output status; receives C0/protected/C4/total errors, valid-context flags and
 *                `MBE_PROCESS_FLAG_PROVOICE`.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeImbe7100x4400Frame(const char imbe_fr[7][24], char imbe_d[88], mbe_process_result* result);
/**
 * @brief Decode a soft IMBE 7100x4400 frame to converted IMBE 4400 parameter bits without synthesis.
 * @param imbe_fr Input soft frame as 7x24 bitplanes; not modified.
 * @param imbe_d  Output hard parameter bits (88), converted to the 7200x4400/IMBE 4400 layout.
 * @param result  Optional output status; also sets `MBE_PROCESS_FLAG_SOFT_INPUT`.
 * @return Corrected error total (`c0_errors + protected_errors`).
 */
MBE_API int mbe_decodeImbe7100x4400SoftFrame(const mbe_soft_bit imbe_fr[7][24], char imbe_d[88],
                                             mbe_process_result* result);
/**
 * @brief Process IMBE 7100x4400 frame into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param result   Optional output status populated by decode and synthesis.
 * @param imbe_fr  Input frame as 7x24 bitplanes.
 * @param imbe_d   Scratch/output parameter bits (88, converted to 7200 layout).
 * @param cur_mp,prev_mp,prev_mp_enhanced Parameter state as per Dataf variant.
 * @return Total error count on success, or a negative `MBE_STATUS_*` code.
 */
MBE_API int mbe_processImbe7100x4400Framef(float* aout_buf, mbe_process_result* result, const char imbe_fr[7][24],
                                           char imbe_d[88], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                           mbe_parms* prev_mp_enhanced);
/** @brief Process IMBE 7100x4400 frame into 16-bit PCM. */
MBE_API int mbe_processImbe7100x4400Frame(short* aout_buf, mbe_process_result* result, const char imbe_fr[7][24],
                                          char imbe_d[88], mbe_parms* cur_mp, mbe_parms* prev_mp,
                                          mbe_parms* prev_mp_enhanced);
/**
 * @brief Process a soft IMBE 7100x4400 frame into float PCM.
 * @param result Optional output status populated by decode and synthesis.
 * @return Final total error count used by the wrapper.
 */
MBE_API int mbe_processImbe7100x4400SoftFramef(float* aout_buf, mbe_process_result* result,
                                               const mbe_soft_bit imbe_fr[7][24], char imbe_d[88], mbe_parms* cur_mp,
                                               mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/** @brief Process a soft IMBE 7100x4400 frame into int16 PCM; see float variant for result semantics. */
MBE_API int mbe_processImbe7100x4400SoftFrame(short* aout_buf, mbe_process_result* result,
                                              const mbe_soft_bit imbe_fr[7][24], char imbe_d[88], mbe_parms* cur_mp,
                                              mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/**
 * @brief Get a pointer to a static NUL-terminated version string.
 *        The returned pointer remains valid for the lifetime of the program.
 */
MBE_API const char* mbe_versionString(void);

/**
 * @brief Set thread-local RNG seeds used by synthesis noise generators.
 *        Applies to comfort noise and unvoiced LCG cold-start state.
 * @param seed Any 32-bit seed value. A zero seed is accepted and remapped to
 *             an internal non-zero state.
 */
MBE_API void mbe_setThreadRngSeed(uint32_t seed);
/**
 * @brief Copy MBE parameter set from one struct to another.
 * @param source_mp Source parameters.
 * @param destination_mp Destination parameters. If either pointer is NULL, this is a no-op.
 */
MBE_API void mbe_moveMbeParms(const mbe_parms* source_mp, mbe_parms* destination_mp);
/**
 * @brief Replace current parameters with the last known parameters.
 * @param cur_mp Destination parameters to fill.
 * @param prev_mp Source parameters from previous frame. If either pointer is NULL, this is a no-op.
 */
MBE_API void mbe_useLastMbeParms(mbe_parms* cur_mp, const mbe_parms* prev_mp);
/**
 * @brief Initialize parameter state for decoding and synthesis.
 * @param cur_mp Output: current parameter state.
 * @param prev_mp Output: previous parameter state (zeroed/reset).
 * @param prev_mp_enhanced Output: enhanced previous parameter state. If any output pointer is NULL, this is a no-op.
 */
MBE_API void mbe_initMbeParms(mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);
/**
 * @brief Apply spectral amplitude enhancement in-place.
 *
 * Invalid parameter state, including an out-of-range harmonic count `L`, is
 * ignored.
 * @param cur_mp In/out parameter set to enhance.
 */
MBE_API void mbe_spectralAmpEnhance(mbe_parms* cur_mp);
/**
 * @brief Synthesize tone frame (AMBE tone indices) into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49).
 * @param cur_mp   Current parameter set (tone synthesis state). NULL or invalid tone inputs synthesize silence.
 */
MBE_API void mbe_synthesizeTonef(float* aout_buf, const char* ambe_d, mbe_parms* cur_mp);
/**
 * @brief Synthesize tone for D-STAR style indices into float PCM.
 * @param aout_buf Output buffer of 160 float samples.
 * @param ambe_d   AMBE parameter bits (49); the tone volume is read from them.
 * @param cur_mp   Current parameter set. NULL, invalid bits or an unknown index synthesize silence.
 * @param ID1      Tone index: single tones 5..122, DTMF 128..143 (128 + 4 * column + row),
 *                 call progress 144..147.
 */
MBE_API void mbe_synthesizeTonefdstar(float* aout_buf, const char* ambe_d, mbe_parms* cur_mp, int ID1);
/** @brief Fill float PCM buffer with 160 samples of silence. */
MBE_API void mbe_synthesizeSilencef(float* aout_buf);
/** @brief Fill 16-bit PCM buffer with 160 samples of silence. */
MBE_API void mbe_synthesizeSilence(short* aout_buf);
/**
 * @brief Synthesize one speech frame into float PCM.
 *
 * If `cur_mp` or `prev_mp` has an out-of-range harmonic count `L`, the output
 * buffer is filled with silence.
 * @param aout_buf Output buffer of 160 float samples.
 * @param cur_mp   Current parameter set.
 * @param prev_mp  Previous parameter set.
 */
MBE_API void mbe_synthesizeSpeechf(float* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp);
/**
 * @brief Synthesize one speech frame into 16-bit PCM.
 *
 * If `cur_mp` or `prev_mp` has an out-of-range harmonic count `L`, the output
 * buffer is filled with silence.
 * @param aout_buf Output buffer of 160 16-bit samples.
 * @param cur_mp   Current parameter set.
 * @param prev_mp  Previous parameter set.
 */
MBE_API void mbe_synthesizeSpeech(short* aout_buf, mbe_parms* cur_mp, mbe_parms* prev_mp);
/**
 * @brief Convert 160 float samples to clipped/scaled 16-bit PCM.
 *
 * Applies the same scaling used by the `short` entry points: a fixed gain of
 * `7.0` and soft clipping at 95% of int16 full-scale before converting to
 * `short`. Non-finite input samples are handled before conversion: NaN becomes
 * zero and infinities clip to the corresponding bound. This makes the output
 * equivalent to calling the corresponding `short` synthesis APIs with the same
 * input state.
 * @param float_buf Input 160 float samples.
 * @param aout_buf  Output 160 16-bit samples.
 */
MBE_API void mbe_floattoshort(const float* float_buf, short* aout_buf);

/* === Frame repeat and muting functions === */

/** Maximum consecutive frame repeats before muting. */
#define MBE_MAX_FRAME_REPEATS     4

/** IMBE muting threshold (8.75% error rate). */
#define MBE_MUTING_THRESHOLD_IMBE 0.0875f

/** AMBE muting threshold (9.6% error rate). */
#define MBE_MUTING_THRESHOLD_AMBE 0.096f

/**
 * @brief Check if frame should be muted due to excessive errors.
 * @param mp Parameter set to check.
 * @return Non-zero if frame should be muted.
 */
MBE_API int mbe_requiresMuting(const mbe_parms* mp);

/**
 * @brief Check if max repeat threshold has been exceeded.
 * @param mp Parameter set to check.
 * @return Non-zero if repeatCount >= MBE_MAX_FRAME_REPEATS.
 */
MBE_API int mbe_isMaxFrameRepeat(const mbe_parms* mp);

/**
 * @brief Generate comfort noise for muted frames.
 * @param aout_buf Output buffer of 160 float samples.
 */
MBE_API void mbe_synthesizeComfortNoisef(float* aout_buf);

/**
 * @brief Generate comfort noise for muted frames (16-bit).
 * @param aout_buf Output buffer of 160 16-bit samples.
 */
MBE_API void mbe_synthesizeComfortNoise(short* aout_buf);

/* === Adaptive smoothing functions === */

/**
 * @brief Apply adaptive smoothing to parameters based on error rates.
 *        Implements JMBE Algorithms #111-116.
 *
 * Invalid parameter state, including an out-of-range harmonic count `L` in
 * either parameter set, is ignored.
 * @param cur_mp Current frame parameters (modified in-place).
 * @param prev_mp Previous frame parameters (for local energy).
 */
MBE_API void mbe_applyAdaptiveSmoothing(mbe_parms* cur_mp, const mbe_parms* prev_mp);

/**
 * @brief Check if adaptive smoothing is required based on error rates.
 * @param mp Parameter set to check.
 * @return Non-zero if smoothing should be applied.
 */
MBE_API int mbe_requiresAdaptiveSmoothing(const mbe_parms* mp);

#ifdef __cplusplus
}
#endif

#endif // MBELIB_NEO_PUBLIC_MBEBELIB_H
