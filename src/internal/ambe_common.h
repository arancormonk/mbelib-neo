// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Internal helpers for AMBE 3600x{2400,2450} ECC and demodulation.
 *
 * Declares common routines used by both AMBE 3600x2400 and 3600x2450
 * implementations to correct C0 with Golay24-compatible behavior
 * (Golay(23,12) + parity bit handling), demodulate C1, and extract the
 * 49-bit parameter vector.
 */

#ifndef MBELIB_NEO_INTERNAL_AMBE_COMMON_H
#define MBELIB_NEO_INTERNAL_AMBE_COMMON_H

#include "mbelib-neo/mbelib.h"

/**
 * @brief Correct C0 for AMBE 3600x{2400,2450} with Golay24 parity behavior.
 *
 * Applies Golay(23,12) decoding to `fr[0][1..23]` and then applies
 * Golay24 parity-bit correction on `fr[0][0]` when the protected
 * 23-bit codeword has zero syndrome.
 *
 * @param fr AMBE frame as 4x24 bitplanes (modified).
 * @return Number of corrected bit errors in C0.
 */
int mbe_eccAmbe3600C0_common(char fr[4][24]);
int mbe_eccAmbe3600C0Soft_common(mbe_soft_bit fr[4][24]);

/**
 * @brief Demodulate AMBE 3600x{2400,2450} C1 in-place.
 *
 * Uses a pseudo-random sequence derived from the C0 payload to remove
 * the interleaving/modulation applied to C1.
 *
 * @param fr AMBE frame as 4x24 bitplanes (modified).
 */
void mbe_demodulateAmbe3600Data_common(char fr[4][24]);
void mbe_demodulateAmbe3600DataSoft_common(mbe_soft_bit fr[4][24]);

/**
 * @brief Extract 49 parameter bits from C0..C3 with ECC.
 *
 * Copies C0, demodulates C1 and applies Golay(23,12), and copies C2/C3
 * into the 49-bit output parameter vector.
 *
 * @param fr     AMBE frame as 4x24 bitplanes (modified by demodulation).
 * @param out49  Output parameter bits (49 entries).
 * @return Number of corrected bit errors in protected fields.
 */
int mbe_eccAmbe3600Data_common(char fr[4][24], char* out49);
int mbe_eccAmbe3600DataSoft_common(mbe_soft_bit fr[4][24], char* out49);

/** AMBE+2 silence-frame fundamental 2*pi/32, TIA-102.BABA-1 4.1 eq 1 (radians/sample). */
#define MBE_AMBE_SILENCE_W0 ((float)0.19634954084936207)
/** AMBE+2 silence-frame harmonic count, TIA-102.BABA-1 4.1 eq 2. */
#define MBE_AMBE_SILENCE_L  14

/*
 * Initial AMBE model. TIA-102.BABA-1 fixes only the initial prediction state:
 * gamma(-1) = 0, L(-1) = 15 and a constant log magnitude (a constant cancels
 * exactly in eq 43, so 0 is equivalent to the spec's 1). The initial
 * fundamental is unspecified; the silence fundamental keeps every initial
 * harmonic below Nyquist (15 * 2*pi/32 < pi). JMBE's W124 default used
 * w0 = (pi/32) * 2*pi, which put harmonics 6..15 above Nyquist.
 */
#define MBE_AMBE_INIT_W0    MBE_AMBE_SILENCE_W0
#define MBE_AMBE_INIT_L     15

/**
 * @brief Set the initial AMBE model: w0, L, K = 0, gamma = 0, unit unvoiced amplitudes.
 *
 * Shared by the AMBE decode paths and the AMBE 2400 encoder, whose silence
 * reset mirrors the decoder's state.
 *
 * @param mp Parameter set to rewrite (phase, error and synthesis state untouched).
 */
void mbe_setAmbeDefaultModel_common(mbe_parms* mp);

/**
 * @brief Initialize AMBE parameter state to the initial AMBE model.
 *
 * Sets the model from mbe_setAmbeDefaultModel_common() with zero phases,
 * default smoothing and error state, the AMBE muting threshold and a cold
 * unvoiced-noise start, and copies it to all three parameter sets.
 *
 * @param cur_mp  Output current parameter state.
 * @param prev_mp Output previous parameter state.
 * @param prev_mp_enhanced Output enhanced previous parameter state.
 */
void mbe_initAmbeParms_common(mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

/**
 * @brief Set AMBE parameters to JMBE-style ERASURE defaults (W120).
 *
 * Mirrors AMBEModelParameters.setDefaults(FrameType.ERASURE): w0=0,
 * L=9, unvoiced bands, unit amplitudes, zero gain term. Keeps phase and
 * unvoiced overlap/noise state from `state_src` to preserve synthesizer
 * continuity behavior.
 *
 * @param mp Target parameter set to rewrite.
 * @param state_src Source state for phase/noise continuity (may be NULL).
 */
void mbe_setAmbeErasureParms_common(mbe_parms* mp, const mbe_parms* state_src);

/**
 * @brief Ensure AMBE parameter state is initialized with AMBE defaults.
 *
 * If the previous state does not appear to be AMBE-initialized, this
 * routine resets all three state structs via mbe_initAmbeParms_common().
 *
 * @param cur_mp  In/out current parameter state.
 * @param prev_mp In/out previous parameter state.
 * @param prev_mp_enhanced In/out enhanced previous parameter state.
 */
void mbe_ensureAmbeDefaults_common(mbe_parms* cur_mp, mbe_parms* prev_mp, mbe_parms* prev_mp_enhanced);

#endif /* MBELIB_NEO_INTERNAL_AMBE_COMMON_H */
