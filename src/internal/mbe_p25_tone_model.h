// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */
/**
 * @file
 * @brief Whether a model is a P25 tone model as the decoder plays one (private API).
 *
 * P25's rendering (mbe_subframe.c) reads a tone model's phases as decoded and
 * its voicing as it stands, where it reads every other voiced model's phases
 * with the unvoiced factor taken out of its unvoiced bands and its voicing by
 * the 500 Hz columns (docs/dvsi-parity.md, "P25 tones sent as voice frames").
 * The previous model it reads is caller state, so the synthesizer asks the
 * codec, whose decode plays the tone model (imbe4400_play_tone() in
 * imbe7200x4400.c), whether a model is one: the answer is checked there with
 * the decoder's own functions (its tone rows and its pitch function, the one
 * body, not inlined, the decode stores the pitch with), as mbe_pulsed.h's code
 * tests are.
 */
#ifndef MBEINT_MBE_P25_TONE_MODEL_H
#define MBEINT_MBE_P25_TONE_MODEL_H

#include "mbelib-neo/mbelib.h"

/**
 * @brief Whether @p mp is a P25 tone model as imbe4400_play_tone() leaves it.
 *
 * Marked MBE_FRAME_MARK_TONE; Vl[0] a tone word (its flag, a pitch code b0
 * with tone rows, and harmonics l1 <= l2 that are those of a row of that
 * code); the model at the pitch (by bit pattern), harmonic count and band
 * count that code decodes to; exactly harmonics l1 and l2 voiced (Vl 1), every
 * other harmonic 1..L unvoiced (Vl 0); the tone harmonics' magnitudes
 * positive, finite and normal, a dual tone's two in the ratio of their droops
 * sinc^2(pi l / 256) the player writes (within 2^-20 of it), and every other
 * magnitude 0 by its bit pattern; log2Ml[1..L] finite. Not constrained: the
 * tone's common level, the finite values of log2Ml, Vl and Ml above L,
 * gamma, the pulse train and the rest of the synthesis state, which are state
 * like any other (changing them changes the next frame). Every model the
 * player marks passes (imbe4400_take_tone() takes a tone only where the two
 * levels it will write are positive normal floats).
 * @param mp A model; any state (its harmonic count is checked).
 * @return 1 for a tone model, 0 otherwise.
 */
int mbe_imbe4400_tone_model(const mbe_parms* mp);

#endif /* MBEINT_MBE_P25_TONE_MODEL_H */
