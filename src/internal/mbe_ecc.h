// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 */

/**
 * @file
 * @brief Internal ECC helpers shared by decode and encode paths.
 *
 * The Golay/Hamming decoders are exposed publicly; the corresponding
 * encoders are only needed inside the library (e.g. by the AMBE
 * 3600x2400 D-STAR encoder) and live here.
 */

#ifndef MBEINT_MBE_ECC_H
#define MBEINT_MBE_ECC_H

/**
 * @brief Encode a 12-bit data word into a (23,12) Golay codeword.
 *
 * Data bits occupy codeword positions 11..22 (MSB at 22), parity bits
 * positions 0..10. This is the exact inverse of mbe_golay2312().
 *
 * @param in12  Input data bits (12), in[0] = MSB.
 * @param out23 Output codeword bits (23).
 */
void mbe_golay2312_encode(const char in12[12], char out23[23]);

/**
 * @brief Encode 11 data bits into a (15,11) Hamming codeword with the IMBE
 *        7200x4400 generator.
 *
 * The data occupy codeword positions 14..4 (in11[0] at 14) and the parity
 * bits positions 3..0, so mbe_hamming1511() returns in11 unchanged.
 *
 * @param in11  Input data bits (11).
 * @param out15 Output codeword bits (15).
 */
void mbe_hamming1511_encode(const char in11[11], char out15[15]);

#endif /* MBEINT_MBE_ECC_H */
