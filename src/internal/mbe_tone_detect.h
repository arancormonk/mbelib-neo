// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 Rhizomatica
 * Author: Rafael Diniz <rafael@rhizomatica.org>
 *
 * Based on Bruce Perens' ham_digital_modes
 * (https://github.com/BrucePerens/hams_open).
 */

/**
 * @file
 * @brief Encoder-side tone detection (private API).
 *
 * C port of ham_digital_modes' float tone_detect.rs: recognizes a DTMF digit
 * or a single sustained tone in one 20 ms frame, matched by its author to the
 * AMBE-3000 chip's own detector, so an encoder can send a tone frame instead
 * of coding the tone as speech.
 */
#ifndef MBELIB_NEO_INTERNAL_MBE_TONE_DETECT_H
#define MBELIB_NEO_INTERNAL_MBE_TONE_DETECT_H

#define MBE_TONE_DETECT_FRAME 160

enum mbe_tone_detect_kind {
    MBE_TONE_DETECT_NONE = 0,
    MBE_TONE_DETECT_DTMF,   /* row 0-3 (697..941 Hz), column 0-3 (1209..1633 Hz) */
    MBE_TONE_DETECT_SINGLE, /* index = round(f / 31.25 Hz) */
};

struct mbe_tone_detection {
    enum mbe_tone_detect_kind kind;
    int row;
    int col;
    int index;
    double hz;
    double amplitude; /* per-tone peak amplitude, 16-bit scale */
};

/* Detect a tone in one frame on the 16-bit scale; returns the detection kind. */
enum mbe_tone_detect_kind mbe_tone_detect(const double frame[MBE_TONE_DETECT_FRAME], struct mbe_tone_detection* out);

#endif /* MBELIB_NEO_INTERNAL_MBE_TONE_DETECT_H */
