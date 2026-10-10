// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Tone detection for the AMBE+2 and D-STAR encoders' tone frames (private API).
 *
 * Recognizes the tones TIA-102.BABA-1 Table 9 can carry in one 20 ms span:
 * single tones (index 5..122, f = 31.25 Hz * index), DTMF, KNOX and
 * call-progress dual tones.
 */
#ifndef MBELIB_NEO_INTERNAL_MBE_TONE_DETECT_H
#define MBELIB_NEO_INTERNAL_MBE_TONE_DETECT_H

#include "mbe_unvoiced_fft.h"

#define MBE_TONE_SPAN   160

/* The encoders centre the span this many samples after the voice analysis,
 * which aligns their decisions at a tone's start and end with DVSI's. DVSI's
 * D-STAR and rate-33 encoders make the same decisions frame for frame. */
#define MBE_TONE_OFFSET 16

struct mbe_tone_detection {
    int id;          /* Table 9 tone index */
    float amplitude; /* peak amplitude per component (16-bit scale); a pair gives the geometric mean */
};

/*
 * Look for a tone in span (MBE_TONE_SPAN DC-filtered samples on the 16-bit
 * scale). Returns 1 with *out filled, 0 when the span is not a supported
 * tone, or a negative MBE_STATUS_* value from the FFT.
 */
int mbe_tone_detect(mbe_fft_plan* fft, const float span[MBE_TONE_SPAN], struct mbe_tone_detection* out);

#endif /* MBELIB_NEO_INTERNAL_MBE_TONE_DETECT_H */
