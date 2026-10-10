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

/*
 * Tone timing, approximating DVSI's encoders: a call-progress tone is
 * sent once it has filled MBE_TONE_CP_CONFIRM consecutive spans (60 ms, which
 * an 80 ms tone always does; DVSI sends none for 40-45 ms bursts), then while
 * it is detected and for up to MBE_TONE_CP_HOLD frames after, as its latest
 * detection. DVSI starts one frame earlier, which needs a detection in the
 * partly filled first span. A DTMF, KNOX or single tone is sent from its first
 * detection and ends that hold, and is not held itself: DVSI's encoders send
 * one for the frames whose span it fills at least 45%, as the detector does.
 * Where one changes directly to another, the span holding both carries the
 * newer if it is a DTMF or KNOX tone that fills the newest half of the span,
 * else the older if it fills the oldest half, as DVSI's encoders send them.
 */
#define MBE_TONE_CP_CONFIRM 3
#define MBE_TONE_CP_HOLD    2

/** Per-stream tone timing state. */
struct mbe_tone_tracker {
    int run_id;                     /* call-progress tone detected in the last frame, or -1 */
    int run;                        /* consecutive frames it has been detected in, at most MBE_TONE_CP_CONFIRM */
    int hold;                       /* frames last may still be sent without a detection */
    int sending;                    /* last is being sent */
    struct mbe_tone_detection last; /* latest detection of the call-progress tone sent */
    int tone_id;                    /* DTMF, KNOX or single tone returned for the last frame, or -1 */
};

void mbe_tone_tracker_reset(struct mbe_tone_tracker* tracker);

/*
 * Detect a tone in span (as mbe_tone_detect()) and decide what this frame
 * sends. Returns 1 with *out filled, 0 for voice, or a negative MBE_STATUS_*
 * value from the FFT, which leaves the tracker unchanged.
 */
int mbe_tone_track(struct mbe_tone_tracker* tracker, mbe_fft_plan* fft, const float span[MBE_TONE_SPAN],
                   struct mbe_tone_detection* out);

#endif /* MBELIB_NEO_INTERNAL_MBE_TONE_DETECT_H */
