// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_QUALITY_DVSI_H
#define MBE_QUALITY_DVSI_H

#include <stddef.h>

/* Developer fixtures for DVSI hard-decision test vectors (local use only).
 * imbe7200: 18-byte P25 frames, MSB first, on-air interleave.
 * ambe2450: 9-byte rate-33 frames, MSB first, DMR interleave.
 * ambe2400: 9-byte D-STAR DV data, LSB first.
 * Returns the frame byte count for a codec, or MBE_STATUS_INVALID_ARGUMENT.
 */
int mbe_quality_dvsi_frame_bytes(const char* codec);

/* Deinterleave one vector frame into a rectangular row-major 0/1 frame
 * (184 bits for imbe7200, 96 for AMBE), keeping channel bits unchanged.
 * Returns the bit count or a negative MBE_STATUS_* value; frame is written
 * only on success.
 */
int mbe_quality_frame_from_dvsi(const char* codec, const unsigned char* bytes, size_t count, char* frame,
                                size_t capacity);

#endif /* MBE_QUALITY_DVSI_H */
