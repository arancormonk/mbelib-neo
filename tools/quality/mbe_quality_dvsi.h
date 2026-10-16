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

/* DVSI 4-bit soft-decision vectors (*_sd.bit, written with -sdbits 4): one
 * value 0..15 per bit of the corresponding hard-decision file, two per byte,
 * high nibble first, in that file's MSB-first bit order; 8..15 decide 1, and
 * the distance from 7.5 is the confidence. Places each value in the cell its
 * hard bit occupies (4 * frame bytes in, the hard import's cell count out);
 * cells without a channel bit read 0. Returns the cell count or a negative
 * MBE_STATUS_* value; frame is written only on success.
 */
int mbe_quality_soft_frame_from_dvsi(const char* codec, const unsigned char* bytes, size_t count, unsigned char* frame,
                                     size_t capacity);

#endif /* MBE_QUALITY_DVSI_H */
