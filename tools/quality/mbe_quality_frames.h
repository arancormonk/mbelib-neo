// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_QUALITY_FRAMES_H
#define MBE_QUALITY_FRAMES_H

#include <stddef.h>

/* Developer fixtures, not a wire-conformance encoder or installed API.
 * Input is numeric 0/1 canonical parameter data (88 IMBE or 49 AMBE bits).
 * Success returns 184 (imbe7200), 168 (imbe7100), or 96 (ambe2450/ambe2400).
 * Invalid arguments/bits return the corresponding MBE_STATUS_* value without
 * modifying frame. Output is numeric 0/1 in rectangular row-major order.
 */
int mbe_quality_frame_from_data(const char* codec, const char* data, size_t count, char* frame, size_t capacity);

#endif /* MBE_QUALITY_FRAMES_H */
