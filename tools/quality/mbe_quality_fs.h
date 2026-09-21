// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Private filesystem checks for the quality tools. */
#ifndef MBE_QUALITY_FS_H
#define MBE_QUALITY_FS_H

#include <stdio.h>

/* Open a binary output, truncating an existing file. New files are created with
 * owner read/write permissions (0600 on POSIX; Windows also applies its ACL).
 * Existing permissions are retained. NULL indicates failure and sets errno. */
FILE* mbe_quality_open_output(const char* path);

/* Return 1 only if both paths open/stat successfully and identify the same file.
 * Missing or inaccessible paths return 0; this is not a path-validation check. */
int mbe_quality_same_file(const char* a, const char* b);

/* Resolve an existing path or an absent leaf under an existing directory.
 * The caller frees the returned string. NULL indicates failure and sets errno:
 * EINVAL for unsupported spellings, ELOOP for a dangling symlink, ENOMEM for
 * allocation failure, EOVERFLOW for an unrepresentable length, EIO when the
 * platform cannot resolve the path, or the underlying filesystem error.
 * Win32 accepts ordinary ASCII drive paths only: non-ASCII or control bytes,
 * UNC, device and extended namespaces, alternate data streams (a ':' after the
 * drive letter), and leaves ending in a dot or space are rejected (the
 * navigation components "." and ".." are fine).
 * A volume without a drive letter resolves to a Volume{GUID} spelling that is
 * only meaningful for comparison. These checks do not protect against
 * concurrent replacement of paths. */
char* mbe_quality_normalized_path(const char* path);

/* Compare path strings, case-insensitively on Windows. Normalize first when
 * comparing different spellings; use same_file as well to detect hard links. */
int mbe_quality_paths_equal(const char* a, const char* b);

#endif /* MBE_QUALITY_FS_H */
