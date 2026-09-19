// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_EXAMPLE_FILE_H
#define MBE_EXAMPLE_FILE_H

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

/* Check before opening the output, including existing symlink/hardlink aliases. */
static inline int
example_paths_differ(const char* input, const char* output) {
    struct stat in_stat, out_stat;
    if (strcmp(input, output) == 0
        || (stat(input, &in_stat) == 0 && stat(output, &out_stat) == 0 && in_stat.st_dev == out_stat.st_dev
            && in_stat.st_ino == out_stat.st_ino)) {
        (void)fprintf(stderr, "input and output must be different files\n");
        return 0;
    }
    return 1;
}

static inline FILE*
example_open_file(const char* path, const char* mode) {
    if (path == NULL || path[0] == '\0' || mode == NULL || (strcmp(mode, "rb") != 0 && strcmp(mode, "wb") != 0)) {
        errno = EINVAL;
        (void)fprintf(stderr, "invalid file path or mode: %s\n", strerror(errno));
        return NULL;
    }
    // nosemgrep: mbelib-neo.no-raw-file-open -- Approved example helper validates arguments and reports open errors.
    FILE* fp = fopen(path, mode);
    if (fp == NULL) {
        (void)fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    }
    return fp;
}

static inline int
example_close_file(FILE* fp) {
    if (fp != NULL && fclose(fp) != 0) {
        (void)fprintf(stderr, "close error: %s\n", strerror(errno));
        return -1;
    }
    return 0;
}

#endif
