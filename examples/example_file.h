// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef MBE_EXAMPLE_FILE_H
#define MBE_EXAMPLE_FILE_H

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

#ifdef _WIN32
static inline HANDLE
example_open_identity_handle(const char* path) {
    HANDLE handle = CreateFileA(path, 0, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING,
                                FILE_FLAG_BACKUP_SEMANTICS, NULL);
    if (handle != INVALID_HANDLE_VALUE || GetLastError() != ERROR_FILENAME_EXCED_RANGE) {
        return handle;
    }
    /* Some ANSI APIs reject extended paths that the CRT can still open. */
    UINT code_page = AreFileApisANSI() != 0 ? CP_ACP : CP_OEMCP;
    int length = MultiByteToWideChar(code_page, 0, path, -1, NULL, 0);
    if (length == 0) {
        return INVALID_HANDLE_VALUE;
    }
    WCHAR* wide_path = (WCHAR*)malloc((size_t)length * sizeof(WCHAR));
    if (wide_path == NULL) {
        SetLastError(ERROR_NOT_ENOUGH_MEMORY);
        return INVALID_HANDLE_VALUE;
    }
    if (MultiByteToWideChar(code_page, 0, path, -1, wide_path, length) != 0) {
        handle = CreateFileW(wide_path, 0, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING,
                             FILE_FLAG_BACKUP_SEMANTICS, NULL);
    }
    DWORD error = GetLastError();
    free(wide_path);
    SetLastError(error);
    return handle;
}

static inline int
example_windows_paths_differ(const char* input, const char* output) {
    /* Windows CRT stat does not provide meaningful inode identities. */
    HANDLE in_handle = example_open_identity_handle(input);
    if (in_handle == INVALID_HANDLE_VALUE) {
        DWORD error = GetLastError();
        if (error == ERROR_FILE_NOT_FOUND || error == ERROR_PATH_NOT_FOUND) {
            return 1; /* Let example_open_file report the input error. */
        }
        (void)fprintf(stderr, "cannot inspect input file: Windows error %lu\n", (unsigned long)error);
        return 0;
    }
    HANDLE out_handle = example_open_identity_handle(output);
    if (out_handle == INVALID_HANDLE_VALUE) {
        DWORD error = GetLastError();
        (void)CloseHandle(in_handle);
        if (error == ERROR_FILE_NOT_FOUND || error == ERROR_PATH_NOT_FOUND) {
            return 1;
        }
        (void)fprintf(stderr, "cannot inspect output file: Windows error %lu\n", (unsigned long)error);
        return 0;
    }
    BY_HANDLE_FILE_INFORMATION in_info, out_info;
    int have_info =
        GetFileInformationByHandle(in_handle, &in_info) != 0 && GetFileInformationByHandle(out_handle, &out_info) != 0;
    DWORD error = GetLastError();
    (void)CloseHandle(in_handle);
    (void)CloseHandle(out_handle);
    if (!have_info) {
        (void)fprintf(stderr, "cannot inspect file identity: Windows error %lu\n", (unsigned long)error);
        return 0;
    }
    if (in_info.dwVolumeSerialNumber == out_info.dwVolumeSerialNumber
        && in_info.nFileIndexHigh == out_info.nFileIndexHigh && in_info.nFileIndexLow == out_info.nFileIndexLow) {
        (void)fprintf(stderr, "input and output must be different files\n");
        return 0;
    }
    return 1;
}
#endif

/* Check file identities before opening the output, including link aliases. */
static inline int
example_paths_differ(const char* input, const char* output) {
    if (strcmp(input, output) == 0) {
        (void)fprintf(stderr, "input and output must be different files\n");
        return 0;
    }
#ifdef _WIN32
    return example_windows_paths_differ(input, output);
#else
    struct stat in_stat, out_stat;
    if (stat(input, &in_stat) == 0 && stat(output, &out_stat) == 0 && in_stat.st_dev == out_stat.st_dev
        && in_stat.st_ino == out_stat.st_ino) {
        (void)fprintf(stderr, "input and output must be different files\n");
        return 0;
    }
    return 1;
#endif
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
