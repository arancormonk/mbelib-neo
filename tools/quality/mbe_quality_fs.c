// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2026 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */
/** @file File identities and resolved paths for the quality tools. */
#include "mbe_quality_fs.h"

#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#if defined(_WIN32)
#include <io.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

FILE*
mbe_quality_open_output(const char* path) {
#if defined(_WIN32)
    int fd = _open(path, _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY, _S_IREAD | _S_IWRITE);
#else
    int fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
#endif
    if (fd < 0) {
        return NULL;
    }
#if defined(_WIN32)
    FILE* stream = _fdopen(fd, "wb");
#else
    FILE* stream = fdopen(fd, "wb");
#endif
    if (!stream) {
        int saved_errno = errno;
#if defined(_WIN32)
        _close(fd);
#else
        close(fd);
#endif
        errno = saved_errno;
    }
    return stream;
}

int
mbe_quality_same_file(const char* a, const char* b) {
#if defined(_WIN32)
    HANDLE ha = CreateFileA(a, 0, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING,
                            FILE_FLAG_BACKUP_SEMANTICS, NULL);
    if (ha == INVALID_HANDLE_VALUE) {
        return 0;
    }
    HANDLE hb = CreateFileA(b, 0, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING,
                            FILE_FLAG_BACKUP_SEMANTICS, NULL);
    if (hb == INVALID_HANDLE_VALUE) {
        CloseHandle(ha);
        return 0;
    }
    /* Volume serial + file index also catches hard-linked inputs. */
    BY_HANDLE_FILE_INFORMATION ia, ib;
    int same = GetFileInformationByHandle(ha, &ia) && GetFileInformationByHandle(hb, &ib)
               && ia.dwVolumeSerialNumber == ib.dwVolumeSerialNumber && ia.nFileIndexHigh == ib.nFileIndexHigh
               && ia.nFileIndexLow == ib.nFileIndexLow;
    CloseHandle(ha);
    CloseHandle(hb);
    return same;
#else
    struct stat sa, sb;
    return stat(a, &sa) == 0 && stat(b, &sb) == 0 && sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
#endif
}

char*
mbe_quality_normalized_path(const char* path) {
    if (!path || !*path) {
        errno = EINVAL;
        return NULL;
    }
#if defined(_WIN32)
    /* Reject UNC, extended and device namespaces before Win32 rewrites them.
     * Forward slashes are separators too. */
    if ((path[0] == '\\' || path[0] == '/') && (path[1] == '\\' || path[1] == '/')) {
        errno = EINVAL;
        return NULL;
    }
    /* The ANSI Win32 APIs and _stricmp do not fold case the way NTFS does for
     * non-ASCII names, so absent outputs differing only in such letters could
     * alias. Keep the supported alphabet ASCII. */
    for (const unsigned char* p = (const unsigned char*)path; *p; ++p) {
        if (*p >= 0x80u || *p < 0x20u) {
            errno = EINVAL;
            return NULL;
        }
    }
    /* A ':' anywhere after the drive letter names an alternate data stream or
     * a device; neither is an ordinary file the alias checks can reason about. */
    const char* drive_end =
        (((path[0] >= 'A' && path[0] <= 'Z') || (path[0] >= 'a' && path[0] <= 'z')) && path[1] == ':') ? path + 2
                                                                                                       : path;
    if (strchr(drive_end, ':')) {
        errno = EINVAL;
        return NULL;
    }
    const char* input_leaf = path;
    for (const char* p = path; *p; ++p) {
        if (*p == '\\' || *p == '/' || *p == ':') {
            input_leaf = p + 1;
        }
    }
    size_t input_length = strlen(input_leaf);
    if (input_length && strcmp(input_leaf, ".") != 0 && strcmp(input_leaf, "..") != 0
        && (input_leaf[input_length - 1] == '.' || input_leaf[input_length - 1] == ' ')) {
        errno = EINVAL;
        return NULL;
    }
    /* Root the path before splitting it: C:foo uses that drive's current
     * directory, while \foo uses the current drive's root. */
    DWORD needed = GetFullPathNameA(path, 0, NULL, NULL);
    if (needed == 0) {
        errno = EIO;
        return NULL;
    }
    char* full = malloc(needed);
    if (!full) {
        errno = ENOMEM;
        return NULL;
    }
    DWORD written = GetFullPathNameA(path, needed, full, NULL);
    if (written == 0 || written >= needed) {
        free(full);
        errno = EIO;
        return NULL;
    }
    /* Reserved device names such as NUL or CON come back as \\.\ spellings. */
    if (full[0] == '\\' && full[1] == '\\') {
        free(full);
        errno = EINVAL;
        return NULL;
    }
    HANDLE handle = CreateFileA(full, 0, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING,
                                FILE_FLAG_BACKUP_SEMANTICS, NULL);
    char* leaf = NULL;
    if (handle == INVALID_HANDLE_VALUE) {
        DWORD error = GetLastError();
        DWORD attributes = GetFileAttributesA(full);
        if (attributes != INVALID_FILE_ATTRIBUTES && (attributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0u) {
            /* An existing reparse point that cannot be opened is dangling. */
            free(full);
            errno = ELOOP;
            return NULL;
        }
        if (error != ERROR_FILE_NOT_FOUND && error != ERROR_PATH_NOT_FOUND) {
            free(full);
            errno = EIO;
            return NULL;
        }
        char* slash = strrchr(full, '\\');
        if (!slash || !slash[1]) {
            free(full);
            errno = EINVAL;
            return NULL;
        }
        size_t leaf_size = strlen(slash + 1) + 1;
        leaf = malloc(leaf_size);
        if (!leaf) {
            free(full);
            errno = ENOMEM;
            return NULL;
        }
        memcpy(leaf, slash + 1, leaf_size);
        /* Keep the separator, especially for a root parent such as C:\. */
        slash[1] = '\0';
        handle = CreateFileA(full, 0, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING,
                             FILE_FLAG_BACKUP_SEMANTICS, NULL);
        BY_HANDLE_FILE_INFORMATION info;
        if (handle == INVALID_HANDLE_VALUE || !GetFileInformationByHandle(handle, &info)
            || (info.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0u) {
            if (handle != INVALID_HANDLE_VALUE) {
                CloseHandle(handle);
            }
            free(leaf);
            free(full);
            errno = EIO;
            return NULL;
        }
    }
    free(full);
    needed = GetFinalPathNameByHandleA(handle, NULL, 0, FILE_NAME_NORMALIZED | VOLUME_NAME_DOS);
    if (needed == 0) {
        CloseHandle(handle);
        free(leaf);
        errno = EIO;
        return NULL;
    }
    /* Windows reports the required size including the terminator; Wine reports
     * it without. Allow one extra byte so both resolve the path. */
    char* parent = malloc((size_t)needed + 1u);
    if (!parent) {
        CloseHandle(handle);
        free(leaf);
        errno = ENOMEM;
        return NULL;
    }
    written = GetFinalPathNameByHandleA(handle, parent, needed + 1u, FILE_NAME_NORMALIZED | VOLUME_NAME_DOS);
    CloseHandle(handle);
    if (written == 0 || written > needed) {
        free(parent);
        free(leaf);
        errno = EIO;
        return NULL;
    }
    /* A junction can resolve onto a network share even from a drive path. */
    if (_strnicmp(parent, "\\\\?\\UNC\\", 8) == 0) {
        free(parent);
        free(leaf);
        errno = EINVAL;
        return NULL;
    }
    /* Drop the extended prefix returned by the handle query for DOS paths. */
    if (strncmp(parent, "\\\\?\\", 4) == 0) {
        memmove(parent, parent + 4, strlen(parent + 4) + 1);
    }
    if (!leaf) {
        return parent;
    }
    const char separator = '\\';
#else
    char* existing = realpath(path, NULL);
    if (existing) {
        return existing;
    }
    if (errno != ENOENT) {
        return NULL;
    }
    struct stat unresolved;
    if (lstat(path, &unresolved) == 0) {
        errno = ELOOP; /* A dangling leaf must not be treated as a new output. */
        return NULL;
    }
    if (errno != ENOENT) {
        return NULL;
    }
    char* copy = strdup(path);
    if (!copy) {
        return NULL;
    }
    char* slash = strrchr(copy, '/');
    const char* leaf = slash ? slash + 1 : copy;
    if (!*leaf) {
        free(copy);
        errno = EINVAL;
        return NULL;
    }
    if (slash) {
        *slash = '\0';
    }
    char* parent = realpath(slash ? (*copy ? copy : "/") : ".", NULL);
    if (!parent) {
        int error = errno;
        free(copy);
        errno = error;
        return NULL;
    }
    const char separator = '/';
#endif
    size_t parent_length = strlen(parent);
    size_t leaf_length = strlen(leaf);
    char* resolved = NULL;
    int overflow = leaf_length > SIZE_MAX - 2 || parent_length > SIZE_MAX - 2 - leaf_length;
    if (!overflow) {
        size_t size = parent_length + leaf_length + 2;
        resolved = malloc(size);
        if (resolved) {
            /* A root already ends in a separator. Avoid two equivalent but
             * lexically different spellings when comparing absent paths. */
            char join[2] = {separator, '\0'};
            if (parent_length && parent[parent_length - 1] == separator) {
                join[0] = '\0';
            }
            snprintf(resolved, size, "%s%s%s", parent, join, leaf);
        }
    }
    free(parent);
#if defined(_WIN32)
    free(leaf);
#else
    free(copy);
#endif
    if (!resolved) {
        errno = overflow ? EOVERFLOW : ENOMEM;
    }
    return resolved;
}

int
mbe_quality_paths_equal(const char* a, const char* b) {
#if defined(_WIN32)
    return _stricmp(a, b) == 0;
#else
    return strcmp(a, b) == 0;
#endif
}
