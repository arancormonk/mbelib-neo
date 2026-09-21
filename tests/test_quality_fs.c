// SPDX-License-Identifier: GPL-2.0-or-later
/** @file Regression tests for quality-tool file identities and path aliases. */
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "mbe_quality_fs.h"

static void
check(int condition, const char* message) {
    if (!condition) {
        fprintf(stderr, "%s (errno=%d)\n", message, errno);
        exit(EXIT_FAILURE);
    }
}

static void
make_file(const char* path) {
    FILE* file = fopen(path, "wb");
    check(file != NULL, "cannot create test file");
    check(fputs("quality filesystem test\n", file) >= 0, "cannot write test file");
    check(fclose(file) == 0, "cannot close test file");
}

static void
assert_normalized_equal(const char* a, const char* b) {
    char* na = mbe_quality_normalized_path(a);
    char* nb = mbe_quality_normalized_path(b);
    check(na != NULL && nb != NULL, "cannot normalize test paths");
    assert(mbe_quality_paths_equal(na, nb));
    free(na);
    free(nb);
}

static void
assert_rejected(const char* path, int expected_errno) {
    errno = 0;
    char* normalized = mbe_quality_normalized_path(path);
    assert(normalized == NULL);
    assert(errno == expected_errno);
    free(normalized);
}

int
main(void) {
#if defined(_WIN32)
    char original[MAX_PATH];
    DWORD length = GetCurrentDirectoryA(sizeof(original), original);
    check(length > 0 && length < sizeof(original), "cannot save working directory");
    char temp[MAX_PATH];
    length = GetTempPathA(sizeof(temp), temp);
    check(length > 0 && length < sizeof(temp), "cannot locate temporary directory");
    char directory[MAX_PATH];
    /* Reserve a unique name, then replace the file with our test directory. */
    check(GetTempFileNameA(temp, "mbq", 0, directory) != 0, "cannot reserve temporary name");
    check(DeleteFileA(directory), "cannot remove temporary placeholder");
    check(CreateDirectoryA(directory, NULL), "cannot create temporary directory");
    check(SetCurrentDirectoryA(directory), "cannot enter temporary directory");
    check(CreateDirectoryA("child", NULL), "cannot create child directory");
#else
    char original[4096];
    check(getcwd(original, sizeof(original)) != NULL, "cannot save working directory");
    const char* tmpdir = getenv("TMPDIR");
    if (!tmpdir || !*tmpdir) {
        tmpdir = "/tmp";
    }
    char directory[4096];
    check(snprintf(directory, sizeof(directory), "%s/mbe-quality-fs-XXXXXX", tmpdir) < (int)sizeof(directory),
          "temporary directory path is too long");
    check(mkdtemp(directory) != NULL, "cannot create temporary directory");
    check(chdir(directory) == 0, "cannot enter temporary directory");
    check(mkdir("child", 0700) == 0, "cannot create child directory");
#endif
    make_file("first");
    make_file("second");
    assert(mbe_quality_same_file("first", "first"));
    assert(mbe_quality_paths_equal("first", "first"));
    assert(!mbe_quality_same_file("first", "second"));
    assert(!mbe_quality_paths_equal("first", "second"));
    assert(!mbe_quality_same_file("first", "absent"));

#if defined(_WIN32)
    int linked = CreateHardLinkA("hardlink", "first", NULL) != 0;
#else
    int linked = link("first", "hardlink") == 0;
#endif
    if (linked) {
        assert(mbe_quality_same_file("first", "hardlink"));
        check(remove("hardlink") == 0, "cannot remove hard link");
    } else {
        puts("SKIP: hard-link creation unavailable");
    }

    char* parent = mbe_quality_normalized_path(".");
    char* absent = mbe_quality_normalized_path("absent");
    check(parent != NULL && absent != NULL, "cannot normalize absent output");
    size_t size = strlen(parent) + sizeof("/absent");
    char* expected = malloc(size);
    check(expected != NULL, "cannot allocate expected path");
#if defined(_WIN32)
    snprintf(expected, size, "%s\\absent", parent);
#else
    snprintf(expected, size, "%s/absent", parent);
#endif
    assert(mbe_quality_paths_equal(expected, absent));
    free(expected);
    free(absent);
    assert_normalized_equal("first", "./child/../first");
    assert_normalized_equal("absent", "./child/../absent");
    assert_normalized_equal(".", "child/..");
    assert_normalized_equal("child", "./child/.");
    assert_rejected("", EINVAL);
#if defined(_WIN32)
    assert_rejected("missing-parent/absent", EIO);
    assert_rejected("first/absent", EIO);
#else
    assert_rejected("missing-parent/absent", ENOENT);
    assert_rejected("first/absent", ENOTDIR);
#endif

#if defined(_WIN32)
    assert(mbe_quality_paths_equal("absent", "ABSENT"));
    assert_normalized_equal("absent", "ABSENT");
    assert_normalized_equal("first", "FIRST");
    check(strlen(parent) >= 3 && parent[1] == ':' && parent[2] == '\\', "expected a rooted drive path");
    char drive_relative[] = "C:absent";
    drive_relative[0] = parent[0];
    assert_normalized_equal(drive_relative, "absent");
    /* Repeat from a subdirectory so the spelling follows the current directory
     * rather than becoming a literal C:absent leaf. */
    check(SetCurrentDirectoryA("child"), "cannot enter child directory");
    assert_normalized_equal(drive_relative, "absent");
    check(SetCurrentDirectoryA(".."), "cannot return to temporary directory");

    /* Derive a root leaf from the unique temporary name so it cannot exist. */
    const char* temp_leaf = strrchr(parent, '\\');
    check(temp_leaf != NULL && temp_leaf[1] != '\0', "cannot derive a unique root leaf");
    temp_leaf++;
    size = sizeof("C:\\") + strlen(temp_leaf) + sizeof("-absent");
    char* root_absolute = malloc(size);
    char* root_rooted = malloc(size);
    check(root_absolute != NULL && root_rooted != NULL, "cannot allocate root spellings");
    snprintf(root_absolute, size, "%c:\\%s-absent", parent[0], temp_leaf);
    snprintf(root_rooted, size, "\\%s-absent", temp_leaf);
    assert_normalized_equal(root_absolute, root_rooted);
    /* Walk to the drive root from this directory using only relative steps. */
    size = strlen(parent) * 3 + strlen(temp_leaf) + sizeof("-absent");
    char* root_relative = calloc(size, 1);
    check(root_relative != NULL, "cannot allocate root-relative spelling");
    size_t offset = 0;
    for (const char* p = parent + 2; *p; ++p) {
        if (*p == '\\' && p[1]) {
            memcpy(root_relative + offset, "..\\", 3);
            offset += 3;
        }
    }
    snprintf(root_relative + offset, size - offset, "%s-absent", temp_leaf);
    assert_normalized_equal(root_absolute, root_relative);
    free(root_relative);
    free(root_rooted);
    free(root_absolute);

    assert_rejected("absent.", EINVAL);
    assert_rejected("absent ", EINVAL);
    assert_rejected("first.", EINVAL);
    assert_rejected("first ", EINVAL);
    assert_rejected("\\\\?\\C:\\absent", EINVAL);
    assert_rejected("\\\\.\\C:\\absent", EINVAL);
    assert_rejected("\\\\server\\share\\absent", EINVAL);
    assert_rejected("//server/share/absent", EINVAL);
    assert_rejected("first:stream", EINVAL);
    assert_rejected("nul", EINVAL);
    assert_rejected("CON", EINVAL);
    assert_rejected("caf\xe9.wav", EINVAL);
    assert_rejected("tab\tname", EINVAL);

    /* Symlink creation may need Developer Mode or elevated privileges. */
    if (CreateSymbolicLinkA("symlink", "first", 0)) {
        assert(mbe_quality_same_file("first", "symlink"));
        assert_normalized_equal("first", "symlink");
        check(DeleteFileA("symlink"), "cannot remove symlink");
        check(CreateSymbolicLinkA("dangling", "absent", 0), "cannot create dangling symlink");
        assert_rejected("dangling", ELOOP);
        check(DeleteFileA("dangling"), "cannot remove dangling symlink");
    } else {
        puts("SKIP: symlink creation unavailable");
    }
#else
    assert(!mbe_quality_paths_equal("absent", "ABSENT"));
    check(symlink("first", "symlink") == 0, "cannot create symlink");
    assert(mbe_quality_same_file("first", "symlink"));
    assert_normalized_equal("first", "symlink");
    check(remove("symlink") == 0, "cannot remove symlink");
    check(symlink("child", "parent-link") == 0, "cannot create parent symlink");
    assert_normalized_equal("child/absent", "parent-link/absent");
    check(remove("parent-link") == 0, "cannot remove parent symlink");
    check(symlink("absent", "dangling") == 0, "cannot create dangling symlink");
    assert_rejected("dangling", ELOOP);
    check(remove("dangling") == 0, "cannot remove dangling symlink");
#endif
    free(parent);
    check(remove("first") == 0 && remove("second") == 0, "cannot remove test files");
#if defined(_WIN32)
    check(RemoveDirectoryA("child"), "cannot remove child directory");
    check(SetCurrentDirectoryA(original), "cannot restore working directory");
    check(RemoveDirectoryA(directory), "cannot remove temporary directory");
#else
    check(rmdir("child") == 0, "cannot remove child directory");
    check(chdir(original) == 0, "cannot restore working directory");
    check(rmdir(directory) == 0, "cannot remove temporary directory");
#endif
    return 0;
}
