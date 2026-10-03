#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Run .githooks/pre-push on a push it must refuse and on one with nothing to check.

The hook hands its checks one path per line, so a pushed name with a newline in it
would split into two paths that do not exist: skipped with a warning for a ref that
is not checked out, or reported as local deletions for the one that is. The hook
must refuse the push and say why, and leave nothing behind in the temporary
directory. Both cases end before the hook runs any check, so the scratch repository
holds only the hook and tools/push_changed_files.sh.
"""

import os
import shutil
import subprocess
import sys
import tempfile


def git(git_exe, repo, *args):
    result = subprocess.run([git_exe, *args], cwd=repo, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def main():
    bash, git_exe, source_root = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
    with tempfile.TemporaryDirectory(prefix="mbe_pre_push_paths_") as scratch:
        repo = os.path.join(scratch, "repo")
        hook_tmp = os.path.join(scratch, "tmp")
        os.makedirs(os.path.join(repo, ".githooks"))
        os.makedirs(os.path.join(repo, "tools"))
        os.makedirs(hook_tmp)
        shutil.copy2(os.path.join(source_root, ".githooks", "pre-push"), os.path.join(repo, ".githooks", "pre-push"))
        shutil.copy2(os.path.join(source_root, "tools", "push_changed_files.sh"),
                     os.path.join(repo, "tools", "push_changed_files.sh"))
        git(git_exe, repo, "init", "-q")
        git(git_exe, repo, "config", "user.name", "pre-push test")
        git(git_exe, repo, "config", "user.email", "pre-push@example.invalid")
        git(git_exe, repo, "config", "commit.gpgsign", "false")
        with open(os.path.join(repo, "base.c"), "w", encoding="utf-8") as handle:
            handle.write("int base;\n")
        git(git_exe, repo, "add", "-A")
        git(git_exe, repo, "commit", "-q", "-m", "base")
        base = git(git_exe, repo, "rev-parse", "HEAD")

        path = os.path.dirname(os.path.abspath(git_exe)) + os.pathsep + os.environ.get("PATH", "")
        env = dict(os.environ, PATH=path, TMPDIR=hook_tmp)

        def run_hook(local_sha, remote_sha):
            line = f"refs/heads/topic {local_sha} refs/heads/topic {remote_sha}\n"
            return subprocess.run([bash, ".githooks/pre-push", "origin"], cwd=repo, env=env, input=line,
                                  capture_output=True, text=True, timeout=120)

        result = run_hook(base, base)
        assert result.returncode == 0, result.stdout + result.stderr
        assert not os.listdir(hook_tmp), os.listdir(hook_tmp)
        print("PASS a push with nothing to check passes and leaves no temporary file")

        with open(os.path.join(repo, "line\nbreak.c"), "w", encoding="utf-8") as handle:
            handle.write("int newline;\n")
        git(git_exe, repo, "add", "-A")
        git(git_exe, repo, "commit", "-q", "-m", "newline")
        result = run_hook(git(git_exe, repo, "rev-parse", "HEAD"), base)
        assert result.returncode != 0, result.stdout + result.stderr
        assert "path with a newline in it" in result.stderr, result.stderr
        assert not os.listdir(hook_tmp), os.listdir(hook_tmp)
        print("PASS a pushed name with a newline is refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
