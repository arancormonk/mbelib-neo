#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the paths tools/push_changed_files.sh lists for .githooks/pre-push.

Computed inline with `git diff --name-only ... || true`, a remote SHA the clone
did not have made the diff fail, which read as "nothing changed", so the hook
exited without checking anything; and git's quoting of a name with non-ASCII
bytes, quotes or tabs kept that file out of every check. This builds a scratch
repository whose history gives each way of finding the comparison a different
answer, and checks each one, the raw names, a remote named like an option, and
that a failing diff fails.
"""

import os
import subprocess
import sys
import tempfile

ZEROS = "0" * 40
MISSING = "1234567890abcdef1234567890abcdef12345678"
ODD_NAMES = ["src/café.c", "src/quote\"name.c", "src/tab\tname.c"]


def git(git_exe, repo, *args):
    result = subprocess.run([git_exe, *args], cwd=repo, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def write(repo, path, text):
    full = os.path.join(repo, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def commit(git_exe, repo, message):
    git(git_exe, repo, "add", "-A")
    git(git_exe, repo, "commit", "-q", "-m", message)
    return git(git_exe, repo, "rev-parse", "HEAD")


def list_paths(bash, git_exe, script, repo, remote, head, remote_sha):
    path = os.path.dirname(os.path.abspath(git_exe)) + os.pathsep + os.environ.get("PATH", "")
    result = subprocess.run(
        [bash, script, remote, "refs/heads/topic", head, remote_sha],
        cwd=repo, env=dict(os.environ, PATH=path), capture_output=True, timeout=60,
    )
    return result.returncode, result.stdout, result.stderr.decode("utf-8", "replace")


def main():
    bash, git_exe, script = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
    with tempfile.TemporaryDirectory(prefix="mbe_push_changed_files_") as repo:
        git(git_exe, repo, "init", "-q")
        git(git_exe, repo, "config", "user.name", "push-changed-files test")
        git(git_exe, repo, "config", "user.email", "push-changed-files@example.invalid")
        git(git_exe, repo, "config", "commit.gpgsign", "false")
        # The quoting case needs git's default quoting, whatever the global config says.
        git(git_exe, repo, "config", "core.quotePath", "true")

        # c0 -> c1 (x.c) -> c2 (y.c) -> c3 (z.c and the odd names) is the pushed branch; side is a commit off c1
        # that only a remote branch has.
        for name in ("keep", "x", "y", "z"):
            write(repo, f"src/{name}.c", f"int {name};\n")
        c0 = commit(git_exe, repo, "c0")
        write(repo, "src/x.c", "int x1;\n")
        c1 = commit(git_exe, repo, "c1")
        write(repo, "src/side.c", "int side;\n")
        side = commit(git_exe, repo, "side")
        git(git_exe, repo, "checkout", "-q", "--detach", c1)
        write(repo, "src/y.c", "int y2;\n")
        commit(git_exe, repo, "c2")
        write(repo, "src/z.c", "int z3;\n")
        for name in ODD_NAMES:
            write(repo, name, "int odd;\n")
        head = commit(git_exe, repo, "c3")

        def expect_paths(label, remote, remote_sha, *paths):
            returncode, stdout, stderr = list_paths(bash, git_exe, script, repo, remote, head, remote_sha)
            assert returncode == 0, (label, stderr)
            got = sorted(name.decode("utf-8") for name in stdout.split(b"\0") if name)
            assert got == sorted([*paths, *ODD_NAMES]), (label, got)
            print(f"PASS {label}")
            return stderr

        git(git_exe, repo, "update-ref", "refs/remotes/origin/main", c0)
        expect_paths("existing ref: compared with the remote's SHA", "origin", c1, "src/y.c", "src/z.c")
        expect_paths("new ref: compared with the remote's default branch", "origin", ZEROS,
                     "src/x.c", "src/y.c", "src/z.c")
        stderr = expect_paths("remote SHA missing from the clone: compared with the default branch", "origin",
                              MISSING, "src/x.c", "src/y.c", "src/z.c")
        assert "is not in this clone" in stderr, stderr
        expect_paths("a remote named -h is a remote, not a request for help", "-h", c1, "src/y.c", "src/z.c")

        git(git_exe, repo, "update-ref", "-d", "refs/remotes/origin/main")
        git(git_exe, repo, "update-ref", "refs/remotes/origin/near", side)
        git(git_exe, repo, "update-ref", "refs/remotes/origin/far", c0)
        expect_paths("no default branch: compared with the nearest merge base", "origin", ZEROS, "src/y.c", "src/z.c")

        git(git_exe, repo, "update-ref", "-d", "refs/remotes/origin/near")
        git(git_exe, repo, "update-ref", "-d", "refs/remotes/origin/far")
        expect_paths("no remote branches: compared with the empty tree", "origin", ZEROS,
                     "src/keep.c", "src/x.c", "src/y.c", "src/z.c")

        git(git_exe, repo, "config", "diff.context", "notanumber")
        returncode, _, stderr = list_paths(bash, git_exe, script, repo, "origin", head, c1)
        assert returncode != 0, stderr
        assert "failed for refs/heads/topic" in stderr, stderr
        print("PASS a failing diff fails instead of listing nothing")
    return 0


if __name__ == "__main__":
    sys.exit(main())
