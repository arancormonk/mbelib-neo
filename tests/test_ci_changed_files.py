#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check which pull-request changes tools/ci_changed_files.sh counts as fuzz targets.

PR fuzzing runs on every pull request so its checks can be required, and skips the
fuzzers when this count is zero. It must count exactly what the workflow's former
`paths:` filter matched: deletions, both sides of a rename and src/external/ included.
"""

import os
import subprocess
import sys
import tempfile

BASE_FILES = {
    "CMakeLists.txt": "project(t C)\n",
    ".clusterfuzzlite/build.sh": "#!/bin/bash\n",
    ".github/workflows/cflite_pr.yml": "name: fuzz\n",
    ".github/workflows/other.yml": "name: other\n",
    "docs/notes.md": "notes\n",
    "fuzz/fuzz_frame.c": "int f;\n",
    "include/mbelib-neo/mbelib.h": "int h;\n",
    "src/core/mbelib.c": "int c;\n",
    "src/external/pffft/pffft.c": "int p;\n",
    "sub/CMakeLists.txt": "add_library(s s.c)\n",
    "tools/helper.sh": "#!/bin/bash\n",
}


def run(git, repo, *args):
    result = subprocess.run([git, *args], cwd=repo, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def write(repo, path, text):
    full = os.path.join(repo, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def commit_all(git, repo, message):
    run(git, repo, "add", "-A")
    run(git, repo, "commit", "-q", "-m", message)
    return run(git, repo, "rev-parse", "HEAD")


def run_helper(bash, script, repo, base, head):
    """Run the helper from the repository and return (returncode, GITHUB_OUTPUT values, stderr)."""
    # Outside the repository, so no case commits another's output.
    scratch = os.path.dirname(repo)
    out_dir = os.path.join(scratch, "ci-out")
    github_output = os.path.join(scratch, "github-output")
    if os.path.exists(github_output):
        os.remove(github_output)
    # The helper runs bare `git`: put the git this test was given first on its PATH.
    path = os.path.dirname(os.path.abspath(GIT)) + os.pathsep + os.environ.get("PATH", "")
    env = dict(os.environ, GITHUB_OUTPUT=github_output, PATH=path)
    result = subprocess.run(
        [bash, script, "--base", base, "--head", head, "--no-header-expansion", "--out-dir", out_dir],
        cwd=repo, env=env, capture_output=True, text=True, timeout=60,
    )
    outputs = {}
    if os.path.exists(github_output):
        with open(github_output, encoding="utf-8") as handle:
            outputs = dict(line.split("=", 1) for line in handle.read().splitlines() if "=" in line)
    return result.returncode, outputs, result.stderr


def listed_targets(repo, outputs):
    with open(os.path.join(os.path.dirname(repo), "ci-out", "fuzz_targets.txt"), encoding="utf-8") as handle:
        listed = [line for line in handle.read().splitlines() if line]
    assert int(outputs["fuzz_targets"]) == len(listed), (outputs, listed)
    return listed


def fuzz_targets(bash, git, script, repo, base, change):
    """Branch from base, apply change(repo), and return (count, listed paths, outputs)."""
    run(git, repo, "checkout", "-q", "--detach", base)
    change(repo)
    commit_all(git, repo, "change")
    returncode, outputs, stderr = run_helper(bash, script, repo, base, "HEAD")
    assert returncode == 0, stderr
    listed = listed_targets(repo, outputs)
    return int(outputs["fuzz_targets"]), listed, outputs


def modify(path):
    return lambda repo: write(repo, path, "changed\n")


def delete(path):
    return lambda repo: os.remove(os.path.join(repo, path))


GIT = "git"


def main():
    global GIT
    bash, git, script = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
    GIT = git
    with tempfile.TemporaryDirectory(prefix="mbe_ci_changed_files_") as scratch:
        repo = os.path.join(scratch, "repo")
        os.mkdir(repo)
        run(git, repo, "init", "-q")
        run(git, repo, "config", "user.name", "ci-changed-files test")
        run(git, repo, "config", "user.email", "ci-changed-files@example.invalid")
        run(git, repo, "config", "commit.gpgsign", "false")
        for path, text in BASE_FILES.items():
            write(repo, path, text)
        base = commit_all(git, repo, "base")

        cases = [
            ("docs only", modify("docs/notes.md"), []),
            ("library source", modify("src/core/mbelib.c"), ["src/core/mbelib.c"]),
            ("public header", modify("include/mbelib-neo/mbelib.h"), ["include/mbelib-neo/mbelib.h"]),
            ("fuzz harness", modify("fuzz/fuzz_frame.c"), ["fuzz/fuzz_frame.c"]),
            ("fuzz build script", modify(".clusterfuzzlite/build.sh"), [".clusterfuzzlite/build.sh"]),
            ("this workflow", modify(".github/workflows/cflite_pr.yml"), [".github/workflows/cflite_pr.yml"]),
            ("another workflow", modify(".github/workflows/other.yml"), []),
            ("top-level CMakeLists.txt", modify("CMakeLists.txt"), ["CMakeLists.txt"]),
            ("nested CMakeLists.txt", modify("sub/CMakeLists.txt"), []),
            ("tool script", modify("tools/helper.sh"), []),
            ("vendored source", modify("src/external/pffft/pffft.c"), ["src/external/pffft/pffft.c"]),
            ("deleted header", delete("include/mbelib-neo/mbelib.h"), ["include/mbelib-neo/mbelib.h"]),
        ]
        for label, change, expected in cases:
            count, listed, _ = fuzz_targets(bash, git, script, repo, base, change)
            assert listed == expected, (label, listed, expected)
            assert count == len(expected), (label, count)
            print(f"PASS {label}: fuzz_targets={count}")

        def rename_out(repo_dir):
            os.makedirs(os.path.join(repo_dir, "lib"), exist_ok=True)
            os.rename(os.path.join(repo_dir, "src/core/mbelib.c"), os.path.join(repo_dir, "lib/mbelib.c"))

        count, listed, outputs = fuzz_targets(bash, git, script, repo, base, rename_out)
        assert listed == ["src/core/mbelib.c"], listed
        print(f"PASS rename out of src/: fuzz_targets={count}")

        # The vendored tree stays out of the other target sets, as before.
        _, _, outputs = fuzz_targets(bash, git, script, repo, base, modify("src/external/pffft/pffft.c"))
        assert outputs["semgrep_targets"] == "0" and outputs["format_files"] == "0", outputs
        print("PASS vendored source stays out of the analysis target sets")

        # The selection follows the pull request's own changes (base...head), as GitHub's paths filter did. An
        # edit main already made the same way stays in scope against the head; against the merge commit it would
        # vanish, which is why cflite_pr.yml passes the pull request's head.
        run(git, repo, "checkout", "-q", "--detach", base)
        write(repo, "src/core/mbelib.c", "same edit\n")
        main_tip = commit_all(git, repo, "main makes the edit")
        run(git, repo, "checkout", "-q", "--detach", base)
        write(repo, "src/core/mbelib.c", "same edit\n")
        write(repo, "docs/notes.md", "and docs\n")
        pr_head = commit_all(git, repo, "pull request makes it too")
        returncode, outputs, stderr = run_helper(bash, script, repo, main_tip, pr_head)
        assert returncode == 0, stderr
        assert listed_targets(repo, outputs) == ["src/core/mbelib.c"], outputs
        run(git, repo, "checkout", "-q", "--detach", main_tip)
        run(git, repo, "merge", "-q", "--no-edit", "--no-ff", pr_head)
        merge = run(git, repo, "rev-parse", "HEAD")
        returncode, outputs, stderr = run_helper(bash, script, repo, main_tip, merge)
        assert returncode == 0 and listed_targets(repo, outputs) == [], (outputs, stderr)
        print("PASS an edit main already has stays in scope against the pull request head")

    return 0


if __name__ == "__main__":
    sys.exit(main())
