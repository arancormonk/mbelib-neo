#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""tools/shell_lint.sh fails when either shellcheck or shfmt does.

The wrapper ran both tools in one group piped through tee, and a group's status is
its last command's, so a script that only shellcheck rejected passed the lint.
"""

import os
import shutil
import subprocess
import sys
import tempfile

SCRIPTS = {
    # Passes both tools.
    "clean.sh": "#!/usr/bin/env bash\nset -euo pipefail\necho \"ok\"\n",
    # shellcheck SC2034 (assigned, never used); shfmt has nothing to change.
    "shellcheck_only.sh": "#!/usr/bin/env bash\nunused_value=1\n",
    # Clean for shellcheck; shfmt wants two-space indentation.
    "shfmt_only.sh": "#!/usr/bin/env bash\nif true; then\n    echo \"ok\"\nfi\n",
}


def main():
    # usage: test_shell_lint.py <tools/shell_lint.sh> [bash git shellcheck shfmt]
    # Tools not given are looked up on PATH.
    wrapper = os.path.abspath(sys.argv[1])
    given = sys.argv[2:6]
    names = ("bash", "git", "shellcheck", "shfmt")
    bash, git, shellcheck, shfmt = (
        given[index] if index < len(given) else shutil.which(name) for index, name in enumerate(names)
    )
    assert all((bash, git, shellcheck, shfmt)), "bash, git, shellcheck and shfmt are required"
    tools = [git, shellcheck, shfmt]
    with tempfile.TemporaryDirectory(prefix="mbe_shell_lint_") as repo:
        subprocess.run([git, "init", "-q"], cwd=repo, check=True, timeout=60)
        for name, text in SCRIPTS.items():
            with open(os.path.join(repo, name), "w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
        # The wrapper finds its tools on PATH: put the ones this test was given first.
        path = os.pathsep.join([*(os.path.dirname(os.path.abspath(tool)) for tool in tools),
                                os.environ.get("PATH", "")])
        env = dict(os.environ, PATH=path)
        for name, expected in (("clean.sh", 0), ("shellcheck_only.sh", 1), ("shfmt_only.sh", 1)):
            result = subprocess.run(
                [bash, wrapper, "--", name], cwd=repo, env=env, capture_output=True, text=True, timeout=120,
            )
            assert result.returncode == expected, (name, result.returncode, result.stdout, result.stderr)
            print(f"PASS {name}: shell lint exit {result.returncode}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
