#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
set -euo pipefail

ROOT_DIR=$(git rev-parse --show-toplevel 2> /dev/null || pwd)
cd "$ROOT_DIR"

export MBE_HOOK_FAIL_ON_MISSING_TOOLS=1
export MBE_HOOK_RUN_SCAN_BUILD=1

tools/preflight_ci.sh "$@"
python3 tools/legal_inventory.py
python3 tests/test_legal_release.py
tools/cmake_format_check.sh
tools/shell_lint.sh
tools/workflow_lint.sh
tools/check_workflow_git_pins.sh
tools/check_workflow_download_pins.sh
tools/zizmor.sh
tools/osv_scan.sh
tools/gitleaks.sh
