#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
# Linux-only controlled loader selection; Python is developer tooling, not a library dependency.
set -euo pipefail
export MBE_QUALITY_RUNNER_ARGV0="$0"
exec python3 "$(dirname -- "${BASH_SOURCE[0]}")/run_quality_ab.py" "$@"
