#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
# Run from the repository root. Python 3 plus git (legacy) or ffmpeg (heldout).
set -euo pipefail

if [[ ! -f tools/quality/fetch_speech_corpus.py ]]; then
  echo "usage: tools/quality/fetch_speech_corpus.sh [--set legacy|heldout] (from repository root)" >&2
  exit 2
fi

exec python3 tools/quality/fetch_speech_corpus.py "$@"
