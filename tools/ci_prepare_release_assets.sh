#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Package/SBOM basenames are preserved by actions/upload-artifact.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python3 "$script_dir/release_assets.py" --prefix "$ARCHIVE_PREFIX"
