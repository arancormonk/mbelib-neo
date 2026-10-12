#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Package/SBOM basenames are preserved by actions/upload-artifact.
set -euo pipefail
mapfile -t packages < <(
  find dist -maxdepth 1 -type f \
    \( -name "${ARCHIVE_PREFIX}-*.tar.gz" -o -name "${ARCHIVE_PREFIX}-*.zip" \) |
    sort
)
if [ "${#packages[@]}" -ne 3 ]; then
  echo "Expected 3 release packages, found ${#packages[@]}." >&2
  find dist -maxdepth 1 -type f -print >&2
  exit 1
fi
for platform in linux macos windows; do
  if ! printf '%s\n' "${packages[@]}" | grep -Fq "/${ARCHIVE_PREFIX}-${platform}-"; then
    echo "Missing ${platform} release package." >&2
    exit 1
  fi
done
for package in "${packages[@]}"; do
  package_name="$(basename "$package")"
  package_stem="${package_name%.tar.gz}"
  package_stem="${package_stem%.zip}"
  case "${package_stem#"${ARCHIVE_PREFIX}-"}" in
    linux-x64 | linux-arm64 | macos-x64 | macos-arm64 | windows-x64 | windows-arm64) ;;
    *)
      echo "Unsupported release package platform or architecture: ${package_name}" >&2
      exit 1
      ;;
  esac
done
mapfile -t sboms < <(
  find dist -maxdepth 1 -type f -name "${ARCHIVE_PREFIX}-*.spdx.json" |
    sort
)
if [ "${#sboms[@]}" -ne 3 ]; then
  echo "Expected 3 release package SBOMs, found ${#sboms[@]}." >&2
  find dist -maxdepth 1 -type f -print >&2
  exit 1
fi
for package in "${packages[@]}"; do
  if [ ! -f "${package}.spdx.json" ]; then
    echo "Missing SBOM for $(basename "$package")." >&2
    exit 1
  fi
done
: > dist/SHA256SUMS.txt
for package in "${packages[@]}"; do
  (
    cd dist
    sha256sum "$(basename "$package")"
  ) >> dist/SHA256SUMS.txt
done
cat dist/SHA256SUMS.txt
