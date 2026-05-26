#!/usr/bin/env bash
set -euo pipefail

artifact=${1:-build/dev-release/libmbe-neo.so}
build_dir=${2:-build/dev-release}
compile_db="${build_dir}/compile_commands.json"

if [[ ! -e "$artifact" ]]; then
  echo "release artifact not found: $artifact" >&2
  exit 2
fi

if [[ -f "$compile_db" ]]; then
  if ! grep -q -- '-D_FORTIFY_SOURCE=2' "$compile_db"; then
    echo "release compile commands are missing -D_FORTIFY_SOURCE=2" >&2
    exit 1
  fi
  if ! grep -q -- '-fstack-protector-strong' "$compile_db"; then
    echo "release compile commands are missing -fstack-protector-strong" >&2
    exit 1
  fi
fi

if [[ "$(uname -s)" != "Linux" ]]; then
  echo "ELF hardening checks are Linux-only; compile flag checks completed."
  exit 0
fi

if ! command -v readelf > /dev/null 2>&1; then
  echo "readelf is required for ELF hardening checks." >&2
  exit 2
fi

if ! readelf -h "$artifact" | grep -Eq 'Type:[[:space:]]*DYN'; then
  echo "release artifact is not position-independent (ELF type DYN expected)." >&2
  exit 1
fi

if ! readelf -l "$artifact" | grep -q 'GNU_RELRO'; then
  echo "release artifact is missing GNU_RELRO." >&2
  exit 1
fi

if ! readelf -d "$artifact" | grep -Eq 'BIND_NOW|FLAGS_1.*NOW'; then
  echo "release artifact is missing BIND_NOW/full RELRO." >&2
  exit 1
fi
