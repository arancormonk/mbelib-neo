#!/usr/bin/env bash
set -euo pipefail

# Run Semgrep for additional SAST checks. Default mode is advisory (non-erroring).
# Use --strict to fail on findings.
# Excludes third-party code under src/external. Strict mode also loads
# project-specific guardrail rules from semgrep/mbelib-neo.yml.

ROOT_DIR=$(git rev-parse --show-toplevel 2> /dev/null || pwd)
cd "$ROOT_DIR"

usage() {
  cat << 'USAGE'
Usage: tools/semgrep.sh [--strict] [--config <config>] [--] [paths...]

Options:
  --strict          Fail on findings (--error).
  --config CONFIG   Semgrep config/rule pack (default: p/default; strict also
                    adds p/c, p/security-audit, and semgrep/mbelib-neo.yml).
                    May be supplied multiple times.

Arguments:
  paths...          Optional paths to scan. Default: src include examples bench tests fuzz tools .github/workflows

Environment:
  MBE_SEMGREP_SARIF_OUT  Optional SARIF output path for GitHub code scanning.
USAGE
}

STRICT=0
CONFIGS=("p/default")
CUSTOM_CONFIGS=0
TARGETS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --strict)
      STRICT=1
      shift
      ;;
    --config)
      if [[ $# -lt 2 ]]; then
        echo "Missing value for --config" >&2
        exit 2
      fi
      if [[ $CUSTOM_CONFIGS -eq 0 ]]; then
        CONFIGS=()
        CUSTOM_CONFIGS=1
      fi
      CONFIGS+=("$2")
      shift 2
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    --)
      shift
      TARGETS+=("$@")
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      TARGETS+=("$1")
      shift
      ;;
  esac
done

if ! command -v semgrep > /dev/null 2>&1; then
  echo "semgrep not found. Install with: pipx install semgrep (or pip install semgrep)." >&2
  exit 1
fi

if [[ ${#TARGETS[@]} -eq 0 ]]; then
  TARGETS=(src include examples bench tests fuzz tools .github/workflows)
fi

LOG_FILE=".semgrep.local.out"

ARGS=(
  --metrics=off
  --disable-version-check
  --no-git-ignore
  --exclude .deps-arch-toolchain
  --exclude .deps-arch-toolchain/**
  --exclude src/external
  --exclude src/external/**
  --exclude build
  --exclude build/**
)
if [[ $STRICT -eq 1 ]]; then
  ARGS+=(--error)
  if [[ $CUSTOM_CONFIGS -eq 0 ]]; then
    CONFIGS+=(p/c p/security-audit semgrep/mbelib-neo.yml)
  fi
fi
for config in "${CONFIGS[@]}"; do
  ARGS+=(--config "$config")
done
if [[ -n "${MBE_SEMGREP_SARIF_OUT:-}" ]]; then
  ARGS+=(--sarif --output "$MBE_SEMGREP_SARIF_OUT")
fi

run_semgrep() {
  local append_log=0
  local rc=0
  local -a tee_args=("$LOG_FILE")

  if [[ "${1:-}" == "--append-log" ]]; then
    append_log=1
    shift
  fi
  if [[ $append_log -eq 1 ]]; then
    tee_args=(-a "$LOG_FILE")
  fi

  set +e
  semgrep "$@" 2>&1 | tee "${tee_args[@]}"
  rc=${PIPESTATUS[0]}
  set -e
  return "$rc"
}

if run_semgrep "${ARGS[@]}" "${TARGETS[@]}"; then
  rc=0
else
  rc=$?
fi

if [[ $rc -ne 0 ]] && grep -q 'io_uring_queue_init' "$LOG_FILE"; then
  echo "Semgrep failed while initializing io_uring; retrying with legacy parallelism." | tee -a "$LOG_FILE" >&2
  if run_semgrep --append-log --x-parmap "${ARGS[@]}" "${TARGETS[@]}"; then
    rc=0
  else
    rc=$?
  fi
fi

if [[ $rc -eq 0 ]]; then
  echo "Semgrep completed. Full output in $LOG_FILE"
else
  echo "Semgrep found issues (or failed). See $LOG_FILE for details." >&2
fi

exit "$rc"
