# Code Quality Guardrails

mbelib-neo relies on layered checks because no single compiler, analyzer, or review pass catches every bad change. These guardrails apply to hand-written, copied, generated, and bulk-edited code equally.

## Review Expectations

- A human owns every submitted change and is responsible for understanding the behavior, failure modes, and module boundaries it touches.
- Broad or risky changes need a second reviewer. Treat codec, DSP, API/ABI, workflow, dependency, and release changes as risky by default.
- Copied, generated, or bulk-written code must be adapted to local conventions before review, not accepted as a black box.
- Static-analysis suppressions must explain why the warning is a false positive or why the local exception is acceptable. Keep suppressions close to the narrowest affected code.
- Do not add unrelated refactors to quality or security fixes. Smaller diffs make analyzer output and review results more reliable.

## High-Risk Change Checklist

Use this checklist when a change touches decoder logic, external input, allocation, dependencies, workflows, or packaging:

- Add or update focused tests under the matching `tests/` area.
- For new public APIs, verify headers live under `include/mbelib-neo/` and are included as `<mbelib-neo/...>`.
- For new dependencies, document the dependency in CMake, README/install notes, and license notices as applicable.
- For workflow changes, keep `GITHUB_TOKEN` permissions least-privilege and avoid interpolating untrusted GitHub context directly into shell scripts.
- For release changes, verify artifact and checksum steps still run on tag builds.
- For dependency or workflow changes, also check the supply-chain policy in `docs/supply-chain-guardrails.md`.

## Required Local Checks

Run the smallest useful set before opening a PR, then broaden it when the change is risky.

- Normal C changes: `cmake --build --preset dev-debug -j` and `ctest --preset dev-debug --output-on-failure`.
- Normal pre-push check: `tools/preflight_ci.sh`.
- Broad or high-risk changes: `tools/quality_preflight.sh`.
- Sanitizer-sensitive code: `ctest --preset asan-ubsan-debug --output-on-failure` after configuring/building the matching preset.
- CMake changes: `tools/cmake_format_check.sh`.
- Workflow changes: `tools/workflow_lint.sh` and `tools/zizmor.sh`.
- Dependency input changes: `tools/osv_scan.sh`.

## Project-Specific Guardrails

The repository intentionally blocks or flags patterns that are easy to reintroduce during large edits:

- Do not execute shells or spawn processes from project-owned C without explicit design review.
- Do not include bundled third-party headers directly outside approved integration points.
- Keep workflow scripts defensive: pass untrusted context through environment variables or action inputs, not direct expression interpolation in `run:` blocks.
- Keep analyzer and linter output actionable. Prefer fixing root causes over widening suppressions.
