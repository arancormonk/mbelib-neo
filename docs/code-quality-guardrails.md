# Code Quality Guardrails

mbelib-neo relies on layered checks because no single compiler, analyzer, or review pass catches every bad change. These guardrails apply to hand-written, copied, generated, and bulk-edited code equally.

## Review Expectations

- A human owns every submitted change and is responsible for understanding the behavior, failure modes, and module boundaries it touches.
- Broad or risky changes should get a second reviewer when practical. For solo-maintainer changes, treat passing guardrails plus documented risk as the minimum review baseline.
- Copied, generated, or bulk-written code must be adapted to local conventions before review, not accepted as a black box.
- Static-analysis suppressions must explain why the warning is a false positive or why the local exception is acceptable. Keep suppressions close to the narrowest affected code.
- Do not add unrelated refactors to quality or security fixes. Smaller diffs make analyzer output and review results more reliable.

## High-Risk Change Checklist

Use this checklist when a change touches decoder logic, external input, allocation, dependencies, workflows, or packaging:

- Add or update focused tests under the matching `tests/` area.
- For new public APIs, verify headers live under `include/mbelib-neo/` and are included as `<mbelib-neo/...>`.
- For new dependencies, document the dependency in CMake, README/install notes, and license notices as applicable.
- For workflow changes, keep `GITHUB_TOKEN` permissions least-privilege and avoid interpolating untrusted GitHub context directly into shell scripts.
- For release changes, verify signed-tag validation, SBOM, artifact, attestation, checksum, and release-hardening checks still run on tag builds.
- For dependency or workflow changes, also check the supply-chain policy in `docs/supply-chain-guardrails.md`.

## Required Local Checks

Run the smallest useful set before opening a PR, then broaden it when the change is risky.

- Normal C changes: `cmake --build --preset dev-debug -j` and `ctest --preset dev-debug --output-on-failure`.
- Normal pre-push check: `tools/preflight_ci.sh`. The hook checks the files `tools/push_changed_files.sh` lists
  for each pushed ref. When the remote's commit is not in the clone it compares with the remote's default branch,
  and when no comparison works the push fails instead of going out unchecked. CI's pull-request jobs take their
  files from `tools/ci_changed_files.sh` the same way: names stay raw (`git diff -z`), and a diff that fails fails
  the job.
- Broad or high-risk changes: `tools/quality_preflight.sh`.
- Sanitizer-sensitive code: `ctest --preset asan-ubsan-debug --output-on-failure` after configuring/building the matching preset.
- CMake changes: `tools/cmake_format_check.sh`.
- Workflow changes: `tools/workflow_lint.sh` and `tools/zizmor.sh`.
- Workflow source/download pin changes: `tools/check_workflow_git_pins.sh` and `tools/check_workflow_download_pins.sh`.
- Release hardening changes: run `tools/check_release_hardening.sh build/dev-release/libmbe-neo.so build/dev-release` on Linux and let CI verify macOS Mach-O and Windows PE release artifacts.
- Dependency input changes: `tools/osv_scan.sh`.

## Required PR Status Checks

Branch protection must require the repository guardrail job for pull requests.
Keep `repository guardrails` required whenever new repository-level security
guardrails are added.

## CI Runner Budget

All runner jobs share 18 fixed, repository-wide concurrency groups. This leaves
headroom below the account's 20-runner limit for dispatch and cancellation. It is
an independent budget for mbelib-neo, not a reservation against dsd-neo's budget:
GitHub concurrency groups do not coordinate between repositories. If both run at
once, the account limit still applies.

Each group runs one job and uses `cancel-in-progress: false` with `queue: max`,
so a new main run does not replace pending PR jobs. GitHub allows up to 100 queued
jobs per group; the queue is bounded, and dispatch order between groups is not
guaranteed. See [GitHub's concurrency documentation](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/control-workflow-concurrency).
Workflow-level concurrency still cancels an older run of the same workflow and
PR/ref. CI additionally separates event types, so manual release recovery cannot
cancel the regular main push checks. Matrix `max-parallel` limits the largest bursts within a run; the fixed
groups enforce the budget across workflows and overlapping PR/main runs.

| Group | Jobs sharing the runner |
| --- | --- |
| `ci-worker-01` | Linux Debug, cppcheck, Clang Release, Linux release package |
| `ci-worker-02` | Linux Release, optimized GCC tests, Linux install/consumer |
| `ci-worker-03` | Linux NOTONES, consumer-defaults GCC tests, warnings build |
| `ci-worker-04` | macOS Debug, macOS install/consumer, macOS release package |
| `ci-worker-05` | macOS Release |
| `ci-worker-06` | macOS NOTONES |
| `ci-worker-07` | Windows Debug, Windows install/consumer, Windows release package |
| `ci-worker-08` | Windows Release, ARM32 default-FPU tests |
| `ci-worker-09` | Windows NOTONES, ARM32 NEON tests |
| `ci-worker-10` | clang-tidy/IWYU/GCC analyzer, C/C++ CodeQL |
| `ci-worker-11` | scan-build, ASAN/UBSAN tests |
| `ci-worker-12` | optimized Clang tests, address-sanitizer fuzzing |
| `ci-worker-13` | AVX-512 link check, undefined-behavior fuzzing |
| `ci-fast-1` | PR classification, Semgrep, Lizard, OSV, Python CodeQL |
| `ci-fast-2` | formatting, other guardrails, Actions CodeQL, dependency review, Scorecard, release publication |
| `ci-control` | closed-PR cleanup, release-tag validation |
| `aur-mbelib-neo-git`, `aur-mbelib-neo` | the existing AUR publication locks |

`tools/workflow_lint.sh` also runs `tools/check_ci_concurrency.py`. Every runner
job must use a fixed group from this table and a non-replacing queue. New matrix
entries must declare `ci_group`; only the budgeted local classifier may be called
as a reusable workflow. actionlint 1.7.12 does not recognize `queue` yet, so its
exact unknown-key diagnostic is suppressed after the budget checker validates
that setting. Other actionlint diagnostics still fail.

PR workflows classify paths before allocating expensive jobs. Only known prose
(`docs/*.md` and the root README, CHANGELOG, CONTRIBUTING, SECURITY and AGENTS
Markdown files) skips native builds and their analysis. Unknown paths, source,
build, tool and workflow inputs retain full coverage, including deletions and
both sides of renames. Classification uses the PR's own base/head diff, even when
main already contains an equivalent change. A failed classifier or missing
output selects full checks. Main, tag, scheduled and manual runs retain their
full checks; release validation, hardening and publication prerequisites remain.

The six existing analysis gates still precede native builds. Standalone jobs can
skip before allocating a runner. Required matrix jobs always expand when those
gates permit them, with an explicit successful step for irrelevant changes, so
every existing check name is reported. Relevant CMake formatting checks scan the
whole tree. CodeQL, repository/script guardrails, Lizard, secret scanning and
dependency review continue on every PR. Fuzzing retains its library/harness/build
selection and also covers changes to the shared classifier.

Closing or merging a PR starts `ci-cancel-closed-pr.yml`. Its small job has
`actions: write` only for cancellation and executes the cleanup helper from the
trusted default branch, never the PR. It snapshots active `pull_request` runs,
matches the repository and PR identity (or source repository, branch and PR
lifetime when GitHub clears the association after merge), and preserves open or
reopened PRs. Push, release, scheduled and unrelated PR runs are excluded. This
hook becomes active after the workflow is merged into the default branch.

The October 2026 review found 39 PR jobs, including a burst of 20 build/test jobs
after analysis. PR #133's CI continued for 37 minutes after merge, using about
21 runner-minutes after closure. Its overlapping main CI took 37 minutes versus
7 minutes for a quiet main run earlier that day. Native test steps were only
5–27 seconds (about a minute under ARM emulation), so CTest remains serial.
Scorecard already runs directly without a polling runner. After rollout, compare
queue time and completion time for overlapping PR/main runs, confirm all 40
required contexts appear, and confirm closed PRs stop consuming runners.

## Project-Specific Guardrails

The repository intentionally blocks or flags patterns that are easy to reintroduce during large edits:

- Do not execute shells or spawn processes from project-owned C without explicit design review.
- Do not include bundled third-party headers directly outside approved integration points.
- Keep public GitHub source fetches in CI pinned through `tools/ci-dependency-pins.env` and `tools/fetch-pinned-git.sh`.
- Keep CI fallback container images pinned by digest and covered by `tools/check_workflow_download_pins.sh`.
- Keep workflow scripts defensive: pass untrusted context through environment variables or action inputs, not direct expression interpolation in `run:` blocks.
- Keep analyzer and linter output actionable. Prefer fixing root causes over widening suppressions.
- Compare computed floating-point values with a domain-appropriate tolerance. Exact
  reproducibility checks should compare object bytes; literal sentinel checks and
  IEEE NaN self-comparisons remain allowed. The Semgrep rule covers same-function
  float/double declarations, parameters, arrays, and pointer parameters, including
  quality tools and tests. CodeQL supplies type-aware coverage across headers.
- Strict default Semgrep runs exercise `semgrep/mbelib-neo.c` regression fixtures
  before scanning. `.semgrepignore` keeps tests in directory scans while excluding
  vendored code, build outputs, and intentional rule fixtures.
