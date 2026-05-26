# Contributing to mbelib-neo

mbelib-neo accepts changes through GitHub issues and pull requests. Keep changes
small, focused, and reviewable.

## Contribution Process

1. Open or reference an issue for user-visible behavior changes, API changes,
   dependency changes, release changes, and security-sensitive work.
2. Create a branch from `main`.
3. Make the smallest coherent change that solves the problem.
4. Add or update tests and documentation when behavior, public API, build
   behavior, release behavior, or security posture changes.
5. Run the relevant local checks before opening a pull request.
6. Open a pull request against `main` and complete the pull request template.
7. Address review comments and wait for required CI checks to pass before merge.

Security vulnerabilities must not be reported through public issues or pull
requests. Use the private reporting process in `SECURITY.md`.

## Legal Certification

Non-trivial contributions must use the Developer Certificate of Origin 1.1
sign-off. Add this line to each commit message:

```text
Signed-off-by: Your Name <your.email@example.com>
```

The sign-off means you certify that you have the right to submit the work under
the project's license. The DCO text is published at:

https://developercertificate.org/

Use `git commit -s` to add the sign-off automatically.

## Acceptable Contribution Requirements

A contribution is acceptable only when it meets the requirements below or the
pull request explains why an exception is appropriate.

- The change is compatible with GPL-2.0-or-later and does not introduce license
  conflicts.
- Project-maintained source files include an SPDX license identifier near the
  top of the file.
- New or changed public API is documented in `include/mbelib-neo/` comments and
  in user-facing documentation when needed.
- Major functionality changes include focused automated tests.
- Bug fixes include regression tests when practical.
- Generated build outputs, compiled binaries, private captures, credentials,
  keys, and machine-specific configuration are not committed.
- Vendored third-party code is placed under `src/external/`, keeps upstream
  license notices, and is documented in dependency records.
- Workflow, dependency, packaging, release, codec, DSP, public API, and security
  changes receive extra review when practical, or document solo-maintainer
  residual risk in the pull request.
- Static-analysis suppressions are narrow and explain why the local exception is
  acceptable.

## Coding Standards

The project uses C99 and CMake.

- C code follows the repository `.clang-format` configuration, based on LLVM
  style with 4-space indentation and a 120-column limit.
- CMake files are formatted with `gersemi` through
  `tools/cmake_format_check.sh`.
- Public symbols use the `mbe_` prefix unless compatibility requires otherwise.
- Keep internal symbols `static` where practical.
- Avoid direct inclusion of bundled third-party headers outside approved
  integration points.
- Do not spawn shells or external processes from project-owned C code without
  explicit design review.

Run formatting with:

```sh
tools/format.sh
```

## Tests and Local Checks

Run the smallest useful check set for the change, then broaden it for risky
changes.

```sh
cmake --preset dev-debug
cmake --build --preset dev-debug -j
ctest --preset dev-debug --output-on-failure
```

Normal pre-push check:

```sh
tools/preflight_ci.sh
```

Broad or high-risk changes:

```sh
tools/quality_preflight.sh
```

Additional focused checks:

- CMake changes: `tools/cmake_format_check.sh`
- Workflow changes: `tools/workflow_lint.sh` and `tools/zizmor.sh`
- Dependency input changes: `tools/osv_scan.sh`
- Sanitizer-sensitive code: configure and build `asan-ubsan-debug`, then run
  `ctest --preset asan-ubsan-debug --output-on-failure`

## Test Policy

Major new functionality must add or update automated tests. For this project,
"major" includes public API changes, codec behavior changes, DSP changes,
dependency updates that affect compiled code, security-sensitive behavior,
workflow changes that gate releases, and packaging/install changes.

The preferred test style is a focused regression test under `tests/` that can
run through CTest. If a behavior cannot be tested directly, the pull request
must explain the reason and include the best practical substitute, such as a
golden hash, build/install consumer check, analyzer rule, or fuzz target update.

## Review Requirements

Every pull request is reviewed for:

- correctness and compatibility with surrounding module design
- security impact, including input handling and dependency changes
- API and ABI compatibility
- licensing and attribution
- test coverage and documentation coverage
- packaging and release impact

Risky or broad changes should receive a second human review before merge.
