# Testing Policy

mbelib-neo uses automated tests, static analysis, sanitizers, and fuzzing to
reduce regression and security risk.

## Test Suites

The CTest suite includes:

- API/version/result helper checks
- ECC tests for hard and soft Golay/Hamming paths
- noise determinism and frame-state determinism checks
- parameter and synthesis behavior checks
- float-to-int16 conversion parity checks
- golden PCM hash regression checks
- SIMD architecture detection checks

Run the default test suite with:

```sh
cmake --preset dev-debug
cmake --build --preset dev-debug -j
ctest --preset dev-debug --output-on-failure
```

## Continuous Integration

GitHub Actions runs tests and quality checks on pull requests and pushes to the
primary branch. Required checks include cross-platform builds, sanitizer tests,
static analysis, workflow linting, dependency review, secret scanning, OSV
scanning, and install/consume checks.

## Regression Test Requirement

At least 50% of bugs fixed in the last six months should include regression
tests. A pull request that fixes a bug should add a regression test unless:

- the behavior cannot be reproduced reliably in automation
- the fix is entirely documentation or packaging metadata
- a better guardrail exists, such as a static-analysis rule or workflow check

When no regression test is added for a bug fix, the pull request must explain
why.

## Major Functionality Test Requirement

Major new functionality must add or update automated tests. Major functionality
includes:

- public API changes
- codec or DSP behavior changes
- external input handling changes
- dependency changes that affect compiled code
- security-sensitive workflow or release changes
- installation and packaging behavior changes

## Coverage Target

The project target is at least 80% statement coverage for project-owned source
files when measured with a FLOSS C coverage tool. Vendored code under
`src/external/` is excluded from project-owned coverage accounting.

Coverage evidence should be generated from an instrumented Debug build and
recorded in a pull request, issue, or CI artifact before claiming coverage-based
badge criteria.

## Dynamic Analysis

For memory-safety-sensitive C changes, run sanitizer tests:

```sh
cmake --preset asan-ubsan-debug
cmake --build --preset asan-ubsan-debug -j
ctest --preset asan-ubsan-debug --output-on-failure
```

Frame-processing paths are also covered by ClusterFuzzLite PR fuzzing with
AddressSanitizer, including fixed-size hard/soft frame decode paths and
parameter synthesis paths.
