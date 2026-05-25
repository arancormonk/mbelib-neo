## Summary

- 

## Review

- [ ] Changes were reviewed against the surrounding module design.
- [ ] Behavior, ownership boundaries, and failure modes were reviewed by a human.
- [ ] Risky or broad changes have a second reviewer.
- [ ] Copied, generated, or bulk-written code has been read, understood, and adapted to mbelib-neo conventions.
- [ ] Non-trivial commits include DCO `Signed-off-by` lines.

## Quality Review

- [ ] New decoder, file input, allocator, threading, dependency, or public API changes include focused tests.
- [ ] Codec, DSP, workflow, dependency, and release changes received extra scrutiny.
- [ ] Static-analysis suppressions include a reason and are limited to the narrowest scope.
- [ ] No new direct bundled-third-party includes, shell execution, broad workflow permissions, or release publishing paths were added without justification.

## Tests And Guardrails

- [ ] `cmake --build --preset dev-debug -j`
- [ ] `ctest --preset dev-debug --output-on-failure`
- [ ] `tools/preflight_ci.sh`
- [ ] `tools/quality_preflight.sh` for broad or high-risk changes
- [ ] `tools/cmake_format_check.sh` for CMake changes
- [ ] `tools/zizmor.sh` for workflow changes
- [ ] `tools/osv_scan.sh` for dependency input changes

## Risk

- Security/dependency impact:
- API/ABI compatibility impact:
- Packaging/release impact:
