# Roadmap

This roadmap covers the next year of mbelib-neo maintenance and development.
It is a planning document, not a promise that every item will ship.

## Current Scope

mbelib-neo provides C library primitives for IMBE/AMBE vocoder frame handling,
error correction, synthesis, packaging, tests, and developer tooling.

The project intends to remain:

- a C library with installable headers
- compatible with CMake and pkg-config consumers
- focused on AMBE/IMBE codec primitives and related synthesis behavior
- testable through CTest and GitHub Actions
- conservative about new runtime dependencies

## 2026 Priorities

- Stabilize the 2.x API and ABI after the v2.0.0 ABI break.
- Keep release packaging working across Linux, macOS, and Windows.
- Improve public API documentation and generated Doxygen output.
- Expand regression tests for soft-decision decode paths, error handling, and
  synthesis edge cases.
- Keep static analysis, fuzzing, dependency review, and secret scanning active
  in CI.
- Improve release verification documentation and release artifact provenance.
- Reduce project continuity risk by recruiting additional reviewers or backup
  maintainers.

## Out of Scope for the Next Year

- Implementing patented codec encoders or patent-clearing downstream use cases.
- Adding network services, daemon modes, or user account management.
- Adding non-FLOSS runtime dependencies.
- Replacing CMake as the primary build system.
- Guaranteeing bit-identical floating-point audio across all compilers,
  architectures, and optimization settings.
- Providing legal advice about patent licensing or jurisdiction-specific use.

## Compatibility Policy

- Patch releases should preserve source and binary compatibility within the same
  major version.
- Minor releases should preserve ABI unless release notes explicitly say
  otherwise.
- Major releases may break ABI and must document migration notes.

## Maintenance Policy

Security fixes are handled on the default branch and latest tagged release line.
Older releases may receive fixes when the maintainer determines that the impact
and user base justify the work.
