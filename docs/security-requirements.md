# Security Requirements

mbelib-neo is a native C library. Its main security goal is to process
potentially malformed codec frames without memory corruption, credential
exposure, or unexpected process execution.

## What Users Can Expect

- Public APIs should not read or write outside documented caller-provided
  buffers.
- Processing functions should be reentrant when each stream uses separate
  `mbe_parms` state.
- The library should not execute shells, spawn processes, open network
  connections, or access credentials during normal audio processing.
- The library should not store authentication credentials or private keys.
- The library should not implement network protocols or TLS.
- Build and release workflows should use least-privilege permissions and avoid
  exposing release credentials to pull-request code.
- Security vulnerabilities should be handled through the private reporting
  process in `SECURITY.md`.

## What Users Should Not Expect

- mbelib-neo does not provide patent licensing advice or patent clearance.
- mbelib-neo is not a sandbox. Applications that process untrusted captures
  should still use process isolation, privilege separation, and normal media
  hardening practices.
- Floating-point synthesis output is not guaranteed to be bit-identical across
  all compiler, architecture, SIMD, and optimization combinations.
- The library does not authenticate frame contents or determine whether an
  upstream radio/capture source is trustworthy.

## Threat Model

Relevant threats include:

- malformed frame data triggering buffer overflows or undefined behavior
- crafted inputs driving excessive CPU work or unstable state transitions
- mistakes in SIMD or architecture-specific code paths
- accidental release of generated binaries, private captures, credentials, or
  signing material
- compromised CI/CD dependencies or workflow changes
- vulnerabilities in vendored third-party code

Out-of-scope threats include:

- physical radio-layer attacks outside the library process
- patent enforcement or licensing disputes
- malicious applications intentionally misusing the library API
- compromise of GitHub or package-hosting infrastructure outside the project's
  control

## Trust Boundaries

- Caller-provided frame, data, state, and output buffers are untrusted until
  checked by the API implementation.
- Vendored code under `src/external/` is treated as external code and reviewed
  separately from project-owned sources.
- CI workflow files are security-sensitive because they control analyzer,
  packaging, and release behavior.
- Release signing keys, GitHub Actions secrets, and AUR credentials are
  sensitive resources and must not be exposed to untrusted pull-request code.

## Secure Design Controls

The project applies the following controls:

- explicit CMake source lists instead of broad source globbing
- public API declarations separated from private implementation headers
- test coverage for API helpers, ECC behavior, frame-state determinism, golden
  PCM regression behavior, parameter behavior, float conversion parity, and SIMD
  architecture detection
- sanitizer CI for AddressSanitizer and UndefinedBehaviorSanitizer
- ClusterFuzzLite PR fuzzing with AddressSanitizer for frame-processing paths
- default-on Release-like compiler/linker hardening for supported Clang/GCC
  targets, with Linux release verification in CI
- static analysis through CodeQL, Semgrep, clang-tidy, cppcheck, scan-build,
  GCC `-fanalyzer`, and include-what-you-use
- Gitleaks and GitHub secret scanning for credential leakage
- OSV-Scanner and dependency review for dependency vulnerabilities
- pinned GitHub Actions, pinned CI source checkouts, and zizmor workflow
  security checks
- branch protection requiring status checks before merge, including repository
  guardrails

## Solo Maintainer Operating Mode

mbelib-neo is allowed to operate as a solo-maintainer project, but the absence
of a routine second reviewer is treated as residual risk, not as implicit
approval. For security-sensitive, workflow, release, parser, dependency, public
API, codec, DSP, or broad runtime changes, the maintainer should keep changes
narrowly scoped, document risk in the pull request or release notes, run the
relevant guardrails, and seek outside review when practical.

## Cryptography and Credentials

The library does not implement cryptographic protocols, password storage, TLS,
or credential management. Criteria about cryptographic algorithm agility,
credential agility, TLS certificate verification, and network protocol security
are not applicable to the library runtime.

Project release and repository operations do use cryptographic mechanisms:

- signed Git tags for releases
- GitHub HTTPS for repository and release distribution
- SHA256 release checksums
- SPDX SBOMs and GitHub artifact attestations for release packages
- GitHub artifact attestations for checksum manifests
- verification instructions in `docs/release-verification.md`

## Vulnerability Handling

Vulnerability reports are handled through `SECURITY.md`. Confirmed
vulnerabilities should result in a fix, release note or advisory, affected
version information, and reporter credit unless anonymity is requested.
