# OpenSSF OSPS Baseline Evidence

This document records mbelib-neo's evidence for OpenSSF OSPS Baseline Level 1.

## Repository and Access Controls

- Sensitive repository resources are hosted on GitHub. Maintainers with access
  to those resources must keep multi-factor authentication or passkeys enabled
  on their GitHub accounts.
- New collaborators are added through GitHub repository permissions. Access is
  intentionally granted by the maintainer, and elevated access is not granted by
  default.
- The primary branch is `main`. GitHub branch protection is required for `main`,
  including required status checks and repository guardrails before merge. Direct
  commits to `main` are blocked by the protected branch workflow.
- Deleting the protected `main` branch is disabled in GitHub branch protection.

## Build, Release, and Secrets Controls

- CI/CD workflows must treat GitHub event metadata, branch names, issue titles,
  pull request fields, and other user-controlled values as untrusted input.
  Workflow scripts must pass untrusted values through action inputs or
  environment variables instead of directly interpolating them into shell
  scripts.
- Pull request workflows use least-privilege `GITHUB_TOKEN` permissions. The
  default GitHub Actions workflow permission for the repository is read-only.
  Release publication elevates `contents: write` only in trusted upstream
  version-tag paths and publishes with the workflow `GITHUB_TOKEN`.
- Official project channels and distribution channels are GitHub HTTPS URLs.
- The project prevents accidental credential storage through contributor policy,
  GitHub secret scanning, Gitleaks CI, and review of workflow and release
  changes.

## Documentation and Project Scope

- `README.md` documents basic usage, build, installation, examples, release
  flow, and public API notes.
- `docs/issue-reporting.md` documents defect reporting through GitHub Issues.
- GitHub Issues, pull requests, and Discussions are the public mechanisms for
  proposed changes and usage obstacles.
- `CONTRIBUTING.md` documents the contribution process and acceptable
  contribution requirements.
- This project has one authoritative source repository:
  `https://github.com/arancormonk/mbelib-neo`.

## Licensing and Repository Contents

- Source code is licensed under GPL-2.0-or-later. Third-party license notices
  are retained under `LICENSES/` and `src/external/pffft/`.
- Release install trees include the project license, GPL text, ISC text, and
  upstream attribution.
- The public Git repository records changes, authors, and timestamps.
- Generated executable artifacts, compiled binaries, release archives, private
  captures, credentials, keys, and machine-specific configuration must not be
  committed.

## Dependency and Vulnerability Contacts

- Direct compiled and tooling dependencies are documented in
  `docs/dependencies.md`, CMake metadata, `.github/requirements/`, workflow
  files, and Dependabot configuration.
- Security contacts and private vulnerability reporting are documented in
  `SECURITY.md`.
