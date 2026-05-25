# Dependency Management

mbelib-neo keeps runtime dependencies small and reviewable.

## Compiled Dependencies

The project has one vendored compiled third-party component:

- PFFFT/FFTPACK under `src/external/pffft/`

The vendored code retains upstream license files and attribution. Updates should
be small, traceable to upstream, and reviewed separately from unrelated changes.

System dependencies used by consumers are:

- C99 compiler
- CMake 3.20 or newer
- platform math library where required by the toolchain

## Tooling Dependencies

CI and local quality tools are tracked through:

- `.github/requirements/*.in`
- `.github/requirements/*.txt`
- `.github/dependabot.yml`
- `.github/workflows/*.yml`
- `tools/*.sh`

Hashed Python requirements are used where Python tooling is installed in CI.
GitHub Actions are pinned to immutable commit SHAs by policy.

## Monitoring

The project monitors dependencies through:

- Dependabot for GitHub Actions and Python requirements
- GitHub dependency review on pull requests
- OSV-Scanner through `tools/osv_scan.sh`
- scheduled guardrail CI
- OpenSSF Scorecard

If a vulnerability is reported in a dependency:

1. Determine whether the vulnerable code is present and reachable.
2. Update or patch the dependency when exploitable.
3. If not exploitable, document the reason in the narrowest available
   suppression or advisory note.
4. Add tests or checks when the issue could regress.

## Update Policy

Dependency updates should:

- preserve license notices
- avoid unrelated refactors
- update documentation when dependency behavior or requirements change
- run focused codec tests and `tools/osv_scan.sh`
- receive human review when they affect compiled code, workflows, packaging, or
  release behavior
