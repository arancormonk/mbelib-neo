# Supply-Chain Guardrails

mbelib-neo keeps CI dependencies explicit so analyzer results and release artifacts are reproducible.

## GitHub Actions

- Workflow security is checked with `tools/workflow_lint.sh`, Semgrep workflow rules, CodeQL Actions analysis, Gitleaks, OpenSSF Scorecard, and `tools/zizmor.sh`.
- The current zizmor policy in `.github/zizmor.yml` requires every external `uses:` action to be pinned to a full commit SHA.
- For new or refreshed third-party actions in release, packaging, signing, attestation, or upload workflows, resolve the upstream tag to the immutable commit and pin that SHA.
- Public GitHub source checkouts in CI must go through `tools/fetch-pinned-git.sh` with SHAs from `tools/ci-dependency-pins.env`; `tools/check_workflow_git_pins.sh` blocks floating `git clone` and `git ls-remote` usage in workflows and CI helper scripts.
- Workflow changes that add secrets, write permissions, artifact publication, release upload, or external actions need human review.

## Pinned CI Source Checkouts

`tools/ci-dependency-pins.env` is the checked-in source of truth for CI-only
GitHub source dependencies such as `include-what-you-use`.

To refresh a pin, resolve the upstream commit outside CI, update the matching
SHA in `tools/ci-dependency-pins.env`, then run
`tools/check_workflow_git_pins.sh` and the affected CI/local build path.

## Vendored Code

Vendored third-party sources live under `src/external/` and are excluded from project-owned analyzer findings where the tool supports path filtering.

- Keep bundled dependency updates small and traceable to an upstream release or commit.
- Preserve upstream license files and attribution.
- Re-run focused codec tests plus `tools/osv_scan.sh` after dependency updates.

## Vulnerability Scanning

- `tools/osv_scan.sh` runs OSV-Scanner source scanning for lockfiles, manifests, and vendored C/C++ dependency fingerprints.
- PR CI runs OSV when dependency inputs change; scheduled guardrails run a full repository scan.
- If OSV reports a vulnerability that is not exploitable in mbelib-neo, add the narrowest possible `osv-scanner.toml` ignore with a reason and expiry date.
