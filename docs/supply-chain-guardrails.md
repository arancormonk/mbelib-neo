# Supply-Chain Guardrails

mbelib-neo keeps CI dependencies explicit so analyzer results and release artifacts are reproducible.

## GitHub Actions

- Workflow security is checked with `tools/workflow_lint.sh`, Semgrep workflow rules, CodeQL Actions analysis, Gitleaks, OpenSSF Scorecard, and `tools/zizmor.sh`.
- The current zizmor policy in `.github/zizmor.yml` requires every external `uses:` action to be pinned to a full commit SHA.
- For new or refreshed third-party actions in release, packaging, signing, attestation, or upload workflows, resolve the upstream tag to the immutable commit and pin that SHA.
- Workflow changes that add secrets, write permissions, artifact publication, release upload, or external actions need human review.

## Vendored Code

Vendored third-party sources live under `src/external/` and are excluded from project-owned analyzer findings where the tool supports path filtering.

- Keep bundled dependency updates small and traceable to an upstream release or commit.
- Preserve upstream license files and attribution.
- Re-run focused codec tests plus `tools/osv_scan.sh` after dependency updates.

## Vulnerability Scanning

- `tools/osv_scan.sh` runs OSV-Scanner source scanning for lockfiles, manifests, and vendored C/C++ dependency fingerprints.
- PR CI runs OSV when dependency inputs change; scheduled guardrails run a full repository scan.
- If OSV reports a vulnerability that is not exploitable in mbelib-neo, add the narrowest possible `osv-scanner.toml` ignore with a reason and expiry date.
