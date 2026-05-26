# Supply-Chain Guardrails

mbelib-neo keeps CI dependencies explicit so analyzer results and release artifacts are reproducible.

## GitHub Actions

- Workflow security is checked with `tools/workflow_lint.sh`, Semgrep workflow rules, CodeQL Actions analysis, Gitleaks, OpenSSF Scorecard, and `tools/zizmor.sh`.
- The current zizmor policy in `.github/zizmor.yml` requires every external `uses:` action to be pinned to a full commit SHA.
- For new or refreshed third-party actions in release, packaging, signing, attestation, or upload workflows, resolve the upstream tag to the immutable commit and pin that SHA.
- Public GitHub source checkouts in CI must go through `tools/fetch-pinned-git.sh` with SHAs from `tools/ci-dependency-pins.env`; `tools/check_workflow_git_pins.sh` blocks floating `git clone` and `git ls-remote` usage in workflows and CI helper scripts.
- Release and CI workflows must not execute mutable downloaded helper binaries or mutable container tags. CI fallback container images for Arch Linux, Gitleaks, and OSV-Scanner are pinned by digest in `tools/ci-dependency-pins.env`; `tools/check_workflow_download_pins.sh` enforces these image and download rules.
- Release packaging workflows must generate package SBOMs and attest both the release packages and their SBOMs before publish jobs upload assets. Release publish jobs keep the checksum manifest and manifest attestation as an additional integrity layer.
- Workflow changes that add secrets, write permissions, artifact publication, release upload, or external actions need human review.

## Pinned CI Source Checkouts

`tools/ci-dependency-pins.env` is the checked-in source of truth for CI-only
GitHub source dependencies such as `include-what-you-use`, plus CI fallback
container image digests.

To refresh a pin, resolve the upstream commit outside CI, update the matching
SHA in `tools/ci-dependency-pins.env`, then run
`tools/check_workflow_git_pins.sh`, `tools/check_workflow_download_pins.sh`,
and the affected CI/local build path.

When refreshing a container image pin, resolve the image digest outside CI,
update the matching `tools/ci-dependency-pins.env` value, and update any
workflow `container.image` literal that GitHub Actions cannot read from the env
file.

## AUR SSH Host Keys

The AUR update workflow currently discovers `aur.archlinux.org`'s SSH host key
at runtime with `ssh-keyscan`, then uses `StrictHostKeyChecking=yes` against
that discovered key. This prevents silent host-key changes after setup within
the same job, but it is still trust-on-first-use for each run and remains
vulnerable if the runner's network path is intercepted during key discovery.

The robust control is to store an out-of-band verified AUR ED25519 host key in
this repository or in a GitHub Actions secret, write that exact key to
`known_hosts`, and fail the job if the live host key differs. Refresh that pin
only after verifying the new key through Arch/AUR-controlled channels outside
the CI network path. Until that pin exists, treat AUR SSH host-key discovery as
an accepted residual risk for the AUR publish jobs.

## Vendored Code

Vendored third-party sources live under `src/external/` and are excluded from project-owned analyzer findings where the tool supports path filtering.

- Keep bundled dependency updates small and traceable to an upstream release or commit.
- Preserve upstream license files and attribution.
- Re-run focused codec tests plus `tools/osv_scan.sh` after dependency updates.

## Vulnerability Scanning

- `tools/osv_scan.sh` runs OSV-Scanner source scanning for lockfiles, manifests, and vendored C/C++ dependency fingerprints.
- PR CI runs OSV when dependency inputs change; scheduled guardrails run a full repository scan.
- If OSV reports a vulnerability that is not exploitable in mbelib-neo, add the narrowest possible `osv-scanner.toml` ignore with a reason and expiry date.
