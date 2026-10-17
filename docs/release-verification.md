<!-- SPDX-License-Identifier: GPL-2.0-or-later -->
# Release Verification

Official releases are published through GitHub Releases:

https://github.com/arancormonk/mbelib-neo/releases

## Release Identifiers

Release tags use the form `vX.Y.Z`. The tag must match the project version
declared by CMake, must be an annotated tag signed by a trusted mbelib-neo
release key, and must point at a commit already contained in `origin/main`, or
the release workflow fails before packaging.

## Release Signing

Release tags are signed with one of the trusted mbelib-neo release keys checked
into the repository:

```sh
gpg --import release-keys/arancormonk-desktop-2026.pgp release-keys/arancormonk-laptop-2026.pgp
```

Verify a release tag:

```sh
git fetch --tags https://github.com/arancormonk/mbelib-neo.git
git tag -v v1.2.7
```

The release workflow pins these exact primary key fingerprints before trusting
the checked-in key files:

```text
5FAF 0C47 C8E1 F95D 33CD 83B1 E42E 43AD D853 F280  (desktop)
3561 9AC5 0DF5 FB9E A053 296E 5C77 FAD4 4C3E 67A7  (laptop)
```

The repository also stores the original upstream mbelib author public key for
historical attribution:

```text
mbelib_Author.pgp
```

That attribution key is:

```text
9E7A 5527 9CDC EBF7 BF1B D772 4F98 E863 EA5E FE2C
```

Do not confuse the upstream attribution key with the current release-signing keys
unless release notes explicitly say otherwise.

GitHub repository rulesets should also protect `v*.*.*` tags, restrict release
tag creation to the maintainer/admin role, block tag deletion and force updates,
and require signed tags where the GitHub ruleset UI supports it. Workflow
verification remains mandatory even when repository rulesets are enabled.

## Asset Integrity

Release assets are distributed over GitHub HTTPS. Each release includes
`SHA256SUMS.txt` with hashes for packaged assets. The tag release workflow
generates SPDX SBOMs and GitHub artifact attestations for the packaged assets,
attests the checksum manifest, then creates or overwrites GitHub Release assets
with the workflow `GITHUB_TOKEN` after packaging succeeds.

Download the release assets and checksum file from the release page, then run:

```sh
sha256sum -c SHA256SUMS.txt
```

On macOS, use:

```sh
shasum -a 256 -c SHA256SUMS.txt
```

Then verify the GitHub artifact attestations. For example:

```sh
gh attestation verify mbelib-neo-<version>-linux-<arch>.tar.gz \
  --repo arancormonk/mbelib-neo
gh attestation verify SHA256SUMS.txt \
  --repo arancormonk/mbelib-neo
```

Review the matching SPDX SBOM file (required for every archive):

```sh
less mbelib-neo-<version>-linux-<arch>.tar.gz.spdx.json
```

## Expected Release Author

The expected release author is the project maintainer, `arancormonk`, using the
public release-signing keys documented above.

## Release Review

Before publication, release packaging must pass the required CI checks for the
tag, including signed-tag validation, build, test, static-analysis, sanitizer,
install/consume, SBOM, attestation, and release validation jobs.

## Complete source and legal review from 2.6.0

A release requires three platform binary archives and one
`mbelib-neo-<version>-source.tar.gz`, each with a matching `.spdx.json`.
`SHA256SUMS.txt` covers all eight files. The Linux packaging job builds source
from the validated tag's tracked tree, includes the original pinned PFFFT
archive, and verifies extraction, build, tests, installation and static/shared
consumers in a network-disabled namespace. Source ordering, ownership, modes
and timestamps are normalized; `SOURCE-IDENTITY.json` records content hashes.

SPDX file checksums and the package verification code retain SHA-1 for
[SPDX 2.3 compatibility](https://spdx.github.io/spdx-spec/v2.3.1/package-information/#79-package-verification-code-field).
They are inventory identifiers, not the release's cryptographic integrity
check. Archive, SBOM and source-identity verification uses SHA-256, together
with signed tags and attestations. The single `spdx_sha1` helper has a narrow,
documented scanner exception; the SHA-1 security rule remains enabled elsewhere.
Source staging compares archived contents directly with the Git blob bytes.

Verify both the source and its SBOM before building, for example:

```sh
gh attestation verify mbelib-neo-2.6.0-source.tar.gz --repo arancormonk/mbelib-neo
gh attestation verify mbelib-neo-2.6.0-source.tar.gz.spdx.json --repo arancormonk/mbelib-neo
```

See [source-build.md](source-build.md), the [component inventory](../legal/components.json),
[licensing status](../legal/status.md) and [patent review](../legal/patents.md).
Release preparation must update their dated research when codec behavior or
evidence changes. The disclosed ETSI and patent questions remain open; passing
these engineering gates is not legal clearance.

Recovery requires four archives and their signed SBOMs for tags >=2.6.0;
older validated tags retain their original three-package format. Recovery
verifies every archive's provenance, exact signed SBOM predicate, and (for new
releases) the SBOM file's provenance before publication. Checksums are rebuilt
only after collection, and the manifest is attested. Existing historical tags
and published release assets are not changed by this remediation.

### Rebuilding packages for an unchanged tag

If the original tag passed every analysis, native, optimized, sanitizer and
consumer check but packaging was skipped or failed, dispatch `CI` from `main`
with `release_tag`, that original `release_run_id`, and `rebuild_packages=true`.
This option supports the component-aware release format from 2.6.0 onward.
This path verifies the signed tag and original checks, reruns the current
workflow gates, and checks out the tag's exact commit for all three package
builds and their tests. It never creates or moves a tag.

Rebuild attestations identify the reviewed workflow commit on `main`; the
signed SBOMs and corresponding-source identity identify the original release
commit. Publication verifies both identities, all package/SBOM attestations,
and exact signed SBOM contents. Release notes link the rebuild run and record
both commits. Using current verification tooling, pass `--builder-run-id RUN`
to `tools/ci_verify_release_recovery.py --assets DIRECTORY` with `RELEASE_TAG`
and `RELEASE_COMMIT` set to the unchanged release. A rebuild whose publication
fails can rerun its failed publish job; ordinary publication-only recovery
continues to require packages from an original tag build.

### Recovering the stable AUR update

Release recovery also runs the stable AUR job using the validated tag's version,
source archive and dependency pins. If GitHub assets are already published and
only AUR needs updating, dispatch `CI` from `main` with `release_tag`, the original
successful tag-check `release_run_id`, and `update_aur_only=true`. This mode
requires the same signed tag and original checks as package rebuilding, reruns
the current native checks, and skips package creation and GitHub publication.
Do not combine it with `rebuild_packages`. Check the AUR job's steps: missing
credentials or AUR maintenance can leave an update skipped despite a green job.

Manual CI and push CI use separate workflow concurrency groups, so starting
release recovery does not cancel the normal `main` push checks.
