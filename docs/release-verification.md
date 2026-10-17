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
