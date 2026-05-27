# Release Verification

Official releases are published through GitHub Releases:

https://github.com/arancormonk/mbelib-neo/releases

## Release Identifiers

Release tags use the form `vX.Y.Z`. The tag must match the project version
declared by CMake, must be an annotated tag signed by the trusted mbelib-neo
release key, and must point at a commit already contained in `origin/main`, or
the release workflow fails before packaging.

## Release Signing

Release tags are signed with the trusted mbelib-neo release key checked into the
repository:

```sh
gpg --import release-keys/arancormonk-2026.pgp
```

Verify a release tag:

```sh
git fetch --tags https://github.com/arancormonk/mbelib-neo.git
git tag -v v1.2.7
```

The release workflow pins this exact primary key fingerprint before trusting the
checked-in key file:

```text
5FAF 0C47 C8E1 F95D 33CD 83B1 E42E 43AD D853 F280
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

Do not confuse the upstream attribution key with the current release-signing key
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

Review the matching SPDX SBOM file when present:

```sh
less mbelib-neo-<version>-linux-<arch>.tar.gz.spdx.json
```

## Expected Release Author

The expected release author is the project maintainer, `arancormonk`, using the
public release-signing key documented above.

## Release Review

Before publication, release packaging must pass the required CI checks for the
tag, including signed-tag validation, build, test, static-analysis, sanitizer,
install/consume, SBOM, attestation, and release validation jobs.
