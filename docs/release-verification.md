# Release Verification

Official releases are published through GitHub Releases:

https://github.com/arancormonk/mbelib-neo/releases

## Release Identifiers

Release tags use the form `vX.Y.Z`. The tag must match the project version
declared by CMake, or the release workflow fails before packaging.

## Release Signing

Release tags are signed with the maintainer release key. Obtain the maintainer's
public keys from GitHub:

```sh
curl -fsSL https://github.com/arancormonk.gpg | gpg --import
```

Verify a release tag:

```sh
git fetch --tags https://github.com/arancormonk/mbelib-neo.git
git tag -v v1.2.7
```

The expected release-signing key fingerprint for v1.2.7 is:

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

## Asset Integrity

Release assets are distributed over GitHub HTTPS. Each release includes
`SHA256SUMS.txt` with hashes for packaged assets. The tag release workflow
creates or overwrites GitHub Release assets with the workflow `GITHUB_TOKEN`
after packaging succeeds.

Download the release assets and checksum file from the release page, then verify
the checksum file signature when `SHA256SUMS.txt.asc` is present:

```sh
gpg --verify SHA256SUMS.txt.asc SHA256SUMS.txt
```

Then run:

```sh
sha256sum -c SHA256SUMS.txt
```

On macOS, use:

```sh
shasum -a 256 -c SHA256SUMS.txt
```

## Expected Release Author

The expected release author is the project maintainer, `arancormonk`, using the
public release-signing key documented above.

## Release Review

Before publication, release packaging must pass the required CI checks for the
tag, including build, test, static-analysis, sanitizer, install/consume, and
release validation jobs.
