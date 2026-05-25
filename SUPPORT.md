# Support Policy

## Supported Versions

Security fixes are handled on the default branch and the latest tagged release
line.

Older tagged releases are generally superseded by the latest release. If an older
release has a significant security impact and users cannot reasonably upgrade,
the maintainer may publish a backport or mitigation note.

## Upgrade Path

Users should upgrade to the newest tagged release unless release notes describe a
reason not to.

For ABI-breaking releases:

- the major version changes
- release notes describe the compatibility impact
- public headers and README API notes describe the new behavior
- older tags remain available in Git

## Getting Support

- Defects and enhancement requests: open a GitHub issue.
- Usage questions: use GitHub Discussions or an issue when discussion needs to
  be tracked.
- Security vulnerabilities: use the private process in `SECURITY.md`.

## Security Update Expectations

Confirmed vulnerabilities are fixed on the default branch first. The release
notes or advisory will identify affected versions, available fixed versions, and
known mitigations when practical.
