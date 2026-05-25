# Governance

mbelib-neo is maintained as a small, maintainer-led open source project.

## Decision Model

The project uses maintainer consensus where possible. For routine changes, a
maintainer may merge a pull request after the required checks pass and review
concerns are addressed.

For broad or risky changes, the maintainer should wait for additional review or
public discussion before merging. Risky changes include:

- public API or ABI changes
- codec, DSP, or external-input handling changes
- dependency, vendored-code, or licensing changes
- CI/CD, release, signing, or packaging changes
- security policy or vulnerability-handling changes

If consensus is not reached, the maintainer makes the final decision and should
document the rationale in the issue or pull request.

## Roles and Responsibilities

### Maintainer

Current maintainer: `@arancormonk`

Responsibilities:

- triage issues and pull requests
- review and merge changes
- maintain CI/CD and release automation
- publish releases
- manage repository settings and branch protection
- handle private vulnerability reports
- maintain project documentation and roadmap
- decide whether proposed work is in scope

### Security Contact

Current security contact: `@arancormonk` through GitHub private vulnerability
reporting.

Responsibilities:

- acknowledge private vulnerability reports
- coordinate fixes and disclosure
- request additional information when needed
- publish advisories or release notes for confirmed vulnerabilities
- credit reporters unless anonymity is requested

### Contributor

Responsibilities:

- follow `CONTRIBUTING.md`
- certify non-trivial contributions with DCO sign-off
- add tests and documentation for meaningful behavior changes
- keep pull requests focused and reviewable
- avoid public disclosure of private security reports

## Sensitive Resources

Sensitive project resources include:

- GitHub repository administration
- GitHub private vulnerability reports and draft advisories
- release signing keys
- package publication credentials
- GitHub Actions secrets
- AUR publication credentials

Access to sensitive resources is limited to the maintainer and any explicitly
approved backup or successor. New access must be intentionally granted and
reviewed before elevation.

## Access Continuity

The project is designed to continue with minimal interruption if the maintainer
is unavailable.

The maintainer must keep an offline emergency access packet for a trusted
successor. The packet should include:

- repository recovery instructions
- release and package publication instructions
- location of release signing public keys
- instructions for rotating or revoking secrets
- authority to appoint or transfer maintainership if the maintainer is
  confirmed unavailable

The trusted successor does not need routine project access, but must be able to
help recover enough access to create and close issues, accept proposed changes,
and publish a release within one week after unavailability is confirmed.

Do not publish private keys, recovery codes, personal emergency contacts, or
secret locations in the repository.

## Bus Factor

The current public bus factor is one maintainer. This is intentionally documented
so users can evaluate continuity risk. The roadmap includes reducing this risk by
recruiting reviewers and potential backup maintainers.
