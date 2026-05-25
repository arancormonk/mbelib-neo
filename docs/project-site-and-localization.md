# Project Site, Accessibility, and Localization

## Project Sites

The official project site, repository, issue tracker, discussions, and release
downloads are hosted on GitHub:

https://github.com/arancormonk/mbelib-neo

The project does not operate a separate website, user database, or password
store. Authentication for repository and discussion access is handled by GitHub.

## Password Security

mbelib-neo does not store passwords for external users. Project sites that allow
login are GitHub-hosted and rely on GitHub's authentication and password storage
controls.

The library runtime also does not store user passwords, authentication tokens,
or private cryptographic keys.

## Accessibility

Project documentation is maintained as Markdown text in the repository so it can
be read through GitHub, cloned locally, rendered by common Markdown tools, or
converted to other formats.

When adding documentation:

- use descriptive headings
- keep tables simple
- avoid conveying important information by color alone
- provide text alternatives for meaningful images if images are added
- keep generated documentation navigable through headings and link text

mbelib-neo is a C library and does not provide a graphical user interface.

## Internationalization

The library does not generate end-user prose, sort human-language text, process
locale-sensitive text, or provide a graphical interface. Runtime output is
limited to version/status strings and caller-controlled audio buffers.

For that reason, software internationalization is currently not applicable to the
library runtime. Documentation and issue discussion are maintained in English.
