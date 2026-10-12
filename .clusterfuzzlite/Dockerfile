# SPDX-License-Identifier: GPL-2.0-or-later
FROM gcr.io/oss-fuzz-base/base-builder@sha256:0d2576e650679ea4e3dd99cf371badafff54412c96a012a11ac61ba3f06ae5bf

COPY . $SRC/mbelib-neo
WORKDIR $SRC/mbelib-neo
RUN . tools/ci-dependency-pins.env && \
    tools/fetch-pinned-git.sh https://github.com/arancormonk/pffft.git "$PFFFT_SHA" "$SRC/pffft"
COPY .clusterfuzzlite/build.sh $SRC/
