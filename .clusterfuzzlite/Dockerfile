# SPDX-License-Identifier: GPL-2.0-or-later
FROM gcr.io/oss-fuzz-base/base-builder@sha256:1f9cf74bc3d7029b528494542f9d9851b256dceed38d3f420989434dc5ae0760

COPY . $SRC/mbelib-neo
WORKDIR $SRC/mbelib-neo
RUN . tools/ci-dependency-pins.env && \
    tools/fetch-pinned-git.sh https://github.com/arancormonk/pffft.git "$PFFFT_SHA" "$SRC/pffft"
COPY .clusterfuzzlite/build.sh $SRC/
