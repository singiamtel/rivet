#! /usr/bin/env bash

RIVET_VERSION=${RIVET_VERSION:-4.0.2}
HERWIG_VERSION=7.3.0
THEPEG_VERSION=2.3.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BUILDFLAGS+=" --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILDFLAGS+=" --build-arg THEPEG_VERSION=${THEPEG_VERSION}"
BUILDFLAGS+=" --build-arg HERWIG_VERSION=${HERWIG_VERSION}"

# TODO: single-source the heavily shared logic below
PKG="hepstore/rivet-herwig"
TAGS="${RIVET_VERSION}-${HERWIG_VERSION} ${RIVET_VERSION}"
test "$LATEST" = 1 && TAGS+=" latest"
for tag in TAGS; do TAGFLAGS="$TAGFLAGS -t $PKG:$tag"; done

dx_build $BUILDFLAGS $TAGFLAGS

if [[ "$PUSH" = 1 ]]; then
    for tag in $TAGS; do
        xdocker push $PKG:$tag
    done
fi
