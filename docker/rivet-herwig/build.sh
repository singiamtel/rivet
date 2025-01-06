#! /usr/bin/env bash

RIVET_VERSION=${RIVET_VERSION:-4.0.2}
# NOTE: have to use branches until Rivet4 is supported in an H7 release
THEPEG_VERSION=release-2-3 #2.3.0
HERWIG_VERSION=herwig-7-3 #7.3.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BUILDFLAGS+=" --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILDFLAGS+=" --build-arg THEPEG_VERSION=${THEPEG_VERSION}"
BUILDFLAGS+=" --build-arg HERWIG_VERSION=${HERWIG_VERSION}"

PKG="hepstore/rivet-herwig"
TAGS="${RIVET_VERSION}-${HERWIG_VERSION} ${RIVET_VERSION}"

source "$SCRIPT_DIR/../docker-build-and-push.sh"
