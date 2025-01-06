#! /usr/bin/env bash

RIVET_VERSION=${RIVET_VERSION:-4.0.2}
SHERPA_VERSION=2.2.16
# TODO: update to Sherpa 3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BUILDFLAGS+=" --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILDFLAGS+=" --build-arg SHERPA_VERSION=${SHERPA_VERSION}"

PKG="hepstore/rivet-sherpa"
TAGS="${RIVET_VERSION}-${SHERPA_VERSION} ${RIVET_VERSION}"

source "$SCRIPT_DIR/../docker-build-and-push.sh"
