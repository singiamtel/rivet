#! /usr/bin/env bash

RIVET_VERSION=${RIVET_VERSION:-4.0.2}
PYTHIA_VERSION=8312

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BUILDFLAGS+=" --build-arg RIVET_VERSION=$RIVET_VERSION"
BUILDFLAGS+=" --build-arg PYTHIA_VERSION=$PYTHIA_VERSION"

PKG="hepstore/rivet-pythia"
TAGS="${RIVET_VERSION}-${PYTHIA_VERSION} ${RIVET_VERSION}"

source "$SCRIPT_DIR/../docker-build-and-push.sh"
