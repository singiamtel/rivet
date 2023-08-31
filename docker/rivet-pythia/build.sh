#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.8
PYTHIA_VERSION=8309

#PLATFLAGS="--platform linux/amd64,linux/arm64"
#BUILD="docker buildx build -f Dockerfile $PLATFLAGS $DOCKERFLAGS"
#if [[ -n "$PLATFLAGS" && "$PUSH" = 1 ]]; then BUILD="$BUILD --push"; fi

BUILD="docker build . -f Dockerfile $DOCKERFLAGS"

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg PYTHIA_VERSION=${PYTHIA_VERSION}"

BUILD="$BUILD -t hepstore/rivet-pythia:${RIVET_VERSION}-${PYTHIA_VERSION}"
BUILD="$BUILD -t hepstore/rivet-pythia:${RIVET_VERSION}"
if [[ "$LATEST" = 1 ]]; then BUILD="$BUILD -t hepstore/rivet-pythia:latest"; fi

echo "Building: $BUILD"
$BUILD

if [[ "$PUSH" = 1 ]]; then
    docker push -a hepstore/rivet-pythia
fi
