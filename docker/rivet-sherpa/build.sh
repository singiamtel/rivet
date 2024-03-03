#! /usr/bin/env bash

set -e

RIVET_VERSION=${RIVET_VERSION:-3.1.10}
SHERPA_VERSION=2.2.15

#PLATFLAGS="--platform linux/amd64,linux/arm64"
#BUILD="docker buildx build -f Dockerfile $PLATFLAGS $DOCKERFLAGS ."
#if [[ -n "$PLATFLAGS" && "$PUSH" = 1 ]]; then BUILD="$BUILD --push"; fi

BUILD="docker build . -f Dockerfile --progress=plain $DOCKERFLAGS"

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg SHERPA_VERSION=${SHERPA_VERSION}"

BUILD="$BUILD -t hepstore/rivet-sherpa:${RIVET_VERSION}-${SHERPA_VERSION}"
BUILD="$BUILD -t hepstore/rivet-sherpa:${RIVET_VERSION}"
if [[ "$LATEST" = 1 ]]; then BUILD="$BUILD -t hepstore/rivet-sherpa:latest"; fi

echo "Building: $BUILD"
$BUILD

if [[ "$PUSH" = 1 ]]; then
    docker push -a hepstore/rivet-sherpa
fi
