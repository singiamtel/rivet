#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.8
HERWIG_VERSION=7.2.3
THEPEG_VERSION=2.2.3

#PLATFLAGS="--platform linux/amd64,linux/arm64"
#BUILD="docker buildx build -f Dockerfile $PLATFLAGS $DOCKERFLAGS ."
#if [[ -n "$PLATFLAGS" && "$PUSH" = 1 ]]; then BUILD="$BUILD --push"; fi

BUILD="docker build . -f Dockerfile $DOCKERFLAGS"

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg THEPEG_VERSION=${THEPEG_VERSION}"
BUILD="$BUILD --build-arg HERWIG_VERSION=${HERWIG_VERSION}"

BUILD="$BUILD -t hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}"
BUILD="$BUILD -t hepstore/rivet-herwig:${RIVET_VERSION}"
if [[ "$LATEST" = 1 ]]; then BUILD="$BUILD -t hepstore/rivet-herwig:latest"; fi

echo "Building: $BUILD"
$BUILD

if [[ "$PUSH" = 1 ]]; then
    docker push -a hepstore/rivet-herwig
fi
