#! /usr/bin/env bash

set -e

RIVET_VERSION=${RIVET_VERSION:-4.0.1}
MG5_VERSION=3.5.5
MG5_URL=https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.5.tar.gz

#PLATFLAGS="--platform linux/amd64,linux/arm64"
#BUILD="docker buildx build -f Dockerfile $PLATFLAGS $DOCKERFLAGS ."
#if [[ -n "$PLATFLAGS" && "$PUSH" = 1 ]]; then BUILD="$BUILD --push"; fi

BUILD="docker build . -f Dockerfile --progress=plain $DOCKERFLAGS"

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg MG5_URL=${MG5_URL}"

BUILD="$BUILD -t hepstore/rivet-mg5amcnlo:${RIVET_VERSION}-${MG5_VERSION}"
BUILD="$BUILD -t hepstore/rivet-mg5amcnlo:${RIVET_VERSION}"
if [[ "$LATEST" = 1 ]]; then BUILD="$BUILD -t hepstore/rivet-mg5amcnlo:latest"; fi

echo "Building: $BUILD"
$BUILD

if [[ "$PUSH" = 1 ]]; then
    docker push -a hepstore/rivet-mg5amcnlo
fi
