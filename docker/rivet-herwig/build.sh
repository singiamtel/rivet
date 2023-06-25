#! /usr/bin/env bash

set -e

RIVET_VERSION=3.1.8
HERWIG_VERSION=7.2.3
THEPEG_VERSION=2.2.3

PLATFLAGS="--platform linux/amd64,linux/arm64"
BUILD="docker buildx build -f Dockerfile $PLATFLAGS $DOCKERFLAGS ."

test "$FORCE" && BUILD="$BUILD --no-cache"

BUILD="$BUILD --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILD="$BUILD --build-arg THEPEG_VERSION=${THEPEG_VERSION}"
BUILD="$BUILD --build-arg HERWIG_VERSION=${HERWIG_VERSION}"
test "$TEST" && BUILD="echo $BUILD"

tag="hepstore/rivet-herwig:${RIVET_VERSION}-${HERWIG_VERSION}"
echo "Building $tag"
$BUILD -f Dockerfile -t $tag

docker tag $tag hepstore/rivet-herwig:$RIVET_VERSION
if [[ "$LATEST" = 1 ]]; then
    docker tag $tag hepstore/rivet-herwig:latest
fi

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    sleep ${SLEEP:-1}m
    docker push hepstore/rivet-herwig:$RIVET_VERSION
    if [[ "$LATEST" = 1 ]]; then
        sleep ${SLEEP:-1}m
        docker push hepstore/rivet-herwig:latest
    fi
fi
