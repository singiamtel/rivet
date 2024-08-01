#! /usr/bin/env bash

set -e

RIVET_VERSION=${RIVET_VERSION:-4.0.1}
PYTHIA_VERSION=${PYTHIA_VERSION:-8311}  # TODO: 8312 changes the example names
MG5_URL=https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.5.tar.gz

#PLATFLAGS="--platform linux/amd64,linux/arm64"
BUILDFLAGS="$PLATFLAGS $DOCKERFLAGS"
function xdocker { echo "docker $@"; docker "$@"; }
if [[ -n "$PLATFLAGS" ]]; then
    if [[ "$PUSH" = 1 ]]; then PUSH="--push"; fi
    function dx_build { xdocker buildx build -f Dockerfile --progress=plain "$@" $PUSH .; }
else
    function dx_build { xdocker build -f Dockerfile --progress=plain "$@" .; }
fi
#test "$FORCE" && BUILD="$BUILD --no-cache"

BUILDFLAGS="$BUILDFLAGS --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILDFLAGS="$BUILDFLAGS --build-arg PYTHIA_VERSION=${PYTHIA_VERSION}"
BUILDFLAGS="$BUILDFLAGS --build-arg MG5_URL=${MG5_URL}"
test "$TEST" && BUILDFLAGS="echo $BUILDFLAGS"

tag="hepstore/rivet-tutorial:${RIVET_VERSION}"
echo "Building $tag"
dx_build $BUILDFLAGS -t $tag

if [[ "$LATEST" = 1 ]]; then
    docker tag $tag hepstore/rivet-tutorial:latest
fi

if [[ "$PUSH" = 1 ]]; then
    docker push $tag
    if [[ "$LATEST" = 1 ]]; then
        sleep ${SLEEP:-1}m
        docker push hepstore/rivet-tutorial:latest
    fi
fi
