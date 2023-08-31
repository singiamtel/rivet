#! /usr/bin/env bash

set -e

#PLATFLAGS="--platform linux/amd64,linux/arm64"
BUILDFLAGS="$PLATFLAGS $DOCKERFLAGS"
function xdocker { echo "docker $@"; docker "$@"; }
if [[ -n "$PLATFLAGS" ]]; then
    if [[ "$PUSH" = 1 ]]; then PUSH="--push"; fi
    function dx_build { xdocker buildx build -f Dockerfile "$@" $PUSH .; }
else
    function dx_build { xdocker build -f Dockerfile "$@" .; }
fi

test "$FORCE" && BUILDFLAGS="$BUILDFLAGS --no-cache"

YODA_BRANCH=yoda-1.9.8
BUILDFLAGS="$BUILDFLAGS --build-arg YODA_BRANCH=$YODA_BRANCH"

## Last branch name -> latest
for RIVET_BRANCH in rivet-3.1.8; do
    RIVET_VERSION=${RIVET_BRANCH#rivet-}

    MSG="Building Rivet $RIVET_VERSION image with architecture ="
    for CC in gcc; do  # clang
        ARCH=ubuntu-$CC-hepmc3-py3$ARM
        echo "@@ Building Rivet $RIVET_VERSION image with architecture = $ARCH"
        tags="hepstore/rivet:$RIVET_VERSION-$arch$ARM hepstore/rivet:$RIVET_VERSION$ARM"
        TAGFLAGS=""; for t in $tags; do TAGFLAGS="$TAGFLAGS -t $t"; done
        dx_build $BUILDFLAGS --build-arg RIVET_BRANCH=$RIVET_BRANCH --build-arg ARCH=$ARCH $TAGFLAGS
        if [[ "$PUSH" = 1 ]]; then
            xdocker push $tag
            sleep ${SLEEP:-1}m
        fi
        echo -e "\n\n\n"
    done
done

## Convenience tags
if [[ "$LATEST" = 1 ]]; then
    dx_build $BUILDFLAGS --build-arg RIVET_BRANCH=$RIVET_BRANCH --build-arg ARCH=$ARCH $TAGFLAGS -t hepstore/rivet:latest
    if [[ "$PUSH" = 1 ]]; then
        xdocker push hepstore/rivet:latest
    fi
fi
