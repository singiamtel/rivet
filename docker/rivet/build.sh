#! /usr/bin/env bash

set -e

YODA_BRANCH=${YODA_BRANCH:-yoda-2.0.2}
RIVET_BRANCHES=${RIVET_BRANCHES:-rivet-4.0.2}


#PLATFLAGS="--platform linux/amd64,linux/arm64"
BUILDFLAGS="$PLATFLAGS $DOCKERFLAGS"
function xdocker { echo "docker $@"; docker "$@"; }
if [[ -n "$PLATFLAGS" ]]; then
    if [[ "$PUSH" = 1 ]]; then PUSH="--push"; fi
    function dx_build { xdocker buildx build -f Dockerfile --progress=plain "$@" $PUSH .; }
else
    function dx_build { xdocker build -f Dockerfile --progress=plain "$@" .; }
fi

test "$FORCE" && BUILDFLAGS="$BUILDFLAGS --no-cache"

BUILDFLAGS="$BUILDFLAGS --build-arg YODA_BRANCH=$YODA_BRANCH"

## Last branch name -> latest
for RIVET_BRANCH in $RIVET_BRANCHES; do
    RIVET_VERSION=${RIVET_BRANCH#rivet-}
    RIVET_VERSION=${RIVET_VERSION#release-}

    MSG="Building Rivet $RIVET_VERSION image with architecture ="
    for CC in gcc; do  # clang
        ARCH=ubuntu-$CC-hepmc3-py3$ARM
        echo "@@ Building Rivet $RIVET_VERSION image with architecture = $ARCH"
        tags="hepstore/rivet:$RIVET_VERSION-$ARCH$ARM hepstore/rivet:$RIVET_VERSION$ARM"
        TAGFLAGS=""; for t in $tags; do TAGFLAGS="$TAGFLAGS -t $t"; done
        dx_build $BUILDFLAGS --build-arg RIVET_BRANCH=$RIVET_BRANCH --build-arg ARCH=$ARCH $TAGFLAGS
        if [[ "$PUSH" = 1 ]]; then
            for tag in $tags; do
                xdocker push $tag
            done
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
