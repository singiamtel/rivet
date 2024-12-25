#! /usr/bin/env bash

YODA_BRANCH=${YODA_BRANCH:-yoda-2.0.2}
RIVET_BRANCHES=${RIVET_BRANCHES:-rivet-4.0.2}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BUILDFLAGS==" --build-arg YODA_BRANCH=$YODA_BRANCH"
PKG="hepstore/rivet"

for RIVET_BRANCH in $RIVET_BRANCHES; do
    RIVET_VERSION=${RIVET_BRANCH#rivet-}

    MSG="Building Rivet $RIVET_VERSION image with architecture ="
    for CC in gcc; do  # clang  #< last architecture in the list will also get the no-arch name
        ARCH=ubuntu-$CC-hepmc3-py3
        echo "@@ Building Rivet $RIVET_VERSION image with architecture = $ARCH"
        TAGS="$PKG:$RIVET_VERSION-$ARCH $PKG:$RIVET_VERSION"
        TAGFLAGS=""; for t in $TAGS; do TAGFLAGS="$TAGFLAGS -t $t"; done
        dx_build $BUILDFLAGS --build-arg RIVET_BRANCH=$RIVET_BRANCH --build-arg ARCH=$ARCH $TAGFLAGS
        if [[ "$PUSH" = 1 ]]; then
            for tag in $TAGS; do
                xdocker push $tag
            done
            sleep ${SLEEP:-1}
        fi
        echo -e "\n\n\n"
    done
done

if [[ "$LATEST" = 1 ]]; then
    dx_build $BUILDFLAGS --build-arg RIVET_BRANCH=$RIVET_BRANCH --build-arg ARCH=$ARCH $TAGFLAGS -t $PKG:latest
    if [[ "$PUSH" = 1 ]]; then
        xdocker push hepstore/rivet:latest
    fi
fi
