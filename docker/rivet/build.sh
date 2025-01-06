#! /usr/bin/env bash

YODA_BRANCH=${YODA_BRANCH:-yoda-2.0.2}
RIVET_BRANCHES=${RIVET_BRANCHES:-rivet-4.0.2}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BUILDFLAGS+=" --build-arg YODA_BRANCH=$YODA_BRANCH"
BUILDFLAGS+=" --build-arg ARCH=ubuntu-gcc-hepmc3-py3"

PKG="hepstore/rivet"

if [[ "$MAIN" = 1 ]]; then RIVET_BRANCHES="main $RIVET_BRANCHES"; fi

for RIVET_BRANCH in $RIVET_BRANCHES; do
    RIVET_VERSION=${RIVET_BRANCH#rivet-}
    MSG="Building $PKG with Rivet=$RIVET_VERSION"

    TAGS="${RIVET_VERSION}"
    test "$LATEST" = 1 && TAGS+=" latest"
    TAGFLAGS=""; for tag in $TAGS; do TAGFLAGS+=" -t $PKG:$tag"; done

    echo "$MSG -> $PKG:${TAGS// /,}"
    dx_build $BUILDFLAGS --build-arg RIVET_BRANCH=$RIVET_BRANCH $TAGFLAGS

    if [[ "$PUSH" = 1 ]]; then
        for tag in $TAGS; do xdocker push $PKG:$tag; done
    fi
    echo -e "\n\n\n"
done

# if [[ "$LATEST" = 1 ]]; then
#     dx_build $BUILDFLAGS --build-arg RIVET_BRANCH=$RIVET_BRANCH --build-arg ARCH=$ARCH $TAGFLAGS -t $PKG:latest
#     if [[ "$PUSH" = 1 ]]; then
#         xdocker push hepstore/rivet:latest
#     fi
# fi
