#! /usr/bin/env bash

RIVET_VERSION=${RIVET_VERSION:-4.0.2}
MG5_VERSION=3.5.7
MG5_URL=https://launchpad.net/mg5amcnlo/3.0/3.6.x/+download/MG5_aMC_v3.5.7.tar.gz

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BUILDFLAGS+=" --build-arg RIVET_VERSION=${RIVET_VERSION}"
BUILDFLAGS+=" --build-arg MG5_URL=${MG5_URL}"

PKG="hepstore/rivet-mg5amcnlo"
TAGS="${RIVET_VERSION}-${MG5_VERSION} ${RIVET_VERSION}"

source "$SCRIPT_DIR/../docker-build-and-push.sh"

# test "$LATEST" = 1 && TAGS+=" latest"
# for tag in $TAGS; do TAGFLAGS+=" -t $PKG:$tag"; done

# dx_build $BUILDFLAGS $TAGFLAGS

# if [[ "$PUSH" = 1 ]]; then
#     for tag in $TAGS; do xdocker push $PKG:$tag; done
# fi
