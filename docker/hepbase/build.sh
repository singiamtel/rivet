#! /usr/bin/env bash

RIVETBS_VERSION=4.0.2
LHAPDF_VERSION=6.5.4
HEPMC_VERSION=3.3.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

BASEARGS="--build-arg RIVETBS_VERSION=$RIVETBS_VERSION"
BASEARGS+=" --build-arg LHAPDF_VERSION=$LHAPDF_VERSION"
BASEARGS+=" --build-arg HEPMC_VERSION=3.3.0"

if [[ -z "$OS" || -x "$TOOLS" ]]; then
    echo "You need to set OS=ubuntu|fedora and TOOLS=gcc|llvm variables"
    exit 1
fi

PKG=hepstore/hepbase-$OS-$TOOLS

echo "Building $PKG image"
dx_build $BUILDFLAGS -f Dockerfile.$OS --build-arg LATEX=0 --build-arg BUILD_TOOLS=$TOOLS -t $PKG:latest
test "$PUSH" = 1 && xdocker push $PKG && sleep $SLEEP

PKG+="-latex"
echo "Building $PKG image"
dx_build $BUILDFLAGS -f Dockerfile.$OS --build-arg LATEX=1 --build-arg BUILD_TOOLS=$TOOLS -t $PKG:latest
test "$PUSH" = 1 && xdocker push $PKG && sleep $SLEEP
