#! /usr/bin/env bash

RIVETBS_VERSION=4.0.2
LHAPDF_VERSION=6.5.5
HEPMC_VERSION=3.3.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

if [[ -z "$OS" || -x "$TOOLS" ]]; then
    echo "You need to set OS=ubuntu|fedora and TOOLS=gcc|llvm variables"
    exit 1
fi

PKG=hepstore/hepbase-$OS-$TOOLS
BASEARGS="--build-arg RIVETBS_VERSION=$RIVETBS_VERSION"
BASEARGS+=" --build-arg LHAPDF_VERSION=$LHAPDF_VERSION"
BASEARGS+=" --build-arg HEPMC_VERSION=$HEPMC_VERSION"
BASEARGS+=" --build-arg BUILD_TOOLS=$TOOLS"

echo "Building Docker $PKG image"
dx_build $BUILDFLAGS -f Dockerfile.$OS $BASEARGS -t $PKG:latest
test "$PUSH" = 1 && xdocker push $PKG
echo $?
