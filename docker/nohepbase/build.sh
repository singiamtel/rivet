#! /usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

if [[ -z "$OS" || -x "$TOOLS" ]]; then
    echo "You need to set OS=ubuntu|fedora and TOOLS=gcc|llvm variables"
    exit 1
fi

PKG=hepstore/nohepbase-$OS-$TOOLS

echo "Building Docker $PKG image"
dx_build $BUILDFLAGS -f Dockerfile.$OS $BASEARGS --build-arg LATEX=0 --build-arg BUILD_TOOLS=$TOOLS -t $PKG:latest
test "$PUSH" = 1 && xdocker push $PKG && sleep $SLEEP

PKG+="-latex"
echo "Building Docker $PKG image"
dx_build $BUILDFLAGS -f Dockerfile.$OS $BASEARGS --build-arg LATEX=1 --build-arg BUILD_TOOLS=$TOOLS -t $PKG:latest
test "$PUSH" = 1 && xdocker push $PKG && sleep $SLEEP
