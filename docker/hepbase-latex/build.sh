#! /usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

if [[ -z "$OS" || -x "$TOOLS" ]]; then
    echo "You need to set OS=ubuntu|fedora and TOOLS=gcc|llvm variables"
    exit 1
fi

PKG=hepstore/hepbase-$OS-$TOOLS-latex
BASEARGS="--build-arg BUILD_TOOLS=$TOOLS"

echo "Building Docker $PKG image"
dx_build $BUILDFLAGS -f Dockerfile.$OS $BASEARGS -t $PKG:latest
test "$PUSH" = 1 && xdocker push $PKG && sleep $SLEEP
