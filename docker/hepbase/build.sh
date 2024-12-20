#! /usr/bin/env bash

RIVETBS_VERSION=4.0.2
LHAPDF_VERSION=6.5.4

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../docker-common.sh"

#PLATFLAGS="--platform linux/amd64,linux/arm64"
#BUILD="docker buildx build $PLATFLAGS $DOCKERFLAGS" # --squash"
BUILD="docker build --progress=plain $PLATFLAGS $DOCKERFLAGS" # --squash"

test "$FORCE" && BUILD="$BUILD --no-cache"

test "$TEST" && BUILD="echo $BUILD"


for vhepmc in 3.2.7; do # 3.3.0
    for tex in 0 1; do

        MSG="Building hepbase image with HepMC=$vhepmc and TeX=$tex"

        BASEARGS="--build-arg RIVETBS_VERSION=$RIVETBS_VERSION"
        BASEARGS="$BASEARGS --build-arg HEPMC_VERSION=$vhepmc"
        BASEARGS="$BASEARGS --build-arg LHAPDF_VERSION=$LHAPDF_VERSION"
        BASEARGS="$BASEARGS --build-arg LATEX=$tex"

        GCCARGS="$BASEARGS --build-arg BUILD_TOOLS=GCC"
        CLANGARGS="$BASEARGS --build-arg BUILD_TOOLS=LLVM"
        INTELARGS="$BASEARGS --build-arg BUILD_TOOLS=Intel"
        if [[ "$tex" = 1 ]]; then TEXSUFFIX="-latex"; else TEXSUFFIX=""; fi

        echo "@@ $MSG on Ubuntu with GCC compilers"
        tag=hepstore/hepbase-ubuntu-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX$ARM
        $BUILD -f Dockerfile.ubuntu $GCCARGS -t $tag .
        test "$PUSH" = 1 && docker push $tag && sleep ${SLEEP:-1}m
        echo -e "\n\n\n"

        echo "@@ $MSG on Ubuntu with clang compilers"
        tag=hepstore/hepbase-ubuntu-clang-hepmc${vhepmc:0:1}-py3$TEXSUFFIX$ARM
        $BUILD -f Dockerfile.ubuntu ${CLANGARGS/gfortran/flang} -t $tag .
        test "$PUSH" = 1 && docker push $tag && sleep ${SLEEP:-1}m
        echo -e "\n\n\n"

        echo "@@ $MSG on Fedora with GCC compilers"
        tag=hepstore/hepbase-fedora-gcc-hepmc${vhepmc:0:1}-py3$TEXSUFFIX$ARM
        $BUILD -f Dockerfile.fedora $GCCARGS -t $tag .
        test "$PUSH" = 1 && docker push $tag && sleep ${SLEEP:-1}m
        echo -e "\n\n\n"

    done
done
