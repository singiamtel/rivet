ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}
LABEL maintainer="rivet-developers@cern.ch"

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y rsync pypy3 \
    && cd /usr/local \
    && git clone https://gitlab.com/openloops/OpenLoops.git \
    && cd OpenLoops && ./scons \
    && ./openloops libinstall ppln ppln_ew pptt \
    && apt-get -y autoremove \
    && apt-get -y autoclean

# TODO: update Sherpa to CMake build with v3

ARG SHERPA_VERSION
ENV PATH="${PATH}:/usr/lib64/openmpi/bin"
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y libtool texinfo texlive-metapost feynmf libsqlite3-dev swig openmpi-bin libopenmpi-dev \
    && mkdir /code && cd /code \
    && git clone -b v${SHERPA_VERSION} https://gitlab.com/sherpa-team/sherpa.git \
    && cd sherpa && autoreconf -fi \
    && ./configure --help \
    && ./configure --prefix=/usr/local \
      --enable-lhapdf=/usr/local --enable-fastjet=/usr/local \
      --enable-hepmc3=/usr/local --enable-openloops=/usr/local/OpenLoops \
      --enable-rivet=/usr/local --enable-analysis \
      --disable-pyext \
      --enable-gzip \
      --enable-mpi CC=mpicc CXX=mpic++ \
    && make -j $(nproc --ignore=1) CXXFLAGS="-O2" AM_CXXFLAGS="-O2" \
    && make install \
    && mkdir /usr/local/Sherpa \
    && mv Examples /usr/local/Sherpa \
    && cd / && rm -r /code \
    && apt-get -y autoremove \
    && apt-get -y autoclean

RUN export DEBIAN_FRONTEND=noninteractive \
    && echo 'export PYTHONPATH="/usr/local/local/lib/python3.11/dist-packages:$PYTHONPATH"' >> /etc/profile.d/05-usrlocal.sh

WORKDIR /work
