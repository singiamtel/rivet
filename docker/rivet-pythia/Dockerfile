ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

ARG PYTHIA_VERSION

RUN export DEBIAN_FRONTEND=noninteractive \
    && export TERM=xterm \
    && mkdir /code && cd /code \
    && wget --no-verbose https://pythia.org/download/pythia83/pythia${PYTHIA_VERSION}.tgz -O- | tar xz \
    && cd pythia*/ && ./configure --{prefix,with-{hepmc3,lhapdf6,rivet}}=/usr/local --with-gzip \
    && make -j $(nproc --ignore=1) && make install \
    && cd examples \
    && for i in 144 164; do make main$i && cp main$i /usr/local/bin/pythia8-main$i; done \
    && rm -r /code

WORKDIR /work
