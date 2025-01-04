ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y rsync libgsl-dev libboost-all-dev \
    && apt-get install -y cython3 python3-six \
    && apt-get -y autoremove \
    && apt-get -y autoclean

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y rsync libgsl-dev time \
    && apt-get -y autoremove \
    && apt-get -y autoclean

# # TODO: fix local/local in packages, remove from here
# RUN export DEBIAN_FRONTEND=noninteractive \
#     && echo 'export PYTHONPATH="/usr/local/local/lib/python3.11/dist-packages:$PYTHONPATH"' >> /etc/profile.d/05-usrlocal.sh

ARG THEPEG_VERSION
RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && if [[ ${THEPEG_VERSION} =~ '[0-9]\.[0-9]\.[0-9]' ]]; then \
         wget https://thepeg.hepforge.org/downloads/ThePEG-${THEPEG_VERSION}.tar.bz2 -O - | tar xj && \
         cd ThePEG-${THEPEG_VERSION}; \
       else \
         apt-get install -y mercurial && \
         hg clone https://phab.hepforge.org/source/thepeghg -r ${THEPEG_VERSION} ThePEG && \
         cd ThePEG && autoreconf -if; \
       fi \
    && ./configure --help \
    && ./configure --prefix=/usr/local --with-{fastjet,rivet,hepmc,lhapdf}=/usr/local --with-hepmcversion=3 \
    && make -j$(nproc --ignore=1) \
    && make check && make install

RUN export DEBIAN_FRONTEND=noninteractive \
    && . /etc/profile \
    && for i in cteq6l1 CT10 CT14lo CT14nlo MSTW2008nlo68cl MMHT2014lo68cl NNPDF23_nlo_as_0119 NNPDF30_nlo_as_0118; do lhapdf install --verbose $i; done

ARG HERWIG_VERSION
RUN export DEBIAN_FRONTEND=noninteractive \
    && . /etc/profile \
    && cd /code \
    && if [[ ${THEPEG_VERSION} =~ '[0-9]\.[0-9]\.[0-9]' ]]; then \
         wget https://herwig.hepforge.org/downloads/Herwig-${HERWIG_VERSION}.tar.bz2 -O - | tar xj && \
         cd Herwig-${HERWIG_VERSION}; \
       else \
         apt-get install -y mercurial && \         
         hg clone https://phab.hepforge.org/source/herwighg -r ${HERWIG_VERSION} Herwig && \
         cd Herwig && autoreconf -if; \
       fi \
    && ./configure --help \
    && ./configure --prefix=/usr/local --with-{thepeg,fastjet}=/usr/local --without-evtgen \
    && make -j$(nproc --ignore=1) \
    && make check && make install

# RUN export DEBIAN_FRONTEND=noninteractive \
#     && . /etc/profile \
#     && cd /code \
#     && cd Herwig-${HERWIG_VERSION} \
#     && sed 's/^read EvtGenDecayer.in/#&/' src/defaults/decayers.in.in -i \
#     && tail src/defaults/decayers.in.in \
#     && ./config.status src/defaults/decayers.in \
#     && tail src/defaults/decayers.in \
#     && make check \
#     && make install

WORKDIR /work
