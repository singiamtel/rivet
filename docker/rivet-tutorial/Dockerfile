ARG RIVET_VERSION
ARG PYTHIA_VERSION
FROM hepstore/rivet-pythia:${RIVET_VERSION}-${PYTHIA_VERSION}
LABEL maintainer="rivet-developers@cern.ch"

CMD /bin/bash

##########

ARG MG5_URL
RUN wget --no-verbose $MG5_URL -O- | tar xz

WORKDIR /work

RUN true \
  && apt-get update -y && apt-get install -y gnuplot-nox && apt-get autoclean -y \
  && ln -s MG5_aMC_v* MG5_aMC \
  && cd MG5_aMC \
  && echo -e "set auto_update 0\ninstall pythia8\n" > mg5setup.txt \
  && bin/mg5_aMC mg5setup.txt \
  && rm mg5setup.txt

RUN lhapdf install NNPDF23_lo_as_0130_qed NNPDF40_nlo_as_01180

RUN export DEBIAN_FRONTEND=noninteractive \
    && echo 'export PYTHONPATH="/usr/local/local/lib/python3.11/dist-packages:$PYTHONPATH"' >> /etc/profile.d/05-usrlocal.sh

##########

WORKDIR /work
ADD . /work
RUN true \
  && cd MG5_aMC/ \
  && patch -p1 < /work/mgpatch.diff \
  && cd /work && rm Dockerfile build.sh* *.diff
