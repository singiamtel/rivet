ARG ARCH=ubuntu-gcc-hepmc3-py3
FROM hepstore/hepbase-${ARCH}-latex
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "--login", "-c"]

ARG YODA_BRANCH
ARG RIVET_BRANCH
RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && texhash \
    && git clone --depth 1 --branch $YODA_BRANCH https://gitlab.com/hepcedar/yoda.git yoda/ \
    && cd yoda/ \
    && autoconf --version \
    && autoreconf -i \
    && ./configure && make -j $(nproc --ignore=1) && make install \
    && cd .. \
    && rm -rf /code

RUN export DEBIAN_FRONTEND=noninteractive \
    && mkdir /code && cd /code \
    && git clone --depth 1 --branch $RIVET_BRANCH https://gitlab.com/hepcedar/rivet.git rivet/ \
    && cd rivet/ \
    && autoconf --version \
    && autoreconf -i \
    && if [[ -d /usr/local/include/HepMC3 ]]; then HEPMC_FLAG="--with-hepmc3=/usr/local"; fi \
    && ./configure $HEPMC_FLAG && make -j $(nproc --ignore=1) && make install \
    && texhash \
    && echo "source /usr/local/etc/bash_completion.d/rivet-completion" > /etc/profile.d/rivet-completion.sh \
    && echo "source /usr/local/etc/bash_completion.d/yoda-completion" > /etc/profile.d/yoda-completion.sh \
    && texhash \
    && rm -rf /code

# TODO: fix local/local in packages
RUN export DEBIAN_FRONTEND=noninteractive \
    && echo 'export PYTHONPATH="/usr/local/lib/python3.10/dist-packages:$PYTHONPATH"' >> /etc/profile.d/05-usrlocal.sh \
    && echo 'export PYTHONPATH="/usr/local/local/lib/python3.10/dist-packages:$PYTHONPATH"' >> /etc/profile.d/05-usrlocal.sh \
    && apt-get install -y cm-super dvipng \
    && apt-get autoremove -y \
    && apt-get autoclean -y

WORKDIR /work
