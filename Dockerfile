FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

COPY . /app/pono/
WORKDIR /app/pono/

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y \
        build-essential \
        cmake \
        curl \
        git \
        libgmp-dev \
        m4 \
        meson \
        ninja-build \
        pkg-config \
        python3 \
        python3-pyparsing \
        wget
RUN apt-get autoremove --purge
RUN apt-get clean -y
RUN ./contrib/setup-bison.sh
RUN ./contrib/setup-flex.sh
RUN ./contrib/setup-btor2tools.sh

## To setup MathSAT dependencies, uncomment the following line:
# RUN ./ci-scripts/setup-msat.sh -y
## Then, add --with-msat to the command below
RUN ./contrib/setup-smt-switch.sh

## To setup IC3ia dependencies, uncomment the following line:
# RUN ./contrib/setup-ic3ia.sh
## To build Pono with MathSAT and IC3ia,
## add --with-msat and --with-msat-ic3ia to the command below, respectively
RUN ./configure.sh --static

RUN cd build/ && make -j$(nproc)

RUN ln -s /app/pono/build/pono /usr/local/bin/pono
ENTRYPOINT ["/usr/local/bin/pono"]
