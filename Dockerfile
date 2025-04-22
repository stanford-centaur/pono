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
        python3-pyparsing
RUN apt-get autoremove --purge
RUN apt-get clean -y
RUN ./contrib/setup-bison.sh
RUN ./contrib/setup-flex.sh
RUN ./contrib/setup-smt-switch.sh
RUN ./contrib/setup-btor2tools.sh
RUN ./configure.sh --static
RUN cd build/ && \
    make -j$(nproc)

RUN ln -s /app/pono/build/pono /usr/local/bin/pono
ENTRYPOINT ["/usr/local/bin/pono"]
