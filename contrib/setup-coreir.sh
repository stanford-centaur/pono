#!/bin/bash

COREIR_VERSION=86925c395e45917fc80a9d843ebc5cd2f74dee8c

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

mkdir -p $DEPS

if [ ! -d "$DEPS/coreir" ]; then
    cd $DEPS
    git clone https://github.com/rdaly525/coreir.git
    cd coreir
    git checkout -f $COREIR_VERSION
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$DEPS/coreir/local
    make -j${NPROC}
    make install
    cd $DIR
else
    echo "$DEPS/coreir already exists. If you want to rebuild, please remove it manually."
fi

if [ -f $DEPS/coreir/local/lib/libcoreir.so ] ; then \
    echo "It appears coreir was successfully built in $DEPS/coreir/local."
    echo "You may now build pono with: ./configure.sh && cd build && make"
else
    echo "Building coreir failed."
    echo "You might be missing some dependencies."
    echo "Please see their github page for installation instructions:https://github.com/rdaly525/coreir"
    exit 1
fi

