#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

mkdir -p $DEPS

if [ ! -d "$DEPS/smt-switch" ]; then
    cd $DEPS
    git clone -b btor-smtcomp19 https://github.com/makaimann/smt-switch
    cd smt-switch
    ./contrib/setup-btor.sh
    ./contrib/setup-msat.sh
    ./configure.sh --btor --msat --prefix=local
    cd build
    make -j$(nproc)
    make test
    make install
    cd $DIR
else
    echo "$DEPS/smt-switch already exists. If you want to rebuild, please remove it manually."
fi

if [ -f $DEPS/smt-switch/local/lib/libsmt-switch.so ] && [ -f $DEPS/smt-switch/local/lib/libsmt-switch-btor.so ] ; then \
    echo "It appears smt-switch with boolector was successfully installed to $DEPS/smt-switch/local."
    echo "You may now build cosa2 with: ./configure.sh && cd build && make"
else
    echo "Building smt-switch failed."
    echo "You might be missing some dependencies."
    echo "Please see the github page for installation instructions: https://github.com/makaimann/smt-switch"
    exit 1
fi
