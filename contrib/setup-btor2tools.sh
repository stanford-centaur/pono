#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

mkdir -p $DEPS

if [ ! -d "$DEPS/btor2tools" ]; then
    cd $DEPS
    git clone --depth 1 https://github.com/Boolector/btor2tools.git btor2tools
    cd btor2tools
    CFLAGS="" ./configure.sh -fPIC
    make -j${NPROC}
    cd $DIR
else
    echo "$DEPS/btor2tools already exists. If you want to rebuild, please remove it manually."
fi

if [ -f $DEPS/btor2tools/build/libbtor2parser.a ] ; then \
    echo "It appears btor2tools was successfully built in $DEPS/btor2tools/build."
    echo "You may now build cosa2 with: ./configure.sh && cd build && make"
else
    echo "Building btor2tools failed."
    echo "You might be missing some dependencies."
    echo "Please see their github page for installation instructions: https://github.com/Boolector/btor2tools.git"
    exit 1
fi

