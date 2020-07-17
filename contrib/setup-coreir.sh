#!/bin/bash

COREIR_VERSION=f49da00d96e41a05ec32befc10e3115de4d71d4d

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

usage () {
    cat <<EOF
Usage: $0 [<option> ...]

Sets up the smt-switch API for interfacing with SMT solvers through a C++ API.

-h, --help              display this message and exit
EOF
    exit 0
}

die () {
    echo "*** configure.sh: $*" 1>&2
    exit 1
}

ENABLE_PYTHON=default

while [ $# -gt 0 ]
do
    case $1 in
        -h|--help) usage;;
        *) die "unexpected argument: $1";;
    esac
    shift
done


mkdir -p $DEPS

if [ ! -d "$DEPS/coreir" ]; then
    cd $DEPS
    git clone https://github.com/rdaly525/coreir.git
    cd coreir
    git checkout -f $COREIR_VERSION
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$DEPS/coreir/local
    make -j2
    make install
    cd $DIR
    cd ../
else
    echo "$DEPS/coreir already exists. If you want to rebuild, please remove it manually."
fi

if [ -f $DEPS/coreir/local/bin/coreir ] ; then \
    echo "It appears coreir was successfully built in $DEPS/coreir/local."
    echo "You may now build pono with: ./configure.sh && cd build && make"
else
    echo "Building coreir failed."
    echo "You might be missing some dependencies."
    echo "Please see their github page for installation instructions:https://github.com/rdaly525/coreir"
    exit 1
fi

