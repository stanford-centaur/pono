#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

usage () {
    cat <<EOF
Usage: $0 [<option> ...]

Sets up the smt-switch API for interfacing with SMT solvers through a C++ API.

-h, --help              display this message and exit
--with-cvc4             include CVC4 (default: off)
--with-msat             include MathSAT which is under a custom non-BSD compliant license (default: off)
-y, --auto-yes          automatically agree to conditions (default: off)
EOF
    exit 0
}

die () {
    echo "*** configure.sh: $*" 1>&2
    exit 1
}

WITH_MSAT=default
WITH_CVC4=default

while [ $# -gt 0 ]
do
    case $1 in
        -h|--help) usage;;
        --with-msat)
            WITH_MSAT=ON
            CONF_OPTS="$CONF_OPTS --msat";;
        --with-cvc4)
            WITH_CVC4=ON
            CONF_OPTS="$CONF_OPTS --cvc4";;
        -y|--auto-yes) MSAT_OPTS=--auto-yes;;
        *) die "unexpected argument: $1";;
    esac
    shift
done

mkdir -p $DEPS

if [ ! -d "$DEPS/smt-switch" ]; then
    cd $DEPS
    git clone https://github.com/makaimann/smt-switch
    cd smt-switch
    ./contrib/setup-btor.sh

    if [[ "$WITH_MSAT" != default ]]; then
        ./contrib/setup-msat.sh $MSAT_OPTS
    fi

    if [[ "$WITH_CVC4" != default ]]; then
        ./contrib/setup-cvc4.sh
    fi

    ./configure.sh --btor $CONF_OPTS --prefix=local --static
    cd build
    make -j$(nproc)
    make test
    make install
    cd $DIR
else
    echo "$DEPS/smt-switch already exists. If you want to rebuild, please remove it manually."
fi

if [ 0 -lt $(ls $DEPS/smt-switch/local/lib/libsmt-switch* 2>/dev/null | wc -w) ]; then
    echo "It appears smt-switch with boolector was successfully installed to $DEPS/smt-switch/local."
    echo "You may now build cosa2 with: ./configure.sh && cd build && make"
else
    echo "Building smt-switch failed."
    echo "You might be missing some dependencies."
    echo "Please see the github page for installation instructions: https://github.com/makaimann/smt-switch"
    exit 1
fi
