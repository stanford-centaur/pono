#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

SMT_SWITCH_VERSION=499621b009ca0e86a65339aed59157bfaa874776

usage () {
    cat <<EOF
Usage: $0 [<option> ...]

Sets up the smt-switch API for interfacing with SMT solvers through a C++ API.

-h, --help              display this message and exit
--with-msat             include MathSAT which is under a custom non-BSD compliant license (default: off)
--cvc5-home             use an already downloaded version of cvc5
--python                build python bindings (default: off)
EOF
    exit 0
}

die () {
    echo "*** configure.sh: $*" 1>&2
    exit 1
}

WITH_MSAT=default
CONF_OPTS=""
WITH_PYTHON=default
cvc5_home=default

while [ $# -gt 0 ]
do
    case $1 in
        -h|--help) usage;;
        --with-msat)
            WITH_MSAT=ON
            CONF_OPTS="$CONF_OPTS --msat --msat-home=../mathsat";;
        --python)
            WITH_PYTHON=YES
            CONF_OPTS="$CONF_OPTS --python";;
        --cvc5-home) die "missing argument to $1 (see -h)" ;;
        --cvc5-home=*)
            cvc5_home=${1##*=}
            # Check if cvc5_home is an absolute path and if not, make it
            # absolute.
            case $cvc5_home in
                /*) ;;                            # absolute path
                *) cvc5_home=$(pwd)/$cvc5_home ;; # make absolute path
            esac
            CONF_OPTS="$CONF_OPTS --cvc5-home=$cvc5_home"
            ;;
        *) die "unexpected argument: $1";;
    esac
    shift
done

mkdir -p $DEPS

if [ ! -d "$DEPS/smt-switch" ]; then
    cd $DEPS
    git clone https://github.com/makaimann/smt-switch
    cd smt-switch
    git checkout -f $SMT_SWITCH_VERSION
    ./contrib/setup-btor.sh
    if [ $cvc5_home = default ]; then
        ./contrib/setup-cvc5.sh
    fi
    # pass bison/flex directories from smt-switch perspective
    ./configure.sh --btor --cvc5 $CONF_OPTS --prefix=local --static --smtlib-reader --bison-dir=../bison/bison-install --flex-dir=../flex/flex-install
    cd build
    make -j$(nproc)
    # TODO put this back
    # temporarily disable due to test-disjointset issue
    # make test
    make install
    cd $DIR
else
    echo "$DEPS/smt-switch already exists. If you want to rebuild, please remove it manually."
fi

if [ 0 -lt $(ls $DEPS/smt-switch/local/lib/libsmt-switch* 2>/dev/null | wc -w) ]; then
    echo "It appears smt-switch with boolector and cvc5 was successfully installed to $DEPS/smt-switch/local."
    echo "You may now build pono with: ./configure.sh && cd build && make"
else
    echo "Building smt-switch failed."
    echo "You might be missing some dependencies."
    echo "Please see the github page for installation instructions: https://github.com/makaimann/smt-switch"
    exit 1
fi
