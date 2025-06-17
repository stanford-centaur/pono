#!/bin/bash
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

SMT_SWITCH_VERSION=80ff7514812a8011f417f82220f8ef1d3a62e248

usage () {
    cat <<EOF
Usage: $0 [<option> ...]

Sets up the smt-switch API for interfacing with SMT solvers through a C++ API.

-h, --help              display this message and exit
--with-btor             include Boolector (default: off)
--with-msat             include MathSAT which is under a custom non-BSD compliant license (default: off)
--with-yices2           include Yices2 which is under a custom non-BSD compliant license (default: off)
--cvc5-home             use an already downloaded version of cvc5
--python                build python bindings (default: off)
EOF
    exit 0
}

die () {
    echo "*** configure.sh: $*" 1>&2
    exit 1
}

WITH_BOOLECOR=default
WITH_MSAT=default
WITH_YICES2=default
CONF_OPTS=""
WITH_PYTHON=default
cvc5_home=default
INCL_SOLVER_STR="bitwuzla and cvc5"

while [ $# -gt 0 ]
do
    case $1 in
        -h|--help) usage;;
        --with-msat)
            WITH_MSAT=ON
            CONF_OPTS="$CONF_OPTS --msat --msat-home=../mathsat"
            INCL_SOLVER_STR="mathsat, $INCL_SOLVER_STR";;
        --with-btor)
            WITH_BOOLECTOR=ON
            CONF_OPTS="$CONF_OPTS --btor"
            INCL_SOLVER_STR="boolector, $INCL_SOLVER_STR";;
        --with-yices2)
            WITH_YICES2=ON
            CONF_OPTS="$CONF_OPTS --yices2"
            INCL_SOLVER_STR="yices2, $INCL_SOLVER_STR";;
        --python)
            WITH_PYTHON=YES
            CONF_OPTS="$CONF_OPTS --python";;
        --cvc5-home) die "missing argument to $1 (see -h)" ;;
        --cvc5-home=*)
            cvc5_home=${1##*=}
            # Check if cvc5_home is an absolute path and
            # if not, make it absolute.
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
    git clone https://github.com/stanford-centaur/smt-switch
    cd smt-switch
    git checkout -f $SMT_SWITCH_VERSION
    ./contrib/setup-bitwuzla.sh
    if [ $cvc5_home = default ]; then
        ./contrib/setup-cvc5.sh
    fi
    if [ $WITH_BOOLECTOR = ON ]; then
        ./contrib/setup-btor.sh
    fi
    # pass bison/flex directories from smt-switch perspective
    ./configure.sh --bitwuzla --cvc5 $CONF_OPTS --prefix=local --static --smtlib-reader --bison-dir=../bison/bison-install --flex-dir=../flex/flex-install
    cd build
    make -j$(nproc)
    make test
    make install
    cd $DIR
else
    echo "$DEPS/smt-switch already exists. If you want to rebuild, please remove it manually."
fi

if [ 0 -lt $(ls $DEPS/smt-switch/local/lib/libsmt-switch* 2>/dev/null | wc -w) ]; then
    echo "It appears smt-switch with $INCL_SOLVER_STR was successfully installed to $DEPS/smt-switch/local."
    echo "You may now build pono with: ./configure.sh && cd build && make"
else
    echo "Building smt-switch failed."
    echo "You might be missing some dependencies."
    echo "Please see the github page for installation instructions: https://github.com/stanford-centaur/smt-switch"
    exit 1
fi
