#!/bin/bash

COREIR_VERSION=2f5035c6712481346c20c0a9956dc3f880cac1d2

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

usage () {
    cat <<EOF
Usage: $0 [<option> ...]

Sets up CoreIR.

-h, --help              display this message and exit
--python                build python bindings (default: off)
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
        --python)
            ENABLE_PYTHON=yes
            ;;
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
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$DEPS/coreir/local
    make -j$(nproc)
    make install
    if [[ "$ENABLE_PYTHON" != default ]]; then
        echo "Pip installing coreir Python bindings"
        git clone https://github.com/leonardt/pycoreir.git
        cd pycoreir
        export PATH=$PATH:$DEPS/coreir/local/bin
        which -a coreir
        if [[ "$?" != 0 ]]; then
            echo "It looks like building CoreIR failed. Unable to find coreir binary in $DEPS/coreir/local/bin."
            echo "Cannot install coreir python bindings without locating binary."
            exit 1
        fi
        pip install -e .
    fi
    cd $DIR
    cd ../
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

