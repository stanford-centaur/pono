#!/bin/bash
set -e

SLANG_VERSION=v7.0

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

mkdir -p $DEPS

if [ ! -d "$DEPS/slang" ]; then
    cd $DEPS
    git clone https://github.com/MikePopoloski/slang.git slang
    cd slang
    git checkout -f $SLANG_VERSION
    cmake -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$DEPS/slang/local" \
        -DSLANG_INCLUDE_TESTS=OFF \
        -DSLANG_INCLUDE_TOOLS=OFF \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    cmake --build build -j$(nproc)
    cmake --install build
    cd $DIR
else
    echo "$DEPS/slang already exists. If you want to rebuild, please remove it manually."
fi

if [ -d "$DEPS/slang/local/include/slang" ] && \
   [ -f "$DEPS/slang/local/lib/libsvlang.a" -o \
     -f "$DEPS/slang/local/lib/libslang.a" -o \
     -f "$DEPS/slang/local/lib64/libsvlang.a" -o \
     -f "$DEPS/slang/local/lib64/libslang.a" ]; then
    echo "It appears slang was successfully built in $DEPS/slang."
    echo "You may now build pono with: ./configure.sh --with-slang && cd build && make"
else
    echo "Building slang failed."
    echo "You might be missing some dependencies."
    echo "Please see: https://github.com/MikePopoloski/slang"
    exit 1
fi
