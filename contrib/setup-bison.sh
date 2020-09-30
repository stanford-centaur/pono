#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

mkdir -p $DEPS

if [ -d "$DEPS/bison" ]; then
    echo "It appears bison has already been downloaded to $DEPS/bison"
    echo "If you'd like to rebuild, please delete it and run this script again"
    exit 1
fi

curl http://ftp.gnu.org/gnu/bison/bison-3.5.tar.gz --output $DEPS/bison-3.5.tar.gz

if [ ! -f "$DEPS/bison-3.5.tar.gz" ]; then
    echo "It seems like downloading bison to $DEPS/bison-3.5.tar.gz failed"
    exit 1
fi

cd $DEPS
tar -xzf bison-3.5.tar.gz
rm bison-3.5.tar.gz
mv ./bison-3.5 ./bison
cd bison
mkdir bison-install
./configure --prefix $DEPS/bison/bison-install --exec-prefix $DEPS/bison/bison-install
make -j$(nproc)
make install
cd $DIR

if [ ! -f "$DEPS/bison/bison-install/bin/bison" ]; then
    echo "It seems like installing bison to $DEPS/bison/bison-install failed"
    exit 1
fi
