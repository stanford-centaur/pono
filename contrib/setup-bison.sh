#!/bin/bash
set -e

BISON_VERSION=3.7.4

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

mkdir -p $DEPS

if [ -d "$DEPS/bison" ]; then
    echo "It appears bison has already been downloaded to $DEPS/bison"
    echo "If you'd like to rebuild, please delete it and run this script again"
    exit 1
fi

curl http://ftp.gnu.org/gnu/bison/bison-$BISON_VERSION.tar.gz --output $DEPS/bison-$BISON_VERSION.tar.gz

if [ ! -f "$DEPS/bison-$BISON_VERSION.tar.gz" ]; then
    echo "It seems like downloading bison to $DEPS/bison-$BISON_VERSION.tar.gz failed"
    exit 1
fi

cd $DEPS
tar -xzf bison-$BISON_VERSION.tar.gz
rm bison-$BISON_VERSION.tar.gz
mv ./bison-$BISON_VERSION ./bison
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
