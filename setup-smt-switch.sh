#!/bin/bash

git clone https://github.com/makaimann/smt-switch
cd smt-switch
./contrib/setup-btor.sh
./configure.sh --btor --prefix=local
cd build
make -j$(nproc)
make test
make install
cd ../../
