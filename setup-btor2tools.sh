#!/bin/bash

git clone --depth 1 https://github.com/Boolector/btor2tools.git btor2tools
cd btor2tools

./configure.sh -fPIC
make -j${NPROC}
cd ../
