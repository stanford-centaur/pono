#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# First argument is the install prefix
prefix=${1:-install}

# Download MathSAT
./ci-scripts/setup-msat.sh -y

# Install dependencies
./contrib/setup-smt-switch.sh --with-msat
./contrib/setup-btor2tools.sh

# Configure
./configure.sh --prefix=$prefix --with-msat --static

# Build
make -C build -j install
