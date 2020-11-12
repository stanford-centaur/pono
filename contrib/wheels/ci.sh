#!/usr/bin/env bash
set -e

cd pono/
mkdir -p build
yum install -y gmp-static libedit-devel
pip install toml setuptools pexpect Cython==0.29
python contrib/wheels/build_wheel.py bdist_wheel
auditwheel show dist/*
auditwheel repair dist/*
