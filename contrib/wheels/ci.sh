#!/usr/bin/env bash
set -e

cd cosa2/
mkdir -p build
pip install Cython==0.29
python contrib/wheels/build_wheel.py bdist_wheel
auditwheel show dist/*
auditwheel repair dist/*
