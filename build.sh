#!/usr/bin/env bash

git submodule update --init --recursive

cmake --preset release
cmake --build build -j"$(nproc)"
