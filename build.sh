#!/usr/bin/env bash

git submodule update --init

cmake --preset release
cmake --build build -j"$(nproc)"
