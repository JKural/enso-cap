#!/usr/bin/env bash

git submodule update --init --recursive
if ! git submodule update --init --recursive; then
  echo "Error: Failed to initialize the submodules"
  echo "Check internet connection"
  exit 1
fi

if ! cmake -B build -DCMAKE_BUILD_TYPE=Release; then
  echo "Error: CMake configuration failed."
  echo "Make sure that:"
  echo "1. All dependencies are installed"
  echo "2. Your CMake version is 3.13 or greater"
  echo "3. Your C++ compiler supports C++17"
  echo ""
  echo "In case of encountering a different error we kindly ask you to send"
  echo "the CMake output logs to jakub.kural@student.uj.edu.pl"
  exit 2
fi

if ! cmake --build build -j"$(nproc)"; then
  echo "Error: Build failed."
  echo "Make sure that:"
  echo "1. All dependencies are installed"
  echo "2. Your CMake version is 3.13 or greater"
  echo "3. Your C++ compiler supports C++17"
  echo ""
  echo "In case of encountering a different error we kindly ask you to send"
  echo "the CMake output logs to jakub.kural@student.uj.edu.pl"
  exit 3
fi
