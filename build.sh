#!/bin/bash

mkdir -p build
rm -rf build/*
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="~/vcpkg/scripts/buildsystems/vcpkg.cmake" -GNinja ..

cd ..

cmake --build build