#!/bin/sh
clear
rm -rf build
mkdir build
cd build
cmake ..
make
./test_ML
cd ..
