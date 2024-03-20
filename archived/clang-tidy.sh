#!/bin/bash

# clang-tidy check
echo '==== clang-tidy ===='
clang-tidy-8 src/*.cpp -- -I./src -I/opt/intel/mkl/include
clang-tidy-8 examples/perovskite/*.cpp -- -I./src -I./examples/perovskite -I/opt/intel/mkl/include
clang-tidy-8 examples/diffusion_2d/*.cpp -- -I./src -I./examples/diffusion_2d -I/opt/intel/mkl/include
clang-tidy-8 examples/robertson/*.cpp -- -I./src -I./examples/robertson -I/opt/intel/mkl/include
echo

# cppcheck
echo '==== cppcheck ===='
echo 'perovskite:'
cppcheck --enable=all --std=c++11 --quiet -I./src/ -I./examples/perovskite -I/opt/intel/mkl/include/ src/*.cpp examples/perovskite/*.cpp
echo 'diffusion_2d:'
cppcheck --enable=all --std=c++11 --quiet -I./src/ -I./examples/diffusion_2d -I/opt/intel/mkl/include/ src/*.cpp examples/diffusion_2d/*.cpp
echo 'robertson:'
cppcheck --enable=all --std=c++11 --quiet -I./src/ -I./examples/robertson -I/opt/intel/mkl/include/ src/*.cpp examples/robertson/*.cpp
