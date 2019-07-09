#!/bin/bash
clang-tidy-6.0 src/*.cpp -- -I/opt/intel/mkl/include
clang-tidy-6.0 examples/perovskite/*.cpp -- -I/opt/intel/mkl/include
clang-tidy-6.0 examples/diffusion_2d/*.cpp -- -I/opt/intel/mkl/include
clang-tidy-6.0 examples/robertson/*.cpp -- -I/opt/intel/mkl/include
