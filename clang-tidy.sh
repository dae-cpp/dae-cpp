#!/bin/bash
clang-tidy-6.0 src/*.cpp -- -I./src -I/opt/intel/mkl/include
clang-tidy-6.0 examples/perovskite/*.cpp -- -I./src -I./examples/perovskite -I/opt/intel/mkl/include
clang-tidy-6.0 examples/diffusion_2d/*.cpp -- -I./src -I./examples/diffusion_2d -I/opt/intel/mkl/include
clang-tidy-6.0 examples/robertson/*.cpp -- -I./src -I./examples/robertson -I/opt/intel/mkl/include
