#!/bin/bash

# Build a static library:
mkdir build_static_lib
cd build_static_lib
g++ -c -O3 -Wall -std=c++11 -m64 -fopenmp ../src/*.cpp -I/opt/intel/mkl/include -I../src -I../src/external
ar rcs libdaecpp.a *.o
rm *.o

# Compile examples
g++ -O3 -Wall -std=c++11 -m64 -fopenmp ../examples/perovskite/*.cpp -o perovskite -I/opt/intel/mkl/include -I../examples/perovskite -I../src -I../src/external -L/opt/intel/mkl/lib/intel64 -L. -ldaecpp -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
g++ -O3 -Wall -std=c++11 -m64 -fopenmp ../examples/diffusion_2d/*.cpp -o diffusion_2d -I/opt/intel/mkl/include -I../examples/diffusion_2d -I../src -I../src/external -L/opt/intel/mkl/lib/intel64 -L. -ldaecpp -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
g++ -O3 -Wall -std=c++11 -m64 -fopenmp ../examples/robertson/*.cpp -o robertson -I/opt/intel/mkl/include -I../examples/robertson -I../src -I../src/external -L/opt/intel/mkl/lib/intel64 -L. -ldaecpp -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
g++ -O3 -Wall -std=c++11 -m64 -fopenmp ../examples/simple_dae/*.cpp -o simple_dae -I/opt/intel/mkl/include -I../examples/simple_dae -I../src -I../src/external -L/opt/intel/mkl/lib/intel64 -L. -ldaecpp -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
g++ -O3 -Wall -std=c++11 -m64 -fopenmp ../examples/two_bodies/*.cpp -o two_bodies -I/opt/intel/mkl/include -I../examples/two_bodies -I../src -I../src/external -L/opt/intel/mkl/lib/intel64 -L. -ldaecpp -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

echo 'How to run an example using 4 cores:'
echo 'cd build_static_lib/'
echo 'source ../set_MKL_env'
echo 'OMP_NUM_THREADS=4 ./perovskite'
