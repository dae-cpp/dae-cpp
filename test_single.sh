rm -r build/
mkdir build
cd build
cmake -DDAE_SINGLE=ON ..
make -j 4
ctest
