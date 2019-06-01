rm -r build_single/
mkdir build_single
cd build_single
cmake -DDAE_SINGLE=ON ..
make -j 4
ctest
cd ..

rm -r build/
mkdir build
cd build
cmake ..
make -j 4
ctest
