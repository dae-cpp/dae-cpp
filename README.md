# dae-cpp

A simple but powerful C++ Differential Algebraic Equation (DAE) solver.

## What is dae-cpp

A set of C++ routines that numerically solve a user-defined system of DAEs. Such system may contain both differential and algebraic equations and can be written in the following matrix-vector form:

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cmathbf%7BM%7D%5Cfrac%7Bd%5Cmathbf%7Bx%7D%7D%7Bdt%7D%3D%5Cmathbf%7Bf%7D%5Cleft%20%28%20%5Cmathbf%7Bx%7D%20%5Cright%20%29">
</p>

where Mass matrix **M** can be singular, and RHS **f**(**x**) is a nonlinear function of a real vector **x** and time *t*.

For the numerical integration the solver uses implicit [BDF](https://en.wikipedia.org/wiki/Backward_differentiation_formula) (Backward Differentiation Formula) method of orders 1-6 (can be defined by a user) with adaptive time stepping. BDF time stepper reduces the original DAE system to a system of nonlinear equations that the solver resolves using iterative Newton root-finding [algorithm](https://en.wikipedia.org/wiki/Newton%27s_method). Each Newton iteration a system of linear algebraic equations is solved using Parallel Direct Sparse Solver ([Intel MKL PARDISO](https://software.intel.com/en-us/mkl-developer-reference-c-intel-mkl-pardiso-parallel-direct-sparse-solver-interface)). The sparse solver performs 3 steps: reordering and symbolic factorisation of Jacobian matrix, then numerical factorisation, and then back substitution + iterative refinement. Finally, depending on the convergence rate of the Newton method and user-defined accuracy, the DAE solver may adjust the time step and initiate a new step in time.

The main features of the solver:

- May resolve DAE systems of 10<sup>8</sup> equations and even more (depending on the machine's RAM).
- A user may provide analytical Jacobian matrix for better performance or use built-in parallel function provided by the solver to estimate numerical Jacobian.
- Utilises all the available cores on the machine for better performance (this can be overridden by a user).
- Allows a user to adjust a lot of parameters related to the solution process in order to achieve better accuracy and performance. On the other hand, this is optional. Default values should work fine in most cases.
- The library provides a simple C++ [interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module for plotting.
- Easy-to-follow examples to kick-start the user's project.

## Installation

This is a cross-platform software that should work well on both Linux (e.g. Ubuntu) and Windows. It should work under macOS as well (but not tested yet). The main library (DAE solver itself) and all examples have only one external dependency: [Intel Math Kernel Library](https://software.intel.com/en-us/mkl), a fast and very well optimised math library. So the first step in the installation process is to download and install Intel MKL: [Linux](https://software.intel.com/en-us/mkl/choose-download/linux), [Windows](https://software.intel.com/en-us/mkl/choose-download/windows), [macOS](https://software.intel.com/en-us/mkl/choose-download/macos).

An alternative way to download and install Intel MKL on Ubuntu (via `apt-get`)is the following:

```bash
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list
sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
sudo apt-get update
sudo apt-get install intel-mkl-2019.3-062
```

In order to uninstall:

```bash
sudo apt-get autoremove intel-mkl-2019.3-062
```

### Linux

On Linux make sure you have `git`, `cmake` and `g++` installed:

```bash
sudo apt-get install g++ cmake git
```

Then download dae-cpp library:

```bash
git clone https://github.com/ikorotkin/dae-cpp.git
```

The easiest way to install the library and compile all examples is just to create a new directory and execute `cmake`:

```bash
cd dae-cpp
mkdir build
cmake -DCMAKE_INSTALL_PREFIX=/install/path .. 
make
make install
```

where `/install/path` is a user-defined path where the package should be installed.

Note that `cmake` will try to find Intel MKL at its default location: `/opt/intel/mkl`. If the installation path is different, please provide it with the following `cmake` option: `-DDAE_MKL_DIR=/path/to/intel/mkl/root/dir`.

TODO:
- Describe more `cmake` options available (DAE_SINGLE, DAE_FORTRAN_STYLE, DAE_TEST, etc.)
- Mention `ccmake`
- Mention about tests (`ctest`)
- Windows installation (Microsoft Visual Studio 2017)
- Plotting with `matplotlib`
- What installation dir contains

## How to use

Please refer to [perovskite.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.cpp) while this section is in progress.
