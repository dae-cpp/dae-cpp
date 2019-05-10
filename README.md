# dae-cpp

A simple but powerful C++ Differential Algebraic Equation (DAE) solver.

## What is dae-cpp

A C++ library that numerically solves a user-defined system of DAEs (an initial value problem). The system may contain both differential and algebraic equations and can be written in the following matrix-vector form:

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cmathbf%7BM%7D%5Cfrac%7Bd%5Cmathbf%7Bx%7D%7D%7Bdt%7D%3D%5Cmathbf%7Bf%7D%5Cleft%20%28%20%5Cmathbf%7Bx%7D%20%5Cright%20%29">
</p>

where Mass matrix **M** can be singular, and the RHS **f**(**x**) is a nonlinear function of a real vector **x** and time *t*.

For the numerical integration the solver uses implicit [BDF](https://en.wikipedia.org/wiki/Backward_differentiation_formula) (Backward Differentiation Formula) method of orders 1-6 (can be defined by a user) with adaptive time stepping.

### How does it work

BDF time stepper reduces the original DAE system to a system of nonlinear equations that the solver resolves using iterative [Newton root-finding algorithm](https://en.wikipedia.org/wiki/Newton%27s_method). Each Newton iteration a system of linear algebraic equations is solved using Parallel Direct Sparse Solver ([Intel MKL PARDISO](https://software.intel.com/en-us/mkl-developer-reference-c-intel-mkl-pardiso-parallel-direct-sparse-solver-interface)). The sparse solver performs 3 steps: reordering and symbolic factorisation of Jacobian matrix, then numerical factorisation, and then back substitution + iterative refinement. Finally, depending on the convergence rate of the Newton method and user-defined accuracy, the DAE solver may adjust the time step and initiate a new iteration in time.

### The main features of the solver

- May resolve DAE systems of 10<sup>8</sup> equations and even more (depending on the machine's RAM).
- A user may provide analytical Jacobian matrix for better performance or use built-in parallel function provided by the solver to estimate numerical Jacobian.
- Utilises all the available cores on the machine for better performance (this can be overridden by a user).
- Allows a user to adjust a lot of parameters related to the solution process in order to achieve better accuracy and performance. On the other hand, this is optional. Default values should work fine in most cases.
- The library provides a simple [C++ interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module for plotting.
- Easy-to-follow examples (see, for example, [perovskite.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.cpp)) to kick-start the user's project.

## Installation

This is a cross-platform software that should work well on both Linux (e.g. Ubuntu) and Windows. It should work under macOS as well (but not tested yet). The main library (DAE solver itself) and all examples have only one external dependency: [Intel Math Kernel Library](https://software.intel.com/en-us/mkl), a fast and very well optimised math library. So the first step in the installation process is to download and install Intel MKL: [Linux](https://software.intel.com/en-us/mkl/choose-download/linux), [Windows](https://software.intel.com/en-us/mkl/choose-download/windows), [macOS](https://software.intel.com/en-us/mkl/choose-download/macos).

An alternative and probably the most convenient way to download and install Intel MKL on Ubuntu (using APT Repository) is the following.

Install the GPG key for the repository:

```bash
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
```

Add the APT Repository:

```bash
sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
```

Update the list of packages and install the library:

```bash
sudo apt-get update
sudo apt-get install intel-mkl-2019.3-062
```

This will install Intel MKL 2019.3. The list of all available versions can be found [here](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo).

### Linux

On Linux make sure you have `git`, `cmake` and `g++` installed:

```bash
sudo apt-get install g++ cmake cmake-curses-gui git
```

Then download dae-cpp library:

```bash
git clone https://github.com/ikorotkin/dae-cpp.git
```

The easiest way to install the library and compile all examples is just to create a new directory and execute `cmake`:

```bash
cd dae-cpp
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/install/path ..
make
make install
```

where `/install/path` is the user-defined path where the package should be installed.

Note that `cmake` will try to find Intel MKL at its default location: `/opt/intel/mkl`. If the installation path is different, please provide it with the following `cmake` option: `-DDAE_MKL_DIR=/path/to/intel/mkl/root/dir`.

Instead of `cmake -DCMAKE_INSTALL_PREFIX=/install/path ..` you might consider to use `ccmake ..`. This is a GUI for `cmake` that will allow you to see all the options available before building the solver.

#### Test

Once installed with `DAE_TEST=ON` (it is ON by default), the solver can perform a quick self test:

```bash
ctest
```

During this test the solver will solve DAE systems from `examples` directory using both analytical and numerical Jacobians, and then compare the results with the reference solutions.

**TODO:**
- Describe more `cmake` options available (`DAE_SINGLE`, `DAE_FORTRAN_STYLE`, `DAE_TEST`, `DAE_LONG_INT`, `DAE_BUILD_EXAMPLES`)
- What installation directory contains

### Windows

Setting up the solver in Microsoft Visual Studio 2017. This has been tested but needs to be described...


## How to use

Please refer to [perovskite.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.cpp) as an example.

The main usage algorithm can be the following. Consider we have a system of DAEs written in a matrix-vector form, with some Mass matrix, RHS, and some initial conditions.

### Step 0. Include dae-cpp into the project

Include the solver's main header to the project. A shortcut to the solver's namespace (`daecpp`) can be added as well:

```cpp
#include "path/to/dae-cpp/include/solver.h"
namespace dae = daecpp;
```

### Step 1. Define the DAE parameters and initial state vector

For example, for *N* equations we should define the state vector with the size *N* and initialise it in accordance with the initial conditions:

```cpp
// State vector
dae::state_type x(N);

// Initial conditions
for(MKL_INT i = 0; i < N; i++)
{
    x[i] = 1.0;
}
```

We can access to each element of the state vector **x** as to `std::vector` from STL. Also note that instead of `int` or any other integer types we should use `MKL_INT` type. This gives us possibility to re-compile the project with `DAE_LONG_INT` option, so the code will work fine even for extremely huge systems (with *N* more than 10<sup>7</sup>).

### Step 2. Set up the RHS

Create MyRHS class that inherits the abstract `dae::RHS` class from dae-cpp library. The parent RHS class contains a pure virtual functor (operator `()`), that must be overridden in the child class. See, for example, [perovskite_RHS.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_RHS.cpp).

Once the RHS class is overridden, we can create an instance of the child class with some user-defined parameter container *p*:

```cpp
MyRHS rhs(p);
```

### Step 3. Set up the Mass matrix

Create MyMassMatrix class that inherits the abstract `dae::MassMatrix` class from dae-cpp library. Similar to the previous step, the parent MassMatrix class contains a pure virtual functor (operator `()`), that must be overridden in the child class. Refer to [perovskite_Mass.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_Mass.cpp) as an example. Note that the matrix should be defined in [three array sparse format](https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format).

Create an instance of the child MyMassMatrix class with the given size *N*:

```cpp
MyMassMatrix mass(N);
```

### Step 4. Set up Jacobian matrix

We can provide analytical Jacobian by overriding `dae::Jacobian` class from the dae-cpp library ([example](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_Jacobian.cpp)) or just use numerically estimated one (this may significantly slow down the computation for large *N*). If provided, analytical Jacobian matrix should be defined in [three array sparse format](https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format) similar to the Mass matrix.

If we don't provide analytical Jacobian we should estimate it with the given tolerance:

```cpp
dae::Jacobian jac(rhs, 1.0e-6);
```

Note that we should pass an instance of the user-defined RHS in order to estimate numerical Jacobian.

### Step 5. Set the solver options

The solver has a lot of options. They all have some default values (defined in [solver_options.h](https://github.com/ikorotkin/dae-cpp/blob/master/src/solver_options.h)) but they can be overridden by a user:

```cpp
// Create an instance of the solver options and update some of the solver
// parameters defined in solver_options.h
dae::SolverOptions opt;

// For example, let's change the absolute tolerance
opt.atol = 1.0e-6;
```

### Step 6. Solve the system

Now we are ready to create an instance of the solver with particular RHS, Mass matrix, Jacobian and the solver options, and then start the solver:

```cpp
dae::Solver solve(rhs, jac, mass, opt, t1);
solve(x);
```

Here *t1* is the integration time (0 < *t* < *t1*).

Solution at time *t1* will be written into vector **x** (initial conditions will be overwritten). That's it!

### Step 7 (optional). Plot results

Solution can be visualised using a simple [C++ interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module. For example, if `python`, `numpy` and `matplotlib` are installed, the [perovskite](https://github.com/ikorotkin/dae-cpp/tree/master/examples/perovskite) example will produce the following plot:

<p align="center">
  <img src="http://korotkin.ru/public/figure.png">
</p>

## Licensing

- dae-cpp is fully open source under [MIT license](https://github.com/ikorotkin/dae-cpp/blob/master/LICENSE).
- Intel MKL is free for use and redistribution under [Intel Simplified Software License](https://software.intel.com/en-us/license/intel-simplified-software-license).
