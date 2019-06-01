# dae-cpp

A simple but powerful C++ solver for Differential Algebraic Equation (DAE) systems.

## What is dae-cpp

A cross-platform, parallel C++ library that numerically solves a user-defined system of DAEs (an initial value problem). The system may contain both differential and algebraic equations and can be written in the following matrix-vector form:

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cmathbf%7BM%7D%5Cfrac%7Bd%5Cmathbf%7Bx%7D%7D%7Bdt%7D%3D%5Cmathbf%7Bf%7D%5Cleft%20%28%20%5Cmathbf%7Bx%7D%20%5Cright%20%29">
</p>

where Mass matrix **M** can be singular, and the RHS **f**(**x**) is a nonlinear function of a real vector **x** and time *t*.

For the numerical integration the solver uses implicit [BDF](https://en.wikipedia.org/wiki/Backward_differentiation_formula) (Backward Differentiation Formula) method of orders 1-6 (can be defined by a user) with adaptive time stepping.

### How does it work

BDF time stepper reduces the original DAE system to a system of nonlinear equations that the solver resolves using iterative [Newton root-finding algorithm](https://en.wikipedia.org/wiki/Newton%27s_method). Each Newton iteration a system of linear algebraic equations is solved using Parallel Direct Sparse Solver ([Intel MKL PARDISO](https://software.intel.com/en-us/mkl-developer-reference-c-intel-mkl-pardiso-parallel-direct-sparse-solver-interface)). The sparse solver performs 3 steps: reordering and symbolic factorisation of Jacobian matrix, then numerical factorisation, and then back substitution + iterative refinement. Finally, depending on the convergence rate of the Newton method, variability of the solution and user-defined accuracy, the DAE solver may adjust the time step and initiate a new iteration in time.

### The main features of the solver

- May resolve DAE systems of 10<sup>8</sup> equations and even more (depending on the machine's RAM).
- A user may provide analytical Jacobian matrix for better performance or use built-in parallel function provided by the solver to estimate numerical Jacobian.
- Utilises all the available cores on the machine for better performance (this can be overridden by a user).
- Allows a user to adjust most of the parameters related to the solution process in order to achieve better accuracy and performance. On the other hand, this is optional. Default values should work fine in most cases.
- The library provides a simple [C++ interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module for plotting.
- Easy-to-follow examples (see, for example, [perovskite.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.cpp) or [diffusion_2d.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/diffusion_2d/diffusion_2d.cpp)) to kick-start the user's project.

## Installation

This is a cross-platform software that should work on both Linux (e.g. Ubuntu) and Windows. It should work under macOS as well (but not tested yet). The main library (DAE solver itself) and all examples have only one external dependency: [Intel Math Kernel Library](https://software.intel.com/en-us/mkl), a fast and very well optimised math library. So the first step in the installation process is to download and install Intel MKL: [Linux](https://software.intel.com/en-us/mkl/choose-download/linux), [Windows](https://software.intel.com/en-us/mkl/choose-download/windows), [macOS](https://software.intel.com/en-us/mkl/choose-download/macos).

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

The easiest way to install the library and compile all examples is just to create a new directory and then execute `cmake` and `make`:

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

Instead of `cmake -DCMAKE_INSTALL_PREFIX=/install/path ..` you might consider using `ccmake ..`, a GUI for `cmake` that will allow you to see all the options available before building the solver.

#### Test the solver

The DAE solver can perform a quick self test. To build the test, dae-cpp should be installed with `DAE_TEST=ON` option (it is ON by default). To start the test, from the build directory execute `ctest`:

```bash
ctest
```

During this test the solver will solve DAE systems from [examples](https://github.com/ikorotkin/dae-cpp/tree/master/examples) directory using analytical (if available) and numerical Jacobians, and then compare the results with the reference solutions.

#### More building options

- `DAE_LONG_INT` - Use long integer representation for huge systems (more than ~10<sup>7</sup> equations). This option is OFF by default. For relatively small systems it is recommended to leave it OFF.
- `DAE_FORTRAN_STYLE` - If ON, the matrices will be defined using FORTRAN style (one-based indexing of columns and rows). By default it is OFF (zero-based indexing).
- `DAE_SINGLE` - If ON, the single precision will be used in the solver instead of double. Single precision may ruin the accuracy. It is highly recommended to leave this option OFF. This option exists for the future compatibility with CUDA implementations of the solver.
- `DAE_BUILD_EXAMPLES` - Build all the examples, ON by default.
- `DAE_TEST` - Build automatic solver test, ON by default. The test can be executed by the command `ctest` from the building directory.

### Windows

Setting up the solver in Microsoft Visual Studio 2017. This has been tested but needs to be described...


## How to use

Please refer to [perovskite.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.cpp) or [diffusion_2d](https://github.com/ikorotkin/dae-cpp/blob/master/examples/diffusion_2d/diffusion_2d.cpp) as an example.

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

We can get access to each element of the state vector **x** as to `std::vector` from STL. Also note that instead of `int` or any other integer types we should use `MKL_INT` type. This gives us possibility to re-compile the project with `DAE_LONG_INT` option, so the code will work fine even for extremely huge systems (with *N* more than 10<sup>7</sup>).

### Step 2. Set up the RHS

Create MyRHS class that inherits the abstract `daecpp::RHS` class from dae-cpp library. The parent RHS class contains a pure virtual functor (operator `()`), that must be overridden in the child class. See, for example, [perovskite_RHS.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_RHS.cpp) or [diffusion_2d_RHS.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/diffusion_2d/diffusion_2d_RHS.cpp).

Once the RHS class is overridden, we can create an instance of the child class with some user-defined parameter container *p*:

```cpp
MyRHS rhs(p);
```

In the child MyRHS class the user can also override `stop_condition` virtual function. By default (if not overridden) the function always returns `false`. The user may override this behaviour and set up one or several stop conditions for the solver depending on the solution **x** at the current time *t*. As soon as the function returns `true`, the solver will finalise the current time step and return the current solution. A trivial example of the stop condition function can be found in [perovskite_RHS.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_RHS.cpp).

### Step 3. Set up the Mass matrix

Create MyMassMatrix class that inherits the abstract `daecpp::MassMatrix` class from dae-cpp library. Similar to the previous step, the parent MassMatrix class contains a pure virtual functor (operator `()`), that must be overridden in the child class. Refer to [perovskite_Mass.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_Mass.cpp) as an example. Note that the matrix should be defined in [three array sparse format](https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format).

Create an instance of the child MyMassMatrix class with the given size *N*:

```cpp
MyMassMatrix mass(N);
```

If the Mass matrix is a simple identity matrix, one can use `daecpp::MassMatrixIdentity` class from dae-cpp library instead of inheriting `daecpp::MassMatrix`. This will create identity Mass matrix in sparse format with the given size *N*:

```cpp
dae::MassMatrixIdentity mass(N);
```

### Step 4. Set up Jacobian matrix

We can provide analytical Jacobian by overriding `daecpp::Jacobian` class from the dae-cpp library ([example](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_Jacobian.cpp)) or just use numerically estimated one (this may significantly slow down the computation for large *N*). If provided, analytical Jacobian matrix should be defined in [three array sparse format](https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format) similar to the Mass matrix.

If we don't provide analytical Jacobian we should estimate it with the given tolerance:

```cpp
dae::Jacobian jac(rhs, 1.0e-6);
```

Note that we should pass an instance of the user-defined RHS in order to estimate numerical Jacobian.

### Step 5. Set the solver options

The solver has lots of options related to the solution process. They all have some default values (defined in [solver_options.h](https://github.com/ikorotkin/dae-cpp/blob/master/src/solver_options.h)) but they can be overridden by a user:

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
dae::Solver solve(rhs, jac, mass, opt);
solve(x, t1);
```

Here *t*<sub>1</sub> is the integration time (0 < *t* < *t*<sub>1</sub>), and **x** is the initial condition vector defined above.

Solution at time *t*<sub>1</sub> will be written into vector **x** (initial conditions will be overwritten). That's it!

Note that in order to get intermediate solutions at times *t*<sub>a</sub>, *t*<sub>b</sub>, *t*<sub>c</sub>, etc. (0 < *t*<sub>a</sub> < *t*<sub>b</sub> < *t*<sub>c</sub> < ... < *t*<sub>1</sub>), for example, for plotting, one can call the solver at the given times:

```cpp
solve(x, t_a);  // solves the system in the interval [0; t_a] and stores the solution in x
solve(x, t_b);  // continues solving in the interval [t_a; t_b], replaces the solution in x
solve(x, t_c);  // continues solving in the interval [t_b; t_c], replaces the solution in x
                // ...
solve(x, t1);   // continues solving in the interval [t_c; t1], stores the final solution at time t1 in x
```

Every call the solver will take the previous solution **x** (if available from the previous call) and overwrite it with a new one at the given time.

But a proper (and more efficient) way to get intermediate results is to override `virtual void observer(daecpp::state_type &x, const double t)` function from `daecpp::Solver` class. This observer function receives solution vector **x** and the current time *t* every time step and allows a user to get access to the solution during the solving process at each time layer.

### Step 7 (optional). Plot results

Solution can be visualised using a simple [C++ interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module. For example, if `python`, `numpy` and `matplotlib` are installed, the [perovskite](https://github.com/ikorotkin/dae-cpp/tree/master/examples/perovskite) example will produce the following plot:

<p align="center">
  <img src="http://korotkin.ru/public/perovskite.png">
</p>

Note that by default the plotting is switched off in the examples, but the plotting-related code can be activated using `#define PLOTTING` at the very beginning of each example. Activating the plotting refers to `matplotlibcpp.h` header located in `src/external/matplotlib-cpp/` directory.

The second example, [diffusion_2d](https://github.com/ikorotkin/dae-cpp/tree/master/examples/diffusion_2d) will produce a two-dimensional Gaussian function, a solution of two-dimensional diffusion problem with an instantaneous point source in the middle of the plane:

<p align="center">
  <img src="http://korotkin.ru/public/diffusion_2d.png">
</p>

## Licensing

- dae-cpp is fully open source under [MIT license](https://github.com/ikorotkin/dae-cpp/blob/master/LICENSE).
- Intel MKL is free for use and redistribution under [Intel Simplified Software License](https://software.intel.com/en-us/license/intel-simplified-software-license).
