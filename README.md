# dae-cpp

![tests](https://github.com/dae-cpp/dae-cpp/actions/workflows/cmake-multi-platform.yml/badge.svg)
![version](https://img.shields.io/badge/version-2.0.0-blue)
[![Static Badge](https://img.shields.io/badge/Documentation-8A2BE2?logo=githubpages&logoColor=fff&style=flat)](https://dae-cpp.github.io/)

**A simple but powerful header-only C++ solver for [systems of Differential-Algebraic Equations](https://en.wikipedia.org/wiki/Differential-algebraic_system_of_equations) (DAE).**

**NOTE:** `dae-cpp` has been redesigned and there were breaking changes between `v1.x` and `v2.x`. If your project still relies on the old `dae-cpp` (`v1.x`), it is archived in the [legacy](https://github.com/dae-cpp/dae-cpp/tree/legacy) branch. For the new version (`v2.x`), see [Documentation](https://dae-cpp.github.io/) and the notes below.

## What is `dae-cpp`

`dae-cpp` is a cross-platform, header-only C++17 library for solving stiff systems of DAEs (an initial value problem). DAE systems can contain both differential and algebraic equations and can be written in the following matrix-vector form:

$$\mathbf{M}(t) \frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = \mathbf{f}(\mathbf{x}, t),$$

to be solved in the interval $`t \in [0, t_\mathrm{end}]`$ with the initial condition $`\mathbf{x}\rvert_{t=0} = \mathbf{x}_0`$. Here $`\mathbf{M}(t)`$ is the mass matrix (can depend on time), $`\mathbf{x}(t)`$ is the state vector, and $`\mathbf{f}(\mathbf{x}, t)`$ is the (nonlinear) vector function of the state vector $`\mathbf{x}`$ and time $t$.

### How does it work

The DAE solver uses implicit Backward Differentiation Formulae (BDF) of orders I-IV with adaptive time stepping. Every time step, the BDF integrator reduces the original DAE system to a system of nonlinear equations, which is solved using iterative [Quasi-Newton](https://en.wikipedia.org/wiki/Quasi-Newton_method) root-finding algorithm. The Quasi-Newton method reduces the problem further down to a system of linear equations, which is solved using [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), a versatile and fast C++ template library for linear algebra.
Eigen's sparse solver performs two steps: factorization (decomposition) of the Jacobian matrix and the linear system solving itself. This gives us the numerical solution of the entire DAE system at the current time step. Finally, depending on the convergence rate of the Quasi-Newton method, variability of the solution, and user-defined accuracy, the DAE solver adjusts the time step size and initiates a new iteration in time.

### The main features of the solver

- Header only, no pre-compilation required.
- Uses [automatic](https://en.wikipedia.org/wiki/Automatic_differentiation) (algorithmic, exact) differentiation ([autodiff](https://autodiff.github.io/) package) to compute the Jacobian matrix, if it is not provided by the user.
- Fourth-order variable-step implicit BDF time integrator that preserves accuracy even when the time step rapidly changes.
- A very flexible and customizable variable time stepping algorithm based on the solution stability and variability.
- Mass matrix can be non-static (can depend on time) and it can be singular.
- The library is extremely easy to use. A simple DAE can be set up using just a few lines of code (see [Quick Start](#quick-start) example below).

## Installation

This library is header only, no need to install, just copy `dae-cpp`, `Eigen`, and `autodiff` folders into your project.

Examples and tests can be compiled using CMake (see [Testing](#testing)).

## Testing

If you already have cloned the project without `--recurse-submodules` option, you can initialize and update `googletest` submodule by running

```bash
git submodule update --init
```

Then build and run the tests:

```bash
mkdir build && cd build
cmake ..
make
ctest
```

## Documentation, examples, and CHANGELOG

- For more information about the solver, please refer to the [Documentation](https://dae-cpp.github.io/) pages.
- Ready to use examples are given in the [examples](https://github.com/dae-cpp/dae-cpp/tree/master/examples) directory of this repository.
- All notable user-facing changes to this project are documented in the [CHANGELOG](https://dae-cpp.github.io/CHANGELOG.html).

## Quick Start

Consider the following (trivial) DAE system as a quick example:

```math
\left\{
    \begin{alignedat}{3}
        \dot x & = y, \\
        y & = \cos(t),
    \end{alignedat}
\right.
```

with the initial condition:

```math
\left\{
    \begin{alignedat}{3}
        x\rvert_{t=0} & = 0, \\
        y\rvert_{t=0} & = 1.
    \end{alignedat}
\right.
```

This system contains one simple differential equation and one algebraic equation. The analytic solution is the following:

```math
\left\{
    \begin{alignedat}{3}
        x(t) & = \sin(t), \\
        y(t) & = \cos(t).
    \end{alignedat}
\right.
```

Below is a simplified procedure of defining and solving the DAE system using `dae-cpp`.

### Step 0. Include `dae-cpp` header into the project

```cpp
#include <dae-cpp/solver.hpp>
```

Optionally, add `daecpp` namespace:

```cpp
using namespace daecpp;
```

### Step 1. Define the mass matrix of the system

Tha mass matrix contains only one non-zero element:

$$
\mathbf{M} =
\begin{vmatrix}
1 & 0 \\
0 & 0
\end{vmatrix}.
$$

```cpp
struct MyMassMatrix
{
    void operator()(sparse_matrix &M, const double t)
    {
        M(0, 0, 1.0); // Row 0, column 0, non-zero element 1.0
    }
};
```

### Step 2. Define the vector function (RHS) of the system

```cpp
struct MyRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];          // y
        f[1] = cos(t) - x[1]; // cos(t) - y
    }
};
```

### Step 3. Set up the DAE system

```cpp
MyMassMatrix mass; // Mass matrix object
MyRHS rhs;         // Vector-function object

System my_system(mass, rhs); // Defines the DAE system object
```

### Step 4. Solve the system

```cpp
state_vector x0{0, 1}; // The initial state vector (initial condition)
double t{1.0};         // The integration interval: t = [0, 1.0]

my_system.solve(x0, t); // Solves the system with initial condition `x0` and time `t`
```

or simply

```cpp
my_system.solve({0, 1}, 1.0);
```

Solution vector of vectors `x` and the corresponding vector of times `t` will be stored in `my_system.sol.x` and `my_system.sol.t`, respectively.

The entire source code is provided in the [Quick Start example](https://github.com/dae-cpp/dae-cpp/blob/master/examples/quick_start/quick_start.cpp).

For more information, refer to the [Documentation](https://dae-cpp.github.io/).

### (Optional) Step 5. Define the Jacobian matrix to boost the computation speed

Differentiating the RHS w.r.t. $x$ and $y$ gives the following Jacobian matrix:

$$
\mathbf{J} =
\begin{vmatrix}
0 & 1 \\
0 & -1
\end{vmatrix}.
$$

This matrix can be defined in the code as

```cpp
struct MyJacobian
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(2);  // Pre-allocates memory for 2 non-zero elements (optional)
        J(0, 1, 1.0);  // Row 0, column 1, non-zero element 1.0
        J(1, 1, -1.0); // Row 1, column 1, non-zero element -1.0
    }
};
```

Then add user-defined Jacobian to the DAE system definition:

```cpp
System my_system(mass, rhs, MyJacobian()); // Defines the DAE system with Jacobian
```

### (Optional) Step 6. Tweak the solver options

For example, restrict the maximum time step:

```cpp
my_system.opt.dt_max = 0.1;   // Update `dt_max`
my_system.solve({0, 1}, 1.0); // Restart the computation
```

## Contribution and Feedback

Thank you for considering contributing to the project! Whether you're an experienced developer or just starting out, your ideas and improvements can make this project even better. No contribution is too small!

### How to contribute

0. Create a [GitHub issue](https://github.com/dae-cpp/dae-cpp/issues) if you want to suggest or discuss your changes.
1. Fork the repository and clone it to your local machine.
2. Create a new branch for your contributions.
3. Make your changes and ensure they adhere to our coding standards.
4. Add tests and examples (if relevant), test your changes thoroughly.
5. Submit a pull request with a clear description of your changes and why they are beneficial.

### Feedback

Feel free to create a [GitHub issue](https://github.com/dae-cpp/dae-cpp/issues) for any questions, suggestions, or feedback you may have.

## Licensing

- [dae-cpp](https://github.com/dae-cpp/dae-cpp) is licensed under the [MIT](https://github.com/dae-cpp/dae-cpp/blob/master/LICENSE) license.
- [autodiff](https://github.com/autodiff/autodiff) is licensed under the [MIT](https://github.com/autodiff/autodiff/blob/main/LICENSE) license.
- [Eigen](https://eigen.tuxfamily.org/) is licensed under the [MPL2](https://www.mozilla.org/en-US/MPL/2.0/) license.
