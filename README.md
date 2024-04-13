# dae-cpp

A simple but powerful header-only C++ solver for systems of Differential-Algebraic Equations (DAE).

## What is dae-cpp

$$\mathbf{M}(t) \frac{\mathrm{d}\boldsymbol{x}}{\mathrm{d}t} = \boldsymbol{f}(\boldsymbol{x}, t),$$

where...

The initial condition $\boldsymbol{x}\rvert_{t=0} = \boldsymbol{x_0}$... Interval $t \in [0, t_\mathrm{end}]$.

### How does it work

TODO: Quasi-Newton method...

### The main features of the solver

## Installation

Header only, no need to install, just copy `dae-cpp`, `Eigen`, and `autodiff` folders into your project.

Examples and tests can be compiled using CMake.

## Testing

If you already cloned the project without `--recurse-submodules`, you can initialize and update `googletest` submodule by running

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

## Quick start

TODO: This is still work in progress.

Trivial DAE system:

```math
\left\{
    \begin{alignedat}{3}
        \dot x & = y, \\
        y & = \cos(t).
    \end{alignedat}
\right.
```

Initial condition:

```math
\left\{
    \begin{alignedat}{3}
        x\rvert_{t=0} & = 0, \\
        y\rvert_{t=0} & = 1.
    \end{alignedat}
\right.
```

Analytic solution:

```math
\left\{
    \begin{alignedat}{3}
        x(t) & = \sin(t), \\
        y(t) & = \cos(t).
    \end{alignedat}
\right.
```

 See [example](https://github.com/dae-cpp/dae-cpp/blob/master/examples/quick_start/quick_start.cpp).

### Step 0. Include dae-cpp header into the project

```cpp
#include <dae-cpp/solver.hpp>
```

### Step 1. Define the mass matrix of the system

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

### Step 2. Define the vector-function (RHS) of the problem

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

my_system.solve(x0, t); // Solves the system with the given initial condition `x0` and time `t`
```

or simply

```cpp
my_system.solve({0, 1}, 1.0);
```

Solution vector of vectors `x` and the corresponding vector of times `t` are stored in `my_system.sol.x` and `my_system.sol.t`, respectively.

### (Optional) Step 5. Define the Jacobian matrix to boost the computation speed

$$
\mathbf{J} =
\begin{vmatrix}
0 & 1 \\
0 & -1
\end{vmatrix}.
$$

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

0. Create a GitHub issue if you want to suggest or discuss your changes.
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
