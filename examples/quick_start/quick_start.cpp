/*
 * Solves a trivial system of differential algebraic equations as a quick example:
 *
 * | x' = y
 * | y  = cos(t)
 *
 * Initial condition (t = 0):
 * | x = 0
 * | y = 1
 *
 * Analytic solution:
 * | x = sin(t)
 * | y = cos(t)
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

// Main dae-cpp header
#include <dae-cpp/solver.hpp>

// dae-cpp namespace
using namespace daecpp;

/*
 * Singular mass matrix in sparse format:
 *
 * M = |1 0|
 *     |0 0|
 */
struct MyMassMatrix
{
    /*
     * Defines the mass matrix `M` of the DAE system `M dx/dt = f`.
     * The mass matrix should be defined in sparse format (non-zero elements only) and can depend on time `t`.
     * Matrix `M` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &M, const double t)
    {
        M(0, 0, 1.0); // Row 0, column 0, non-zero element 1.0
    }
};

/*
 * The vector-function (RHS) of the problem:
 *
 * f(x,y,t) = | y
 *            | cos(t) - y
 */
struct MyRHS
{
    /*
     * Defines the RHS (vector function) `f` of the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the RHS vector `f`.
     * Vector `f` is already pre-allocated with `f.size() == x.size()`.
     *
     * For the given DAE system,
     * | x = x[0]
     * | y = x[1]
     */
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];          // y
        f[1] = cos(t) - x[1]; // cos(t) - y
    }
};

/*
 * (OPTIONAL) Analytic Jacobian in sparse format.
 * The DAE solver will use automatic (algorithmic) differentiation if analytic Jacobian is not provided.
 * However, providing the analytic Jacobian can significantly speed up the computation for large systems.
 *
 * Differentiating the RHS w.r.t. x[0] and x[1] gives the following Jacobian matrix:
 *
 * J = |0  1|
 *     |0 -1|
 */
struct MyJacobian
{
    /*
     * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the Jacobian matrix `J`.
     * Matrix `J` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(2);  // Pre-allocates memory for 2 non-zero elements (optional)
        J(0, 1, 1.0);  // Row 0, column 1, non-zero element 1.0
        J(1, 1, -1.0); // Row 1, column 1, non-zero element -1.0
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 */
int main()
{
    MyMassMatrix mass; // Mass matrix object
    MyRHS rhs;         // The vector-function object

    System my_system(mass, rhs); // Defines the DAE system object

    // Alternatively:
    // System my_system((MyMassMatrix()), (MyRHS())); // Note the parentheses

    // Class `System` is a wrapper for a more flexible function daecpp::solve(...).

    state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
    double t_end{1.0};     // Solution interval: t = [0, t_end]

    // To update the solver options, for example, to enable solver output during computation, use `opt` property:
    // my_system.opt.verbosity = verbosity::extra;

    my_system.solve(x0, t_end); // Solves the DAE system `my_system` with the given initial condition `x0` and time `t_end`

    // Alternatively:
    // my_system.solve(x0, t_end, MyJacobian()); // Add analytic Jacobian to speed up the computation (for big systems)
    // or
    // int status = my_system.solve({0, 1}, 1.0, MyJacobian()); // Note that `solve` can return the exit code

    my_system.sol.print(); // Prints solution on screen

    // Print the absolute error of the numerical solution `x` and `y` at time `t_end` using
    // `my_system.sol.x` - a vector that contains numerical solution `x` and `y` at times `my_system.sol.t`.
    std::cout << "Abs. error: " << my_system.sol.x.back()[0] - sin(t_end) << '\t' << my_system.sol.x.back()[1] - cos(t_end) << '\n';

    return 0;
}
