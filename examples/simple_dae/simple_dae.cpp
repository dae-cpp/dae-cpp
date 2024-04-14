/*
 * Solves another simple system of differential algebraic equations.
 * This example introduces Solution Manager, that will work as solution observer.
 *
 * The system:
 * | x' = y
 * | x*x + y*y = 1
 *
 * Initial condition (t = 0):
 * | x = 0
 * | y = 1
 *
 * Analytic solution for t = [0, pi/2]:
 * | x = sin(t)
 * | y = cos(t)
 *
 * There is also another, trivial solution x = 1, y = 0, t > pi/2,
 * so we will be solving the system in the interval t = [0, pi/2].
 *
 * Every time step, we will check that:
 *   1. x*x + y*y = 1
 *   2. x(t) = sin(t)
 * and print the absolute errors.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

// Main dae-cpp header
#include <dae-cpp/solver.hpp>

// dae-cpp namespace
using namespace daecpp;

constexpr double pi{3.14159265358979323846};

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
 *            | x*x + y*y - 1
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
        f[0] = x[1];                          // y
        f[1] = x[0] * x[0] + x[1] * x[1] - 1; // x*x + y*y - 1
    }
};

/*
 * (OPTIONAL) Analytic Jacobian in sparse format.
 * The DAE solver will use automatic (algorithmic) differentiation if analytic Jacobian is not provided.
 * However, providing the analytic Jacobian can significantly speed up the computation for large systems.
 *
 * Differentiating the RHS w.r.t. x[0] and x[1] gives the following Jacobian matrix:
 *
 * J = |0  1 |
 *     |2x 2y|
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
        J.reserve(3);        // Pre-allocates memory for 3 non-zero elements (optional)
        J(0, 1, 1.0);        // Row 0, column 1, non-zero element `1.0`
        J(1, 0, 2.0 * x[0]); // Row 1, column 0, non-zero element `2x`
        J(1, 1, 2.0 * x[1]); // Row 1, column 1, non-zero element `2y`
    }
};

/*
 * User-defined Solution Manager to post-process solution every time step.
 * In this example, it works as a passive observer that saves solution error during computation.
 */
class MyObserver
{
    state_vector &_err1, &_err2; // Vectors of the solution errors

public:
    MyObserver(state_vector &err1, state_vector &err2) : _err1(err1), _err2(err2) {}

    /*
     * Solution Manager functor will be called every time step providing the time `t` and
     * the corresponding solution `x` for further post-processing.
     * If the functor returns an integer != 0 (`true`), the computation will immediately stop.
     */
    virtual int operator()(const state_vector &x, const double t)
    {
        double e1 = std::abs(x[0] * x[0] + x[1] * x[1] - 1.0); // Error 1
        double e2 = std::abs(std::sin(t) - x[0]);              // Error 2

        // Every time step, print time t, x(t), y(t), Error 1, Error 2
        // std::cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << e1 << '\t' << e2 << '\n';

        _err1.push_back(e1);
        _err2.push_back(e2);

        return 0;
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 */
int main()
{
    state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
    double t_end{pi / 2};  // Solution interval: t = [0, t_end]

    state_vector error1, error2; // Absolute errors

    int status = solve(MyMassMatrix(), MyRHS(), x0, t_end, MyObserver(error1, error2)); // Solves the DAE system

    // To add analytic Jacobian:
    // int status = solve(MyMassMatrix(), MyRHS(), MyJacobian(), x0, t_end, MyObserver(error1, error2));

    // To add user-defined solver options:
    // SolverOptions opt;
    // opt.verbosity = verbosity::extra;
    // opt.dt_init = 0.001;
    // int status = solve(MyMassMatrix(), MyRHS(), MyJacobian(), x0, t_end, MyObserver(error1, error2), opt);

    // Print both errors at the very end of computation
    std::cout << "Abs. error: " << error1.back() << '\t' << error2.back() << '\n';

    return 0;
}
