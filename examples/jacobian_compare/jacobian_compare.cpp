/*
 * In this example, we do not solve any DAEs. Instead, we define a simple vector function
 * (for demonstration purposes), the corresponding Jacobian matrix to boost the calculation performance,
 * and then we use a built-in helper class `JacobianCompare` to compare our manually derived Jacobian
 * with the one computed algorithmically from the vector function.
 *
 * We will make a few mistakes in the analytic Jacobian on purpose to see what information
 * `JacobianComapre` can provide.
 *
 * Note that `JacobianCompare` can work with Jacobians derived from the given Jacobian shapes as well.
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

/*
 * The vector-function (RHS) of the problem.
 */
struct MyRHS
{
    /*
     * Defines the RHS (vector function) `f` of the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the RHS vector `f`.
     * Vector `f` is already pre-allocated with `f.size() == x.size()`.
     */
    void operator()(state_type &f, const state_type &x, const double t)
    {
        // `dual_type` for automatic differentiation
        dual_type u = x[0];
        dual_type v = x[1];
        dual_type w = x[2];
        dual_type z = x[3];

        f[0] = z * 0.5;
        f[1] = u * u + v * v + w * w;
        f[2] = sin(u);
        f[3] = exp(-v) + t * z;
    }
};

/*
 * Analytic Jacobian in sparse format.
 * Obtained by differentiation the RHS w.r.t. vector `x`.
 * Contains a few errors:
 *     - one missing element
 *     - one element computed incorrectly
 *     - one element is in wrong place (incorrect `i` and `j` indices)
 */
struct MyJacobianBad
{
    /*
     * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the Jacobian matrix `J`.
     * Matrix `J` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        float_type u = x[0];
        float_type v = x[1];
        float_type w = x[2];
        float_type z = x[3];

        J.reserve(6); // Pre-allocates memory for 6 non-zero elements

        J(0, 3, 1.0);      // Row 0, column 3 - error, should be `0.5`
        J(1, 0, 2.0 * u);  // Row 1, column 0
        J(1, 1, 2.0 * v);  // Row 1, column 1
        J(1, 2, 2.0 * w);  // Row 1, column 2
        J(2, 1, cos(u));   // Row 2, column 1 - error, should be column 0
        J(3, 1, -exp(-v)); // Row 3, column 1
        // Forgot to differentiate f[3] w.r.t. `z` - missing element (3, 3)
    }
};

/*
 * Correct analytic Jacobian in sparse format.
 * Obtained by differentiation the RHS w.r.t. vector `x`.
 */
struct MyJacobianGood
{
    /*
     * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the Jacobian matrix `J`.
     * Matrix `J` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        float_type u = x[0];
        float_type v = x[1];
        float_type w = x[2];
        float_type z = x[3];

        J.reserve(7); // Pre-allocates memory for 7 non-zero elements

        J(0, 3, 0.5);      // Row 0, column 3
        J(1, 0, 2.0 * u);  // Row 1, column 0
        J(1, 1, 2.0 * v);  // Row 1, column 1
        J(1, 2, 2.0 * w);  // Row 1, column 2
        J(2, 0, cos(u));   // Row 2, column 1
        J(3, 1, -exp(-v)); // Row 3, column 1
        J(3, 3, t);        // Row 3, column 3
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 */
int main()
{
    // Fill vectors x at which the Jacobian matrix will be tested
    state_vector x0 = {1.0, 2.0, 3.0, 4.0};
    state_vector x1 = {-0.5, -0.1, 0.5, 0.1};

    // Check "bad" Jacobian at `x0` and `x1` and times t = 1 and 2
    {
        auto jac_comparison = JacobianCompare(MyJacobianBad(), MyRHS());

        auto N_errors_1 = jac_comparison(x0, 1.0);
        auto N_errors_2 = jac_comparison(x1, 2.0);

        std::cout << "Tested \"bad\" Jacobian at x0, found " << N_errors_1 << " errors.\n";
        std::cout << "Tested \"bad\" Jacobian at x1, found " << N_errors_2 << " errors.\n\n";
    }

    // This will show the following summary:
    //
    // Jacobian matrix comparison summary at time t = 1:
    // -- Found 4 difference(s) compared to the automatic (reference) Jacobian:
    // ----------------------------------------------------------------------------------------
    //  Row        | Column     | Reference value    | User-defined value | Absolute error
    // ------------+------------+--------------------+--------------------+--------------------
    //  2          | 0          | 0.540302           | 0                  | -0.540302
    //  2          | 1          | 0                  | 0.540302           | 0.540302
    //  0          | 3          | 0.5                | 1                  | 0.5
    //  3          | 3          | 1                  | 0                  | -1
    // ----------------------------------------------------------------------------------------
    //
    // Here we can see that:
    //     - element (2, 1) should be (2, 0),
    //     - element (0, 3) is incorrect,
    //     - element (3, 3) is not defined at all

    // Check "good" Jacobian at `x0` and `x1` and times t = 1 and 2
    {
        auto jac_comparison = JacobianCompare(MyJacobianGood(), MyRHS());

        auto N_errors_1 = jac_comparison(x0, 1.0);
        auto N_errors_2 = jac_comparison(x1, 2.0);

        std::cout << "Tested \"good\" Jacobian at x0, found " << N_errors_1 << " errors.\n";
        std::cout << "Tested \"good\" Jacobian at x1, found " << N_errors_2 << " errors.\n";
    }

    // Note that Jacobian comparisons are element-by-element and hence can be slow.
    // These comparisons should be removed from the production runs.

    return 0;
}
