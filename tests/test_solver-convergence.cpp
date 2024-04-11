/*
 * Testing:
 * Order of convergence and accuracy of the BDF integrators used in the DAE solver.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#include <dae-cpp/solver.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

/*
 * Mass matrix
 */
class TestMassMatrix
{
public:
    void operator()(sparse_matrix &M, const double t)
    {
        M(0, 0, 1.0);
    }
};

/*
 * RHS (vector function)
 */
class TestRHS
{
    const double _param;

public:
    TestRHS(const double p) : _param(p) {}

    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1] * _param;
        f[1] = x[0] + x[1];
    }
};

/*
 * Analytic Jacobian
 */
class TestJacobian
{
    const double _param;

public:
    TestJacobian(const double p) : _param(p) {}

    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(3);
        J(0, 1, _param);
        J(1, 0, 1.0);
        J(1, 1, 1.0);
    }
};

// Absolute errors
constexpr double conv_err{0.2};
constexpr double max_abs_err{1e-4};

/*
 * Tests BDF order vs expected (estimated) order of convergence on uniform grid
 */
TEST(SolverConvergence, ConvergenceUniform)
{
    double param{1e-3};

    double t = 10.0 / param;

    // Initial condition
    state_vector x0 = {1.0, -1.0};

    // Exact solution (`x[0]`) at time `t`
    double x_exact = exp(-param * t);

    System sys((TestMassMatrix()), (TestRHS(param)));

    sys.opt.atol = 1e-10;
    sys.opt.rtol = 1e-10;

    for (unsigned order = 1; order <= DAECPP_MAX_ORDER; order++)
    {
        sys.opt.BDF_order = order;
        sys.opt.dt_max = t / 20.0; // Relatively big time step

        ASSERT_EQ(sys.solve(x0, t, TestJacobian(param)), 0);
        ASSERT_EQ(sys.status, 0);

        // Absolute error
        double first = abs(sys.sol.x.back()[0] - x_exact);

        if (order == 1)
        {
            EXPECT_LT(first, max_abs_err * 5); // First order method with big time step is innacurate
        }
        else
        {
            EXPECT_LT(first, max_abs_err);
        }

        sys.opt.dt_max /= 10; // Reduce time step 10 times

        ASSERT_EQ(sys.solve(x0, t, TestJacobian(param)), 0);
        ASSERT_EQ(sys.status, 0);

        // Absolute error after reducing the time step
        double second = abs(sys.sol.x.back()[0] - x_exact);

        EXPECT_LT(second, max_abs_err);

        // Order of convergence estimation
        double order_est = log10(first / second);

        EXPECT_NEAR(order_est, static_cast<double>(order) + 0.1, conv_err);
    }
}

/*
 * Tests BDF order vs expected (estimated) order of convergence on non-uniform grid
 */
TEST(SolverConvergence, NonUniform)
{
    double param{1e-3};

    // List of output times that breaks the time step uniformity
    std::vector<double> t_list =
        {77, 166, 499, 999, 1450, 1990, 2222, 2888, 3333, 3850, 4444, 4930, 5555, 6060, 6500, 7070, 7777, 8333, 9090, 9500, 9090, 10000};

    // Initial condition
    state_vector x0 = {1.0, -1.0};

    // Exact solution (`x[0]`) at time `t`
    double x_exact = exp(-param * 10000.0);

    System sys((TestMassMatrix()), (TestRHS(param)));

    sys.opt.atol = 1e-10;
    sys.opt.rtol = 1e-10;

    for (unsigned order = 1; order <= DAECPP_MAX_ORDER; order++)
    {
        sys.opt.BDF_order = order;
        sys.opt.dt_max = 10000.0 / 20.0; // Relatively big time step

        ASSERT_EQ(sys.solve(x0, t_list, TestJacobian(param)), 0);
        ASSERT_EQ(sys.status, 0);

        // Absolute error
        double first = abs(sys.sol.x.back()[0] - x_exact);

        if (order == 1)
        {
            EXPECT_LT(first, max_abs_err * 5); // First order method with big time step is innacurate
        }
        else
        {
            EXPECT_LT(first, max_abs_err);
        }

        sys.opt.dt_max /= 10; // Reduce time step 10 times

        ASSERT_EQ(sys.solve(x0, t_list, TestJacobian(param)), 0);
        ASSERT_EQ(sys.status, 0);

        // Absolute error after reducing the time step
        double second = abs(sys.sol.x.back()[0] - x_exact);

        EXPECT_LT(second, max_abs_err);

        // Order of convergence estimation
        double order_est = log10(first / second);

        if (order == 1)
        {
            EXPECT_NEAR(order_est, static_cast<double>(order) + 0.1, conv_err);
        }
        else
        {
            EXPECT_GT(order_est, static_cast<double>(order) - 1.0); // Higher-order methods lose up to 1 order of accuracy
        }
    }
}

} // namespace
