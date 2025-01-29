/*
 * Testing:
 * Order of convergence and accuracy of the BDF integrators used in the DAE solver.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#include <dae-cpp/solver.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

/*
 * Singular mass matrix in sparse format:
 *
 * M = |1 0|
 *     |0 0|
 */
struct TestMassMatrix
{
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
struct TestRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];          // y
        f[1] = cos(t) - x[1]; // cos(t) - y
    }
};

/*
 * Analytic Jacobian in sparse format.
 *
 * J = |0  1|
 *     |0 -1|
 */
struct TestJacobian
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(2);  // Pre-allocates memory for 2 non-zero elements (optional)
        J(0, 1, 1.0);  // Row 0, column 1, non-zero element 1.0
        J(1, 1, -1.0); // Row 1, column 1, non-zero element -1.0
    }
};

// Absolute errors
constexpr double conv_err{0.05};           // I.e., for example, for 4th order BDF, the estimated order is in the interval [3.95, 4.05],
                                           // which is quite tight!
constexpr double conv_err_nonuniform{0.1}; // Absolute error for non-uniform case.
                                           // E.g., 4th order BDF is expected to be in the [3.8, 4.0] interval.

/*
 * Tests BDF order vs expected (estimated) order of convergence on uniform grid with automatic Jacobian
 */
TEST(SolverConvergence, Uniform)
{
    System sys((TestMassMatrix()), (TestRHS()));

    state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
    double t_end{10.0};    // Solution interval: t = [0, t_end] -- several periods of sin(t)

    sys.opt.atol = 1e-10;
    sys.opt.rtol = 1e-10;
    sys.opt.dt_init = 1e-6;
    sys.opt.variability_threshold_low = 1.0;
    sys.opt.variability_threshold_high = 1.0;

    std::vector<double> norm_dt_1(DAECPP_MAX_ORDER + 1);
    std::vector<double> norm_dt_2(DAECPP_MAX_ORDER + 1);

    for (unsigned order = 1; order <= DAECPP_MAX_ORDER; order++)
    {
        sys.sol.x.clear();
        sys.sol.t.clear();

        sys.opt.BDF_order = order;
        sys.opt.dt_max = 0.1;

        ASSERT_EQ(sys.solve(x0, t_end), 0);
        ASSERT_EQ(sys.status, 0);

        // Computes Norm 1 for dt1
        for (std::size_t i = 0; i < sys.sol.t.size(); ++i)
        {
            double err = std::abs(sys.sol.x[i][0] - sin(sys.sol.t[i]));
            if (err > norm_dt_1[order])
            {
                norm_dt_1[order] = err;
            }
        }

        sys.opt.dt_max /= 10;

        sys.sol.x.clear();
        sys.sol.t.clear();

        ASSERT_EQ(sys.solve(x0, t_end), 0);
        ASSERT_EQ(sys.status, 0);

        // Computes Norm 1 for dt2
        for (std::size_t i = 0; i < sys.sol.t.size(); ++i)
        {
            double err = abs(sys.sol.x[i][0] - sin(sys.sol.t[i]));
            if (err > norm_dt_2[order])
            {
                norm_dt_2[order] = err;
            }
        }

        double order_estimated = log10(norm_dt_1[order] / norm_dt_2[order]);

        // Estimated order of convergence
        EXPECT_NEAR(order_estimated, static_cast<double>(order), conv_err);

        // The higher order, the less error
        if (order > 1)
        {
            EXPECT_LT(norm_dt_1[order], norm_dt_1[order - 1]);
            EXPECT_LT(norm_dt_2[order], norm_dt_2[order - 1]);
        }
    }
}

/*
 * Tests BDF order vs expected (estimated) order of convergence on uniform grid with analytic Jacobian
 */
TEST(SolverConvergence, UniformWithJacobian)
{
    System sys((TestMassMatrix()), (TestRHS()));

    state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
    double t_end{10.0};    // Solution interval: t = [0, t_end] -- several periods of sin(t)

    sys.opt.atol = 1e-10;
    sys.opt.rtol = 1e-10;
    sys.opt.dt_init = 1e-6;
    sys.opt.variability_threshold_low = 1.0;
    sys.opt.variability_threshold_high = 1.0;

    std::vector<double> norm_dt_1(DAECPP_MAX_ORDER + 1);
    std::vector<double> norm_dt_2(DAECPP_MAX_ORDER + 1);

    for (unsigned order = 1; order <= DAECPP_MAX_ORDER; order++)
    {
        sys.sol.x.clear();
        sys.sol.t.clear();

        sys.opt.BDF_order = order;
        sys.opt.dt_max = 0.1;

        ASSERT_EQ(sys.solve(x0, t_end, TestJacobian()), 0);
        ASSERT_EQ(sys.status, 0);

        // Computes Norm 1 for dt1
        for (std::size_t i = 0; i < sys.sol.t.size(); ++i)
        {
            double err = std::abs(sys.sol.x[i][0] - sin(sys.sol.t[i]));
            if (err > norm_dt_1[order])
            {
                norm_dt_1[order] = err;
            }
        }

        sys.opt.dt_max /= 10;

        sys.sol.x.clear();
        sys.sol.t.clear();

        ASSERT_EQ(sys.solve(x0, t_end, TestJacobian()), 0);
        ASSERT_EQ(sys.status, 0);

        // Computes Norm 1 for dt2
        for (std::size_t i = 0; i < sys.sol.t.size(); ++i)
        {
            double err = abs(sys.sol.x[i][0] - sin(sys.sol.t[i]));
            if (err > norm_dt_2[order])
            {
                norm_dt_2[order] = err;
            }
        }

        double order_estimated = log10(norm_dt_1[order] / norm_dt_2[order]);

        // Estimated order of convergence
        EXPECT_NEAR(order_estimated, static_cast<double>(order), conv_err);

        // The higher order, the less error
        if (order > 1)
        {
            EXPECT_LT(norm_dt_1[order], norm_dt_1[order - 1]);
            EXPECT_LT(norm_dt_2[order], norm_dt_2[order - 1]);
        }
    }
}

/*
 * Tests BDF order vs expected (estimated) order of convergence on non-uniform grid
 */
TEST(SolverConvergence, NonUniform)
{
    System sys((TestMassMatrix()), (TestRHS()));

    state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
    double t_end{10.0};    // Solution interval: t = [0, t_end] -- several periods of sin(t)

    // Lists of output times that break the time step uniformity
    std::vector<double> t_list_1;
    std::vector<double> t_list_2;

    // Populate list 1 for bigger time step
    double t{0.123};
    while (t < t_end)
    {
        t_list_1.push_back(t);
        t += 0.123;
    }
    t_list_1.push_back(t_end);

    // Populate list 2 for smaller time step
    t = 0.0134;
    while (t < t_end)
    {
        t_list_2.push_back(t);
        t += 0.0134;
    }
    t_list_2.push_back(t_end);

    sys.opt.atol = 1e-10;
    sys.opt.rtol = 1e-10;
    sys.opt.dt_init = 1e-6;
    sys.opt.variability_threshold_low = 1.0;
    sys.opt.variability_threshold_high = 1.0;

    std::vector<double> norm_dt_1(DAECPP_MAX_ORDER + 1);
    std::vector<double> norm_dt_2(DAECPP_MAX_ORDER + 1);

    for (unsigned order = 1; order <= DAECPP_MAX_ORDER; order++)
    {
        sys.sol.x.clear();
        sys.sol.t.clear();

        sys.opt.BDF_order = order;
        sys.opt.dt_max = 0.1;

        ASSERT_EQ(sys.solve(x0, t_list_1), 0);
        ASSERT_EQ(sys.status, 0);

        // Computes Norm 1 for dt1
        for (std::size_t i = 0; i < sys.sol.t.size(); ++i)
        {
            double err = std::abs(sys.sol.x[i][0] - sin(sys.sol.t[i]));
            if (err > norm_dt_1[order])
            {
                norm_dt_1[order] = err;
            }
        }

        sys.opt.dt_max /= 10;

        sys.sol.x.clear();
        sys.sol.t.clear();

        ASSERT_EQ(sys.solve(x0, t_list_2), 0);
        ASSERT_EQ(sys.status, 0);

        // Computes Norm 1 for dt2
        for (std::size_t i = 0; i < sys.sol.t.size(); ++i)
        {
            double err = abs(sys.sol.x[i][0] - sin(sys.sol.t[i]));
            if (err > norm_dt_2[order])
            {
                norm_dt_2[order] = err;
            }
        }

        double order_estimated = log10(norm_dt_1[order] / norm_dt_2[order]);

        // Estimated order of convergence
        EXPECT_NEAR(order_estimated, static_cast<double>(order) - 0.1, conv_err_nonuniform);

        // The higher order, the less error
        if (order > 1)
        {
            EXPECT_LT(norm_dt_1[order], norm_dt_1[order - 1]);
            EXPECT_LT(norm_dt_2[order], norm_dt_2[order - 1]);
        }
    }
}

} // namespace
