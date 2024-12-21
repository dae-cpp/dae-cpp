/*
 * Integration test of the solver.
 * Based on the `perovskite_model` example.
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

struct MyParams
{
    int_type N{4000};   // Number of discretization points
    double L{1.0};      // Length of the domain
    double lambda{1.0}; // Lambda parameter
};

class MyMassMatrix
{
    MyParams p; // Parameters

public:
    explicit MyMassMatrix(MyParams &params) : p(params) {}

    void operator()(sparse_matrix &M, const double t)
    {
        M.reserve(p.N); // Reserves memory for `N` non-zero elements

        for (int_type i = 0; i < p.N; ++i)
        {
            M(i, i, 1.0);
        }
    }
};

class MyRHS
{
    MyParams p; // Parameters

public:
    explicit MyRHS(MyParams &params) : p(params) {}

    void operator()(state_type &f, const state_type &x, const double t)
    {
        // Derived parameters
        const double h = p.L / (double)(p.N - 1);           // Cell size
        const double invh2 = 1.0 / (h * h);                 // Inverse cell size squared
        const double invlam2 = 1.0 / (p.lambda * p.lambda); // Inverse lambda squared

        // An alias for p.N
        const int_type N = p.N;

        // RHS for the ion concentration P = dFlux/dx
        for (int_type i = 1; i < N - 1; ++i)
        {
            f[i] = (x[i + 1] - 2.0 * x[i] + x[i - 1] + 0.5 * ((x[i + 1] + x[i]) * (x[N + i + 1] - x[N + i]) - (x[i] + x[i - 1]) * (x[N + i] - x[N + i - 1]))) * invh2;
        }
        f[0] = (x[1] - x[0] + 0.5 * (x[1] + x[0]) * (x[N + 1] - x[N])) * invh2;                                  // Left BC
        f[N - 1] = -(x[N - 1] - x[N - 2] + 0.5 * (x[N - 1] + x[N - 2]) * (x[2 * N - 1] - x[2 * N - 2])) * invh2; // Right BC

        // RHS for the potential Phi
        for (int_type i = 1; i < N - 1; i++)
        {
            f[i + N] = (x[i + 1 + N] - 2.0 * x[i + N] + x[i - 1 + N]) * invh2 - (1.0 - x[i]) * invlam2;
        }
        f[N] = x[N] + t;                 // Left BC
        f[2 * N - 1] = x[2 * N - 1] - t; // Right BC
    }
};

class MyJacobian
{
    MyParams p; // Parameters

public:
    explicit MyJacobian(MyParams &params) : p(params) {}

    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        // Derived parameters
        const double h = p.L / (double)(p.N - 1);           // Cell size
        const double invh2 = 1.0 / (h * h);                 // Inverse cell size squared
        const double invlam2 = 1.0 / (p.lambda * p.lambda); // Inverse lambda squared

        // An alias for p.N
        const int_type N = p.N;

        J.reserve(12 * N); // Overestimating but it's better than underestimate

        // Fills Jacobian row by row
        for (int_type i = 0; i < 2 * N; ++i)
        {
            if (i == 0)
            {
                J(i, 0, (-1.0 + 0.5 * (x[N + 1] - x[N])) * invh2);
                J(i, 1, (1.0 + 0.5 * (x[N + 1] - x[N])) * invh2);
                J(i, N, -0.5 * (x[0] + x[1]) * invh2);
                J(i, N + 1, 0.5 * (x[0] + x[1]) * invh2);
            }
            else if (i < N - 1)
            {
                J(i, i - 1, (1.0 - 0.5 * (x[N + i] - x[N + i - 1])) * invh2);
                J(i, i, (-2.0 + 0.5 * (x[N + i + 1] - 2.0 * x[N + i] + x[N + i - 1])) * invh2);
                J(i, i + 1, (1.0 + 0.5 * (x[N + i + 1] - x[N + i])) * invh2);
                J(i, N + i - 1, 0.5 * (x[i] + x[i - 1]) * invh2);
                J(i, N + i, -0.5 * (x[i + 1] + 2.0 * x[i] + x[i - 1]) * invh2);
                J(i, N + i + 1, 0.5 * (x[i + 1] + x[i]) * invh2);
            }
            else if (i == N - 1)
            {
                J(i, i - 1, (1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2);
                J(i, i, (-1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2);
                J(i, N + i - 1, 0.5 * (x[N - 1] + x[N - 2]) * invh2);
                J(i, N + i, -0.5 * (x[N - 1] + x[N - 2]) * invh2);
            }
            else if (i == N)
            {
                J(i, N, 1.0);
            }
            else if (i < 2 * N - 1)
            {
                J(i, i - N, invlam2);
                J(i, i - 1, invh2);
                J(i, i, -2.0 * invh2);
                J(i, i + 1, invh2);
            }
            else if (i == 2 * N - 1)
            {
                J(i, 2 * N - 1, 1.0);
            }
            else
            {
                ERROR("Perovskite model: index i is out of boundaries.");
            }
        }
    }
};

class MySolutionManager
{
    // A reference to the solution holder object
    SolutionHolder &m_sol;

public:
    explicit MySolutionManager(SolutionHolder &sol) : m_sol(sol) {}

    virtual int operator()(const state_vector &x, const double t)
    {
        m_sol.x.emplace_back(x);
        m_sol.t.emplace_back(t);

        if (x[0] < 0)
        {
            return -1; // if x[0] is less than 0 (should never happen), then the solver will stop
        }
        else
        {
            return 0;
        }
    }
};

// Absolute errors
constexpr double abs_err{0.0};
constexpr double abs_err_fine{1e-14};

TEST(Integration, PerovskiteModel)
{
    // Solution parameters
    MyParams params;
    params.N = 4000;

    // Initial condition
    state_vector x0(2 * params.N);
    for (int_type i = 0; i < params.N; ++i)
    {
        x0[i] = 1.0;
    }

    // Solution interval: t = [0, t_end]
    double t_end{10.0};

    // Solution holder
    SolutionHolder sol;

    // Solver options
    SolverOptions opt;
    opt.verbosity = verbosity::off;           // Prints computation time and basic info
    opt.solution_variability_control = false; // Switches off solution variability control for better performance

    ASSERT_EQ(JacobianCompare(MyJacobian(params), MyRHS(params))(x0, 0.0), 0);

    // Solve the DAE system
    int status = solve(MyMassMatrix(params), MyRHS(params), MyJacobian(params),
                       x0, t_end,
                       MySolutionManager(sol), opt);

    ASSERT_EQ(status, 0);

    // Soluton vs time `t` is in the `sol` object.
    ASSERT_GT(sol.x.size(), 0);
    ASSERT_GT(sol.t.size(), 0);

    ASSERT_DOUBLE_EQ(sol.t.back(), t_end);

    // Comparison with MATLAB numerical solution
    auto N = params.N;
    EXPECT_NEAR(sol.x.back()[0], 19.9949, 0.05);
    EXPECT_NEAR(sol.x.back()[(N - 1) / 10] * 0.1 + sol.x.back()[(N - 1) / 10 + 1] * 0.9, 2.72523, 1e-2);
    EXPECT_NEAR(sol.x.back()[(N - 1) / 5] * 0.2 + sol.x.back()[(N - 1) / 5 + 1] * 0.8, 0.382148, 1e-3);
    EXPECT_NEAR(sol.x.back()[N], -10.0, abs_err_fine);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 1] * 0.2 + sol.x.back()[N + (N - 1) / 5 * 1 + 1] * 0.8, -6.04056, 1e-4);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 2] * 0.4 + sol.x.back()[N + (N - 1) / 5 * 2 + 1] * 0.6, -2.08970, 0.006);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 3] * 0.6 + sol.x.back()[N + (N - 1) / 5 * 3 + 1] * 0.4, 1.90021, 0.015);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 4] * 0.8 + sol.x.back()[N + (N - 1) / 5 * 4 + 1] * 0.2, 5.93011, 0.02);
    EXPECT_NEAR(sol.x.back().back(), 10.0, abs_err_fine);
}

TEST(Integration, PerovskiteModelDefault)
{
    // Solution parameters
    MyParams params;
    params.N = 4000;

    // Initial condition
    state_vector x0(2 * params.N);
    for (int_type i = 0; i < params.N; ++i)
    {
        x0[i] = 1.0;
    }

    // Solution interval: t = [0, t_end]
    double t_end{10.0};

    // Solution holder
    SolutionHolder sol;

    // Solve the DAE system
    int status = solve(MyMassMatrix(params), MyRHS(params), MyJacobian(params),
                       x0, t_end,
                       MySolutionManager(sol));

    ASSERT_EQ(status, 0);

    // Soluton vs time `t` is in the `sol` object.
    ASSERT_GT(sol.x.size(), 0);
    ASSERT_GT(sol.t.size(), 0);

    ASSERT_DOUBLE_EQ(sol.t.back(), t_end);

    // Comparison with MATLAB numerical solution
    auto N = params.N;
    EXPECT_NEAR(sol.x.back()[0], 19.9949, 0.05);
    EXPECT_NEAR(sol.x.back()[(N - 1) / 10] * 0.1 + sol.x.back()[(N - 1) / 10 + 1] * 0.9, 2.72523, 1e-2);
    EXPECT_NEAR(sol.x.back()[(N - 1) / 5] * 0.2 + sol.x.back()[(N - 1) / 5 + 1] * 0.8, 0.382148, 1e-3);
    EXPECT_NEAR(sol.x.back()[N], -10.0, abs_err_fine);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 1] * 0.2 + sol.x.back()[N + (N - 1) / 5 * 1 + 1] * 0.8, -6.04056, 1e-4);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 2] * 0.4 + sol.x.back()[N + (N - 1) / 5 * 2 + 1] * 0.6, -2.08970, 0.006);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 3] * 0.6 + sol.x.back()[N + (N - 1) / 5 * 3 + 1] * 0.4, 1.90021, 0.015);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 4] * 0.8 + sol.x.back()[N + (N - 1) / 5 * 4 + 1] * 0.2, 5.93011, 0.02);
    EXPECT_NEAR(sol.x.back().back(), 10.0, abs_err_fine);
}

} // namespace
