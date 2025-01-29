/*
 * Integration test of the solver.
 * Based on the `jacobian_shape` example.
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

struct MyParams
{
    int_type N{4000};   // Number of discretization points
    double L{1.0};      // Length of the domain
    double lambda{1.0}; // Lambda parameter
};

class MyMassMatrix
{
    const MyParams &p; // Parameters

public:
    explicit MyMassMatrix(const MyParams &params) : p(params) {}

    /*
     * Defines the mass matrix `M` of the DAE system `M dx/dt = f`.
     * The mass matrix should be defined in sparse format (non-zero elements only) and can depend on time `t`.
     * Matrix `M` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &M, const double t)
    {
        M.reserve(p.N); // Reserves memory for `N` non-zero elements

        for (int_type i = 0; i < p.N; ++i)
        {
            M(i, i, 1.0);
        }
    }
};

class MyRHS : public VectorFunctionElements
{
    const MyParams &p; // Parameters

    // Derived parameters
    const double h = p.L / (double)(p.N - 1);           // Cell size
    const double invh2 = 1.0 / (h * h);                 // Inverse cell size squared
    const double invlam2 = 1.0 / (p.lambda * p.lambda); // Inverse lambda squared

    /*
     * RHS for the ion concentration P = dFlux/dx:
     * dP/dt = d/dx(dP/dx + P * dPhi/dx)
     */
    state_value eq1(const state_type &x, const double t, const int_type i_global) const
    {
        const state_value *P = x.data();
        const state_value *Phi = x.data() + p.N;

        const int_type i = i_global;

        if (i == 0)
            return (P[i + 1] - P[i] + 0.5 * (P[i + 1] + P[i]) * (Phi[i + 1] - Phi[i])) * invh2; // Left BC
        else if (i == p.N - 1)
            return -(P[i] - P[i - 1] + 0.5 * (P[i] + P[i - 1]) * (Phi[i] - Phi[i - 1])) * invh2; // Right BC
        else
            return (P[i + 1] - 2.0 * P[i] + P[i - 1] + 0.5 * ((P[i + 1] + P[i]) * (Phi[i + 1] - Phi[i]) - (P[i] + P[i - 1]) * (Phi[i] - Phi[i - 1]))) * invh2;
    }

    /*
     * RHS for the potential Phi:
     * d^2(Phi)/dx^2 - (1 - P)/lambda^2 = 0
     */
    state_value eq2(const state_type &x, const double t, const int_type i_global) const
    {
        const state_value *P = x.data();
        const state_value *Phi = x.data() + p.N;

        const int_type i = i_global - p.N;

        if (i == 0)
            return Phi[i] + t; // Left BC
        else if (i == p.N - 1)
            return Phi[i] - t; // Right BC
        else
            return (Phi[i + 1] - 2.0 * Phi[i] + Phi[i - 1]) * invh2 - (1.0 - P[i]) * invlam2;
    }

public:
    explicit MyRHS(const MyParams &params) : p(params) {}

    /*
     * All equations combined.
     * This function returns the i-th component of the vector function.
     */
    state_value equations(const state_type &x, const double t, const int_type i) const
    {
        if (i < p.N)
            return eq1(x, t, i);
        else if (i < 2 * p.N)
            return eq2(x, t, i);
        else
        {
            ERROR("Equation system: index i is out of boundaries.");
        }
    }
};

class MyJacobianShape : public JacobianMatrixShape<MyRHS>
{
    const int_type N{0}; // The number of discretization points for each equation

    /*
     * Defines all non-zero elements of the Jacobian matrix row by row.
     * Here we use `add_element()` helper method from the `JacobianMatrixShape` class.
     */
    void m_define_Jacobian_shape()
    {
        reserve(10 * N); // Reserve memory for (6N + 4N) non-zero elements

        for (int_type i = 0; i < 2 * N; ++i)
        {
            if (i == 0)
            {
                add_element(i, {i + 1, i, N + i + 1, N + i}); // The order does not matter
            }
            else if (i < N - 1)
            {
                add_element(i, {i - 1, i, i + 1, N + i - 1, N + i, N + i + 1});
            }
            else if (i == N - 1)
            {
                add_element(i, {i - 1, i, N + i - 1, N + i});
            }
            else if (i == N)
            {
                add_element(i, N); // Can pass one index instead of a vector
            }
            else if (i < 2 * N - 1)
            {
                add_element(i, {i - N, i - 1, i, i + 1});
            }
            else if (i == 2 * N - 1)
            {
                add_element(i, 2 * N - 1);
            }
            else
            {
                ERROR("Jacobian shape: index i is out of boundaries.");
            }
        }
    }

public:
    MyJacobianShape(MyRHS &rhs, const int_type N) : JacobianMatrixShape(rhs), N(N)
    {
        m_define_Jacobian_shape();
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
constexpr double abs_err{1e-14};

TEST(Integration, JacobianShape)
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

    // The vector function and Jacobian objects
    MyRHS rhs = MyRHS(params);
    MyJacobianShape jac = MyJacobianShape(rhs, params.N);

    ASSERT_EQ(JacobianCompare(jac, rhs)(x0, 0.0), 0);

    // Solve the DAE system using automatic Jacobian computed from the user-defined shape
    int status = solve(MyMassMatrix(params), rhs, jac,
                       x0, t_end,
                       MySolutionManager(sol), opt);

    ASSERT_EQ(status, 0);

    // Soluton vs time `t` is in the `sol` object
    ASSERT_GT(sol.x.size(), 0);
    ASSERT_GT(sol.t.size(), 0);

    ASSERT_DOUBLE_EQ(sol.t.back(), t_end);

    // Comparison with MATLAB numerical solution
    auto N = params.N;
    EXPECT_NEAR(sol.x.back()[0], 19.9949, 0.05);
    EXPECT_NEAR(sol.x.back()[(N - 1) / 10] * 0.1 + sol.x.back()[(N - 1) / 10 + 1] * 0.9, 2.72523, 1e-2);
    EXPECT_NEAR(sol.x.back()[(N - 1) / 5] * 0.2 + sol.x.back()[(N - 1) / 5 + 1] * 0.8, 0.382148, 1e-3);
    EXPECT_NEAR(sol.x.back()[N], -10.0, abs_err);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 1] * 0.2 + sol.x.back()[N + (N - 1) / 5 * 1 + 1] * 0.8, -6.04056, 1e-4);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 2] * 0.4 + sol.x.back()[N + (N - 1) / 5 * 2 + 1] * 0.6, -2.08970, 0.006);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 3] * 0.6 + sol.x.back()[N + (N - 1) / 5 * 3 + 1] * 0.4, 1.90021, 0.015);
    EXPECT_NEAR(sol.x.back()[N + (N - 1) / 5 * 4] * 0.8 + sol.x.back()[N + (N - 1) / 5 * 4 + 1] * 0.2, 5.93011, 0.02);
    EXPECT_NEAR(sol.x.back().back(), 10.0, abs_err);
}

} // namespace
