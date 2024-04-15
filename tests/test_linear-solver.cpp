/*
 * Testing:
 * Linear solver used in the DAE solver.
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

// Absolute errors
constexpr double abs_err{1e-15};
constexpr double abs_err_big_system{1e-7};

TEST(LinearSolver, SolutionCheck)
{
    // Linear solver
    Eigen::SparseLU<core::eimat> linsolver;

    core::eimat Jb(3, 3); // Linear system matrix
    core::eivec b(3);     // The RHS of the linear system

    Jb.coeffRef(0, 0) = 1.0;
    Jb.coeffRef(0, 1) = 2.0;
    Jb.coeffRef(0, 2) = 3.0;
    Jb.coeffRef(1, 1) = 4.0;
    Jb.coeffRef(1, 2) = 5.0;
    Jb.coeffRef(2, 0) = 6.0;
    Jb.coeffRef(2, 1) = 7.0;

    b[0] = 0.5;
    b[1] = 1.5;
    b[2] = 2.5;

    linsolver.compute(Jb); // Decomposition

    ASSERT_EQ(linsolver.info(), Eigen::Success);

    core::eivec x = linsolver.solve(b); // Solution

    ASSERT_EQ(linsolver.info(), Eigen::Success);

    core::eivec xb = Jb * x;

    b -= xb; // Should be 0

    EXPECT_NEAR(b[0], 0.0, abs_err);
    EXPECT_NEAR(b[1], 0.0, abs_err);
    EXPECT_NEAR(b[2], 0.0, abs_err);
}

struct Params
{
    int_type N{10000};  // Number of discretization points
    double L{1000.0};   // Length of the domain
    double lambda{0.1}; // Lambda parameter
};

class TestJacobian
{
    Params p; // Parameters

public:
    explicit TestJacobian(Params &params) : p(params) {}

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
                J(i, 0, (-1.0 + 0.5 * (x[N + 1] - x[N])) * invh2 - 1.0);
                J(i, 1, (1.0 + 0.5 * (x[N + 1] - x[N])) * invh2);
                J(i, N, -0.5 * (x[0] + x[1]) * invh2);
                J(i, N + 1, 0.5 * (x[0] + x[1]) * invh2);
            }
            else if (i < N - 1)
            {
                J(i, i - 1, (1.0 - 0.5 * (x[N + i] - x[N + i - 1])) * invh2);
                J(i, i, (-2.0 + 0.5 * (x[N + i + 1] - 2.0 * x[N + i] + x[N + i - 1])) * invh2 - 1.0);
                J(i, i + 1, (1.0 + 0.5 * (x[N + i + 1] - x[N + i])) * invh2);
                J(i, N + i - 1, 0.5 * (x[i] + x[i - 1]) * invh2);
                J(i, N + i, -0.5 * (x[i + 1] + 2.0 * x[i] + x[i - 1]) * invh2);
                J(i, N + i + 1, 0.5 * (x[i + 1] + x[i]) * invh2);
            }
            else if (i == N - 1)
            {
                J(i, i - 1, (1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2);
                J(i, i, (-1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2 - 1.0);
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

TEST(LinearSolver, BigSystem)
{
    // Linear solver
    Eigen::SparseLU<core::eimat> linsolver;

    Params p;

    state_vector x0(2 * p.N);
    sparse_matrix J;

    for (int_type i = 0; i < p.N; ++i)
    {
        x0[i] = 1.0;
    }

    TestJacobian jac(p);

    jac(J, x0, 0.0);

    core::eimat Jb = J.convert(2 * p.N); // Linear system matrix
    core::eivec b(2 * p.N);              // The RHS of the linear system

    for (auto &b_ : b)
    {
        b_ = 1.24;
    }

    linsolver.compute(Jb); // Decomposition

    ASSERT_EQ(linsolver.info(), Eigen::Success);

    core::eivec x = linsolver.solve(b); // Solution

    ASSERT_EQ(linsolver.info(), Eigen::Success);

    core::eivec xb = Jb * x;

    b -= xb; // Should be 0

    for (std::size_t i = 0; i < b.size(); ++i)
    {
        ASSERT_NEAR(b[i], 0.0, abs_err_big_system) << "i=" << i;
    }
}

} // namespace
