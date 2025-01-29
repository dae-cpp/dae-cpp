/*
 * Integration test of the solver.
 * Solves purely algebraic system (zero mass matrix).
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

struct MyRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[0] - sin(t);            // x = sin(t)
        f[1] = 1.0 - x[0] * x[0] - x[1]; // y = 1 - x*x
    }
};

struct MyJacobian
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J(0, 0, 1.0);
        J(1, 0, -2.0 * x[0]);
        J(1, 1, -1.0);
    }
};

// Absolute errors
constexpr double abs_err_0{1e-14};
constexpr double abs_err_1{1e-14};

/*
 * Time-dependent algebraic system
 */
TEST(Integration, Algebraic)
{
    MyRHS rhs; // The vector-function object

    state_vector x0{0.1, 0.9}; // Initial condition: x = 0, y = 1 -- set it slightly off deliberatly
    double t_end{1.0};         // Solution interval: t = [0, t_end]

    System my_system(MassMatrixZero(), rhs);

    ASSERT_EQ(my_system.solve(x0, t_end), 0);
    ASSERT_EQ(my_system.status, 0);

    ASSERT_GT(my_system.sol.x.size(), 0);
    ASSERT_GT(my_system.sol.t.size(), 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - pow(cos(t_end), 2.0)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);

    ASSERT_EQ(my_system.solve(x0, t_end, MyJacobian()), 0);
    ASSERT_EQ(my_system.status, 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - pow(cos(t_end), 2.0)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);
}

struct MyRHS_static
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = 2.0 * x[0] - x[1];  // y = 2x
        f[1] = x[0] * x[0] - x[1]; // y = x*x
    }
};

/*
 * Time-independent (static) algebraic system
 */
TEST(Integration, AlgebraicStatic)
{
    MyRHS_static rhs; // The vector-function object

    state_vector x0{3.0, 3.0}; // Initial guess
    double t_end{1.0};         // Solution interval: t = [0, t_end]

    System my_system(MassMatrixZero(), rhs);

    ASSERT_EQ(my_system.solve(x0, t_end), 0);
    ASSERT_EQ(my_system.status, 0);

    ASSERT_GT(my_system.sol.x.size(), 0);
    ASSERT_GT(my_system.sol.t.size(), 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - 2.0), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - 4.0), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);
}

} // namespace
