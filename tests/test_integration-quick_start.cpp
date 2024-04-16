/*
 * Integration test of the solver.
 * Based on the `quick_start` example.
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

struct MyMassMatrix
{
    void operator()(sparse_matrix &M, const double t)
    {
        M(0, 0, 1.0); // Row 0, column 0, non-zero element 1.0
    }
};

struct MyRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];          // y
        f[1] = cos(t) - x[1]; // cos(t) - y
    }
};

struct MyJacobian
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(2);  // Pre-allocates memory for 2 non-zero elements (optional)
        J(0, 1, 1.0);  // Row 0, column 1, non-zero element 1.0
        J(1, 1, -1.0); // Row 1, column 1, non-zero element -1.0
    }
};

// Absolute errors
constexpr double abs_err_0{1e-5};
constexpr double abs_err_1{1e-14};

TEST(Integration, QuickStart)
{
    MyMassMatrix mass; // Mass matrix object
    MyRHS rhs;         // The vector-function object

    state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
    double t_end{1.0};     // Solution interval: t = [0, t_end]

    System my_system(mass, rhs); // System 1

    ASSERT_EQ(my_system.solve(x0, t_end), 0);
    ASSERT_EQ(my_system.status, 0);

    ASSERT_GT(my_system.sol.x.size(), 0);
    ASSERT_GT(my_system.sol.t.size(), 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - cos(t_end)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);

    System my_system2((MyMassMatrix()), (MyRHS())); // System 2

    ASSERT_EQ(my_system2.solve(x0, t_end), 0);
    ASSERT_EQ(my_system2.status, 0);

    ASSERT_GT(my_system2.sol.x.size(), 0);
    ASSERT_GT(my_system2.sol.t.size(), 0);

    EXPECT_LT(std::abs(my_system2.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system2.sol.x.back()[1] - cos(t_end)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system2.sol.t.back(), t_end);

    my_system.opt.verbosity = verbosity::off;
    my_system2.opt.verbosity = verbosity::off;

    ASSERT_EQ(my_system.solve(x0, t_end, MyJacobian()), 0);
    ASSERT_EQ(my_system2.solve(x0, t_end, MyJacobian()), 0);
    ASSERT_EQ(my_system.status, 0);
    ASSERT_EQ(my_system2.status, 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - cos(t_end)), abs_err_1);
    EXPECT_LT(std::abs(my_system2.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system2.sol.x.back()[1] - cos(t_end)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);
    EXPECT_DOUBLE_EQ(my_system2.sol.t.back(), t_end);

    ASSERT_EQ(my_system.solve({0, 1}, 1.0, MyJacobian()), 0);
    ASSERT_EQ(my_system2.solve({0, 1}, 1.0, MyJacobian()), 0);
    ASSERT_EQ(my_system.status, 0);
    ASSERT_EQ(my_system2.status, 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - cos(t_end)), abs_err_1);
    EXPECT_LT(std::abs(my_system2.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system2.sol.x.back()[1] - cos(t_end)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);
    EXPECT_DOUBLE_EQ(my_system2.sol.t.back(), t_end);
}

} // namespace
