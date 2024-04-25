/*
 * Integration test of the solver.
 * Solves system of ODEs (identity mass matrix).
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

struct MyRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = cos(t);                         // dx/dt = cos(t)
        f[1] = x[0] * x[0] + pow(cos(t), 2.0); // dy/dt = x*x + cos^2(t)
    }
};

struct MyJacobian
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J(1, 0, 2.0 * x[0]);
    }
};

// Absolute errors
constexpr double abs_err_0{1e-5};
constexpr double abs_err_1{1e-5};

/*
 * Solves ODE system
 */
TEST(Integration, ODE)
{
    MyRHS rhs; // The vector-function object

    state_vector x0{0.0, 1.0}; // Initial condition: x = 0, y = 1
    double t_end{1.0};         // Solution interval: t = [0, t_end]

    System my_system(MassMatrixIdentity(x0.size()), rhs);

    ASSERT_EQ(my_system.solve(x0, t_end), 0);
    ASSERT_EQ(my_system.status, 0);

    ASSERT_GT(my_system.sol.x.size(), 0);
    ASSERT_GT(my_system.sol.t.size(), 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - (t_end + 1.0)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);

    ASSERT_EQ(my_system.solve(x0, t_end, MyJacobian()), 0);
    ASSERT_EQ(my_system.status, 0);

    EXPECT_LT(std::abs(my_system.sol.x.back()[0] - sin(t_end)), abs_err_0);
    EXPECT_LT(std::abs(my_system.sol.x.back()[1] - (t_end + 1.0)), abs_err_1);
    EXPECT_DOUBLE_EQ(my_system.sol.t.back(), t_end);
}

} // namespace
