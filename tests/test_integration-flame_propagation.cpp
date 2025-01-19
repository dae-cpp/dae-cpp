/*
 * Integration test of the solver.
 * Based on the `flame_propagation` example.
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
        state_value y = x[0];
        f[0] = y * y - y * y * y;
    }
};

// Absolute error
constexpr double abs_err_sol{2e-3};

TEST(Integration, FlamePropagation)
{
    MyRHS rhs; // The vector-function object

    System my_system(MassMatrixIdentity(1), rhs); // Defines the DAE system object

    const double delta = 1e-4; // Parameter `delta`: 0.01 - not very stiff problem, 1e-4 - stiff problem

    state_vector x0{delta};    // Initial condition: y(0) = delta
    double t_end{2.0 / delta}; // Solution interval: 0 < t < (2 / delta)

    // Update the solver options
    my_system.opt.dt_init = t_end / 100.0;           // Increase the initial time step
    my_system.opt.dt_increase_threshold_delta = -2;  // Try to increase the time step less often
    my_system.opt.variability_threshold_low = 0.10;  // Do not increase time step if the solution changes more than 10%
    my_system.opt.variability_threshold_high = 0.10; // Decrease the time step if the solution changes more than 10%

    // Solves the DAE system `my_system` with the given initial condition `x0` and time `t_end`
    ASSERT_EQ(my_system.solve(x0, t_end), 0);
    ASSERT_EQ(my_system.status, 0);

    ASSERT_GT(my_system.sol.x.size(), 0);
    ASSERT_GT(my_system.sol.t.size(), 0);

    for (std::size_t i = 0; i < my_system.sol.t.size(); i++)
    {
        if (my_system.sol.t[i] > (1.0 / delta * 1.001))
        {
            EXPECT_NEAR(my_system.sol.x[i][0], 1.0, abs_err_sol) << "t = " << my_system.sol.t[i];
        }
    }
}

} // namespace
