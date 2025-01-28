/*
 * Testing:
 * class SolutionManager, SolutionHolder, Solution
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

TEST(SolutionManager, Default)
{
    SolutionManager mgr;

    state_vector x(16);

    EXPECT_EQ(mgr(x, 10.0), 0);
}

TEST(SolutionManager, Definition)
{
    struct TestSolutionManager : SolutionManager
    {
        int operator()(const state_vector &x, const double t)
        {
            return 42;
        }
    };

    TestSolutionManager mgr;

    state_vector x(16);

    EXPECT_EQ(mgr(x, 10.0), 42);
}

TEST(SolutionManager, SolutionHolder)
{
    SolutionHolder sol;

    ASSERT_EQ(sol.x.size(), 0);
    ASSERT_EQ(sol.t.size(), 0);

    sol.x.push_back({1.0, 2.0, 3.0});
    sol.x.push_back({4.0, 5.0, 6.0});
    sol.t.push_back(10.0);
    sol.t.push_back(11.0);

    ASSERT_EQ(sol.x.size(), 2);
    ASSERT_EQ(sol.t.size(), 2);

    EXPECT_DOUBLE_EQ(sol.x[0][0], 1.0);
    EXPECT_DOUBLE_EQ(sol.x[0][1], 2.0);
    EXPECT_DOUBLE_EQ(sol.x[0][2], 3.0);
    EXPECT_DOUBLE_EQ(sol.x[1][0], 4.0);
    EXPECT_DOUBLE_EQ(sol.x[1][1], 5.0);
    EXPECT_DOUBLE_EQ(sol.x[1][2], 6.0);
    EXPECT_DOUBLE_EQ(sol.t[0], 10.0);
    EXPECT_DOUBLE_EQ(sol.t[1], 11.0);
}

TEST(SolutionManager, SolutionClass)
{
    SolutionHolder sol;
    Solution sol_obj(sol);

    state_vector x = {1.5, 2.5};

    sol_obj(x, 10.0);

    x = {11.5, 12.5};

    sol_obj(x, 11.0);

    ASSERT_EQ(sol.x.size(), 2);
    ASSERT_EQ(sol.t.size(), 2);

    EXPECT_DOUBLE_EQ(sol.x[0][0], 1.5);
    EXPECT_DOUBLE_EQ(sol.x[0][1], 2.5);
    EXPECT_DOUBLE_EQ(sol.x[1][0], 11.5);
    EXPECT_DOUBLE_EQ(sol.x[1][1], 12.5);
    EXPECT_DOUBLE_EQ(sol.t[0], 10.0);
    EXPECT_DOUBLE_EQ(sol.t[1], 11.0);
}

TEST(SolutionManager, SolutionClassList)
{
    SolutionHolder sol;
    Solution sol_obj(sol, {1.0, 2.0, 12.0, 1.0, 0.0, -5.0, 11.0, 11.0, 1e9}); // Duplicates and out of range values

    state_vector x = {1.5, 2.5};
    sol_obj(x, 10.0);

    x = {11.5, 12.5};
    sol_obj(x, 11.0);

    x = {25.0, 35.0};
    sol_obj(x, 12.0);

    // Only two values should be written
    ASSERT_EQ(sol.x.size(), 2);
    ASSERT_EQ(sol.t.size(), 2);

    EXPECT_DOUBLE_EQ(sol.x[0][0], 11.5);
    EXPECT_DOUBLE_EQ(sol.x[0][1], 12.5);
    EXPECT_DOUBLE_EQ(sol.x[1][0], 25.0);
    EXPECT_DOUBLE_EQ(sol.x[1][1], 35.0);
    EXPECT_DOUBLE_EQ(sol.t[0], 11.0);
    EXPECT_DOUBLE_EQ(sol.t[1], 12.0);
}

struct MyRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = -1; // dx/dt = -1
    }
};

// Absolute error
constexpr double abs_err{1e-6};

class MySolutionManager
{
    SolutionHolder &m_sol;

    bool m_keep_reducing_time_step{false};

    void m_save_solution(const state_vector &x, const double t)
    {
        m_sol.x.emplace_back(x);
        m_sol.t.emplace_back(t);
    }

public:
    MySolutionManager(SolutionHolder &sol) : m_sol(sol) {}

    /*
     * Solution Manager functor will be called every time step providing the time `t` and
     * the corresponding solution `x` for further post-processing.
     */
    int operator()(const state_vector &x, const double t)
    {
        if (std::abs(x[0] - 1.0) < abs_err)
        {
            m_save_solution(x, t);
            return solver_command::stop_intergration;
        }

        if (x[0] < 1.0)
        {
            m_keep_reducing_time_step = true;
            return solver_command::decrease_time_step_and_redo;
        }

        m_save_solution(x, t);

        if (m_keep_reducing_time_step)
        {
            return solver_command::decrease_time_step;
        }

        return 0;
    }
};

TEST(SolutionManager, SolverCommands)
{
    MyRHS rhs; // The vector-function object

    state_vector x0{2.0}; // Initial condition: x = 2
    double t_end{100.0};  // Solution interval: t = [0, t_end] - should stop earlier

    SolutionHolder sol;

    auto status = solve(MassMatrixIdentity(x0.size()), rhs, x0, t_end, MySolutionManager(sol));

    ASSERT_EQ(status, 0);

    ASSERT_GT(sol.x.size(), 0);
    ASSERT_GT(sol.t.size(), 0);
    EXPECT_GT(sol.t.back(), 0.0);

    EXPECT_NEAR(sol.x.back()[0], 1.0, abs_err); // Should stop at x = 1.0
}

} // namespace
