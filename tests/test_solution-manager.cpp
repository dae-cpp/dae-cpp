/*
 * Testing:
 * class SolutionManager, SolutionHolder, Solution
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#include <dae-cpp/solution-manager.hpp>

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

} // namespace
