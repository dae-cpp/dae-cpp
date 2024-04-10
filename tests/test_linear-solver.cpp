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

// Absolute error
constexpr double abs_err{1e-15};

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

} // namespace
