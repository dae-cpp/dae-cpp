/*
 * Testing:
 * struct SolverOptions
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#include <dae-cpp/solver-options.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

/*
 * To make sure default options pass the checks, can be adjusted, and are accessible
 */
TEST(SolverOptions, DefaultCheck)
{
    SolverOptions opt;

    opt.check();

    opt.atol = 0.0123;

    EXPECT_DOUBLE_EQ(opt.atol, 0.0123);
}

} // namespace
