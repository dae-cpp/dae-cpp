/*
 * Testing:
 * class VectorFunction, VectorFunctionElements
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#include <dae-cpp/vector-function.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

TEST(VectorFunction, Definition)
{
    struct TestVectorFunction : VectorFunction
    {
        void operator()(state_type &f, const state_type &x, const double t) const
        {
            ASSERT_EQ(x.size(), 2);

            f[0] = x[0];
            f[1] = x[1] * t;
        }
    };

    TestVectorFunction rhs;

    state_type x(2);

    x[0] = 4.0;
    x[1] = 6.0;

    state_type f(2);

    constexpr double t{10.0};

    rhs(f, x, t);

    EXPECT_DOUBLE_EQ(f[0].val(), 4.0);
    EXPECT_DOUBLE_EQ(f[1].val(), 6.0 * t);

    EXPECT_EQ(f.size(), 2);
}

TEST(VectorFunctionElements, Definition)
{
    struct TestVectorFunction : VectorFunctionElements
    {
        state_value equations(const state_type &x, const double t, const int_type i) const
        {
            ASSERT(x.size() == 2, "Incorrect size of x: " << x.size());

            if (i == 0)
                return x[0];
            else if (i == 1)
                return x[1] * t;
            else
            {
                ERROR("Index i in TestVectorFunction is out of boundaries, i = " << i);
            }
        }
    };

    TestVectorFunction rhs;

    state_type x(2);

    x[0] = 4.0;
    x[1] = 6.0;

    state_type f(2);

    constexpr double t{10.0};

    rhs(f, x, t);

    EXPECT_DOUBLE_EQ(f[0].val(), 4.0);
    EXPECT_DOUBLE_EQ(f[1].val(), 6.0 * t);

    EXPECT_EQ(f.size(), 2);
}

} // namespace
