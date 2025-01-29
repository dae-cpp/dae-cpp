/*
 * Testing:
 * class JacobianMatrix, JacobianAutomatic
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#include <dae-cpp/jacobian-matrix.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

TEST(JacobianMatrix, Definition)
{
    struct TestJacobian : JacobianMatrix
    {
        void operator()(sparse_matrix &J, const state_vector &x, const double t) const
        {
            EXPECT_EQ(J.N_elements(), 0);

            J.reserve(3);
            J(0, 1, 1.0);
            J(1, 0, 2.0 * t);
            J(1, 1, 3.0 * x[1]);
        }
    };

    TestJacobian jac;
    sparse_matrix J;

    state_vector x{4.0, 6.0, 8.0};

    constexpr double t{10.0};

    jac(J, x, t);

    J.check();

    EXPECT_DOUBLE_EQ(J.A[0], 1.0);
    EXPECT_DOUBLE_EQ(J.A[1], 2.0 * t);
    EXPECT_DOUBLE_EQ(J.A[2], 3.0 * x[1]);

    EXPECT_EQ(J.i[0], 0);
    EXPECT_EQ(J.i[1], 1);
    EXPECT_EQ(J.i[2], 1);

    EXPECT_EQ(J.j[0], 1);
    EXPECT_EQ(J.j[1], 0);
    EXPECT_EQ(J.j[2], 1);

    EXPECT_EQ(J.N_elements(), 3);
}

TEST(JacobianMatrix, Automatic)
{
    struct TestVectorFunction
    {
        void operator()(state_type &f, const state_type &x, const double t)
        {
            f[0] = x[0] * x[0] + 25.0;
            f[1] = x[1] * t + 4.0 * x[0] * x[1] + 42.0 * t;
        }
    };

    TestVectorFunction rhs;
    JacobianAutomatic jac(rhs);
    sparse_matrix J;
    state_vector x(2);

    x[0] = 4.0;
    x[1] = 7.0;

    constexpr double t{10.0};

    jac(J, x, t);

    J.check();

    ASSERT_EQ(J.N_elements(), 3); // One of the elements is 0

    EXPECT_DOUBLE_EQ(J.dense(2).coeffRef(0, 0), 2.0 * x[0]);
    EXPECT_DOUBLE_EQ(J.dense(2).coeffRef(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(J.dense(2).coeffRef(1, 0), 4.0 * x[1]);
    EXPECT_DOUBLE_EQ(J.dense(2).coeffRef(1, 1), t + 4.0 * x[0]);
}

TEST(JacobianMatrix, AutomaticSparse)
{
    struct TestVectorFunction
    {
        void operator()(state_type &f, const state_type &x, const double t)
        {
            f[0] = 2.4 * sin(3.0 * x[0]) * t + 15.0;
            f[1] = 2.8 * exp(x[1] * x[1]) * t + 42.0 * t;
            f[2] = 3.2 * exp(x[0] * x[2] * t);
        }
    };

    TestVectorFunction rhs;
    JacobianAutomatic jac(rhs);
    sparse_matrix J;
    state_vector x(3);

    x[0] = 4.0;
    x[1] = 1.5;
    x[2] = 0.8;

    constexpr double t{10.0};

    jac(J, x, t);

    J.check();

    ASSERT_EQ(J.N_elements(), 4);

    // Sparse matrix fills up column by column
    EXPECT_DOUBLE_EQ(J.A[0], 2.4 * t * 3.0 * cos(3.0 * x[0]));
    EXPECT_DOUBLE_EQ(J.A[1], 3.2 * t * x[2] * exp(x[0] * x[2] * t));
    EXPECT_DOUBLE_EQ(J.A[2], 2.8 * t * 2.0 * x[1] * exp(x[1] * x[1]));
    EXPECT_DOUBLE_EQ(J.A[3], 3.2 * t * x[0] * exp(x[0] * x[2] * t));

    EXPECT_EQ(J.i[0], 0);
    EXPECT_EQ(J.i[1], 2);
    EXPECT_EQ(J.i[2], 1);
    EXPECT_EQ(J.i[3], 2);

    EXPECT_EQ(J.j[0], 0);
    EXPECT_EQ(J.j[1], 0);
    EXPECT_EQ(J.j[2], 1);
    EXPECT_EQ(J.j[3], 2);
}

} // namespace
