/*
 * Testing:
 * class JacobianMatrixShape
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#include <dae-cpp/jacobian-matrix.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

struct TestRHS : VectorFunctionElements
{
    state_value equations(const state_type &x, const double t, const int_type i) const
    {
        if (i == 0)
        {
            return x[0] * x[1];
        }
        else if (i == 1)
        {
            return x[0] - x[1];
        }
        else
        {
            ERROR("Incorrect index i: " << i);
        }
    }
};

TEST(JacobianMatrixShape, Definition)
{
    struct TestJacobian : JacobianMatrixShape<TestRHS>
    {
        explicit TestJacobian(TestRHS rhs) : JacobianMatrixShape(rhs)
        {
            add_element(0, {0, 1});
            add_element(1, 1); // Missing element (1, 0) on purpose
        }
    };

    TestJacobian jac(TestRHS{});

    sparse_matrix J;
    state_vector x{4.0, 6.0};

    constexpr double t{10.0};

    jac(J, x, t);

    J.check();

    EXPECT_EQ(J.N_elements(), 3);

    auto Jd = J.dense(2);

    EXPECT_DOUBLE_EQ(Jd(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(Jd(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(Jd(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(Jd(1, 1), -1.0);
}

TEST(JacobianMatrixShape, DefinitionInline)
{
    JacobianMatrixShape jac(TestRHS{});

    jac.clear();
    jac.reserve(2);

    jac.add_element(0, {0, 1});
    jac.add_element(1, 1); // Missing element (1, 0) on purpose

    sparse_matrix J;
    state_vector x{4.0, 6.0};

    constexpr double t{10.0};

    jac(J, x, t);

    J.check();

    EXPECT_EQ(J.N_elements(), 3);

    auto Jd = J.dense(2);

    EXPECT_DOUBLE_EQ(Jd(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(Jd(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(Jd(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(Jd(1, 1), -1.0);
}

TEST(JacobianMatrixShape, AddElement)
{
    JacobianMatrixShape jac(TestRHS{});

    jac.add_element(0, {0, 1}); // Add as a vector
    jac.add_element(1, 1);      // Add as a single element

    sparse_matrix J;
    state_vector x{4.0, 6.0};

    constexpr double t{10.0};

    jac(J, x, t);

    J.check();

    EXPECT_EQ(J.N_elements(), 3);

    auto Jd = J.dense(2);

    EXPECT_DOUBLE_EQ(Jd(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(Jd(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(Jd(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(Jd(1, 1), -1.0);
}

TEST(JacobianMatrixShape, Clear)
{
    JacobianMatrixShape jac(TestRHS{});

    sparse_matrix J;
    state_vector x{4.0, 6.0};

    constexpr double t{10.0};

    jac.add_element(0, {0, 1}); // Add as a vector
    jac.add_element(1, 0);      // Add as a single element
    jac.clear();

    jac(J, x, t);

    J.check();

    EXPECT_EQ(J.N_elements(), 0);

    jac.add_element(0, {0, 1}); // Add as a vector
    jac.add_element(1, 1);      // Add as a single element (a different one)

    jac(J, x, t);

    J.check();

    EXPECT_EQ(J.N_elements(), 3);

    auto Jd = J.dense(2);

    EXPECT_DOUBLE_EQ(Jd(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(Jd(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(Jd(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(Jd(1, 1), -1.0);
}

} // namespace
