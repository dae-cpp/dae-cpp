/*
 * Testing:
 * class JacobianCompare
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#include <dae-cpp/jacobian-matrix.hpp>
#include <dae-cpp/vector-function.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

struct MyRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        // `dual_type` for automatic differentiation
        dual_type u = x[0];
        dual_type v = x[1];
        dual_type w = x[2];
        dual_type z = x[3];

        f[0] = z * 0.5;
        f[1] = u * u + v * v + w * w;
        f[2] = sin(u);
        f[3] = exp(-v) + t * z;
    }
};

struct MyJacobianBad
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        float_type u = x[0];
        float_type v = x[1];
        float_type w = x[2];
        float_type z = x[3];

        J.reserve(6);

        J(0, 3, 1.0);      // Row 0, column 3 - error, should be `0.5`
        J(1, 0, 2.0 * u);  // Row 1, column 0
        J(1, 1, 2.0 * v);  // Row 1, column 1
        J(1, 2, 2.0 * w);  // Row 1, column 2
        J(2, 1, cos(u));   // Row 2, column 1 - error, should be column 0
        J(3, 1, -exp(-v)); // Row 3, column 1
        // Forgot to differentiate f[3] w.r.t. `z` - missing element (3, 3)
    }
};

struct MyJacobianGood
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        float_type u = x[0];
        float_type v = x[1];
        float_type w = x[2];
        float_type z = x[3];

        J.reserve(7);

        J(0, 3, 0.5);      // Row 0, column 3
        J(1, 0, 2.0 * u);  // Row 1, column 0
        J(1, 1, 2.0 * v);  // Row 1, column 1
        J(1, 2, 2.0 * w);  // Row 1, column 2
        J(2, 0, cos(u));   // Row 2, column 1
        J(3, 1, -exp(-v)); // Row 3, column 1
        J(3, 3, t);        // Row 3, column 3
    }
};

class MyRHSShape : public VectorFunctionElements
{
public:
    dual_type equations(const state_type &x, const double t, const int_type i) const
    {
        // `dual_type` for automatic differentiation
        dual_type u = x[0];
        dual_type v = x[1];
        dual_type w = x[2];
        dual_type z = x[3];

        if (i == 0)
            return z * 0.5;
        else if (i == 1)
            return u * u + v * v + w * w;
        else if (i == 2)
            return sin(u);
        else if (i == 3)
            return exp(-v) + t * z;
        else
        {
            ERROR("Index i in MyRHSShape is out of boundaries, i = " << i);
        }
    }
};

class MyJacobianShapeBad : public JacobianMatrixShape<MyRHSShape>
{
public:
    MyJacobianShapeBad(MyRHSShape &rhs) : JacobianMatrixShape(rhs)
    {
        add_element(0, 3);
        add_element(1, {0, 1, 2, 3}); // (1, 3) is not necessary (not an error)
        add_element(2, 1);            // Should be (2, 0)
        add_element(3, {1, 2});       // Should be {1, 3}
    }
};

class MyJacobianShapeGood : public JacobianMatrixShape<MyRHSShape>
{
public:
    MyJacobianShapeGood(MyRHSShape &rhs) : JacobianMatrixShape(rhs)
    {
        add_element(0, 3);
        add_element(1, {0, 1, 2, 3}); // (1, 3) is not necessary (not an error)
        add_element(2, 0);
        add_element(3, {1, 3});
    }
};

// Fill vectors x at which the Jacobian matrix will be tested
state_vector x0 = {1.0, 2.0, 3.0, 4.0};
state_vector x1 = {-0.5, -0.1, 0.5, 0.1};

TEST(JacobianCompare, AnalyticJacobianBad)
{
    auto jac_comparison = JacobianCompare(MyJacobianBad(), MyRHS());

    ASSERT_EQ(jac_comparison(x0, 1.0), 4);
    ASSERT_EQ(jac_comparison(x1, 2.0), 4);
}

TEST(JacobianCompare, AnalyticJacobianGood)
{
    auto jac_comparison = JacobianCompare(MyJacobianGood(), MyRHS());

    ASSERT_EQ(jac_comparison(x0, 1.0), 0);
    ASSERT_EQ(jac_comparison(x1, 2.0), 0);
}

TEST(JacobianCompare, JacobianShapeBad)
{
    auto rhs = MyRHSShape();

    auto jac_comparison = JacobianCompare(MyJacobianShapeBad(rhs), rhs);

    ASSERT_EQ(jac_comparison(x0, 1.0), 2);
    ASSERT_EQ(jac_comparison(x1, 2.0), 2);
}

TEST(JacobianCompare, JacobianShapeGood)
{
    auto rhs = MyRHSShape();

    auto jac_comparison = JacobianCompare(MyJacobianShapeGood(rhs), rhs);

    ASSERT_EQ(jac_comparison(x0, 1.0), 0);
    ASSERT_EQ(jac_comparison(x1, 2.0), 0);
}

} // namespace
