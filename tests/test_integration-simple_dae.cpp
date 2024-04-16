/*
 * Integration test of the solver.
 * Based on the `simple_dae` example.
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

constexpr double pi{3.14159265358979323846};

struct MyMassMatrix
{
    void operator()(sparse_matrix &M, const double t)
    {
        M(0, 0, 1.0); // Row 0, column 0, non-zero element 1.0
    }
};

struct MyRHS
{
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];                          // y
        f[1] = x[0] * x[0] + x[1] * x[1] - 1; // x*x + y*y - 1
    }
};

struct MyJacobian
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(3);        // Pre-allocates memory for 3 non-zero elements (optional)
        J(0, 1, 1.0);        // Row 0, column 1, non-zero element `1.0`
        J(1, 0, 2.0 * x[0]); // Row 1, column 0, non-zero element `2x`
        J(1, 1, 2.0 * x[1]); // Row 1, column 1, non-zero element `2y`
    }
};

class MyObserver
{
    state_vector &_err1, &_err2; // Vectors of the solution errors

public:
    MyObserver(state_vector &err1, state_vector &err2) : _err1(err1), _err2(err2) {}

    virtual int operator()(const state_vector &x, const double t)
    {
        double e1 = std::abs(x[0] * x[0] + x[1] * x[1] - 1.0); // Error 1
        double e2 = std::abs(std::sin(t) - x[0]);              // Error 2

        _err1.push_back(e1);
        _err2.push_back(e2);

        return 0;
    }
};

// Absolute errors
constexpr double abs_err_0{1e-14};
constexpr double abs_err_1{2e-8};

TEST(Integration, SimpleDAE)
{
    state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
    double t_end{pi / 2};  // Solution interval: t = [0, t_end]

    state_vector error1, error2; // Absolute errors

    int status = solve(MyMassMatrix(), MyRHS(), x0, t_end, MyObserver(error1, error2)); // Solves the DAE system

    ASSERT_EQ(status, 0);

    ASSERT_GT(error1.size(), 0);
    ASSERT_GT(error2.size(), 0);
    EXPECT_LT(error1.back(), abs_err_0);
    EXPECT_LT(error2.back(), abs_err_1);

    // With Jacobian
    status = solve(MyMassMatrix(), MyRHS(), MyJacobian(), x0, t_end, MyObserver(error1, error2));

    ASSERT_EQ(status, 0);

    EXPECT_LT(error1.back(), abs_err_0);
    EXPECT_LT(error2.back(), abs_err_1);

    // With user-defined solver options
    SolverOptions opt;
    opt.verbosity = verbosity::off;
    status = solve(MyMassMatrix(), MyRHS(), MyJacobian(), x0, t_end, MyObserver(error1, error2), opt);

    ASSERT_EQ(status, 0);

    EXPECT_LT(error1.back(), abs_err_0);
    EXPECT_LT(error2.back(), abs_err_1);
}

} // namespace
