/*
 * Testing:
 * class System and solver() definitions
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#include <utility>

#include <dae-cpp/solver.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

/*
 * Mass matrix
 */
class TestMassMatrix
{
    const double _param;

public:
    TestMassMatrix(const double p) : _param(p) {}

    void operator()(sparse_matrix &M, const double t)
    {
        M.reserve(1);
        M(0, 0, 1.0 / _param);
    }
};

/*
 * RHS (vector function)
 */
class TestRHS
{
    const double _param;

public:
    TestRHS(const double p) : _param(p) {}

    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];
        f[1] = x[0] * _param + x[1];
    }
};

/*
 * Analytic Jacobian
 */
class TestJacobian
{
    const double _param;

public:
    TestJacobian(const double p) : _param(p) {}

    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(3);
        J(0, 1, 1.0);
        J(1, 0, _param);
        J(1, 1, 1.0);
    }
};

/*
 * Custom Solution Manager
 */
class TestManager
{
    typedef std::vector<std::pair<float_type, float_type>> pair_type;

    pair_type &_x;
    std::vector<double> &_t;

public:
    TestManager(pair_type &x, std::vector<double> &t) : _x(x), _t(t) {}

    int operator()(const state_vector &x, const double t)
    {
        _x.push_back(std::make_pair(x[0], x[1]));
        _t.push_back(t);
        return 0;
    }
};

/*
 * Analytic solution
 */
auto x_exact(double p, double t)
{
    double e = exp(-p * p * t);
    return std::make_pair(e, -p * e);
}

// Absolute error
constexpr double abs_err{1e-3};

/*
 * Tests all possible parameter combinations to call `solve(...)` function
 */
TEST(SolverDefinition, SolverFunctions)
{
    double param{0.5};
    double time{1.0};

    std::vector<double> time_vec = {0.1, 0.2, 0.5, -0.1, 0.5, 0.8, 1.0, 0.0}; // Unsorted with duplicates and negatives

    TestMassMatrix mass(param);
    TestRHS rhs(param);
    TestJacobian jac(param);
    SolverOptions opt;

    // Initial condition
    state_vector x0(2);
    x0[0] = 1.0;
    x0[1] = -param * x0[0];

    // Numerical solution
    std::vector<std::pair<float_type, float_type>> x;
    std::vector<double> t;

    TestManager mgr(x, t);

    // Minimal
    ASSERT_EQ(solve(mass, rhs, x0, time, mgr), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac
    ASSERT_EQ(solve(mass, rhs, jac, x0, time, mgr), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac + opt
    ASSERT_EQ(solve(mass, rhs, jac, x0, time, mgr, opt), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + opt
    ASSERT_EQ(solve(mass, rhs, x0, time, mgr, opt), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // Minimal with vector of times
    ASSERT_EQ(solve(mass, rhs, x0, time_vec, mgr), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac
    ASSERT_EQ(solve(mass, rhs, jac, x0, time_vec, mgr), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac + opt
    ASSERT_EQ(solve(mass, rhs, jac, x0, time_vec, mgr, opt), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + opt
    ASSERT_EQ(solve(mass, rhs, x0, time_vec, mgr, opt), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // Automatic Jacobian
    JacobianAutomatic jac_auto(rhs);
    ASSERT_EQ(solve(mass, rhs, jac_auto, x0, time, mgr, opt), 0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);
}

/*
 * Tests all possible parameter combinations to call `solve(...)` function,
 * creating temporary objects "in place", in the function call.
 */
TEST(SolverDefinition, SolverFunctionsInPlace)
{
    double param{0.5};
    double time{1.0};

    SolverOptions opt;

    // Numerical solution
    std::vector<std::pair<float_type, float_type>> x;
    std::vector<double> t;

    // Minimal
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    {1.0, -0.5},
                    1.0,
                    TestManager(x, t)),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    TestJacobian(param),
                    {1.0, -0.5},
                    1.0,
                    TestManager(x, t)),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac + opt
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    TestJacobian(param),
                    {1.0, -0.5},
                    1.0,
                    TestManager(x, t),
                    opt),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + opt
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    {1.0, -0.5},
                    1.0,
                    TestManager(x, t),
                    opt),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // Minimal with vector of times
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    {1.0, -0.5},
                    {0.1, 0.2, 0.5, -0.1, 0.5, 0.8, 1.0, 0.0},
                    TestManager(x, t)),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    TestJacobian(param),
                    {1.0, -0.5},
                    {0.1, 0.2, 0.5, -0.1, 0.5, 0.8, 1.0, 0.0},
                    TestManager(x, t)),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + jac + opt
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    TestJacobian(param),
                    {1.0, -0.5},
                    {0.1, 0.2, 0.5, -0.1, 0.5, 0.8, 1.0, 0.0},
                    TestManager(x, t),
                    opt),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // + opt
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    {1.0, -0.5},
                    {0.1, 0.2, 0.5, -0.1, 0.5, 0.8, 1.0, 0.0},
                    TestManager(x, t),
                    opt),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);

    // Explicit automatic Jacobian
    ASSERT_EQ(solve(TestMassMatrix(param),
                    TestRHS(param),
                    JacobianAutomatic(TestRHS(param)),
                    {1.0, -0.5},
                    1.0,
                    TestManager(x, t),
                    opt),
              0);
    EXPECT_NEAR(x.back().first, x_exact(param, time).first, abs_err);
    EXPECT_NEAR(x.back().second, x_exact(param, time).second, abs_err);
}

} // namespace
