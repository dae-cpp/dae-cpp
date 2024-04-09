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
        M(0, 0, 2.0 * _param);
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
    state_vector &_x;
    state_vector &_y;
    std::vector<double> &_t;

public:
    TestManager(state_vector &x, state_vector &y, std::vector<double> &t) : _x(x), _y(y), _t(t) {}

    int operator()(const state_vector &x, const double t)
    {
        _x.push_back(x[0]);
        _y.push_back(x[1]);
        _t.push_back(t);
        return 0;
    }
};

TEST(SolverDefinition, SolverFunctions)
{
    double param{0.5};
    double time{1.0};

    TestMassMatrix mass(param);
    TestRHS rhs(param);

    // Initial condition
    state_vector x0(2);
    x0[0] = 1.0;
    x0[1] = -param * x0[0];

    // Solution
    state_vector x, y;
    std::vector<double> t;

    TestManager mgr(x, y, t);

    EXPECT_EQ(solve(mass, rhs, x0, time, mgr), 0);
}

} // namespace
