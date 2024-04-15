/*
 * Testing:
 * The backward differentiation formula (BDF) used for the numerical integration in the DAE solver.
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

using namespace daecpp::core;
using namespace daecpp::core::detail;

// Time derivative absolute error
constexpr double abs_err{1e-7};

/*
 * Tests BDF-I time derivative approximation and its corresponding derivative w.r.t. xk
 */
TEST(TimeIntegrator, SchemeExact)
{
    constexpr std::size_t size{6};

    eivec dxdt(size);
    rvec xk(size);
    SolverState state(size);

    double val{};

    xk = {0.0, -1.5, 10, 4000, 4566.6, 42};
    state.x[0] = {0.01, 1.5, 4000, 10, 4567, 42};

    state.dt[0] = 1.23;
    val = time_derivative_approx(dxdt, xk, state, size);
    for (std::size_t i = 0; i < size; ++i)
    {
        EXPECT_DOUBLE_EQ(dxdt[i], -(state.x[0][i] - xk[i]) / state.dt[0]);
    }
    EXPECT_DOUBLE_EQ(val, 1.0 / state.dt[0]);

    state.dt[0] = 1.4e-6;
    val = time_derivative_approx(dxdt, xk, state, size);
    for (std::size_t i = 0; i < size; ++i)
    {
        EXPECT_DOUBLE_EQ(dxdt[i], -(state.x[0][i] - xk[i]) / state.dt[0]);
    }
    EXPECT_DOUBLE_EQ(val, 1.0 / state.dt[0]);

    state.dt[0] = 2.4e6;
    val = time_derivative_approx(dxdt, xk, state, size);
    for (std::size_t i = 0; i < size; ++i)
    {
        EXPECT_DOUBLE_EQ(dxdt[i], -(state.x[0][i] - xk[i]) / state.dt[0]);
    }
    EXPECT_DOUBLE_EQ(val, 1.0 / state.dt[0]);
}

/*
 * Tests BDF-I time derivative approximation
 */
TEST(TimeIntegrator, SchemeI)
{
    constexpr std::size_t size{3};

    eivec dxdt(size);
    rvec xk(size);
    SolverState state(size);

    state.order = 1;

    double val{};
    double h0{};

    // Equation 1:
    auto f1 = [](const double t)
    {
        return -7.5 * t + 4.2;
    };

    // Equation 2:
    auto f2 = [](const double t)
    {
        return 6.5 * t + 2.2;
    };

    // Equation 3:
    auto f3 = [](const double t)
    {
        return -1.33 * t - 1.1;
    };

    // Test matrix
    std::vector<double> dt_list = {1.23e-8, 1.23e-6, 1.23e-4, 1.23e-2, 1.23, 1.23e+2};

    for (const auto dt : dt_list)
    {
        state.dt[0] = dt;

        h0 = state.dt[0];

        xk = {f1(0), f2(0), f3(0)};
        state.x[0] = {f1(h0), f2(h0), f3(h0)};

        val = time_derivative_approx(dxdt, xk, state, size);

        EXPECT_NEAR(dxdt[0], 7.5, abs_err) << "Eq.1; dt=" << dt;
        EXPECT_NEAR(dxdt[1], -6.5, abs_err) << "Eq.2; dt=" << dt;
        EXPECT_NEAR(dxdt[2], 1.33, abs_err) << "Eq.3; dt=" << dt;
    }
}

/*
 * Tests BDF-II time derivative approximation
 */
TEST(TimeIntegrator, SchemeII)
{
    constexpr std::size_t size{3};

    eivec dxdt(size);
    rvec xk(size);
    SolverState state(size);

    state.order = 2;

    double val{};
    double h0{}, h1{}, h01{};

    // Equation 1:
    auto f1 = [](const double t)
    {
        return 1.45 * t * t - 7.5 * t + 4.2;
    };

    // Equation 2:
    auto f2 = [](const double t)
    {
        return -1.45 * t * t + 6.5 * t + 2.2;
    };

    // Equation 3:
    auto f3 = [](const double t)
    {
        return 0.11 * t * t - 1.33 * t - 1.1;
    };

    // Test matrix
    std::vector<double> dt_list = {1.23e-6, 1.23e-5, 1.23e-4, 1.23e-3, 1.23e-2, 1.23e-1, 1.23, 12.3};

    for (const auto dt0 : dt_list)
    {
        for (const auto dt1 : dt_list)
        {
            state.dt[0] = dt0;
            state.dt[1] = dt1;

            h0 = state.dt[0];
            h1 = state.dt[1];
            h01 = h0 + h1;

            xk = {f1(0), f2(0), f3(0)};
            state.x[0] = {f1(h0), f2(h0), f3(h0)};
            state.x[1] = {f1(h01), f2(h01), f3(h01)};

            val = time_derivative_approx(dxdt, xk, state, size);

            EXPECT_NEAR(dxdt[0], 7.5, abs_err) << "Eq.1; dt0=" << dt0 << "; dt1=" << dt1;
            EXPECT_NEAR(dxdt[1], -6.5, abs_err) << "Eq.2; dt0=" << dt0 << "; dt1=" << dt1;
            EXPECT_NEAR(dxdt[2], 1.33, abs_err) << "Eq.3; dt0=" << dt0 << "; dt1=" << dt1;
        }
    }
}

/*
 * Tests BDF-III time derivative approximation
 */
TEST(TimeIntegrator, SchemeIII)
{
    constexpr std::size_t size{3};

    eivec dxdt(size);
    rvec xk(size);
    SolverState state(size);

    state.order = 3;

    double val{};
    double h0{}, h1{}, h2{}, h01{}, h012{};

    // Equation 1:
    auto f1 = [](const double t)
    {
        return 1.45 * t * t - 7.5 * t + 4.2;
    };

    // Equation 2:
    auto f2 = [](const double t)
    {
        return -1.45 * t * t + 6.5 * t + 2.2;
    };

    // Equation 3:
    auto f3 = [](const double t)
    {
        return 0.11 * t * t - 1.33 * t - 1.1;
    };

    // Test matrix
    std::vector<double> dt_list = {1.23e-4, 1.23e-3, 1.23e-2, 1.23e-1, 1.23};

    for (const auto dt0 : dt_list)
    {
        for (const auto dt1 : dt_list)
        {
            for (const auto dt2 : dt_list)
            {
                state.dt[0] = dt0;
                state.dt[1] = dt1;
                state.dt[2] = dt2;

                h0 = state.dt[0];
                h1 = state.dt[1];
                h2 = state.dt[2];
                h01 = h0 + h1;
                h012 = h01 + h2;

                xk = {f1(0), f2(0), f3(0)};
                state.x[0] = {f1(h0), f2(h0), f3(h0)};
                state.x[1] = {f1(h01), f2(h01), f3(h01)};
                state.x[2] = {f1(h012), f2(h012), f3(h012)};

                val = time_derivative_approx(dxdt, xk, state, size);

                EXPECT_NEAR(dxdt[0], 7.5, abs_err) << "Eq.1; dt0=" << dt0 << "; dt1=" << dt1 << "; dt2=" << dt2;
                EXPECT_NEAR(dxdt[1], -6.5, abs_err) << "Eq.2; dt0=" << dt0 << "; dt1=" << dt1 << "; dt2=" << dt2;
                EXPECT_NEAR(dxdt[2], 1.33, abs_err) << "Eq.3; dt0=" << dt0 << "; dt1=" << dt1 << "; dt2=" << dt2;
            }
        }
    }
}

/*
 * Tests BDF-IV time derivative approximation
 */
TEST(TimeIntegrator, SchemeIV)
{
    constexpr std::size_t size{3};

    eivec dxdt(size);
    rvec xk(size);
    SolverState state(size);

    state.order = 4;

    double val{};
    double h0{}, h1{}, h2{}, h3{}, h01{}, h012{}, h0123{};

    // Equation 1:
    auto f1 = [](const double t)
    {
        return 1.45 * t * t - 7.5 * t + 4.2;
    };

    // Equation 2:
    auto f2 = [](const double t)
    {
        return -1.45 * t * t + 6.5 * t + 2.2;
    };

    // Equation 3:
    auto f3 = [](const double t)
    {
        return 0.11 * t * t - 1.33 * t - 1.1;
    };

    // Test matrix
    std::vector<double> dt_list = {1.23e-3, 1.23e-2, 1.23e-1, 1.23};

    for (const auto dt0 : dt_list)
    {
        for (const auto dt1 : dt_list)
        {
            for (const auto dt2 : dt_list)
            {
                for (const auto dt3 : dt_list)
                {
                    state.dt[0] = dt0;
                    state.dt[1] = dt1;
                    state.dt[2] = dt2;
                    state.dt[3] = dt3;

                    h0 = state.dt[0];
                    h1 = state.dt[1];
                    h2 = state.dt[2];
                    h3 = state.dt[3];
                    h01 = h0 + h1;
                    h012 = h01 + h2;
                    h0123 = h012 + h3;

                    xk = {f1(0), f2(0), f3(0)};
                    state.x[0] = {f1(h0), f2(h0), f3(h0)};
                    state.x[1] = {f1(h01), f2(h01), f3(h01)};
                    state.x[2] = {f1(h012), f2(h012), f3(h012)};
                    state.x[3] = {f1(h0123), f2(h0123), f3(h0123)};

                    val = time_derivative_approx(dxdt, xk, state, size);

                    // TODO: The scheme does not like relatively large h0, hence 3.0 * abs_err
                    EXPECT_NEAR(dxdt[0], 7.5, 3.0 * abs_err) << "Eq.1; dt0=" << dt0 << "; dt1=" << dt1 << "; dt2=" << dt2 << "; dt3=" << dt3;
                    EXPECT_NEAR(dxdt[1], -6.5, 3.0 * abs_err) << "Eq.2; dt0=" << dt0 << "; dt1=" << dt1 << "; dt2=" << dt2 << "; dt3=" << dt3;
                    EXPECT_NEAR(dxdt[2], 1.33, 3.0 * abs_err) << "Eq.3; dt0=" << dt0 << "; dt1=" << dt1 << "; dt2=" << dt2 << "; dt3=" << dt3;
                }
            }
        }
    }
}

} // namespace
