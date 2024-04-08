/*
 * Testing:
 * class MassMatrix, MassMatrixIdentity, MassMatrixZero
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#include <dae-cpp/mass-matrix.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

TEST(MassMatrix, Definition)
{
    struct TestMassMatrix : MassMatrix
    {
        void operator()(sparse_matrix &M, const double t) const
        {
            EXPECT_EQ(M.N_elements(), 0);

            M.reserve(2);
            M(0, 0, 1.0);
            M(1, 1, 2.0 * t);
        }
    };

    TestMassMatrix mass;
    sparse_matrix M;

    constexpr double t{10.0};

    mass(M, t);

    M.check();

    EXPECT_DOUBLE_EQ(M.A[0], 1.0);
    EXPECT_DOUBLE_EQ(M.A[1], 2.0 * t);

    EXPECT_EQ(M.i[0], 0);
    EXPECT_EQ(M.i[1], 1);

    EXPECT_EQ(M.j[0], 0);
    EXPECT_EQ(M.j[1], 1);

    EXPECT_EQ(M.N_elements(), 2);
}

TEST(MassMatrix, Identity)
{
    constexpr double N{1000};
    constexpr double t{10.0};

    MassMatrixIdentity mass(N);
    sparse_matrix M;

    mass(M, t);

    M.check();

    for (std::size_t i = 0; i < N; ++i)
    {
        EXPECT_DOUBLE_EQ(M.A[i], 1.0);
        EXPECT_EQ(M.i[i], i);
        EXPECT_EQ(M.j[i], i);
    }

    EXPECT_EQ(M.N_elements(), N);
}

TEST(MassMatrix, Zero)
{
    constexpr double N{1000};
    constexpr double t{10.0};

    MassMatrixZero mass;
    sparse_matrix M;

    mass(M, t);

    M.check();

    EXPECT_EQ(M.N_elements(), 0);
}

} // namespace
