/*
 * Testing:
 * struct sparse_matrix
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#include <dae-cpp/sparse-matrix.hpp>

#include "gtest/gtest.h"

namespace
{

using namespace daecpp;

TEST(SparseMatrix, AddElement)
{
    sparse_matrix M;

    M.add_element(1, 2, 1.5);
    M.add_element(3, 4, 2.5);
    M.add_element(5, 6, 3.5);

    EXPECT_DOUBLE_EQ(M.A[0], 1.5);
    EXPECT_DOUBLE_EQ(M.A[1], 2.5);
    EXPECT_DOUBLE_EQ(M.A[2], 3.5);

    EXPECT_EQ(M.i[0], 1);
    EXPECT_EQ(M.i[1], 3);
    EXPECT_EQ(M.i[2], 5);

    EXPECT_EQ(M.j[0], 2);
    EXPECT_EQ(M.j[1], 4);
    EXPECT_EQ(M.j[2], 6);

    EXPECT_EQ(M.A.size(), 3);
    EXPECT_EQ(M.i.size(), 3);
    EXPECT_EQ(M.j.size(), 3);
}

TEST(SparseMatrix, AddElementOperator)
{
    sparse_matrix M;

    // Using operator () instead of add_element
    M(1, 2, 1.5);
    M(3, 4, 2.5);
    M(5, 6, 3.5);
    M(5, 6, 3.5); // Duplicate

    EXPECT_DOUBLE_EQ(M.A[0], 1.5);
    EXPECT_DOUBLE_EQ(M.A[1], 2.5);
    EXPECT_DOUBLE_EQ(M.A[2], 3.5);
    EXPECT_DOUBLE_EQ(M.A[3], 3.5);

    EXPECT_EQ(M.i[0], 1);
    EXPECT_EQ(M.i[1], 3);
    EXPECT_EQ(M.i[2], 5);
    EXPECT_EQ(M.i[3], 5);

    EXPECT_EQ(M.j[0], 2);
    EXPECT_EQ(M.j[1], 4);
    EXPECT_EQ(M.j[2], 6);
    EXPECT_EQ(M.j[3], 6);

    EXPECT_EQ(M.A.size(), 4);
    EXPECT_EQ(M.i.size(), 4);
    EXPECT_EQ(M.j.size(), 4);
}

TEST(SparseMatrix, Reserve)
{
    sparse_matrix M;

    M.reserve(2);

    M.add_element(1, 2, 1.5);
    M.add_element(3, 4, 2.5);

    EXPECT_DOUBLE_EQ(M.A[0], 1.5);
    EXPECT_DOUBLE_EQ(M.A[1], 2.5);

    EXPECT_EQ(M.i[0], 1);
    EXPECT_EQ(M.i[1], 3);

    EXPECT_EQ(M.j[0], 2);
    EXPECT_EQ(M.j[1], 4);

    EXPECT_EQ(M.A.size(), 2);
    EXPECT_EQ(M.i.size(), 2);
    EXPECT_EQ(M.j.size(), 2);
}

TEST(SparseMatrix, Clear)
{
    sparse_matrix M;

    M.reserve(2);

    M.add_element(1, 2, 1.5);
    M.add_element(3, 4, 2.5);

    M.clear();

    EXPECT_EQ(M.A.size(), 0);
    EXPECT_EQ(M.i.size(), 0);
    EXPECT_EQ(M.j.size(), 0);

    M.add_element(1, 1, 12.5);

    EXPECT_EQ(M.A.size(), 1);
    EXPECT_EQ(M.i.size(), 1);
    EXPECT_EQ(M.j.size(), 1);
}

TEST(SparseMatrix, Check)
{
    sparse_matrix M;

    M.reserve(3);

    M.add_element(0, 1, 1e-6);
    M.add_element(10, 0, 1e3);
    M.add_element(10, 0, 2e3); // Duplicate element is OK

    EXPECT_DOUBLE_EQ(M.A[0], 1e-6);
    EXPECT_DOUBLE_EQ(M.A[1], 1000.0);
    EXPECT_DOUBLE_EQ(M.A[2], 2000.0);

    EXPECT_EQ(M.i[0], 0);
    EXPECT_EQ(M.i[1], 10);
    EXPECT_EQ(M.i[2], 10);

    EXPECT_EQ(M.j[0], 1);
    EXPECT_EQ(M.j[1], 0);
    EXPECT_EQ(M.j[2], 0);

    M.check();
}

TEST(SparseMatrix, CheckZero)
{
    sparse_matrix M; // Zero matrix

    M.check();

    EXPECT_EQ(M.A.size(), 0);
    EXPECT_EQ(M.i.size(), 0);
    EXPECT_EQ(M.j.size(), 0);
}

TEST(SparseMatrix, Dense)
{
    sparse_matrix M;

    M.reserve(3);
    M.add_element(0, 0, 1.0);
    M.add_element(0, 1, 2.0);
    M.add_element(1, 0, 3.0);

    M.check();

    auto A = M.dense(2);

    EXPECT_DOUBLE_EQ(A.coeffRef(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A.coeffRef(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(A.coeffRef(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(A.coeffRef(1, 1), 0.0); // It wasn't defined, should be 0

    // Alternatively
    EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(A(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(A(1, 1), 0.0);
}

TEST(SparseMatrix, Convert)
{
    sparse_matrix M1, M2, M3;

    M1.reserve(4);
    M1.add_element(0, 0, 1.0);
    M1.add_element(0, 1, 2.0);
    M1.add_element(1, 0, 3.0);
    M1.add_element(1, 1, 4.0);

    M2.reserve(3);
    M2.add_element(0, 0, -5.0);
    M2.add_element(1, 1, 2.5);
    M2.add_element(1, 1, 2.5); // Duplicate element -- should be summed up

    M1.check();
    M2.check();
    M3.check(); // Empty (zero) matrix

    auto A1 = M1.convert(2);
    auto A2 = M2.convert(2);
    auto A3 = M3.convert(2);

    EXPECT_DOUBLE_EQ(A1.toDense().coeffRef(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A1.toDense().coeffRef(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(A1.toDense().coeffRef(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(A1.toDense().coeffRef(1, 1), 4.0);

    EXPECT_DOUBLE_EQ(A2.toDense().coeffRef(0, 0), -5.0);
    EXPECT_DOUBLE_EQ(A2.toDense().coeffRef(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(A2.toDense().coeffRef(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(A2.toDense().coeffRef(1, 1), 5.0);

    EXPECT_DOUBLE_EQ(A3.toDense().coeffRef(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(A3.toDense().coeffRef(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(A3.toDense().coeffRef(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(A3.toDense().coeffRef(1, 1), 0.0);

    auto SUM = 2.0 * A1 + A2 + 10.0 * A3;

    EXPECT_DOUBLE_EQ(SUM.toDense().coeffRef(0, 0), -3.0);
    EXPECT_DOUBLE_EQ(SUM.toDense().coeffRef(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(SUM.toDense().coeffRef(1, 0), 6.0);
    EXPECT_DOUBLE_EQ(SUM.toDense().coeffRef(1, 1), 13.0);

    EXPECT_DOUBLE_EQ(SUM.sum(), 20);
}

TEST(SparseMatrix, NonZeros)
{
    sparse_matrix M;

    M.reserve(4);
    M.add_element(0, 0, 1.0);
    M.add_element(0, 1, 2.0);
    M.add_element(1, 0, 3.0);
    M.add_element(1, 1, 4.0);

    M.check();

    EXPECT_EQ(M.N_elements(), 4);

    sparse_matrix Z; // Empty matrix

    Z.check();

    EXPECT_EQ(Z.N_elements(), 0);
}

} // namespace
