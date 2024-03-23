#include <dae-cpp/sparse-matrix.hpp>

#include "gtest/gtest.h"

// Testing:
// struct sparse_matrix

namespace
{

using namespace daecpp;

TEST(SparseMatrix, AddElement)
{
    sparse_matrix M;

    M.add_element(1.5, 1, 2);
    M.add_element(2.5, 3, 4);
    M.add_element(3.5, 5, 6);

    EXPECT_DOUBLE_EQ(M.A[0], 1.5);
    EXPECT_DOUBLE_EQ(M.A[1], 2.5);
    EXPECT_DOUBLE_EQ(M.A[2], 3.5);

    EXPECT_EQ(M.i[0], 1);
    EXPECT_EQ(M.i[1], 3);
    EXPECT_EQ(M.i[2], 5);

    EXPECT_EQ(M.j[0], 2);
    EXPECT_EQ(M.j[1], 4);
    EXPECT_EQ(M.j[2], 6);
}

TEST(SparseMatrix, Reserve)
{
    sparse_matrix M;

    M.reserve(2);

    M.add_element(1.5, 1, 2);
    M.add_element(2.5, 3, 4);

    EXPECT_DOUBLE_EQ(M.A[0], 1.5);
    EXPECT_DOUBLE_EQ(M.A[1], 2.5);

    EXPECT_EQ(M.i[0], 1);
    EXPECT_EQ(M.i[1], 3);

    EXPECT_EQ(M.j[0], 2);
    EXPECT_EQ(M.j[1], 4);
}

TEST(SparseMatrix, Check)
{
    sparse_matrix M;

    M.reserve(2);

    M.add_element(1e-6, 0, 1);
    M.add_element(1e3, 10, 0);

    EXPECT_DOUBLE_EQ(M.A[0], 1e-6);
    EXPECT_DOUBLE_EQ(M.A[1], 1000.0);

    EXPECT_EQ(M.i[0], 0);
    EXPECT_EQ(M.i[1], 10);

    EXPECT_EQ(M.j[0], 1);
    EXPECT_EQ(M.j[1], 0);

    M.check();
}

TEST(SparseMatrix, CheckZero)
{
    sparse_matrix M;

    M.check();
}

TEST(SparseMatrix, Convert)
{
    sparse_matrix M1, M2;

    M1.reserve(4);
    M1.add_element(1.0, 0, 0);
    M1.add_element(2.0, 0, 1);
    M1.add_element(3.0, 1, 0);
    M1.add_element(4.0, 1, 1);

    M2.reserve(2);
    M2.add_element(-5.0, 0, 0);
    M2.add_element(5.0, 1, 1);

    M1.check();
    M2.check();

    auto A1 = M1.convert(2);
    auto A2 = M2.convert(2);

    auto A3 = A1 + A2;
    auto A4 = 2.0 * A1 + A2;

    EXPECT_DOUBLE_EQ(A3.sum(), 10);
    EXPECT_DOUBLE_EQ(A4.sum(), 20);
}

TEST(SparseMatrix, Dense)
{
    sparse_matrix M1, M2;

    M1.reserve(4);
    M1.add_element(1.0, 0, 0);
    M1.add_element(2.0, 0, 1);
    M1.add_element(3.0, 1, 0);
    M1.add_element(4.0, 1, 1);

    M2.reserve(2);
    M2.add_element(-5.0, 0, 0);
    M2.add_element(5.0, 1, 1);

    M1.check();
    M2.check();

    auto A1 = M1.dense(2);
    auto A2 = M2.dense(2);

    auto A3 = A1 + A2;
    auto A4 = 2.0 * A1 + A2;

    EXPECT_DOUBLE_EQ(A3.sum(), 10);
    EXPECT_DOUBLE_EQ(A4.sum(), 20);
}

} // namespace
