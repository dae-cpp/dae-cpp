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

    EXPECT_EQ(M.A.size(), 3);
    EXPECT_EQ(M.i.size(), 3);
    EXPECT_EQ(M.j.size(), 3);
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

    EXPECT_EQ(M.A.size(), 2);
    EXPECT_EQ(M.i.size(), 2);
    EXPECT_EQ(M.j.size(), 2);
}

TEST(SparseMatrix, Clear)
{
    sparse_matrix M;

    M.reserve(2);

    M.add_element(1.5, 1, 2);
    M.add_element(2.5, 3, 4);

    M.clear();

    EXPECT_EQ(M.A.size(), 0);
    EXPECT_EQ(M.i.size(), 0);
    EXPECT_EQ(M.j.size(), 0);

    M.add_element(12.5, 1, 1);

    EXPECT_EQ(M.A.size(), 1);
    EXPECT_EQ(M.i.size(), 1);
    EXPECT_EQ(M.j.size(), 1);
}

TEST(SparseMatrix, Check)
{
    sparse_matrix M;

    M.reserve(3);

    M.add_element(1e-6, 0, 1);
    M.add_element(1e3, 10, 0);
    M.add_element(2e3, 10, 0); // Duplicate element is OK

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

    EXPECT_EQ(M.A.size(), 0);
    EXPECT_EQ(M.i.size(), 0);
    EXPECT_EQ(M.j.size(), 0);
}

TEST(SparseMatrix, Dense)
{
    sparse_matrix M;

    M.reserve(4);
    M.add_element(1.0, 0, 0);
    M.add_element(2.0, 0, 1);
    M.add_element(3.0, 1, 0);
    M.add_element(4.0, 1, 1);

    M.check();

    auto A = M.dense(2);

    EXPECT_DOUBLE_EQ(A.coeffRef(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A.coeffRef(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(A.coeffRef(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(A.coeffRef(1, 1), 4.0);
}

TEST(SparseMatrix, Convert)
{
    sparse_matrix M1, M2, M3;

    M1.reserve(4);
    M1.add_element(1.0, 0, 0);
    M1.add_element(2.0, 0, 1);
    M1.add_element(3.0, 1, 0);
    M1.add_element(4.0, 1, 1);

    M2.reserve(3);
    M2.add_element(-5.0, 0, 0);
    M2.add_element(2.5, 1, 1);
    M2.add_element(2.5, 1, 1); // Duplicate element -- should be summed up

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
    M.add_element(1.0, 0, 0);
    M.add_element(2.0, 0, 1);
    M.add_element(3.0, 1, 0);
    M.add_element(4.0, 1, 1);

    M.check();

    EXPECT_EQ(M.N_elements(), 4);

    sparse_matrix Z; // Empty matrix

    Z.check();

    EXPECT_EQ(Z.N_elements(), 0);
}

} // namespace
