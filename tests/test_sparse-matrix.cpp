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

} // namespace
