#include <dae-cpp/jacobian-matrix.hpp>

#include "gtest/gtest.h"

// Testing:
// class JacobianMatrix

namespace
{

using namespace daecpp;

TEST(JacobianMatrix, Definition)
{
    struct TestJacobian : JacobianMatrix
    {
        void operator()(sparse_matrix &J, const state_vector &x, const double t) const
        {
            EXPECT_EQ(J.N_elements(), 0);

            J.reserve(3);
            J(0, 1, 1.0);
            J(1, 0, 2.0 * t);
            J(1, 1, 3.0 * x[1]);
        }
    };

    TestJacobian jac;
    sparse_matrix J;

    state_vector x{4.0, 6.0, 8.0};

    constexpr double t{10.0};

    jac(J, x, t);

    J.check();

    EXPECT_DOUBLE_EQ(J.A[0], 1.0);
    EXPECT_DOUBLE_EQ(J.A[1], 2.0 * t);
    EXPECT_DOUBLE_EQ(J.A[2], 3.0 * x[1]);

    EXPECT_EQ(J.i[0], 0);
    EXPECT_EQ(J.i[1], 1);
    EXPECT_EQ(J.i[2], 1);

    EXPECT_EQ(J.j[0], 1);
    EXPECT_EQ(J.j[1], 0);
    EXPECT_EQ(J.j[2], 1);

    EXPECT_EQ(J.N_elements(), 3);
}

} // namespace
