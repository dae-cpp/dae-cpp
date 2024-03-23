#include <dae-cpp/rhs.hpp>

#include "gtest/gtest.h"

// Testing:
// class RHS, RHS_empty

namespace
{

using namespace daecpp;

TEST(RHS, Definition)
{
    struct TestRHS : public RHS
    {
        void operator()(const state_type &x, state_type &f, const double t) const
        {
            EXPECT_EQ(x.size(), 2);

            f[0] = x[0];
            f[1] = x[1] * t;
        }
    };

    TestRHS rhs;

    state_type x{4.0, 6.0};
    state_type f(2);

    constexpr double t{10.0};

    rhs(x, f, t);

    EXPECT_DOUBLE_EQ(f[0], 4.0);
    EXPECT_DOUBLE_EQ(f[1], 60.0);

    EXPECT_EQ(f.size(), 2);
}

TEST(RHS, EmptyRHS)
{
    core::RHS_empty rhs;
}

} // namespace
