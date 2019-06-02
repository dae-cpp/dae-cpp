/*
 * The RHS of the system. Finite Volume method.
 */

#include "diffusion_2d_RHS.h"

void MyRHS::operator()(const daecpp::state_type &x, daecpp::state_type &f,
                       const double t)
{
    // Number of cells along x and y axis
    const MKL_INT N = m_N;

    // Using `daecpp::float_type` here instead of `double` to support single
    // precision (only if DAE_SINGLE=ON). But in most cases this is not
    // necessary, one can use `double` only since `daecpp::float_type` is
    // equivalent to `double` by default (DAE_SINGLE=OFF).
    const daecpp::float_type dinvh2 = m_D * (double)(N) * (double)(N);

    // Flux value on the boundaries
    const daecpp::float_type bflux = 0.0;

    // clang-format off
    for(MKL_INT i = 0; i < N; i++)
    {
        for(MKL_INT j = 0; j < N; j++)
        {
            MKL_INT ind = j + i*N;

            // Neumann boundary conditions.
            // A more efficient (performance-wise) way to implement BCs would be
            // moving them into a separate loop over all boundary cells in order
            // to avoid if-conditions within the loop. But the way suggested
            // bellow is less messy and much easier to read/modify for a human.
            const daecpp::float_type Fr = (j != N-1) ? (x[ind+1] - x[ind]) : bflux;  // Right
            const daecpp::float_type Fu = (i != N-1) ? (x[ind+N] - x[ind]) : bflux;  // Up
            const daecpp::float_type Fl = (j != 0)   ? (x[ind] - x[ind-1]) : bflux;  // Left
            const daecpp::float_type Fd = (i != 0)   ? (x[ind] - x[ind-N]) : bflux;  // Down

            f[ind] = (Fr - Fl + Fu - Fd)*dinvh2;
        }
    }
    // clang-format on
}
