/*
 * The RHS of the system. Finite Volume method.
 */

#include "diffusion_2d_RHS.h"

void MyRHS::operator()(const daecpp::state_type &x, daecpp::state_type &f,
                       const double t)
{
    // Locals
    const MKL_INT N     = m_N;
    const double  invh2 = (double)N * (double)N;

    // clang-format off
    for(MKL_INT i = 0; i < N; i++)
    {
        for(MKL_INT j = 0; j < N; j++)
        {
            MKL_INT ind = j + i*N;

            // Neumann boundary conditions
            double Fr = (j != N-1) ? (x[ind+1] - x[ind]) : 0.0;
            double Fu = (i != N-1) ? (x[ind+N] - x[ind]) : 0.0;
            double Fl = (j != 0)   ? (x[ind] - x[ind-1]) : 0.0;
            double Fd = (i != 0)   ? (x[ind] - x[ind-N]) : 0.0;

            f[ind] = m_D*(Fr - Fl + Fu - Fd)*invh2;
        }
    }
    // clang-format on
}
