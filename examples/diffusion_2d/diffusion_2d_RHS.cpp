/*
 * The RHS of the system. Finite Volume method.
 */

#include "diffusion_2d_RHS.h"

void MyRHS::operator()(const daecpp::state_type &x, daecpp::state_type &f,
                       const double t)
{
    // Locals
    const MKL_INT N     = m_N;
    const double  invh2 = (double)(N * N);

    daecpp::float_type Fr, Fl, Fu, Fd;

    // clang-format off

    for(MKL_INT i = 0; i < N; i++)
    {
        for(MKL_INT j = 0; j < N; j++)
        {
            MKL_INT ind = j + i*N;

            if(j != N-1)
                Fr = (x[ind+1] - x[ind]);
            else
                Fr = 0;

            if(j != 0)
                Fl = (x[ind] - x[ind-1]);
            else
                Fl = 0;

            if(i != N-1)
                Fu = (x[ind+N] - x[ind]);
            else
                Fu = 0;

            if(i != 0)
                Fd = (x[ind] - x[ind-N]);
            else
                Fd = 0;

            f[ind] = m_D*(Fr - Fl + Fu - Fd)*invh2;
        }
    }

    // clang-format on
}
