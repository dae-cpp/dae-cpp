/*
 * TODO: Description of the class
 */

#include "perovskite_RHS.h"

void MyRHS::operator()(const daecpp::state_type &x, daecpp::state_type &f,
                       const double t)
{
    // Locals
    const int    N       = m_p.N;  // const int N = (int)(x.size()) / 2;
    const double invh2   = m_p.invh * m_p.invh;
    const double invlam2 = 1.0 / (m_p.lambda * m_p.lambda);

    // RHS for P = dFlux/dx
#pragma omp parallel for
    for(int i = 1; i < N - 1; i++)
    {
        // clang-format off
        f[i] = (x[i+1] - 2.0*x[i] + x[i-1] + 0.5*((x[i+1] + x[i])*(x[N+i+1] - x[N+i]) - (x[i] + x[i-1])*(x[N+i] - x[N+i-1])))*invh2;
        // clang-format on
    }
    // clang-format off
    f[0] = (x[1] - x[0] + 0.5*(x[1] + x[0])*(x[N+1] - x[N]))*invh2;  // Left BC
    f[N-1] = -(x[N-1] - x[N-2] + 0.5*(x[N-1] + x[N-2])*(x[2*N-1] - x[2*N-2]))*invh2;  // Right BC
    // clang-format on

    // RHS for Phi
#pragma omp parallel for
    for(int i = 1; i < N - 1; i++)
    {
        // clang-format off
        f[i+N] = (x[i+1+N] - 2.0*x[i+N] + x[i-1+N])*invh2 - (1.0 - x[i])*invlam2;
        // clang-format on
    }
    // clang-format off
    f[N] = (x[N+1] - 3.0*x[N] - 2.0*t)*invh2 - (1.0 - x[0])*invlam2;  // Left BC
    f[2*N-1] = (2.0*t - 3.0*x[2*N-1] + x[2*N-2])*invh2 - (1.0 - x[N-1])*invlam2;  // Right BC
    // clang-format on
}
