/*
 * The RHS of the system. Finite Difference method.
 */

#include "perovskite_RHS.h"

void MyRHS::operator()(const daecpp::state_type &x, daecpp::state_type &f,
                       const double t)
{
    // Locals
    const MKL_INT N       = m_p.N;
    const double  invh2   = m_p.invh * m_p.invh;
    const double  invlam2 = 1.0 / (m_p.lambda * m_p.lambda);

    // clang-format off

    // RHS for the ion concentration P = dFlux/dx
    for(MKL_INT i = 1; i < N-1; i++)
    {
        f[i] = (x[i+1] - 2.0*x[i] + x[i-1] + 0.5*((x[i+1] + x[i])*(x[N+i+1] - x[N+i]) - (x[i] + x[i-1])*(x[N+i] - x[N+i-1])))*invh2;
    }
    f[0] = (x[1] - x[0] + 0.5*(x[1] + x[0])*(x[N+1] - x[N]))*invh2;                   // Left BC
    f[N-1] = -(x[N-1] - x[N-2] + 0.5*(x[N-1] + x[N-2])*(x[2*N-1] - x[2*N-2]))*invh2;  // Right BC

    // RHS for the potential Phi
    for(MKL_INT i = 1; i < N-1; i++)
    {
        f[i+N] = (x[i+1+N] - 2.0*x[i+N] + x[i-1+N])*invh2 - (1.0 - x[i])*invlam2;
    }
    f[N] = x[N] + t;          // Left BC
    f[2*N-1] = x[2*N-1] - t;  // Right BC

    // clang-format on
}

/*
 * This is for the demonstration purpose.
 * We can force the solver to stop earlier by overriding
 * stop_condition(const daecpp::state_type &x, const double t) function.
 */
bool MyRHS::stop_condition(const daecpp::state_type &x, const double t)
{
    if(x[0] < 0)
    {
        return true;  // if x[0] is less than 0 (should never happen), then the
                      // solver should stop
    }
    else
    {
        return false;
    }
}
