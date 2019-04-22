/*
 * TODO: Description
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

void TimeIntegrator::operator()(state_type &x, state_type x_prev[], double &t)
{
    const int size = (int)(x.size());

    const double invdt = 1.0 / m_dt;

    state_type dxdt(size);

    if(m_scheme == 1)  // first order BDF
    {
        for(int i = 0; i < size; i++)
        {
            dxdt[i] = (x[i] - x_prev[i][0]) * invdt;
        }
    }
    else
    {
        // internal error
    }

    state_type M(size * size);  // Overkill! TODO: this should flexible

    vector_type_int ia(size + 1), ja(size * size);  // Same here

    m_mass(M, ia, ja, t);

    // vector b = sparse matrix M * vector dxdt

    // vector x -= vector b

    // sparse matrix Mdt = sparse matrix M * invdt element-wise

    // sparse matrix J = sparse matrix J - sparse matrix Mdt
}

}  // namespace daecpp_namespace_name
