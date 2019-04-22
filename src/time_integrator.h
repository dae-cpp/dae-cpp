/*
 * TODO: Description of the class
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"
#include "jacobian.h"
#include "mass_matrix.h"
#include "solver_options.h"

namespace daecpp_namespace_name
{

class TimeIntegrator
{
    RHS &m_rhs;

    Jacobian &m_jac;

    MassMatrix &m_mass;

    SolverOptions &m_opt;

    int m_scheme = 1;

    double m_dt;

public:
    /*
     * TODO: Description
     */
    TimeIntegrator(RHS &rhs, Jacobian &jac, MassMatrix &mass,
                   SolverOptions &opt)
        : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt), m_dt(opt.dt_init)
    {
    }

    /*
     * x_prev is a C-style matrix containing history of the previous states
     */
    void operator()(state_type &x, state_type x_prev[], double &t);
};

}  // namespace daecpp_namespace_name
