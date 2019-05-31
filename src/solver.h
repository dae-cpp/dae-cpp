/*
 * The main solver class
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"
#include "jacobian.h"
#include "mass_matrix.h"
#include "solver_options.h"
#include "time_integrator.h"

namespace daecpp_namespace_name
{

class Solver
{
    RHS &m_rhs;  // RHS

    Jacobian &m_jac;  // Jacobian matrix

    MassMatrix &m_mass;  // Mass matrix

    SolverOptions &m_opt;  // Solver options

    void check_pardiso_error(MKL_INT err);

public:
    /*
     * Receives user-defined RHS, Jacobian, Mass matrix and solver options
     */
    Solver(RHS &rhs, Jacobian &jac, MassMatrix &mass, SolverOptions &opt)
        : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
    {
    }

    /*
     * Integrates the system of DAEs on the interval t = [t0; t1] and returns
     * result in the array x. Parameter t0 can be overriden in the solver
     * options (t0 = 0 by default).
     * The data stored in x (initial conditions) will be overwritten.
     */
    void operator()(state_type &x, const double t1);
};

}  // namespace daecpp_namespace_name
