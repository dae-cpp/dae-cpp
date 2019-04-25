/*
 * TODO: Description of the class
 */

#pragma once

#include <iostream>
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

    // Mass matrix container
    sparse_matrix_holder m_M;

public:
    /*
     * TODO: Description
     */
    TimeIntegrator(RHS &rhs, Jacobian &jac, MassMatrix &mass,
                   SolverOptions &opt, const MKL_INT N)
        : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
    {
        // assert m_opt.N > 0
        m_M.A.reserve(m_opt.max_size_mult * N);
        m_M.ia.reserve(N + 1);
        m_M.ja.reserve(m_opt.max_size_mult * N);

        state_type::size_type sz = m_M.A.capacity();

        m_mass(m_M);

        // MKL sparse matrix check

        // not important for mass matrix - can be deleted
        if(sz != m_M.A.capacity())
        {
            // test this
            double mult = (double)m_M.A.capacity() / (double)sz;
            std::cout << "\nWarning: Mass matrix capacity changed. Suggested "
                         "multiplyer max_size_mult should be at least "
                      << mult * m_opt.max_size_mult << '\n';
        }
    }

    /*
     * x_prev is a C-style matrix containing history of the previous states
     */
    void operator()(sparse_matrix_holder &Jt, state_type &b, state_type &x,
                    state_type x_prev[], double &t, const double dt);
};

}  // namespace daecpp_namespace_name
