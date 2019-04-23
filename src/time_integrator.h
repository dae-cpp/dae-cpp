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
    state_type m_M;
    vector_type_int m_ia;
    vector_type_int m_ja;

public:
    /*
     * TODO: Description
     */
    TimeIntegrator(RHS &rhs, Jacobian &jac, MassMatrix &mass,
                   SolverOptions &opt, const int N)
        : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
    {
        // assert m_opt.N > 0
        m_M.reserve(m_opt.max_size_mult * N);
        m_ia.reserve(N + 1);
        m_ja.reserve(m_opt.max_size_mult * N);

        state_type::size_type sz = m_M.capacity();

        m_mass(m_M, m_ia, m_ja);

        // MKL sparse matrix check

        // not important for mass matrix - can be deleted
        if(sz != m_M.capacity()) {
            // test this
            double mult = (double)m_M.capacity() / (double)sz;
            std::cout << "\nWarning: Mass matrix capacity changed. Suggested multiplyer max_size_mult should be at least " << mult*m_opt.max_size_mult << '\n';
        }
    }

    /*
     * x_prev is a C-style matrix containing history of the previous states
     */
    void operator()(state_type &Jt, vector_type_int &Jt_ia, vector_type_int &Jt_ja, state_type &x, state_type x_prev[], double &t, const double dt);
};

}  // namespace daecpp_namespace_name
