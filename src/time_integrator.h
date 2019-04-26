/*
 * TODO: Description of the class
 */

#pragma once

#include <iostream>

#include <mkl_spblas.h>

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

    // Descriptor of sparse mass matrix properties
    struct matrix_descr m_descrA;

    // Structure with sparse mass matrix stored in CSR format
    sparse_matrix_t m_csrA;

public:
    /*
     * TODO: Description
     */
    TimeIntegrator(RHS &rhs, Jacobian &jac, MassMatrix &mass,
                   SolverOptions &opt, const MKL_INT size)
        : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
    {
        // Reserve memory for at least 3-diagonal mass matrix
        m_M.A.reserve(3 * size);
        m_M.ja.reserve(3 * size);
        m_M.ia.reserve(size + 1);

        // Get static mass matrix
        m_mass(m_M);

        // TODO: MKL sparse matrix check

        m_descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
        m_descrA.mode = SPARSE_FILL_MODE_UPPER;
        m_descrA.diag = SPARSE_DIAG_NON_UNIT;

        // TODO: check output
        // double precision
        mkl_sparse_d_create_csr(&m_csrA, SPARSE_INDEX_BASE_ONE, size, size,
                                m_M.ia.data(), m_M.ia.data() + 1, m_M.ja.data(),
                                m_M.A.data());

        // TODO: Analyze sparse matrix; choose proper kernels and workload
        // balancing strategy mkl_sparse_optimize(csrA);
    }

    ~TimeIntegrator() { mkl_sparse_destroy(m_csrA); }

    /*
     * x_prev is a C-style matrix containing history of the previous states
     */
    void operator()(sparse_matrix_holder &Jt, state_type &b, state_type &x,
                    const state_type x_prev[], const double t, const double dt);
};

}  // namespace daecpp_namespace_name
