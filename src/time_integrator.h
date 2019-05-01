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

    /*
    * Sparse matrix checker
    */
    int matrix_checker(sparse_matrix_holder &A, MKL_INT size);

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

        // User defined sparse matrix check
        if(matrix_checker(m_M, size))
            exit(1);

        // Init mass matrix descriptor for sparse solver
        m_descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
        m_descrA.mode = SPARSE_FILL_MODE_UPPER;
        m_descrA.diag = SPARSE_DIAG_NON_UNIT;

        // Create sparse matrix descriptor
        sparse_status_t sp_status = mkl_sparse_d_create_csr(&m_csrA, SPARSE_INDEX_BASE_ONE, size, size,
                                m_M.ia.data(), m_M.ia.data() + 1, m_M.ja.data(),
                                m_M.A.data());

        if(sp_status != SPARSE_STATUS_SUCCESS)
        {
            std::cout << "ERROR: Could not create sparse matrix descriptor for Mass matrix.\n";
            exit(1);
        }

        // Analyze sparse matrix, choose proper kernels and workload
        // balancing strategy
        sparse_status_t opt_status = mkl_sparse_optimize(m_csrA);

        if(opt_status != SPARSE_STATUS_SUCCESS)
        {
            std::cout << "WARNING: Sparse matrix analysis and optimisation failed.\n";
        }
    }

    ~TimeIntegrator() { mkl_sparse_destroy(m_csrA); }

    /*
     * x_prev is a C-style matrix containing history of the previous states
     */
    void operator()(sparse_matrix_holder &Jt, state_type &b, state_type &x,
                    const state_type x_prev[], const double t, const double dt);
};

}  // namespace daecpp_namespace_name
