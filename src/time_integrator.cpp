/*
 * Numerical time integrator. BDF method.
 */

#include <iostream>
#include <chrono>

#include "time_integrator.h"

namespace daecpp_namespace_name
{

TimeIntegrator::TimeIntegrator(RHS &rhs, Jacobian &jac, MassMatrix &mass,
                               SolverOptions &opt)
    : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
{
    // Get static mass matrix
    m_mass(m_M);

    // Convert it to Intel MKL three-array format if necessary
    m_matrix_converter(m_M);

    // Extract the mass matrix size
    const MKL_INT size = m_M.ia.size() - 1;

    // User defined sparse matrix check
    if(m_matrix_checker(m_M, size))
    {
        std::cout << "Error in Mass matrix. This error is fatal.\n";
        exit(11);
    }

    // Init mass matrix descriptor for sparse solver
    m_descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    m_descrA.mode = SPARSE_FILL_MODE_UPPER;
    m_descrA.diag = SPARSE_DIAG_NON_UNIT;

#ifdef DAE_FORTRAN_STYLE
    sparse_index_base_t style = SPARSE_INDEX_BASE_ONE;
#else
    sparse_index_base_t style = SPARSE_INDEX_BASE_ZERO;
#endif

    // Create sparse matrix descriptor
#ifdef DAE_SINGLE
    sparse_status_t sp_status =
        mkl_sparse_s_create_csr(&m_csrA, style, size, size, m_M.ia.data(),
                                m_M.ia.data() + 1, m_M.ja.data(), m_M.A.data());
#else
    sparse_status_t sp_status =
        mkl_sparse_d_create_csr(&m_csrA, style, size, size, m_M.ia.data(),
                                m_M.ia.data() + 1, m_M.ja.data(), m_M.A.data());
#endif

    if(sp_status != SPARSE_STATUS_SUCCESS)
    {
        std::cout << "ERROR: Could not create sparse matrix descriptor for "
                     "Mass matrix. This error is fatal.\n";
        exit(12);
    }

    // Analyze sparse matrix, choose proper kernels and workload
    // balancing strategy
    sparse_status_t opt_status = mkl_sparse_optimize(m_csrA);

    if(opt_status != SPARSE_STATUS_SUCCESS)
    {
        std::cout
            << "WARNING: Sparse matrix analysis and optimisation failed.\n";
    }
}

void TimeIntegrator::integrate(sparse_matrix_holder &J, state_type &b,
                               const state_type &x,
                               const state_type_matrix &x_prev, const double t,
                               const double dt[], const bool do_jac)
{
    const MKL_INT size = (MKL_INT)(x.size());

    double alpha  = 0.0;
    int    scheme = m_scheme - 1;

    state_type dxdt(size);

    // Variable time stepper for BDF-2
    if(m_scheme == 2)
    {
        for(MKL_INT i = 0; i < size; i++)
        {
            dxdt[i] = (2.0 * dt[0] + dt[1]) / (dt[0] * (dt[0] + dt[1])) * x[i] -
                      (dt[0] + dt[1]) / (dt[0] * dt[1]) * x_prev[0][i] +
                      dt[0] / (dt[1] * (dt[0] + dt[1])) * x_prev[1][i];
        }

        alpha = (2.0 * dt[0] + dt[1]) / (dt[0] * (dt[0] + dt[1]));
    }
    else  // Constant time stepper
    {
        const double invdt = 1.0 / dt[0];

        for(MKL_INT i = 0; i < size; i++)
        {
            double val = BDF_COEF[scheme][0] * x[i];
            for(int d = 1; d <= m_scheme; d++)
            {
                val += BDF_COEF[scheme][d] * x_prev[d - 1][i];
            }
            val *= invdt / BDF_COEF[scheme][7];
            dxdt[i] = val;
        }

        alpha = invdt * ALPHA_COEF[scheme];
    }

    // Initialise clock
    using clock     = std::chrono::high_resolution_clock;
    using time_unit = std::chrono::microseconds;

    // Calculate the RHS
    {
        auto tic0 = clock::now();
        m_rhs(x, b, t);
        auto tic1 = clock::now();

        // Update the RHS timer
        m_rhs_time +=
            std::chrono::duration_cast<time_unit>(tic1 - tic0).count() * 1e-6;
    }

    // b := -M * dxdt + b
#ifdef DAE_SINGLE
    mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, m_csrA, m_descrA,
                    dxdt.data(), 1.0, b.data());
#else
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, m_csrA, m_descrA,
                    dxdt.data(), 1.0, b.data());
#endif

    if(do_jac)
    {
        // Clear temporary Jacobian
        m_J.A.clear();
        m_J.ia.clear();
        m_J.ja.clear();

        // Reserve memory for at least 3-diagonal Jacobian
        m_J.A.reserve(3 * size);
        m_J.ia.reserve(size + 1);
        m_J.ja.reserve(3 * size);

        // Calculate Jacobian
        auto tic0 = clock::now();
        m_jac(m_J, x, t);
        m_matrix_converter(m_J);  // Converts Jacobian if necessary
        auto tic1 = clock::now();

        // Update Jacobian timer
        m_jac_time +=
            std::chrono::duration_cast<time_unit>(tic1 - tic0).count() * 1e-6;

        // Sparse matrix check
        if(m_matrix_checker(m_J, size))
        {
            std::cout << "Error in Jacobian matrix. This error is fatal.\n";
            exit(13);
        }

        // Clear previous Jacobian matrix
        J.A.clear();
        J.ia.clear();
        J.ja.clear();

        std::size_t nzmax = m_M.A.size() + m_J.A.size();

        // If new size is greater than the current capacity,
        // a reallocation happens
        J.A.reserve(nzmax);
        J.ia.reserve(size + 1);
        J.ja.reserve(nzmax);

        // Replaces deprecated mkl_dcsradd()
        // J: = m_J - M*alpha
        m_matrix_add(-alpha, m_M, m_J, J);
    }
}

}  // namespace daecpp_namespace_name
