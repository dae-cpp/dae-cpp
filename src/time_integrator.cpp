/*
 * TODO: Description
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

TimeIntegrator::TimeIntegrator(RHS &rhs, Jacobian &jac, MassMatrix &mass,
                               SolverOptions &opt, const MKL_INT size)
    : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
{
    // Reserve memory for at least 1-diagonal mass matrix
    m_M.A.reserve(size);
    m_M.ja.reserve(size);
    m_M.ia.reserve(size + 1);

    // Get static mass matrix
    m_mass(m_M);

    // User defined sparse matrix check
    if(matrix_checker(m_M, size))
    {
        std::cout << "Error in Mass matrix.\n";
        exit(1);
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
                     "Mass matrix.\n";
        exit(1);
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

void TimeIntegrator::operator()(sparse_matrix_holder &Jt, state_type &b,
                                sparse_matrix_holder &J, state_type &x,
                                const state_type_matrix &x_prev, const double t,
                                const double dt)
{
    const MKL_INT size = (MKL_INT)(x.size());

    const double invdt = 1.0 / dt;

    double alpha  = 0;
    int    scheme = m_scheme - 1;

    state_type dxdt(size);

    for(MKL_INT i = 0; i < size; i++)
    {
        double val = 0;
        val += BDF_COEF[scheme][0] * x[i];
        for(int d = 1; d <= m_scheme; d++)
        {
            val += BDF_COEF[scheme][d] * x_prev[d - 1][i];
        }
        val *= invdt / BDF_COEF[scheme][7];
        dxdt[i] = val;
    }

    alpha = invdt * ALPHA_COEF[scheme];

    // Calculate RHS
    m_rhs(x, b, t);

    // TODO: check output
    // b := -M * dxdt + b
#ifdef DAE_SINGLE
    mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, m_csrA, m_descrA,
                    dxdt.data(), 1.0, b.data());
#else
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, m_csrA, m_descrA,
                    dxdt.data(), 1.0, b.data());
#endif

    J.A.clear();
    J.ia.clear();
    J.ja.clear();

    // Calculate Jacobian
    m_jac(J, x, t);

    // Sparse matrix check
    if(matrix_checker(J, size))
    {
        std::cout << "Error in Jacobian matrix.\n";
        exit(1);
    }

    size_t nzmax = m_M.A.size() + J.A.size();

    if(Jt.A.size() != nzmax)
        Jt.A.resize(nzmax);
    if(Jt.ia.size() != (size_t)(size) + 1)
        Jt.ia.resize(size + 1);
    if(Jt.ja.size() != nzmax)
        Jt.ja.resize(nzmax);

    // Replaces deprecated mkl_dcsradd()
    // Jt: = J - M*invdt
    matrix_add(-alpha, m_M, J, Jt);
}

}  // namespace daecpp_namespace_name
