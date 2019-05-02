/*
 * TODO: Description
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

void TimeIntegrator::operator()(sparse_matrix_holder &Jt, state_type &b,
                                sparse_matrix_holder &J, state_type &x,
                                const state_type x_prev[], const double t,
                                const double dt)
{
    const int size = (int)(x.size());

    const double invdt = 1.0 / dt;

    state_type dxdt(size);

    if(m_scheme == 1)  // first order BDF
    {
        for(int i = 0; i < size; i++)
        {
            dxdt[i] = (x[i] - x_prev[0][i]) * invdt;
        }
    }
    else
    {
        // internal error
    }

    // Calculate RHS
    m_rhs(x, b, t);

    // TODO: check output
    // b := -M * dxdt + b
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, m_csrA, m_descrA,
                    dxdt.data(), 1.0, b.data());

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
    matrix_add(-invdt, m_M, J, Jt);
}

}  // namespace daecpp_namespace_name
