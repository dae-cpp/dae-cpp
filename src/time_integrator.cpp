/*
 * TODO: Description
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

void TimeIntegrator::operator()(sparse_matrix_holder &Jt, state_type &b, sparse_matrix_holder &J,
                                state_type &x, const state_type x_prev[],
                                const double t, const double dt)
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

    // Calculate full RHS
    // ==================
    // b := -M * dxdt + b

    m_rhs(x, b, t);

    // TODO: check output
    // double precision
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, m_csrA, m_descrA,
                    dxdt.data(), 1.0, b.data());

    // Calculate full Jacobian
    // =======================
    // Jt: = J - M*invdt

    J.A.clear();
    J.ia.clear();
    J.ja.clear();

    m_jac(J, x, t);

    // Sparse matrix check
    if(matrix_checker(J, size))
        exit(1);

    int request = 0;
    int sort    = 0;
    int nzmax   = m_M.A.size() + J.A.size();
    int info;

    double beta = -invdt;

    if(Jt.A.size() != nzmax)
        Jt.A.resize(nzmax);
    if(Jt.ia.size() != size + 1)
        Jt.ia.resize(size + 1);
    if(Jt.ja.size() != nzmax)
        Jt.ja.resize(nzmax);

    // https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-46768951-3369-4425-AD16-643C0E445373.htm

    /*void*/ mkl_dcsradd("N", &request, &sort, &size, &size, J.A.data(),
                         J.ja.data(), J.ia.data(), &beta, m_M.A.data(),
                         m_M.ja.data(), m_M.ia.data(), Jt.A.data(),
                         Jt.ja.data(), Jt.ia.data(), &nzmax, &info);  // double

    // OUTPUT: new Jt and new x
}

}  // namespace daecpp_namespace_name
