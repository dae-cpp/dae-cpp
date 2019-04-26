/*
 * TODO: Description
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

void TimeIntegrator::operator()(sparse_matrix_holder &Jt, state_type &b,
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

    sparse_matrix_holder J;

    J.A.reserve(m_opt.max_size_mult * size);
    J.ia.reserve(size + 1);
    J.ja.reserve(m_opt.max_size_mult * size);

    state_type::size_type sza = J.A.capacity();
    state_type::size_type szj = J.ja.capacity();

    m_jac(J, x, t);

    // MKL sparse matrix check

    if(sza != J.A.capacity() || szj != J.ja.capacity())
    {
        // test this
        double multa = (double)J.A.capacity() / (double)sza;
        double multj = (double)J.ja.capacity() / (double)szj;
        double mult  = multa > multj ? multa : multj;
        std::cout << "\nWarning: Jacobian matrix capacity is bigger than "
                     "expected. This may lead to performance degradation. "
                     "Suggested multiplier max_size_mult should be at least "
                  << mult * m_opt.max_size_mult << '\n';
        m_opt.max_size_mult *= (int)(mult + 1.001);
    }

    int request = 0;
    int sort    = 0;
    int nzmax   = m_M.A.size() + J.A.size();
    int info;

    double beta = -invdt;

    Jt.A.resize(nzmax);
    Jt.ia.resize(size + 1);
    Jt.ja.resize(nzmax);

    // https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-46768951-3369-4425-AD16-643C0E445373.htm

    /*void*/ mkl_dcsradd("N", &request, &sort, &size, &size, J.A.data(),
                         J.ja.data(), J.ia.data(), &beta, m_M.A.data(),
                         m_M.ja.data(), m_M.ia.data(), Jt.A.data(),
                         Jt.ja.data(), Jt.ia.data(), &nzmax, &info);  // double

    // OUTPUT: new Jt and new x
}

}  // namespace daecpp_namespace_name
