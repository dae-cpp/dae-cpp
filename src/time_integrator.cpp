/*
 * TODO: Description
 */

#include <mkl_types.h>
#include <mkl_spblas.h>
#include "time_integrator.h"

namespace daecpp_namespace_name
{

// t may be const
void TimeIntegrator::operator()(sparse_matrix_holder &Jt, state_type &b, state_type &x,
                                state_type x_prev[], double &t, const double dt)
{
    const int size = (int)(x.size());

    const double invdt = 1.0 / dt;

    state_type dxdt(size);

    if(m_scheme == 1)  // first order BDF
    {
        for(int i = 0; i < size; i++)
        {
            dxdt[i] = (x[i] - x_prev[i][0]) * invdt;
        }
    }
    else
    {
        // internal error
    }

    m_rhs(x, b, t);

    // vector b = sparse matrix M * vector dxdt
    // vector x -= vector b

    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;

    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;

    // Structure with sparse matrix stored in CSR format
    sparse_matrix_t csrA;

    // check output
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, (MKL_INT)size,
                            (MKL_INT)size, m_M.ia.data(), m_M.ia.data() + 1,
                            m_M.ja.data(),
                            m_M.A.data());  // double precision, BASE_ZERO

    // Analyze sparse matrix; choose proper kernels and workload balancing
    // strategy
    // mkl_sparse_optimize ( csrA );

    // check output
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csrA, descrA,
                    dxdt.data(), 1.0, b.data());  // double precision
    // x := -M * dxdt + x

    mkl_sparse_destroy(csrA);

    // sparse matrix Mdt = sparse matrix M * invdt element-wise
    // sparse matrix J = sparse matrix J - sparse matrix Mdt
    // J: = J - M*invdt

    // https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-46768951-3369-4425-AD16-643C0E445373.htm

    sparse_matrix_holder J;

    // assert m_opt.N > 0
    J.A.reserve(m_opt.max_size_mult * size);
    J.ia.reserve(size + 1);
    J.ja.reserve(m_opt.max_size_mult * size);

    state_type::size_type sz = J.A.capacity();

    m_jac(J, x, t);

    // MKL sparse matrix check

    if(sz != J.A.capacity())
    {
        // test this
        double mult = (double)J.A.capacity() / (double)sz;
        std::cout << "\nWarning: Jacobian matrix capacity is bigger than "
                     "expected. This may lead to performance degradation. "
                     "Suggested multiplier max_size_mult should be at least "
                  << mult * m_opt.max_size_mult << '\n';
        // increase max_size_mult here!
    }

    int *  sort = nullptr;
    double beta = -invdt;
    int    info;
    int    nzmax = m_M.A.size() + J.A.size();

    Jt.A.resize(nzmax);
    Jt.ia.resize(size + 1);
    Jt.ja.resize(nzmax);

    /*void*/ mkl_dcsradd("N", 0, sort, &size, &size, J.A.data(), J.ja.data(),
                         J.ia.data(), &beta, m_M.A.data(), m_M.ja.data(),
                         m_M.ia.data(), Jt.A.data(), Jt.ja.data(), Jt.ia.data(),
                         &nzmax, &info);  // double

    // OUTPUT: new Jt and new x
}

}  // namespace daecpp_namespace_name
