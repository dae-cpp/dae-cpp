/*
 * TODO: Description
 */

#include <mkl_types.h>
#include <mkl_spblas.h>
#include "time_integrator.h"

namespace daecpp_namespace_name
{

void TimeIntegrator::operator()(state_type &Jt, vector_type_int &Jt_ia, vector_type_int &Jt_ja, state_type &x, state_type x_prev[], double &t, const double dt)
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

    // vector b = sparse matrix M * vector dxdt
    // vector x -= vector b

    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;

    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;

    // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       csrA;

    // check output
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, (MKL_INT)size, (MKL_INT)size, m_ia.data(), m_ia.data()+1, m_ja.data(), m_M.data());  // double precision, BASE_ZERO
    
    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    //mkl_sparse_optimize ( csrA );

    // check output
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csrA, descrA, dxdt.data(), 1.0, x.data());  // double precision
    // x := -M * dxdt + x

    mkl_sparse_destroy(csrA);


    // sparse matrix Mdt = sparse matrix M * invdt element-wise
    // sparse matrix J = sparse matrix J - sparse matrix Mdt
    // J: = J - M*invdt

    // https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-46768951-3369-4425-AD16-643C0E445373.htm
    
    state_type J;
    vector_type_int J_ia;
    vector_type_int J_ja;
    
    // assert m_opt.N > 0
    J.reserve(m_opt.max_size_mult * size);
    J_ia.reserve(size + 1);
    J_ja.reserve(m_opt.max_size_mult * size);

    state_type::size_type sz = J.capacity();

    m_jac(J, J_ia, J_ja, x, t);

    // MKL sparse matrix check

    if(sz != J.capacity()) {
        // test this
        double mult = (double)J.capacity() / (double)sz;
        std::cout << "\nWarning: Jacobian matrix capacity is bigger than expected. This may lead to performance degradation. Suggested multiplier max_size_mult should be at least " << mult*m_opt.max_size_mult << '\n';
        // increase max_size_mult here!
    }

    int *sort = nullptr;
    double beta = -invdt;
    int info;
    int nzmax = m_M.size() + J.size();

    Jt.resize(nzmax);
    Jt_ia.resize(size + 1);
    Jt_ja.resize(nzmax);

    /*void*/ mkl_dcsradd('N', 0, sort, &size, &size,
        J.data(), J_ja.data(), J_ia.data(), &beta, m_M.data(), m_ja.data(), m_ia.data(),
        Jt, Jt_ja, Jt_ia, &nzmax, &info);  // double

    // time lapse
    t += dt;
    
    // OUTPUT: new Jt and new x, new t
}

}  // namespace daecpp_namespace_name
