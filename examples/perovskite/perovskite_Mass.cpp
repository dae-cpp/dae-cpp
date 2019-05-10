/*
 * Mass matrix in three array sparse format:
 * https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format
 */

#include "perovskite_Mass.h"

void MyMassMatrix::operator()(daecpp::sparse_matrix_holder &M)
{
    const MKL_INT N = m_N;

    for(MKL_INT i = 0; i < 2 * N; i++)
    {
        if(i < N)
        {
            M.A.push_back(1.0);
        }
        else
        {
            M.A.push_back(0.0);
        }
        M.ja.push_back(i);
        M.ia.push_back(i);
    }

    M.ia.push_back(2 * N);
}
