/*
 * Converts a matrix from simple three-array format (A, i, j), where A
 * is the non-zero element, i is the column index of A, j is the row index of A,
 * to Intel MKL sparse three-array format described here:
 * https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

/*
 * Input: matrix holder M with simple three-array format
 * Output: matrix holder M with Intel MKL three-array format
 */
void matrix_converter(daecpp::sparse_matrix_holder &M)
{
    // Check if we need to convert the matrix at all.
    // The matrix will be converted if the size of all three arrays (A, ia, ja)
    // is the same. That means each non-zero element has its column and row
    // index.
    const std::size_t size = M.A.size();
    if((size != M.ia.size()) || (size != M.ja.size()))
        return;  // The sizes are different, do nothing

    // The sizes are equal. Do conversion. Create a temporary matrix holder:
    daecpp::sparse_matrix_holder O;

    // We only need to update ia indexes. The first index is 0:
    O.ia.push_back(0);

    // Looking for the other ia indexes
    double row = 0;  // current row
    for(std::size_t i = 0; i < size; i++)
    {
        if(M.ia[i] != row)
        {
            row++;
            O.ia.push_back(i);
        }
    }

    // Finalise ia
    O.ia.push_back(size);

    // Copy ia index to M
    M.ia.clear();
    M.ia = O.ia;
}

void TimeIntegrator::m_matrix_converter(daecpp::sparse_matrix_holder &M)
{
    matrix_converter(M);
}

void MassMatrix::m_matrix_converter(daecpp::sparse_matrix_holder &M)
{
    matrix_converter(M);
}

void Jacobian::m_matrix_converter(daecpp::sparse_matrix_holder &M)
{
    matrix_converter(M);
}

}  // namespace daecpp_namespace_name
