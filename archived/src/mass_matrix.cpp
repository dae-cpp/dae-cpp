/*
 * Helper function to create identity Mass matrix
 */

#include "mass_matrix.h"

namespace daecpp_namespace_name
{

void MassMatrixIdentity::operator()(daecpp::sparse_matrix_holder &M)
{
    M.A.resize(m_N, 1);
    M.ja.resize(m_N);
    M.ia.resize(m_N + 1);

    for(MKL_INT i = 0; i < m_N; i++)
    {
        M.ja[i] = i;
        M.ia[i] = i;
    }

    M.ia[m_N] = m_N;
}

}  // namespace daecpp_namespace_name
