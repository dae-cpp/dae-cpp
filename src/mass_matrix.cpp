/*
 * Helper function to create identity Mass matrix
 */

#include "mass_matrix.h"

namespace daecpp_namespace_name
{

void MassMatrixIdentity::operator()(daecpp::sparse_matrix_holder &M)
{
    for(MKL_INT i = 0; i < m_N; i++)
    {
        M.A.push_back(1.0);
        M.ja.push_back(i);
        M.ia.push_back(i);
    }

    M.ia.push_back(m_N);
}

}  // namespace daecpp_namespace_name
