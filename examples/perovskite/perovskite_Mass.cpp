/*
 * TODO: Description of the class
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
        M.ja.push_back(i + 1);
        M.ia.push_back(i + 1);
    }

    M.ia.push_back(2 * N + 1);
}
