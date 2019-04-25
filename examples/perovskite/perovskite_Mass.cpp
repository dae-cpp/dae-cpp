/*
 * TODO: Description of the class
 */

#include "perovskite_Mass.h"

void MyMassMatrix::operator()(daecpp::sparse_matrix_holder &M)
{
    const int N = ((int)(M.ia.size()) - 1) / 2;  // check this

    for(int i = 0; i < N; i++)
    {
        M.A[i]  = 1.0;
        M.ja[i] = i + 1;
        M.ia[i] = i + 1;
    }
    M.ia[N] = N + 1;
}
