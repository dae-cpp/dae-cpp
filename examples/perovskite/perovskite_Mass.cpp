/*
 * TODO: Description of the class
 */

#include "perovskite_Mass.h"

#include "mass_matrix.h"

void Mass_Matrix::operator()(state_type &M, vector_type_int &ia, vector_type_int &ja)
{
    const int N = m_N; // number of points

    for (int i = 0; i < N; i++)
    {
        M[i] = 1.0;
        ja[i] = i + 1;
        ia[i] = i + 1;
    }
    ia[N] = N + 1;
}

void Mass_Matrix::operator()(state_type &M)
{
    const int N = m_N; // number of points

    for (int i = 0; i < N; i++)
    {
        M[i] = 1.0;
        M[i+N] = 0.0;
    }
}
