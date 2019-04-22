/*
 * TODO: Description of the class
 */

#include <cmath>
#include "jacobian.h"

namespace daecpp_namespace_name
{

// make sure m_tol is defined!

/*
void Jacobian::operator() (const state_type &x, const state_type &x_prev,
state_type &J, const double t)
{
    const int size = (int)(x.size());

    state_type f0(size);

    m_rhs(x, f0, t);

#pragma omp parallel for num_threads(6)
    for(int j = 0; j < size; j++)
    {
        state_type f1(size);
        state_type x1(size);

        for(int i = 0; i < size; i++)
        {
            x1[i] = x[i];
        }
        x1[j] += m_tol;

        m_rhs(x1, f1, t);

        for(int i = 0; i < size; i++)
        {
            //J(i, j) = (f1[i] - f0[i])/m_tol;
            J[j*size + i] = (f1[i] - f0[i])/m_tol;
        }
    }
}
*/

// TODO: it is not general at the moment!!!
void Jacobian::operator()(state_type &J, vector_type_int &ia,
                          vector_type_int &ja, state_type &x, const double t)
{
    const int size   = (int)(x.size());
    const int size05 = size / 2;

    state_type f0(size);
    state_type f1(size);
    int        k[9];

    m_rhs(x, f0, t);

    int c = 0, c1 = 0;
    for(int i = 0; i < size; i++)
    {
        x[i] += m_tol;

        m_rhs(x, f1, t);

        int first = 1;

        k[0] = i - size05 - 1;
        k[1] = i - size05;
        k[2] = i - size05 + 1;
        k[3] = i - 1;
        k[4] = i;
        k[5] = i + 1;
        k[6] = i + size05 - 1;
        k[7] = i + size05;
        k[8] = i + size05 + 1;

        for(int j = 0; j < 9; j++)
        {
            if(k[j] < 0 || k[j] >= size)
                continue;

            double val = (f1[k[j]] - f0[k[j]]) / m_tol;  // a[i * 2 * N + j];
            if(std::abs(val) < 1.0e-6)
                continue;
            else
            {
                J[c]  = val;
                ja[c] = k[j] + 1;
                c++;
                if(first)
                {
                    first  = 0;
                    ia[c1] = c;
                    c1++;
                }
            }
        }

        x[i] -= m_tol;
    }
    ia[size] = c + 1;
}

}  // namespace daecpp_namespace_name
