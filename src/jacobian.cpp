/*
 * TODO: Description
 */

#include <cmath>
#include "jacobian.h"

namespace daecpp_namespace_name
{

/*
 * Numerical Jacobian
 * Calls rhs N times, hence O(N^2) operations
 */
void Jacobian::operator()(sparse_matrix_holder &J, state_type &x, const double t)
{
    const int size = (int)(x.size());

    state_type f0(size);
    state_type f1(size);

    m_rhs(x, f0, t);

    int ci = 0, cj = 0;

    for(int j = 0; j < size; j++)
    {
        x[j] += m_tol;

        m_rhs(x, f1, t);

        bool first = true;  // first element in a row

        for(int i = 0; i < size; i++)  // loop over columns
        {
            double der = (f1[i] - f0[i]) / m_tol;

            if(std::abs(der) < m_tol)  // skip zero element
            {
                continue;
            }
            else
            {
                J.A[cj]  = der;    // write derivative
                J.ja[cj] = i + 1;  // write column number -- FORTRAN style here
                cj++;

                if(first)
                {
                    J.ia[ci] = cj;  // write ID of the first element in a row --
                                  // FORTRAN style here
                    first = false;
                    ci++;
                }
            }
        }

        x[j] -= m_tol;
    }

    J.ia[size] = cj + 1;  // FORTRAN style here
}

}  // namespace daecpp_namespace_name
