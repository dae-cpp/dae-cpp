/*
 * TODO: Description
 */

#include <cmath>
#include <mkl_spblas.h>

#include "jacobian.h"

namespace daecpp_namespace_name
{

/*
 * Numerical Jacobian.
 * Calls rhs N times, hence O(N^2) operations.
 * Probably should hire MKL jacobi routine
 * https://software.intel.com/en-us/mkl-developer-reference-c-jacobian-matrix-calculation-routines
 */
void Jacobian::operator()(sparse_matrix_holder &J, state_type &x,
                          const double t)
{
    const int size = (int)(x.size());

    state_type f0(size);
    state_type f1(size);

    m_rhs(x, f0, t);

    int cj = 0;

    for(int j = 0; j < size; j++)
    {
        x[j] += m_tol;

        m_rhs(x, f1, t);

        bool is_first = true;  // first element in a row

        for(int i = 0; i < size; i++)  // loop over columns
        {
            double der = (f1[i] - f0[i]) / m_tol;

            if(std::abs(der) < m_tol)  // skip zero element
            {
                continue;
            }
            else
            {
                J.A.push_back(der);     // write derivative
                J.ja.push_back(i + 1);  // write column number -- FORTRAN style
                cj++;

                if(is_first)
                {
                    J.ia.push_back(cj);  // write ID of the first element in a
                                         // row -- FORTRAN style here
                    is_first = false;
                }
            }
        }

        x[j] -= m_tol;
    }

    J.ia.push_back(cj + 1);  // FORTRAN style here

    // The code above produces transposed Jacobian.
    // The only way to transpose it I found so far is
    // using mkl_dcsradd: add zero matrix and transpose.
    // This needs to be improved...

    // Create zero matrix A
    sparse_matrix_holder A;
    for(int i = 0; i < size; i++)
    {
        A.A.push_back(0.0);
        A.ja.push_back(i + 1);
        A.ia.push_back(i + 1);
    }
    A.ia.push_back(size + 1);

    // Init mkl_dcsradd and perform matrix addition
    int request = 0;
    int sort = 0;
    int nzmax = J.A.size();
    int info;

    double beta = 1.0;

    sparse_matrix_holder Jt;
    Jt.A.resize(nzmax);
    Jt.ia.resize(size + 1);
    Jt.ja.resize(nzmax);

    // https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-46768951-3369-4425-AD16-643C0E445373.htm

    mkl_dcsradd("T", &request, &sort, &size, &size, A.A.data(),
        A.ja.data(), A.ia.data(), &beta, J.A.data(),
        J.ja.data(), J.ia.data(), Jt.A.data(),
        Jt.ja.data(), Jt.ia.data(), &nzmax, &info);  // double

    J.A = Jt.A;
    J.ia = Jt.ia;
    J.ja = Jt.ja;
}

}  // namespace daecpp_namespace_name
