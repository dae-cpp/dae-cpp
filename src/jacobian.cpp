/*
 * TODO: Description
 */

#include <cmath>
#include <thread>
#include <mkl_spblas.h>

#include "jacobian.h"

namespace daecpp_namespace_name
{

/*
 * Numerical Jacobian. Central difference scheme. Parallel version.
 * Calls rhs 2*N times, hence O(N^2) operations.
 *//*
void Jacobian::operator()(sparse_matrix_holder &J, const state_type &x,
                          const double t)
{
    const int size = (int)(x.size());

    // Get max number of threads
    unsigned int nth = std::thread::hardware_concurrency();

    // Transposed Jacobian holder
    sparse_matrix_holder Jt;
    Jt.A.resize(J.A.capacity());
    Jt.ia.resize(size + 1);
    Jt.ja.resize(J.ja.capacity());

    int work_thread = 0;

#pragma omp parallel for num_threads(nth) schedule(static)
    for(int th = 0; th < nth; th++)
    {
        int n = size;
        int start_th, end_th;

        start_th = (n * th) / nth;
        end_th   = (n * (th + 1)) / nth;

        state_type f0(n);
        state_type f1(n);
        state_type x1(n);

        x1 = x;

        double jacd;
        double tol  = m_tol;
        double tol2 = 2.0 * tol;

        sparse_matrix_holder jac;

        jac.A.reserve(J.A.capacity());
        jac.ia.reserve(size + 1);
        jac.ja.reserve(J.ja.capacity());

        int ci = 0;

        for(int j = start_th; j < end_th; j++)
        {
            x1[j] -= tol;

            m_rhs(x1, f0, t);

            x1[j] += tol2;

            m_rhs(x1, f1, t);

            bool is_first = true;

            for(int i = 0; i < n; i++)
            {
                jacd = (f1[i] - f0[i]) / tol2;

                if(std::abs(jacd) < tol)
                    continue;

                jac.A.push_back(jacd);
                jac.ja.push_back(i + 1);

                ci++;

                if(is_first)
                {
                    jac.ia.push_back(ci);
                    is_first = false;
                }
            }

            x1[j] -= tol;
        }

        while(true)
        {
            if(th == work_thread)
            {
                Jt.A.insert(Jt.A.end(), jac.A.begin(), jac.A.end());
                Jt.ia.insert(Jt.ia.end(), jac.ia.begin(), jac.ia.end());
                Jt.ja.insert(Jt.ja.end(), jac.ja.begin(), jac.ja.end());

                work_thread++;
            }
            if(work_thread == nth)
                break;
        }
    }

    // Transpose the matrix using mkl_dcsradd
    // https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-46768951-3369-4425-AD16-643C0E445373.htm

    sparse_matrix_holder A;
    for(int i = 0; i < size; i++)
    {
        A.A.push_back(0.0);
        A.ja.push_back(i + 1);
        A.ia.push_back(i + 1);
    }
    A.ia.push_back(size + 1);

    int request = 0;
    int sort    = 0;
    int nzmax   = Jt.A.size();
    int info;

    double beta = 1.0;

    J.A.resize(nzmax);
    J.ia.resize(size + 1);
    J.ja.resize(nzmax);

    mkl_dcsradd("T", &request, &sort, &size, &size, A.A.data(), A.ja.data(),
                A.ia.data(), &beta, Jt.A.data(), Jt.ja.data(), Jt.ia.data(),
                J.A.data(), J.ja.data(), J.ia.data(), &nzmax,
                &info);  // double precision
}
*/


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
    int sort    = 0;
    int nzmax   = J.A.size();
    int info;

    double beta = 1.0;

    sparse_matrix_holder Jt;
    Jt.A.resize(nzmax);
    Jt.ia.resize(size + 1);
    Jt.ja.resize(nzmax);

    // https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-46768951-3369-4425-AD16-643C0E445373.htm

    mkl_dcsradd("T", &request, &sort, &size, &size, A.A.data(), A.ja.data(),
                A.ia.data(), &beta, J.A.data(), J.ja.data(), J.ia.data(),
                Jt.A.data(), Jt.ja.data(), Jt.ia.data(), &nzmax,
                &info);  // double

    J.A  = Jt.A;
    J.ia = Jt.ia;
    J.ja = Jt.ja;
}

}  // namespace daecpp_namespace_name
