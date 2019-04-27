/*
 * TODO: Description
 */

#include <omp.h>
#include <cmath>
#include <mkl_spblas.h>

#include "jacobian.h"

// clang-format off
#define JACOBIAN_SCHEME 1  // 0 - Central differences: (f(x+h) - f(x-h)) / (2*h)
                           // 1 - Almost twice faster scheme: (f(x+h) - f(x)) / h
// clang-format on

namespace daecpp_namespace_name
{

/*
 * Numerical Jacobian. Parallel version.
 * Calls rhs 2*N times, hence O(N^2) operations.
 */
void Jacobian::operator()(sparse_matrix_holder &J, const state_type &x,
                          const double t)
{
    const int size = (int)(x.size());

    // Get max number of threads.
    // This can be set from the terminal using export OMP_NUM_THREADS=N,
    // where N is the number of threads
    int nth = omp_get_max_threads();

    // Transposed Jacobian holder
    sparse_matrix_holder Jt;
    Jt.A.reserve(J.A.capacity());
    Jt.ia.reserve(size + 1);
    Jt.ja.reserve(J.ja.capacity());

    int work_thread    = 0;
    int work_thread_ia = 0;

#pragma omp parallel for num_threads(nth) schedule(static, 1)
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

#if JACOBIAN_SCHEME == 1
        m_rhs(x1, f0, t);
#endif

        int ci = 0;

        for(int j = start_th; j < end_th; j++)
        {
#if JACOBIAN_SCHEME == 0
            x1[j] -= tol;
            m_rhs(x1, f0, t);
            x1[j] += tol2;
#else
            x1[j] += tol;
#endif
            m_rhs(x1, f1, t);

            bool is_first = true;

            for(int i = 0; i < n; i++)
            {
#if JACOBIAN_SCHEME == 0
                jacd = (f1[i] - f0[i]) / tol2;
#else
                jacd = (f1[i] - f0[i]) / tol;
#endif
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

        int ia_shift;

        while(true)
        {
#pragma omp critical
            if(th == work_thread)
            {
                ia_shift = Jt.A.size();

                Jt.A.insert(Jt.A.end(), jac.A.begin(), jac.A.end());
                Jt.ja.insert(Jt.ja.end(), jac.ja.begin(), jac.ja.end());

                work_thread++;
            }
            if(th < work_thread)
                break;
        }

        for(vector_type_int::iterator it = jac.ia.begin(); it != jac.ia.end();
            ++it)
        {
            *it += ia_shift;
        }

        while(true)
        {
#pragma omp critical
            if(th == work_thread_ia)
            {
                Jt.ia.insert(Jt.ia.end(), jac.ia.begin(), jac.ia.end());

                work_thread_ia++;
            }
            if(th < work_thread_ia)
                break;
        }
    }  // Parallel section end

    Jt.ia.push_back(Jt.A.size() + 1);

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

}  // namespace daecpp_namespace_name
