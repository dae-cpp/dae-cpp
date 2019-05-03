/*
 * TODO: Description
 */

#include <cmath>      // fabs()
#include <algorithm>  // for std::copy

#if defined(_OPENMP)
#include <omp.h>      // to catch omp_get_max_threads() and OpenMP locks
#endif

#include "jacobian.h"

// clang-format off
#define JACOBIAN_SCHEME 0  // 0 - Central differences: (f(x+h) - f(x-h)) / (2*h)
                           // 1 - Faster but less accurate scheme: (f(x+h) - f(x)) / h
// clang-format on

namespace daecpp_namespace_name
{

/*
 * Numerical Jacobian. Parallel version.
 * Calls RHS up to 2*N times, hence O(N^2) operations.
 */
void Jacobian::operator()(sparse_matrix_holder &J, const state_type &x,
                          const double t)
{
    const int size = (int)(x.size());

    // Get max number of threads.
    // This can be defined using "export OMP_NUM_THREADS=N",
    // where N is the number of threads
#if defined(_OPENMP)
    const int nth = omp_get_max_threads();
#else
    const int nth = 1;
#endif

    state_type_matrix J_values(size, state_type(0));

    std::vector<std::vector<int>> J_rows(size, std::vector<int>(0));

    std::vector<std::vector<int>> sizes_local(size, std::vector<int>(nth));
    std::vector<std::vector<int>> shift_local(size, std::vector<int>(nth));

#if defined(_OPENMP)
    int th_barrier1 = 0;
    int th_barrier2 = 0;

    omp_lock_t writelock1, writelock2;

    omp_init_lock(&writelock1);
    omp_init_lock(&writelock2);
#endif

#pragma omp parallel for num_threads(nth) schedule(static, 1)
    for(int th = 0; th < nth; th++)
    {
        int start_th, end_th;

        start_th = (size * th) / nth;
        end_th   = (size * (th + 1)) / nth;

        state_type f0(size);
        state_type f1(size);
        state_type x1(size);

        x1 = x;

        std::vector<std::vector<float_type>> values_local(
            size, std::vector<float_type>(0));
        std::vector<std::vector<int>> rows_local(size, std::vector<int>(0));

#if JACOBIAN_SCHEME == 1
        m_rhs(x1, f0, t);
#endif

        for(int j = start_th; j < end_th; j++)
        {

#if JACOBIAN_SCHEME == 0
            x1[j] -= m_tol;
            m_rhs(x1, f0, t);
            x1[j] += 2.0 * m_tol;
#else
            x1[j] += m_tol;
#endif

            m_rhs(x1, f1, t);

            for(int i = 0; i < size; i++)
            {
                double jacd;

#if JACOBIAN_SCHEME == 0
                jacd = (f1[i] - f0[i]) / (2.0 * m_tol);
#else
                jacd = (f1[i] - f0[i]) / m_tol;
#endif

                if(std::fabs(jacd) < m_tol)
                    continue;

                values_local[i].push_back(jacd);
                rows_local[i].push_back(j + 1);
                sizes_local[i][th]++;
            }

            x1[j] -= m_tol;
        }

        // Synchronise the threads and assemble local arrays into global
        // Jacobian matrix

#if defined(_OPENMP)
#pragma omp atomic
        th_barrier1++;
        while(true)  // Wait for all threads
        {
            omp_set_lock(&writelock1);
            if(th_barrier1 == nth)
            {
                omp_unset_lock(&writelock1);
                break;
            }
            else
            {
                omp_unset_lock(&writelock1);
            }
        }
#endif

        for(int i = 0; i < size; i++)
        {
            if(th == 0)
            {
                int size_total = 0;
                for(int k = 0; k < nth; k++)
                {
                    size_total += sizes_local[i][k];
                }
                J_values[i].resize(size_total);
                J_rows[i].resize(size_total);
            }
            else
            {
                for(int k = 0; k < th; k++)
                {
                    shift_local[i][th] += sizes_local[i][k];
                }
            }
        }

#if defined(_OPENMP)
        // Make sure the master thread resized the global arrays
        if(th == 0)
        {
#pragma omp atomic
            th_barrier2++;
        }
        while(true)  // Wait for master thread
        {
            omp_set_lock(&writelock2);
            if(th_barrier2 > 0)
            {
                omp_unset_lock(&writelock2);
                break;
            }
            else
            {
                omp_unset_lock(&writelock2);
            }
        }
#endif

        for(int i = 0; i < size; i++)
        {
            std::copy(values_local[i].begin(), values_local[i].end(),
                      J_values[i].begin() + shift_local[i][th]);
            std::copy(rows_local[i].begin(), rows_local[i].end(),
                      J_rows[i].begin() + shift_local[i][th]);
        }

    }  // Parallel section end

    // Unroll global array to CSR format and transpose Jacobian

    int ci = 1;

    for(int i = 0; i < size; i++)
    {
        J.A.insert(J.A.end(), J_values[i].begin(), J_values[i].end());
        J.ja.insert(J.ja.end(), J_rows[i].begin(), J_rows[i].end());
        J.ia.push_back(ci);

        ci += J_values[i].size();
    }

    J.ia.push_back(ci);
}

}  // namespace daecpp_namespace_name
