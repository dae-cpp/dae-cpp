/*
 * Performs numerical differentiation of the RHS with the given tolerance to
 * estimate numerical Jacobian matrix
 */
#include <iostream>   // std::cout
#include <iomanip>    // std::setw etc.
#include <cmath>      // std::abs
#include <algorithm>  // std::copy

#if defined(_OPENMP)
#include <omp.h>  // to catch omp_get_max_threads()
#endif

#include "jacobian.h"

// clang-format off
#define JACOBIAN_SCHEME 1  // 0 - Central differences: (f(x+h) - f(x-h)) / (2*h)
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
    const MKL_INT size   = (MKL_INT)(x.size());
    const double  tol    = m_tol;
    const double  tol2   = tol * tol;
    const double  invtol = 1.0 / tol;

    // Get max number of threads.
    // This can be defined using "export OMP_NUM_THREADS=N",
    // where N is the number of threads
#if defined(_OPENMP)
    const int nth = omp_get_max_threads();
#else
    const int nth = 1;
#endif

    state_type_matrix J_values(size, state_type(0));

    std::vector<std::vector<MKL_INT>> J_rows(size, std::vector<MKL_INT>(0));

    std::vector<std::vector<MKL_INT>> sizes_local(size,
                                                  std::vector<MKL_INT>(nth));
    std::vector<std::vector<MKL_INT>> shift_local(size,
                                                  std::vector<MKL_INT>(nth));

    // Parallel section
#pragma omp parallel num_threads(nth)
    {
#if defined(_OPENMP)
        int th = omp_get_thread_num();
#else
        int th = 0;
#endif

        int start_th = (size * th) / nth;
        int end_th   = (size * (th + 1)) / nth;

        state_type f0(size);
        state_type f1(size);
        state_type x1(size);

        x1 = x;

        std::vector<std::vector<float_type>> values_local(
            size, std::vector<float_type>(0));
        std::vector<std::vector<MKL_INT>> rows_local(size,
            std::vector<MKL_INT>(0));

#if JACOBIAN_SCHEME == 1
        m_rhs(x1, f0, t);
#endif

        for(MKL_INT j = start_th; j < end_th; j++)
        {
            float_type x1_backup = x1[j];

#if JACOBIAN_SCHEME == 0
            x1[j] -= tol;
            m_rhs(x1, f0, t);
            x1[j] = x1_backup + tol;
#else
            x1[j] += tol;
#endif

            m_rhs(x1, f1, t);

            for(MKL_INT i = 0; i < size; i++)
            {
                double diff = f1[i] - f0[i];

                if(std::abs(diff) < tol2)
                    continue;

#if JACOBIAN_SCHEME == 0
                double jacd = diff * invtol * 0.5;
#else
                double jacd = diff * invtol;
#endif

                values_local[i].push_back((float_type)jacd);
#ifdef DAE_FORTRAN_STYLE
                rows_local[i].push_back(j + 1);
#else
                rows_local[i].push_back(j);
#endif
                sizes_local[i][th]++;
            }

            x1[j] = x1_backup;
        }

#if defined(_OPENMP)
#pragma omp barrier  // Synchronise the threads
#endif

        // Resize the global arrays
        for(MKL_INT i = 0; i < size; i++)
        {
            if(th == 0)
            {
                MKL_INT size_total = 0;
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
#pragma omp barrier  // Make sure the master thread resized the global arrays
#endif

        // Assemble local arrays into global Jacobian matrix
        for(MKL_INT i = 0; i < size; i++)
        {
            std::copy(values_local[i].begin(), values_local[i].end(),
                      J_values[i].begin() + shift_local[i][th]);
            std::copy(rows_local[i].begin(), rows_local[i].end(),
                      J_rows[i].begin() + shift_local[i][th]);
        }

    }  // Parallel section end

#ifdef DAE_FORTRAN_STYLE
    MKL_INT ci = 1;
#else
    MKL_INT ci = 0;
#endif

    // Unroll global array to CSR format and transpose Jacobian
    for(MKL_INT i = 0; i < size; i++)
    {
        J.A.insert(J.A.end(), J_values[i].begin(), J_values[i].end());
        J.ja.insert(J.ja.end(), J_rows[i].begin(), J_rows[i].end());
        J.ia.push_back(ci);

        ci += J_values[i].size();
    }

    J.ia.push_back(ci);
}

/*
 * Helper function to show Jacobian structure
 */
void Jacobian::print(const state_type &x, const double t)
{
    if(x.size() > 1000)
    {
        std::cout << "\nJacobian::print -- too much output. Skipped.\n";
        return;
    }

    sparse_matrix_holder J;
    this->operator()(J, x, t);

    std::cout << std::right;
    std::cout << "\nJacobian matrix at time t = " << t << ':';
    std::cout << "\n-----------------------------------------\n";
    std::cout << std::setw(7) << "i" << std::setw(16) << "J.A |"
              << std::setw(10) << "J.ja |" << std::setw(8) << "J.ia";
    std::cout << "\n-----------------------------------------\n";

    for(std::size_t i = 0; i < J.A.size(); i++)
    {
        std::cout << std::setw(7) << i << ": " << std::setw(12) << J.A[i]
                  << " | ";

        if(i < J.ja.size())
            std::cout << std::setw(7) << J.ja[i] << " | ";
        else
            std::cout << std::setw(7) << "???" << " | ";  // Error in Jacobian
                                                          // matrix structure
        if(i < J.ia.size())
            std::cout << std::setw(7) << J.ia[i];
        std::cout << std::endl;
    }
}

}  // namespace daecpp_namespace_name
