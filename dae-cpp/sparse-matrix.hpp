/*
 * Defines sparse matrix holder in 3-array format.
 *
 * There are two reasons behind this class:
 *  - We can easily switch the linear solver, if needed, by overriding the matrix converter
 *    without changing the main RHS/Mass/Jacobian classes. So even if the linear solver under
 *    the hood changes, the user-defined matrices will still be compatible with the DAE solver.
 *  - The class provides a set of helper tools for the user to define sparse matrices
 *    in a very straightforward way, e.g.:
 *
 *        sparse_matrix M;            // Creates empty sparse matrix M
 *        M.add_element(1, 2, 42.0);  // Fills row 1, column 2 with value 42.0
 *
 * Cons:
 *  - Overhead: the class has to convert matrices to (currently) Eigen format.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#ifndef DAECPP_SPARSE_MATRIX_H
#define DAECPP_SPARSE_MATRIX_H

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Defines sparse matrix holder in 3-array format
 */
struct sparse_matrix
{
    core::rvec A; // Non-zero element A_{ij}
    core::ivec i; // Row index (i) of the element A_{ij}
    core::ivec j; // Column index (j) of the element A_{ij}

    /*
     * Adds non-zero element.
     * Duplicated elements will be summed up.
     *
     * Parameters:
     *     ind_i - row index of the element (int_type)
     *     ind_j - column index of the element (int_type)
     *     A_ij - non-zero element (float_type)
     */
    inline void operator()(const int_type ind_i, const int_type ind_j, const float_type A_ij)
    {
        A.emplace_back(A_ij);
        i.emplace_back(ind_i);
        j.emplace_back(ind_j);
    }

    /*
     * Adds non-zero element.
     * Duplicated elements will be summed up.
     *
     * Parameters:
     *     ind_i - row index of the element (int_type)
     *     ind_j - column index of the element (int_type)
     *     A_ij - non-zero element (float_type)
     */
    inline void add_element(const int_type ind_i, const int_type ind_j, const float_type A_ij)
    {
        operator()(ind_i, ind_j, A_ij);
    }

    /*
     * Reserves memory for `N_elements` non-zero elements
     */
    inline void reserve(const int_type N_elements)
    {
        A.reserve(N_elements);
        i.reserve(N_elements);
        j.reserve(N_elements);
    }

    /*
     * Clears the matrix
     */
    inline void clear()
    {
        A.clear();
        i.clear();
        j.clear();
    }

    /*
     * Performs basic checks of the matrix structure.
     * If it fails, exits with error code (it is not possible to recover).
     */
    inline void check() const noexcept
    {
        constexpr char msg[] = "Three-array sparse matrix check failed. Inconsistent array size.";
        ASSERT(A.size() == i.size(), msg);
        ASSERT(A.size() == j.size(), msg);
    }

    /*
     * Represents the matrix in dense format.
     * Suitable for printing using std::cout.
     *
     * Parameter:
     *     N - matrix size (square matrix N x N) (int_type)
     */
    inline auto dense(const int_type N) const
    {
        return convert(N).toDense();
    }

    /*
     * Converts matrix from 3-array format to Eigen::SparseMatrix format.
     *
     * Parameter:
     *     N - matrix size (square matrix N x N) (int_type)
     *
     * Returns:
     *     Eigen::SparseMatrix<float_type> (core::eimat) sparse matrix
     */
    inline core::eimat convert(const int_type N) const
    {
        ASSERT(N > 0, "daecpp::sparse_matrix.convert(const int_type N): N must be positive.");

        const std::size_t N_el = N_elements(); // Number of non-zero elements

        Eigen::VectorXi nonzeros_per_col = Eigen::VectorXi::Constant(N, 0); // Number of non-zero elements per column

        // Counts non-zero elements per column and checks the indices
        for (std::size_t k = 0; k < N_el; ++k)
        {
            ASSERT(i[k] < N, "Sparse matrix index i is out of boundaries: i = " << i[k] << ", matrix size is " << N << "x" << N);
            ASSERT(j[k] < N, "Sparse matrix index j is out of boundaries: j = " << j[k] << ", matrix size is " << N << "x" << N);
            nonzeros_per_col[j[k]]++;
        }

        core::eimat M(N, N); // Eigen::SparseMatrix<float_type> sparse matrix

        M.reserve(nonzeros_per_col);

        for (std::size_t k = 0; k < N_el; ++k)
        {
            M.coeffRef(i[k], j[k]) += A[k];
        }

        // M.makeCompressed(); // It is already compressed

        return M;
    }

    /*
     * Returns the number of non-zero elements in the matrix
     */
    inline std::size_t N_elements() const noexcept
    {
        return A.size();
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SPARSE_MATRIX_H
