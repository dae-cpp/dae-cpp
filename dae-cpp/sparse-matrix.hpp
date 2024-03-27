/*
 * Defines sparse matrix holder in 3-array format.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
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
    fvec A; // Non-zero element A_{ij}
    ivec i; // Row index (i) of the element A_{ij}
    ivec j; // Column index (j) of the element A_{ij}

    /*
     * Adds non-zero element.
     * Duplicated elements will be summed up.
     *
     * Parameters:
     *     A_ij - non-zero element (float_type)
     *     ind_i - row index of the element (int_type)
     *     ind_j - column index of the element (int_type)
     */
    inline void operator()(const float_type A_ij, const int_type ind_i, const int_type ind_j)
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
     *     A_ij - non-zero element (float_type)
     *     ind_i - row index of the element (int_type)
     *     ind_j - column index of the element (int_type)
     */
    inline void add_element(const float_type A_ij, const int_type ind_i, const int_type ind_j)
    {
        operator()(A_ij, ind_i, ind_j);
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
        const char *msg = "Three-array sparse matrix check failed. Inconsistent array size.";
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
     *     Eigen::SparseMatrix<float_type> sparse matrix
     */
    Eigen::SparseMatrix<float_type> convert(const int_type N) const
    {
        typedef Eigen::Triplet<float_type> T;

        std::vector<T> triplet;

        const std::size_t N_el = N_elements(); // Number of non-zero elements

        triplet.reserve(N_el);

        for (std::size_t k = 0; k < N_el; ++k)
        {
            ASSERT(i[k] < N, "Sparse matrix index i is out of boundaries: i = " << i[k] << ", matrix size is " << N << "x" << N);
            ASSERT(j[k] < N, "Sparse matrix index j is out of boundaries: j = " << j[k] << ", matrix size is " << N << "x" << N);
            triplet.emplace_back(T(i[k], j[k], A[k]));
        }

        Eigen::SparseMatrix<float_type> M(N, N);

        M.setFromTriplets(triplet.begin(), triplet.end());

        return M; // TODO: Check if Eigen::SparseMatrix has move-semantics
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
