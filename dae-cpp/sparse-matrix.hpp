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

#include <Eigen/Sparse>

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Defines sparse matrix holder in 3-array format.
 */
struct sparse_matrix
{
    fvec A; // Non-zero element A_{ij}
    ivec i; // Row index (i) of the element A_{ij}
    ivec j; // Column index (j) of the element A_{ij}

    // Adds non-zero element. Duplicated elements will be summed up.
    inline void operator()(const float_type A, const int_type i, const int_type j);

    // Adds non-zero element. Duplicated elements will be summed up.
    inline void add_element(const float_type A, const int_type i, const int_type j);

    // Reserves memory for `N_elements` non-zero elements
    inline void reserve(const int_type N_elements);

    // Performs basic checks of the matrix structure
    inline void check() const noexcept;

    // Represents the matrix in dense format. Suitable for printing using std::cout.
    inline auto dense(const int_type N) const;

    // Converts matrix from 3-array format to Eigen::SparseMatrix format
    Eigen::SparseMatrix<float_type> convert(const int_type N) const;

private:
    // Returns the number of non-zero elements in the matrix
    inline std::size_t N_elements() const noexcept;
};

/*
 * Adds non-zero element.
 * Duplicated elements will be summed up.
 *
 * Parameters:
 *     A_ij - non-zero element (float-type)
 *     ind_i - row index of the element (int_type)
 *     ind_j - column index of the element (int_type)
 */
inline void sparse_matrix::operator()(float_type A_ij, int_type ind_i, int_type ind_j)
{
    A.push_back(A_ij);
    i.push_back(ind_i);
    j.push_back(ind_j);
}

/*
 * Adds non-zero element.
 * Duplicated elements will be summed up.
 *
 * Parameters:
 *     A_ij - non-zero element (float-type)
 *     ind_i - row index of the element (int_type)
 *     ind_j - column index of the element (int_type)
 */
inline void sparse_matrix::add_element(const float_type A, const int_type i, const int_type j)
{
    operator()(A, i, j);
}

/*
 * Reserves memory for `N_elements` non-zero elements.
 */
inline void sparse_matrix::reserve(int_type N_elements)
{
    A.reserve(N_elements);
    i.reserve(N_elements);
    j.reserve(N_elements);
}

/*
 * Performs basic checks of the matrix structure.
 * If it fails, exits with error code (it is not possible to recover).
 */
inline void sparse_matrix::check() const noexcept
{
    ASSERT(A.size() == i.size(), "");
    ASSERT(A.size() == j.size(), "");
}

/*
 * Represents the matrix in dense format.
 * Suitable for printing using std::cout.
 *
 * Parameter:
 *     N - matrix size (square matrix N x N) (int_type)
 */
inline auto sparse_matrix::dense(const int_type N) const
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
Eigen::SparseMatrix<float_type> sparse_matrix::convert(const int_type N) const
{
    typedef Eigen::Triplet<float_type> T;

    std::vector<T> triplet;

    const std::size_t N_el = N_elements(); // Number of non-zero elements

    triplet.reserve(N_el);

    for (std::size_t k = 0; k < N_el; ++k)
    {
        ASSERT(i[k] < N, "Sparse matrix index i is out of boundaries: i = " << i[k] << ", matrix size is " << N << "x" << N);
        ASSERT(j[k] < N, "Sparse matrix index j is out of boundaries: j = " << j[k] << ", matrix size is " << N << "x" << N);
        triplet.push_back(T(i[k], j[k], A[k]));
    }

    Eigen::SparseMatrix<float_type> M(N, N);

    M.setFromTriplets(triplet.begin(), triplet.end());

    return M; // TODO: Check if Eigen::SparseMatrix has move-semantics
}

/*
 * Returns the number of non-zero elements in the matrix.
 */
inline std::size_t sparse_matrix::N_elements() const noexcept
{
    return A.size();
}

} // namespace daecpp_namespace_name

#endif // DAECPP_SPARSE_MATRIX_H