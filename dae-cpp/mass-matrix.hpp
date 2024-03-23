/*
 * Mass matrix classes.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_MASS_MATRIX_H
#define DAECPP_MASS_MATRIX_H

#include "sparse-matrix.hpp"

namespace daecpp_namespace_name
{

/*
 * Parent Mass Matrix class. This class is abstract and must be inherited.
 */
class MassMatrix
{
    // TODO: Matrix converter
    // /*
    //  * Sparse matrix converter from simple three-array format to Intel MKL
    //  * three array format.
    //  * Input: matrix holder M with simple three-array format
    //  * Output: matrix holder M with Intel MKL three-array format
    //  */
    // void m_matrix_converter(daecpp::sparse_matrix_holder &M);

public:
    /*
     * The matrix should be defined in sparse format,
     * see three array sparse format decription on
     * https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format
     *
     * The Mass matrix is static, i.e. it does not depend on time t and
     * the vector x
     *
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(sparse_matrix &M, const double t) const = 0;
};

/*
 * Helper child class to create identity Mass matrix of size N
 */
class MassMatrixIdentity : public MassMatrix
{
    const int_type _N{0}; // Indentity matrix size

public:
    MassMatrixIdentity(const int_type N) : MassMatrix(), _N(N) {}

    // TODO: Can we setup identity matrix in Eigen straight away?
    void operator()(sparse_matrix &M, const double t) const
    {
        M.A.resize(_N, 1.0);
        M.i.resize(_N); // Resize and then overwrite in a loop worked faster than reserve and push_back
        M.j.resize(_N);

        for (int_type i = 0; i < _N; ++i)
        {
            M.i[i] = i;
            M.j[i] = i;
        }
    }
};

// TODO: Zero matrix
class MassMatrixZero : public MassMatrix
{
    // const int_type _N{0}; // Zero matrix size

public:
    MassMatrixZero() : MassMatrix() {}

    // TODO: Can we setup zero matrix in Eigen straight away?
    void operator()(sparse_matrix &M, const double t) const
    {
        M.A.clear();
        M.i.clear();
        M.j.clear();
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_MASS_MATRIX_H
