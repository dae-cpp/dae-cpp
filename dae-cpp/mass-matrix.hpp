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

#include "typedefs.hpp"

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
    virtual void operator()(sparse_matrix &M) const = 0;
};

/*
 * Helper child class to create identity Mass matrix of size N
 */
class MassMatrixIdentity : public MassMatrix
{
    const int_type _N{0}; // Indentity matrix size

public:
    MassMatrixIdentity(const int_type N) : MassMatrix(), _N(N) {}

    // TODO: Redo
    void operator()(daecpp::sparse_matrix &M) const
    {
        M.A.resize(_N, 1);
        M.j.resize(_N);
        M.i.resize(_N + 1);

        for (int_type i = 0; i < _N; ++i)
        {
            M.i[i] = i;
            M.j[i] = i;
        }

        M.i[_N] = _N;
    }
};

// TODO: Zero matrix

} // namespace daecpp_namespace_name

#endif // DAECPP_MASS_MATRIX_H
