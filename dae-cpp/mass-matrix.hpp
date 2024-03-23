/*
 * Mass matrix class.
 * Defines the Mass matrix `M` of the DAE system `M dx/dt = f`.
 * This class is abstract and must be inherited.
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
 * Mass matrix class.
 * This class is abstract and must be inherited.
 */
class MassMatrix
{
public:
    /*
     * Defines the Mass matrix `M` of the DAE system `M dx/dt = f`.
     * The Mass matrix should be defined in 3-array sparse format and can depend on time t.
     * Matrix M is empty and should be filled with non-zero elements.
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(sparse_matrix &M, const double t) const = 0;
};

/*
 * Helper class to create identity Mass matrix of size N
 */
class MassMatrixIdentity : public MassMatrix
{
    const std::size_t _N{0}; // Indentity matrix size

public:
    MassMatrixIdentity(const std::size_t N) : MassMatrix(), _N(N) {}

    void operator()(sparse_matrix &M, const double t) const
    {
        M.A.resize(_N, 1.0);
        M.i.resize(_N); // Resize and then overwrite in a loop worked faster than reserve and push_back
        M.j.resize(_N);

        for (std::size_t i = 0; i < _N; ++i)
        {
            M.i[i] = i;
            M.j[i] = i;
        }
    }
};

/*
 * Helper class to create zero Mass matrix
 */
class MassMatrixZero : public MassMatrix
{
public:
    MassMatrixZero() : MassMatrix() {}

    void operator()(sparse_matrix &M, const double t) const
    {
        M.A.clear();
        M.i.clear();
        M.j.clear();
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_MASS_MATRIX_H
