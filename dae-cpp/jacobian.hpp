/*
 * Jacobian matrix class.
 * If not overridden by the user, performs automatic (algorithmic) differentiation of the RHS
 * using autodiff package to estimate Jacobian matrix numerically.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_JACOBIAN_H
#define DAECPP_JACOBIAN_H

#include "rhs.hpp"
#include "sparse-matrix.hpp"

namespace daecpp_namespace_name
{

class Jacobian
{
    const RHS &_rhs;

    RHS_empty _rhs_empty;

    // TODO: Matrix converter
    // /*
    //  * Sparse matrix converter from simple three-array format to Intel MKL
    //  * three array format.
    //  * Input: matrix holder M with simple three-array format
    //  * Output: matrix holder M with Intel MKL three-array format
    //  */
    // void m_matrix_converter(daecpp::sparse_matrix &M);

public:
    explicit Jacobian(const RHS &rhs) : _rhs(rhs) {}
    Jacobian() : _rhs(_rhs_empty) {}

    /*
     * Can be overriden to provide analytical Jacobian
     */
    virtual void operator()(sparse_matrix &J, const state_type &x, const double t) const;

    // /*
    //  * Helper function to show Jacobian structure on screen (in sparse format)
    //  */
    // void print(const state_type &x, const double t);

    // /*
    //  * Helper function to write Jacobian matrix to a file (in dense format)
    //  */
    // void dump(const state_type &x, const double t);

    // /*
    //  * Helper function to compare two Jacobians and write the differences.
    //  * Comparison will be made with the external Jacobian jac (usually,
    //  * numerical Jacobian) using vector x at time t with the given tolerance.
    //  * Returns the number of differences found.
    //  */
    // int compare(Jacobian jac, const state_type &x, const double t,
    //             const double tol);
};

/*
 * TODO: Numerical Jacobian. Parallel version.
 * Calls RHS up to 2*N times, hence O(N^2) operations.
 */
void Jacobian::operator()(sparse_matrix &J, const state_type &x, const double t) const
{
}

} // namespace daecpp_namespace_name

#endif // DAECPP_JACOBIAN_H
