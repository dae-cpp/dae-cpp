/*
 * Jacobian matrix class.
 * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
 * This class is abstract and must be inherited.
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

// #include <autodiff/forward/real.hpp>
// #include <autodiff/forward/real/eigen.hpp>

// using namespace autodiff;

namespace daecpp_namespace_name
{

/*
 * Jacobian matrix class.
 * This class is abstract and must be inherited.
 */
class Jacobian
{
public:
    /*
     * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
     * Takes vector x and time t and returns the Jacobian matrix J.
     * Matrix J is empty and should be filled with non-zero elements.
     * This function is pure virtual and must be overriden to provide analytical Jacobian.
     */
    virtual void operator()(sparse_matrix &J, const state_type &x, const double t) const = 0;
};

/*
 * TODO: Numerical Jacobian.
 * Performs automatic (algorithmic) differentiation of the RHS using `autodiff` package.
 */
class JacobianNumerical : public Jacobian
{
    const RHS &_rhs; // The RHS for differentiation

    // The vector function with parameters for which the Jacobian is needed
    // VectorXreal f(const VectorXreal& f_, const state_type &x, const double t)
    // {
    // }

public:
    explicit JacobianNumerical(const RHS &rhs) : Jacobian(), _rhs(rhs) {}

    /*
     * TODO: Numerical Jacobian.
     * Performs automatic (algorithmic) differentiation of the RHS using `autodiff` package.
     */
    void operator()(sparse_matrix &J, const state_type &x, const double t) const
    {
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_JACOBIAN_H
