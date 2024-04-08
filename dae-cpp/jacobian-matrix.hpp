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

#ifndef DAECPP_JACOBIAN_MATRIX_H
#define DAECPP_JACOBIAN_MATRIX_H

#include "sparse-matrix.hpp"
#include "vector-function.hpp"

namespace daecpp_namespace_name
{

/*
 * Jacobian matrix class.
 * This class is abstract and must be inherited.
 */
class JacobianMatrix
{
public:
    /*
     * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
     * Takes vector x and time t and returns the Jacobian matrix J.
     * Matrix J is empty and should be filled with non-zero elements.
     * This function is pure virtual and must be overriden to provide analytical Jacobian.
     */
    virtual void operator()(sparse_matrix &J, const state_vector &x, const double t) const = 0;

    virtual ~JacobianMatrix() {}
};

/*
 * Automatic (algorithmic) Jacobian class.
 * Performs algorithmic differentiation of the RHS using `autodiff` package.
 */
template <class RHS>
class JacobianAutomatic : public JacobianMatrix
{
    RHS &_rhs; // The RHS for differentiation (it is not moved because we need it in the solver)

public:
    explicit JacobianAutomatic(RHS &rhs) : JacobianMatrix(), _rhs(rhs) {}

    /*
     * Automatic (algorithmic) Jacobian.
     * Performs algorithmic differentiation of the RHS using `autodiff` package.
     */
    void operator()(sparse_matrix &J, const state_vector &x, const double t) const
    {
        std::size_t size = x.size(); // System size

        state_type x_(size); // Vectors of `dual` numbers are defined with `_` suffix

        // Conversion to dual numbers for automatic differentiation
        for (std::size_t k = 0; k < size; ++k)
        {
            x_[k] = x[k];
        }

        // The vector lambda-function with parameters for which the Jacobian is needed
        auto f = [&_rhs = _rhs, size](const state_type &x_, const double t)
        {
            state_type f_(size);
            _rhs(f_, x_, t);
            return f_;
        };

        // Dense Jacobian matrix
        Eigen::MatrixXd jac = autodiff::jacobian(f, wrt(x_), at(x_, t));

        // Convert dense matrix to sparse format
        for (std::size_t j = 0; j < size; ++j)
        {
            for (std::size_t i = 0; i < size; ++i)
            {
                const double val = jac(i, j);

                if (std::abs(val) > DAECPP_SPARSE_MATRIX_ELEMENT_TOLERANCE)
                {
                    J(i, j, val);
                }
            }
        }
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_JACOBIAN_MATRIX_H
