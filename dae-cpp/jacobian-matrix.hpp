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

#include "timer.hpp" // TODO: temporary here

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
};

/*
 * Numerical Jacobian.
 * Performs automatic (algorithmic) differentiation of the RHS using `autodiff` package.
 */
class JacobianNumerical : public JacobianMatrix
{
    const VectorFunction &_rhs; // The RHS for differentiation

public:
    explicit JacobianNumerical(const VectorFunction &rhs) : JacobianMatrix(), _rhs(rhs) {}

    /*
     * Numerical Jacobian.
     * Performs automatic (algorithmic) differentiation of the RHS using `autodiff` package.
     */
    void operator()(sparse_matrix &J, const state_vector &x, const double t) const
    {
        std::size_t size = x.size();

        state_type x__(size);

        for (std::size_t k = 0; k < size; ++k)
        {
            x__[k] = x[k];
        }

        // state_type f__(size);

        // The vector lambda-function with parameters for which the Jacobian is needed
        // auto f = [this](const state_type &x, const double &t, state_type &f)
        auto f = [&_rhs = _rhs](const state_type &x, const double t, const std::size_t size)
        {
            state_type f(size);
            // this->_rhs(f, x, t);
            _rhs(f, x, t);
            return f;
        };

        Eigen::MatrixXd jac;

        static double jac_numerical{0.0};
        {
            Timer timer(&jac_numerical);
            // jac = autodiff::jacobian(f, wrt(x__), at(x__, t, f__));
            jac = autodiff::jacobian(f, wrt(x__), at(x__, t, size));
        }
        std::cout << '|' << jac_numerical / 1000. << ' ';

        // Convert matrix to sparse format
        static double jac_numerical_2{0.0};
        {
            Timer timer(&jac_numerical_2);
            for (std::size_t j = 0; j < size; ++j)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    // double val = jac.coeffRef(i, j);
                    const double val = jac(i, j);

                    if (std::abs(val) > 1e-14)
                    {
                        J(val, i, j);
                    }
                }
            }
        }
        std::cout << jac_numerical_2 / 1000. << '|';
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_JACOBIAN_MATRIX_H
