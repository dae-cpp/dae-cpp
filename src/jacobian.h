/*
 * TODO: Description of the class
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"

namespace daecpp_namespace_name
{

class Jacobian
{

    RHS &m_rhs;

    double m_tol = 1.0e-6;

public:
    /*
     * TODO: Description
     */
    Jacobian(RHS &rhs) : m_rhs(rhs) {}
    Jacobian(RHS &rhs, const double tol) : m_rhs(rhs), m_tol(tol) {}

    /*
     * TODO: Numerical differentiation of RHS with the given tolerance.
     * Can be overriden to provide analytical Jacobian.
     * x is not a constant here to avoid unnecessary array copying.
     */
    virtual void operator()(sparse_matrix_holder &J, const state_type &x, const double t);

    /*
     * TODO: Helper function to print Jacobian
     */
};

}  // namespace daecpp_namespace_name
