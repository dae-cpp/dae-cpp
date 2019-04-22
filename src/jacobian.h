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

    double m_tol;

public:
    /*
     * TODO:
     */
    Jacobian(RHS &rhs) : m_rhs(rhs) {}
    Jacobian(RHS &rhs, const double tol) : m_rhs(rhs), m_tol(tol) {}

    /*
     * TODO: Numerical differentiation of RHS with the give tolerance
     * Can be overriden to provide analytical Jacobian
     * x is not a constant here to avoid unnecessary array copying
     */
    void operator()(state_type &J, vector_type_int &ia, vector_type_int &ja,
                    state_type &x, const double t);

    /*
     * TODO: Helper function to print Jacobian
     */
};

}  // namespace daecpp_namespace_name
