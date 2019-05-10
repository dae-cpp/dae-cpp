/*
 * Calculates Jacobian matrix.
 * If not overridden by a user, performs numerical differentiation of the RHS
 * with the given tolerance to estimate numerical Jacobian matrix.
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"

namespace daecpp_namespace_name
{

class Jacobian
{

    RHS &m_rhs;

#ifdef DAE_SINGLE
    double m_tol = 1.0e-3;
#else
    double m_tol = 1.0e-6;
#endif

public:
    Jacobian(RHS &rhs) : m_rhs(rhs) {}
    Jacobian(RHS &rhs, const double tol) : m_rhs(rhs), m_tol(tol) {}

    /*
     * Can be overriden to provide analytical Jacobian.
     */
    virtual void operator()(sparse_matrix_holder &J, const state_type &x,
                            const double t);

    /*
     * TODO: Helper function to show Jacobian structure
     */
};

}  // namespace daecpp_namespace_name
