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
    const double m_tol = 1.0e-3;  // Default tolerance
    const double m_eps = 1.0e-6;  // The order of the rounding unit
#else
    const double m_tol = 1.0e-6;   // Default tolerance
    const double m_eps = 1.0e-13;  // The order of the rounding unit
#endif

public:
    explicit Jacobian(RHS &rhs) : m_rhs(rhs) {}
    Jacobian(RHS &rhs, const double tol) : m_rhs(rhs), m_tol(tol)
    {
        // TODO: Check user's tol parameter. Too small tol may lead to crash.
    }

    /*
     * Can be overriden to provide analytical Jacobian
     */
    virtual void operator()(sparse_matrix_holder &J, const state_type &x,
                            const double t);

    /*
     * Helper function to show Jacobian structure
     */
    void print(const state_type &x, const double t);
};

}  // namespace daecpp_namespace_name
