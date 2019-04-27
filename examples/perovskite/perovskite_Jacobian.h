/*
 * TODO: Description of the class
 */

#pragma once

#include "../../src/jacobian.h"
#include "perovskite_parameters.h"

class MyJacobian : public daecpp::Jacobian
{
    MyParams m_p;

public:
    MyJacobian(daecpp::RHS &rhs, MyParams p) : daecpp::Jacobian(rhs), m_p(p) {}
    void operator()(daecpp::sparse_matrix_holder &J, const daecpp::state_type &x,
                    const double t);
};
