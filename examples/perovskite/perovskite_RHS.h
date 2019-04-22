/*
 * TODO: Description of the class
 */

#pragma once

#include "../../src/RHS.h"
#include "perovskite_parameters.h"

class MyRHS : public daecpp::RHS
{
    MyParams m_p;

public:
    MyRHS(MyParams p) : RHS(), m_p(p) {}
    void operator()(const daecpp::state_type& x, daecpp::state_type& f,
                    const double t);
};
