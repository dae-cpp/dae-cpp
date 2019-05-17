/*
 * The RHS definition
 */

#pragma once

#include "../../src/RHS.h"
#include "perovskite_parameters.h"

class MyRHS : public daecpp::RHS
{
    MyParams m_p;

public:
    MyRHS(MyParams p) : daecpp::RHS(), m_p(p) {}

    void operator()(const daecpp::state_type &x, daecpp::state_type &f,
                    const double t);

    // Optional: override solver stop condition (it is always FALSE by default)
    bool stop_condition(const daecpp::state_type &x, const double t);
};
