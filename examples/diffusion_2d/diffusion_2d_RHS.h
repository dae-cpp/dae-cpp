/*
 * The RHS definition
 */

#pragma once

#include "../../src/RHS.h"

class MyRHS : public daecpp::RHS
{
    MKL_INT m_N;
    double  m_D;

public:
    MyRHS(MKL_INT N, double D) : daecpp::RHS(), m_N(N), m_D(D) {}

    void operator()(const daecpp::state_type &x, daecpp::state_type &f,
                    const double t);
};
