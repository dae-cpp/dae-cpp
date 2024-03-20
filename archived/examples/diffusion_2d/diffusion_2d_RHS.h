/*
 * The RHS definition
 */

#pragma once

#include "../../src/RHS.h"

class MyRHS : public daecpp::RHS
{
    const MKL_INT m_N;  // Number of cells along axis
    const double  m_D;  // Diffusion coefficient (dimensionless)

public:
    MyRHS(MKL_INT N, double D) : daecpp::RHS(), m_N(N), m_D(D) {}

    void operator()(const daecpp::state_type &x, daecpp::state_type &f,
                    const double t);
};
