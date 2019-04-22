/*
 * TODO: Description of the class
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"
#include "mass_matrix.h"

namespace daecpp_namespace_name
{

class TimeIntegrator
{

    int m_scheme;

public:
    /*
     * TODO: Description
     */
    TimeIntegrator(int scheme) : m_scheme(scheme) {}

    void operator()(RHS &rhs, MassMatrix &mass, state_type &x);
};

}  // namespace daecpp_namespace_name
