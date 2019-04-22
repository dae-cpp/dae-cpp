/*
 * TODO: Description
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

void TimeIntegrator::operator()(RHS &rhs, MassMatrix &mass, state_type &x)
{
    // mass matrix * dt (element-wise)
    // * rhs vector

    // x = rhs + (x-x_prev)/dt * m
}

}  // namespace daecpp_namespace_name
