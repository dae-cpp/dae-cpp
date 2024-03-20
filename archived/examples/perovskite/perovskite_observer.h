/*
 * User-defined observer for perovskite problem.
 * Just prints out current time t and potential phi on the right boundary.
 * According to the boundary conditions, phi on the right boundary should be
 * equal to the current time exactly. This is what we check every time step.
 */

#pragma once

#include <iostream>

#include "../../src/solver.h"

class MySolver : public daecpp::Solver
{
public:
    MySolver(daecpp::RHS &rhs, daecpp::Jacobian &jac, daecpp::MassMatrix &mass,
             daecpp::SolverOptions &opt)
        : daecpp::Solver(rhs, jac, mass, opt)
    {
    }

    /*
     * Overloaded observer.
     * Receives current solution vector and the current time every time step.
     * Prints current time t and potential phi on the right boundary.
     */
    void observer(daecpp::state_type &x, const double t)
    {
        const MKL_INT size = (MKL_INT)(x.size());
        std::cout << " | t = " << t << " == phi(x=1) = " << x[size - 1];
    }
};
