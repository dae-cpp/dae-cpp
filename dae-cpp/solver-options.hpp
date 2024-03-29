/*
 * Defines DAE solver options.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_SOLVER_OPTIONS_H
#define DAECPP_SOLVER_OPTIONS_H

#include "typedefs.hpp"

#define DAECPP_MAX_ORDER 4
#define DAECPP_TIMESTEP_ROUNDING_ERROR 1e-14

namespace daecpp_namespace_name
{

namespace verbosity
{
enum VerbosityLevel
{
    off,
    info,
    extra
};
} // namespace verbosity

/*
 * Defines DAE solver options
 */
struct SolverOptions
{
    // Verbosity level.
    // Can be `off` (no output), `info` (prints basic information), or `extra` (full output).
    // Default value is `off`.
    verbosity::VerbosityLevel verbosity{verbosity::off};

    // Initial time step.
    // Should be relatively small, since the first time step is always first order.
    double dt_init{0.01};

    // Minimum time step
    double dt_min = 1e-10;

    // Maximum time step
    double dt_max = 1e10;

#ifdef DAECPP_SINGLE
    double atol = 1.0e-3;      // Absolute tolerance for the Newton algorithm
    double rtol = 1.0e-6;      // Relative tolerance for the Newton algorithm
    double max_value = 1.0e20; // Solution shouldn't be higher than this, otherwise, it is considered as diverged
#else
    double atol = 1.0e-6;       // Absolute tolerance for the Newton algorithm
    double rtol = 1.0e-6;       // Relative tolerance for the Newton algorithm
    double max_value = 1.0e100; // Solution shouldn't be higher than this
#endif

    // Order of BDF implicit numerical integration method:
    // 1 - first order BDF, 2 - BDF-2, ..., 6 - BDF-6
    // Default is BDF-2 since it fully supports variable time stepping
    int BDF_order = 4;

    //
    int linear_solver = 0;

    // Perform Jacobian update, Reordering, Symbolic and Numerical Factorization
    // every Newton iteration. Changing to 'false' can increase speed but also
    // can lead to instability.
    // bool fact_every_iter = true;

    // If fact_every_iter = false, update Jacobian every fact_iter Newton
    // iterations
    // int fact_iter = 15;

    // Maximum number of Newton iterations. If the Newton method fails to
    // converge after max_Newton_iter iterations, the solver reduces time step
    // and tries to make the current step again.
    unsigned int max_Newton_iter{15};

    // If Newton method fails to converge within 'max_Newton_iter' iterations
    // 'newton_failed_attempts' times in a row, the solver will try to update
    // Jacobian every single iteration next time step.
    // int max_Newton_failed_attempts = 3;

    // Simple Adaptive Time Stepping options
    int dt_increase_threshold = 2;   // Time step amplification threshold
    int dt_decrease_threshold = 7;   // Time step reduction threshold
    double dt_increase_factor = 2.0; // Time step amplification factor
    double dt_decrease_factor = 2.0; // Time step reduction factor

    int num_threads = 1;

    /*
     * Checks the user-defined parameters
     */
    void check() const
    {
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_OPTIONS_H
