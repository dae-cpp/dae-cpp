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
    // Should be relatively small because the first time step is always first order.
    double dt_init{0.01};

    // Minimum time step.
    // If the solver has to reduce the time step below `dt_min` but still fails to converge, the simulation will stop.
    double dt_min{1e-10};

    // Maximum time step.
    // The solver will not increase the time step above `dt_max`.
    double dt_max{1e10};

    // Absolute tolerance for the Newton algorithm.
    double atol{1.0e-6};

    // Relative tolerance for the Newton algorithm.
    double rtol{1.0e-6};

    // If the absolute error estimation exceeds `max_err_abs`, the Newton iterations will be considered as diverged.
    double max_err_abs{1e100};

    // Order of BDF implicit numerical integration method:
    // 1 - first order BDF, 2 - BDF-2, 3 - BDF-3, 4 - BDF-4 (default).
    unsigned int BDF_order{4};

    // Non-linear solver algorithm:
    // 0 - Classic Newton method (usually the most stable but the slowest),
    // 1 - Quasi-Newton method I (balanced with focus on stability, updates Jacobian and performs factorization every 2nd iteration),
    // 2 - Quasi-Newton method II (balanced with focus on speed, updates Jacobian and performs factorization every 3rd iteration),
    // 3 - Quasi-Newton method III (can be the fastest but less stable, may require tweaking the time step increase/decrease thresholds,
    //     updates Jacobian and performs factorization every 4th iteration).
    unsigned int Newton_scheme{1};

    // If `true`, the mass matrix will be updated only once, at the beginning of computation.
    // Setting this option as `true` slightly speeds up the computation if the mass matrix is static (i.e., does not depend on time).
    bool is_mass_matrix_static{false};

    // Maximum number of Newton iterations.
    // If the algorithm fails to converge after `max_Newton_iter` iterations, the Newton iterations will be considered as diverged.
    // The solver will try to roll back and decrease the time step. TODO: Conflicts with dynamic thresholds.
    unsigned int max_Newton_iter{15};

    // If Newton method fails to converge within 'max_Newton_iter' iterations
    // 'newton_failed_attempts' times in a row, the solver will try to update
    // Jacobian every single iteration next time step.
    // int max_Newton_failed_attempts = 3;

    // Time step amplification threshold. Can be negative or positive
    unsigned int dt_increase_threshold_delta{0};

    // Time step reduction threshold
    unsigned int dt_decrease_threshold_delta{0};

    // Time step amplification factor
    double dt_increase_factor{2.0};

    // Time step reduction factor
    double dt_decrease_factor{2.0};

    // Number of threads
    int num_threads{1};

    /*
     * TODO: Checks the user-defined parameters
     */
    void check() const
    {
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_OPTIONS_H
