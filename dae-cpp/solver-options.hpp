/*
 * Defines DAE solver options.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#ifndef DAECPP_SOLVER_OPTIONS_H
#define DAECPP_SOLVER_OPTIONS_H

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

// Solver verbosity levels
namespace verbosity
{
enum VerbosityLevel
{
    off,
    normal,
    extra
};
} // namespace verbosity

/*
 * Defines DAE solver options
 */
struct SolverOptions
{
    // Verbosity level.
    // Can be `verbosity::off` (no output), `verbosity::normal` (prints basic information), or `verbosity::extra` (full output).
    // Default value is `verbosity::off`.
    verbosity::VerbosityLevel verbosity{verbosity::off};

    // Initial time step (0.01 by default).
    // Should be relatively small because the first time step is always first order.
    double dt_init{0.01};

    // Minimum time step (1e-10 by default).
    // If the solver has to reduce the time step below `dt_min` but still fails to converge, the simulation will stop.
    double dt_min{1e-10};

    // Maximum time step (1e10 by default).
    // The solver will not increase the time step above `dt_max`.
    double dt_max{1e10};

    // Absolute tolerance for the Newton algorithm. Default value is 1.0e-6.
    double atol{1e-6};

    // Relative tolerance for the Newton algorithm. Default value is 1.0e-6.
    double rtol{1e-6};

    // Maximum absolute error.
    // If the absolute error estimation exceeds `max_err_abs`, the Newton iterations will be considered as diverged.
    double max_err_abs{1e100};

    // Order of the BDF implicit numerical integration method:
    // 1 - first order BDF, 2 - BDF-2, 3 - BDF-3, 4 - BDF-4 (default).
    unsigned int BDF_order{4};

    // Non-linear solver algorithm:
    // 0 - Classic Newton method (usually the most stable but the slowest),
    // 1 (default) - Quasi-Newton method I (balanced with focus on stability, updates Jacobian and performs factorization every 2nd iteration),
    // 2 - Quasi-Newton method II (balanced with focus on speed, updates Jacobian and performs factorization every 3rd iteration),
    // 3 - Quasi-Newton method III (can be the fastest but less stable, may require tweaking the time step increase/decrease thresholds,
    //     updates Jacobian and performs factorization every 4th iteration).
    unsigned int Newton_scheme{1};

    // If `true`, the mass matrix will be updated only once (at the beginning of computation).
    // Setting this option as `true` slightly speeds up the computation, but the mass matrix must be static (i.e., it must be independent on time).
    // Default value is `false`.
    bool is_mass_matrix_static{false};

    // Maximum number of the Jacobian matrix updates and factorizations per time step.
    // If the algorithm fails to converge after `max_Jacobian_updates` Jacobian matrix updates and factorizations per time step,
    // the Newton iterations will be considered as diverged. The solver will try to roll back and decrease the time step.
    // Default value is 8.
    int max_Jacobian_updates{8};

    // If the Newton method fails to converge 'max_Newton_failed_attempts' times in a row,
    // the solver will try to update Jacobian matrix every single iteration next time step (for Quasi-Newton methods).
    // Default value is 3.
    unsigned int max_Newton_failed_attempts{3};

    // If `true`, the solver will try to recover from the linear solver failure by rolling back and decreasing the time step.
    // Default value is `true`.
    bool recover_from_linsolver_failure{true};

    // Time step amplification threshold delta (0 by default). Can be negative or positive integer number.
    // The solver increases the time step if the number of successful Newton iterations per time step is less than the threshold.
    // Example: If the solver increases the time step too often (too fast), decrease the time step amplification threshold by
    // setting `dt_increase_threshold_delta` to -1.
    int dt_increase_threshold_delta{0};

    // Time step reduction threshold delta (0 by default). Can be negative or positive integer number.
    // The solver decreases the time step if the number of successful Newton iterations per time step is greater than the threshold.
    // Example: If the solver struggles to converge but keeps doing many Newton iterations without reducing the time step,
    // decrease the time step reduction threshold by setting `dt_decrease_threshold_delta` to -1 (or -2 and less).
    int dt_decrease_threshold_delta{0};

    // Turns ON/OFF solution variability control.
    // Solution variability control tightens up the adaptive time stepping algorithm, and it is ON by default.
    // Switching it OFF can lead to a significant speed boost for big systems, but it can also lead to instability.
    bool solution_variability_control{true};

    // Solution variability tolerance.
    // Solution variability is the maximum relative change in the solution every time step.
    // If the absolute value of the solution is very close to zero, the relative change can be very high, leading to unnecessary small time steps.
    // To prevent this, the solver will ignore the values below the solution variability tolerance.
    // The default value is 1e-4. Consider increasing it if the solver reduces the time step too much.
    double variability_tolerance{1e-4};

    // Solution variability lower threshold.
    // The higher the value, the more likely the time step can be increased.
    // If the maximum relative change in the solution is above the lower threshold, the time step will NOT be increased.
    // Default value is 0.15.
    double variability_threshold_low{0.15};

    // Solution variability higher threshold.
    // The higher the value, the less likely the time step will be reduced.
    // If the maximum relative change in the solution is above the higher threshold, the time step will be reduced.
    // Default value is 0.5.
    double variability_threshold_high{0.5};

    // Time step amplification factor.
    // If the solution is stable enough, the solver will increase the time step by a factor of `dt_increase_factor`.
    // Default value is 2.0.
    double dt_increase_factor{2.0};

    // Time step reduction factor.
    // If the solution is not stable enough, the solver will decrease the time step by a factor of `dt_decrease_factor`.
    // Default value is 2.0.
    double dt_decrease_factor{2.0};

    // Number of threads
    unsigned int num_threads{1}; // TODO: Not used yet

    /*
     * Checks the user-defined parameters
     */
    void check() const
    {
        CHECK((verbosity >= 0) && (verbosity < 3), "Unknown dae-cpp verbosity level: " << verbosity);
        ASSERT(dt_init > 0.0, "Initial time step `dt_init` must be positive.");
        ASSERT(dt_min > 0.0, "Minimum time step `dt_min` must be positive.");
        ASSERT(dt_max > 0.0, "Maximum time step `dt_max` must be positive.");
        ASSERT(dt_min < dt_init, "Minimum time step `dt_min` must be less than the initial time step `dt_init`.");
        ASSERT(dt_min > 0.0, "Minimum time step `dt_min` should be positive.");
        ASSERT(dt_max >= dt_init, "Maximum time step `dt_max` must be greater than or equal to the initial time step `dt_init`.");
        ASSERT(atol > 0.0, "Absolute tolerance `atol` for the Newton algorithm must be positive.");
        ASSERT(rtol > 0.0, "Relative tolerance `rtol` for the Newton algorithm must be positive.");
        ASSERT(max_err_abs > 0.0, "Maximum absolute error `max_err_abs` must be positive. In most cases it should be a very big number.");
        ASSERT((BDF_order >= 1) && (BDF_order <= DAECPP_MAX_ORDER), "Unknown order of the BDF implicit numerical integration method: " << BDF_order << ".\n`BDF_order` must be 1, 2, 3, or 4.");
        ASSERT(Newton_scheme <= 10, "Non-linear solver algorithm (defined by `Newton_scheme` option) can be 0 (classic Newton method), 1, 2, or 3 (Quasi-Newton methods).");
        ASSERT(max_Jacobian_updates >= 4, "Maximum number of the Jacobian matrix updates and factorizations per time step should be at least 4.");
        ASSERT(max_Newton_failed_attempts > 0, "`max_Newton_failed_attempts` should be positive.");
        ASSERT(dt_increase_factor >= 1.0, "Time step amplification factor `dt_increase_factor` should be greater than 1.");
        ASSERT(dt_decrease_factor >= 1.0, "Time step reduction factor `dt_decrease_factor` should be greater than 1.");
        ASSERT(variability_tolerance >= 0.0, "Solution variability tolerance cannot be negative.");
        ASSERT(variability_threshold_low >= 0.0, "Solution variability lower threshold cannot be negative.");
        ASSERT(variability_threshold_high >= 0.0, "Solution variability higher threshold cannot be negative.");
        ASSERT(variability_threshold_low <= variability_threshold_high, "Solution variability lower threshold should be less than the higher threshold.");
        ASSERT(num_threads > 0, "Number of threads `num_threads` should be 1 or more.");
        CHECK(num_threads <= 32, "Using more than 32 threads can lead to a significant performance degradation.");

        if (verbosity > 0)
        {
            CHECK(Newton_scheme <= 3, "Quasi-Newton method " << Newton_scheme << " can lead to instability and inaccurate solution.");

            if (is_mass_matrix_static)
            {
                NOTE("Using static mass matrix...");
            }
        }
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_OPTIONS_H
