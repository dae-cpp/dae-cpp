/*
 * Timer class -- measures time in ms.
 *
 * Usage example:
 *
 *     double time = 0.0;
 *     {
 *         Timer timer(&time);
 *         << TASK >>
 *     }
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_TIMER_H
#define DAECPP_TIMER_H

#include <chrono>

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Main timer class
 */
class Timer
{
    // Sets up clock
    using clock = std::chrono::steady_clock;
    using time_unit = std::chrono::microseconds;

    // Starts timer
    std::chrono::time_point<clock> tic = clock::now();

    // Points to a specific timer
    double *_t;

public:
    Timer(double *t) : _t(t) {}

    ~Timer()
    {
        *_t += std::chrono::duration_cast<time_unit>(clock::now() - tic).count() * 1e-3; // ms
    }
};

namespace core
{

/*
 * Specific timers
 */
struct Time
{
    double total{0.0}; // Total time spent by the DAE solver

    // TODO: Might be worth combining them into something like std::map?
    double init{0.0};            // Initialization time
    double time_derivative{0.0}; // Time to compute the time derivative approximation
    double rhs{0.0};             // Time to compute and convert the RHS
    double mass{0.0};            // Time to compute and convert the Mass matrix
    double jacobian{0.0};        // Time to compute and convert the Jacobian matrix
    double linear_algebra{0.0};  // Time to perform linear algebra operations
    double factorization{0.0};   // Time to perform sparse matrix factorization
    double linear_solver{0.0};   // Time spent by the linear solver
    double error_check{0.0};     // Time to perform error checks
    double history{0.0};         // Time to update the solution vector history

    /*
     * Returns time spent for initialization and other calculations not covered by the specific timers
     */
    double other() const
    {
        double sum = init + time_derivative + rhs + mass + jacobian + linear_algebra + factorization + linear_solver + error_check + history;
        return total - sum;
    }
};

} // namespace core
} // namespace daecpp_namespace_name

#endif // DAECPP_TIMER_H
