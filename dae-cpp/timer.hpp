/*
 * Timer class -- measures time in ms.
 *
 * Usage example:
 *
 *     double time{};
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
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#ifndef DAECPP_TIMER_H
#define DAECPP_TIMER_H

#include <chrono>

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Main timer class. Measures time in ms.
 */
class Timer
{
    // Sets up clock
    using clock = std::chrono::steady_clock;
    using time_unit = std::chrono::microseconds;

    // Starts timer
    std::chrono::time_point<clock> m_tic = clock::now();

    // Points to a specific timer
    double *m_t;

public:
    explicit Timer(double *t) : m_t(t) {}

    ~Timer()
    {
        *m_t += std::chrono::duration_cast<time_unit>(clock::now() - m_tic).count() * 1e-3; // ms
    }
};

namespace core
{
namespace timer
{

enum specific_timers_enum
{
    init,            // Initialization time
    time_derivative, // Time to compute the time derivative approximation
    rhs,             // Time to compute and convert the RHS
    mass,            // Time to compute and convert the Mass matrix
    jacobian,        // Time to compute and convert the Jacobian matrix
    linear_algebra,  // Time to perform linear algebra operations
    factorization,   // Time to perform sparse matrix factorization
    linear_solver,   // Time spent by the linear solver
    error_check,     // Time to perform error checks
    manager          // Time spent to call Solution Manager
};

/*
 * Specific timers
 */
class Time
{
    const int m_N{10}; // Total number of the specific timers

public:
    double total{0.0}; // Total time spent by the DAE solver

    std::vector<double> timers; // A vector of specific timers

    Time() : timers(std::vector<double>(m_N, 0.0)) {}

    /*
     * Returns time spent for initialization and other calculations not covered by the specific timers
     */
    double other() const
    {
        double sum{0.0};
        for (const auto t : timers)
        {
            sum += t;
        }
        return total - sum;
    }
};

} // namespace timer
} // namespace core
} // namespace daecpp_namespace_name

#endif // DAECPP_TIMER_H
