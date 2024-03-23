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
namespace core
{

/*
 * Specific timers (a singleton)
 */
class Timers
{
public:
    double total{0.0}; // Total time

    static Timers &get()
    {
        static Timers instance;
        volatile int dummy{}; // So the function call will not be optimised away
        return instance;
    }

private:
    Timers() = default;
    ~Timers() = default;
    Timers(const Timers &) = delete;
    Timers &operator=(const Timers &) = delete;
};

} // namespace core

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

} // namespace daecpp_namespace_name

#endif // DAECPP_TIMER_H
