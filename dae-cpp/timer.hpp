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

#pragma once

#include <chrono>

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Specific timers (a singleton)
 */
class Timers
{
    static Timers *_instance;

    Timers() {}

public:
    double total{0.0}; // Total time

    static Timers *getInstance()
    {
        if (_instance == nullptr)
        {
            _instance = new Timers();
        }
        return _instance;
    }
};

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

    // Stores time
    double *_t;

public:
    Timer(double *t) : _t(t) {}

    ~Timer()
    {
        *_t += std::chrono::duration_cast<time_unit>(clock::now() - tic).count() * 1e-3; // ms
    }
};

} // namespace daecpp_namespace_name
