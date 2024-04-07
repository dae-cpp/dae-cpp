/*
 * Observer class.
 * Observer functor will be called every time step providing the time `t` and
 * the corresponding solution `x` for further post-processing.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_OBSERVER_H
#define DAECPP_OBSERVER_H

#include <algorithm>

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Observer class
 */
struct DefaultObserver
{
    /*
     * Observer functor will be called every time step providing the time `t` and
     * the corresponding solution `x` for further post-processing.
     * It does nothing by default.
     */
    virtual void operator()(const state_vector &x, const double t) {}

    virtual ~DefaultObserver() {}
};

/*
 * Saves solution vector `x` and time `t` every time step into vectors `x` and `t`
 */
class Solution : public DefaultObserver
{
    std::vector<double> _t_out;

public:
    // Vector of solution vectors `x` written every time step
    std::vector<state_vector> x;

    // Vector of the corresponding solution times `t` written every time step
    std::vector<double> t;

    /*
     * Saves solution vector `x` and time `t` every time step into vectors `x` and `t`.
     *
     * Parameter:
     *     `t_output` - (optional) a vector of output times for writing (`std::vector<double>`)
     */
    Solution(const std::vector<double> &t_output = {}) : _t_out(t_output)
    {
        if (_t_out.size() > 0)
        {
            std::sort(_t_out.begin(), _t_out.end());
        }
    }

    void operator()(const state_vector &x_, const double t_)
    {
        if (_t_out.size() > 0)
        {
            if (!std::binary_search(_t_out.begin(), _t_out.end(), t_))
            {
                return;
            }
        }
        x.emplace_back(x_);
        t.emplace_back(t_);
    }

    // TODO: Print `x` and `t`
    // void print()

    // TODO: Save to a file
    // void save(const char delimiter = ',')
};

} // namespace daecpp_namespace_name

#endif // DAECPP_OBSERVER_H
