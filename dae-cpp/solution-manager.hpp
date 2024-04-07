/*
 * Solution Manager class.
 * Solution Manager functor will be called every time step providing the time `t` and
 * the corresponding solution `x` for further post-processing.
 * If the functor returns an integer != 0 (`true`), the computation will immediately stop.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_SOLUTION_MANAGER_H
#define DAECPP_SOLUTION_MANAGER_H

#include <algorithm>

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Default Solution Manager class
 */
struct SolutionManager
{
    /*
     * Solution Manager functor will be called every time step providing the time `t` and
     * the corresponding solution `x` for further post-processing.
     * If the functor returns an integer != 0 (`true`), the computation will immediately stop.
     * It does nothing by default and returns 0.
     */
    virtual int operator()(const state_vector &x, const double t)
    {
        return 0;
    }

    virtual ~SolutionManager() {}
};

/*
 * Solution holder class.
 * Holds solution vectors `x` and the corresponding times `t`.
 */
struct SolutionHolder
{
    // Vector of solution vectors `x`
    std::vector<state_vector> x;

    // Vector of the corresponding solution times `t`
    std::vector<double> t;

    // TODO: Print `x` and `t`
    // void print()

    // TODO: Save to a file
    // void save(const char delimiter = ',')
};

/*
 * Solution class.
 * Writes solution vector `x` and time `t` every time step or every specific time
 * from `t_output` vector into vectors `x` and `t`.
 */
class Solution : public SolutionManager
{
    // Output times
    std::vector<double> _t_out; // Need a copy for sorting

    // A reference to the solution holder object
    SolutionHolder &_sol;

public:
    /*
     * Writes solution vector `x` and time `t` every time step or every specific time
     * from `t_output` vector into vectors `x` and `t`.
     *
     * Parameters:
     *     `sol` - Solution holder object (SolutionHolder)
     *     `t_output` - (optional) a vector of output times for writing (`std::vector<double>`)
     */
    explicit Solution(SolutionHolder &sol, const std::vector<double> &t_output = {}) : _sol(sol), _t_out(t_output)
    {
        if (_t_out.size() > 0)
        {
            std::sort(_t_out.begin(), _t_out.end());
        }
    }

    /*
     * Solution Manager functor will be called every time step providing the time `t` and
     * the corresponding solution `x` for further post-processing.
     */
    int operator()(const state_vector &x_, const double t_)
    {
        if (_t_out.size() > 0)
        {
            if (!std::binary_search(_t_out.begin(), _t_out.end(), t_))
            {
                return 0;
            }
        }

        _sol.x.emplace_back(x_);
        _sol.t.emplace_back(t_);

        return 0;
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLUTION_MANAGER_H
