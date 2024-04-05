/*
 * Observer class.
 * Observer operator () will be called every time step and provide the time `t`
 * and the corresponding solution `x`.
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

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Observer class
 */
class Observer
{
public:
    /*
     * Observer operator () will be called every time step and provide the time `t`
     * and the corresponding solution `x`.
     * It does nothing by default.
     */
    virtual void operator()(const state_vector &x, const double t) {}

    virtual ~Observer() {}
};

/*
 * Saves solution vector `x` and time `t` every time step into vectors
 */
class Solution : public Observer
{
    std::vector<state_vector> x_sol;
    std::vector<double> t_sol;

public:
    /*
     * Saves solution vector `x` and time `t` every time step into vectors.
     *
     * Parameters:
     *     `x_sol` - vector of solution vectors written every time step (`std::vector<daecpp::state_vector>`)
     *     `t_sol` - vector of the corresponding solution times (`std::vector<double>`)
     */
    Solution(std::vector<state_vector> &x_sol, std::vector<double> &t_sol) : x_sol(x_sol), t_sol(t_sol) {}

    void operator()(const state_vector &x, const double t)
    {
        x_sol.emplace_back(x);
        t_sol.emplace_back(t);
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_OBSERVER_H
