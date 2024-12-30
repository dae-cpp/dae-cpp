/*
 * Solves equation of flame propagation:
 *
 * | y' = y^2 - y^3
 *
 * Initial condition (t = 0):
 * | y = delta
 *
 * Solution interval:
 * 0 < t < (2 / delta)
 *
 * When you light a match, the ball of flame grows rapidly until it reaches a critical size.
 * Then it remains at that size because the amount of oxygen being consumed by the combustion
 * in the interior of the ball balances the amount of oxigen available through the surface.
 * The scalar variable y(t) represents the radius of the ball.
 *
 * If parameter `delta` is small (1e-4), the problem is stiff.
 *
 * See, for example,
 * https://uk.mathworks.com/company/technical-articles/stiff-differential-equations.html
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

// Main dae-cpp header
#include <dae-cpp/solver.hpp>

// dae-cpp namespace
using namespace daecpp;

/*
 * The vector-function (RHS) of the problem:
 *
 * f(y,t) = y^2 - y^3
 */
struct MyRHS
{
    /*
     * Defines the RHS (vector function) `f` of the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the RHS vector `f`.
     * Vector `f` is already pre-allocated with `f.size() == x.size()`.
     */
    void operator()(state_type &f, const state_type &x, const double t)
    {
        state_value y = x[0]; // Note `state_value` (not `double`) because the solver will automatically differentiate `f` w.r.t. `y`
        f[0] = y * y - y * y * y;
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 */
int main()
{
    MyRHS rhs; // The vector-function object

    System my_system(MassMatrixIdentity(1), rhs); // Defines the DAE system object, note the identity 1x1 mass matrix

    const double delta = 1e-4; // Parameter `delta`: 0.01 - not very stiff problem, 1e-4 - stiff problem

    state_vector y0{delta};    // Initial condition: y(0) = delta
    double t_end{2.0 / delta}; // Solution interval: 0 < t < (2 / delta)

    // Update the solver options
    my_system.opt.dt_init = t_end / 100.0;           // Increase the initial time step (but not too much, because the first time step is of 1st order)
    my_system.opt.dt_increase_threshold_delta = -2;  // Try to increase the time step less often
    my_system.opt.variability_threshold_low = 0.10;  // Do not increase the time step if the solution changes more than 10%
    my_system.opt.variability_threshold_high = 0.10; // Decrease the time step if the solution changes more than 10%
    // my_system.opt.verbosity = verbosity::extra;   // Prints the time step information on screen

    // Solves the DAE system `my_system` with the given initial condition `y0` and time `t_end`
    my_system.solve(y0, t_end);

    // Prints solution on screen
    my_system.sol.print();

    return my_system.status;
}
