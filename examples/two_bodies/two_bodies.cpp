/*
 * Two bodies with with masses m and 10*m head to each other.
 * The friction force F is constant and depends on the sign of velocity v:
 * F = -sign(v) * F0.
 * After collision, the body with mass m bounces in the other direction
 * (and hence flips the sign of the friction force).
 *
 * The system of equations is the following:
 *
 * dv1/dt = -sign(v1) * F0 / m
 * dv2/dt = -sign(v2) * F0 / (10 * m)
 * dx1/dt = v1
 * dx2/dt = v2
 *
 * Here v1 and x1 is the velocity and the coordinate of the first body,
 * v2 and x2 is the velocity and the coordinate of the second body.
 *
 * When x1 reaches x2, the velocity direction of the body with mass m changes to
 * the opposite.
 *
 * As an example, let us consider the following initial conditions:
 * v1 = 10, v2 = -2, x1 = 0, x2 = 5 for t = 0.
 */

#include <iostream>
#include <cmath>

#include "../../src/solver.h"  // the main header of dae-cpp library solver

using namespace daecpp;

// python3 + numpy + matplotlib should be installed in order to enable plotting
// #define PLOTTING

#ifdef PLOTTING
#include "../../src/external/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

/*
 * RHS of the problem
 * =============================================================================
 */
class MyRHS : public RHS
{
    const double f0 = 10.0;
    const double m  = 1.0;

public:
    /*
     * Receives the current solution vector x and the current time t.
     * Defines the RHS f.
     */
    void operator()(const state_type &x, state_type &f, const double t)
    {
        double F1 = -std::copysign(1.0, x[0]) * f0;
        double F2 = -std::copysign(1.0, x[1]) * f0;

        f[0] = F1 / m;
        f[1] = F2 / (10.0 * m);
        f[2] = x[0];
        f[3] = x[1];
    }
};

/*
 * (Optional) Observer
 * =============================================================================
 * Every time step checks that
 * (1) x*x + y*y = 1, and
 * (2) x(t) - sin(t) = 0 for t <= pi/2, x(t) = 1 for t > pi/2
 * and prints solution and errors to console.
 */
class MySolver : public Solver
{
public:
    MySolver(RHS &rhs, Jacobian &jac, MassMatrix &mass, SolverOptions &opt)
        : Solver(rhs, jac, mass, opt)
    {
    }

#ifdef PLOTTING
    state_type x_axis, v1, v2, x1, x2;  // For plotting
#endif

    /*
     * Overloaded observer.
     * Receives current solution vector and the current time every time step.
     */
    void observer(state_type &x, const double t)
    {
        std::cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << '\t'
                  << x[3] << '\n';

        // When x1 reaches x2, change v1 velocity sign to negative
        if(x[2] > x[3])
        {
            x[0] = -std::abs(x[0]);
        }

#ifdef PLOTTING
        // Save data for plotting
        x_axis.push_back(t);
        v1.push_back(x[0]);
        v2.push_back(x[1]);
        x1.push_back(x[2]);
        x2.push_back(x[3]);
#endif
    }
};

/*
 * (Optional) Analytical Jacobian in simplified 3-array sparse format
 * =============================================================================
 */
class MyJacobian : public Jacobian
{
public:
    explicit MyJacobian(RHS &rhs) : Jacobian(rhs) {}

    /*
     * Receives the current solution vector x and the current time t. Defines
     * the analytical Jacobian matrix J.
     */
    void operator()(sparse_matrix_holder &J, const state_type &x,
                    const double t)
    {
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 * Returns '0' if solution comparison is OK or '1' if solution error is above
 * the acceptable tolerances.
 */
int main()
{
    // Solution time 0 <= t <= t1
    double t1 = 1.0;

    // Define the state vector for 4 equations
    state_type x(4);

    // Initial conditions
    x[0] = 10;  // Initial velocity v1
    x[1] = -2;  // Initial velocity v2
    x[2] = 0;   // Initial coordinate x1
    x[3] = 5;   // Initial coordinate x2

    // Set up the RHS of the problem.
    // Class MyRHS inherits abstract RHS class from dae-cpp library.
    MyRHS rhs;

    // Set up the Mass Matrix of the problem. In this case this matrix is
    // identity, so we can use a helper class provided by the library.
    MassMatrixIdentity mass(4);

    // Create an instance of the solver options and update some of the solver
    // parameters defined in solver_options.h
    SolverOptions opt;

    opt.time_stepping = 1;  // Use simple stability-based adaptive time stepping
                            // algorithm
    opt.bdf_order = 1;      // Use BDF-1
    opt.verbosity = 0;      // Suppress output to screen (we have our own output
                            // defined in Observer function above)
    opt.dt_init = 0.01;     // Change the initial time step
    opt.dt_max  = 0.01;     // Restrict the maximum time step

    // We can override Jacobian class from dae-cpp library and provide
    // analytical Jacobian.
    // MyJacobian jac(rhs);

    // Or we can use numerically estimated Jacobian with the given tolerance.
    Jacobian jac(rhs, 1e-8);

    // Create an instance of the solver with particular RHS, Mass matrix,
    // Jacobian and solver options
    MySolver solve(rhs, jac, mass, opt);

    // Now we are ready to solve the set of DAEs
    std::cout << "\nStarting DAE solver...\n";
    std::cout << "time\tv1\tv2\tx1\tx2\n";

    // Solve the system
    int status = solve(x, t1);

    // Plot the solution
#ifdef PLOTTING
    // plt::figure();
    // plt::figure_size(640, 480);
    // plt::named_semilogx("x", solve.x_axis, solve.x0);
    // plt::named_semilogx("y", solve.x_axis, solve.x1);
    // plt::xlabel("time");
    // plt::title("Two bodies");
    // plt::grid(true);
    // plt::legend();

    // // Save figure
    // const char *filename = "two_bodies.png";
    // std::cout << "Saving result to " << filename << "...\n";
    // plt::save(filename);
#endif

    if(status)
        std::cout << "...Test FAILED\n\n";
    else
        std::cout << "...done\n\n";

    return status;
}
