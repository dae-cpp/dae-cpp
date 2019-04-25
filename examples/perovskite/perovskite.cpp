/*
 * TODO: Description of the problem
 * Keywords: perovskite solar cell, potential, ion concentration,
 * singular mass matrix
 *
 */

#include <chrono>
#include <iostream>

#include "../../src/solver.h"  // dae-cpp library solver

#include "perovskite_RHS.h"         // RHS of the problem
#include "perovskite_Mass.h"        // Mass Matrix definition
#include "perovskite_Jacobian.h"    // Jacobian of the problem
#include "perovskite_parameters.h"  // Parameter container of the problem

namespace dae = daecpp;  // A shortcut to dae-cpp library namespace

// python3 + numpy + matplotlib should be installed in order to enable plotting
//#define PLOTTING

#ifdef PLOTTING
#include "../../src/external/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

int main()
{
    // These parameters can be obtained from a parameter file or
    // as command line options. Here for simplicity we define them as constants.
    const int    N      = 4000;   // Number of cells
    const double L      = 1.0;    // Space interval length
    const double lambda = 1.0;    // Lambda parameter
    const double t1     = 100.0;  // Integration time (0 < t < t1)

    // Pass the parameters to the container
    MyParams p(N, L, lambda, t1);

    std::cout << "N = " << p.N << "; lambda = " << p.lambda << "; t = " << p.t1
              << '\n';

    using clock     = std::chrono::high_resolution_clock;
    using time_unit = std::chrono::milliseconds;

    // Define state vectors. Here 2*N is the total number of the equations.
    // We are going to carry out two independent simulations: with analytical
    // Jacobian and with numerically estimated one, hence two vectors.
    dae::state_type x1(2 * N);
    dae::state_type x2(2 * N);

    // Initial conditions
    for(int i = 0; i < N; i++)
    {
        x1[i]     = 1.0;  // for P - ion concentration
        x1[i + N] = 0.0;  // for Phi - potential
    }
    x2 = x1;  // x1 and x2 will be overwritten by the solver

    // Set up the RHS of the problem.
    // Class MyRHS inherits abstract RHS class from dae-cpp library.
    MyRHS rhs(p);

    // We can override Jacobian class from dae-cpp library
    // and provide analytical Jacobian
    MyJacobian jac(rhs, p);

    // Set up the Mass Matrix of the problem.
    // MyMassMatrix inherits abstract MassMatrix class from dae-cpp library.
    MyMassMatrix mass(N);

    // Update some of the solver options
    dae::SolverOptions opt;
    opt.atol = 1.0e-6;  // Absolute tolerance

    // Create an instance of the solver
    dae::Solver solve(rhs, jac, mass, opt, p.t1);

    // Now we are ready to solve the set of DAEs
    std::cout << "Starting DAE solver...\n";

    {
        auto tic0 = clock::now();
        solve(x1);
        auto tic1 = clock::now();

        std::cout
            << "Solver execution time: "
            << std::chrono::duration_cast<time_unit>(tic1 - tic0).count() /
                   1000.0
            << " sec." << '\n';
    }

    std::cout << "Solution check:\n";
    std::cout << x1[0] << ' ' << x1[N - 1] << '\n';
    std::cout << x1[N] << ' ' << x1[2 * N - 1] << '\n';

    // If we don't provide analytical Jacobian we need to estimate it
    // with a given tolerance:
    dae::Jacobian jac_est(rhs, 1.0e-6);

    // Create a new instance of the solver for estimated Jacobian
    dae::Solver solve_slow(rhs, jac_est, mass, opt, p.t1);

    // Solve the set of DAEs again
    std::cout << "\nStarting DAE solver with estimated Jacobian...\n";

    {
        auto tic0 = clock::now();
        solve_slow(x2);
        auto tic1 = clock::now();

        std::cout
            << "Solver execution time: "
            << std::chrono::duration_cast<time_unit>(tic1 - tic0).count() /
                   1000.0
            << " sec." << '\n';
    }

    std::cout << "Solution check:\n";
    std::cout << x2[0] << ' ' << x2[N - 1] << '\n';
    std::cout << x2[N] << ' ' << x2[2 * N - 1] << '\n';

#ifdef PLOTTING
    state_type x_axis(N), p(N), phi(N), d2phi(N);
    for(int i = 0; i < N; i++)
    {
        x_axis[i] = (0.5 + i) / N;
        p[i]      = x[i];
        phi[i]    = x[i + N];
        if(i > 0 && i < (N - 1))
            d2phi[i] =
                1.0 - lambda * lambda *
                          (x[i + 1 + N] - 2.0 * x[i + N] + x[i - 1 + N]) *
                          ((double)(N) * (double)(N));
    }
    d2phi[0] = 1.0 - lambda * lambda * (x[N + 1] - 3.0 * x[N] - 2.0 * t1) *
                         ((double)(N) * (double)(N));
    d2phi[N - 1] = 1.0 - lambda * lambda *
                             (2.0 * t1 - 3.0 * x[2 * N - 1] + x[2 * N - 2]) *
                             ((double)(N) * (double)(N));

    plt::figure();
    plt::figure_size(800, 600);
    plt::named_plot("1 - lambda^2 * d2Phi/dx2", x_axis, d2phi, "k.");
    plt::named_plot("P(x)", x_axis, p, "b-");
    plt::named_plot("Phi(x)", x_axis, phi, "r-");
    plt::xlabel("x");
    plt::ylabel("P and Phi");
    plt::xlim(0.0, 1.0);
    plt::grid(true);
    plt::legend();

    // Save figure
    const char *filename = "figure.png";
    std::cout << "Saving result to " << filename << "...\n";
    plt::save(filename);
#endif

    std::cout << "...done\n\n";

    return 0;
}
