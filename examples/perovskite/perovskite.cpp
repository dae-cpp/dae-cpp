/*
 * TODO: Description of the problem
 * Keywords: perovskite solar cell, potential, ion concentration,
 * singular mass matrix
 *
 */

#include <iostream>
#include <chrono>
#include <cmath>

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

void solution_check(dae::state_type &x);

int main()
{
    // These parameters can be obtained from a parameter file or
    // as command line options. Here for simplicity we define them as constants.
    const int    N      = 4000;  // Number of cells
    const double L      = 1.0;   // Space interval length
    const double lambda = 1.0;   // Lambda parameter
    const double t1     = 10.0;  // Integration time (0 < t < t1)

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
    std::cout << "\nStarting DAE solver...\n";

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

    // Compare solution with an alternative solver (e.g. MATLAB)
    solution_check(x1);

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

    // Compare solution
    solution_check(x2);

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

void solution_check(dae::state_type &x)
{
    std::cout << "Solution check:\n";

    const int N     = (int)(x.size()) / 2;
    const int N_sol = 9;

    double sol[N_sol];

    // MATLAB ode15s solution (Finite Elements), N = 4001 points.
    // Note that Finite Volume method has to extrapolate values on the
    // boundaries and interpolate solution to the nodes, that's why comparison
    // FE - FV is not 100% accurate. Given here for reference.
    const double ode15s_MATLAB[N_sol] = {19.9949,    2.72523, 0.382148,
                                         0.00101573, -10.0,   -6.04056,
                                         -2.08970,   1.90021, 5.93011};
    // clang-format off
    sol[0] = x[0];                                           // P(x = 0)
    sol[1] = x[(N-1)/10] * 0.1 + x[(N-1)/10+1] * 0.9;        // P(x = 0.1)
    sol[2] = x[(N-1)/5] * 0.2 + x[(N-1)/5+1] * 0.8;          // P(x = 0.2)
    sol[3] = (x[N/2-1] + x[N/2]) / 2.0;                      // P(x = 0.5)
    sol[4] = x[N];                                           // Phi(x = 0)
    sol[5] = x[N+(N-1)/5*1] * 0.2 + x[N+(N-1)/5*1+1] * 0.8;  // Phi(x = 0.2)
    sol[6] = x[N+(N-1)/5*2] * 0.4 + x[N+(N-1)/5*2+1] * 0.6;  // Phi(x = 0.4)
    sol[7] = x[N+(N-1)/5*3] * 0.6 + x[N+(N-1)/5*3+1] * 0.4;  // Phi(x = 0.6)
    sol[8] = x[N+(N-1)/5*4] * 0.8 + x[N+(N-1)/5*4+1] * 0.2;  // Phi(x = 0.8)
    // clang-format on

    std::cout << "  MATLAB ode15s\t<->  dae-cpp\t(rel. error)\n";

    double err_max = 0;

    for(int i = 0; i < N_sol; i++)
    {
        double error = (sol[i] - ode15s_MATLAB[i]) / ode15s_MATLAB[i] * 100.0;
        if(std::fabs(error) > err_max)
            err_max = std::fabs(error);
        std::cout << "     " << ode15s_MATLAB[i] << "\t<->  " << sol[i] << "\t("
                  << error << "%)\n";
    }
    std::cout << "Maximum relative error: " << err_max << "%\n";
}
