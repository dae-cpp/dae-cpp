/*
 * This example solves the system of ODEs that describes diffusion in 2D plane:
 *
 * dC/dt = D*(d/dx(dC/dx) + d/dy(dC/dy)),
 *
 * where C is the concentration (dimensionless) on the square unit domain,
 * 0 <= x <= 1 and 0 <= y <= 1. D is the diffusion coefficient.
 *
 * Initial condition is an instantaneous point source in two dimensions:
 * C(x,y,t=0) = delta_function(x-1/2,y-1/2).
 * Boundary conditions: dC/dx = dC/dy = 0 on the boundaries.
 *
 * The system will be resolved using Finite Volume approach in the time
 * interval 0 <= t <= 0.01, and compared with the analytical solution.
 * Since the scheme is conservative and there are no sources in the domain,
 * the total concentration in the domain should be constant and equal to 1.
 *
 * Keywords: diffusion equation, 2D, finite volume method.
 */

#include <iostream>
#include <chrono>
#include <cmath>

#include "../../src/solver.h"  // the main header of dae-cpp library solver

#include "diffusion_2d_RHS.h"  // RHS of the problem

namespace dae = daecpp;  // A shortcut to dae-cpp library namespace

// python3 + numpy + matplotlib should be installed in order to enable plotting
// #define PLOTTING

#ifdef PLOTTING
#include "../../src/external/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

// To compare dae-cpp solution with the analytical solution
int solution_check(dae::state_type &x, MKL_INT N, double t, double D);

/*
 * MAIN FUNCTION
 * =============================================================================
 * Returns '0' if solution comparison is OK or '1' if solution error is above
 * acceptable tolerance.
 */
int main()
{
    // These parameters can be obtained from a parameter file or as command line
    // options. Here for simplicity we define them as constants.
    const MKL_INT N  = 51;    // Number of cells along one axis
    const double  D  = 1.0;   // Diffusion coefficient
    const double  t1 = 0.01;  // Integration time (0 < t < t1)

    std::cout << "N = " << N << "; D = " << D << "; t = " << t1 << '\n';

    using clock     = std::chrono::high_resolution_clock;
    using time_unit = std::chrono::milliseconds;

    // Define state vectors. Here N*N is the total number of the equations.
    dae::state_type x(N * N);

    // Initial conditions
    for(MKL_INT i = 0; i < N * N; i++)
    {
        x[i] = 0.0;
    }
    x[N * N / 2] = N * N;  // 1/(h*h) -- numerical delta-function

    // Set up the RHS of the problem.
    // Class MyRHS inherits abstract RHS class from dae-cpp library.
    MyRHS rhs(N, D);

    // Set up the Mass Matrix of the problem. In this case this matrix is
    // identity, so we can use a helper class provided by dae-cpp library.
    dae::MassMatrixIdentity mass(N * N);

    // Create an instance of the solver options and update some of the solver
    // parameters defined in solver_options.h
    dae::SolverOptions opt;

    opt.dt_init         = 5.0e-5;  // Change initial time step
    opt.fact_every_iter = false;   // Gain some speed (delay the update
                                   // of Jacobian and the matrix factorisation)

    // We can override Jacobian class from dae-cpp library and provide
    // analytical Jacobian. But we will use numerically estimated one.
    dae::Jacobian jac_est(rhs);

    // Create an instance of the solver with particular RHS, Mass matrix,
    // Jacobian and solver options
    dae::Solver solve(rhs, jac_est, mass, opt);

    // Now we are ready to solve the set of DAEs
    std::cout << "\nStarting DAE solver...\n";

    {
        auto tic0 = clock::now();
        solve(x, t1 / 4);  // This line can be removed. It is given here just as
                           // an example. Here we produce an intermediate
                           // solution at time t = (t1 / 4). This solution
                           // will be stored in the vector x. Note that a better
                           // way to get intermediate results is to override
                           // observer function from daecpp::Solver class.
        solve(x, t1);      // Reuse vector x as an initial condition and get
                           // the final solution at time t = t1.
        auto tic1 = clock::now();

        std::cout
            << "Solver execution time: "
            << std::chrono::duration_cast<time_unit>(tic1 - tic0).count() /
                   1000.0
            << " sec." << '\n';
    }

    // Compare result with the analytical solution
    int check_result = solution_check(x, N, t1, D);

    // Plot the solution
#ifdef PLOTTING
    const double h = 1.0 / (double)N;

    dae::state_type_matrix x_axis, y_axis, z_axis;

    for(MKL_INT i = 0; i < N; i++)
    {
        dae::state_type x_row, y_row, z_row;

        for(MKL_INT j = 0; j < N; j++)
        {
            x_row.push_back((double)j * h + h * 0.5);
            y_row.push_back((double)i * h + h * 0.5);
            z_row.push_back(x[j + i * N]);
        }

        x_axis.push_back(x_row);
        y_axis.push_back(y_row);
        z_axis.push_back(z_row);
    }

    plt::figure();
    plt::figure_size(800, 600);
    plt::plot_surface(x_axis, y_axis, z_axis);

    // Save figure
    const char *filename = "diffusion_2d.png";
    std::cout << "Saving result to " << filename << "...\n";
    plt::save(filename);
#endif

    if(check_result)
        std::cout << "...Test FAILED\n\n";
    else
        std::cout << "...done\n\n";

    return check_result;
}

/*
 * Returns analytical solution
 */
double analyt(double x, double y, double t, double D)
{
    double Dt4 = 1.0 / (D * t * 4.0);
    return Dt4 / 3.1415926535897932384626433832795 *
           exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) * Dt4);
}

/*
 * Returns '0' if solution comparison is OK or '1' if the error is above
 * acceptable tolerance
 */
int solution_check(dae::state_type &x, MKL_INT N, double t, double D)
{
    std::cout << "Solution check:\n";

    const double h = 1.0 / (double)N;

    double total_C = 0;
    double err_max = 0;

    for(MKL_INT i = 0; i < N; i++)
    {
        for(MKL_INT j = 0; j < N; j++)
        {
            MKL_INT ind = j + i * N;

            total_C += x[ind];

            double xi = (double)j * h + h * 0.5;
            double yi = (double)i * h + h * 0.5;
            double an = analyt(xi, yi, t, D);

            double error;

            if(an > 1.0)
            {
                error = (x[ind] - an) / an * 100.0;  // relative error

                if(fabs(error) > err_max)
                {
                    err_max = fabs(error);
                }
            }
            else
            {
                // error = (x[ind] - an);  // absolute error
            }
        }
    }

    total_C *= h * h;

    double err_conc = fabs(total_C - 1.0) * 100;

    std::cout << "Total concentration:    " << total_C << " (" << err_conc
              << "% deviation from the analytical value)\n";
    std::cout << "Maximum relative error: " << err_max << "%\n";

#ifdef DAE_SINGLE
    if(err_max < 1.0 && err_conc < 2.0e-5)
#else
    if(err_max < 1.0 && err_conc < 1.0e-10)
#endif
        return 0;
    else
        return 1;
}
