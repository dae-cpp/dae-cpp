/*
 * TODO: Description
 *
 * Keywords: Robertson problem, stiff DAE.
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

// To compare dae-cpp solution with the analytical solution
// int solution_check(state_type &x, MKL_INT N, double t, double D);

class MyMassMatrix : public MassMatrix
{
public:
    void operator()(daecpp::sparse_matrix_holder &M)
    {
        M.A.resize(3);
        M.ja.resize(3);
        M.ia.resize(4);

        M.A[0] = 1;
        M.A[1] = 1;
        M.A[2] = 0;

        M.ja[0] = 0;
        M.ja[1] = 1;
        M.ja[2] = 2;

        M.ia[0] = 0;
        M.ia[1] = 1;
        M.ia[2] = 2;
        M.ia[3] = 3;
    }
};

/*
 * RHS of the problem
 * =============================================================================
 */
class MyRHS : public RHS
{
public:
    void operator()(const daecpp::state_type &x, daecpp::state_type &f,
                    const double t)
    {
        f[0] = -0.04 * x[0] + 1.0e4 * x[1] * x[2];
        f[1] = 0.04 * x[0] - 1.0e4 * x[1] * x[2] - 3.0e7 * x[1] * x[1];
        f[2] = x[0] + x[1] + x[2] - 1;
    }
};

class MySolver : public Solver
{
public:
    MySolver(daecpp::RHS &rhs, daecpp::Jacobian &jac, daecpp::MassMatrix &mass,
             daecpp::SolverOptions &opt)
        : daecpp::Solver(rhs, jac, mass, opt)
    {
    }

    /*
     * Overloaded observer.
     * Receives current solution vector and the current time every time step.
     * Prints current time t and potential phi on the right boundary.
     */
    void observer(daecpp::state_type &x, const double t)
    {
        std::cout << " | " << x[0] << ' ' << 1e4 * x[1] << ' ' << x[2]
                  << " == " << x[0] + x[1] + x[2] - 1.0;
    }
};

class MyJacobian : public Jacobian
{
public:
    MyJacobian(daecpp::RHS &rhs) : daecpp::Jacobian(rhs) {}

    void operator()(daecpp::sparse_matrix_holder &J,
                    const daecpp::state_type &x, const double t)
    {
        J.A.resize(9);
        J.ja.resize(9);
        J.ia.resize(4);

        J.A[0] = -0.04;
        J.A[1] = 1.0e4 * x[2];
        J.A[2] = 1.0e4 * x[1];
        J.A[3] = 0.04;
        J.A[4] = -1.0e4 * x[2] - 6.0e7 * x[1];
        J.A[5] = -1.0e4 * x[1];
        J.A[6] = 1.0;
        J.A[7] = 1.0;
        J.A[8] = 1.0;

        J.ja[0] = 0;
        J.ja[1] = 1;
        J.ja[2] = 2;
        J.ja[3] = 0;
        J.ja[4] = 1;
        J.ja[5] = 2;
        J.ja[6] = 0;
        J.ja[7] = 1;
        J.ja[8] = 2;

        J.ia[0] = 0;
        J.ia[1] = 3;
        J.ia[2] = 6;
        J.ia[3] = 9;
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 * Returns '0' if solution comparison is OK or '1' if solution error is above
 * acceptable tolerance.
 */
int main()
{
    const double t1 = 4.0e6;

    // Define state vector
    state_type x(3);

    // Initial conditions
    // Use inconsistent initial condition to test initialization
    x[0] = 1;
    x[1] = 0;
    x[2] = 1e-3;  // Should be 0

    // Set up the RHS of the problem.

    // Class MyRHS inherits abstract RHS class from dae-cpp library.
    MyRHS rhs;

    // Set up the Mass Matrix of the problem. In this case this matrix is
    // identity, so we can use a helper class provided by dae-cpp library.
    MyMassMatrix mass;

    // Create an instance of the solver options and update some of the solver
    // parameters defined in solver_options.h
    SolverOptions opt;

    opt.dt_init = 1.0e-6;  // Change initial time step
    // opt.fact_every_iter = false;   // Gain some speed (delay the update
    // of Jacobian and the matrix factorisation)

    opt.verbosity             = 2;
    opt.dt_max                = t1 / 100;
    opt.time_stepping         = 1;
    opt.dt_increase_threshold = 2;
    // opt.dt_decrease_threshold = 6;
    // opt.atol = 1e-7;
    // opt.bdf_order = 6;

    // We can override Jacobian class from dae-cpp library and provide
    // analytical Jacobian. But we will use numerically estimated one.
    Jacobian jac_est(rhs, 1e-10);
    jac_est.print(x, 0);

    MyJacobian jac(rhs);
    jac.print(x, 0);

    // Create an instance of the solver with particular RHS, Mass matrix,
    // Jacobian and solver options
    MySolver solve(rhs, jac, mass, opt);

    // Now we are ready to solve the set of DAEs
    std::cout << "\nStarting DAE solver...\n";

    solve(x, t1);

    std::cout << " | " << x[0] << ' ' << 1e4 * x[1] << ' ' << x[2]
              << " == " << x[0] + x[1] + x[2] << '\n';

    // Compare result with the analytical solution
    const double x_ref[3]     = {0.00051675, 2.068e-9, 0.99948324};
    const double conservation = std::abs(x[0] + x[1] + x[2] - 1);
    double       result       = 0.0;
    for(int i = 0; i < 3; i++)
        result += std::abs(x[i] - x_ref[i]) / x_ref[i] * 100;

    std::cout << result << "% " << conservation << '\n';
    //    int check_result = solution_check(x, N, t1, D);

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

    const bool check_result = (result > 1.0 || conservation > 1e-14);

    if(check_result)
        std::cout << "...Test FAILED\n\n";
    else
        std::cout << "...done\n\n";

    return check_result;
}

/*
 * Returns '0' if solution comparison is OK or '1' if the error is above
 * acceptable tolerance
 */
/*
int solution_check(state_type &x, MKL_INT N, double t, double D)
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

            if(an > 1.0)
            {
                double error = (x[ind] - an) / an * 100.0;  // relative error

                if(std::abs(error) > err_max)
                {
                    err_max = std::abs(error);
                }
            }
        }
    }

    total_C *= h * h;

    double err_conc = std::abs(total_C - 1.0) * 100;

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
*/