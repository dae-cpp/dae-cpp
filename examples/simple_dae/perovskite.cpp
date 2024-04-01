// g++ -O3 perovskite.cpp -o perovskite.exe -I/home/ivan/workspace/dae-cpp

#include <iostream>
// #include <cmath>
// #include <algorithm>  // for std::max_element

// #include <dae-cpp/Eigen/Dense>
#include <dae-cpp/solver.hpp> // the main header of dae-cpp library solver

const int N0 = 4000;       // Number of points
const double L = 1.0;      // Space interval length
const double lambda = 1.0; // Lambda parameter
const double t1 = 10.0;    // Integration time (0 < t < t1)

// Derived parameters
const double h = L / (double)(N0 - 1); // cell size
const double invh = 1.0 / h;           // inverse cell size

using namespace daecpp;

struct MyMassMatrix : MassMatrix
{
    void operator()(sparse_matrix &M, const double t) const
    {
        M.reserve(2 * N0);
        for (int i = 0; i < N0; ++i)
            M(1.0, i, i);
    }
};

/*
 * RHS of the problem
 * =============================================================================
 */
class MyRHS : public RHS
{
public:
    /*
     * Receives current solution vector x and the current time t.
     * Defines the RHS f.
     */
    void operator()(state_type &f, const state_type &x, const double t) const
    {
        // Locals
        const int N = N0;
        const double invh2 = invh * invh;
        const double invlam2 = 1.0 / (lambda * lambda);

        // RHS for the ion concentration P = dFlux/dx
        for (int i = 1; i < N - 1; i++)
        {
            f[i] = (x[i + 1] - 2.0 * x[i] + x[i - 1] + 0.5 * ((x[i + 1] + x[i]) * (x[N + i + 1] - x[N + i]) - (x[i] + x[i - 1]) * (x[N + i] - x[N + i - 1]))) * invh2;
        }
        f[0] = (x[1] - x[0] + 0.5 * (x[1] + x[0]) * (x[N + 1] - x[N])) * invh2;                                  // Left BC
        f[N - 1] = -(x[N - 1] - x[N - 2] + 0.5 * (x[N - 1] + x[N - 2]) * (x[2 * N - 1] - x[2 * N - 2])) * invh2; // Right BC

        // RHS for the potential Phi
        for (int i = 1; i < N - 1; i++)
        {
            f[i + N] = (x[i + 1 + N] - 2.0 * x[i + N] + x[i - 1 + N]) * invh2 - (1.0 - x[i]) * invlam2;
        }
        f[N] = x[N] + t;                 // Left BC
        f[2 * N - 1] = x[2 * N - 1] - t; // Right BC
    }
};

/*
 * (Optional) Analytical Jacobian in simplified 3-array sparse format
 * =============================================================================
 */

struct MyJacobian : Jacobian
{
    void operator()(sparse_matrix &J, const state_type &x, const double t) const
    {
        // J(1.0, 0, 1);

        // Locals
        const int N = N0;
        const int size = (int)(x.size());
        const double invh2 = invh * invh;
        const double invlam2 = 1.0 / (lambda * lambda);

        J.reserve(6 * size); // Overestimating but it's better than underestimate

        for (int i = 0; i < size; i++)
        {
            if (i == 0)
            {
                J((-1.0 + 0.5 * (x[N + 1] - x[N])) * invh2, i, 0);
                J((1.0 + 0.5 * (x[N + 1] - x[N])) * invh2, i, 1);
                J(-0.5 * (x[0] + x[1]) * invh2, i, N);
                J(0.5 * (x[0] + x[1]) * invh2, i, N + 1);
            }
            else if (i < N - 1)
            {
                J((1.0 - 0.5 * (x[N + i] - x[N + i - 1])) * invh2, i, i - 1);
                J((-2.0 + 0.5 * (x[N + i + 1] - 2.0 * x[N + i] + x[N + i - 1])) * invh2, i, i);
                J((1.0 + 0.5 * (x[N + i + 1] - x[N + i])) * invh2, i, i + 1);
                J(0.5 * (x[i] + x[i - 1]) * invh2, i, N + i - 1);
                J(-0.5 * (x[i + 1] + 2.0 * x[i] + x[i - 1]) * invh2, i, N + i);
                J(0.5 * (x[i + 1] + x[i]) * invh2, i, N + i + 1);
            }
            else if (i == N - 1)
            {
                J((1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2, i, i - 1);
                J((-1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2, i, i);
                J(0.5 * (x[N - 1] + x[N - 2]) * invh2, i, N + i - 1);
                J(-0.5 * (x[N - 1] + x[N - 2]) * invh2, i, N + i);
            }
            else if (i == N)
            {
                J(1.0, i, N);
            }
            else if (i < 2 * N - 1)
            {
                J(invlam2, i, i - N);
                J(invh2, i, i - 1);
                J(-2.0 * invh2, i, i);
                J(invh2, i, i + 1);
            }
            else // i == 2*N-1
            {
                J(1.0, i, 2 * N - 1);
            }
        }
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
    double time = 0.0;
    {
        Timer timer(&time);

        // Solution time 0 <= t <= pi
        double t{10.0};

        // Define the state vector
        state_type x(2 * N0);

        // Initial conditions
        for (int i = 0; i < N0; i++)
        {
            x[i] = 1.0;      // for P - ion concentration
            x[i + N0] = 0.0; // for Phi - potential
        }

        // Set up the RHS of the problem.
        // Class MyRHS inherits abstract RHS class from dae-cpp library.
        MyRHS rhs;

        // Set up the Mass Matrix of the problem.
        // MyMassMatrix inherits abstract MassMatrix class from dae-cpp library.
        MyMassMatrix mass;
        // MassMatrixIdentity mass(2);

        // Create an instance of the solver options and update some of the solver
        // parameters defined in solver_options.h
        SolverOptions opt;

        opt.BDF_order = 4;
        opt.atol = 1e-6;
        opt.rtol = 1e-6;
        // opt.dt_max = 0.05;
        opt.is_mass_matrix_static = true;
        opt.Newton_scheme = 1;
        opt.dt_increase_threshold_delta = 0;

        opt.verbosity = verbosity::extra; // Suppress output to screen (we have our own output
        //                         // defined in Observer function above)

        // We can override Jacobian class from dae-cpp library and provide
        // analytical Jacobian.
        MyJacobian jac;

        // Or we can use numerically estimated Jacobian with the given tolerance.
        // Jacobian jac_est(rhs, 1e-6);

        // Create an instance of the solver with particular RHS, Mass matrix,
        // Jacobian and solver options
        // MySolver solve(rhs, jac, mass);
        System simple_dae(mass, rhs, jac, opt);

        // Now we are ready to solve the set of DAEs
        // std::cout << "Starting DAE solver...\n";
        // std::cout << "time\tx\ty\terror1\terror2\n";

        // Solve the system

        // t_output will be move-constructed if the user provides an
        // r-value (either directly from a temporary, or by moving from an lvalue)

        // Keep a copy:
        // std::vector<string> items = { "1", "2", "3" };
        // Test t;
        // t.someFunction(items); // pass items directly - we keep a copy

        // Don't keep a copy:
        // std::vector<string> items = { "1", "2", "3" };
        // Test t;
        // t.someFunction(std::move(items)); // move items - we don't keep a copy
        // or t.someFunction({ "1", "2", "3" });

        // std::vector<double> t_out{1, 2, 3, 4, 5, 2, 5, 0, -5};
        int status = simple_dae.solve(x, t);
        // int status = simple_dae.solve(x, t, {1, 2, 3, 4, 5});
        // int status = simple_dae.solve(x, t, std::move(t_out));
        // std::cout << "t_out size: " << t_out.size() << '\n';

        // using Eigen::MatrixXd;
        std::cout << "x = " << x[0] << " " << x[N0] << " " << x[N0 + (N0 - 1) / 5 * 3] * 0.6 + x[N0 + (N0 - 1) / 5 * 3 + 1] * 0.4 << " " << x[2 * N0 - 1] << " ";
        // std::cout << "Time: " << t << "\t" << x[0] << "\t" << x[1] << "\t" << std::exp(-10.0) << "\t" << -std::exp(-10.0) << '\n';

        // MatrixXd m(2, 2);
        // m(0, 0) = 3;
        // m(1, 0) = 2.5;
        // m(0, 1) = -1;
        // m(1, 1) = m(1, 0) + m(0, 1);
        // std::cout << m << std::endl;

        // Check errors
        //     double max_err1 = *std::max_element(solve.err1.begin(), solve.err1.end());
        //     double max_err2 = *std::max_element(solve.err2.begin(), solve.err2.end());

        //     std::cout << "\nMaximum absoulte error (1) x*x + y*y = 1: " << max_err1
        //               << '\n';
        //     std::cout << "Maximum absolute error (2) x(t) - sin(t) = 0 for t <= pi/2 "
        //                  "or x(t) = 1 for t > pi/2: "
        //               << max_err2 << '\n';

        // #ifdef DAE_SINGLE
        //     const bool check_result = (max_err1 > 1e-6 || max_err2 > 1e-6 || status);
        // #else
        //     const bool check_result = (max_err1 > 1e-15 || max_err2 > 1e-6 || status);
        // #endif

        //     if(check_result)
        //         std::cout << "...Test FAILED\n\n";
        //     else
    }

    std::cout << "...done. Time = " << time << "\n\n";

    return 0;
}
