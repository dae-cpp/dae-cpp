/*
 * Solves a very simple system of differential algebraic equation as a test:
 *
 * x' = y
 * y  = -x  // 0  = x*x + y*y - 1
 *
 * Initial conditions are: x = 0, y = 1 for t = 0.
 *
 * The solution of this system is
 *
 * x = sin(t), y = cos(t), 0 <= t <= pi/2;
 * x = 1, y = 0, t > pi/2.
 *
 * Each time step we will check that
 * (1) x*x + y*y = 1 for any t, and
 * (2) x(t) = sin(t) for t <= pi/2, x(t) = 1 for t > pi/2
 *
 * with the absolute tolerance at least 1e-6.
 */

// g++ simple_dae.cpp -o simple_dae.exe -I/home/ivan/workspace/dae-cpp

#include <iostream>
// #include <cmath>
// #include <algorithm>  // for std::max_element

// #include <dae-cpp/Eigen/Dense>
#include <dae-cpp/solver.hpp> // the main header of dae-cpp library solver

using namespace daecpp;

/*
 * Singular mass matrix in simplified 3-array sparse format
 * =============================================================================
 * The matrix has the following form:
 * M = |1 0|
 *     |0 0|
 */
// class MyMassMatrix : public MassMatrix
// {
// public:
//     void operator()(sparse_matrix &M) const
//     {
//         // M.A.resize(1); // Number of non-zero elements
//         // M.i.resize(1); // Number of non-zero elements
//         // M.j.resize(1); // Number of non-zero elements

//         // // Non-zero elements
//         // M.A[0] = 1;
//         // // M.A[1] = 0;

//         // // Column index of each element given above
//         // M.j[0] = 0;
//         // // M.j[1] = 1;

//         // // Row index of each element in M.A:
//         // M.i[0] = 0;
//         // // M.i[1] = 1;
//         M.reserve(1);
//         M.add_element(1.0, 0, 0);
//     }
// };

struct MyMassMatrix : MassMatrix
{
    void operator()(sparse_matrix &M, const double t) const
    {
        M.reserve(1);
        M(1.0, 0, 0);
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
        f[0] = x[1];
        f[1] = x[0] + x[1]; // x[0] * x[0] + x[1] * x[1] - 1.0;
    }
};

// /*
//  * (Optional) Observer
//  * =============================================================================
//  * Every time step checks that
//  * (1) x*x + y*y = 1, and
//  * (2) x(t) - sin(t) = 0 for t <= pi/2, x(t) = 1 for t > pi/2
//  * and prints solution and errors to console.
//  */
// class MySolver : public Solver
// {
// public:
//     MySolver(RHS &rhs, Jacobian &jac, MassMatrix &mass, SolverOptions &opt)
//         : Solver(rhs, jac, mass, opt)
//     {
//     }

// #ifdef PLOTTING
//     state_type x_axis, x0, x1;  // For plotting
// #endif
//     state_type err1, err2;  // To check errors

//     /*
//      * Overloaded observer.
//      * Receives current solution vector and the current time every time step.
//      */
//     void observer(state_type &x, const double t)
//     {
//         double e1 = std::abs(x[0] * x[0] + x[1] * x[1] - 1.0);
//         double e2 = 0;

//         if(t <= 1.5707963)
//             e2 = std::abs(std::sin(t) - x[0]);
//         else
//             e2 = std::abs(x[0] - 1.0);

//         std::cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << e1 << '\t'
//                   << e2 << '\n';

//         err1.push_back(e1);
//         err2.push_back(e2);

// #ifdef PLOTTING
//         // Save data for plotting
//         x_axis.push_back(t);
//         x0.push_back(x[0]);
//         x1.push_back(x[1]);
// #endif
//     }
// };

/*
 * (Optional) Analytical Jacobian in simplified 3-array sparse format
 * =============================================================================
 */

struct MyJacobian : Jacobian
{
    void operator()(sparse_matrix &J, const state_type &x, const double t) const
    {
        J.reserve(3);
        J(1.0, 0, 1);
        J(1.0, 1, 0);
        J(1.0, 1, 1);
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
    double time=0.0;
    {
        Timer timer(&time);

    // Solution time 0 <= t <= pi
    double t{1.0};

    // Define the state vector
    state_type x(2);

    // Initial conditions
    x[0] = 1;
    x[1] = -1;

    // Set up the RHS of the problem.
    // Class MyRHS inherits abstract RHS class from dae-cpp library.
    MyRHS rhs;

    // Set up the Mass Matrix of the problem.
    // MyMassMatrix inherits abstract MassMatrix class from dae-cpp library.
    MyMassMatrix mass;
    // MassMatrixIdentity mass(2);

    // Create an instance of the solver options and update some of the solver
    // parameters defined in solver_options.h
    // SolverOptions opt;

    // opt.dt_init = 1.0e-2;   // Change the initial time step.
    //                         // It should be relatively small, because the first
    //                         // step in time is first order accuracy.
    //                         // Reducing dt_init decreases the error (2)
    // opt.time_stepping = 1;  // Use simple stability-based adaptive time stepping
    //                         // algorithm.
    // opt.bdf_order = 6;      // Use BDF-6
    // opt.verbosity = 0;      // Suppress output to screen (we have our own output
    //                         // defined in Observer function above)

    // We can override Jacobian class from dae-cpp library and provide
    // analytical Jacobian.
    MyJacobian jac;

    // Or we can use numerically estimated Jacobian with the given tolerance.
    // Jacobian jac_est(rhs, 1e-6);

    // Create an instance of the solver with particular RHS, Mass matrix,
    // Jacobian and solver options
    // MySolver solve(rhs, jac, mass);
    Solver solve(mass, rhs, jac);

    // Now we are ready to solve the set of DAEs
    std::cout << "Starting DAE solver...\n";
    // std::cout << "time\tx\ty\terror1\terror2\n";

    // Solve the system
    int status = solve(x, t);

    // using Eigen::MatrixXd;

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
