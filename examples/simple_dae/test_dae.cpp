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

// g++ -O3 -std=c++17 -Wall test_dae.cpp -o test_dae.exe -I/home/ivan/workspace/dae-cpp

#include <iostream>
// #include <cmath>
// #include <algorithm>  // for std::max_element

// #include <dae-cpp/Eigen/Dense>
#include <dae-cpp/solver.hpp> // the main header of dae-cpp library solver

using namespace daecpp;

struct MyMassMatrix
{
    void operator()(sparse_matrix &M, const double t) const // const is optional
    {
        M(0, 0, 1.0);
        // M(0, 0, 1.0);
        // M(1, 1, 1.0);
    }
};

/*
 * RHS of the problem
 * =============================================================================
 */
class MyRHS // : public RHS
{
public:
    /*
     * Receives current solution vector x and the current time t.
     * Defines the RHS f.
     */
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];
        f[1] = cos(t) - x[1];
        // f[0] = sin(t) - x[1];
        // f[1] = sin(t) - x[2];
        // f[2] = cos(t) - x[0];
    }
};

/*
 * (Optional) Analytical Jacobian in simplified 3-array sparse format
 * =============================================================================
 */

struct MyJacobian //: JacobianMatrix
{
    void operator()(sparse_matrix &J, const state_vector &x, const double t) // const
    {
        J.reserve(3);
        J(0, 1, 1.0);
        J(1, 1, -1.0);
        // J(0, 1, -1.0);
        // J(1, 2, -1.0);
        // J(2, 0, -1.0);
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
        double t{100.0};

        // Define the state vector
        state_vector x(2);

        // Initial conditions
        x[0] = 0;
        x[1] = 1;
        // x[0] = 1;
        // x[1] = 0;
        // x[2] = -2;

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

        // opt.BDF_order = 1;
        opt.atol = 1e-12;
        opt.rtol = 1e-12;
        opt.Newton_scheme = 1;
        // opt.dt_init = 0.001;
        // opt.dt_max = 1024;

        // opt.verbosity = verbosity::extra; // Suppress output to screen (we have our own output
        //                         // defined in Observer function above)

        // We can override Jacobian class from dae-cpp library and provide
        // analytical Jacobian.
        MyJacobian jac;

        // Or we can use numerically estimated Jacobian with the given tolerance.
        // Jacobian jac_est(rhs, 1e-6);

        // Create an instance of the solver with particular RHS, Mass matrix,
        // Jacobian and solver options
        // MySolver solve(rhs, jac, mass);
        // System test_dae(mass, rhs, jac, opt);
        // System test_dae(mass, rhs, opt);

        // daecpp::solve(mass, rhs, jac, x, t, opt);
        // daecpp::solve(mass, rhs, {1.0, -1.0}, 1.0, {0.5}, opt);
        // daecpp::solve(MyMassMatrix(), MyRHS(), x, 10.0, opt);

        // daecpp::solve(MyMassMatrix(), MyRHS(), jac, {0, -1}, 3.14, opt);

        // std::vector<double> t_list = {77, 166, 499, 999, 1450, 1990, 2222, 2888, 3333, 3850, 4444, 4930, 5555, 6060, 6500, 7070, 7777, 8333, 9090, 9500, 9090, 10000};
        // std::vector<double> t_list;
        // double val = 0.12345;
        // while (val < t)
        // {
        //     t_list.push_back(val);
        //     val += 34.5678; //12.3567;
        // }
        // t_list.push_back(t);

        System sys(MyMassMatrix(), rhs);

        sys.opt = opt;

        // double exact = sin(t) - 2.0 * cos(t); // z
        // double exact = 2.0 * sin(t); // y
        double exact = sin(t); //x

        for (int order = 1; order <= 4; order++)
        {

            sys.opt.BDF_order = order;
            sys.opt.dt_max = 1.0;

            sys.solve(x, t, jac);
        
            double first = abs(sys.sol.x.back()[0] - exact);
            std::cout << std::setprecision(14) << first << '\n';
            // std::cout << "Status: " << sys.status << '\n'
            //           << '\n'
            //           << '\n';

            // while (val < t)
            // {
            //     t_list.push_back(val);
            //     val += 3.45678; // 12.3567;
            // }

            // for (int i = 0; i < 6; ++i)
            // {
            sys.opt.dt_max /= 2;
            sys.solve(x, t, jac);
            double second = abs(sys.sol.x.back()[0] - exact);
            std::cout << second << '\n';
            // std::cout << "Status: " << sys.status << '\n'
            //           << '\n'
            //           << '\n';
            // }
            std::cout << "======== " << log2(second / first) << '\n';
        }

        // sys.sol.print();

        // std::cout << jac << '\n';
        // std::cout <<<< '\n';

        // daecpp::solve(MassMatrixZero(), MyRHS(), {0, 0}, opt.dt_init, opt);

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
        // int status = test_dae.solve(x, t);
        // int status = test_dae.solve(x, t, {1, 2, 3, 4, 5});
        // int status = test_dae.solve(x, t, std::move(t_out));
        // std::cout << "t_out size: " << t_out.size() << '\n';

        // using Eigen::MatrixXd;
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
