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

// g++ -O3 -std=c++17 -Wall simple_dae.cpp -o simple_dae.exe -I/home/ivan/workspace/dae-cpp

#include <iostream>
// #include <cmath>
// #include <algorithm>  // for std::max_element

// #include <dae-cpp/Eigen/Dense>
#include <dae-cpp/solver.hpp> // the main header of dae-cpp library solver

using namespace daecpp;

/*
 * Singular mass matrix in sparse format:
 *
 * M = |1 0|
 *     |0 0|
 */
struct MyMassMatrix
{
    /*
     * Defines the mass matrix `M` of the DAE system `M dx/dt = f`.
     * The mass matrix should be defined in sparse format (non-zero elements only) and can depend on time `t`.
     * Matrix `M` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &M, const double t)
    {
        M(0, 0, 1.0); // Row 0, column 0, non-zero element 1.0
    }
};

/*
 * The vector-function (RHS) of the problem:
 *
 * f(x,y,t) = | y
 *            | cos(t) - y
 */
struct MyRHS
{
    /*
     * Defines the RHS (vector function) `f` of the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the RHS vector `f`.
     * Vector `f` is already pre-allocated with `f.size() == x.size()`.
     *
     * For the given DAE system,
     * | x = x[0]
     * | y = x[1]
     */
    void operator()(state_type &f, const state_type &x, const double t)
    {
        f[0] = x[1];          // y
        f[1] = cos(t) - x[1]; // cos(t) - y
    }
};

/*
 * (OPTIONAL) Analytic Jacobian in sparse format.
 * The DAE solver will use automatic (algorithmic) differentiation if analytic Jacobian is not provided.
 * However, providing the analytic Jacobian can significantly speed up the computation for large systems.
 *
 * Differentiating the RHS w.r.t. x[0] and x[1] gives the following Jacobian matrix:
 *
 * J = |0  1|
 *     |0 -1|
 */
struct MyJacobian
{
    /*
     * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the Jacobian matrix `J`.
     * Matrix `J` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        J.reserve(2);  // Pre-allocates memory for 2 non-zero elements (optional)
        J(0, 1, 1.0);  // Row 0, column 1, non-zero element 1.0
        J(1, 1, -1.0); // Row 1, column 1, non-zero element -1.0
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

        MyMassMatrix mass; // Mass matrix object
        MyRHS rhs;         // The vector-function object

        System my_system(mass, rhs); // Defines the DAE system object

        state_vector x0{0, 1}; // Initial condition: x = 0, y = 1
        double t_end{10.0};     // Solution interval: t = [0, t_end]

    std::vector<double> t_list; // List of output times that breaks the time step uniformity

    double t{0.123};
    while (t < t_end)
    {
        t_list.push_back(t);
        t += 0.123;
    }
    t_list.push_back(t_end);

        my_system.opt.verbosity = verbosity::extra;
        my_system.opt.atol = 1e-10;
        my_system.opt.rtol = 1e-10;
        my_system.opt.dt_init = 0.000001;
        my_system.opt.variability_threshold_high = 1.0;
        my_system.opt.variability_threshold_low = 1.0;

        // my_system.solve(x0, t_end); // Solves the DAE system `my_system` with the given initial condition `x0` and time `t_end`

        // std::cout << "Abs. error: " << my_system.sol.x.back()[0] - sin(t_end) << '\t' << my_system.sol.x.back()[1] - cos(t_end) << '\n';

        double exact = sin(t_end); // x

        for (int order = 1; order <= 4; order++)
        {
            std::cout << "ORDER: " << order << '\n';

            my_system.sol.x.clear();
            my_system.sol.t.clear();

            my_system.opt.BDF_order = order;
            // my_system.opt.dt_init = 0.1;
            my_system.opt.dt_max = 0.1;

            my_system.solve(x0, t_list);

            // norm 1
            double norm1{0.0};
            for (std::size_t i = 0; i < my_system.sol.t.size(); ++i)
            {
                double err = abs(my_system.sol.x[i][0] - sin(my_system.sol.t[i]));
                if (err > norm1)
                {
                    norm1 = err;
                }
            }
            std::cout << "First: " << norm1 << '\n';

            break;
            // my_system.opt.dt_init /= 2;
            my_system.opt.dt_max /= 10;

            my_system.sol.x.clear();
            my_system.sol.t.clear();

            my_system.solve(x0, t_list);

            // norm 1
            double norm2{0.0};
            for (std::size_t i = 0; i < my_system.sol.t.size(); ++i)
            {
                double err = abs(my_system.sol.x[i][0] - sin(my_system.sol.t[i]));
                if (err > norm2)
                {
                    norm2 = err;
                }
            }
            std::cout << "second: " << norm2 << '\n';

            std::cout << "======== " << log10(norm1 / norm2) << '\n';
        }

        // Or we can use numerically estimated Jacobian with the given tolerance.
        // Jacobian jac_est(rhs, 1e-6);

        // Create an instance of the solver with particular RHS, Mass matrix,
        // Jacobian and solver options
        // MySolver solve(rhs, jac, mass);
        // System simple_dae(mass, rhs, jac, opt);
        // System simple_dae(mass, rhs, opt);

        // daecpp::solve(mass, rhs, jac, x, t, opt);
        // daecpp::solve(mass, rhs, {1.0, -1.0}, 1.0, {0.5}, opt);
        // daecpp::solve(MyMassMatrix(), MyRHS(), x, 10.0, opt);

        // daecpp::solve(MyMassMatrix(), MyRHS(), jac, {0, -1}, 3.14, opt);

        // System sys((MyMassMatrix()), (MyRHS()));

        // sys.opt = opt;

        // sys.solve({0, 1}, 3.14, jac);

        // sys.sol.print();

        // std::cout << "Status: " << sys.status << '\n';
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
        // int status = simple_dae.solve(x, t);
        // int status = simple_dae.solve(x, t, {1, 2, 3, 4, 5});
        // int status = simple_dae.solve(x, t, std::move(t_out));
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
