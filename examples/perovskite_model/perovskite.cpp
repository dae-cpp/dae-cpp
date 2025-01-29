/*
 * This example solves a big DAE system that describes potential distribution
 * and ion concentration in a perovskite solar cell:
 *
 * | dP/dt = d/dx(dP/dx + P * dPhi/dx),
 * | d^2(Phi)/dx^2 = (1 - P)/lambda^2,
 *
 * where P is the ion concentration (dimensionless),
 * Phi is the potential (dimensionless) along axis x, 0 <= x <= 1.
 * Lambda parameter is a constant.
 *
 * Initial condition (t = 0):
 * | P = 1
 * | Phi = 0
 *
 * Boundary conditions:
 * | (dP/dx + P * dPhi/dx) = 0 for x = 0 and x = 1,
 * | Phi(t,x=0) = -t,
 * | Phi(t,x=1) = t.
 *
 * The system will be solved using Finite Difference approach in the time interval t = [0, 10].
 *
 * This example introduces Solution Manager, that will work as solution observer and a simple event function.
 * It also shows how parameters, such as the number of discretization points, can be passed to the mass matrix,
 * vector function, and Jacobian matrix.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

// Main dae-cpp header
#include <dae-cpp/solver.hpp>

// dae-cpp namespace
using namespace daecpp;

/*
 * DAE system parameters
 */
struct MyParams
{
    int_type N{4000};   // Number of discretization points
    double L{1.0};      // Length of the domain
    double lambda{1.0}; // Lambda parameter
};

/*
 * Singular mass matrix in sparse format
 */
class MyMassMatrix
{
    MyParams p; // Parameters

public:
    explicit MyMassMatrix(MyParams &params) : p(params) {}

    /*
     * Defines the mass matrix `M` of the DAE system `M dx/dt = f`.
     * The mass matrix should be defined in sparse format (non-zero elements only) and can depend on time `t`.
     * Matrix `M` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &M, const double t)
    {
        M.reserve(p.N); // Reserves memory for `N` non-zero elements

        for (int_type i = 0; i < p.N; ++i)
        {
            M(i, i, 1.0);
        }
    }
};

/*
 * The vector-function (RHS) of the problem
 */
class MyRHS
{
    MyParams p; // Parameters

public:
    explicit MyRHS(MyParams &params) : p(params) {}

    /*
     * Defines the RHS (vector function) `f` of the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the RHS vector `f`.
     * Vector `f` is already pre-allocated with `f.size() == x.size()`.
     */
    void operator()(state_type &f, const state_type &x, const double t)
    {
        // Derived parameters
        const double h = p.L / (double)(p.N - 1);           // Cell size
        const double invh2 = 1.0 / (h * h);                 // Inverse cell size squared
        const double invlam2 = 1.0 / (p.lambda * p.lambda); // Inverse lambda squared

        // An alias for p.N
        const int_type N = p.N;

        // RHS for the ion concentration P = dFlux/dx
        for (int_type i = 1; i < N - 1; ++i)
        {
            f[i] = (x[i + 1] - 2.0 * x[i] + x[i - 1] + 0.5 * ((x[i + 1] + x[i]) * (x[N + i + 1] - x[N + i]) - (x[i] + x[i - 1]) * (x[N + i] - x[N + i - 1]))) * invh2;
        }
        f[0] = (x[1] - x[0] + 0.5 * (x[1] + x[0]) * (x[N + 1] - x[N])) * invh2;                                  // Left BC
        f[N - 1] = -(x[N - 1] - x[N - 2] + 0.5 * (x[N - 1] + x[N - 2]) * (x[2 * N - 1] - x[2 * N - 2])) * invh2; // Right BC

        // RHS for the potential Phi
        for (int_type i = 1; i < N - 1; i++)
        {
            f[i + N] = (x[i + 1 + N] - 2.0 * x[i + N] + x[i - 1 + N]) * invh2 - (1.0 - x[i]) * invlam2;
        }
        f[N] = x[N] + t;                 // Left BC
        f[2 * N - 1] = x[2 * N - 1] - t; // Right BC
    }
};

/*
 * Analytic Jacobian in sparse format.
 * The DAE solver will use automatic (algorithmic) differentiation if analytic Jacobian is not provided.
 * However, providing the analytic Jacobian can significantly speed up the computation for large systems.
 */
class MyJacobian
{
    MyParams p; // Parameters

public:
    explicit MyJacobian(MyParams &params) : p(params) {}

    /*
     * Defines the Jacobian matrix (matrix of the RHS derivatives) for the DAE system `M dx/dt = f`.
     * Takes vector `x` and time `t` and returns the Jacobian matrix `J`.
     * Matrix `J` is empty and should be filled with non-zero elements.
     */
    void operator()(sparse_matrix &J, const state_vector &x, const double t)
    {
        // Derived parameters
        const double h = p.L / (double)(p.N - 1);           // Cell size
        const double invh2 = 1.0 / (h * h);                 // Inverse cell size squared
        const double invlam2 = 1.0 / (p.lambda * p.lambda); // Inverse lambda squared

        // An alias for p.N
        const int_type N = p.N;

        J.reserve(12 * N); // Overestimating but it's better than underestimate

        // Fills Jacobian row by row
        for (int_type i = 0; i < 2 * N; ++i)
        {
            if (i == 0)
            {
                J(i, 0, (-1.0 + 0.5 * (x[N + 1] - x[N])) * invh2);
                J(i, 1, (1.0 + 0.5 * (x[N + 1] - x[N])) * invh2);
                J(i, N, -0.5 * (x[0] + x[1]) * invh2);
                J(i, N + 1, 0.5 * (x[0] + x[1]) * invh2);
            }
            else if (i < N - 1)
            {
                J(i, i - 1, (1.0 - 0.5 * (x[N + i] - x[N + i - 1])) * invh2);
                J(i, i, (-2.0 + 0.5 * (x[N + i + 1] - 2.0 * x[N + i] + x[N + i - 1])) * invh2);
                J(i, i + 1, (1.0 + 0.5 * (x[N + i + 1] - x[N + i])) * invh2);
                J(i, N + i - 1, 0.5 * (x[i] + x[i - 1]) * invh2);
                J(i, N + i, -0.5 * (x[i + 1] + 2.0 * x[i] + x[i - 1]) * invh2);
                J(i, N + i + 1, 0.5 * (x[i + 1] + x[i]) * invh2);
            }
            else if (i == N - 1)
            {
                J(i, i - 1, (1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2);
                J(i, i, (-1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2);
                J(i, N + i - 1, 0.5 * (x[N - 1] + x[N - 2]) * invh2);
                J(i, N + i, -0.5 * (x[N - 1] + x[N - 2]) * invh2);
            }
            else if (i == N)
            {
                J(i, N, 1.0);
            }
            else if (i < 2 * N - 1)
            {
                J(i, i - N, invlam2);
                J(i, i - 1, invh2);
                J(i, i, -2.0 * invh2);
                J(i, i + 1, invh2);
            }
            else if (i == 2 * N - 1)
            {
                J(i, 2 * N - 1, 1.0);
            }
            else
            {
                ERROR("Perovskite model: index i is out of boundaries.");
            }
        }
    }
};

/*
 * User-defined Solution Manager to post-process solution every time step.
 * In this example, it works as a passive observer and as an event function.
 */
class MySolutionManager
{
    // A reference to the solution holder object
    SolutionHolder &m_sol;

public:
    explicit MySolutionManager(SolutionHolder &sol) : m_sol(sol) {}

    /*
     * Solution Manager functor will be called every time step providing the time `t` and
     * the corresponding solution `x` for further post-processing.
     * If the functor returns an integer != 0 (`true`), the computation will immediately stop.
     */
    int operator()(const state_vector &x, const double t)
    {
        m_sol.x.emplace_back(x);
        m_sol.t.emplace_back(t);

        if (x[0] < 0)
        {
            return -1; // if x[0] is less than 0 (should never happen), then the solver will stop
        }
        else
        {
            return 0;
        }
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 */
int main()
{
    // Solution parameters
    MyParams params;
    params.N = 4000;

    // Initial condition
    state_vector x0(2 * params.N);
    for (int_type i = 0; i < params.N; ++i)
    {
        x0[i] = 1.0;
    }

    // Solution interval: t = [0, t_end]
    double t_end{10.0};

    // Solution holder
    SolutionHolder sol;

    // Solver options
    SolverOptions opt;
    opt.verbosity = verbosity::normal;        // Prints computation time and basic info
    opt.solution_variability_control = false; // Switches off solution variability control for better performance

    // Solve the DAE system
    int status = solve(MyMassMatrix(params), MyRHS(params), MyJacobian(params),
                       x0, t_end,
                       MySolutionManager(sol), opt);

    // Soluton vs time `t` is in the `sol` object.
    // Print P at the left boundary, and Phi at the right boundary of the domain.
    std::cout << "MySolutionManager: "
              << "t = " << sol.t.back()
              << ", P_left = " << sol.x.back()[0]
              << ", Phi_right = " << sol.x.back().back()
              << '\n';

    return status;
}
