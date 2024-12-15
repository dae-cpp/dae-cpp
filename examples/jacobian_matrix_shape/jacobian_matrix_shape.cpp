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
 * Copyright (c) 2024 Ivan Korotkin
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
 * The vector-function (RHS) of the problem.
 * Based on `VectorFunctionJacobianShape` class.
 */
class JacShapeRHS : public VectorFunctionJacobianShape
{
    const MyParams p; // Parameters

    // Derived parameters
    const double h = p.L / (double)(p.N - 1);           // Cell size
    const double invh2 = 1.0 / (h * h);                 // Inverse cell size squared
    const double invlam2 = 1.0 / (p.lambda * p.lambda); // Inverse lambda squared

    /*
     * RHS for the ion concentration P = dFlux/dx:
     * dP/dt = d/dx(dP/dx + P * dPhi/dx)
     */
    dual_type eq1(const state_type &x, const double t, const int_type i_global) const
    {
        const dual_type *P = x.data();
        const dual_type *Phi = x.data() + p.N;

        const int_type i = i_global;

        if (i == 0)
            return (P[i + 1] - P[i] + 0.5 * (P[i + 1] + P[i]) * (Phi[i + 1] - Phi[i])) * invh2; // Left BC
        else if (i == p.N - 1)
            return -(P[i] - P[i - 1] + 0.5 * (P[i] + P[i - 1]) * (Phi[i] - Phi[i - 1])) * invh2; // Right BC
        else
            return (P[i + 1] - 2.0 * P[i] + P[i - 1] + 0.5 * ((P[i + 1] + P[i]) * (Phi[i + 1] - Phi[i]) - (P[i] + P[i - 1]) * (Phi[i] - Phi[i - 1]))) * invh2;
    }

    /*
     * RHS for the potential Phi:
     * d^2(Phi)/dx^2 - (1 - P)/lambda^2 = 0
     */
    dual_type eq2(const state_type &x, const double t, const int_type i_global) const
    {
        const dual_type *P = x.data();
        const dual_type *Phi = x.data() + p.N;

        const int_type i = i_global - p.N;

        if (i == 0)
            return Phi[i] + t; // Left BC
        else if (i == p.N - 1)
            return Phi[i] - t; // Right BC
        else
            return (Phi[i + 1] - 2.0 * Phi[i] + Phi[i - 1]) * invh2 - (1.0 - P[i]) * invlam2;
    }

public:
    explicit JacShapeRHS(MyParams &params) : p(params) {}

    /*
     * All equations combined
     */
    dual_type equations(const state_type &x, const double t, const int_type i) const
    {
        if (i < p.N)
            return eq1(x, t, i);
        else if (i < 2 * p.N)
            return eq2(x, t, i);
        else
        {
            ERROR("Perovskite model: index i is out of boundaries.");
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

    // Fill initial condition
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

    {
        std::cout << "-- Automatic Jacobian derived from the user-defined shape:\n\n";

        JacobianMatrixShape jac = JacobianMatrixShape(JacShapeRHS(params));

        // An alias for p.N - the number of discretization points for each equation
        const int_type N = params.N;

        // Defines all non-zero elements of the Jacobian matrix row by row
        for (int_type i = 0; i < 2 * N; ++i)
        {
            if (i == 0)
            {
                jac.add_element(i, {i + 1, i, N + i + 1, N + i}); // The order does not matter
            }
            else if (i < N - 1)
            {
                jac.add_element(i, {i - 1, i, i + 1, N + i - 1, N + i, N + i + 1});
            }
            else if (i == N - 1)
            {
                jac.add_element(i, {i - 1, i, N + i - 1, N + i});
            }
            else if (i == N)
            {
                jac.add_element(i, N); // Can pass one index instead of a vector
            }
            else if (i < 2 * N - 1)
            {
                jac.add_element(i, {i - N, i - 1, i, i + 1});
            }
            else if (i == 2 * N - 1)
            {
                jac.add_element(i, 2 * N - 1);
            }
            else
            {
                ERROR("Perovskite model: index i is out of boundaries.");
            }
        }

        // Solve the DAE system using automatic Jacobian computed from the user-defined shape
        solve(MyMassMatrix(params), JacShapeRHS(params), jac,
              x0, t_end,
              MySolutionManager(sol), opt);
    }

    return 0;
}
