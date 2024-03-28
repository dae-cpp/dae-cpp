/*
 * The main solver class definition.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_SOLVER_H
#define DAECPP_SOLVER_H

#include <iomanip>

#include "jacobian.hpp"
#include "mass-matrix.hpp"
#include "rhs.hpp"
#include "solver-options.hpp"
#include "timer.hpp"

namespace daecpp_namespace_name
{
namespace core
{

/*
 * Stores the solver state
 */
struct SolverState
{
    double t[MAX_ORDER];  // Current and previous integration times
    double dt[MAX_ORDER]; // Current and previous time steps

    std::array<state_type, MAX_ORDER> x; // Current and previous states

    int order{1}; // Current integration order (always starts from 1)

    SolverState(const std::size_t size)
    {
        for (int i = 0; i < MAX_ORDER; ++i)
        {
            t[i] = 0.0;
            dt[i] = 0.0;
            x[i].resize(size);
        }
    }
};

} // namespace core

/*
 * The main solver class
 */
class System
{
    const RHS &_rhs;           // RHS
    const Jacobian &_jac;      // Jacobian matrix
    const MassMatrix &_mass;   // Mass matrix
    const SolverOptions &_opt; // Solver options

    std::vector<double> _t_out; // Vector of output times

public:
    System(const MassMatrix &mass, const RHS &rhs)
        : _mass(mass), _rhs(rhs), _jac(JacobianNumerical(rhs)), _opt(SolverOptions()) {}

    System(const MassMatrix &mass, const RHS &rhs, const Jacobian &jac)
        : _mass(mass), _rhs(rhs), _jac(jac), _opt(SolverOptions()) {}

    System(const MassMatrix &mass, const RHS &rhs, const SolverOptions &opt)
        : _mass(mass), _rhs(rhs), _jac(JacobianNumerical(rhs)), _opt(opt) {}

    System(const MassMatrix &mass, const RHS &rhs, const Jacobian &jac, const SolverOptions &opt)
        : _mass(mass), _rhs(rhs), _jac(jac), _opt(opt) {}

    /*
     * Integrates the system of DAEs on the interval t = [t0; t1] and returns
     * result in the array x. Parameter t0 can be overriden in the solver
     * options (t0 = 0 by default).
     * The data stored in x (initial conditions) will be overwritten.
     * Returns 0 in case of success or error code if integration is failed.
     */
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
    int solve(state_type &x, double &t_end, const std::vector<double> t_output = {})
    {
        core::Time time;
        // Timer
        {
            // Measures total time
            // Timer timer(&core::Timers::get().total_time);
            Timer timer(&time.total);

            // Initial output
            PRINT(_opt.verbosity >= 1, "Starting dae-cpp solver...");

            // Initialize output times
            _t_out = std::move(t_output);

            // Sort vector of output times and erase duplicates
            _t_out.push_back(t_end);
            std::sort(_t_out.begin(), _t_out.end());
            _t_out.erase(std::unique(_t_out.begin(), _t_out.end()), _t_out.end());

            // Throw an error if target time t < 0
            if (_t_out.back() < 0.0)
            {
                ERROR("Target time t cannot be negative. The solver integrates from 0 to t.");
            }

            // Check user-defined solver options
            _opt.check();

            // System size
            std::size_t size = x.size();
            ASSERT(size > 0, "Initial condition vector is empty.");

            // Solver state
            core::SolverState state(size);

            // Set initial time step
            state.dt[0] = _opt.dt_init;

            // Copy initial state
            state.x[0] = x;

            // Solution vector used for Newton iterations
            state_type xk = x;

            // The RHS vector
            state_type f(size);

            // Mass matrix holder
            sparse_matrix M;

            // Jacobian matrix holder
            sparse_matrix J;

            // Time derivative approximation
            core::eivec dxdt(size);

            // Counts number of time steps
            u_int64_t steps{0};

            // Counts total linear algebra solver calls
            u_int64_t calls{0};

            // Counts how many times the Newton iterator failed to converge within
            // max_Newton_iter iterations in a row.
            // unsigned int n_iter_failed = 0;

            // Output after initialization
            PRINT(_opt.verbosity >= 2, "Float size:      " << 8 * sizeof(float_type) << " bit");
            PRINT(_opt.verbosity >= 2, "Integer size:    " << 8 * sizeof(int_type) << " bit");
            PRINT(_opt.verbosity >= 1, "DAE system size: " << size << " equations");
            PRINT(_opt.verbosity >= 1, "Calculating...");

            /*
             * Output time loop
             */
            for (const auto &t1 : _t_out)
            {
                PRINT(_opt.verbosity >= 2, "-- Integration time t = " << t1 << ":");

                if (t1 < 0.0)
                {
                    WARNING("Negative integration time t = " << t1 << ". Skipped.");
                    continue;
                }

                // TODO: Check initial time step is less than t
                // TODO: For t=0 do at least one iteration

                /*
                 * Time loop
                 */
                while (state.t[0] < t1)
                {
                    state.t[0] += state.dt[0]; // Time step lapse

                    steps++; // Number of time steps

                    if (_opt.verbosity >= 1)
                    {
                        std::cout << std::left
                                  << "Step " << std::setw(7) << steps
                                  << " :: t = " << std::setw(12) << state.t[0]
                                  << " :: ";
                    }

                    if (_opt.verbosity >= 2)
                    {
                        std::cout << "BDF-" << state.order
                                  << ": dt=" << state.dt[0]
                                  << " :: ";
                    }

                    // // Can be set to true by the solver if it fails to converge
                    // bool fact_every_iter = (n_iter_failed >= m_opt.newton_failed_attempts)
                    //                            ? true
                    //                            : m_opt.fact_every_iter;

                    u_int32_t iter; // Newton iteration loop index - we will need this value later

                    /*
                     * Newton iteration loop
                     */
                    for (iter = 0; iter < _opt.max_Newton_iter; ++iter)
                    {
                        // Reordering, Symbolic and Numerical Factorization (slow) + recalculate Jac
                        // if(fact_every_iter || iter == 0 || !(iter % m_opt.fact_iter))

                        // Time integrator - find full Jb and b (from x,t,dt - state)
                        // Returns time derivative approximation and its corresponding derivative w.r.t. xk
                        double alpha = time_derivative_approx(dxdt, xk, state, size);

                        // Get RHS
                        _rhs(f, xk, state.t[0]);
                        ASSERT(f.size() == size, "The RHS vector size (" << f.size() << ") does not match the initial condition vector size (" << size << ").");
                        Eigen::Map<core::eivec> f_(f.data(), f.size()); // Does not copy data... check

                        // If above crashes due to alignment:
                        // core::eivec f_ = Eigen::Map<core::eivec, Eigen::Unaligned>(f.data(), f.size()); // Makes a copy

                        // Get Mass Matrix
                        M.clear();
                        _mass(M, state.t[0]);
                        M.check();
                        auto M_ = M.convert(size);

                        // Get Jac
                        J.clear();
                        _jac(J, xk, state.t[0]);
                        J.check();
                        auto Jb = J.convert(size);

                        // b = M(t) * dxdt - f // Note it's with '-'
                        core::eivec b = M_ * dxdt;
                        b -= f_;

                        // Jb = J - d/dx (M * [dx/dt])
                        Jb -= M_ * alpha; /// state.dt[0];

                        // Solve linear system Jb dx = b
                        Eigen::SparseLU<Eigen::SparseMatrix<float_type>> lu(Jb); // LU method
                        auto dx = lu.solve(b);
                        calls++;
                        std::putchar('#');

                        // x_k+1 = x_k - delta_x
                        for (int_type i = 0; i < size; ++i)
                            xk[i] += dx[i];

                        // Check error and convergence -- loop for i
                        // Checks NaN, atol, rtol, and x[i] -= dx[i]

                        // break if convereged

                    } // Newton iteration loop

                    // Updates state history
                    for (std::size_t i = 0; i < size; ++i)
                    {
                        // Manually unrolled loop
                        state.x[3][i] = state.x[2][i];
                        state.x[2][i] = state.x[1][i];
                        state.x[1][i] = state.x[0][i];
                        state.x[0][i] = xk[i];
                    }

                    // Updates time step history
                    for (int k = core::MAX_ORDER - 1; k > 0; --k)
                    {
                        state.dt[k] = state.dt[k - 1];
                    }

                    state.dt[0] = 0.09; // New time step

                    // Updates time integration order
                    if (state.order < core::MAX_ORDER)
                    {
                        state.order++;
                    }

                    // Newton iteration finished
                    std::putchar('\n');

                    // Adjust the last time step if needed
                    if (t1 - state.t[0] < state.dt[0])
                    {
                        state.dt[0] = t1 - state.t[0];
                    }

                    // We may already reached the target time
                    if (state.dt[0] < core::TIMESTEP_ROUNDING_ERROR)
                    {
                        break;
                    }

                } // Time loop

            } // for (const auto &t1 : _t_out)

            // Return solution and the final time
            x = state.x[0];
            t_end = state.t[0];

        } // Global timer

        // Final output
        PRINT(_opt.verbosity >= 1, "Total time: " << time.total << " ms");

        return 0;
    }

private:
    inline double time_derivative_approx(core::eivec &dxdt, const state_type &xk, const core::SolverState &state, const std::size_t size) const
    {
        double alpha{0.0}; // Derivative w.r.t. xk

        const double h0 = state.dt[0];
        const double h1 = state.dt[1];
        const double h2 = state.dt[2];
        const double h3 = state.dt[3];

        const double h01 = h0 + h1;
        const double h12 = h1 + h2;
        const double h23 = h2 + h3;
        const double h012 = h0 + h1 + h2;
        const double h123 = h1 + h2 + h3;
        const double h0123 = h0 + h1 + h2 + h3;

        switch (state.order)
        {
        case 1:
            alpha = 1.0 / h0;
            for (std::size_t i = 0; i < size; ++i)
            {
                dxdt[i] = (xk[i] - state.x[0][i]) * alpha;
            }
            break;

        case 2:
            alpha = (2.0 * h0 + h1) / (h0 * h01);
            for (std::size_t i = 0; i < size; ++i)
            {
                dxdt[i] = alpha * xk[i] -
                          h01 / (h0 * h1) * state.x[0][i] +
                          h0 / (h1 * h01) * state.x[1][i];
            }
            break;

        case 3:
            alpha = (3.0 * h0 * h0 + h1 * h12 + 2.0 * h0 * (2.0 * h1 + h2)) / (h0 * h01 * h012);
            for (std::size_t i = 0; i < size; ++i)
            {
                dxdt[i] = alpha * xk[i] -
                          h01 * h012 / (h0 * h1 * h12) * state.x[0][i] +
                          h0 * h012 / (h1 * h2 * h01) * state.x[1][i] -
                          h0 * h01 / (h2 * h12 * h012) * state.x[2][i];
            }
            break;

        case 4:
            alpha = (4.0 * h0 * h0 * h0 +
                     h1 * h12 * h123 +
                     3.0 * h0 * h0 * (3.0 * h1 + 2.0 * h2 + h3) +
                     2.0 * h0 * (3.0 * h1 * h1 + h2 * h23 + 2.0 * h1 * (2.0 * h2 + h3))) /
                    (h0 * h01 * h012 * h0123);
            for (std::size_t i = 0; i < size; ++i)
            {
                dxdt[i] = alpha * xk[i] -
                          h01 * h012 * h0123 / (h0 * h1 * h12 * h123) * state.x[0][i] +
                          h0 * h012 * h0123 / (h1 * h01 * h2 * h23) * state.x[1][i] -
                          h0 * h01 * h0123 / (h2 * h12 * h012 * h3) * state.x[2][i] +
                          h0 * h01 * h012 / (h3 * h23 * h123 * h0123) * state.x[3][i];
            }
            break;

        default:
            ERROR("Unsupported time integration order.");
        }

        return alpha;
    }

    // /*
    //  * Virtual Observer. Called by the solver every time step.
    //  * Receives current solution vector and the current time.
    //  * Can be overriden by a user to get intermediate results (for example,
    //  * for plotting).
    //  */
    // virtual void observer(state_type &x, const double t)
    // {
    //     return;  // It does nothing by default
    // }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_H
