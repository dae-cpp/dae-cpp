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

    std::vector<double> _t_out; // Output times

    std::size_t _steps{0}; // Counts number of time steps
    std::size_t _calls{0}; // Counts total linear algebra solver calls

    core::Time _time;

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
        // Timer
        {
            // Measures total time
            // Timer timer(&core::Timers::get().total_time);
            Timer timer(&_time.total);

            // Initial output
            PRINT(_opt.verbosity >= 1, "Starting dae-cpp solver...");

            // Initialize output times
            _t_out = std::move(t_output);

            // Solver state
            core::SolverState state;

            // System size
            std::size_t size = x.size();
            ASSERT(size > 0, "Initial condition vector is empty.");

            // Reserve memory for the solution history
            for (auto &x_ : state.x)
            {
                x_.reserve(size);
            }

            // Initial time
            state.t[0] = 0.0;
            state.t[1] = 0.0;
            state.t[2] = 0.0;

            // Initial time step
            state.dt[0] = _opt.dt_init;
            state.dt[1] = 0.0;
            state.dt[2] = 0.0;

            // Copy initial state
            state.x[0] = x;
            state.x[1].resize(size);
            state.x[2].resize(size);

            // TODO: Check user-defined solver options
            // _opt.check();

            // Full Jacobian matrix holder
            // sparse_matrix_holder J;

            // Full RHS vector
            // state_type b(m_size);

            // Solution vector used for Newton iterations
            // state_type xk(m_size);

            // Counts how many times the Newton iterator failed to converge within
            // max_Newton_iter iterations in a row.
            // int n_iter_failed = 0;

            // Sort vector of output times and erase duplicates
            _t_out.push_back(t_end);
            std::sort(_t_out.begin(), _t_out.end());
            _t_out.erase(std::unique(_t_out.begin(), _t_out.end()), _t_out.end());

            // Throw an error if target time t < 0
            if (_t_out.back() < 0.0)
            {
                ERROR("Target time t cannot be negative. The solver integrates from 0 to t.");
            }

            // Output after initialization
            PRINT(_opt.verbosity >= 2, "Float size:      " << 8 * sizeof(float_type) << " bit");
            PRINT(_opt.verbosity >= 2, "Integer size:    " << 8 * sizeof(int_type) << " bit");
            PRINT(_opt.verbosity >= 1, "DAE system size: " << size << " equations");
            PRINT(_opt.verbosity >= 1, "Calculating...");

            state_type xk = x; // TODO:
            state_type f(size);
            core::eivec dxdt(size); // TODO: Can be single

            // state_type b

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
                while (state.t[0] < (t1 + state.dt[0] * 0.5))
                {
                    std::cout << state.t[0] << "\t" << state.x[0][0] << "\t" << state.x[0][1] << "\t" << std::exp(-state.t[0]) << "\t" << -std::exp(-state.t[0]) << '\n';

                    state.t[0] += state.dt[0]; // Time step lapse

                    _steps++; // Number of time steps

                    // if(m_opt.verbosity > 1)
                    // {
                    //     std::cout << std::left;
                    //     std::cout << "\nStep " << std::setw(7) << m_steps
                    //               << " :: t = " << std::setw(12) << m_iterator_state.t
                    //               << " :: ";
                    //     // std::cout.flush();  // This degrades performance in some cases
                    // }

                    // if(m_opt.verbosity > 2)
                    // {
                    //     std::cout << "BDF-" << m_iterator_state.current_scheme << ": ";
                    //     std::cout << "dt=" << m_iterator_state.dt[0]
                    //               << ", dt_prev=" << m_iterator_state.dt[1] << ": ";
                    // }

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
                        // std::cout << "  -- "<< iter << "\t" << state.x[0][0] << "\t" << state.x[0][1] << '\n';

                        // Reordering, Symbolic and Numerical Factorization (slow) + recalculate Jac
                        // if(fact_every_iter || iter == 0 || !(iter % m_opt.fact_iter))

                        // Time integrator - find full Jb and b (from x,t,dt - state)
                        // Returns time derivative approximation and its corresponding derivative w.r.t. xk
                        double alpha = time_derivative_approx(dxdt, xk, state, size);

                        // Get RHS
                        _rhs(f, xk, state.t[0]);
                        Eigen::Map<core::eivec> f_(f.data(), f.size()); // Does not copy data... check

                        // If above crashes due to alignment:
                        // core::eivec f_ = Eigen::Map<core::eivec, Eigen::Unaligned>(f.data(), f.size()); // Makes a copy

                        // Get Mass Matrix
                        sparse_matrix M; // Create in advance then reset
                        _mass(M, state.t[0]);
                        M.check();
                        auto M_ = M.convert(size);

                        // Get Jac
                        sparse_matrix J; // Create in advance then reset
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
                        auto xk_ = lu.solve(b);

                        // x_k+1 = x_k - delta_x
                        for (int_type i = 0; i < size; ++i)
                            xk[i] += xk_[i];

                        // std::cout << "  -- "<< iter << "\t" << state.x[0][0] << "\t" << state.x[0][1] << '\n';

                        // break;

                        // calls++;

                        // Check error and convergence -- loop for i
                        // Checks NaN, atol, rtol, and x[i] -= dx[i]

                        // Print # - one char

                        // break if convereged
                    }

                    for (int_type i = 0; i < size; ++i)
                    {
                        state.x[2][i] = state.x[1][i];
                        state.x[1][i] = state.x[0][i];
                        state.x[0][i] = xk[i];
                    }

                    state.order++;

                    state.dt[2] = state.dt[1];
                    state.dt[1] = state.dt[0];
                    state.dt[0] = 0.09;

                    // break;

                } //
            }     // for (const auto &t1 : _t_out)
        }

        // std::cout << "Total time: " << core::Timers::get().total_time << '\n';
        std::cout << "Total time: " << _time.total << '\n';
        return 0;
    }

private:
    inline double time_derivative_approx(core::eivec &dxdt, const state_type &xk, const core::SolverState &state, const std::size_t size)
    {
        double alpha{0.0}; // Derivative w.r.t. xk
        const double *dt = state.dt;
        if (state.order == 1)
        {
            alpha = 1.0 / dt[0];
            for (int_type i = 0; i < size; ++i)
            {
                dxdt[i] = (xk[i] - state.x[0][i]) * alpha;
            }
        }
        else if (state.order == 2)
        {
            const double dt01 = dt[0] + dt[1];
            alpha = (2.0 * dt[0] + dt[1]) / (dt[0] * dt01);
            for (int_type i = 0; i < size; ++i)
            {
                dxdt[i] = alpha * xk[i] -
                          dt01 / (dt[0] * dt[1]) * state.x[0][i] +
                          dt[0] / (dt[1] * dt01) * state.x[1][i];
            }
        }
        else
        {
            const double dt01 = dt[0] + dt[1];
            const double dt12 = dt[1] + dt[2];
            const double dt02 = dt[0] + dt[2];
            const double dt012 = dt[0] + dt[1] + dt[2];
            alpha = (3.0 * dt[0] * dt[0] + dt[1] * dt12 + 2.0 * dt[0] * (dt[1] + dt12)) / (dt[0] * dt01 * dt012);
            for (int_type i = 0; i < size; ++i)
            {
                dxdt[i] = alpha * xk[i] -
                          dt01 * dt012 / (dt[0] * dt[1] * dt12) * state.x[0][i] +
                          dt[0] * dt012 / (dt[1] * dt[2] * dt01) * state.x[1][i] -
                          dt[0] * dt01 / (dt[2] * dt12 * dt012) * state.x[2][i];
            }
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
