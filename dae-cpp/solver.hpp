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

#include <Eigen/SparseLU>

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
    double dt[DAECPP_MAX_ORDER]; // Current and previous time steps

    double t{0.0};      // Current integration time
    double t_prev{0.0}; // Previous integration time

    std::array<state_type, DAECPP_MAX_ORDER> x; // Current and previous states

    int order{1}; // Current integration order (always starts from 1)

    /*
     * Resets solver state and allocates memory for solution history
     */
    SolverState(const std::size_t size)
    {
        for (int i = 0; i < DAECPP_MAX_ORDER; ++i)
        {
            dt[i] = 0.0;

            try
            {
                x[i].resize(size);
            }
            catch (const std::exception &e)
            {
                ERROR("Failed to allocate memory for the solver state.\n"
                      << e.what());
            }
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
     * Integrates the system of DAEs in the interval `t = [0; t_end]`.
     * Returns solution `x` at time `t_end` (the initial condition stored in `x` will be overwritten).
     * If the solver stops earlier, `t_end` will be overwritten with the actual time.
     *
     * Parameters:
     *     x - initial condition (daecpp::state_type)
     *     t_end - integration interval `t = [0; t_end]` (double)
     *     t_output - (optional) a vector of output times (std::vector<double>)
     *
     * Returns:
     *     0 if integration is successful or error code if integration is failed (int)
     */
    int solve(state_type &x, double &t_end, const std::vector<double> t_output = {})
    {
        // Specific timers
        core::Time time;

        // Counts total linear algebra solver calls
        u_int64_t n_calls{0};

        // Global timer
        {
            // Measures total time
            Timer timer_global(&time.total);

            // Initial output
            PRINT(_opt.verbosity >= 1, "Starting dae-cpp solver...");

            // Vector of output times
            std::vector<double> t_out = std::move(t_output);

            // Sort vector of output times and erase duplicates
            t_out.push_back(t_end);
            std::sort(t_out.begin(), t_out.end());
            t_out.erase(std::unique(t_out.begin(), t_out.end()), t_out.end());

            // Throw an error if target time t < 0
            ASSERT(t_out.back() >= 0.0, "Target time t_end cannot be negative. The solver integrates from 0 to t_end.");

            // Check user-defined solver options
            _opt.check();

            // System size
            auto size = x.size();
            ASSERT(size > 0, "Initial condition vector x is empty.");

            // Solver state
            core::SolverState state(size);

            // Alias for the current time step
            double &dt = state.dt[0];

            // Set initial time step
            dt = _opt.dt_init;

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

            // Linear solver
            Eigen::SparseLU<core::eimat> linsolver;

            // Eigen::SparseMatrix<float_type> matrices
            core::eimat M_; // Mass matrix (converted)
            core::eimat Jb; // Linear system matrix

            // Eigen::VectorX vectors
            core::eivec f_; // The RHS vector (converted)
            core::eivec b;  // The RHS of the linear system
            core::eivec dx; // Linear system solution

            // Counts number of time steps
            u_int64_t n_steps{0};

            // Counts how many times the Newton iterator failed to converge in a row
            u_int32_t n_iter_failed{0};

            // Output after initialization
            PRINT(_opt.verbosity >= 2, "Float size:      " << 8 * sizeof(float_type) << " bit");
            PRINT(_opt.verbosity >= 2, "Integer size:    " << 8 * sizeof(int_type) << " bit");
            PRINT(_opt.verbosity >= 1, "DAE system size: " << size << " equations");
            PRINT(_opt.verbosity >= 2, "Calculating...");

            /*
             * Output time loop
             */
            for (const auto &t1 : t_out)
            {
                PRINT(_opt.verbosity >= 2, "\n-- Integration time t = " << t1 << ":");

                if (t1 < 0.0)
                {
                    WARNING("Negative integration time t = " << t1 << ". Skipped.");
                    continue;
                }

                // Adjust the initial time step if needed
                if (dt > t1 - state.t)
                {
                    dt = t1 - state.t;
                    ASSERT(dt > 0.0, "Negative or zero time step: dt = "
                                         << dt << ".\n"
                                         << "This assertion triggered after adjusting the time step to match the output time t_output = "
                                         << t1 << ".");
                }

                /*
                 * Time loop
                 */
                while (state.t < t1)
                {
                    state.t_prev = state.t; // Save the previous time
                    state.t += dt;          // Time step lapse

                    n_steps++; // Number of time steps

                    if (_opt.verbosity >= 1)
                    {
                        std::cout << std::left
                                  << "Step " << std::setw(8) << n_steps
                                  << " :: t = " << std::setw(11) << state.t
                                  << " :: ";
                    }

                    if (_opt.verbosity >= 2)
                    {
                        std::cout << "BDF-" << state.order
                                  << ": dt=" << std::setw(11) << dt
                                  << " :: ";
                    }

                    bool is_diverged{false}; // True if Newton iterations diverged

                    u_int32_t iter{}; // Newton iteration loop index - we will need this value later

                    /*
                     * Newton iteration loop
                     */
                    for (iter = 0; iter < _opt.max_Newton_iter; ++iter)
                    {
                        // Returns time derivative approximation dxdt and its corresponding derivative w.r.t. xk
                        double alpha{}; // Derivative w.r.t. xk
                        {
                            Timer timer(&time.time_derivative);
                            alpha = time_derivative_approx(dxdt, xk, state, size);
                        }

                        // Get and convert the RHS
                        {
                            Timer timer(&time.rhs);
                            _rhs(f, xk, state.t);
                            ASSERT(f.size() == size, "The RHS vector size (" << f.size() << ") does not match the initial condition vector size (" << size << ").");
                            f_ = Eigen::Map<core::eivec, Eigen::Unaligned>(f.data(), f.size());
                        }

                        // Get and convert the Mass matrix
                        if (!_opt.is_mass_matrix_static || (_opt.is_mass_matrix_static && !iter && !n_calls))
                        {
                            Timer timer(&time.mass);
                            M.clear();
                            _mass(M, state.t);
                            M.check();
                            M_ = M.convert(size);
                        }

                        // Get and convert the Jacobian matrix
                        if (iter % 2 == 0) // TODO: Add solver option
                        {
                            Timer timer(&time.jacobian);
                            J.clear();
                            _jac(J, xk, state.t);
                            J.check();
                            Jb = J.convert(size);
                        }

                        // Matrix-vector operations
                        {
                            Timer timer(&time.linear_algebra);

                            // b = M(t) * [dx/dt] - f
                            b = M_ * dxdt;
                            b -= f_;

                            // Jb = J - d/dxk (M(t) * [dx/dt])
                            if (iter % 2 == 0) // TODO: Add solver option
                            {
                                Jb -= M_ * alpha;
                            }
                        }

                        // Factorization
                        if (iter % 2 == 0) // TODO: Add solver option
                        {
                            Timer timer(&time.factorization);
                            linsolver.compute(Jb);
                            if (linsolver.info() != Eigen::Success)
                            {
                                ERROR("Decomposition failed."); // TODO: Try to save
                            }
                        }

                        // Solve linear system Jb dx = b
                        {
                            Timer timer(&time.linear_solver);
                            dx = linsolver.solve(b);
                            if (linsolver.info() != Eigen::Success)
                            {
                                ERROR("Solving failed."); // TODO: Try to save
                            }
                            n_calls++;
                            print_char(_opt.verbosity >= 1, '#');
                        }

                        bool is_converged = true; // Assume the iterations converged

                        // Check convergence/divergence and update xk
                        {
                            Timer timer(&time.error_check);

                            for (std::size_t i = 0; i < size; ++i)
                            {
                                // Absolute error
                                double err_abs = std::abs(dx[i]);

                                // Solution diverged. Roll back to the previous state and redo with reduced time step.
                                if (err_abs > _opt.max_value || std::isnan(dx[i]))
                                {
                                    PRINT(_opt.verbosity >= 1, " <- diverged");

                                    // Trying to roll back and reduce the time step
                                    is_diverged = true;
                                    n_iter_failed++;
                                    state.t = state.t_prev;
                                    n_steps--;
                                    dt /= _opt.dt_decrease_factor;
                                    if (dt < _opt.dt_min)
                                    {
                                        PRINT(_opt.verbosity >= 1, "The time step was reduced to `t_min` but the scheme failed to converge.");
                                        goto result; // Abort all loops and go straight to the results.
                                                     // Using goto here is much more clear than using a sequence of `break` statements.
                                    }

                                    break;
                                }

                                // Relative error check
                                if (state.x[0][i] != 0.0)
                                {
                                    double err_rel = err_abs / std::abs(state.x[0][i]);
                                    if (err_rel > _opt.rtol)
                                    {
                                        is_converged = false;
                                    }
                                }

                                // Absolute error check
                                if (err_abs > _opt.atol)
                                {
                                    is_converged = false;
                                }

                                // x_{k+1} = x_{k} - delta_x
                                xk[i] += dx[i];
                            }
                        }

                        // break if convereged (or diverged)
                        if (is_converged || is_diverged)
                        {
                            break;
                        }

                    } // Newton iteration loop

                    if (is_diverged)
                    {
                        continue;
                    }

                    {
                        Timer timer(&time.history);

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
                        for (int k = DAECPP_MAX_ORDER - 1; k > 0; --k)
                        {
                            state.dt[k] = state.dt[k - 1];
                        }
                    }

                    // Make decision about new time step
                    if (iter >= _opt.dt_decrease_threshold)
                    {
                        dt /= _opt.dt_decrease_factor;
                        print_char(_opt.verbosity >= 1, '<');
                        if (dt < _opt.dt_min)
                        {
                            PRINT(_opt.verbosity >= 1, " <- reached dt_min");
                            PRINT(_opt.verbosity >= 1, "The time step was reduced to `t_min` but the scheme failed to converge.");
                            goto result; // Abort all loops and go straight to the results.
                                         // Using goto here is much more clear than using a sequence of `break` statements.
                        }
                    }
                    if (iter <= _opt.dt_increase_threshold - 1)
                    {
                        dt *= _opt.dt_increase_factor;
                        if (dt > _opt.dt_max)
                        {
                            dt = _opt.dt_max;
                        }
                    }

                    // Adjust the last time step if needed
                    if (dt > t1 - state.t)
                    {
                        dt = t1 - state.t;
                    }

                    // Updates time integration order
                    if (state.order < _opt.BDF_order)
                    {
                        state.order++;
                    }

                    // Newton iteration finished
                    print_char(_opt.verbosity >= 1, '\n');

                    // We may already reached the target time
                    if (dt < DAECPP_TIMESTEP_ROUNDING_ERROR)
                    {
                        dt = state.dt[1];
                        break;
                    }

                } // Time loop

            } // for (const auto &t1 : t_out) - Loop over all output times

        result:

            // Return solution and the final time
            x = state.x[0];
            t_end = state.t;

        } // Global timer

        // Final output TODO: Move to a function
        // TODO: if time.total quite big, print time in seconds
        PRINT(_opt.verbosity >= 1, "\nTime spent by:");
        PRINT(_opt.verbosity >= 1, "  Time derivative approximation: " << time.time_derivative << " ms");
        PRINT(_opt.verbosity >= 1, "  RHS computation:               " << time.rhs << " ms");
        PRINT(_opt.verbosity >= 1, "  Mass matrix computation:       " << time.mass << " ms");
        PRINT(_opt.verbosity >= 1, "  Jacobian matrix computation:   " << time.jacobian << " ms");
        PRINT(_opt.verbosity >= 1, "  Linear algebra operations:     " << time.linear_algebra << " ms");
        PRINT(_opt.verbosity >= 1, "  Matrix factorization:          " << time.factorization << " ms");
        PRINT(_opt.verbosity >= 1, "  Linear solver:                 " << time.linear_solver << " ms  <--  " << n_calls << " linear solver calls");
        PRINT(_opt.verbosity >= 1, "  Error control:                 " << time.error_check << " ms");
        PRINT(_opt.verbosity >= 1, "  Solution history update:       " << time.history << " ms");
        PRINT(_opt.verbosity >= 1, "  Initialization and other:      " << time.other() << " ms");
        PRINT(_opt.verbosity >= 1, "Total time:                      " << time.total << " ms");

        return 0;
    }

private:
    /*
     * Prints a single character if condition is true
     */
    inline void print_char(bool condition, char ch) const noexcept
    {
        if (condition)
        {
            std::putchar(ch);
        }
    }

    /*
     * Returns time derivative approximation and its corresponding derivative w.r.t. xk
     */
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
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_H
