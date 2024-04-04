/*
 * The main solver implementation.
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

#include "jacobian-matrix.hpp"
#include "mass-matrix.hpp"
#include "solver-options.hpp"
#include "timer.hpp"
#include "vector-function.hpp"

namespace daecpp_namespace_name
{
namespace core
{
namespace internal
{

/*
 * Stores current and previous states of the solver
 */
struct SolverState
{
    double dt[DAECPP_MAX_ORDER]; // Current and previous time steps

    double t{0.0};      // Current integration time
    double t_prev{0.0}; // Previous integration time

    std::array<rvec, DAECPP_MAX_ORDER> x; // Current and previous states

    unsigned int order{1}; // Current integration order (always starts from 1)

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

/*
 * Counts interations, function calls, etc.
 */
struct Counters
{
    u_int64_t n_fact_calls{0}; // Number of Jacobian matrix updates and factorizations
    u_int64_t n_lin_calls{0};  // Number of linear solver calls
};

/*
 * Prints a single character if condition is true
 */
inline void print_char(bool condition, char ch) noexcept
{
    if (condition)
    {
        std::putchar(ch);
    }
}

/*
 * Converts and formats simulation time, estimates percentage
 */
inline std::string print_time(double val, double t_total)
{
    double _t_coef{1.0};       // Default conversion coefficient for output
    std::string _t_unit{"ms"}; // Default time units for output

    if (t_total > 1e4)
    {
        _t_unit = "s";  // seconds
        _t_coef = 1e-3; // ms -> s
    }

    std::streamsize ss = std::cout.precision();
    std::ostringstream oss;
    oss << std::setw(11) << val * _t_coef << ' ' << _t_unit;
    oss << std::setprecision(3);
    oss << "  (" << val / t_total * 100.0 << "%)";
    oss << std::setprecision(ss);
    std::string var = oss.str();
    return var;
}

/*
 * Final output
 */
inline void finalize(const Time &t, int v, Counters c)
{
    PRINT(v >= 1, "\nComputation time:" << std::right);
    PRINT(v >= 1, "------------------------------------------------------------");
    PRINT(v >= 1, "  Initialization:           " << print_time(t.init, t.total));
    PRINT(v >= 1, "  Time derivative:          " << print_time(t.time_derivative, t.total));
    PRINT(v >= 1, "  RHS:                      " << print_time(t.rhs, t.total));
    PRINT(v >= 1, "  Mass matrix:              " << print_time(t.mass, t.total));
    PRINT(v >= 1, "  Jacobian matrix:          " << print_time(t.jacobian, t.total));
    PRINT(v >= 1, "  Linear algebra:           " << print_time(t.linear_algebra, t.total));
    PRINT(v >= 1, "  Matrix factorization:     " << print_time(t.factorization, t.total) << "  <--  " << c.n_fact_calls << " calls");
    PRINT(v >= 1, "  Linear solver:            " << print_time(t.linear_solver, t.total) << "  <--  " << c.n_lin_calls << " calls");
    PRINT(v >= 1, "  Error control:            " << print_time(t.error_check, t.total));
    PRINT(v >= 1, "  Solution history update:  " << print_time(t.history, t.total));
    PRINT(v >= 1, "  Other calculations:       " << print_time(t.other(), t.total));
    PRINT(v >= 1, "------------------------------------------------------------");
    PRINT(v >= 1, "Total time:                 " << print_time(t.total, t.total) << std::left);
    PRINT(v >= 1, "------------------------------------------------------------\n");
}

/*
 * Returns time derivative approximation and its corresponding derivative w.r.t. xk
 */
inline double time_derivative_approx(eivec &dxdt, const rvec &xk, const SolverState &state, const std::size_t size)
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

/*
 * Integrates the system of DAEs in the interval `t = [0; t_end]`. TODO: Update description.
 * Mass, RHS, Jacobian will be moved from temporary or if called with std::move()
 * Returns solution `x` at time `t_end` (the initial condition stored in `x` will be overwritten).
 * If the solver stops earlier, `t_end` will be overwritten with the actual time.
 *
 * Parameters:
 *     x - initial condition (state_vector)
 *     t_end - integration interval `t = [0; t_end]` (double)
 *     t_output - (optional) a vector of output times (std::vector<double>)
 *
 * Returns:
 *     0 if integration is successful or error code if integration is failed (int)
 */
template <class Mass, class RHS, class Jacobian>
int solve(Mass mass, RHS rhs, Jacobian jac, const state_vector &x, const double t_end, const std::vector<double> t_output, const SolverOptions &opt, bool is_jac_auto)
{
    // Specific counters
    Counters c;

    // Specific timers
    Time time;

    // Global timer
    {
        // Measures total time
        Timer timer_global(&time.total);

        // Measures initialization time
        Timer *timer_init = new Timer(&time.init);

        // Initial output
        PRINT(opt.verbosity >= 1, "Starting dae-cpp solver...");
        PRINT((opt.verbosity >= 1) && is_jac_auto, "NOTE: Using automatic Jacobian...");

        // Vector of output times
        std::vector<double> t_out = std::move(t_output);

        // Sort vector of output times and erase duplicates
        t_out.push_back(t_end);
        std::sort(t_out.begin(), t_out.end());
        t_out.erase(std::unique(t_out.begin(), t_out.end()), t_out.end());

        // Throw an error if target time t < 0
        ASSERT(t_out.back() >= 0.0, "Target time t_end cannot be negative. The solver integrates from 0 to t_end.");

        // Check user-defined solver options
        opt.check();

        // Time step amplification threshold
        unsigned int dt_increase_threshold = 2 * (opt.Newton_scheme + 1) + opt.dt_increase_threshold_delta;
        ASSERT(dt_increase_threshold > 0, "Too small delta `dt_increase_threshold_delta`: " << opt.dt_increase_threshold_delta);

        // Time step reduction threshold
        unsigned int dt_decrease_threshold = 4 * (opt.Newton_scheme + 1) + opt.dt_decrease_threshold_delta;
        ASSERT(dt_decrease_threshold > 0, "Too small delta `dt_decrease_threshold_delta`: " << opt.dt_decrease_threshold_delta);
        ASSERT(dt_decrease_threshold > dt_increase_threshold, "Adaptive time stepping thresholds are not consistent with each other.");

        // Maximum number of iteration per time step
        unsigned int max_Newton_iter = opt.max_Jacobian_updates * (opt.Newton_scheme + 1);

        // System size
        auto size = x.size();
        ASSERT(size > 0, "Initial condition vector x is empty.");

        // Solver state
        SolverState state(size);

        // Alias for the current time step
        double &dt = state.dt[0];

        // Set initial time step
        dt = opt.dt_init;

        // Copy initial state
        state.x[0] = x;

        // Solution vector used for Newton iterations
        state_vector xk = x;

        // The RHS vector
        state_vector f(size);

        // Mass matrix holder
        sparse_matrix M;

        // Jacobian matrix holder
        sparse_matrix J;

        // Time derivative approximation
        eivec dxdt(size);

        // Linear solver
        Eigen::SparseLU<eimat> linsolver;

        // Eigen::SparseMatrix<float_type> matrices
        eimat M_; // Mass matrix (converted)
        eimat Jb; // Linear system matrix

        // Eigen::VectorX vectors
        eivec f_(size); // The RHS vector (converted)
        eivec b;        // The RHS of the linear system
        eivec dx;       // Linear system solution

        // Counts number of time steps
        u_int64_t n_steps{0};

        // Counts how many times the Newton iterator failed to converge in a row
        u_int32_t n_iter_failed{0};

        // Output after initialization
        PRINT(opt.verbosity >= 2, "Float size:      " << 8 * sizeof(float_type) << " bit");
        PRINT(opt.verbosity >= 2, "Integer size:    " << 8 * sizeof(int_type) << " bit");
        PRINT(opt.verbosity >= 2, "BDF max order:   " << opt.BDF_order);
        PRINT(opt.verbosity >= 2, "Newton scheme:   " << opt.Newton_scheme);
        PRINT(opt.verbosity >= 2, "Max time step:   " << opt.dt_max);
        PRINT(opt.verbosity >= 1, "DAE system size: " << size << " equations");
        PRINT(opt.verbosity >= 1, "Calculating...");

        // End of initialization. Stop the timer.
        delete timer_init;

        /*
         * Output time loop
         */
        for (const auto &t1 : t_out)
        {
            PRINT(opt.verbosity >= 2, "\n-- Integration time t = " << t1 << ":");

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

                if (opt.verbosity >= 2)
                {
                    std::cout << std::left
                              << "Step " << std::setw(8) << n_steps
                              << " :: t = " << std::setw(11) << state.t
                              << " :: ";
                }

                if (opt.verbosity >= 2)
                {
                    std::cout << "BDF-" << state.order
                              << ": dt=" << std::setw(11) << dt
                              << " :: ";
                }

                bool is_diverged{false}; // True if Newton iterations diverged

                unsigned int iter{}; // Newton iteration loop index - we will need this value later

                /*
                 * Newton iteration loop
                 */
                for (iter = 0; iter < max_Newton_iter; ++iter)
                {
                    double alpha{}; // Derivative w.r.t. xk

                    // Returns time derivative approximation dxdt and its corresponding derivative w.r.t. xk
                    {
                        Timer timer(&time.time_derivative);
                        alpha = time_derivative_approx(dxdt, xk, state, size);
                    }

                    // Get and convert the RHS
                    {
                        Timer timer(&time.rhs);

                        state_type f__(size), xk__(size);
                        for (std::size_t k = 0; k < size; ++k)
                        {
                            xk__[k] = xk[k];
                        }
                        rhs(f__, xk__, state.t);
                        ASSERT(f.size() == size, "The RHS vector size (" << f.size() << ") does not match the initial condition vector size (" << size << ").");
                        // f_ = Eigen::Map<eivec, Eigen::Unaligned>(f.data(), f.size());
                        for (std::size_t k = 0; k < size; ++k)
                        {
                            f_[k] = f__[k].val();
                        }
                    }

                    // Get and convert the Mass matrix
                    if (!opt.is_mass_matrix_static || (opt.is_mass_matrix_static && !iter && !c.n_lin_calls))
                    {
                        Timer timer(&time.mass);
                        M.clear();
                        mass(M, state.t);
                        M.check();
                        M_ = M.convert(size);
                    }

                    bool is_fact_enabled = !(iter % (opt.Newton_scheme + 1));

                    // Get and convert the Jacobian matrix
                    if (is_fact_enabled)
                    {
                        Timer timer(&time.jacobian);
                        J.clear();
                        jac(J, xk, state.t);
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
                        if (is_fact_enabled)
                        {
                            Jb -= M_ * alpha;
                        }
                    }

                    // Factorization
                    if (is_fact_enabled)
                    {
                        Timer timer(&time.factorization);

                        linsolver.compute(Jb);

                        if (linsolver.info() != Eigen::Success)
                        {
                            ERROR("Decomposition failed."); // TODO: Try to recover
                        }

                        c.n_fact_calls++;
                    }

                    // Solve linear system Jb dx = b
                    {
                        Timer timer(&time.linear_solver);

                        dx = linsolver.solve(b);

                        if (linsolver.info() != Eigen::Success)
                        {
                            ERROR("Solving failed."); // TODO: Try to recover
                        }

                        c.n_lin_calls++;

                        if (is_fact_enabled)
                        {
                            print_char(opt.verbosity >= 2, '#');
                        }
                        else
                        {
                            print_char(opt.verbosity >= 2, '*');
                        }
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
                            if (err_abs > opt.max_err_abs || std::isnan(dx[i]))
                            {
                                PRINT(opt.verbosity >= 2, " <- diverged");

                                // Trying to roll back and reduce the time step
                                is_diverged = true;
                                n_iter_failed++;
                                state.t = state.t_prev;
                                n_steps--;
                                dt /= opt.dt_decrease_factor;
                                if (dt < opt.dt_min)
                                {
                                    PRINT(opt.verbosity >= 1, "The time step was reduced to `t_min` but the scheme failed to converge.");
                                    goto result; // Abort all loops and go straight to the results.
                                                 // Using goto here is much more clear than using a sequence of `break` statements.
                                }

                                break;
                            }

                            // Relative error check
                            if (state.x[0][i] != 0.0)
                            {
                                double err_rel = err_abs / std::abs(state.x[0][i]);
                                if (err_rel > opt.rtol)
                                {
                                    is_converged = false;
                                }
                            }

                            // Absolute error check
                            if (err_abs > opt.atol)
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
                if (iter >= dt_decrease_threshold)
                {
                    dt /= opt.dt_decrease_factor;
                    print_char(opt.verbosity >= 2, '<');
                    if (dt < opt.dt_min)
                    {
                        PRINT(opt.verbosity >= 2, " <- reached dt_min");
                        PRINT(opt.verbosity >= 1, "The time step was reduced to `t_min` but the scheme failed to converge.");
                        goto result; // Abort all loops and go straight to the results.
                                     // Using goto here is much more clear than using a sequence of `break` statements.
                    }
                }
                if (iter <= dt_increase_threshold - 1)
                {
                    dt *= opt.dt_increase_factor;
                    if (dt > opt.dt_max)
                    {
                        dt = opt.dt_max;
                        print_char(opt.verbosity >= 2, '|');
                    }
                    else
                    {
                        print_char(opt.verbosity >= 2, '>');
                    }
                }

                // Adjust the last time step if needed
                if (dt > t1 - state.t)
                {
                    dt = t1 - state.t;
                }

                // Updates time integration order
                if (state.order < opt.BDF_order)
                {
                    state.order++;
                }

                // Newton iteration finished
                print_char(opt.verbosity >= 2, '\n');

                // We may already reached the target time
                if (dt < DAECPP_TIMESTEP_ROUNDING_ERROR)
                {
                    if (state.dt[2] > DAECPP_TIMESTEP_ROUNDING_ERROR)
                    {
                        dt = state.dt[2];
                    }
                    else if (state.dt[1] > DAECPP_TIMESTEP_ROUNDING_ERROR)
                    {
                        dt = state.dt[1];
                    }
                    else
                    {
                        dt = opt.dt_init;
                    }

                    break;
                }

            } // Time loop

        } // for (const auto &t1 : t_out) - Loop over all output times

    result:
        // Return solution and the final time
        // x = state.x[0];
        // t_end = state.t;

        PRINT(opt.verbosity >= 1, "...done");

    } // Global timer

    // Final output
    finalize(time, opt.verbosity, c);

    // Success
    return 0;
}

} // namespace internal
} // namespace core

template <class Mass, class RHS, class Jacobian>
int solve(Mass mass, RHS rhs, Jacobian jac, const state_vector &x, const double t_end, const std::vector<double> t_output = {}, const SolverOptions &opt = SolverOptions())
{
    return core::internal::solve(mass, rhs, jac, x, t_end, t_output, opt, false);
}

template <class Mass, class RHS>
int solve(Mass mass, RHS rhs, const state_vector &x, const double t_end, const std::vector<double> t_output = {}, const SolverOptions &opt = SolverOptions())
{
    return core::internal::solve(mass, rhs, JacobianAutomatic(rhs), x, t_end, t_output, opt, true);
}

template <class Mass, class RHS>
int solve(Mass mass, RHS rhs, const state_vector &x, const double t_end, const SolverOptions &opt = SolverOptions())
{
    return core::internal::solve(mass, rhs, JacobianAutomatic(rhs), x, t_end, {}, opt, true);
}

template <class Mass, class RHS, class Jacobian>
int solve(Mass mass, RHS rhs, Jacobian jac, const state_vector &x, const double t_end, const SolverOptions &opt = SolverOptions())
{
    return core::internal::solve(mass, rhs, jac, x, t_end, {}, opt, false);
}

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_H
