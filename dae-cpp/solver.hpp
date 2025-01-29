/*
 * The main solver implementation.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#ifndef DAECPP_SOLVER_H
#define DAECPP_SOLVER_H

#include <algorithm>
#include <iomanip>

#include <Eigen/SparseLU>

#include "jacobian-matrix.hpp"
#include "mass-matrix.hpp"
#include "solution-manager.hpp"
#include "solver-options.hpp"
#include "timer.hpp"
#include "vector-function.hpp"
#include "version.hpp"

namespace daecpp_namespace_name
{
namespace core
{
namespace detail
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
    explicit SolverState(const std::size_t size)
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
    uint64_t n_fact_calls{0}; // Number of Jacobian matrix updates and factorizations
    uint64_t n_lin_calls{0};  // Number of linear solver calls
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
    double t_coef{1.0};       // Default conversion coefficient for output
    std::string t_unit{"ms"}; // Default time units for output

    if (t_total > 1e4)
    {
        t_coef = 1e-3; // ms -> s
        t_unit = "s";  // seconds
    }

    std::ostringstream oss;
    oss << std::setw(11) << val * t_coef << ' ' << t_unit;
    oss << std::setprecision(3);
    oss << "  (" << (int)(val / t_total * 1e5) / 1000.0 << "%)";
    std::string var = oss.str();
    return var;
}

/*
 * Final output
 */
inline void finalize(const timer::Time &time, const int v, const Counters c)
{
    auto &t = time.timers;

    PRINT(v >= 1, "\nComputation time:" << std::right);
    PRINT(v >= 1, std::string(64, '-'));
    PRINT(v >= 1, "  Initialization:       " << print_time(t[timer::init], time.total));
    PRINT(v >= 1, "  Time derivative:      " << print_time(t[timer::time_derivative], time.total));
    PRINT(v >= 1, "  RHS:                  " << print_time(t[timer::rhs], time.total));
    PRINT(v >= 1, "  Mass matrix:          " << print_time(t[timer::mass], time.total));
    PRINT(v >= 1, "  Jacobian matrix:      " << print_time(t[timer::jacobian], time.total));
    PRINT(v >= 1, "  Linear algebra:       " << print_time(t[timer::linear_algebra], time.total));
    PRINT(v >= 1, "  Matrix factorization: " << print_time(t[timer::factorization], time.total) << "  <--  " << c.n_fact_calls << " calls");
    PRINT(v >= 1, "  Linear solver:        " << print_time(t[timer::linear_solver], time.total) << "  <--  " << c.n_lin_calls << " calls");
    PRINT(v >= 1, "  Error control:        " << print_time(t[timer::error_check], time.total));
    PRINT(v >= 1, "  Solution Manager:     " << print_time(t[timer::manager], time.total));
    PRINT(v >= 1, "  Other calculations:   " << print_time(time.other(), time.total));
    PRINT(v >= 1, std::string(64, '-'));
    PRINT(v >= 1, "Total time:             " << print_time(time.total, time.total) << std::left);
    PRINT(v >= 1, std::string(64, '-') << '\n');
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
    {
        alpha = (4.0 * h0 * h0 * h0 +
                 h1 * h12 * h123 +
                 3.0 * h0 * h0 * (3.0 * h1 + 2.0 * h2 + h3) +
                 2.0 * h0 * (3.0 * h1 * h1 + h2 * h23 + 2.0 * h1 * (2.0 * h2 + h3))) /
                (h0 * h01 * h012 * h0123);

        const double term0 = h01 * h012 * h0123 / (h0 * h1 * h12 * h123);
        const double term1 = h0 * h012 * h0123 / (h1 * h01 * h2 * h23);
        const double term2 = h0 * h01 * h0123 / (h2 * h12 * h012 * h3);
        const double term3 = h0 * h01 * h012 / (h3 * h23 * h123 * h0123);

        for (std::size_t i = 0; i < size; ++i)
        {
            dxdt[i] = alpha * xk[i] -
                      term0 * state.x[0][i] +
                      term1 * state.x[1][i] -
                      term2 * state.x[2][i] +
                      term3 * state.x[3][i];
        }
    }
    break;

    default:
        ERROR("Unsupported time integration order.");
    }

    return alpha;
}

/*
 * The main solver.
 * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
 *
 * Parameters:
 *     `mass` - Mass matrix (Mass matrix object)
 *     `rhs` - the Right-Hand Side (vector function) of the DAE system (Vector function object)
 *     `jac` - Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
 *     `mgr` - Solution Manager object
 *     `x0` - initial condition (`state_vector`)
 *     `t_end` - integration interval `t = [0; t_end]` (`double`)
 *     `t_output` - a vector of output times (`std::vector<double>`)
 *     `opt` - solver options (`SolverOptions` object)
 *     `is_jac_auto` - `true` if Jacobian is computed automatically, `false` otherwise (`bool`)
 *
 * Returns:
 *     `daecpp::exit_code::success` if integration is successful or error code if integration is failed (`int`)
 */
template <class Mass, class RHS, class Jacobian, class Manager>
inline exit_code::status solve(Mass mass, RHS rhs, Jacobian jac, Manager mgr, const state_vector &x0, const double t_end, const std::vector<double> &t_output, const SolverOptions &opt, bool is_jac_auto)
{
    // Specific counters
    Counters c;

    // Specific timers
    timer::Time time;

    // An alias for the specific timers array
    auto &t = time.timers;

    // Solution outcome (success or error code)
    exit_code::status error_msg{exit_code::unknown};

    // Global timer
    {
        // Measures total time
        Timer timer_global(&time.total);

        // Measures initialization time
        Timer *timer_init = new Timer(&t[timer::init]);

        // Initial output
        PRINT(opt.verbosity >= 1, "Starting dae-cpp version " << version_major << '.' << version_minor << '.' << version_patch << "...");
        PRINT((opt.verbosity >= 1) && is_jac_auto, "NOTE: Using automatic Jacobian...");

        // A copy of the vector of output times
        std::vector<double> t_out = t_output;

        // Sort vector of output times and erase duplicates
        if (t_out.size())
        {
            std::sort(t_out.begin(), t_out.end());
            t_out.erase(std::unique(t_out.begin(), t_out.end()), t_out.end());
        }
        else
        {
            t_out.push_back(t_end);
        }

        // Throw an error if target time t < 0
        ASSERT(t_out.back() >= 0.0, "Target time t_end cannot be negative. The solver integrates from 0 to t_end.");

        // Check user-defined solver options
        opt.check();

        // Time step amplification threshold
        int dt_increase_threshold = 2 * (opt.Newton_scheme + 1) + opt.dt_increase_threshold_delta;
        ASSERT(dt_increase_threshold > 0, "Too small delta `dt_increase_threshold_delta`: " << opt.dt_increase_threshold_delta);

        // Time step reduction threshold
        int dt_decrease_threshold = 4 * (opt.Newton_scheme + 1) + opt.dt_decrease_threshold_delta;
        ASSERT(dt_decrease_threshold > 0, "Too small delta `dt_decrease_threshold_delta`: " << opt.dt_decrease_threshold_delta);
        ASSERT(dt_decrease_threshold > dt_increase_threshold, "Adaptive time stepping thresholds are not consistent with each other.");

        // Maximum number of iteration per time step
        int max_Newton_iter = opt.max_Jacobian_updates * (opt.Newton_scheme + 1);

        // System size
        auto size = x0.size();
        ASSERT(size > 0, "Initial condition vector x is empty.");

        // Solver state
        SolverState state(size);

        // Alias for the current time step
        double &dt = state.dt[0];

        // Set initial time step
        dt = opt.dt_init;

        // Copy initial state
        state.x[0] = x0;

        // Solution vector used for Newton iterations
        state_vector xk = x0;

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
        uint64_t n_steps{0};

        // Counts how many times the Newton iterator failed to converge in a row
        uint32_t n_iter_failed{0};

        // Output after initialization
        PRINT(opt.verbosity >= 2, "Float size:      " << 8 * sizeof(float_type) << " bit");
        PRINT(opt.verbosity >= 2, "Integer size:    " << 8 * sizeof(int_type) << " bit");
        PRINT(opt.verbosity >= 2, "BDF max order:   " << opt.BDF_order);
        PRINT(opt.verbosity >= 2, "Newton scheme:   " << opt.Newton_scheme);
        PRINT(opt.verbosity >= 2, "Max time step:   " << opt.dt_max);
        PRINT(opt.verbosity >= 1, "DAE system size: " << size << " equations");
        PRINT(opt.verbosity >= 1, "Calculating...");

        // Call Solution Manager functor with the initial condition
        try
        {
            Timer timer(&t[timer::manager]);
            if (mgr(x0, 0.0) == solver_command::stop_intergration)
            {
                PRINT(opt.verbosity >= 1, "Stop event in Solution Manager triggered.");
                error_msg = exit_code::success;
                goto result;
            }
        }
        catch (const std::exception &e)
        {
            ERROR("Solution Manager functor call failed.\n"
                  << e.what());
        }

        // End of initialization. Stop the timer.
        delete timer_init;

        /*
         * Output time loop
         */
        for (const auto &t1 : t_out)
        {
            PRINT(opt.verbosity >= 2, "\n-- Integration time t = " << t1 << ":");

            if (t1 <= 0.0)
            {
                WARNING("Negative or zero integration time t = " << t1 << ". Skipped.");
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

            bool delay_timestep_inc{false}; // If true, delays increasing the time step for one time step

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

                int iter{}; // Newton iteration loop index - we will need this value later

                /*
                 * Newton iteration loop
                 */
                for (iter = 0; iter < max_Newton_iter; ++iter)
                {
                    double alpha{}; // Derivative w.r.t. xk

                    // Returns time derivative approximation dxdt and its corresponding derivative w.r.t. xk
                    try
                    {
                        Timer timer(&t[timer::time_derivative]);
                        alpha = time_derivative_approx(dxdt, xk, state, size);
                    }
                    catch (const std::exception &e)
                    {
                        ERROR("Failed to compute the time derivative approximation.\n"
                              << e.what());
                    }

                    // Get and convert the RHS
                    try
                    {
                        Timer timer(&t[timer::rhs]);

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
                    catch (const std::exception &e)
                    {
                        ERROR("Failed to compute and convert the vector function (RHS).\n"
                              << e.what());
                    }

                    // Get and convert the Mass matrix
                    if (!opt.is_mass_matrix_static || (opt.is_mass_matrix_static && !iter && !c.n_lin_calls))
                    {
                        try
                        {
                            Timer timer(&t[timer::mass]);
                            M.clear();
                            mass(M, state.t);
                            M.check();
                            M_ = M.convert(static_cast<int_type>(size));
                        }
                        catch (const std::exception &e)
                        {
                            ERROR("Failed to compute and convert the Mass matrix.\n"
                                  << e.what());
                        }
                    }

                    // Enable/disable factorization depending on the number of failed attempts and the current Newton scheme
                    bool is_fact_enabled = (n_iter_failed < opt.max_Newton_failed_attempts) ? !(iter % (opt.Newton_scheme + 1)) : true;

                    // Get and convert the Jacobian matrix
                    if (is_fact_enabled)
                    {
                        try
                        {
                            Timer timer(&t[timer::jacobian]);
                            J.clear();
                            jac(J, xk, state.t);
                            J.check();
                            Jb = J.convert(static_cast<int_type>(size));
                        }
                        catch (const std::exception &e)
                        {
                            ERROR("Failed to compute and convert the Jacobian matrix.\n"
                                  << e.what());
                        }
                    }

                    // Matrix-vector operations
                    try
                    {
                        Timer timer(&t[timer::linear_algebra]);

                        // b = M(t) * [dx/dt] - f
                        b = M_ * dxdt;
                        b -= f_;

                        // Jb = J - d/dxk (M(t) * [dx/dt])
                        if (is_fact_enabled)
                        {
                            Jb -= M_ * alpha;
                        }
                    }
                    catch (const std::exception &e)
                    {
                        ERROR("Failed to perform matrix-vector operations with Eigen.\n"
                              << e.what());
                    }

                    // Factorization
                    if (is_fact_enabled)
                    {
                        Timer timer(&t[timer::factorization]);

                        linsolver.compute(Jb);

                        c.n_fact_calls++;

                        if (linsolver.info() != Eigen::Success)
                        {
                            PRINT(opt.verbosity >= 2, " <- decomposition failed");
                            if(opt.recover_from_linsolver_failure)
                            {
                                is_diverged = true;
                                break;
                            }
                            error_msg = exit_code::linsolver_failed_decomposition;
                            goto result; // Abort all loops and go straight to the results
                        }
                    }

                    // Solve linear system Jb dx = b
                    {
                        Timer timer(&t[timer::linear_solver]);

                        dx = linsolver.solve(b);

                        c.n_lin_calls++;

                        if (linsolver.info() != Eigen::Success)
                        {
                            PRINT(opt.verbosity >= 2, " <- linear solver failed");
                            if(opt.recover_from_linsolver_failure)
                            {
                                is_diverged = true;
                                break;
                            }
                            error_msg = exit_code::linsolver_failed_solving;
                            goto result; // Abort all loops and go straight to the results
                        }

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
                        Timer timer(&t[timer::error_check]);

                        for (std::size_t i = 0; i < size; ++i)
                        {
                            // Absolute error
                            auto err_abs = std::abs(dx[i]);

                            // Solution diverged. Roll back to the previous state and redo with reduced time step.
                            if (err_abs > opt.max_err_abs || std::isnan(dx[i]))
                            {
                                PRINT(opt.verbosity >= 2, " <- diverged");
                                is_diverged = true;
                                break;
                            }

                            // Relative error check
                            auto x_abs = std::abs(state.x[0][i]);
                            if ((x_abs > DAECPP_FLOAT_TOLERANCE) && ((err_abs / x_abs) > opt.rtol))
                            {
                                is_converged = false;
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

                if (iter == max_Newton_iter)
                {
                    PRINT(opt.verbosity >= 2, " <- couldn't converge");
                    is_diverged = true;
                }

                if (is_diverged)
                {
                    // Trying to roll back and reduce the time step
                    n_iter_failed++;
                    delay_timestep_inc = true;
                    state.t = state.t_prev;
                    n_steps--;
                    dt /= opt.dt_decrease_factor;
                    xk = state.x[0];
                    if (dt < opt.dt_min)
                    {
                        PRINT(opt.verbosity >= 1, "The time step was reduced to `t_min` but the scheme failed to converge.");
                        error_msg = exit_code::diverged;
                        goto result; // Abort all loops and go straight to the results
                    }
                    continue;
                }

                bool decrease_time_step{false}; // If true, decrease the time step

                // Call Solution Manager functor with the current solution and time
                try
                {
                    Timer timer(&t[timer::manager]);
                    auto command = mgr(xk, state.t);
                    if (command == solver_command::decrease_time_step)
                    {
                        decrease_time_step = true;
                    }
                    else if (command == solver_command::decrease_time_step_and_redo)
                    {
                        PRINT(opt.verbosity >= 2, " <- decrease_time_step_and_redo");
                        delay_timestep_inc = true;
                        state.t = state.t_prev;
                        n_steps--;
                        dt /= opt.dt_decrease_factor;
                        xk = state.x[0];
                        if (dt < opt.dt_min)
                        {
                            PRINT(opt.verbosity >= 1, "The time step was reduced to `t_min` but the scheme failed to converge.");
                            error_msg = exit_code::diverged;
                            goto result; // Abort all loops and go straight to the results
                        }
                        continue;
                    }
                    else if (command) // solver_command::stop_intergration
                    {
                        print_char(opt.verbosity >= 2, '\n');
                        PRINT(opt.verbosity >= 1, "Stop event in Solution Manager triggered.");
                        error_msg = exit_code::success;
                        goto result;
                    }
                }
                catch (const std::exception &e)
                {
                    ERROR("Solution Manager functor call failed.\n"
                          << e.what());
                }

                n_iter_failed = 0;

                double variability{0.0}; // Maximum relative variability of the solution

                {
                    // Updates state history
                    for (std::size_t i = 0; i < size; ++i)
                    {
                        // Manually unrolled loop
                        state.x[3][i] = state.x[2][i];
                        state.x[2][i] = state.x[1][i];
                        state.x[1][i] = state.x[0][i];
                        state.x[0][i] = xk[i];

                        // Finds maximum relative variability
                        if (opt.solution_variability_control)
                        {
                            if ((std::abs(state.x[1][i]) > opt.variability_tolerance) &&
                                (std::abs(state.x[0][i]) > opt.variability_tolerance))
                            {
                                double rel_change = std::abs((state.x[1][i] - state.x[0][i]) / state.x[1][i]);
                                if (rel_change > variability)
                                {
                                    variability = rel_change;
                                }
                            }
                        }
                    }

                    // Updates time step history
                    for (int k = DAECPP_MAX_ORDER - 1; k > 0; --k)
                    {
                        state.dt[k] = state.dt[k - 1];
                    }
                }

                // Make decision about new time step
                if ((iter >= dt_decrease_threshold) ||
                    (variability > opt.variability_threshold_high) ||
                    decrease_time_step)
                {
                    dt /= opt.dt_decrease_factor;
                    decrease_time_step = false;
                    print_char(opt.verbosity >= 2, '<');
                    if (dt < opt.dt_min)
                    {
                        PRINT(opt.verbosity >= 2, " <- reached dt_min");
                        PRINT(opt.verbosity >= 1, "The time step was reduced to `t_min` but the scheme failed to converge.");
                        error_msg = exit_code::diverged;
                        goto result; // Abort all loops and go straight to the results
                    }
                }
                else if ((iter <= dt_increase_threshold - 1) &&
                    !delay_timestep_inc &&
                    (variability <= opt.variability_threshold_low))
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

                delay_timestep_inc = false;

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
                if (dt < opt.dt_min)
                {
                    if (state.dt[2] > opt.dt_min)
                    {
                        dt = state.dt[2];
                    }
                    else if (state.dt[1] > opt.dt_min)
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

        error_msg = exit_code::success;

    result: // Using goto here is much more clear than using a sequence of `break` statements

        PRINT(opt.verbosity >= 1, "...done");

    } // Global timer

    // Final output
    finalize(time, opt.verbosity, c);

    // Success
    return error_msg;
}

} // namespace detail
} // namespace core

/*
 * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
 *
 * Parameters:
 *     `mass` - Mass matrix (Mass matrix object)
 *     `rhs` - the Right-Hand Side (vector function) of the DAE system (Vector function object)
 *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
 *     `x0` - initial condition (`state_vector`)
 *     `t_end` - integration interval `t = [0; t_end]` (`double`)
 *     `mgr` - Solution Manager object
 *     `opt` - (optional) solver options (`SolverOptions` object)
 *
 * Returns:
 *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
 */
template <class Mass, class RHS, class Jacobian, class Manager = SolutionManager>
inline exit_code::status solve(Mass mass, RHS rhs, Jacobian jac, const state_vector &x0, const double t_end, Manager mgr = SolutionManager(), const SolverOptions &opt = SolverOptions())
{
    return core::detail::solve(mass, rhs, jac, mgr, x0, t_end, {}, opt, false);
}

/*
 * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
 *
 * Parameters:
 *     `mass` - Mass matrix (Mass matrix object)
 *     `rhs` - the Right-Hand Side (vector function) of the DAE system (Vector function object)
 *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
 *     `x0` - initial condition (`state_vector`)
 *     `t_end` - integration interval `t = [0; t_end]` (`double`)
 *     `mgr` - Solution Manager object
 *     `opt` - (optional) solver options (`SolverOptions` object)
 *
 * Returns:
 *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
 */
template <class Mass, class RHS, class Manager = SolutionManager>
inline exit_code::status solve(Mass mass, RHS rhs, const state_vector &x0, const double t_end, Manager mgr = SolutionManager(), const SolverOptions &opt = SolverOptions())
{
    return core::detail::solve(mass, rhs, JacobianAutomatic(rhs), mgr, x0, t_end, {}, opt, true);
}

/*
 * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
 *
 * Parameters:
 *     `mass` - Mass matrix (Mass matrix object)
 *     `rhs` - the Right-Hand Side (vector function) of the DAE system (Vector function object)
 *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
 *     `x0` - initial condition (`state_vector`)
 *     `t_output` - a vector of output times (`std::vector<double>`)
 *     `mgr` - Solution Manager object
 *     `opt` - (optional) solver options (`SolverOptions` object)
 *
 * Returns:
 *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
 */
template <class Mass, class RHS, class Jacobian, class Manager = SolutionManager>
inline exit_code::status solve(Mass mass, RHS rhs, Jacobian jac, const state_vector &x0, const std::vector<double> &t_output, Manager mgr = SolutionManager(), const SolverOptions &opt = SolverOptions())
{
    return core::detail::solve(mass, rhs, jac, mgr, x0, 0.0, t_output, opt, false);
}

/*
 * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
 *
 * Parameters:
 *     `mass` - Mass matrix (Mass matrix object)
 *     `rhs` - the Right-Hand Side (vector function) of the DAE system (Vector function object)
 *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
 *     `x0` - initial condition (`state_vector`)
 *     `t_output` - a vector of output times (`std::vector<double>`)
 *     `mgr` - Solution Manager object
 *     `opt` - (optional) solver options (`SolverOptions` object)
 *
 * Returns:
 *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
 */
template <class Mass, class RHS, class Manager = SolutionManager>
inline exit_code::status solve(Mass mass, RHS rhs, const state_vector &x0, const std::vector<double> &t_output, Manager mgr = SolutionManager(), const SolverOptions &opt = SolverOptions())
{
    return core::detail::solve(mass, rhs, JacobianAutomatic(rhs), mgr, x0, 0.0, t_output, opt, true);
}

/*
 * DAE System class.
 * Serves as a wrapper for lower level `solve(...)` function calls.
 * Contains Solver Options object and Solution Holder object.
 */
template <class Mass, class RHS>
class System
{
    Mass m_mass; // Mass matrix holder
    RHS m_rhs;   // The system RHS holder

public:
    SolverOptions opt = SolverOptions(); // Solver options

    SolutionHolder sol; // Solution holder

    exit_code::status status{exit_code::unknown}; // Solver status (exit code)

    /*
     * DAE System class serves as a wrapper for lower level `solve(...)` function calls.
     * Contains Solver Options object and Solution Holder object.
     *
     * Parameters:
     *     `mass` - Mass matrix (Mass matrix object)
     *     `rhs` - the Right-Hand Side (vector function) of the DAE system (Vector function object)
     */
    System(Mass mass, RHS rhs) : m_mass(mass), m_rhs(rhs) {}

    /*
     * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
     *
     * Parameters:
     *     `x0` - initial condition (`state_vector`)
     *     `t_end` - integration interval `t = [0; t_end]` (`double`)
     *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
     *
     * Returns:
     *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
     */
    template <class Jacobian>
    exit_code::status solve(const state_vector &x0, const double t_end, Jacobian jac)
    {
        status = core::detail::solve(m_mass, m_rhs, jac, Solution(sol), x0, t_end, {}, opt, false);
        return status;
    }

    /*
     * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
     *
     * Parameters:
     *     `x0` - initial condition (`state_vector`)
     *     `t_end` - integration interval `t = [0; t_end]` (`double`)
     *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
     *
     * Returns:
     *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
     */
    exit_code::status solve(const state_vector &x0, const double t_end)
    {
        status = core::detail::solve(m_mass, m_rhs, JacobianAutomatic(m_rhs), Solution(sol), x0, t_end, {}, opt, true);
        return status;
    }

    /*
     * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
     *
     * Parameters:
     *     `x0` - initial condition (`state_vector`)
     *     `t_output` - a vector of output times (`std::vector<double>`)
     *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
     *
     * Returns:
     *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
     */
    template <class Jacobian>
    exit_code::status solve(const state_vector &x0, const std::vector<double> &t_output, Jacobian jac)
    {
        status = core::detail::solve(m_mass, m_rhs, jac, Solution(sol, t_output), x0, 0.0, t_output, opt, false);
        return status;
    }

    /*
     * Integrates the system of DAEs in the interval `t = [0; t_end]` with the initial condition `x0`.
     *
     * Parameters:
     *     `x0` - initial condition (`state_vector`)
     *     `t_output` - a vector of output times (`std::vector<double>`)
     *     `jac` - (optional) Jacobian matrix (matrix of the RHS derivatives) (Jacobian matrix object)
     *
     * Returns:
     *     `daecpp::exit_code::success` (0) if integration is successful or error code if integration is failed (`int`)
     */
    exit_code::status solve(const state_vector &x0, const std::vector<double> &t_output)
    {
        status = core::detail::solve(m_mass, m_rhs, JacobianAutomatic(m_rhs), Solution(sol, t_output), x0, 0.0, t_output, opt, true);
        return status;
    }
};

} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_H
