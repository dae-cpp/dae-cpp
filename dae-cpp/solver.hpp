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

    core::SolverState _state; // Solver state

    std::vector<double> _t_out; // Output times

    std::size_t _size{0};  // System size
    std::size_t _steps{0}; // Counts number of time steps
    std::size_t _calls{0}; // Counts total linear algebra solver calls

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
    int solve(state_type &x, double &t, const std::vector<double> t_output = {})
    {
        // Timer
        {
            // Measures total time
            Timer timer(&core::Timers::get().total_time);

            // Initial output
            PRINT(_opt.verbosity >= 1, "Starting dae-cpp solver...");

            // Initialize output times
            _t_out = std::move(t_output);

            // System size
            _size = x.size();
            ASSERT(_size > 0, "Initial condition vector is empty.");

            // Initial time
            _state.t[0] = 0.0;
            _state.t[1] = 0.0;

            // Initial time step
            _state.dt[0] = _opt.dt_init;
            _state.dt[1] = 0.0;

            // TODO: Reserve memory for the solution history

            // Copy initial state
            _state.x[0] = x;

            // TODO: Check user-defined solver options
            // _opt.check();

            // We don't need to do anything if t == 0. Throw an error if t < 0.
            if (t == 0.0)
            {
                NOTE("Target time t = 0. The initial condition is the solution.");
                return 0;
            }
            else if (t < 0.0)
            {
                ERROR("Target time t cannot be negative. The solver integrates from 0 to t.");
            }

            // TODO: Check initial time step is less than t

            // Full Jacobian matrix holder
            // sparse_matrix_holder J;

            // Full RHS vector
            // state_type b(m_size);

            // Solution vector used for Newton iterations
            // state_type xk(m_size);

            // Counts how many times the Newton iterator failed to converge within
            // max_Newton_iter iterations in a row.
            // int n_iter_failed = 0;

            // Output of initialization
            PRINT(_opt.verbosity >= 2, "Float size:      " << 8 * sizeof(float_type) << " bit");
            PRINT(_opt.verbosity >= 2, "Integer size:    " << 8 * sizeof(int_type) << " bit");
            PRINT(_opt.verbosity >= 1, "DAE system size: " << _size << " equations");
            PRINT(_opt.verbosity >= 1, "Calculating...");

            /*
             * Time loop
             */

            std::cout << "\nt_out:" << _t_out.size() << '\n';

            sparse_matrix J;
            _jac(J, x, t);
            J.check();

            sparse_matrix M;
            _mass(M, t);
            M.check();

            auto jac = J.convert(x.size());
            auto mass = M.convert(x.size());

            std::cout << "Mass:\n"
                      << mass.toDense() << '\n';
            std::cout << "Jac:\n"
                      << jac.toDense() << '\n';

            Eigen::SparseMatrix<double> sum(x.size(), x.size());

            sum = mass;
            sum += jac;
            // sum = mass + jac; // TODO: check performance

            std::cout << "Sum:\n"
                      << sum.toDense() << '\n';

            std::cout << "Print Jac:\n"
                      << J.dense(x.size()) << '\n';
        }

        std::cout << "Total time: " << core::Timers::get().total_time << '\n';
        return 0;
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
