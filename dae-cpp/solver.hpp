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
#include "timer.hpp"

namespace daecpp_namespace_name
{

class Solver
{
    const RHS &_rhs;         // RHS
    const Jacobian &_jac;    // Jacobian matrix
    const MassMatrix &_mass; // Mass matrix

    // TimeIntegrator *m_ti;  // Pointer to the time integrator

    // struct m_iterator_state_struct  // Keeps the current time layer state
    // {
    //     double t;                   // current time
    //     double dt[2];               // current and previous time steps
    //     double dt_eval;             // new time step
    //     int    current_scheme;      // current BDF order
    //     int    step_counter_local;  // local time step counter
    //     bool   final_time_step;     // do final time step
    // } m_iterator_state;

    // MKL_INT m_size;  // System size

    // std::size_t m_steps = 0;  // Total time iteration counter
    // std::size_t m_calls = 0;  // Total linear algebra solver calls counter

    // // Count the number of the DAE solver calls (for output)
    // std::size_t m_dae_solver_calls = 0;

    // // Timers
    // double m_timer_lin = 0;
    // double m_timer_rhs = 0;
    // double m_timer_jac = 0;
    // double m_timer_tot = 0;

    // // Contains a few latest successful time steps for the time integrator
    // state_type_matrix m_x_prev;

    // // Intel MKL PARDISO control parameters

    // // Simple yet efficient Adaptive Time Stepping
    // int m_adaptive_time_stepping(state_type &x, const state_type_matrix &x_prev,
    //                              int iter);

    // // Scrapes the current time iteration and decreases the time step
    // // Return -1 in case the time step is below dt_min
    // int m_reset_ti_state(state_type &x, const state_type_matrix &x_prev);

    // // Updates time integrator scheme when the time step changes
    // int m_reset_ti_scheme();

    // // Increases the time step
    // void m_increase_dt();

    // // Decreases the time step
    // void m_decrease_dt();

    // // Checks if dt is within the interval defined in solver_options.h
    // int m_check_dt();

    // // Checks PARDISO solver error messages
    // void m_check_pardiso_error(MKL_INT err);

    // protected:
    //     /*
    //      * Expose solver options for the observer in the children classes
    //      */
    //     SolverOptions &m_opt;  // Solver options

public:
    /*
     * Receives user-defined RHS, Jacobian, Mass matrix and solver options.
     * Defined in solver.cpp
     */
    Solver(const MassMatrix &mass, const RHS &rhs, const Jacobian &jac) : _mass(mass), _rhs(rhs), _jac(jac) {}
    // Solver(const MassMatrix &mass, const RHS &rhs, const Jacobian &jac, const SolverOptions &opt) : _mass(mass), _rhs(rhs), _jac(jac), _opt(opt)

    // /*
    //  * Releases memory. Defined in solver.cpp.
    //  */
    // virtual ~Solver();

    /*
     * Integrates the system of DAEs on the interval t = [t0; t1] and returns
     * result in the array x. Parameter t0 can be overriden in the solver
     * options (t0 = 0 by default).
     * The data stored in x (initial conditions) will be overwritten.
     * Returns 0 in case of success or error code if integration is failed.
     */
    int operator()(state_type &x, double &t)
    {
        // Timers::reset();
        // MySingleton::getInstance();
        {
            Timer timer(&core::Timers::get().total);

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

        std::cout << "Time: " << core::Timers::get().total << '\n';
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
