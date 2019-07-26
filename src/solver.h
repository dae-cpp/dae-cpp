/*
 * The main solver class definition
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"
#include "jacobian.h"
#include "mass_matrix.h"
#include "solver_options.h"
#include "time_integrator.h"

namespace daecpp_namespace_name
{

class Solver
{
    RHS &m_rhs;  // RHS

    Jacobian &m_jac;  // Jacobian matrix

    MassMatrix &m_mass;  // Mass matrix

    SolverOptions &m_opt;  // Solver options

    TimeIntegrator *m_ti;  // Pointer to the time integrator

    struct m_iterator_state_struct  // Keeps the current time layer state
    {
        double t;                   // current time
        double dt[2];               // current and previous time steps
        double dt_eval;             // new time step
        int    current_scheme;      // current BDF order
        int    step_counter_local;  // local time step counter
        bool   final_time_step;     // do final time step
    } m_iterator_state;

    MKL_INT m_size;  // System size

    size_t m_steps = 0;  // Total time iteration counter
    size_t m_calls = 0;  // Total linear algebra solver calls counter

    // Count the number of the DAE solver calls (for output)
    size_t m_dae_solver_calls = 0;

    // Timers
    double m_timer_lin = 0;
    double m_timer_rhs = 0;
    double m_timer_jac = 0;
    double m_timer_tot = 0;

    // Contains a few latest successful time steps for the time integrator
    state_type_matrix m_x_prev;

    // Intel MKL PARDISO control parameters
    MKL_INT m_phase;        // Current phase of the solver
    MKL_INT m_maxfct = 1;   // Maximum number of numerical factorizations
    MKL_INT m_mnum   = 1;   // Which factorization to use
    MKL_INT m_mtype  = 11;  // Real unsymmetric matrix
    MKL_INT m_nrhs   = 1;   // Number of right hand sides
    MKL_INT m_msglvl = 0;   // Print statistical information
    MKL_INT m_error  = 0;   // Error flag

    // Intel MKL PARDISO sparse matrix indeces
    MKL_INT *m_ia = nullptr;
    MKL_INT *m_ja = nullptr;

    // Intel MKL PARDISO vectors and sparse matrix non-zero elements
    float_type *m_mkl_a = nullptr;
    float_type *m_mkl_b = nullptr;
    float_type *m_mkl_x = nullptr;

    // Intel MKL PARDISO internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void *m_pt[64];

    // Intel MKL PARDISO auxiliary variables
    double  m_ddum;  // Double dummy
    MKL_INT m_idum;  // Integer dummy

    // Intel MKL PARDISO iparm parameter
    MKL_INT m_iparm[64];

    // Simple yet efficient Adaptive Time Stepping
    int m_adaptive_time_stepping(state_type &x, const state_type_matrix &x_prev,
                                 int iter);

    // Scrapes the current time iteration and decreases the time step
    // Return -1 in case the time step is below dt_min
    int m_reset_ti_state(state_type &x, const state_type_matrix &x_prev);

    // Updates time integrator scheme when the time step changes
    int m_reset_ti_scheme();

    // Increases the time step
    void m_increase_dt();

    // Decreases the time step
    void m_decrease_dt();

    // Checks if dt is within the interval defined in solver_options.h
    int m_check_dt();

    // Checks PARDISO solver error messages
    void m_check_pardiso_error(MKL_INT err);

public:
    /*
     * Receives user-defined RHS, Jacobian, Mass matrix and solver options.
     * Defined in solver.cpp
     */
    Solver(RHS &rhs, Jacobian &jac, MassMatrix &mass, SolverOptions &opt);

    /*
     * Releases memory. Defined in solver.cpp.
     */
    virtual ~Solver();

    /*
     * Integrates the system of DAEs on the interval t = [t0; t1] and returns
     * result in the array x. Parameter t0 can be overriden in the solver
     * options (t0 = 0 by default).
     * The data stored in x (initial conditions) will be overwritten.
     * Returns 0 in case of success or error code if integration is failed.
     */
    int operator()(state_type &x, double &t1);

    /*
     * Virtual Observer. Called by the solver every time step.
     * Receives current solution vector and the current time.
     * Can be overriden by a user to get intermediate results (for example,
     * for plotting).
     */
    virtual void observer(state_type &x, const double t)
    {
        return;  // It does nothing by default
    }
};

}  // namespace daecpp_namespace_name
