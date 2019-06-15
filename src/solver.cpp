/*
 * The main solver
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <mkl_types.h>
#include <mkl_pardiso.h>

#include "solver.h"
#include "time_integrator.h"

namespace daecpp_namespace_name
{

/*
 * The main solver
 * =============================================================================
 */
int Solver::operator()(state_type &x, const double t1)
{
    // Set system size
    m_size = (MKL_INT)(x.size());

    // Check user-defined solver options
    m_opt.check_options();

    // Assert t1 > t0
    if(t1 <= m_opt.t0)
    {
        std::cout << "ERROR: Integration time t1 = " << t1
                  << " cannot be less than the initial time t0 = " << m_opt.t0
                  << std::endl;
        return 1;
    }

    // Check the initial time step
    m_opt.dt_init =
        (m_opt.dt_init > (t1 - m_opt.t0)) ? (t1 - m_opt.t0) : m_opt.dt_init;

    // Initial time
    m_iterator_state.t = m_opt.t0;

    // Initial time step
    m_iterator_state.dt[0] = m_opt.dt_init;
    m_iterator_state.dt[1] = m_iterator_state.dt[0];

    // Initialise the time integrator state structure.
    // Solver starts the first time step using BDF-1 method
    // (since it doesn't have enough history yet)
    m_iterator_state.current_scheme     = 1;
    m_iterator_state.step_counter_local = 0;
    m_iterator_state.final_time_step    = false;

    // Initialise time integrator
    TimeIntegrator ti(m_rhs, m_jac, m_mass, m_opt, m_size);

    // Initial output
    if(m_opt.verbosity > 1)
    {
        std::cout << "Number of equations: " << m_size << std::endl;
        std::cout << "Float precision:     " << 8 * sizeof(float_type)
                  << " bit\n";
        std::cout << "Integer precision:   " << 8 * sizeof(MKL_INT) << " bit\n";
        std::cout << "Numerical algorithm: BDF-" << m_opt.bdf_order
                  << std::endl;
    }

    // Contains a few latest successful time steps for Time Integrator
    state_type_matrix x_prev(m_opt.bdf_order, state_type(m_size));

    // Full Jacobian matrix holder
    sparse_matrix_holder J;

    // Full RHS vector
    state_type b(m_size);

    // Solution vector used for Newton iterations
    state_type xk(m_size);

    // Copy current state vector into the history vector
    x_prev[0] = x;

    // Reset PARDISO pointers
    m_mkl_b = b.data();
    m_mkl_x = xk.data();

    // Load Intel MKL PARDISO iparm parameter from solver_options class
    m_opt.set_iparm_for_pardiso(m_iparm);

    // Memory control variables
    int peak_mem1 = 0, peak_mem2 = 0, peak_mem3 = 0;

    /*
     * Start the solver
     * =========================================================================
     */

    // TODO: Start timer here

    m_iterator_state.t += m_iterator_state.dt[0];

    while(m_iterator_state.t < (t1 + m_iterator_state.dt[0] * 0.5))
    {
        m_iterator_state.step_counter_local++;
        m_steps++;

        if(m_opt.verbosity > 0)
        {
            std::cout << std::left;
            std::cout << "\nStep " << std::setw(7) << m_steps
                      << " :: t = " << std::setw(12) << m_iterator_state.t
                      << " :: ";
            std::cout.flush();
        }

        if(m_opt.verbosity > 1)
        {
            std::cout << "BDF-" << m_iterator_state.current_scheme << ": ";
        }

        ti.set_scheme(m_iterator_state.current_scheme);

        if(m_iterator_state.current_scheme < m_opt.bdf_order)
        {
            m_iterator_state.current_scheme++;
        }

        int iter;  // We need this value later

        for(iter = 0; iter < m_opt.max_Newton_iter; iter++)
        {
            // Reordering, Symbolic and Numerical Factorization
            if(m_opt.fact_every_iter || iter == 0)
            {
                // Time Integrator with updated Jacobian
                ti(J, b, x, x_prev, m_iterator_state.t, m_iterator_state.dt, true);

                // Jacobian can change its size and can be re-allocated.
                // Catch up new array addresses.
                m_mkl_a = J.A.data();
                m_ia    = J.ia.data();
                m_ja    = J.ja.data();

                // PHASE 1.
                // Reordering and Symbolic Factorization. This step also
                // allocates all memory that is necessary for the factorization
                m_phase = 11;
                PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size,
                        m_mkl_a, m_ia, m_ja, &m_idum, &m_nrhs, m_iparm,
                        &m_msglvl, &m_ddum, &m_ddum, &m_error);

                if(m_opt.verbosity > 1)
                {
                    if(m_iparm[14] > peak_mem1 || m_iparm[15] > peak_mem2 ||
                       m_iparm[16] > peak_mem3)
                    {
                        peak_mem1 = m_iparm[14];
                        peak_mem2 = m_iparm[15];
                        peak_mem3 = m_iparm[16];

                        std::cout << "\nPeak memory on symbolic factorization: "
                                  << "              "
                                  << (double)peak_mem1 / 1024.0 << " Mb";
                        std::cout
                            << "\nPermanent memory on symbolic factorization: "
                            << "         "
                            << (double)peak_mem2 / 1024.0 << " Mb";
                        std::cout << "\nPeak memory on numerical factorization "
                                     "and solution: "
                                  << (double)peak_mem3 / 1024.0 << " Mb"
                                  << std::endl;
                    }
                }

                if(m_error != 0)
                {
                    std::cout << "\nERROR during symbolic factorization...\n";
                    m_check_pardiso_error(m_error);
                    return 11;
                }

                // PHASE 2.
                // Numerical factorization
                m_phase = 22;
                PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size,
                        m_mkl_a, m_ia, m_ja, &m_idum, &m_nrhs, m_iparm,
                        &m_msglvl, &m_ddum, &m_ddum, &m_error);

                if(m_error != 0)
                {
                    std::cout << "\nERROR during numerical factorization...\n";
                    m_check_pardiso_error(m_error);
                    return 22;
                }
            }
            else
            {
                // Time Integrator with the previous Jacobian
                ti(J, b, x, x_prev, m_iterator_state.t, m_iterator_state.dt, false);
            }

            // PHASE 3.
            // Back substitution and iterative refinement
            m_phase = 33;
            PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size,
                    m_mkl_a, m_ia, m_ja, &m_idum, &m_nrhs, m_iparm, &m_msglvl,
                    m_mkl_b, m_mkl_x, &m_error);

            if(m_error != 0)
            {
                std::cout << "\nERROR during solution...\n";
                m_check_pardiso_error(m_error);
                return 33;
            }

            m_calls++;

            double tol = 0.0;

            for(MKL_INT i = 0; i < m_size; i++)
            {
                double adiff = std::abs(m_mkl_x[i]);

                if(adiff > m_opt.value_max || std::isnan(m_mkl_x[i]))
                {
                    std::cout << "\nERROR: Newton iterations diverged. "
                              << "Review the tolerances and/or adaptive time "
                                 "stepping.\n";
                    return 2;
                }

                if(adiff > tol)
                {
                    tol = adiff;
                }

                x[i] -= m_mkl_x[i];
            }

            if(m_opt.verbosity > 0)
            {
                std::cout << "#";
                std::cout.flush();
            }

            if(tol < m_opt.atol)
            {
                break;
            }

        }  // for iter

        // Newton iterator failed to converge within max_Newton_iter iterations
        if(iter == m_opt.max_Newton_iter)
        {
            if(m_opt.verbosity > 0)
                std::cout << " <- redo";
            if(m_reset_ti_state(x, x_prev))
                return 3;  // Newton method failed to converge
            continue;
        }

        // The solver has reached the target time t1 or the stop condition
        // triggered.
        if(m_iterator_state.final_time_step)
        {
            break;
        }
        else if(m_rhs.stop_condition(x, m_iterator_state.t))
        {
            m_iterator_state.dt[1] = m_iterator_state.dt[0];
            break;
        }

        // Adaptive time stepping algorithm
        int status = adaptive_time_stepping(x, x_prev, iter);
        if(status < 0)
            return 4;  // The algorithm failed to converge
        else if(status > 0)
            continue;  // Re-run the current time step

        // Looks like the solver has reached the target time t1
        if(m_iterator_state.t + m_iterator_state.dt[0] >= t1)
        {
            m_iterator_state.final_time_step = true;
            // Adjust the last time step size
            m_iterator_state.dt[1] = m_iterator_state.dt[0];
            m_iterator_state.dt[0] = t1 - m_iterator_state.t;
        }

        // Rewrite solution history
        for(int d = m_opt.bdf_order - 1; d > 0; d--)
        {
            x_prev[d] = x_prev[d - 1];
        }
        x_prev[0] = x;

        // Call Observer to provide a user with intermediate results
        observer(x, m_iterator_state.t);

        m_iterator_state.t += m_iterator_state.dt[0];  // Time step lapse

    }  // while t

    m_opt.t0      = m_iterator_state.t;
    m_opt.dt_init = m_iterator_state.dt[1];

    if(m_opt.verbosity > 0)
        std::cout << "\nLinear algebra solver calls: " << m_calls << '\n';

    // Success
    return 0;
}

/*
 * Releases memory
 */
Solver::~Solver()
{
    m_phase = -1;  // Termination and release of memory
    PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size, &m_ddum,
            m_ia, m_ja, &m_idum, &m_nrhs, m_iparm, &m_msglvl, &m_ddum, &m_ddum,
            &m_error);
}

/*
 * Checks PARDISO solver error messages
 */
void Solver::m_check_pardiso_error(MKL_INT err)
{
    if(!err)
    {
        return;  // All good
    }
    else
    {
        std::cout << "Linear algebra solver ERROR: ";

        switch(err)
        {
        case -1:
            std::cout << "input inconsistent.\n";
            break;
        case -2:
            std::cout << "not enough memory.\n";
            break;
        case -3:
            std::cout << "reordering problem.\n";
            break;
        case -4:
            std::cout << "zero pivot, numerical factorization or iterative "
                         "refinement problem.\n";
            break;
        case -5:
            std::cout << "unclassified (internal) error.\n";
            break;
        case -6:
            std::cout << "reordering failed(matrix types 11 and 13 only).\n";
            break;
        case -7:
            std::cout << "diagonal matrix is singular.\n";
            break;
        case -8:
            std::cout << "32-bit integer overflow problem.\n";
            break;
        case -9:
            std::cout << "not enough memory for OOC.\n";
            break;
        case -10:
            std::cout << "error opening OOC files.\n";
            break;
        case -11:
            std::cout << "read/write error with OOC files.\n";
            break;
        case -12:
            std::cout << "pardiso_64 called from 32-bit library.\n";
            break;
        case -13:
            std::cout
                << "interrupted by the (user-defined) mkl_progress function.\n";
            break;
        case -15:
            std::cout << "internal error which can appear for iparm[23] "
                         "(parallel_fact_control) = 10 and iparm[12] = 1. Try "
                         "switch matching off (set iparm[12] = 0 and rerun).\n";
            break;
        default:
            std::cout << "Unknown error.\n";
            break;
        }

        std::cout << std::endl;
    }
}

}  // namespace daecpp_namespace_name
