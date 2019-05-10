/*
 * The main solver
 */

#include <iostream>
#include <cmath>

#include <mkl_types.h>
#include <mkl_spblas.h>
#include <mkl_pardiso.h>

#include "solver.h"

namespace daecpp_namespace_name
{

void Solver::operator()(state_type &x)
{
    // Matrix size
    MKL_INT size = (MKL_INT)(x.size());

    // Check user-defined solver options
    m_opt.check_options();

    // Initialise time integrator
    TimeIntegrator ti(m_rhs, m_jac, m_mass, m_opt, size);

    // Initial time and time step
    double t  = 0.0;
    double dt = m_opt.dt_init;

    // Initial output
    if(m_opt.verbosity > 0)
    {
        std::cout << "Number of equations: " << size << std::endl;
        std::cout << "Float precision:     " << 8 * sizeof(float_type)
                  << " bit\n";
        std::cout << "Integer precision:   " << 8 * sizeof(MKL_INT) << " bit\n";
        std::cout << "Numerical algorithm: BDF-" << m_opt.bdf_order
                  << std::endl;
    }

    // Solver starts the first time step using BDF-1 method
    // (since it doesn't have enough history yet)
    int current_scheme = 1;

    // Contains a few latest successful time steps for Time Integrator
    state_type_matrix x_prev(m_opt.bdf_order, state_type(size));

    // Full Jacobian matrix holder
    sparse_matrix_holder J;

    // Temporary Jacobian matrix holder
    sparse_matrix_holder J_tmp;

    // Reserve memory for at least 3-diagonal temporary Jacobian
    J_tmp.A.reserve(3 * size);
    J_tmp.ia.reserve(size + 1);
    J_tmp.ja.reserve(3 * size);

    // Full RHS vector
    state_type b(size);

    // Solution vector used for Newton iterations
    state_type xk(size);

    MKL_INT *ia = nullptr;
    MKL_INT *ja = nullptr;

    float_type *mkl_a;
    float_type *mkl_b = b.data();
    float_type *mkl_x = xk.data();

    // Copy current state vector into the history vector
    x_prev[0] = x;

    int calls = 0;  // Counts linear algebra solver calls

    // PARDISO control parameters
    MKL_INT phase;        // Current phase of the solver
    MKL_INT maxfct = 1;   // Maximum number of numerical factorizations
    MKL_INT mnum   = 1;   // Which factorization to use
    MKL_INT mtype  = 11;  // Real unsymmetric matrix
    MKL_INT nrhs   = 1;   // Number of right hand sides
    MKL_INT msglvl = 0;   // Print statistical information
    MKL_INT error  = 0;   // Initialize error flag

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void *pt[64];

    // Initialise the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver.
    for(MKL_INT i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    // Auxiliary variables
    double  ddum;  // Double dummy
    MKL_INT idum;  // Integer dummy

    // Intel MKL PARDISO iparm parameter
    MKL_INT iparm[64];

    // Load iparm from solver_options class
    m_opt.set_iparm_for_pardiso(iparm);

    /*
     * Start the solver
     * =========================================================================
     */

    // TODO: Start timer here

    int peak_mem1 = 0, peak_mem2 = 0, peak_mem3 = 0;

    bool final_time_step = false;
    int  step_counter    = 0;

    while(t < (m_t1 + dt * 0.5))
    {
        t += dt;  // Time step lapse

        step_counter++;

        if(m_opt.verbosity > 0)
        {
            std::cout << "\nStep " << step_counter << ": \tt = " << t
                      << "   \t:: ";
            std::cout.flush();
        }

        int iter;

        ti.set_scheme(current_scheme);

        if(current_scheme < m_opt.bdf_order)
        {
            current_scheme++;
        }

        for(iter = 0; iter < m_opt.max_Newton_iter; iter++)
        {
            // Time Integrator
            ti(J, b, J_tmp, x, x_prev, t, dt);

            // Jacobian can change its size and can be re-allocated.
            // Catch up new array addresses.
            mkl_a = J.A.data();
            ia    = J.ia.data();
            ja    = J.ja.data();

            if(m_opt.fact_every_iter || iter == 0)
            {
                // PHASE 1.
                // Reordering and Symbolic Factorization. This step also
                // allocates all memory that is necessary for the factorization
                phase = 11;
                PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia,
                        ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
                if(error != 0)
                {
                    std::cout << "\nERROR during symbolic factorization...\n";
                    check_pardiso_error(error);
                    exit(1);
                }

                if(m_opt.verbosity > 1)
                {
                    if(iparm[14] > peak_mem1 || iparm[15] > peak_mem2 ||
                       iparm[16] > peak_mem3)
                    {
                        peak_mem1 = iparm[14];
                        peak_mem2 = iparm[15];
                        peak_mem3 = iparm[16];

                        std::cout << "\nPeak memory on symbolic factorization: "
                                  << (double)peak_mem1 / 1024.0 << " Mb";
                        std::cout
                            << "\nPermanent memory on symbolic factorization: "
                            << (double)peak_mem2 / 1024.0 << " Mb";
                        std::cout << "\nPeak memory on numerical factorization "
                                     "and solution: "
                                  << (double)peak_mem3 / 1024.0 << " Mb"
                                  << std::endl;
                    }
                }

                // PHASE 2.
                // Numerical factorization
                phase = 22;
                PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia,
                        ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
                if(error != 0)
                {
                    std::cout << "\nERROR during numerical factorization...\n";
                    check_pardiso_error(error);
                    exit(2);
                }
            }

            // PHASE 3.
            // Back substitution and iterative refinement
            phase = 33;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia, ja,
                    &idum, &nrhs, iparm, &msglvl, mkl_b, mkl_x, &error);
            if(error != 0)
            {
                std::cout << "\nERROR during solution...\n";
                check_pardiso_error(error);
                exit(3);
            }

            calls++;

            double tol = 0.0;

            for(MKL_INT i = 0; i < size; i++)
            {
                double adiff = std::fabs(mkl_x[i]);
                if(adiff > tol)
                    tol = adiff;
                x[i] -= mkl_x[i];
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

            // Decrease the time step, scrape the current time iteration and
            // carry out it again.
            t -= dt;
            step_counter--;
            final_time_step = false;
            dt /= m_opt.dt_decrease_factor;
            x = x_prev[0];
            continue;
        }

        if(final_time_step)
        {
            break;
        }

        // Simple yet efficient adaptive time stepping
        if(iter < m_opt.dt_increase_threshold)
        {
            dt *= m_opt.dt_increase_factor;
            if(m_opt.verbosity > 0)
                std::cout << '>';
        }
        else if(iter >= m_opt.dt_decrease_threshold - 1)
        {
            dt /= m_opt.dt_decrease_factor;
            if(m_opt.verbosity > 0)
                std::cout << '<';
        }

        if(t + dt > m_t1)
        {
            final_time_step = true;
            // Adjust the last time step size
            dt = m_t1 - t;
        }

        // Rewrite solution history
        for(int d = m_opt.bdf_order - 1; d > 0; d--)
        {
            x_prev[d] = x_prev[d - 1];
        }
        x_prev[0] = x;

    }  // while t

    if(m_opt.verbosity > 0)
        std::cout << "\nLinear algebra solver calls: " << calls << '\n';

    // Termination and release of memory
    phase = -1;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, &ddum, ia, ja, &idum,
            &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
}

void Solver::check_pardiso_error(MKL_INT err)
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
