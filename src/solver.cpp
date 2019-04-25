/*
 * TODO: Description
 */

#include <iostream>
#include <cmath>

#include <mkl_types.h>
#include <mkl_spblas.h>
#include <mkl_pardiso.h>

#include "solver.h"

#define BDF_MAX_ORDER 1  // Max BDF scheme order currently implemented

namespace daecpp_namespace_name
{

void Solver::operator()(state_type &x)
{
    // Matrix size
    MKL_INT size = (MKL_INT)(x.size());

    // Initialise time integrator
    TimeIntegrator ti(m_rhs, m_jac, m_mass, m_opt, size);

    // Initial time and time step
    double t  = 0.0;
    double dt = m_opt.dt_init;

    // Contains a few latest successful time steps for Time Integrator
    state_type x_prev[BDF_MAX_ORDER];

    // Full Jacobian matrix holder
    sparse_matrix_holder J;

    // Full RHS and solution vector
    state_type b(size), xk(size);

    MKL_INT *ia;
    MKL_INT *ja;

    float_type *mkl_a;
    float_type *mkl_b = b.data();
    float_type *mkl_x = xk.data();

    // Copy current state vector into the history vector
    x_prev[0] = x;

    int calls = 0;  // Counts linear algebra solver calls

    /*
     * Pardiso control parameters
     * ==========================
     */

    MKL_INT phase;  // Current phase of the solver
    MKL_INT mtype  = m_opt.matrix_type;
    MKL_INT nrhs   = 1;  // Number of right hand sides
    MKL_INT maxfct = 1;  // Maximum number of numerical factorizations
    MKL_INT mnum   = 1;  // Which factorization to use
    MKL_INT msglvl = 0;  // Print statistical information
    MKL_INT error  = 0;  // Initialize error flag

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void *pt[64];

    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver.
    for(MKL_INT i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    // Auxiliary variables
    double  ddum;  // Double dummy -- shouldn't be 'float_type'?
    MKL_INT idum;  // Integer dummy

    MKL_INT iparm[64];

    for(MKL_INT i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }

    iparm[0] = 1;  // No solver default, must supply all values in components
                   // iparm[1] - iparm[63].
    iparm[1] = 3;  // The parallel (OpenMP) version of the nested dissection
                   // algorithm. It can decrease the time of computations on
                   // multi-core computers, especially when Intel MKL PARDISO
                   // Phase 1 takes significant time.

    iparm[3] = m_opt.preconditioned_CGS;  // Controls preconditioned CGS/CG

    iparm[4] = 0;  // No user fill-in reducing permutation
    iparm[5] = 0;  // Write solution into x
    iparm[6] = 0;  // Not in use

    iparm[7] = m_opt.refinement_steps;  // Number of iterative refinement steps

    iparm[8]  = 0;   // Not in use
    iparm[9]  = 13;  // Perturb the pivot elements with 1.0E-13
    iparm[10] = 1;   // Use nonsymmetric permutation and scaling MPS

    iparm[11] = 0;   // Conjugate transposed/transpose solve
    iparm[12] = 1;   // Maximum weighted matching algorithm is switched-on
    iparm[13] = 0;   // Output: Number of perturbed pivots
    iparm[14] = 0;   // Not in use
    iparm[15] = 0;   // Not in use
    iparm[16] = 0;   // Not in use
    iparm[17] = -1;  // Output: Number of nonzeros in the factor LU
    iparm[18] = -1;  // Output: Mflops for LU factorization
    iparm[19] = 0;   // Output: Numbers of CG Iterations

    // iparm[23] = 10;
    // iparm[24] = 2;

    /*
     * Start solver
     * ============
     */

    // Start timer here

    bool final_time_step = false;

    for(double t = dt; t <= (m_t1 + dt * 0.5); t += dt)
    {
        std::cout << "\nt = " << t << ": ";

        int iter;

        // remove 20 to options
        for(iter = 0; iter < 20; iter++)
        {
            // estimate new numerical J only for inter == 0?
            ti(J, b, x, x_prev, t, dt);

            // Jacobian can change its size and can be re-allocated.
            // Catch up new array addresses.
            mkl_a = J.A.data();
            ia    = J.ia.data();
            ja    = J.ja.data();

            if(iter == 0)
            {
                // PHASE 1.
                // Reordering and Symbolic Factorization. This step also
                // allocates all memory that is necessary for the factorization
                phase = 11;
                PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia,
                        ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
                if(error != 0)
                {
                    printf("\nERROR during symbolic factorization: %d", error);
                    exit(1);
                }
                // printf("\nReordering completed ... ");
                // printf("\nNumber of nonzeros in factors = %d", iparm[17]);
                // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

                // PHASE 2.
                // Numerical factorization
                phase = 22;
                PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia,
                        ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
                if(error != 0)
                {
                    printf("\nERROR during numerical factorization: %d", error);
                    exit(2);
                }
                // printf("\nFactorization completed ... ");
            }

            // PHASE 3.
            // Back substitution and iterative refinement
            phase = 33;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia, ja,
                    &idum, &nrhs, iparm, &msglvl, mkl_b, mkl_x, &error);
            if(error != 0)
            {
                printf("\nERROR during solution: %d", error);
                exit(3);
            }

            calls++;

            double tol = 0.0;
            for(int i = 0; i < size; i++)
            {
                double adiff = std::abs(mkl_x[i]);
                if(adiff > tol)
                    tol = adiff;
                x[i] -= mkl_x[i];
            }

            std::cout << iter << ':' << tol << ' ';

            if(tol < m_opt.atol)
                break;
        }  // for iter

        std::cout << '=' << iter << "= dt: " << dt;

        if(final_time_step)
            break;

        // Temporary adaptive time stepping
        if(iter <= 3)
            dt *= 1.4;
        else if(iter > 6)
            dt /= 1.4;

        if(t + dt > m_t1)
        {
            final_time_step = true;
            // Adjust the last time step size
            dt = m_t1 - t;
        }

        x_prev[0] = x;

    }  // for t

    std::cout << '\n' << calls << " calls\n";

    // Termination and release of memory
    phase = -1;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, &ddum, ia, ja, &idum,
            &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
}

}  // namespace daecpp_namespace_name
