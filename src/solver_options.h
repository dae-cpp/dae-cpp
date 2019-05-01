/*
 * TODO: Description of the class
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

class SolverOptions
{
public:
    /*
     * List of public solver options
     */

    // You can control the parallel execution of the solver by explicitly
    // setting the MKL_NUM_THREADS environment variable. If fewer OpenMP threads
    // are available than specified, the execution may slow down instead of
    // speeding up. If MKL_NUM_THREADS is not defined, then the solver uses all
    // available processors.
    
    // Reserve memory for N*max_size_mult sparse matrix elements. Usually 10*N
    // elements is more than enough, but if Mass or Jacobian matrices are not
    // very sparse a warning may be given. In this case
    // max_size_mult should be increased.
    int max_size_mult = 10;  // TODO: Get rid of it!

    // Solver constants

    bool fact_every_iter = true;  // Perform symbolic and numerical factorisation every Newton iteration (changing to 'false' can increase speed but also can lead to instability)

    int max_Newton_iter = 20;  // Maximum number of Newton iterations
    
    double atol    = 1.0e-6;  // Absolute tolerance
    double dt_init = 0.1;     // Initial time step

    // Intel MKL PARDISO parameters (iparam). More about iparam:
    // https://software.intel.com/en-us/mkl-developer-reference-c-pardiso-iparm-parameter

    // iparam[3]: Controls preconditioned CGS.
    // 0  - The factorization is always computed as required.
    // 31 - LU-preconditioned CGS iteration with a stopping criterion of 1.0E-3.
    // 61 - LU-preconditioned CGS iteration with a stopping criterion of 1.0E-6.
    MKL_INT preconditioned_CGS = 0;

    // iparam[7]: Maximum number of iterative refinement steps.
    // 0  - The solver automatically performs two steps of iterative refinement.
    // >0 - Maximum number of iterative refinement steps that the solver
    //      performs.
    MKL_INT refinement_steps = 2;

    // iparm[23]: Parallel factorization control.
    // 0  - Intel MKL PARDISO uses the classic algorithm for factorization.
    // 1  - Two-level factorization algorithm. This algorithm generally improves
    //      scalability in case of parallel factorization on many OpenMP threads
    //      (more than eight).
    // 10 - Improved two-level factorization algorithm.
    MKL_INT parallel_fact_control = 10;

    // iparm[24]: Parallel forward/backward solve control.
    // 0 - Intel MKL PARDISO uses the sequential forward and backward solve.
    // 1 - Parallel algorithm for the solve step (in-core mode only).
    // 2 - If there is only one right hand side, Intel MKL PARDISO performs
    //     parallelization. If there are multiple right hand sides, PARDISO
    //     performs parallel forward and backward substitution via matrix.
    MKL_INT parallel_solve_control = 2;

    SolverOptions() = default;

    /*
     * Intel MKL PARDISO iparm parameters:
     * https://software.intel.com/en-us/mkl-developer-reference-c-pardiso-iparm-parameter
     */
    void set_iparm_for_pardiso(MKL_INT *iparm)
    {
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

        iparm[3] = preconditioned_CGS;  // Controls preconditioned CGS/CG

        iparm[4] = 0;  // No user fill-in reducing permutation
        iparm[5] = 0;  // Write solution into x
        iparm[6] = 0;  // Reports the number of iterative refinement steps performed

        iparm[7] = refinement_steps;  // Number of iterative refinement steps

        iparm[9]  = 13;  // Perturb the pivot elements with 1.0E-13 (default)
        iparm[10] = 1;   // Enable scaling. Default for nonsymmetric matrices.
        iparm[11] = 0;   // Conjugate transposed/transpose solve
        iparm[12] = 1;   // Maximum weighted matching algorithm is switched-on
        iparm[13] = 0;   // Output: Number of perturbed pivots
        iparm[14] = 0;   // Output: Peak memory on symbolic factorization
        iparm[15] = 0;   // Output: Permanent memory on symbolic factorization
        iparm[16] = 0;   // Output: Size of factors/Peak memory on numerical factorization and solution
        iparm[17] = 0;   // Report the number of nonzeros in the factor LU
        iparm[18] = 0;   // Report Mflops for LU factorization
        iparm[19] = 0;   // Output: Numbers of CG Iterations

        iparm[23] = parallel_fact_control;  // Parallel factorization control

        if(iparm[23] == 1)
        {
            iparm[10] = 0;  // Disable scaling
            iparm[12] = 0;  // Disable matching
        }

        iparm[24] = parallel_solve_control;  // Parallel forward/backward solve control
        
        iparm[26] = 0;  // Matrix checker
        iparm[27] = 0;  // Single or double precision Intel MKL PARDISO:
                        // 0 - double, 1 - single
        iparm[29] = 0;  // Number of zero or negative pivots
        iparm[30] = 0;  // Partial solve and computing selected components
        iparm[33] = 0;  // Because iparm[1] = 3;
        iparm[34] = 0;  // One- or zero-based indexing of columns and rows:
                        // 0 - one-based, 1 - zero-based
        iparm[59] = 0;  // Intel MKL PARDISO mode - in-core
    }
};

}  // namespace daecpp_namespace_name
