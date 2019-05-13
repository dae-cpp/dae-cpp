/*
 * Contains public solver parameters
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

#define BDF_MAX_ORDER 6  // Max BDF scheme order currently implemented

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

    // Perform symbolic and numerical factorisation every Newton iteration
    // (changing to 'false' can increase speed but also can lead to instability)
    bool fact_every_iter = true;

    // Order of BDF implicit numerical integration method:
    // 1 - first order BDF, 2 - BDF-2, ..., 6 - BDF-6
    int bdf_order = 6;

    // Time stepping algorithm:
    // 1 - Stability-based Simple Adaptive Time Stepping (S-SATS),
    // 2 - Accuracy-based Simple Adaptive Time Stepping (A-SATS) from
    // https://www.sciencedirect.com/science/article/pii/S0377042705005534
    int time_stepping = 1;

    // Maximum number of Newton iterations. If the Newton method fails to
    // converge after max_Newton_iter iterations, the solver reduces time step
    // and tries to make the current step again.
    int max_Newton_iter = 15;

    // Absolute tolerance for the Newton algorithm
#ifdef DAE_SINGLE
    double atol = 1.0e-3;  // Absolute tolerance for single precision
#else
    double atol = 1.0e-6;  // Absolute tolerance for double precision
#endif

#ifdef DAE_SINGLE
    double dt_eps_m =
        1.0e-5;  // The order of the rounding unit for single precision
#else
    double dt_eps_m =
        1.0e-10;  // The order of the rounding unit for double precision
#endif

    // Initial time step
    double dt_init = 0.1;

    // Minimum and maximum time steps
    double dt_min = dt_eps_m;
    double dt_max = 100.0;

    // Verbosity level of the solver:
    // 0 - be silent, 1 - prints some basic information, 2 - chatterbox
    int verbosity = 1;

    // Simple Adaptive Time Stepping options
    int dt_increase_threshold = 3;    // Time step amplification threshold
                                      // (S-SATS only)
    int dt_decrease_threshold = 7;    // Time stepreduction threshold
                                      // (S-SATS only)
    double dt_increase_factor = 1.4;  // Time step amplification factor
    double dt_decrease_factor = 1.4;  // Time step reduction factor
    double dt_eta_min = 0.01;  // Monitor function lower threshold (A-SATS only)
    double dt_eta_max = 0.1;  // Monitor function higher threshold (A-SATS only)

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
    MKL_INT parallel_fact_control = 0;

    // iparm[24]: Parallel forward/backward solve control.
    // 0 - Intel MKL PARDISO uses the sequential forward and backward solve.
    // 1 - Parallel algorithm for the solve step (in-core mode only).
    // 2 - If there is only one right hand side, Intel MKL PARDISO performs
    //     parallelization. If there are multiple right hand sides, PARDISO
    //     performs parallel forward and backward substitution via matrix.
    MKL_INT parallel_solve_control = 0;

    SolverOptions() = default;

    // Initialises Intel MKL PARDISO parameters (iparam) array
    void set_iparm_for_pardiso(MKL_INT *iparm);

    // Checks correctness of the solver parameters
    void check_options();
};

}  // namespace daecpp_namespace_name
