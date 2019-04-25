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
    // very sparse a warning may be given (in verbose mode). In this case
    // max_size_mult should be increased.
    int max_size_mult = 10;

    double atol    = 1.0e-6;  // Absolute tolerance
    double dt_init = 0.1;     // Initial time step

    // MKL PARDISO parameters (iparam). More about iparam:
    // https://software.intel.com/en-us/mkl-developer-reference-c-pardiso-iparm-parameter

    // iparam[3]: Controls preconditioned CGS.
    // 0 by default, 61 - LU-preconditioned CGS iteration with a stopping
    // criterion of 1.0E-6 for nonsymmetric matrices.
    MKL_INT preconditioned_CGS = 0;

    // iparam[7]: Max number of iterative refinement steps.
    // 0 - The solver automatically performs two steps of iterative refinement
    // when perturbed pivots are obtained during the numerical factorization.
    MKL_INT refinement_steps = 2;

    // Solver constants
    const MKL_INT matrix_type = 11;  // Real unsymmetric matrix

    SolverOptions() = default;
};

}  // namespace daecpp_namespace_name
