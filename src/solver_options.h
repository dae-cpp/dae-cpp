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
    int max_size_mult = 10;  // Reserve memory for N*max_size_mult sparse matrix elements. Usually 10*N elements is more than enough, but if Mass or Jacobian matrices are not very sparse a warning may be given (in verbose mode). In this case max_size_mult should be increased.

    double atol = 1.0e-6;  // Absolute tolerance
    double dt_init = 0.1;  // Initial time stem

    SolverOptions() = default;
};

}  // namespace daecpp_namespace_name
