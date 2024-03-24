/*
 * Defines DAE solver options.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_SOLVER_OPTIONS_H
#define DAECPP_SOLVER_OPTIONS_H

#include "typedefs.hpp"

namespace daecpp_namespace_name
{
namespace verbosity
{

enum VerbosityLevel
{
    silent,
    medium,
    loud
};

} // namespace verbosity

/*
 * Defines DAE solver options.
 */
struct SolverOptions
{
    verbosity::VerbosityLevel verbosity{verbosity::silent};

    double dt_init{0.001}; // Initial time step
};

namespace core
{

constexpr int MAX_ORDER{2};

} // namespace core
} // namespace daecpp_namespace_name

#endif // DAECPP_SOLVER_OPTIONS_H
