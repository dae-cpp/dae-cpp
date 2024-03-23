/*
 * The RHS class.
 * Defines the RHS `f` of the DAE system `M dx/dt = f`.
 * This class is abstract and must be inherited.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#ifndef DAECPP_RHS_H
#define DAECPP_RHS_H

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * The RHS class.
 * This class is abstract and must be inherited.
 */
class RHS
{
public:
    /*
     * Defines the RHS `f` of the DAE system `M dx/dt = f`.
     * Takes vector x and time t and returns the RHS vector f.
     * Vector f is already pre-allocated with f.size() == x.size() and filled with 0.
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(state_type &f, const state_type &x, const double t) const = 0;
};

} // namespace daecpp_namespace_name

#endif // DAECPP_RHS_H
