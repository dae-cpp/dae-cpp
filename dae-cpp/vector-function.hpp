/*
 * Vector function class.
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

#ifndef DAECPP_VECTOR_FUNCTION_H
#define DAECPP_VECTOR_FUNCTION_H

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

/*
 * Vector function class.
 * This class is abstract and must be inherited.
 */
class VectorFunction
{
public:
    /*
     * Defines the RHS (vector function) `f` of the DAE system `M dx/dt = f`.
     * Takes vector x and time t and returns the RHS vector f.
     * Vector f is already pre-allocated with f.size() == x.size().
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(state_type &f, const state_type &x, const double t) const = 0;

    virtual ~VectorFunction() {}
};

} // namespace daecpp_namespace_name

#endif // DAECPP_VECTOR_FUNCTION_H
