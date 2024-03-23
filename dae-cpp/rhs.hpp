/*
 * The RHS class.
 * Defines the RHS `f` of `M x = f`.
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
 * Defines the RHS `f` of `M x = f`.
 * This class is abstract and must be inherited.
 */
class RHS
{
public:
    /*
     * Takes vector x and time t and returns the RHS vector f.
     * Vector f is already pre-allocated with f.size() == x.size() and filled with 0.
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(const state_type &x, state_type &f, const double t) const = 0;
};

namespace core
{

/*
 * Defines empty RHS for internal purposes
 */
class RHS_empty : public RHS
{
public:
    void operator()(const state_type &x, state_type &f, const double t) const
    {
        // Generally, empty RHS should never be called
        ERROR("Empty RHS function has been called.");
    }
};

} // namespace core
} // namespace daecpp_namespace_name

#endif // DAECPP_RHS_H
