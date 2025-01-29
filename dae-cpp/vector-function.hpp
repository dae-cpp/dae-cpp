/*
 * Vector function classes.
 * Define the RHS `f` of the DAE system `M dx/dt = f`.
 * This class is abstract and must be inherited.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
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
     * Takes vector `x` and time `t` and returns the RHS vector `f`.
     * Vector `f` is already pre-allocated with f.size() == x.size().
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(state_type &f, const state_type &x, const double t) const = 0;

    virtual ~VectorFunction() {}
};

/*
 * Vector function class used for automatic Jacobian computed from the user-defined Jacobian shape.
 * Used to define the vector function (the RHS) element-by-element (equation-by-equation).
 * Should be used together with `JacobianMatrixShape` class for the Jacobian matrix.
 * This class is abstract and must be inherited.
 */
class VectorFunctionElements
{
public:
    /*
     * Defines the RHS (vector function) `f` of the DAE system `M dx/dt = f`.
     * Vector `f` is already pre-allocated with f.size() == x.size().
     */
    void operator()(state_type &f, const state_type &x, const double t) const
    {
        const int_type size = static_cast<int_type>(x.size()); // System size

        for (int_type i = 0; i < size; ++i)
        {
            f[i] = equations(x, t, i);
        }
    }

    /*
     * All RHS functions `f_i` for each equation `i` in the system.
     * Takes vector `x`, time `t`, index `i`, and returns the RHS value `f_i` for the given equation `i`.
     * I.e., it returns the i-th element of the vector function.
     * This function is pure virtual and must be overriden.
     */
    virtual state_value equations(const state_type &x, const double t, const int_type i) const = 0;

    virtual ~VectorFunctionElements() {}
};

} // namespace daecpp_namespace_name

#endif // DAECPP_VECTOR_FUNCTION_H
