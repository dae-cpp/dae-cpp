/*
 * The RHS class.
 * This class is abstract and must be inherited.
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

class RHS
{
public:
    /*
     * Takes vector x and time t and returns vector f.
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(const state_type &x, state_type &f,
                            const double t) = 0;
};

}  // namespace daecpp_namespace_name
