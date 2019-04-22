/*
 * TODO: Description of the class
 * This class is abstract
 *
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

class RHS
{
public:
    /*
     * TODO: Description
     */
    virtual void operator()(const state_type &x, state_type &f,
                            const double t) = 0;
};

}  // namespace daecpp_namespace_name
