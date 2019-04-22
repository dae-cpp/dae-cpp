/*
 * TODO: Description of the class
 * This class is abstract
 *
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

class MassMatrix
{
public:
    /*
     * TODO: Operator description
     * Sparse matrix format
     * May depend on t but does not depend on x at the moment!
     */
    virtual void operator()(state_type &M, vector_type_int &ia,
                            vector_type_int &ja, const double t) = 0;

    /*
     * TODO: Helper function(s) to create mass matrices in sparse format
     */
};

}  // namespace daecpp_namespace_name
