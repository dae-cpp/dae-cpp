/*
 * Mass matrix class.
 * This class is abstract and must be inherited.
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

class MassMatrix
{
public:
    /*
     * The matrix should be defined in sparse format,
     * see three array sparse format decription on
     * https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format
     *
     * The Mass matrix is static, i.e. it does not depend on time t and
     * the vector x
     *
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(sparse_matrix_holder &M) = 0;

    /*
     * TODO: Helper function(s) to create mass matrices in sparse format
     */
};

}  // namespace daecpp_namespace_name
