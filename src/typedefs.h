/*
 * TODO: Description
 */

#pragma once

#include <mkl_types.h>
#include <vector>

#include "cmake_config.h"

#define daecpp_namespace_name daecpp

namespace daecpp_namespace_name
{

#ifdef DAE_SINGLE
typedef float float_type;
#else
typedef double float_type;
#endif

typedef std::vector<float_type> state_type;
typedef std::vector<MKL_INT>    vector_type_int;

typedef std::vector<std::vector<float_type>> state_type_matrix;

struct sparse_matrix_holder
{
    state_type      A;   // Non-zero elements of the sparse matrix A
    vector_type_int ia;  // Points to the first column index of the given row
                         // in the array ja
    vector_type_int ja;  // Contains column indices of the sparse matrix A
};

}  // namespace daecpp_namespace_name
