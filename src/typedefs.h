/*
 * TODO: Description
 */

#pragma once

#include <mkl_types.h>
#include <vector>

#define daecpp_namespace_name daecpp

namespace daecpp_namespace_name
{

typedef double                  float_type;
typedef std::vector<float_type> state_type;
typedef std::vector<MKL_INT>    vector_type_int;

}  // namespace daecpp_namespace_name
