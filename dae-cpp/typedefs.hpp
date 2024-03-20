/*
 * Defines custom types and shortcuts used in the project.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024 Ivan Korotkin
 */

#pragma once

#include <vector>

#define daecpp_namespace_name daecpp

namespace daecpp_namespace_name
{

#ifdef DAECPP_SINGLE
typedef float float_type;
#else
typedef double float_type;
#endif

// Integer vector
typedef std::vector<int> ivec;

// Floating point (double or single precision) vector
typedef std::vector<float_type> fvec;

// Matrix structure in 3-array format
struct sparse_matrix_holder
{
    fvec A; // Non-zero element A
    ivec i; // Row index of the element A
    ivec j; // Column index of the element A

    void check(); // Checks the matrix structure
};

} // namespace daecpp_namespace_name
