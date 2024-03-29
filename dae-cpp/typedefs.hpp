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

#ifndef DAECPP_TYPEDEFS_H
#define DAECPP_TYPEDEFS_H

#include <vector>

#include <Eigen/Sparse>

#define daecpp_namespace_name daecpp

#include "assert-custom.hpp"

namespace daecpp_namespace_name
{

// dae-cpp library major version
constexpr u_int16_t version_major{2};

// dae-cpp library minor version
constexpr u_int16_t version_major{0};

// Integer type
#ifdef DAECPP_LONG
typedef u_int64_t int_type;
#else
typedef u_int32_t int_type;
#endif

// Floating point scalar
#ifdef DAECPP_SINGLE
typedef float float_type;
#else
typedef double float_type;
#endif

// Integer vector
typedef std::vector<int_type> ivec;

// Floating point (double or single precision) vector
typedef std::vector<float_type> fvec;

// State vector
typedef fvec state_type;

namespace core
{

// Eigen vector type
#ifdef DAECPP_SINGLE
typedef Eigen::VectorXf eivec;
#else
typedef Eigen::VectorXd eivec;
#endif

} // namespace core
} // namespace daecpp_namespace_name

#endif // DAECPP_TYPEDEFS_H
