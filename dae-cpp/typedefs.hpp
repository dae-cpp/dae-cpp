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

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

#define daecpp_namespace_name daecpp

// dae-cpp version
#define DAECPP_VERSION_MAJOR 2
#define DAECPP_VERSION_MINOR 1
#define DAECPP_VERSION_PATCH 0

// Internal constants
#define DAECPP_MAX_ORDER 4
#define DAECPP_SPARSE_MATRIX_ELEMENT_TOLERANCE 1e-14 // Used in automatic (algorithmic) Jacobian

#include "assert-custom.hpp"

namespace daecpp_namespace_name
{

// Solver exit codes
namespace exit_code
{
enum status
{
    success = 0,
    diverged = 1,
    linsolver_failed_decomposition = 2,
    linsolver_failed_solving = 3,
    unknown = 10
};
} // namespace exit_code

// Unsigned integer type
#ifdef DAECPP_LONG
typedef uint64_t int_type;
#else
typedef uint32_t int_type;
#endif

// Floating point scalar
#ifdef DAECPP_SINGLE
typedef float float_type;
#else
typedef double float_type;
#endif

// Floating point dual number for automatic differentiation
typedef autodiff::real1st dual_type;

namespace core
{

// Unsigned integer vector
typedef std::vector<int_type> ivec;

// Floating point (double or single precision) vector
typedef std::vector<float_type> rvec;

// Eigen vector type
#ifdef DAECPP_SINGLE
typedef Eigen::VectorXf eivec;
#else
typedef Eigen::VectorXd eivec;
#endif

// Eigen sparse matrix type
typedef Eigen::SparseMatrix<float_type> eimat;

} // namespace core

// State vector
typedef core::rvec state_vector;

// State vector for automatic (algorithmic) differentiation
typedef autodiff::VectorXreal state_type;

} // namespace daecpp_namespace_name

#endif // DAECPP_TYPEDEFS_H
