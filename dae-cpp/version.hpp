/*
 * Defines major, minor, and patch version of the DAE solver.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#ifndef DAECPP_VERSION_H
#define DAECPP_VERSION_H

#include "typedefs.hpp"

namespace daecpp_namespace_name
{

// dae-cpp library major version
static constexpr uint16_t version_major{DAECPP_VERSION_MAJOR};

// dae-cpp library minor version
static constexpr uint16_t version_minor{DAECPP_VERSION_MINOR};

// dae-cpp library patch version
static constexpr uint16_t version_patch{DAECPP_VERSION_PATCH};

} // namespace daecpp_namespace_name

#endif // DAECPP_VERSION_H
