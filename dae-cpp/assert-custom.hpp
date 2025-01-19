/*
 * Custom header-only run-time testing and printing.
 *
 * This file is part of dae-cpp.
 *
 * dae-cpp is licensed under the MIT license.
 * A copy of the license can be found in the LICENSE file.
 *
 * Copyright (c) 2024-2025 Ivan Korotkin
 */

#ifndef DAECPP_ASSERT_CUSTOM_H
#define DAECPP_ASSERT_CUSTOM_H

#include <iostream>

namespace daecpp_namespace_name
{

/*
 * ERROR(error message)
 *
 * Prints error message and aborts with error code -1.
 *
 * Example:
 * ERROR("v = " << v);
 */
#define ERROR(msg)                                                                \
    {                                                                             \
        std::cerr << "\nERROR: " << msg << "\nThis error is fatal." << std::endl; \
        exit(-1);                                                                 \
    }

/*
 * WARNING(warning message)
 *
 * Prints warning message.
 *
 * Example:
 * WARNING("v = " << v);
 */
#ifdef TESTING
#define WARNING(msg)
#else
#define WARNING(msg) \
    std::cerr << "WARNING: " << msg << std::endl;
#endif

/*
 * NOTE(note message)
 *
 * Prints a note.
 *
 * Example:
 * NOTE("v = " << v);
 */
#ifdef TESTING
#define NOTE(msg)
#else
#define NOTE(msg) \
    std::cerr << "NOTE: " << msg << std::endl;
#endif

/*
 * PRINT(condition, message)
 *
 * Prints a message if condition (e.g., verbosity level) is true and flushes the output.
 *
 * Example:
 * NOTE(verbosity > 1, "v = " << v);
 */
#ifdef TESTING
#define PRINT(condition, msg)
#else
#define PRINT(condition, msg)          \
    if (condition)                     \
    {                                  \
        std::cout << msg << std::endl; \
    }
#endif

/*
 * ASSERT(expression, error message)
 *
 * Prints error message and aborts with error code -1 if expression is false.
 *
 * Example:
 * ASSERT(v > 0, "Variable v is negative or zero: v = " << v);
 */
#define ASSERT(expression, msg) \
    if (!(expression))          \
    {                           \
        ERROR(msg);             \
    }

/*
 * CHECK(expression, warning message)
 *
 * Prints a warning message if expression is false.
 *
 * Example:
 * CHECK(v > 0, "Variable v is negative or zero: v = " << v);
 */
#define CHECK(expression, msg) \
    if (!(expression))         \
    {                          \
        WARNING(msg);          \
    }

#ifdef DEBUG

/*
 * Same as ERROR(message) but works only if DEBUG is defined.
 */
#define DEBUG_ERROR(msg) ERROR(msg)

/*
 * Same as WARNING(message) but works only if DEBUG is defined.
 */
#define DEBUG_WARNING(msg) WARNING(msg)

/*
 * Same as NOTE(message) but works only if DEBUG is defined.
 */
#define DEBUG_NOTE(msg) NOTE(msg)

/*
 * Same as ASSERT(expression, error message) but works only if DEBUG is defined.
 */
#define DEBUG_ASSERT(expression, msg) ASSERT(expression, msg)

/*
 * Same as CHECK(expression, warning message) but works only if DEBUG is defined.
 */
#define DEBUG_CHECK(expression, msg) CHECK(expression, msg)

#else // #ifdef DEBUG

/*
 * DEBUG is not defined. This macro does nothing.
 */
#define DEBUG_ERROR(msg)

/*
 * DEBUG is not defined. This macro does nothing.
 */
#define DEBUG_WARNING(msg)

/*
 * DEBUG is not defined. This macro does nothing.
 */
#define DEBUG_NOTE(msg)

/*
 * DEBUG is not defined. This macro does nothing.
 */
#define DEBUG_ASSERT(expression, msg)

/*
 * DEBUG is not defined. This macro does nothing.
 */
#define DEBUG_CHECK(expression, msg)

#endif // #ifdef DEBUG

} // namespace daecpp_namespace_name

#endif // DAECPP_ASSERT_CUSTOM_H
