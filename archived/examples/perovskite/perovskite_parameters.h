/*
 * Parameter container for the perovskite problem.
 */

#pragma once

struct MyParams
{
    const MKL_INT N      = 4000;  // Number of points
    const double  L      = 1.0;   // Space interval length
    const double  lambda = 1.0;   // Lambda parameter
    const double  t1     = 10.0;  // Integration time (0 < t < t1)

    // Derived parameters
    const double h    = L / (double)(N - 1);  // cell size
    const double invh = 1.0 / h;              // inverse cell size
};
