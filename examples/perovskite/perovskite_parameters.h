/*
 * Parameter container for the perovskite problem.
 * For simplicity all parameters in this example are public. Generally, they
 * should be private with getters/setters.
 */

#pragma once

class MyParams
{
public:
    MKL_INT N;       // Number of cells
    double  L;       // Space interval length
    double  lambda;  // Lambda parameter
    double  t1;      // Integration time (0 < t < t1)

    // Derived parameters
    const double h    = L / (double)(N - 1);  // cell size
    const double invh = 1.0 / h;              // inverse cell size

    MyParams(MKL_INT N, double L, double lambda, double t1)
        : N(N), L(L), lambda(lambda), t1(t1)
    {
    }
};
