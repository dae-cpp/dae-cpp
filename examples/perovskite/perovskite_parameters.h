/*
 * TODO: Description of the class
 */

#pragma once

class MyParams
{
public:
    // For simplicity in this example all parameters are public.
    // Generally, they should be private with getters/setters.
    int    N;
    double L;
    double lambda;
    double t1;

    MyParams(int N, double L, double lambda, double t1)
        : N(N), L(L), lambda(lambda), t1(t1)
    {
    }
}
