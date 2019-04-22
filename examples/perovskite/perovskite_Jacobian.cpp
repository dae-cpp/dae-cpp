/*
 * TODO: Description of the class
 */

#pragma once

#include "perovskite_Jacobian.h"

void MyJacobian::operator() (state_type &x, state_type &J, vector_type_int &ia, vector_type_int &ja, const double t, const double dt)
{
    const int size = (int)(x.size());
    const double m_lambda = 1.0;

    const int N = size/2; // number of points
    const double L = 1.0; // box length
    const double h = L / (double)N; // cell size
    const double invh = 1.0 / h;
    const double invh2 = invh * invh;
    const double invlam2 = 1.0 / (m_lambda*m_lambda);
    const double dtinvh2 = dt * invh2;

    int c = 0, c1 = 0;
    for (int i = 0; i < size; i++)
    {
        ia[c1++] = c + 1;
        
        if(i == 0)
        {
            ja[c] = 1;
            J[c++] = (-1.0 + 0.5*(x[N+1] - x[N]))*dtinvh2 - 1.0;

            ja[c] = 2;
            J[c++] = (1.0 + 0.5*(x[N+1] - x[N]))*dtinvh2;

            ja[c] = N + 1;
            J[c++] = -0.5*(x[0] + x[1])*dtinvh2;

            ja[c] = N + 2;
            J[c++] = 0.5*(x[0] + x[1])*dtinvh2;
        } 
        else if(i < N-1)
        {
            ja[c] = i;
            J[c++] = (1.0 - 0.5*(x[N+i] - x[N+i-1]))*dtinvh2; // dx[i-1]

            ja[c] = i + 1;
            J[c++] = (-2.0 + 0.5*(x[N+i+1] - 2.0*x[N+i] + x[N+i-1]))*dtinvh2 - 1.0; // dx[i]

            ja[c] = i + 2;
            J[c++] = (1.0 + 0.5*(x[N+i+1] - x[N+i]))*dtinvh2; // dx[i+1]

            ja[c] = N + i;
            J[c++] = 0.5*(x[i] + x[i-1])*dtinvh2; // dx[N+i-1]

            ja[c] = N + i + 1;
            J[c++] = -0.5*(x[i+1] + 2.0*x[i] + x[i-1])*dtinvh2; // dx[N+i]

            ja[c] = N + i + 2;
            J[c++] = 0.5*(x[i+1] + x[i])*dtinvh2; // dx[N+i+1]
        }
        else if(i == N-1)
        {
            ja[c] = i;
            J[c++] = (1.0 - 0.5*(x[2*N-1] - x[2*N-2]))*dtinvh2;

            ja[c] = i + 1;
            J[c++] = (-1.0 - 0.5*(x[2*N-1] - x[2*N-2]))*dtinvh2 - 1.0;

            ja[c] = N + i;
            J[c++] = 0.5*(x[N-1] + x[N-2])*dtinvh2;

            ja[c] = N + i + 1;
            J[c++] = -0.5*(x[N-1] + x[N-2])*dtinvh2;
        }
        else if(i == N)
        {
            ja[c] = 1;
            J[c++] = invlam2;

            ja[c] = N + 1;
            J[c++] = -3.0*invh2;

            ja[c] = N + 2;
            J[c++] = invh2;
        }
        else if(i < 2*N-1)
        {
            ja[c] = i - N + 1;
            J[c++] = invlam2; // dx[i]

            ja[c] = i;
            J[c++] = invh2; // dx[N+i-1]

            ja[c] = i + 1;
            J[c++] = -2.0*invh2; // dx[N+i]

            ja[c] = i + 2;
            J[c++] = invh2; // dx[N+i+1]
        }
        else // i == 2*N-1
        {
            ja[c] = N;
            J[c++] = invlam2;

            ja[c] = 2*N - 1;
            J[c++] = invh2;

            ja[c] = 2*N;
            J[c++] = -3.0*invh2;
        }
    }

    ia[c1] = c + 1;
/*
    for (int i = 0; i < 5 * size; i++)
    {
        if (i <= size)
        {
            std::cout << i << ": " << J[i] << '\t' << ja[i] << '\t' << ia[i] << '\n';
        }
        else
        {
            std::cout << i << ": " << J[i] << '\t' << ja[i] << '\n';
        }
    }
    exit(0);
    */
}
