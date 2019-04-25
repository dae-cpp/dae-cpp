/*
 * TODO: Analytical Jacobian description
 */

#pragma once

#include "perovskite_Jacobian.h"

void MyJacobian::operator()(daecpp::state_type &x,
                            daecpp::sparse_matrix_holder &J, const double t,
                            const double dt)
{
    // Locals
    const int    N       = m_p.N;
    const int    size    = (int)(x.size());
    const double invh2   = m_p.invh * m_p.invh;
    const double invlam2 = 1.0 / (m_p.lambda * m_p.lambda);
    const double dtinvh2 = dt * invh2;

    int c = 0, c1 = 0;
    for(int i = 0; i < size; i++)
    {
        J.ia[c1++] = c + 1;

        if(i == 0)
        {
            J.ja[c]  = 1;
            J.A[c++] = (-1.0 + 0.5 * (x[N + 1] - x[N])) * dtinvh2 - 1.0;

            J.ja[c]  = 2;
            J.A[c++] = (1.0 + 0.5 * (x[N + 1] - x[N])) * dtinvh2;

            J.ja[c]  = N + 1;
            J.A[c++] = -0.5 * (x[0] + x[1]) * dtinvh2;

            J.ja[c]  = N + 2;
            J.A[c++] = 0.5 * (x[0] + x[1]) * dtinvh2;
        }
        else if(i < N - 1)
        {
            J.ja[c] = i;
            J.A[c++] =
                (1.0 - 0.5 * (x[N + i] - x[N + i - 1])) * dtinvh2;  // dx[i-1]

            J.ja[c] = i + 1;
            J.A[c++] =
                (-2.0 + 0.5 * (x[N + i + 1] - 2.0 * x[N + i] + x[N + i - 1])) *
                    dtinvh2 -
                1.0;  // dx[i]

            J.ja[c] = i + 2;
            J.A[c++] =
                (1.0 + 0.5 * (x[N + i + 1] - x[N + i])) * dtinvh2;  // dx[i+1]

            J.ja[c]  = N + i;
            J.A[c++] = 0.5 * (x[i] + x[i - 1]) * dtinvh2;  // dx[N+i-1]

            J.ja[c] = N + i + 1;
            J.A[c++] =
                -0.5 * (x[i + 1] + 2.0 * x[i] + x[i - 1]) * dtinvh2;  // dx[N+i]

            J.ja[c]  = N + i + 2;
            J.A[c++] = 0.5 * (x[i + 1] + x[i]) * dtinvh2;  // dx[N+i+1]
        }
        else if(i == N - 1)
        {
            J.ja[c]  = i;
            J.A[c++] = (1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * dtinvh2;

            J.ja[c] = i + 1;
            J.A[c++] =
                (-1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * dtinvh2 - 1.0;

            J.ja[c]  = N + i;
            J.A[c++] = 0.5 * (x[N - 1] + x[N - 2]) * dtinvh2;

            J.ja[c]  = N + i + 1;
            J.A[c++] = -0.5 * (x[N - 1] + x[N - 2]) * dtinvh2;
        }
        else if(i == N)
        {
            J.ja[c]  = 1;
            J.A[c++] = invlam2;

            J.ja[c]  = N + 1;
            J.A[c++] = -3.0 * invh2;

            J.ja[c]  = N + 2;
            J.A[c++] = invh2;
        }
        else if(i < 2 * N - 1)
        {
            J.ja[c]  = i - N + 1;
            J.A[c++] = invlam2;  // dx[i]

            J.ja[c]  = i;
            J.A[c++] = invh2;  // dx[N+i-1]

            J.ja[c]  = i + 1;
            J.A[c++] = -2.0 * invh2;  // dx[N+i]

            J.ja[c]  = i + 2;
            J.A[c++] = invh2;  // dx[N+i+1]
        }
        else  // i == 2*N-1
        {
            J.ja[c]  = N;
            J.A[c++] = invlam2;

            J.ja[c]  = 2 * N - 1;
            J.A[c++] = invh2;

            J.ja[c]  = 2 * N;
            J.A[c++] = -3.0 * invh2;
        }
    }

    J.ia[c1] = c + 1;
    /*
        for(int i = 0; i < 5 * size; i++)
        {
            if(i <= size)
            {
                std::cout << i << ": " << J[i] << '\t' << ja[i] << '\t' << ia[i]
       << '\n';
            }
            else
            {
                std::cout << i << ": " << J[i] << '\t' << ja[i] << '\n';
            }
        }
        exit(0);
    */
}
