/*
 * TODO: Analytical Jacobian description
 */

#include "perovskite_Jacobian.h"

void MyJacobian::operator()(daecpp::sparse_matrix_holder &J,
                            daecpp::state_type &x, const double t)
{
    // Locals
    const int    N       = m_p.N;
    const int    size    = (int)(x.size());
    const double invh2   = m_p.invh * m_p.invh;
    const double invlam2 = 1.0 / (m_p.lambda * m_p.lambda);

    int c = 0;

    for(int i = 0; i < size; i++)
    {
        J.ia.push_back(c + 1);

        if(i == 0)
        {
            J.ja.push_back(1);
            J.A.push_back((-1.0 + 0.5 * (x[N + 1] - x[N])) * invh2);

            J.ja.push_back(2);
            J.A.push_back((1.0 + 0.5 * (x[N + 1] - x[N])) * invh2);

            J.ja.push_back(N + 1);
            J.A.push_back(-0.5 * (x[0] + x[1]) * invh2);

            J.ja.push_back(N + 2);
            J.A.push_back(0.5 * (x[0] + x[1]) * invh2);

            c += 4;
        }
        else if(i < N - 1)
        {
            J.ja.push_back(i);
            J.A.push_back((1.0 - 0.5 * (x[N + i] - x[N + i - 1])) *
                          invh2);  // dx[i-1]

            J.ja.push_back(i + 1);
            J.A.push_back(
                (-2.0 + 0.5 * (x[N + i + 1] - 2.0 * x[N + i] + x[N + i - 1])) *
                invh2);  // dx[i]

            J.ja.push_back(i + 2);
            J.A.push_back((1.0 + 0.5 * (x[N + i + 1] - x[N + i])) *
                          invh2);  // dx[i+1]

            J.ja.push_back(N + i);
            J.A.push_back(0.5 * (x[i] + x[i - 1]) * invh2);  // dx[N+i-1]

            J.ja.push_back(N + i + 1);
            J.A.push_back(-0.5 * (x[i + 1] + 2.0 * x[i] + x[i - 1]) *
                          invh2);  // dx[N+i]

            J.ja.push_back(N + i + 2);
            J.A.push_back(0.5 * (x[i + 1] + x[i]) * invh2);  // dx[N+i+1]

            c += 6;
        }
        else if(i == N - 1)
        {
            J.ja.push_back(i);
            J.A.push_back((1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2);

            J.ja.push_back(i + 1);
            J.A.push_back((-1.0 - 0.5 * (x[2 * N - 1] - x[2 * N - 2])) * invh2);

            J.ja.push_back(N + i);
            J.A.push_back(0.5 * (x[N - 1] + x[N - 2]) * invh2);

            J.ja.push_back(N + i + 1);
            J.A.push_back(-0.5 * (x[N - 1] + x[N - 2]) * invh2);

            c += 4;
        }
        else if(i == N)
        {
            J.ja.push_back(1);
            J.A.push_back(invlam2);

            J.ja.push_back(N + 1);
            J.A.push_back(-3.0 * invh2);

            J.ja.push_back(N + 2);
            J.A.push_back(invh2);

            c += 3;
        }
        else if(i < 2 * N - 1)
        {
            J.ja.push_back(i - N + 1);
            J.A.push_back(invlam2);  // dx[i]

            J.ja.push_back(i);
            J.A.push_back(invh2);  // dx[N+i-1]

            J.ja.push_back(i + 1);
            J.A.push_back(-2.0 * invh2);  // dx[N+i]

            J.ja.push_back(i + 2);
            J.A.push_back(invh2);  // dx[N+i+1]

            c += 4;
        }
        else  // i == 2*N-1
        {
            J.ja.push_back(N);
            J.A.push_back(invlam2);

            J.ja.push_back(2 * N - 1);
            J.A.push_back(invh2);

            J.ja.push_back(2 * N);
            J.A.push_back(-3.0 * invh2);

            c += 3;
        }
    }

    J.ia[size] = c + 1;
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
