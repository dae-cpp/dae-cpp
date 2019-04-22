/*
 * TODO: Description
 */

#include <iostream>
#include <cmath>

#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl_spblas.h>

#include "solver.h"

namespace daecpp_namespace_name
{

void Solver::operator()(state_type &x)
{
    MKL_INT size = (MKL_INT)(x.size());

    state_type      a(5 * size);
    state_type      b(size);
    state_type      x_prev(size);
    vector_type_int jia(size + 1);
    vector_type_int jja(5 * size);

    float_type *mkl_a = a.data();  // A
    float_type *mkl_b = b.data();  // b

    MKL_INT *ia = jia.data();
    MKL_INT *ja = jja.data();

    double dt = 0.1;  // time step

    int count = 0;
    int calls = 0;

    for(int i = 0; i < size; i++)
    {
        x_prev[i] = x[i];
    }

    MKL_INT mtype = 11; /* Real unsymmetric matrix */

    /* RHS and solution vectors. */
    float_type mkl_x[size];

    MKL_INT nrhs = 1; /* Number of right hand sides. */

    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];

    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;

    /* Auxiliary variables. */
    double  ddum; /* Double dummy */
    MKL_INT idum; /* Integer dummy. */

    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for(MKL_INT i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0]  = 1;  /* No solver default */
    iparm[1]  = 3;  /* Fill-in reordering from METIS ---- 3 better perf */
    iparm[3]  = 0;  /* No iterative-direct algorithm */
    iparm[4]  = 0;  /* No user fill-in reducing permutation */
    iparm[5]  = 0;  /* Write solution into x */
    iparm[6]  = 0;  /* Not in use */
    iparm[7]  = 2;  /* Max numbers of iterative refinement steps */
    iparm[8]  = 0;  /* Not in use */
    iparm[9]  = 10; /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;  /* Conjugate transposed/transpose solve */
    iparm[12] = 1;  /* Maximum weighted matching algorithm is switched-on
                       (default for non-symmetric) */
    iparm[13] = 0;  /* Output: Number of perturbed pivots */
    iparm[14] = 0;  /* Not in use */
    iparm[15] = 0;  /* Not in use */
    iparm[16] = 0;  /* Not in use */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[19] = 0;  /* Output: Numbers of CG Iterations */

    // iparm[23] = 10;
    // iparm[24] = 2;

    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum   = 1; /* Which factorization to use. */
    msglvl = 0; /* Print statistical information  */
    error  = 0; /* Initialize error flag */

    /* -------------------------------------------------------------------- */
    /* Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for(MKL_INT i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    // Start timer
    // boost::timer::auto_cpu_timer timer;

    int exitf = 0;

    for(double t = dt; t <= (m_t1 + dt * 0.5); t += dt)
    {
        std::cout << "\nt = " << t << ": ";

        int iter;

        for(iter = 0; iter < 20; iter++)
        {

            m_rhs(x, x_prev, b, t, dt);
            /*if(iter == 0)*/ m_jac(x, x_prev, a, jia, jja, t, dt);

            // Solve the equations A*X = B
            // info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, mkl_a, lda, ipiv,
            // mkl_b, ldb); Solution is in b, details of LU factorization are in
            // a

            if(iter == 0)
            {

                /* --------------------------------------------------------------------
                 */
                /* .. Reordering and Symbolic Factorization. This step also
                 * allocates */
                /* all memory that is necessary for the factorization. */
                /* --------------------------------------------------------------------
                 */
                phase = 11;
                PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia,
                        ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
                if(error != 0)
                {
                    printf("\nERROR during symbolic factorization: %d", error);
                    exit(1);
                }
                // printf("\nReordering completed ... ");
                // printf("\nNumber of nonzeros in factors = %d", iparm[17]);
                // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

                /* --------------------------------------------------------------------
                 */
                /* .. Numerical factorization. */
                /* --------------------------------------------------------------------
                 */
                phase = 22;
                PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia,
                        ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
                if(error != 0)
                {
                    printf("\nERROR during numerical factorization: %d", error);
                    exit(2);
                }
                // printf("\nFactorization completed ... ");
            }

            /* --------------------------------------------------------------------
             */
            /* .. Back substitution and iterative refinement. */
            /* --------------------------------------------------------------------
             */
            phase = 33;

            iparm[11] = 0; /* Conjugate transposed/transpose solve */

            PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, mkl_a, ia, ja,
                    &idum, &nrhs, iparm, &msglvl, mkl_b, mkl_x, &error);
            if(error != 0)
            {
                printf("\nERROR during solution: %d", error);
                exit(3);
            }

            calls++;

            double tol = 0.0;
            for(int i = 0; i < size; i++)
            {
                double adiff = std::abs(mkl_x[i]);
                if(adiff > tol)
                    tol = adiff;
                x[i] -= mkl_x[i];
            }

            std::cout << iter << ':' << tol << ' ';
            // std::cout << iter << ':';
            if(tol < m_tol)
                break;

        }  // for iter

        std::cout << '=' << iter << "= dt: " << dt;

        if(exitf)
            break;

        if(iter <= 3)
            dt *= 1.4;
        else if(iter > 6)
            dt /= 1.4;

        if(t + dt > m_t1)
        {
            dt    = m_t1 - t;
            exitf = 1;
        }

        for(int i = 0; i < size; i++)
        {
            x_prev[i] = x[i];
        }

        count++;

    }  // for t

    std::cout << '\n' << calls << " calls\n";

    // Check for the exact singularity
    /*        if(info > 0)
            {
                std::cout << "The diagonal element " << info << " of the
       triangular factor of A is zero,\n"; std::cout << "so that A is singular.
       The solution could not be computed.\n\n"; exit(1);
            }*/

    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, &ddum, ia, ja, &idum,
            &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
}

}  // namespace daecpp_namespace_name
