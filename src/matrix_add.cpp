/*
 * Performs matrix-matrix addition: C = alpha*A + B.
 * Replaces deprecated Intel MKL mkl_dcsradd() function.
 */

#include "time_integrator.h"

namespace daecpp_namespace_name
{

void TimeIntegrator::matrix_add(const float_type alpha,
                                const sparse_matrix_holder &A,
                                const sparse_matrix_holder &B,
                                sparse_matrix_holder &C)
{
// For compatibility with mkl_dcsradd()
#ifdef DAE_FORTRAN_STYLE
#define FORTRAN_STYLE 1
#else
#define FORTRAN_STYLE 0
#endif

    MKL_INT ia    = 1;
    MKL_INT ja    = 0;
    MKL_INT ib    = 1;
    MKL_INT jb    = 0;
    MKL_INT a_row = 0;
    MKL_INT b_row = 0;
    MKL_INT ic    = 0;
    MKL_INT jc    = 0;
    MKL_INT c_row = -1;

    const MKL_INT sizeA = A.ja.size();
    const MKL_INT sizeB = B.ja.size();

    while(ja < sizeA && jb < sizeB)
    {
        // ib and jb point to the element before
        if(b_row < a_row || (b_row == a_row && B.ja[jb] < A.ja[ja]))
        {
            C.A[jc]  = B.A[jb];
            C.ja[jc] = B.ja[jb];

            if(c_row < b_row)
            {
                C.ia[ic++] = jc + FORTRAN_STYLE;
                c_row++;
            }

            jb++;
            jc++;

            if(B.ia[ib] == jb + FORTRAN_STYLE)
            {
                b_row++;
                ib++;
            }

            continue;
        }

        // ib and jb point to the same position as ia and ja
        if(b_row == a_row && B.ja[jb] == A.ja[ja])
        {
            C.A[jc]  = alpha * A.A[ja] + B.A[jb];
            C.ja[jc] = A.ja[ja];

            if(c_row < a_row)
            {
                C.ia[ic++] = jc + FORTRAN_STYLE;
                c_row++;
            }

            ja++;
            jb++;
            jc++;

            if(A.ia[ia] == ja + FORTRAN_STYLE)
            {
                a_row++;
                ia++;
            }

            if(B.ia[ib] == jb + FORTRAN_STYLE)
            {
                b_row++;
                ib++;
            }

            continue;
        }

        // ib and jb point to the element after
        if(b_row > a_row || (b_row == a_row && B.ja[jb] > A.ja[ja]))
        {
            C.A[jc]  = alpha * A.A[ja];
            C.ja[jc] = A.ja[ja];

            if(c_row < a_row)
            {
                C.ia[ic++] = jc + FORTRAN_STYLE;
                c_row++;
            }

            ja++;
            jc++;

            if(A.ia[ia] == ja + FORTRAN_STYLE)
            {
                a_row++;
                ia++;
            }

            continue;
        }
    }

    while(ja < sizeA)
    {
        C.A[jc]  = alpha * A.A[ja];
        C.ja[jc] = A.ja[ja];

        if(c_row < a_row)
        {
            C.ia[ic++] = jc + FORTRAN_STYLE;
            c_row++;
        }

        ja++;
        jc++;
    }

    while(jb < sizeB)
    {
        C.A[jc]  = B.A[jb];
        C.ja[jc] = B.ja[jb];

        if(c_row < b_row)
        {
            C.ia[ic++] = jc + FORTRAN_STYLE;
            c_row++;
        }

        jb++;
        jc++;
    }

    C.ia[ic] = jc + FORTRAN_STYLE;
}

}  // namespace daecpp_namespace_name
