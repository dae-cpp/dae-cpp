/*
 * Checks a user-defined array used to store a sparse matrix in order to detect
 * issues which could cause problems in routines that require sparse input
 * matrices, such as Intel MKL PARDISO.
 * Based on Intel MKL sparse_matrix_checker() routine.
 */

#include <iostream>
#include <mkl_types.h>
#include <mkl_sparse_handle.h>

#include "time_integrator.h"

namespace daecpp_namespace_name
{

int TimeIntegrator::matrix_checker(sparse_matrix_holder &A, MKL_INT size)
{
    sparse_checker_error_values check_err_val;
    sparse_struct pt;

    sparse_matrix_checker_init(&pt);

    pt.n = size;
    pt.csr_ia = A.ia.data();
    pt.csr_ja = A.ja.data();
    pt.indexing = MKL_ONE_BASED;
    pt.matrix_structure = MKL_GENERAL_STRUCTURE;
    pt.matrix_format = MKL_CSR;
    pt.print_style = MKL_C_STYLE;  // MKL_C_STYLE or MKL_FORTRAN_STYLE
    pt.message_level = MKL_PRINT;  // MKL_PRINT or MKL_NO_PRINT

    check_err_val = sparse_matrix_checker(&pt);

    switch(check_err_val)
    {
    case MKL_SPARSE_CHECKER_SUCCESS:
        // std::cout << "The input array successfully passed all checks.\n";
        return 0;
        break;
    case MKL_SPARSE_CHECKER_NON_MONOTONIC:
        std::cout << "The input array is not 0 or 1 based (ia[0] is not 0 or 1) or elements of ia are not in non-decreasing order as required.\n";
        break;
    case MKL_SPARSE_CHECKER_OUT_OF_RANGE:
        std::cout << "The value of the ja array is lower than the number of the first column or greater than the number of the last column.\n";
        break;
    case MKL_SPARSE_CHECKER_NONORDERED:
        std::cout << "The elements of the ja array are not in non-decreasing order in each row as required.\n";
        break;
    default:
        std::cout << "Unknown error in ia and/or ja arrays.\n";
        break;
    }

    std::cout << "Matrix check details: (" << pt.check_result[0] << ", " << pt.check_result[1] << ", " << pt.check_result[2] << ")\n";

    return 1;
}

}
