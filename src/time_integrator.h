/*
 * Numerical time integrator class definition
 */

#pragma once

#include <mkl_spblas.h>

#include "typedefs.h"
#include "RHS.h"
#include "jacobian.h"
#include "mass_matrix.h"
#include "solver_options.h"

namespace daecpp_namespace_name
{

class TimeIntegrator
{
    RHS &m_rhs;

    Jacobian &m_jac;

    MassMatrix &m_mass;

    SolverOptions &m_opt;

    // clang-format off
    const double BDF_COEF[6][8] = {
        {  1.0,   -1.0,   0.0,    0.0,   0.0,   0.0,  0.0,  1.0},
        {  3.0,   -4.0,   1.0,    0.0,   0.0,   0.0,  0.0,  2.0},
        { 11.0,  -18.0,   9.0,   -2.0,   0.0,   0.0,  0.0,  6.0},
        { 25.0,  -48.0,  36.0,  -16.0,   3.0,   0.0,  0.0, 12.0},
        {137.0, -300.0, 300.0, -200.0,  75.0, -12.0,  0.0, 60.0},
        {147.0, -360.0, 450.0, -400.0, 225.0, -72.0, 10.0, 60.0}};

    const double ALPHA_COEF[6] =
        {1.0, 1.5, 11.0/6.0, 25.0/12.0, 137.0/60.0, 2.45};
    // clang-format on

    // The first time step will be performed using BDF-1
    int m_scheme = 1;

    // Temporary Jacobian matrix holder
    sparse_matrix_holder m_J;

    // Mass matrix container
    sparse_matrix_holder m_M;

    // Descriptor of sparse mass matrix properties
    struct matrix_descr m_descrA;

    // Structure with sparse mass matrix stored in CSR format
    sparse_matrix_t m_csrA;

    /*
     * Sparse matrix checker
     */
    int m_matrix_checker(sparse_matrix_holder &A, MKL_INT size);

    /*
     * Performs matrix-matrix addition: C = alpha*A + B.
     * Replaces deprecated Intel MKL mkl_dcsradd() function.
     */
    void m_matrix_add(const float_type alpha, const sparse_matrix_holder &A,
                      const sparse_matrix_holder &B, sparse_matrix_holder &C);

public:
    TimeIntegrator(RHS &rhs, Jacobian &jac, MassMatrix &mass,
                   SolverOptions &opt);

    ~TimeIntegrator() { mkl_sparse_destroy(m_csrA); }

    void set_scheme(int scheme) { m_scheme = scheme; }

    void integrate(sparse_matrix_holder &J, state_type &b, const state_type &x,
                   const state_type_matrix &x_prev, const double t,
                   const double dt[], const bool do_jac);
};

}  // namespace daecpp_namespace_name
