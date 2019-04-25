/*
 * TODO: Description of the class
 */

#pragma once

#include "../../src/mass_matrix.h"

class MyMassMatrix : public daecpp::MassMatrix
{
    const int m_N;

public:
    MyMassMatrix(const int N) : daecpp::MassMatrix(), m_N(N) {}
    void operator()(daecpp::sparse_matrix_holder &M);
};
