/*
 * TODO: Description of the class
 */

#pragma once

#include "../../src/mass_matrix.h"
//#include "perovskite_parameters.h"

class MyMassMatrix : public MassMatrix
{
    //MyParams m_p;

public:
    //MyMassMatrix(MyParams p) : RHS(), m_p(p) {}
    void operator()(const state_type& x, state_type& f, const double t);
};
