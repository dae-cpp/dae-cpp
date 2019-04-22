/*
 * TODO: Description of the class
 */

#pragma once

#include "../../src/jacobian.h"

class MyJacobian : public Jacobian
{
//    MyParams m_p;

public:
//    MyRHS(MyParams p) : RHS(), m_p(p) {}
    //void operator()(const state_type& x, state_type& f, const double t);
    MyJacobian(RHS &rhs) : Jacobian(rhs) {}
    void operator() (state_type &x, state_type &J, vector_type_int &ia, vector_type_int &ja, const double t, const double dt);
};
