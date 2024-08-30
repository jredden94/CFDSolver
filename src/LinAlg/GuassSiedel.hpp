#pragma once

#include "SolMatrix.hpp"
#include "SolVector.hpp"
#include "BLAS.hpp"
#include "../Common/AD.hpp"

class GSS {
    public:
        GSS();
        ~GSS();

        void Solve(const SolMatrix &mat, const SolVector &b, SolVector &x, const unsigned short nSweeps = 15);
        void LUSolve(const zdouble *mat, const zdouble *b, zdouble *x);
    private:
        unsigned short vec_len, mat_len, nEqn;
        void LUDecomposition(const zdouble *block);
        zdouble *lu_upper, *lu_lower, *b, *y, *x;
};
