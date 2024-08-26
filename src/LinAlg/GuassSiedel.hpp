#pragma once

#include "SolMatrix.hpp"
#include "SolVector.hpp"
#include "BLAS.hpp"

class GSS {
    public:
        GSS();
        ~GSS();

        void Solve(const SolMatrix &mat, const SolVector &b, SolVector &x, const unsigned short nSweeps = 15);
        void LUSolve(const double *mat, const double *b, double *x);
    private:
        unsigned short vec_len, mat_len, nEqn;
        void LUDecomposition(const double *block);
        double *lu_upper, *lu_lower, *b, *y, *x;
};
