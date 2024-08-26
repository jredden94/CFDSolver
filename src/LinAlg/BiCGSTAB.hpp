#pragma once

#include "../Common/Config.hpp"
#include "SolMatrix.hpp"
#include "SolVector.hpp"
#include "Precond.hpp"

class BiCGSTAB {
    public:
        BiCGSTAB();
        ~BiCGSTAB();

        void Solve(const SolMatrix &mat, SolVector &b, SolVector &x, const unique_ptr<Precond> &precond);
    private:
        SolVector Ax, r0, r, p, v, z;
        unsigned short max_iter = 50;
        double tol = 0.001;
};
