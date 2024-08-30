#pragma once

#include "SolMatrix.hpp"
#include "SolVector.hpp"
#include "Precond.hpp"
#include "../Common/Config.hpp"
#include "../Common/AD.hpp"

class GMRES {
    public:
        GMRES();
        ~GMRES();

        void Solve(const SolMatrix &mat, SolVector &b, SolVector &x, const unique_ptr<Precond> &precond);

    private:
        zdouble tol;
        unsigned long max_iter;
        zdouble *g, *sn, *cs, *y, *Hess;
        vector<SolVector> w, z;

        void GramSchmidt(const unsigned long i);
        void ApplyGivens(zdouble s, zdouble c, zdouble &h1, zdouble &h2);
        void GenerateGivens(zdouble &dx, zdouble &dy, zdouble &s, zdouble &c);
        void SolveReduced(unsigned long i, const zdouble *hess, const zdouble *rhs, zdouble *x);
        void Reset(void);
};
