#pragma once

#include "SolMatrix.hpp"
#include "SolVector.hpp"
#include "Precond.hpp"
#include "../Common/Config.hpp"

class GMRES {
    public:
        GMRES();
        ~GMRES();

        void Solve(const SolMatrix &mat, SolVector &b, SolVector &x, const unique_ptr<Precond> &precond);

    private:
        double tol;
        unsigned long max_iter;
        double *g, *sn, *cs, *y, *Hess;
        vector<SolVector> w, z;

        void GramSchmidt(const unsigned long i);
        void ApplyGivens(double s, double c, double &h1, double &h2);
        void GenerateGivens(double &dx, double &dy, double &s, double &c);
        void SolveReduced(unsigned long i, const double *hess, const double *rhs, double *x);
        void Reset(void);
};
