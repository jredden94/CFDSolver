#include "GMRES.hpp"

GMRES::GMRES() { 
    Config *config = &Config::GetConfig();
    max_iter = config->GetLinearSolverMaxIterations();
    tol = config->GetLinearSolverTolerance();
    w.resize(max_iter + 1);
    z.resize(max_iter + 1);
    
    g = new zdouble[max_iter+1];
    sn = new zdouble[max_iter+1];
    cs = new zdouble[max_iter+1];
    y = new zdouble[max_iter];
    Hess = new zdouble [(max_iter+1) * max_iter];

}
GMRES::~GMRES() { 
    delete[] g;
    delete[] sn;
    delete[] cs;
    delete[] y;
    delete[] Hess;
}

void GMRES::Solve(const SolMatrix &mat, SolVector &b, SolVector &x, const unique_ptr<Precond> &precond) { 
    Reset();

    zdouble norm0 = b.Norm();

    // Assuming x0 = 0
    //w[0].SetZeroes();
    w[0] -= b;

    zdouble beta = w[0].Norm();
    w[0] /= -beta;
    g[0] = beta;

    for (auto i = 0ul; i < max_iter; i++) {
        if (beta < tol * norm0) break;

        precond->Apply(w[i], z[i]);
        mat.MatrixVecMult(z[i], w[i+1]);

        GramSchmidt(i);

#define H(I,J) Hess[I * max_iter + J]
        for (auto k = 0ul; k < i; k++)
            ApplyGivens(sn[k], cs[k], H(k,i), H(k+1,i));

        GenerateGivens(H(i,i), H(i+1,i), sn[i], cs[i]);
        ApplyGivens(sn[i], cs[i], g[i], g[i+1]);
#undef H

        beta = fabs(g[i+1]);

        SolveReduced(i, Hess, g, y);

        for (auto k = 0ul; k < i; k++) x.PlusEqualScalarMult(z[k], y[k]);
    }
}

void GMRES::GramSchmidt(unsigned long i) { 
    const zdouble reorth = 0.98;
    zdouble nrm = w[i+1].SquaredNorm();
    zdouble thr = nrm * reorth;

    for (int k = 0; k < i+1; k++) {
        zdouble prod = w[i+1].DotProd(w[k]);
        zdouble h_ki = prod;
        w[i+1].MinusEqualScalarMult(w[k], prod);
        
        if (prod * prod > thr) {
            prod = w[i+1].DotProd(w[k]);
            h_ki += prod;
            w[i+1].MinusEqualScalarMult(w[k], prod);
        }
        Hess[k * max_iter + i] = h_ki;

        nrm -= h_ki * h_ki;
        nrm = max(nrm, 0.0);
        thr = nrm * reorth;
    }

    nrm = w[i+1].Norm();
    Hess[max_iter * (i+1) + i] = nrm;

    w[i+1] /= nrm;
}

void GMRES::ApplyGivens(zdouble s, zdouble c, zdouble &h1, zdouble &h2) {
    zdouble temp = c * h1 + s * h2;
    h2 = c * h2 - s * h1;
    h1 = temp;
}

zdouble Sign(zdouble x, zdouble y) {
    if (y == 0.0) return 0.0;
    return fabs(x) * (y < 0.0 ? -1.0 : 1.0);
}

void GMRES::GenerateGivens(zdouble &dx, zdouble &dy, zdouble &s, zdouble &c) { 
    if (dx == 0 && dy == 0) {
        c = 1.0;
        s = 0.0;
    }
    else if (fabs(dy) > fabs(dx)) {
        zdouble temp = dx / dy;
        dx = sqrt(1.0 + temp * temp);
        s = Sign(1.0 / dx, dy);
        c = temp * s;
    }
    else if (fabs(dy) <= fabs(dx)) {
        zdouble temp = dy / dx;
        dy = sqrt(1.0 + temp * temp);
        c = Sign(1.0 / dy, dx);
        s = temp * c;
    }
    else {
        dx = 0.0;
        dy = 0.0;
        c = 1.0;
        s = 0.0;
    }

    dx = fabs(dx * dy);
    dy = 0.0;
}

void GMRES::SolveReduced(unsigned long n, const zdouble *hess, const zdouble *rhs, zdouble *x) {
    for (int i = 0; i < n; i++) x[i] = rhs[i];

    for (int i = n-1; i >=0; i--) {
        x[i] /= hess[i * max_iter + i];
        for (int j = i - 1; j >= 0; j--) {
            x[j] -= Hess[j * max_iter + i] * x[i];
        }
    }
}

void GMRES::Reset(void) {
    for (auto i = 0ul; i < max_iter + 1; i++) {
        z[i].SetZeroes();
        w[i].SetZeroes();
        g[i] = 0;
        sn[i] = 0;
        cs[i] = 0;
        for (auto j = 0ul; j < max_iter; j++) Hess[i * max_iter + j] = 0;
    }
    for (auto i = 0ul; i < max_iter; i++) y[i] = 0;
}
