#include "BiCGSTAB.hpp"
#include "SolVector.hpp"

BiCGSTAB::BiCGSTAB() { }

BiCGSTAB::~BiCGSTAB() { }
void BiCGSTAB::Solve(const SolMatrix &mat, SolVector &b, SolVector &x, const unique_ptr<Precond> &precond) {

    // Assuming x is initially zero vector
    r = b;
    r0 = r;

    double alpha = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
    p.SetZeroes();
    z.SetZeroes();
    x.SetZeroes();
    Ax.SetZeroes();
    v.SetZeroes();
    double norm0 = b.Norm();

    for (auto i = 0ul; i < max_iter; i++) {
        rho_prime = rho;
        rho = r.DotProd(r0);

        double beta = (rho / rho_prime) * (alpha / omega);
        p.MinusEqualScalarMult(v, omega);
        p *= beta;
        p += r;

        precond->Apply(p, z);
        mat.MatrixVecMult(z, v);

        double r_o_v = r0.DotProd(v);
        alpha = rho / r_o_v;

        x.PlusEqualScalarMult(z, alpha);
        r.MinusEqualScalarMult(v, alpha);

        precond->Apply(r, z);
        mat.MatrixVecMult(z, Ax);

        omega = Ax.SquaredNorm();
        if (omega == 0) break;
        omega = Ax.DotProd(r) / omega;

        x.PlusEqualScalarMult(z, omega);
        r.MinusEqualScalarMult(Ax, omega);

        double r_norm = r.Norm();
        if (r_norm < tol * norm0) break;
    }
}
