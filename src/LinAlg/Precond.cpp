#include "Precond.hpp"

/* Preconditioner Base Class */
Precond::Precond() { }
Precond::~Precond() { }

std::unique_ptr<Precond> Precond::GetPreconditioner() {
    Config *config = &Config::GetConfig();
    Config::Preconditioner type = config->GetPreconditioner();
    switch (type) {
        case (Config::Preconditioner::Jacobi) : return make_unique<Jacobi>(); break;
    }
}

/* Jacobi Block Preconditioner */
Jacobi::Jacobi() { }
Jacobi::~Jacobi() { }

void Jacobi::Build(const SolMatrix &jac_mat) {
    Config *config = &Config::GetConfig();
    unsigned long nEqn = config->GetNumEqn();
    unsigned long nVar = config->GetNumVars();

    for (auto i = 0ul; i < nEqn; i++) {
        zdouble *jacobi_blk = jacobi.GetBlock(i);
        const zdouble *jac_diag = jac_mat.GetBlock(i, i);
        blas::InverseMatrix(jac_diag, jacobi_blk, nVar);
    }
}

void Jacobi::Apply(const SolVector &vec, SolVector &prod) const {
    jacobi.MatrixVecMult(vec, prod);
}

