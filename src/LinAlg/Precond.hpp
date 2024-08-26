#pragma once

#include "SolMatrix.hpp"
#include "SolVector.hpp"
#include "DiagonalMatrix.hpp"

class Precond {
    public:
        enum class Preconditioners { Jacobi };
        Precond();
        virtual ~Precond();

        virtual void Build(const SolMatrix &jac_mat) = 0;
        virtual void Apply(const SolVector &vec, SolVector &prod) const = 0;
        static std::unique_ptr<Precond> GetPreconditioner(void);
};

class Jacobi : public Precond {
    public:
        Jacobi();
        ~Jacobi() override;

        void Build(const SolMatrix &jac_mat) override;
        void Apply(const SolVector &vec, SolVector &prod) const override;

    private:
        DiagonalMatrix jacobi;
};
