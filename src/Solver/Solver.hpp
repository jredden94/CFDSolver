#pragma once

#include <memory>
#include "../Common/Config.hpp"
#include "../Common/AD.hpp"
#include "../DualGrid/Grid.hpp"
#include "../DualGrid/Edge.hpp"
#include "../DualGrid/Node.hpp"
#include "../Numerics/Convection/Convection.hpp"
#include "../LinAlg/SolVector.hpp"
#include "../LinAlg/SolMatrix.hpp"
#include "../LinAlg/BiCGSTAB.hpp"
#include "../LinAlg/GMRES.hpp"
#include "../LinAlg/GuassSiedel.hpp"
#include "../Gradient/LeastSquaresGrad.hpp"
#include "../LinAlg/Precond.hpp"
#include "../Numerics/Viscous/Viscous.hpp"

using namespace std;

class Solver {
    public:
        Solver(const Solver&) = delete;
        virtual ~Solver();
        static std::unique_ptr<Solver> CreateSolver(Grid*);
        virtual void Solve() = 0;

        Solver& operator=(const Solver&) = delete;
        void Init(void);

        const SolVector& GetPrimitives(void) const;

    protected:
        Solver(Config::SolverType);
        LeastSquaresGrad lsq_grad;
        zdouble *v_i, *v_j;
        void UpdatePrimitiveVars(void);
        void PrintResiduals(const unsigned long current_iter, const zdouble cL = 0.0, const zdouble cD = 0.0) const;
        void AdaptCFL(void);

        Config *config;
        mutable bool converged = false;

        bool implicit, muscl;
        unsigned short nVar;
        unsigned short nDim;
        unsigned short nEqn;

        SolVector delU, residual, conVar, primVar;
        SolMatrix jacMat;

        zdouble cfl, gamma;
        vector<zdouble> waveSpeed, dt;
        Grid *grid;
        unique_ptr<Convection> conv_flux;
        bool viscous;
        Viscous visc_flux;

        zdouble tol;
        unique_ptr<Precond> precond;
        BiCGSTAB sysSolver;
        GSS gsSolver;
        GMRES gmres;
        mutable vector<zdouble> res_norm;
        vector<zdouble> nonlinear_res;
        unsigned long nonlinear_res_counter = 0, nonlinear_res_max_count = 25;

};

class Euler : public Solver {
    public:
        Euler();
        ~Euler() override;

        void Solve() override;
        vector<zdouble> GetBoundaryState(const Boundary::BCType &type, 
                const zdouble *leftState, const vector<zdouble> &norm) const;
};
