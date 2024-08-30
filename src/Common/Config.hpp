#pragma once

#include <vector>
#include <mutex>
#include <iostream>
#include "../Limiters/Limiter.hpp"
#include "AD.hpp"

using namespace std;

class Config {
    public:
        enum class ConvFlux { Roe, HLLC };
        enum class SolverType { Euler };
        enum class OutputFile { VTK, TecPlot };
        enum class Preconditioner { Jacobi };

        static Config& GetConfig(void);

        ~Config();
        
        void SetDefault(void);

        const string& GetGridFilename(void) const;
        const ConvFlux& GetConvFluxScheme(void) const;
        const SolverType& GetSolverType(void) const;
        const vector<zdouble>& GetFreestreamVelocity(void) const;
        const zdouble& GetFreestreamDensity(void) const;
        const zdouble& GetFreestreamPressure(void) const;
        const zdouble& GetFreestreamMach(void) const;
        const zdouble& GetAngleOfAttack(void) const;
        const zdouble& GetGamma(void) const;
        const zdouble& GetGasConstant(void) const;
        const zdouble& GetCFL(void) const;
        const zdouble& GetMinCFL(void) const;
        const zdouble& GetMaxCFL(void) const;
        const zdouble& GetCFL_FactorUp(void) const;
        const zdouble& GetCFL_FactorDown(void) const;
        const bool& AdaptiveCFL(void) const;
        const bool& IsSteady(void) const;
        const unsigned short& GetNumDims(void) const;
        const unsigned short& GetNumVars(void) const;
        const unsigned short& GetNumEqn(void) const;
        const zdouble& GetConvergenceTolerance(void) const;
        const unsigned long& GetLinearSolverMaxIterations(void) const;
        const unsigned long& GetMaxIterations(void) const;
        const unsigned long& GetMinIterations(void) const;
        const zdouble& GetLinearSolverTolerance(void) const;
        const Limiter::Type& GetLimiterType(void) const;
        const zdouble& GetLimiterCoefficient(void) const;
        const bool IsImplicit(void) const;
        const bool MUSCL(void) const;
        const bool IsLimited(void) const;
        const Preconditioner GetPreconditioner(void) const;
        const zdouble& GetViscosity(void) const;
        const bool& IsViscous(void) const;
        const zdouble& GetRoeDissipationCoefficient(void) const;
        const zdouble& GetReferenceArea(void) const;

        void SetNumEqn(const unsigned short&);

        Config(const Config&) = delete;
        Config& operator+(const Config&) = delete;

    private:
        Config();
        static void Init();
        static Config *instance;
        static once_flag onceFlag;

        string grid_file;

        ConvFlux convFlux;
        SolverType solver;
        Preconditioner precond;

        Limiter::Type limiterType;
        zdouble venkat_lim_coeff;

        vector<zdouble> vel_inf;
        zdouble rho_inf, p_inf, M_inf;
        zdouble aoa; // angle of attack
        zdouble gamma, gamma_m1, gasConstant;
        bool steady, implicit, muscl, isLimited;
        unsigned short nDim, nVar, nEqn;

        zdouble viscosity;
        bool isViscous;

        bool adaptCFL;
        zdouble cfl_start, cfl_min, cfl_max, cfl_factor_up, cfl_factor_down;
        zdouble convergence_tol, lin_solver_tol;
        unsigned long max_lin_solver_iter, max_iter, min_iter;

        zdouble roe_diss_coeff;

        zdouble ref_area;
};
