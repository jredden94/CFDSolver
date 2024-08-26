#pragma once

#include <vector>
#include <mutex>
#include <iostream>
#include "../Limiters/Limiter.hpp"

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
        const vector<double>& GetFreestreamVelocity(void) const;
        const double& GetFreestreamDensity(void) const;
        const double& GetFreestreamPressure(void) const;
        const double& GetFreestreamMach(void) const;
        const double& GetAngleOfAttack(void) const;
        const double& GetGamma(void) const;
        const double& GetGasConstant(void) const;
        const double& GetCFL(void) const;
        const double& GetMinCFL(void) const;
        const double& GetMaxCFL(void) const;
        const double& GetCFL_FactorUp(void) const;
        const double& GetCFL_FactorDown(void) const;
        const bool& AdaptiveCFL(void) const;
        const bool& IsSteady(void) const;
        const unsigned short& GetNumDims(void) const;
        const unsigned short& GetNumVars(void) const;
        const unsigned short& GetNumEqn(void) const;
        const double& GetConvergenceTolerance(void) const;
        const unsigned long& GetLinearSolverMaxIterations(void) const;
        const unsigned long& GetMaxIterations(void) const;
        const unsigned long& GetMinIterations(void) const;
        const double& GetLinearSolverTolerance(void) const;
        const Limiter::Type& GetLimiterType(void) const;
        const double& GetLimiterCoefficient(void) const;
        const bool IsImplicit(void) const;
        const bool MUSCL(void) const;
        const bool IsLimited(void) const;
        const Preconditioner GetPreconditioner(void) const;
        const double& GetViscosity(void) const;
        const bool& IsViscous(void) const;
        const double& GetRoeDissipationCoefficient(void) const;

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
        double venkat_lim_coeff;

        vector<double> vel_inf;
        double rho_inf, p_inf, M_inf;
        double aoa; // angle of attack
        double gamma, gamma_m1, gasConstant;
        bool steady, implicit, muscl, isLimited;
        unsigned short nDim, nVar, nEqn;

        double viscosity;
        bool isViscous;

        bool adaptCFL;
        double cfl_start, cfl_min, cfl_max, cfl_factor_up, cfl_factor_down;
        double convergence_tol, lin_solver_tol;
        unsigned long max_lin_solver_iter, max_iter, min_iter;

        double roe_diss_coeff;
};
