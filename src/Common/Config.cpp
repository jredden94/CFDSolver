#include "Config.hpp"

Config* Config::instance = nullptr;
once_flag Config::onceFlag;

Config::Config() { }
Config::~Config() { delete instance; }

Config& Config::GetConfig() { 
    call_once(onceFlag, &Config::Init);
    return *instance;
}

void Config::Init() {
    instance = new Config;
}

void Config::SetDefault() {
    grid_file = "flatplate.su2";
    solver = Config::SolverType::Euler;
    convFlux = Config::ConvFlux::Roe;
    venkat_lim_coeff = 0.0001;
    nDim = 2;
    nVar = 4;
    rho_inf = 1.13831;
    vel_inf.resize(3);
    aoa = 1.0;
    vel_inf[0] = 69.1645; //0.8 * cosf(aoa * M_PI / 180.0);
    vel_inf[1] = 0.0; //0.8 * sinf(aoa * M_PI / 180.0);
    vel_inf[2] = 0.0;
    gamma = 1.4;
    p_inf = 97250.; //rho_inf / gamma; //1.0 / gamma; //100000; //1.0 / gamma; // nondimensionalization such that u_inf = M
    gamma_m1 = gamma - 1.0;

    roe_diss_coeff = 1.0;

    cfl_min = 100.0;
    cfl_start = 100.0;
    cfl_max = 10000.0;
    cfl_factor_up = 1.2;
    cfl_factor_down = 0.9;
    adaptCFL = true;

    convergence_tol = 1e-11;
    max_iter = 5000;
    min_iter = 100;
    max_lin_solver_iter = 8;
    lin_solver_tol = 1e-12;
    implicit = true;
    limiterType = Limiter::Type::Venkatakrishnan;
    muscl = true;
    isLimited = true;
    precond = Config::Preconditioner::Jacobi;

    isViscous = true;
    viscosity = 1.83751e-05;
}

const string& Config::GetGridFilename(void) const { return grid_file; }
const Config::ConvFlux& Config::GetConvFluxScheme(void) const { return convFlux; }
const Config::SolverType& Config::GetSolverType(void) const { return solver; }
const vector<double>& Config::GetFreestreamVelocity(void) const { return vel_inf; }
const double& Config::GetFreestreamDensity(void) const { return rho_inf; }
const double& Config::GetFreestreamPressure(void) const { return p_inf; }
const double& Config::GetFreestreamMach(void) const { return M_inf; }
const double& Config::GetAngleOfAttack(void) const { return aoa; }
const double& Config::GetGamma(void) const { return gamma; }
const double& Config::GetGasConstant(void) const { return gasConstant; }
const double& Config::GetCFL(void) const { return cfl_start; }
const double& Config::GetMinCFL(void) const { return cfl_min; }
const double& Config::GetMaxCFL(void) const { return cfl_max; }
const double& Config::GetCFL_FactorUp(void) const { return cfl_factor_up; }
const double& Config::GetCFL_FactorDown(void) const { return cfl_factor_down; }
const bool& Config::AdaptiveCFL(void) const { return adaptCFL; }
const bool& Config::IsSteady(void) const { return steady; }
const unsigned short& Config::GetNumDims(void) const { return nDim; }
const unsigned short& Config::GetNumVars(void) const { return nVar; }
const unsigned short& Config::GetNumEqn(void) const { return nEqn; }
const double& Config::GetConvergenceTolerance(void) const { return convergence_tol; }
const unsigned long& Config::GetLinearSolverMaxIterations(void) const { return max_lin_solver_iter; }
const unsigned long& Config::GetMaxIterations(void) const { return max_iter; }
const unsigned long& Config::GetMinIterations(void) const { return min_iter; }
const double& Config::GetLinearSolverTolerance(void) const { return lin_solver_tol; };
const bool Config::IsImplicit(void) const { return implicit; }
const bool Config::MUSCL(void) const { return muscl; }
const bool Config::IsLimited(void) const { return isLimited; }
const Limiter::Type& Config::GetLimiterType(void) const { return limiterType; }
const double& Config::GetLimiterCoefficient(void) const { return venkat_lim_coeff; }
const Config::Preconditioner Config::GetPreconditioner(void) const { return precond; }
const double& Config::GetViscosity(void) const { return viscosity; }
const bool& Config::IsViscous(void) const { return isViscous; }
const double& Config::GetRoeDissipationCoefficient(void) const { return roe_diss_coeff; }

void Config::SetNumEqn(const unsigned short& nEqn) { this->nEqn = nEqn; }
