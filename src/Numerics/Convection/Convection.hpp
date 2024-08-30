#pragma once

#include <memory>
#include <vector>
#include <cmath>
#include "../../Common/Config.hpp"
#include "../../Common/AD.hpp"

using namespace std;

class Convection {
    public:
        Convection();
        Convection(Config::ConvFlux);
        virtual ~Convection();

        virtual void ComputeFlux(void) = 0;

        void SetStates(const vector<zdouble>& v_i, const vector<zdouble>& v_j, 
                const vector<zdouble>& norm, const zdouble &area_mag);

        void SetStates(const zdouble* v_i, const zdouble* v_j, 
                const zdouble* areaVector, const zdouble &area_mag);

        const vector<zdouble>& Flux(void);
        const zdouble* JacI(void);
        const zdouble* JacJ(void);
        const zdouble& MaxWaveSpeed(void);

        static unique_ptr<Convection> CreateConvFlux(Config *config);

        void PrintIJacobianInfo(void);
        void PrintJJacobianInfo(void);

        void ComputeJacobian(const zdouble *vel, const zdouble &vel_sqr, 
                const zdouble &proj_vel, const zdouble *norm, 
                const zdouble &energy, zdouble *jac, const zdouble scale = 0.5);

    protected:

        void ComputeFlux(const zdouble &rho, const zdouble *vel, const zdouble &enthalpy,
        const zdouble *norm, const zdouble *flux);
        void ComputePTensor(const zdouble &rho, const zdouble *vel,
                const zdouble &soundSpeed, const zdouble *norm, zdouble *p_tensor);
        void ComputeInversePTensor(const zdouble &rho, const zdouble *vel, 
                const zdouble &soundspeed, const zdouble *norm, zdouble *inv_p);

        Config *config;
        unsigned short nVar, nDim;
        bool implicit;

        vector<zdouble> norm, flux, eigen; 
        zdouble maxWaveSpeed, entropy_fix_coeff; 
        zdouble area;

        zdouble gamma, gamma_minus_one;

        unsigned short iDim, jDim, kDim, iVar, jVar, kVar;

        // Left state values
        zdouble *v_i, *u_i;
        vector<zdouble> vel_i, flux_i; 
        zdouble rho_i, p_i, enthalpy_i, energy_i;
        zdouble vel_sqr_i, proj_vel_i, soundSpeed_i;

        // Right state values
        zdouble *v_j, *u_j;
        vector<zdouble> vel_j, flux_j; 
        zdouble rho_j, p_j, enthalpy_j, energy_j;
        zdouble vel_sqr_j, proj_vel_j, soundSpeed_j;

        // Roe Averages
        vector<zdouble> vel_roe;
        zdouble rhoR, rho, p, enthalpy, energy;
        zdouble proj_vel, vel_sqr_roe, soundSpeed, soundSpeed2;

        // Jacobians
        zdouble *jac_i, *jac_j;
};

class Roe : public Convection {
    public:
        Roe();
        ~Roe() override;

        void ComputeFlux() override;

    private:
        zdouble *p_tensor, *p_inv, *del_u;
        zdouble p_lam_pinv;
        zdouble diss_coeff;
};

class HLLC : public Convection {
    public:
        HLLC();
        ~HLLC() override;

        void ComputeFlux() override;

    private:
        zdouble sL, sR, sM;
        zdouble rho_m, pStar, rhoSL, rhoSR;

        zdouble *inter_state;
        zdouble eStar, omega, omegaSM;
        zdouble *dp_du_i, *dsM_du, *drhoStar_du, *dpStar_du, *deStar_du;
};
