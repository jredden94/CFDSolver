#pragma once

#include "../../Common/Config.hpp"
#include "../../Common/AD.hpp"

class Viscous {
    public:
        Viscous(void);
        ~Viscous(void);

        void SetStates(const zdouble *v_i, const zdouble *v_j, const zdouble *norm,
                const zdouble area, const zdouble edge_len, const zdouble *vel_grad_i, const zdouble *vel_grad_j);
        void ComputeResidual(void);

        const zdouble* Flux(void) const;
        const zdouble* JacI(void) const;
        const zdouble* JacJ(void) const;
        const zdouble* StressTensor(void) const;
        void ComputeStressTensor(const zdouble *vel_grad, const zdouble viscosity, zdouble *tau );

    private:
        void ComputeViscousFlux(void);
        void ComputeTauJacobian(void);
        void ComputeJacobians(void);

        unsigned short nVar, nDim;

        zdouble viscosity;
        zdouble *tau, *visc_flux;
        zdouble *v_i, *v_j, *v_mean, *norm, *coord_i, *coord_j;
        zdouble area, edge_len;
        zdouble *vel_grad_i, *vel_grad_j, *vel_grad_mean;


        bool implicit;
        zdouble *tau_jac, *jac_i, *jac_j;
};
