#include "Viscous.hpp"

Viscous::Viscous() { 
    Config *config = &Config::GetConfig(); 
    viscosity = config->GetViscosity();
    nVar = config->GetNumVars();
    nDim = config->GetNumDims();
    implicit = config->IsImplicit();

    v_i = new zdouble[nVar];
    v_j = new zdouble[nVar];
    v_mean = new zdouble[nVar];
    norm = new zdouble[nDim];
    vel_grad_i = new zdouble[nDim * nDim];
    vel_grad_j = new zdouble[nDim * nDim];
    vel_grad_mean = new zdouble[nDim * nDim];
    visc_flux = new zdouble[nVar];
    tau = new zdouble[nDim * nDim];

    if (implicit) {
        tau_jac = new zdouble[nDim * nVar];
        jac_i = new zdouble[nVar * nVar];
        jac_j = new zdouble[nVar * nVar];
    }
}

Viscous::~Viscous() { 
    delete[] v_i;
    delete[] v_j;
    delete[] v_mean;
    delete[] norm;
    delete[] vel_grad_i;
    delete[] vel_grad_j;
    delete[] vel_grad_mean;
    delete[] visc_flux;
    delete[] tau;
    if (implicit) {
        delete[] tau_jac;
        delete[] jac_i;
        delete[] jac_j;
    }
}

const zdouble* Viscous::Flux(void) const { return visc_flux; }
const zdouble* Viscous::JacI(void) const { return jac_i; }
const zdouble* Viscous::JacJ(void) const { return jac_j; }
const zdouble* Viscous::StressTensor(void) const { return tau; }

void Viscous::SetStates(const zdouble *v_i, const zdouble *v_j, const zdouble *norm,
        const zdouble area, const zdouble edge_len, const zdouble *vel_grad_i, const zdouble *vel_grad_j) { 
    for (auto iVar = 0ul; iVar < nVar; iVar++) {
        this->v_i[iVar] = v_i[iVar];
        this->v_j[iVar] = v_j[iVar];
    }

    for (auto iDim = 0ul; iDim < nDim; iDim++) {
        this->norm[iDim] = norm[iDim];
        for (auto jDim = 0ul; jDim < nDim; jDim++) {
            this->vel_grad_i[iDim * nDim + jDim] = vel_grad_i[iDim * nDim + jDim];
            this->vel_grad_j[iDim * nDim + jDim] = vel_grad_j[iDim * nDim + jDim];
        }
    }

    this->area = area;
    this->edge_len = edge_len;
}

void Viscous::ComputeResidual() { 
    for (auto iVar = 0ul; iVar < nVar; iVar++) {
        v_mean[iVar] = 0.5 * (v_i[iVar] + v_j[iVar]);
    }

    zdouble t_mean = v_mean[nVar-1] / (v_mean[0] * 287.0);
    zdouble t_ref = 273.15;
    zdouble suth = 110.4;
    zdouble mu_ref = 1.716e-05;
    zdouble t_nondim = t_mean / t_ref;
    //viscosity = mu_ref * t_nondim * sqrtf(t_nondim * (t_ref + suth) / (t_mean + suth));//(t_nondim * ((t_ref + suth) / (t_mean _ suth)));


    for (auto iDim = 0ul; iDim < nDim; iDim++) {
        for (auto jDim = 0ul; jDim < nDim; jDim++) {
            vel_grad_mean[iDim * nDim + jDim] = 0.5 * (vel_grad_i[iDim * nDim + jDim]
                    + vel_grad_j[iDim * nDim + jDim]);
        }
    }

    ComputeStressTensor(vel_grad_mean, viscosity, tau);
    ComputeViscousFlux();

    if (implicit) {
        ComputeTauJacobian();
        ComputeJacobians();
    }
}

void Viscous::ComputeStressTensor(const zdouble *vel_grad, const zdouble viscosity, zdouble *tau ) { 
    for (auto i = 0ul; i < nDim * nDim; i++) tau[i] = 0.0;

    zdouble vel_div = 0.0;
    for (auto iDim = 0ul; iDim < nDim; iDim++)
        vel_div += vel_grad[iDim * nDim + iDim];

    zdouble p = 2.0 / 3.0 * vel_div * viscosity;

    for (auto iDim = 0ul; iDim < nDim; iDim++) {
        for (auto jDim = 0ul; jDim < nDim; jDim++) {
            tau[iDim * nDim + jDim] = viscosity 
                * (vel_grad[iDim * nDim + jDim] + vel_grad[jDim * nDim + iDim]);
        }
        tau[iDim * nDim + iDim] -= p;
    }
}

void Viscous::ComputeViscousFlux(void) {
    zdouble flux_tensor[5][3];

#define TAU(I,J) tau[I * nDim + J]
    if (nDim == 3) {
    }
    else {
        flux_tensor[0][0] = 0.0;
        flux_tensor[1][0] = TAU(0,0);
        flux_tensor[2][0] = TAU(0,1);
        flux_tensor[3][0] = TAU(0,0) * v_mean[1] + TAU(0,1) * v_mean[2];

        flux_tensor[0][1] = 0.0;
        flux_tensor[1][1] = TAU(1,0);
        flux_tensor[2][1] = TAU(1,1);
        flux_tensor[3][1] = TAU(1,0) * v_mean[1] + TAU(1,1) * v_mean[2];
    }

#undef TAU

    for (auto iVar = 0ul; iVar < nVar; iVar++) {
        visc_flux[iVar] = 0.0;
        for (auto iDim = 0ul; iDim < nDim; iDim++) {
            visc_flux[iVar] += flux_tensor[iVar][iDim] * norm[iDim] * area;
        }
    }
}

void Viscous::ComputeTauJacobian(void) {
    zdouble rho = v_mean[0];
    zdouble xi = viscosity / (rho * edge_len);

#define DEL(I,J) I == J ? zdouble{1.0} : zdouble{0.0}
#define TAU(I,J) tau_jac[I * nVar + J]

    for (auto iDim = 0ul; iDim < nDim; iDim++) {
        for (auto jDim = 0ul; jDim < nDim; jDim++) {
            TAU(iDim, jDim+1) = -xi * (DEL(iDim, jDim) + norm[iDim] * norm[jDim] / zdouble{3.0});
        }

        TAU(iDim, 0) = 0.0;
        for (auto jDim = 0ul; jDim < nDim; jDim++) {
            TAU(iDim, 0) -= TAU(iDim, jDim+1) * v_mean[jDim+1];
        }

        TAU(iDim, nVar-1) = 0.0;
    }

#undef DEL
#undef TAU
}

/* Currently no heat fluxes */
void Viscous::ComputeJacobians(void) {
    zdouble rho = v_mean[0];
    zdouble factor = 0.5 / rho;

    if (nDim == 3) {
    }
    else {
        jac_i[0] = 0.0;
        jac_i[1] = 0.0;
        jac_i[2] = 0.0;
        jac_i[3] = 0.0;
        jac_i[4] = area * tau_jac[0];
        jac_i[5] = area * tau_jac[1];
        jac_i[6] = area * tau_jac[2];
        jac_i[7] = area * tau_jac[3];
        jac_i[8] = area * tau_jac[4];
        jac_i[9] = area * tau_jac[5];
        jac_i[10] = area * tau_jac[6];
        jac_i[11] = area * tau_jac[7];
        zdouble contraction = tau_jac[0] * v_mean[1] + tau_jac[4] * v_mean[2];
        jac_i[12] = area * contraction;
        jac_i[13] = -area * tau_jac[0];
        jac_i[14] = -area * tau_jac[4];
        jac_i[15] = 0.0;

        for (auto i = 0ul; i < nVar * nVar; i++) jac_j[i] = -jac_i[i];

        zdouble proj_visc_flux_vel = visc_flux[1] * v_mean[1] + visc_flux[2] * v_mean[2];

        jac_i[12] -= factor * proj_visc_flux_vel;
        jac_i[13] += factor * visc_flux[1];
        jac_i[14] += factor * visc_flux[2];

        jac_j[12] -= factor * proj_visc_flux_vel;
        jac_j[13] += factor * visc_flux[1];
        jac_j[14] += factor * visc_flux[2];
    }
}
