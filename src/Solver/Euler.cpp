#include "Solver.hpp"
#include <fstream>

Euler::Euler() : Solver(Config::SolverType::Euler) { }
Euler::~Euler() { }

void Euler::Solve() { 
    const vector<unique_ptr<Boundary>> &boundaries = grid->Boundaries();
    const vector<Edge> &edges = grid->Edges();

    // bump
    /*
    boundaries[0]->SetType(Boundary::BCType::WallInviscid);
    boundaries[1]->SetType(Boundary::BCType::Farfield);
    boundaries[2]->SetType(Boundary::BCType::Farfield);
    boundaries[3]->SetType(Boundary::BCType::Farfield);
    */

    // cyl
    /*
    boundaries[0]->SetType(Boundary::BCType::WallInviscid);
    boundaries[1]->SetType(Boundary::BCType::Farfield);
    boundaries[2]->SetType(Boundary::BCType::Farfield);
    boundaries[3]->SetType(Boundary::BCType::Farfield);
    boundaries[4]->SetType(Boundary::BCType::Farfield);
    boundaries[5]->SetType(Boundary::BCType::Farfield);
    boundaries[6]->SetType(Boundary::BCType::Farfield);
    */

    // ramp
    /*
    boundaries[0]->SetType(Boundary::BCType::Farfield);
    boundaries[1]->SetType(Boundary::BCType::WallInviscid);
    */

    // su2 wedge 
    /*
    boundaries[0]->SetType(Boundary::BCType::SupersonicInlet);
    boundaries[1]->SetType(Boundary::BCType::WallInviscid);
    boundaries[2]->SetType(Boundary::BCType::SupersonicOutlet);
    boundaries[3]->SetType(Boundary::BCType::WallInviscid);
    */

    // su2 ramp
    /*
    boundaries[0]->SetType(Boundary::BCType::SupersonicInlet);
    boundaries[1]->SetType(Boundary::BCType::WallInviscid);
    boundaries[2]->SetType(Boundary::BCType::SubsonicOutlet);
    boundaries[3]->SetType(Boundary::BCType::WallInviscid);
    */

    // su2 naca airfoil
    /*
    boundaries[0]->SetType(Boundary::BCType::WallViscous);
    boundaries[1]->SetType(Boundary::BCType::Farfield);
    */

    // su2 flat plate
    boundaries[0]->SetType(Boundary::BCType::Farfield);
    boundaries[1]->SetType(Boundary::BCType::Inlet);
    boundaries[2]->SetType(Boundary::BCType::Outlet);
    boundaries[3]->SetType(Boundary::BCType::WallInviscid);
    boundaries[4]->SetType(Boundary::BCType::WallViscous);

    // airfoil
    /*
    boundaries[0]->SetType(Boundary::BCType::WallViscous);
    boundaries[1]->SetType(Boundary::BCType::Farfield);
    */

    // cube
    /*
    boundaries[0]->SetType(Boundary::BCType::Farfield);
    boundaries[1]->SetType(Boundary::BCType::Farfield);
    boundaries[2]->SetType(Boundary::BCType::WallInviscid);
    boundaries[3]->SetType(Boundary::BCType::WallInviscid);
    boundaries[4]->SetType(Boundary::BCType::Farfield);
    boundaries[5]->SetType(Boundary::BCType::Farfield);
    */

    // wing
    /*
    boundaries[0]->SetType(Boundary::BCType::WallInviscid);
    boundaries[1]->SetType(Boundary::BCType::WallInviscid);
    boundaries[2]->SetType(Boundary::BCType::Farfield);
    */

    // sphere 
    /*
    boundaries[0]->SetType(Boundary::BCType::WallInviscid);
    boundaries[1]->SetType(Boundary::BCType::Farfield);
    boundaries[2]->SetType(Boundary::BCType::Farfield);
    */


    unsigned long maxIter = config->GetMaxIterations();
    unsigned long min_iter = config->GetMinIterations();
    unsigned long currentIter = 0;

    while ((!converged && currentIter < min_iter) || currentIter < maxIter) {

        /*Reset wave speeds, residuals and jacobians. Recompute gradients */
        for (unsigned long i = 0; i < nEqn; i++) waveSpeed[i] = 0;
        residual.SetZeroes();
        if (implicit) {
            jacMat.Zeroes();
            delU.SetZeroes();
        }
        if (muscl || viscous) lsq_grad.ComputeGradients(grid, &primVar);

        /* Interior Fluxes */
        for (const Edge &e : edges) {
            const unsigned long &iNodeI = e.Node1();
            const unsigned long &iNodeJ = e.Node2();
            const double edge_len = e.Length();
            const vector<double> &unit_norm = e.AreaVector();
            const double &area = e.Area();

            primVar.CopyBlock(iNodeI, v_i);
            primVar.CopyBlock(iNodeJ, v_j);

            if (muscl) {
                lsq_grad.Reconstruct(v_i, v_j, iNodeI, iNodeJ, edge_len, unit_norm.data()) ;
            }

            conv_flux->SetStates(v_i, v_j, unit_norm.data(), area);
            conv_flux->ComputeFlux();

            waveSpeed[iNodeI] += conv_flux->MaxWaveSpeed();
            waveSpeed[iNodeJ] += conv_flux->MaxWaveSpeed();
            residual.SubtractBlock(iNodeI, conv_flux->Flux());
            residual.AddBlock(iNodeJ, conv_flux->Flux());
            if (implicit) {
                jacMat.AddBlock(iNodeI, iNodeI, conv_flux->JacI());
                jacMat.AddBlock(iNodeI, iNodeJ, conv_flux->JacI());
                jacMat.SubtractBlock(iNodeJ, iNodeJ, conv_flux->JacJ());
                jacMat.SubtractBlock(iNodeJ, iNodeI, conv_flux->JacJ());
            }

            if (viscous) {
                const double* vel_grad_i = lsq_grad.GetGradientConst(iNodeI, 1);
                const double *vel_grad_j = lsq_grad.GetGradientConst(iNodeJ, 1);
                visc_flux.SetStates(v_i, v_j, unit_norm.data(), area, edge_len, vel_grad_i, vel_grad_j);
                visc_flux.ComputeResidual();
                residual.AddBlock(iNodeI, visc_flux.Flux());
                residual.SubtractBlock(iNodeJ, visc_flux.Flux());
                if (implicit) {
                    jacMat.SubtractBlock(iNodeI, iNodeI, visc_flux.JacI());
                    jacMat.SubtractBlock(iNodeI, iNodeJ, visc_flux.JacI());
                    jacMat.AddBlock(iNodeJ, iNodeJ, visc_flux.JacJ());
                    jacMat.AddBlock(iNodeJ, iNodeI, visc_flux.JacJ());
                }
            }
        }

        /* Boundary Fluxes */
        const vector<unique_ptr<Boundary>> &boundaries = grid->Boundaries();
        for (const auto &bndry : boundaries) {
            const vector<unsigned long> &iNodes = bndry->Nodes();
            const vector<vector<double>> &norms = bndry->NodeNorms();
            const vector<double> &area = bndry->DualAreas();
            unsigned long nNodes = iNodes.size();

            for (auto i = 0ul; i < nNodes; i++) {
                const unsigned long &iNode = iNodes[i];
                primVar.CopyBlock(iNode, v_i);
                conVar.CopyBlock(iNode, v_j);

                if (bndry->Type() == Boundary::BCType::WallViscous) {
                    double *res = residual.GetBlock(iNode);
                    double *jac = jacMat.GetBlock(iNode, iNode);
                    for (size_t iDim = 0; iDim < nDim; iDim++) {
                        v_i[iDim+1] = 0.0;
                        v_j[iDim+1] = 0.0;
                        res[iDim+1] = 0.0;
                    }
                    for (size_t iVar = 1; iVar < nVar-1; iVar++) {
                        for (size_t jVar = 0; jVar < nVar; jVar++) {
                            jac[iVar * nVar + jVar] = iVar == jVar ? 1.0 : 0.0;
                        }
                    }
                    primVar.SetBlock(iNode, v_i);
                    conVar.SetBlock(iNode, v_j);

                }
                else {

                    vector<double> v_j = GetBoundaryState(bndry->Type(), v_i, norms[i]);

                    conv_flux->SetStates(v_i, v_j.data(), norms[i].data(), area[i]);
                    conv_flux->ComputeFlux();

                    waveSpeed[iNode] += conv_flux->MaxWaveSpeed();
                    residual.SubtractBlock(iNode, conv_flux->Flux());
                    if (implicit) jacMat.AddBlock(iNode, iNode, conv_flux->JacI());
                }
            }
        }

        /* Add psuedo time term or update explicit */
        const vector<Node> &nodes = grid->Nodes();
        for (auto i = 0ul; i < nEqn; i++) {
            const double &vol = nodes[i].DualVolume();
            dt[i] = cfl * vol / waveSpeed[i]; 

            if (implicit) jacMat.AddToDiag(i, i, vol / dt[i]);
            else {
                const double * const res = residual.GetBlock(i);
                conVar.AddBlock(i, res, dt[i] / vol);
            }
        }

        /* Solve and Update */
        if (implicit) {
            precond->Build(jacMat);
            gmres.Solve(jacMat, residual, delU, precond);
            conVar.AddVector(delU);
        }

        PrintResiduals(currentIter++);

        UpdatePrimitiveVars();

        if (config->AdaptiveCFL()) AdaptCFL();
    }
}

vector<double> Euler::GetBoundaryState(const Boundary::BCType &type, 
        const double *leftState, const vector<double> &norm) const {
    vector<double> bState(nVar, 0);
    if (type == Boundary::BCType::Farfield) {
        const vector<double> &vel_inf = config->GetFreestreamVelocity();
        const double &rho_inf = config->GetFreestreamDensity();
        const double &p_inf = config->GetFreestreamPressure();

        bState[0] = rho_inf;
        for (auto i = 0ul; i < nDim; i++) {
            bState[i+1] = vel_inf[i];
        }
        bState[nVar-1] = p_inf;

    }
    else if (type == Boundary::BCType::WallInviscid) {
        bState[0] = leftState[0];

        double vel_norm = 0;
        for (size_t i = 0; i < nDim; i++) { 
            bState[i+1] = leftState[i+1];
            vel_norm += leftState[i+1] * norm[i]; 
        }
        for (size_t i = 0; i < nDim; i++) 
            bState[i+1] -= vel_norm * norm[i];

        bState[nVar-1] = leftState[nVar-1];
    }
    else if (type == Boundary::BCType::WallViscous) {
        bState[0] = leftState[0];
        for (size_t i = 0; i < nDim; i++)
            bState[i+1] = 0.0;

        bState[nVar-1] = leftState[nVar-1];
    }
    else if (type == Boundary::BCType::SupersonicInlet) {
        const vector<double> &vel_inf = config->GetFreestreamVelocity();

        bState[0] = config->GetFreestreamDensity();
        for (size_t i = 0; i < nDim; i++)
            bState[i+1] = vel_inf[i];
        bState[nVar-1] = config->GetFreestreamPressure();
    }
    else if (type == Boundary::BCType::SupersonicOutlet) {
        for (unsigned short i = 0; i < nVar; i++) bState[i] = leftState[i];
    }
    else if (type == Boundary::BCType::Inlet) {

        double v_2 = leftState[1] * leftState[1] + leftState[2] * leftState[2];
        double soundspeed = sqrt(gamma * leftState[nVar-1] / leftState[0]);
        double machspeed = sqrt(v_2) / soundspeed;
        const vector<double> &vel_inf = config->GetFreestreamVelocity();

        bState[0] = config->GetFreestreamDensity();
        for (size_t i = 0; i < nDim; i++) bState[i+1] = vel_inf[i];

        if (machspeed >= 1.0) {
            bState[nVar-1] = config->GetFreestreamPressure();
        }
        else {
            bState[nVar-1] = leftState[nVar-1];
        }
    }
    else if (type == Boundary::BCType::Outlet) {
        double v_2 = leftState[1] * leftState[1] + leftState[2] * leftState[2];
        double soundspeed = sqrt(gamma * leftState[nVar-1] / leftState[0]);
        double machspeed = sqrt(v_2) / soundspeed;

        if (machspeed >= 1.0) {
            for (unsigned short i = 0; i < nVar; i++) bState[i] = leftState[i];
        }
        else {
            for (unsigned short i = 0; i < nVar; i++) bState[i] = leftState[i];
            bState[nVar-1] = config->GetFreestreamPressure();
        }
    }
    else {
    }

    return bState;
}
