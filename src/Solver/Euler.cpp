#include "Solver.hpp"
#include "../Common/AD.hpp"
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
    boundaries[0]->SetType(Boundary::BCType::WallViscous);
    boundaries[1]->SetType(Boundary::BCType::Farfield);
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
    boundaries[2]->SetType(Boundary::BCType::SupersonicOutlet);
    boundaries[3]->SetType(Boundary::BCType::WallInviscid);
    */

    // su2 naca airfoil
    boundaries[0]->SetType(Boundary::BCType::WallViscous);
    boundaries[1]->SetType(Boundary::BCType::Farfield);

    // su2 flat plate
    /*
    boundaries[0]->SetType(Boundary::BCType::Farfield);
    boundaries[1]->SetType(Boundary::BCType::Inlet);
    boundaries[2]->SetType(Boundary::BCType::Outlet);
    boundaries[3]->SetType(Boundary::BCType::WallInviscid);
    boundaries[4]->SetType(Boundary::BCType::WallViscous);
    */

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

    while ((!converged || currentIter < min_iter) && currentIter < maxIter) {

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
            const zdouble edge_len = e.Length();
            const vector<zdouble> &area_norm = e.AreaVector();
            const vector<zdouble> &edge_norm = e.EdgeVector();
            const zdouble &area = e.Area();

            primVar.CopyBlock(iNodeI, v_i);
            primVar.CopyBlock(iNodeJ, v_j);

            if (muscl) {
                lsq_grad.Reconstruct(v_i, v_j, iNodeI, iNodeJ, edge_len, edge_norm.data()) ;
            }

            conv_flux->SetStates(v_i, v_j, area_norm.data(), area);
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
                const zdouble* vel_grad_i = lsq_grad.GetGradientConst(iNodeI, 1);
                const zdouble *vel_grad_j = lsq_grad.GetGradientConst(iNodeJ, 1);
                visc_flux.SetStates(v_i, v_j, area_norm.data(), area, edge_len, vel_grad_i, vel_grad_j);
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
            const vector<vector<zdouble>> &norms = bndry->NodeNorms();
            const vector<zdouble> &area = bndry->DualAreas();
            unsigned long nNodes = iNodes.size();

            for (auto i = 0ul; i < nNodes; i++) {
                const unsigned long &iNode = iNodes[i];
                primVar.CopyBlock(iNode, v_i);
                conVar.CopyBlock(iNode, v_j);

                if (bndry->Type() == Boundary::BCType::WallViscous) {
                    zdouble *res = residual.GetBlock(iNode);
                    zdouble *jac = jacMat.GetBlock(iNode, iNode);
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

                    vector<zdouble> v_j = GetBoundaryState(bndry->Type(), v_i, norms[i]);

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
            const zdouble &vol = nodes[i].DualVolume();
            dt[i] = cfl * vol / waveSpeed[i]; 

            if (implicit) jacMat.AddToDiag(i, i, vol / dt[i]);
            else {
                const zdouble * const res = residual.GetBlock(i);
                conVar.AddBlock(i, res, dt[i] / vol);
            }
        }

        /* Solve and Update */
        if (implicit) {
            precond->Build(jacMat);
            gmres.Solve(jacMat, residual, delU, precond);
            conVar.AddVector(delU);
        }


        UpdatePrimitiveVars();

        zdouble F_x = 0.0, F_y = 0.0;
        zdouble p_inf = config->GetFreestreamPressure();
        for (const auto &bound : boundaries) {
            if (bound->Type() == Boundary::BCType::WallViscous) {
                auto wall_nodes = bound->Nodes();
                auto wall_node_norm = bound->NodeNorms();
                for (auto i = 0ul; i < wall_nodes.size(); i++) {

                    unsigned long iNode = bound->Nodes()[i];
                    vector<zdouble> &norm = wall_node_norm[i];
                    zdouble area = bound->DualAreas()[i];
                    const zdouble *v_i = primVar.GetBlock(iNode);
                    F_x += -(v_i[nVar-1] - p_inf) * area * norm[0];
                    F_y += -(v_i[nVar-1] - p_inf) * area * norm[1];
                    if (viscous) {
                        zdouble tau[nDim * nDim];
                        const zdouble *vel_grad = lsq_grad.GetGradientConst(iNode, 1);

                        zdouble t_mean = v_i[nVar-1] / (v_i[0] * 287.0);
                        zdouble t_ref = 273.15;
                        zdouble suth = 110.4;
                        zdouble mu_ref = 1.716e-05;
                        zdouble t_nondim = t_mean / t_ref;
                        zdouble viscosity = config->GetViscosity(); // mu_ref * t_nondim * sqrtf(t_nondim * (t_ref + suth) / (t_mean + suth));//(t_nondim * ((t_ref + suth) / (t_mean _ suth)));

                        visc_flux.ComputeStressTensor(vel_grad, viscosity, tau);
                        zdouble tau_elm[nDim];
                        for (size_t iDim = 0; iDim < nDim; iDim++) {
                            tau_elm[iDim] = 0.0;
                            for (size_t jDim = 0; jDim < nDim; jDim++) {
                                tau_elm[iDim] += tau[iDim * nDim + jDim] * norm[jDim];
                            }
                        }

                        F_x += tau_elm[0] * area;
                        F_y += tau_elm[1] * area;
                        //drag -= area * (tau[1] * norm[1]);
                        //lift -= area * (tau[1] * norm[0]);
                    }
                }
            }
            else continue;
        }

        const vector<zdouble> &free_vel = config->GetFreestreamVelocity();
        const zdouble rho_inf = config->GetFreestreamDensity();
        zdouble vel_2 = 0.0;
        for (size_t iDim = 0; iDim < nDim; iDim++) {
            vel_2 += free_vel[iDim] * free_vel[iDim];
        }

        zdouble factor = 1. / (0.5 * rho_inf * vel_2 * config->GetReferenceArea());

        zdouble cL, cD, aoa;
        aoa = config->GetAngleOfAttack();
        aoa = M_PI * aoa / 180.;
        zdouble ref_area = config->GetReferenceArea();
        zdouble drag = factor * (F_x * cosf(aoa) + F_y * sinf(aoa));
        zdouble lift = factor * (-F_x * sinf(aoa) + F_y * cosf(aoa));

        //cout << "Lift: " << lift << "\tDrag: " << drag << endl;
        PrintResiduals(currentIter++, lift, drag);
        if (config->AdaptiveCFL()) AdaptCFL();
    }

    ofstream f("wall_values.txt");
    const auto &nodes = grid->Nodes();
    f << "X\tY\tPressure\tn_x\tn_y\tarea\n";
    for (const auto &bound : boundaries) {
        if (bound->Type() == Boundary::BCType::WallViscous) {
            auto wall_nodes = bound->Nodes();
            auto wall_node_norm = bound->NodeNorms();
            for (auto i = 0ul; i < wall_nodes.size(); i++) {
                vector<zdouble> &norm = wall_node_norm[i];
                zdouble area = bound->DualAreas()[i];
                unsigned long iNode = bound->Nodes()[i];
                const zdouble *v_i = primVar.GetBlock(iNode);
                const Node &node = nodes[iNode];

                f << node.X() << "\t" << node.Y() << "\t" << v_i[nVar-1] << "\t" << norm[0] << "\t" << norm[1] << "\t" << area << endl;
            }
        }
        else continue;
    }
    f.close();
}

vector<zdouble> Euler::GetBoundaryState(const Boundary::BCType &type, 
        const zdouble *leftState, const vector<zdouble> &norm) const {
    vector<zdouble> bState(nVar, 0);
    if (type == Boundary::BCType::Farfield) {
        const vector<zdouble> &vel_inf = config->GetFreestreamVelocity();
        const zdouble &rho_inf = config->GetFreestreamDensity();
        const zdouble &p_inf = config->GetFreestreamPressure();

        bState[0] = rho_inf;
        for (auto i = 0ul; i < nDim; i++) {
            bState[i+1] = vel_inf[i];
        }
        bState[nVar-1] = p_inf;

    }
    else if (type == Boundary::BCType::WallInviscid) {
        bState[0] = leftState[0];

        zdouble vel_norm = 0;
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
        const vector<zdouble> &vel_inf = config->GetFreestreamVelocity();

        bState[0] = config->GetFreestreamDensity();
        for (size_t i = 0; i < nDim; i++)
            bState[i+1] = vel_inf[i];
        bState[nVar-1] = config->GetFreestreamPressure();
    }
    else if (type == Boundary::BCType::SupersonicOutlet) {
        for (unsigned short i = 0; i < nVar; i++) bState[i] = leftState[i];
    }
    else if (type == Boundary::BCType::Inlet) {

        zdouble v_2 = leftState[1] * leftState[1] + leftState[2] * leftState[2];
        zdouble soundspeed = sqrt(gamma * leftState[nVar-1] / leftState[0]);
        zdouble machspeed = sqrt(v_2) / soundspeed;
        const vector<zdouble> &vel_inf = config->GetFreestreamVelocity();

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
        zdouble v_2 = leftState[1] * leftState[1] + leftState[2] * leftState[2];
        zdouble soundspeed = sqrt(gamma * leftState[nVar-1] / leftState[0]);
        zdouble machspeed = sqrt(v_2) / soundspeed;

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
