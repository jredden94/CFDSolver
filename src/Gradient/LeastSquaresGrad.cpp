#include "LeastSquaresGrad.hpp"

LeastSquaresGrad::LeastSquaresGrad() { 
    Config *config = &Config::GetConfig();
    limiter.SetLimiter(config->GetLimiterType());
    is_limited = config->IsLimited();
    nEqn = config->GetNumEqn();
    nDim = config->GetNumDims();
    nVar = config->GetNumVars();
    blk_len = nVar * nDim;
    arr_len = nEqn * blk_len;
    grad_array = new double[arr_len];
    epsilon = limiter.GetEpsilon();
}

LeastSquaresGrad::~LeastSquaresGrad() { }

void LeastSquaresGrad::Reset() {
    for (auto i = 0ul; i < arr_len; i++) grad_array[i] = 0;
}

void LeastSquaresGrad::ComputeGradients(const Grid *grid, const SolVector *primVar) {
    Reset();
    const vector<Node> &nodes = grid->Nodes();
    unsigned long nNodes = nodes.size();

    if (nDim == 3) {
        for (auto iNode = 0ul; iNode < nNodes; iNode++) {
            const Node &node_i = nodes[iNode];
            SolveLeastSquares3D(iNode, node_i, nodes, primVar);
        }
    }
    else {
        for (auto iNode = 0ul; iNode < nNodes; iNode++) {
            const Node &node_i = nodes[iNode];
            SolveLeastSquares2D(iNode, node_i, nodes, primVar);
        }
    }
}

void LeastSquaresGrad::SolveLeastSquares2D(const unsigned long iNode, const Node &node_i, 
        const vector<Node> &nodes, const SolVector *primVar) {
    const double *v_i = primVar->GetBlock(iNode);
    const vector<unsigned long> &nbrs = node_i.Neighbors();
    const unsigned short nNbrs = nbrs.size();
    a11 = a12 = a22 = 0;

    for (auto iNbr : nbrs) {
        const Node &node_j = nodes[iNbr]; 
        const double *v_j = primVar->GetBlock(iNbr);
        /* A^T A */
        dx = node_j.X() - node_i.X();
        dy = node_j.Y() - node_i.Y();
        a11 += dx * dx;
        a12 += dx * dy;
        a22 += dy * dy;

        /* A^T b */
        for (auto iVar = 0ul; iVar < nVar; iVar++) {
            double *grad = GetGradient(iNode, iVar);
            delta = v_j[iVar] - v_i[iVar];
            grad[0] += dx * delta;
            grad[1] += dy * delta;
        }
    }

    /* inv(A^T A) */
    det_inv = 1.0 / (a11 * a22 - a12 * a12 + epsilon);
    temp = a11;
    a11 = a22 * det_inv;
    a22 = temp * det_inv;
    a12 = -a12 * det_inv;

    /* grad(phi) = inv(A^T A) * (A^T b)  */
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        double *grad = GetGradient(iNode, iVar);
        temp = grad[0];
        grad[0] = a11 * grad[0] + a12 * grad[1];
        grad[1] = a12 * temp + a22 * grad[1];
    }
}

void LeastSquaresGrad::SolveLeastSquares3D(const unsigned long iNode, const Node &node_i, 
        const vector<Node> &nodes, const SolVector *primVar) {
    const double *v_i = primVar->GetBlock(iNode);
    const vector<unsigned long> &nbrs = node_i.Neighbors();
    const unsigned short nNbrs = nbrs.size();
    a11 = a12 = a13 = a22 = a23 = a33 = 0;

    for (auto iNbr : nbrs) {
        const Node &node_j = nodes[iNbr]; 
        const double *v_j = primVar->GetBlock(iNbr);
        /* A^T A */
        dx = node_j.X() - node_i.X();
        dy = node_j.Y() - node_i.Y();
        dz = node_j.Z() - node_i.Z();
        a11 += dx * dx;
        a12 += dx * dy;
        a13 += dx * dz;
        a22 += dy * dy;
        a23 += dy * dz;
        a33 += dz * dz;

        /* A^T b */
        for (auto iVar = 0ul; iVar < nVar; iVar++) {
            double *grad = GetGradient(iNode, iVar);
            delta = v_j[iVar] - v_i[iVar];
            grad[0] += dx * delta;
            grad[1] += dy * delta;
            grad[2] += dz * delta;
        }
    }

    /* inv(A^T A) */
    cofac11 = a22 * a33 - a23 * a23;
    cofac12 = -(a12 * a33 - a23 * a13);
    cofac13 = a12 * a23 - a22 * a13;
    cofac22 = a11 * a33 - a13 * a13;
    cofac23 = -(a11 * a23 - a12 * a13);
    cofac33 = a11 * a22 - a12 * a12;
    det_inv = 1 / (a11 * cofac11 + a12 * cofac12 + a13 * cofac13 + epsilon);
    a11 = cofac11 * det_inv;
    a12 = cofac12 * det_inv;
    a13 = cofac13 * det_inv;
    a22 = cofac22 * det_inv;
    a23 = cofac23 * det_inv;
    a33 = cofac33 * det_inv;

    /* grad(phi) = inv(A^T A) * (A^T b)  */
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        double *grad = GetGradient(iNode, iVar);
        temp = grad[0];
        temp2 = grad[1];
        grad[0] = a11 * grad[0] + a12 * grad[1] + a13 * grad[2];
        grad[1] = a12 * temp + a22 * grad[1] + a23 * grad[2];
        grad[2] = a13 * temp + a23 * temp2 + a33 * grad[2];
    }
}

const double* LeastSquaresGrad::GetGradientConst(const unsigned long iNode, const unsigned short iVar) const { 
    return &grad_array[iNode * blk_len + nDim * iVar];
}

double * LeastSquaresGrad::GetGradient(const unsigned long iNode, const unsigned short iVar) { 
    return &grad_array[iNode * blk_len + nDim * iVar];
}

void LeastSquaresGrad::Reconstruct(double *v_i, double *v_j, const unsigned long iNodeI, 
        const unsigned long iNodeJ, const double edge_len, const double *unit_norm) const {
    const double &epsilon = limiter.GetEpsilon();

    double tmp_vi[nVar], tmp_vj[nVar];
    for (size_t i = 0; i < nVar; i++) {
        tmp_vi[i] = v_i[i];
        tmp_vj[i] = v_j[i];
    }

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        const double *grad_i = GetGradientConst(iNodeI, iVar);
        const double *grad_j = GetGradientConst(iNodeJ, iVar);
        double proj_grad_i = 0, proj_grad_j = 0;
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
            proj_grad_i += grad_i[iDim] * unit_norm[iDim];
            proj_grad_j += grad_j[iDim] * unit_norm[iDim];
        }
        double del_v = v_j[iVar] - v_i[iVar];

        double lim_i = 1.0;
        double lim_j = 1.0;
        if (is_limited) {
            lim_i = limiter.Limit(proj_grad_i, del_v, epsilon);
            lim_j = limiter.Limit(proj_grad_j, del_v, epsilon);
        }

        v_i[iVar] += 0.5 * edge_len * proj_grad_i * lim_i;
        v_j[iVar] -= 0.5 * edge_len * proj_grad_j * lim_j;
    }

    if (v_i[0] <= 0.0 || v_i[nVar-1] <= 0.0 || v_j[0] <= 0.0 || v_j[nVar-1] <= 0.0) {
        for (size_t i = 0; i < nVar; i++) {
            v_i[i] = tmp_vi[i];
            v_j[i] = tmp_vj[i];
        }
    }
}
