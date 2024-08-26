#pragma once

#include "../DualGrid/Grid.hpp"
#include "../LinAlg/SolVector.hpp"
#include "../Common/Config.hpp"
#include "../Limiters/Limiter.hpp"

class LeastSquaresGrad {
    public:
        LeastSquaresGrad();
        ~LeastSquaresGrad();

        void ComputeGradients(const Grid *grid, const SolVector *primVar);

        const double* GetGradientConst(const unsigned long iNode, const unsigned short iVar) const;

        void Reconstruct(double *v_i, double *v_j, const unsigned long iNodeI, const unsigned long iNodeJ, 
                const double edge_len, const double *unit_norm) const;

        void Reset(void);

    private:
        /* Solve (A^T A) grad(phi) = A^T b */
        void SolveLeastSquares2D(const unsigned long iNode, const Node &node_i, const vector<Node> &nodes, const SolVector *primVar);
        void SolveLeastSquares3D(const unsigned long iNode, const Node &node_i, const vector<Node> &nodes, const SolVector *primVar);

        double* GetGradient(const unsigned long iNode, const unsigned short iVar);
        Limiter limiter;
        bool is_limited;

        unsigned long nEqn, arr_len;
        unsigned short nDim, nVar, blk_len;

        double *grad_array;
        double epsilon;
        double dx, dy, dz, det_inv, delta;

        /* (A^T A) and its inverse are symmetric, so only upper triangular entries are computed */
        double a11, a12, a13, a22, a23, a33, temp, temp2;
        double cofac11, cofac12, cofac13, cofac22, cofac23, cofac33;
};
