#pragma once

#include "../DualGrid/Grid.hpp"
#include "../LinAlg/SolVector.hpp"
#include "../Common/Config.hpp"
#include "../Limiters/Limiter.hpp"
#include "../Common/AD.hpp"

class LeastSquaresGrad {
    public:
        LeastSquaresGrad();
        ~LeastSquaresGrad();

        void ComputeGradients(const Grid *grid, const SolVector *primVar);

        const zdouble* GetGradientConst(const unsigned long iNode, const unsigned short iVar) const;

        void Reconstruct(zdouble *v_i, zdouble *v_j, const unsigned long iNodeI, const unsigned long iNodeJ, 
                const zdouble edge_len, const zdouble *unit_norm) const;

    private:
        /* Solve (A^T A) grad(phi) = A^T b */
        void SolveLeastSquares2D(const unsigned long iNode, const Node &node_i, const vector<Node> &nodes, const SolVector *primVar);
        void SolveLeastSquares3D(const unsigned long iNode, const Node &node_i, const vector<Node> &nodes, const SolVector *primVar);

        void Reset(void);
        zdouble* GetGradient(const unsigned long iNode, const unsigned short iVar);
        Limiter limiter;
        bool is_limited;

        unsigned long nEqn, arr_len;
        unsigned short nDim, nVar, blk_len;

        zdouble *grad_array;
        zdouble epsilon, machine_eps;
        zdouble dx, dy, dz, det_inv, delta;

        /* (A^T A) and its inverse are symmetric, so only upper triangular entries are computed */
        zdouble a11, a12, a13, a22, a23, a33, temp, temp2;
        zdouble cofac11, cofac12, cofac13, cofac22, cofac23, cofac33;
};
