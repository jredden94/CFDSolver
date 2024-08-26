#pragma once

#include "SolVector.hpp"
#include "../Common/Config.hpp"
#include "BLAS.hpp"

class DiagonalMatrix{
    public:
        DiagonalMatrix();
        ~DiagonalMatrix();

        void Zeroes(void);
        double* GetBlock(const unsigned long iBlock);
        const double* GetBlock(const unsigned long iBlock) const;
        void SetBlock(const double *block, const unsigned long iBlock);
        void MatrixVecMult(const SolVector &vec, SolVector &prod) const;
    private:
        double *val;
        unsigned long len, nBlk, blk_len, nVar;
};
